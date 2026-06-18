    // ============================================================================
    // File: nps_sim_smearing_new_try.C
    // Purpose: NPS π0 simulation smearing analysis with per-section optimization
    //
    // This code performs chi2 minimization to find optimal smearing parameters
    // (mu, sigma, sigma_pos) for each section of the NPS calorimeter by comparing
    // data and simulation distributions (M_γγ, M_miss, and (p_target+γγ)^2).
    //
    // Uses physics calculations exactly as in nps_analysis.C
    //
    // INPUT REQUIREMENTS:
    //   - Simulation data should be pre-filtered by simc_pi0_analysis.C
    //   - Pre-filtering ensures HMS cuts and event selection are consistent
    //   - This script does NOT apply additional HMS or selection cuts
    //
    // Features:
    //   - Per-section optimization of energy scale (mu), resolution (sigma), and position smearing (sigma_pos)
    //   - Shared smearing implementation across objective, diagnostics, and combined summary pages
    //   - Optional final global p_e_scale fit from M_miss-only objective
    //   - Bilinear interpolation for smooth parameter variation across calorimeter
    //   - Outputs: discrete section parameters (CSV) + interpolated 2D maps (ROOT)
    //
    // Compile:
    //   g++ nps_sim_smearing_new_try.C `root-config --cflags --libs` -lMathMore -O2 -std=c++17 -fopenmp -I../src -o nps_sim_smearing_new_try
    //
    // Usage example: 
    //   ./nps_sim_smearing_new_try /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/combined_branches_LH2.root physics /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/simc_pi0_analysis_output.root simulation out_smear.root 11 16 -24 28 -34 34 0.2 50
    //
    //   ./nps_sim_smearing_new_try /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/production_wfpi0/combined_branches_LH2_wfpi0.root physics /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/production_wfpi0/simc_pi0_analysis_output.root simulation out_smear.root 11 16 -24 28 -34 34 0.2 50
    // 
    // Output files:
    //   - out.root: Per-section histograms and fit parameters
    //   - out_interpolated.root: 2D maps (h_mu_interp, h_sigma_interp, h_sigma_pos_interp)
    //   - section_map.csv: Tabulated calibration parameters
    //   - chi2_scans.pdf: Chi^2 scan plots for each section
    //
    // ============================================================================

    // Include standard library headers FIRST
    #include <iostream>
    #include <vector>
    #include <string>
    #include <fstream>
    #include <cmath>
    #include <limits>
    #include <algorithm>
    #include <sstream>
    #include <utility>
    #include <iomanip>
    #include <set>
    #include <ctime>
    #include <chrono>
    #include <cstdlib>
    #include <memory>
    #include <array>
    #include <map>
    #include <cctype>
    #include <sys/stat.h>
    #include <omp.h>

    // Bring standard library types into global namespace for analysis helpers
    using std::vector;
    using std::string;
    using std::pair;
    using std::cout;
    using std::endl;
    using std::cerr;
    using std::ifstream;
    using std::ofstream;
    using std::ostringstream;
    using std::istringstream;

    // Include ROOT headers SECOND
    #include "TFile.h"
    #include "TTree.h"
    #include "TH1.h"
    #include "TH1D.h"
    #include "TH2D.h"
    #include "TF1.h"
    #include "TCanvas.h"
    #include "TLatex.h"
    #include "TLegend.h"
    #include "TLine.h"
    #include "TPaveText.h"
    #include "TStyle.h"
    #include "TRandom3.h"
    #include "TLorentzVector.h"
    #include "TVector3.h"
    #include "TNamed.h"
    #include "TMultiGraph.h"
    #include "TDirectory.h"
    #include "TSystem.h"
    #include "TFitResult.h"
    #include "TGraph.h"
    #include "TGraph2D.h"
    #include "TMarker.h"
    #include "TObject.h"
    #include "Math/Factory.h"
    #include "Math/Functor.h"
    #include "Math/Minimizer.h"
    #include "Math/QuasiRandom.h"

    // Include project headers LAST (these depend on std and ROOT)
    #include "../src/utils.C"
    #include "../src/nps_helper.h"
    #include "../src/nps_time_bg.h"
    #include "../src/nps_comb_bg.h"
    #include "../src/nps_mmiss_cor.h"
    #include "../src/nps_hms_track_eff.h"
    #include "../src/nps_yao_database_reader.h"

    // Suppress empty-body warning in physics_var.h
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wempty-body"
    #include "../src/nps_physics_var.h"
    #pragma GCC diagnostic pop

    using namespace std;

    // ============================================================================
    // ============================================================================
    //                         USER CONFIGURATION SECTION
    // ============================================================================
    // All user-configurable parameters are in this section for easy access.
    // Modify these values to customize the analysis behavior.
    //
    // RECOMMENDED CONFIGURATION (for most analyses):
    //   - ENERGY_SMEARING_HISTOGRAM = HIST_BOTH (Use BOTH M_γγ and M_miss)
    //   - ENABLE_POSITION_SMEARING = true       (Calibrate angular resolution)
    //   - W_MPI0 = 1.0, W_MMISS = 1.0           (Equal weights)
    //   - USE_SIMPLE_STOCHASTIC_MODEL = true    (σ_E = σ√E model)
    //
    // NOTE: Current defaults in this file are tuned for production_wfpi0 and may
    //       differ from the above recommended profile.
    //
    // CAUTION: Using ENERGY_SMEARING_HISTOGRAM = 2 (M_miss only) will give
    //          POOR M_γγ agreement. See MISSING_MASS_FITTING_ISSUE.txt for details.
    // ============================================================================

    namespace Config {
        // ========================================================================
        // PHYSICS SETTINGS
        // ========================================================================
        
        // Beam energy for missing mass calculation (GeV)
        const double BEAM_ENERGY = 10.538;  // Hall C typical beam energy
        
        // Y mispointing offset (cm) - must match simc_pi0_analysis.C
        const double Y_MISPOINT = 0.103665;  // SIMC infile for x60_4b spec%e%offset%y
        
        // Clamp value for non-positive energy smearing draws (GeV)
        const double NONPOSITIVE_CLAMP = 1e-6;

        // Optional linear pair-energy correction versus reconstructed m_gg:
        //   E_corr = E_smear * f(m_gg), where
        //   f(m_gg) = 1 / (1 + slope * (m_gg - pivot))  when use_inverse=true
        //           =     (1 + slope * (m_gg - pivot))  when use_inverse=false
        // This is applied equally to both photons after per-photon smearing.
        const bool ENABLE_MGG_LINEAR_ENERGY_CORRECTION = false;
        const double MGG_LINEAR_SLOPE = -40.0;          // GeV^-1
        const double MGG_LINEAR_PIVOT_GEV = 0.135;    // GeV
        const bool MGG_LINEAR_USE_INVERSE = false;
        const double MGG_LINEAR_FACTOR_MIN = 0.85;
        const double MGG_LINEAR_FACTOR_MAX = 1.15;
        
        // ========================================================================
        // ENERGY RESOLUTION MODEL
        // ========================================================================
        // Energy-smearing PDF shape for per-photon response tails.
        //   0: Gaussian  -> sigma interpreted as standard-deviation coefficient
        //   1: Landau    -> sigma interpreted as FWHM coefficient
        constexpr int SMEAR_SHAPE_GAUSSIAN = 0;
        constexpr int SMEAR_SHAPE_LANDAU = 1;
        const int ENERGY_SMEAR_SHAPE = SMEAR_SHAPE_GAUSSIAN;  // RECOMMENDED: Landau for more realistic tails, Gaussian for simpler model and faster computation

        // Conversion used in Landau mode:
        //   FWHM ~= LANDAU_FWHM_TO_SCALE * (Landau scale parameter)
        // The exact factor depends on parameterization; 4.0 is a robust practical default.
        const double LANDAU_FWHM_TO_SCALE = 4.0;

        // Best-practice Landau sampling controls:
        //  - Resample a finite number of times to avoid non-physical tails dominating the fit.
        //  - Enforce E > NONPOSITIVE_CLAMP and (optionally) an upper tail cap.
        // The cap is expressed in units of FWHM above the scaled energy mean:
        //   E_max = E_scaled + LANDAU_MAX_FWHM_ABOVE_MPV * FWHM(E_scaled)
        const bool ENABLE_TRUNCATED_LANDAU = true;
        const int LANDAU_MAX_REDRAWS = 64;
        const double LANDAU_MAX_FWHM_ABOVE_MPV = 8.0;

        // RNG seed used for deterministic smearing in every chi2 evaluation.
        // Resetting the RNG to this seed at the start of each evaluation ensures
        // that the chi2 surface is smooth and reproducible — identical random draws
        // are used for every parameter-point evaluation, eliminating stochastic
        // noise that can trap the grid-search optimizer.
        const unsigned int SMEAR_DETERMINISTIC_SEED = 42;

        // Choose between two resolution models:
        //
        // MODEL 1 (Simple Stochastic): σ_E = σ × √E
        //   - Data-driven, 2 parameters per section (μ, σ)
        //   - σ typically 0.02-0.04 GeV^(1/2) for calorimeters
        //   - Fast, good for empirical calibration
        //   - RECOMMENDED for most users
        //
        // MODEL 2 (Full 3-term): σ_E = σ × E × √((A/√E)² + (B/E)² + C²)
        //   - Physics-motivated: stochastic + noise + constant terms
        //   - Requires A,B,C from external measurements or fit
        //   - σ ≈ 1.0 if A,B,C are correct for your detector
        //   - Use only if you know your detector's resolution parameters
        //
        const bool USE_SIMPLE_STOCHASTIC_MODEL = true;  // true = Model 1, false = Model 2
        
        // Model 2 constants (only used if USE_SIMPLE_STOCHASTIC_MODEL = false)
        // Reference values from Wasim's NIM paper draft (E in GeV, ⊕ means quadrature sum):
        const double RESOLUTION_A_DEFAULT = 0.97;  // Stochastic term starting value (0.97%)
        const double RESOLUTION_B_DEFAULT = 1.1;   // Noise term starting value (1.1%)
        const double RESOLUTION_C_DEFAULT = 1.14;  // Constant term starting value (1.14%)
        
        // ========================================================================
        // ========================================================================
        // OPTIMIZATION STRATEGY SELECTION
        // ========================================================================
        // Default optimizer:
        //   deterministic low-discrepancy bounded seeds -> keep best N seeds
        //   -> Minuit2 MIGRAD from each seed -> HESSE/profile diagnostics.
        
        // CHI-SQUARED OBSERVABLE WEIGHTS
        // ========================================================================
        // Relative weights for combined chi-squared:
        // χ²_total = W_MPI0 × χ²_M_γγ + W_MMISS × χ²_M_miss + W_MPGG2 × χ²_(p+γγ)^2
        //
        // Equal weights (1.0, 1.0) treat both observables equally.
        // Current defaults below are intentionally tuned to emphasize M_γγ shape.
        // Adjust based on your calibration strategy and data quality.
        //
        // Advanced usage: Adjust only if you have specific physics reasons
        //   - Higher W_MMISS: If you trust M_miss calibration more
        //   - Higher W_MPI0: If M_γγ has better statistics or precision
        //
        const double W_MPI0 = 1.0;   // Weight for invariant mass chi2 (M_γγ)
        const double W_MMISS = 1.5;  // Weight for missing mass chi2 (M_miss)
        const double W_MPGG2 = 0.0;  // Weight for (p_target + γγ)^2 chi2
        const double W_MPGG2_ENERGY = 1.0;  // Energy scaling factor for mpgg2 calculation: E -> w_mpgg2_energy * E
        
        // ========================================================================
        // HISTOGRAM SELECTION FOR ENERGY SMEARING
        // ========================================================================
        constexpr int HIST_BOTH = 0;
        constexpr int HIST_MPI0_ONLY = 1;
        constexpr int HIST_MMISS_ONLY = 2;
        constexpr int HIST_MPGG2_ONLY = 3;

        // Choose which histogram(s) to use for energy smearing optimization:
        //
        //   0: BOTH histograms (M_γγ AND M_miss) [RECOMMENDED]
        //      - M_miss constrains total energy: E₁+E₂
        //      - M_γγ constrains energy product: E₁×E₂ and opening angle θ
        //      - Together they fully constrain the kinematics
        //      - Uses combined chi2 with weights W_MPI0 and W_MMISS
        //
        //   1: INVARIANT MASS ONLY (M_γγ)
        //      - Good for angular and energy product calibration
        //      - May not match M_miss peak perfectly
        //
        //   2: MISSING MASS ONLY (M_miss) [NOT RECOMMENDED]
        //      - Only constrains total energy E₁+E₂
        //      - Leaves energy asymmetry distribution unconstrained
        //      - Will give POOR M_γγ agreement (see MISSING_MASS_FITTING_ISSUE.txt)
        //
        //   3: TARGET+2γ MASS^2 ONLY ((p_target + γγ)^2)
        //      - Uses: m_p^2 + 4E1E2 sin^2(θ/2) + 2m_p(E1+E2)
        //      - Strongly constrains total neutral energy scale
        //
        // NOTE: section sweeps fit energy response + energy resolution first,
        //       then fit sigma_pos with energy response fixed.
        //
        const int ENERGY_SMEARING_HISTOGRAM = HIST_BOTH;  // RECOMMENDED: HIST_BOTH
        
        // ========================================================================
        // FIT SCAN RANGES AND RESOLUTION
        // ========================================================================
        
        // Energy scale (mu) scan range and resolution
        // μ = 1.0 means no scaling, μ > 1.0 increases energies, μ < 1.0 decreases
        // Typical detector calibration: μ = 0.95-1.05 (2-5% corrections expected)
        // IMPORTANT: Set MU_MIN < 1.0 to allow downward energy scale corrections.
        //   μ < 1.0 is physically meaningful when simulation overestimates photon
        //   energies (common if the SIMC calorimeter response model is slightly off).
        //   Restricting MU_MIN = 1.0 traps the optimizer at the boundary — the fit
        //   then compensates with sigma instead, producing a wrong/poor result.
        const double MU_MIN = 0.95;
        const double MU_MAX = 1.12;
        const int MU_NSTEPS = 100;  // Fine sampling for precise energy scale

        // Energy-dependent photon energy-scale model:
        //   mu_eff(E) = a + b*E + c*ln(E)
        //
        // Enabled  (ENABLE_ENERGY_DEPENDENT_MU = true):
        //   All three coefficients (a, b, c) are per-section fit parameters.
        //   mu_eff(E) gives the scaled energy directly: E_sc = mu_eff(E).
        //
        // Disabled (ENABLE_ENERGY_DEPENDENT_MU = false):
        //   a = 0 and c = 0 are forced, only b is fitted.
        //   This reduces to E_sc = b * E in scalar-mu mode.
        //
        const bool ENABLE_ENERGY_DEPENDENT_MU = true;  // RECOMMENDED: true for best calibration fidelity, false for simpler model and faster computation
        const double MU_ENERGY_MIN_GEV = 0.2;  // energy floor for ln(E) evaluation

        // Initial / seed values for (a, b, c) when ENABLE_ENERGY_DEPENDENT_MU = true.
        // In disabled mode these are irrelevant; only b=1 is used as the mu seed.
        const double MU_ENERGY_A_INIT = 0.0;   // constant offset (GeV)
        const double MU_ENERGY_B_INIT = 1.0;   // linear coefficient (dimensionless)
        const double MU_ENERGY_C_INIT = 0.0;   // logarithmic coefficient (GeV)

        // Fit ranges for a and c (enabled mode only).
        // b always uses MU_MIN / MU_MAX defined above.
        const double MU_A_MIN = -0.20;  const double MU_A_MAX = 0.20;
        const double MU_C_MIN = -0.20;  const double MU_C_MAX = 0.20;
        const double MIGRAD_STEP_MU_A  = 0.0001;
        const double MIGRAD_STEP_MU_C  = 0.0001;

        // Guard against artificial a=c=0 results.
        // Coarse/fine grid scans seed only b and sigma; these deterministic
        // nonzero (a,c) starts test whether MIGRAD is trapped by the zero seed.
        const bool ENABLE_MU_AC_MULTISTART = true;
        const double MU_A_MULTISTART_SPAN = 0.05;  // GeV
        const double MU_C_MULTISTART_SPAN = 0.05;  // GeV

        
        // Resolution parameter (sigma) scan range
        // Units depend on model (see USE_SIMPLE_STOCHASTIC_MODEL above):
        //   Simple model: σ in GeV^(1/2), typically 0.02-0.05 for NPS PbWO4
        //   3-term model: σ dimensionless, ~1.0 if A,B,C are correct
        const double SIGMA_MIN = USE_SIMPLE_STOCHASTIC_MODEL ? 0.01 : 0.5;
        const double SIGMA_MAX = USE_SIMPLE_STOCHASTIC_MODEL ? 0.08 : 2.0;
        const int SIGMA_NSTEPS = 80;  // Fine sampling for resolution

        // ========================================================================
        // VISUALIZATION PERFORMANCE CONTROLS
        // ========================================================================
        // These settings affect only chi2 visualization scans stored for plots.
        // They do NOT change the optimization itself.
        // Lower values give much faster chi2_scans.pdf generation.
        const int VIS_MU_SCAN_POINTS = 28;       // 2D scan points in mu direction
        const int VIS_SIGMA_SCAN_POINTS = 28;    // 2D scan points in sigma direction
        const int VIS_SLICE_POINTS = 24;         // 1D scan points for mu/sigma slices
        const int VIS_SIGMA_POS_POINTS = 24;     // 1D sigma_pos points (if enabled)
        const double RESPONSE_CURVE_E_MIN_GEV = 0.5;
        const double RESPONSE_CURVE_E_MAX_GEV = 6.0;
        
        // ========================================================================
        // POSITION SMEARING CALIBRATION
        // ========================================================================
        // Enable/disable position smearing optimization
        // 
        // PHYSICS: Position smearing affects M_γγ through the opening angle:
        //          M_γγ² = 2E₁E₂(1 - cos θ)  where θ is the angle between photons
        //          
        //          Hit position uncertainty (x,y) → direction uncertainty → θ smearing
        //          This affects the WIDTH and TAIL of the M_γγ distribution!
        //
        // IMPORTANCE: Even if your position resolution is good (~1 mm), enabling this
        //             provides an important cross-check and accounts for residual
        //             angular resolution effects not captured by energy smearing alone.
        //
        // RECOMMENDATION: ENABLE (true) for complete detector calibration
        //                 The fit will find sigma_pos ≈ 0 if position effects are negligible
        //
        const bool ENABLE_POSITION_SMEARING = true;    // RECOMMENDED: true
        const double SIGMA_POS_MIN = 0.0;   // cm (start from zero)
        const double SIGMA_POS_MAX = 0.8;   // cm (NPS: ~0.5-1.5 cm typical, but allow range)
        const int SIGMA_POS_NSTEPS = 100;    // Sufficient sampling for clear minimum

        // Optional energy-dependent position smearing:
        //   sigma_pos_eff(E) = sigma_pos_ref * sqrt(E0 / E)
        // where sigma_pos_ref is the fitted sigma_pos value at reference energy E0.
        const bool ENABLE_ENERGY_DEPENDENT_SIGMA_POS = false;
        const double SIGMA_POS_ENERGY_E0_GEV = 2.0;
        const double SIGMA_POS_ENERGY_MIN_GEV = 0.2;

        // ========================================================================
        // OPTIONAL ELECTRON MOMENTUM SCALING (Missing-mass only)
        // ========================================================================
        // Applies a global scale to scattered-electron momentum in M_miss calculation:
        //   p_e' = p_e_scale * p_e
        // and recomputes electron energy from momentum magnitude and m_e.
        //
        // This parameter is global over the full calorimeter face (never section-wise).
        // Set ENABLE_ELECTRON_MOMENTUM_SCALING = false to keep electron momentum unchanged.
        const bool ENABLE_ELECTRON_MOMENTUM_SCALING = false;  // RECOMMENDED: false unless you have specific reasons to suspect electron momentum scale issues
        const bool ENABLE_FINAL_GLOBAL_PE_SCALE_STAGE4 = false;  // If false, keep per-section p_e_scale from coupled sweeps and skip final global p_e_scale fit
        const double GLOBAL_PE_SCALE_DEFAULT = 1.0;              // initial seed for global p_e_scale fit

        // Initial seeds used before global two-step fit:
        //   Step 1: fit (mu, sigma) with selected combined objective
        //   Step 2: final global p_e_scale fit with M_miss only, keeping (mu, sigma, sigma_pos) fixed
        const double GLOBAL_PE_FIT_MU = 1.0;
        const double GLOBAL_PE_FIT_SIGMA = 0.02;
        const double GLOBAL_PE_FIT_SIGMA_POS = 0.0;

        const double PE_SCALE_MIN = 0.98;
        const double PE_SCALE_MAX = 1.02;
        const int PE_SCALE_NSTEPS = 100;

        
        // Global multistart controls. Workflow:
        //   ROOT/GSL Sobol bounded seed scan -> keep best N seeds
        //   -> Minuit2 MIGRAD each retained seed -> HESSE correlation diagnostics
        //   -> limited profile scans for suspicious strongly correlated parameters.
        const int GLOBAL_MULTISTART_SEEDS = 256;
        const int GLOBAL_MULTISTART_KEEP_BEST = 8;
        const bool ENABLE_HESSE_DIAGNOSTICS = true;
        const bool ENABLE_PROFILE_SCANS = true;
        const bool ENABLE_PARALLEL_SECTION_FITS = true;
        const int PROFILE_SCAN_POINTS = 5;
        const int PROFILE_MAX_PARAMETERS = 1;
        const double PROFILE_CORR_THRESHOLD = 0.90;
        const double PROFILE_RANGE_FRACTION = 0.10;

        // Diagnostic grid densities used by visualization scans
        const int COARSE_GRID_DIVISOR = 3;         // Use Nsteps/3 for coarse search
        const int MAX_REFINEMENT_ITERATIONS = 1000;  // Maximum fine-grid refinement cycles

        // Optional Minuit2/MIGRAD local refinement after coarse/grid seeding.
        // Keeps current robust scan workflow while using a derivative-based
        // optimizer to tighten the final minimum.
        const bool USE_MIGRAD_REFINEMENT = true;
        const int MIGRAD_MAX_FUNCTION_CALLS = 40000;
        const int MIGRAD_MAX_ITERATIONS = 40000;
        const double MIGRAD_TOLERANCE = 1e-4;
        const int MIGRAD_STRATEGY = 1;
        const double MIGRAD_STEP_MU = 0.0001;
        const double MIGRAD_STEP_SIGMA = 0.0001;
        const double MIGRAD_STEP_SIGMA_POS = 0.001;
        const double MIGRAD_STEP_PE_SCALE = 5e-4;
        
        // ========================================================================
        // HISTOGRAM BINNING FOR CHI-SQUARED CALCULATION
        // ========================================================================
        
        // Invariant mass (M_γγ) - focused on π⁰ peak region (135 MeV/c²)
        const double MGGAMMA_MIN = 0.11;    // GeV/c²
        const double MGGAMMA_MAX = 0.15;    // GeV/c²
        const int MGGAMMA_NBINS = 50;      // ~0.4 MeV/bin for fine shape comparison
        
        // Missing mass (M_miss) - foc  used on proton mass region (938 MeV/c²)
        const double MMISS_MIN = 0.6;       // GeV/c²
        const double MMISS_MAX = 1.3;       // GeV/c²
        const int MMISS_NBINS = 80;        // ~4 MeV/bin for detailed peak comparison

        // (p_target + γγ)^2 in GeV²
        const double MPGG2_MIN = 4.0;
        const double MPGG2_MAX = 14.0;
        const int MPGG2_NBINS = 50;
        
        // ========================================================================
        // STATISTICS AND QUALITY CONTROL
        // ========================================================================

        // Chi-squared statistic for histogram comparison:
        //   false: Pearson chi2:  sum (s-d)^2 / (sigma_s^2 + sigma_d^2)
        //          Standard weighted-histogram comparison; uses both sim and data errors.
        //   true:  Baker-Cousins (log-likelihood ratio):  sum 2*(s - d + d*ln(d/s))
        //          Preferred in HEP for template fitting — unbiased for low-count bins
        //          and follows a chi-squared distribution more closely.  Does not use
        //          sim error bars, so it is immune to stochastic smearing noise inflation.
        const bool USE_BAKER_COUSINS_CHI2 = true;
        const double BAKER_COUSINS_EMPTY_SIM_FLOOR_ABS = 1e-9;
        const double BAKER_COUSINS_EMPTY_SIM_FLOOR_FRAC = 1e-12;

        // Optional exclusivity gating in event weights:
        //   false (default): use all events
        //   true: multiply data/sim weights by their configured exclusivity branch
        const bool APPLY_IS_EXCLUSIVE_SELECTION = true;
        const char* DATA_EXCLUSIVITY_BRANCH = "is_exclusive_ellipse_combined";
        const char* SIM_EXCLUSIVITY_BRANCH  = "is_exclusive_ellipse";

        // SIMC de-modeling: remove the generator cross-section model from event
        // weights so the smearing fit uses the accepted phase-space basis.
        const char* SIM_MODEL_XSEC_BRANCH = "sigcm";
        const double SIM_MODEL_XSEC_MIN_ABS = 1e-20;
        
        // Minimum events required per section (both data AND simulation)
        // Sections with insufficient statistics will be skipped with a warning
        const int MIN_EVENTS_PER_SECTION = 100;
        
        // Fit quality threshold: maximum acceptable chi2/ndf
        // Values > 2.0 indicate poor agreement, suggesting:
        //   - Insufficient smearing model complexity
        //   - Systematic differences between data and simulation
        //   - Statistical fluctuations (check if events > MIN_EVENTS_PER_SECTION)
        const double MAX_CHI2_PER_NDF = 2.0;  // Warning threshold
        const bool SKIP_BAD_FITS = false;     // If true, exclude bad sections from output

        // Runtime controls for optimizer only. Final chi2 recompute, final section
        // plots, section maps, and all-section summaries still use the full event
        // buffers and the command-line Nsmear value.
        const int OPTIMIZATION_NSMEAR = 80;
        const int OPT_MAX_SIM_EVENTS_PER_SECTION = 25000;
        const int OPT_MAX_SIM_EVENTS_GLOBAL_PREFIT = 60000;
        const int OPT_SUBSET_MGG_BINS = 32;
        const double OPT_SUBSET_MGG_MIN = 0.05;
        const double OPT_SUBSET_MGG_MAX = 0.25;

        // Section fit orchestration.
        // Iterative coupled sweep model:
        //   prefit: fit one global all-calorimeter response with no section split.
        //   sweep 1: each section starts from global prefit parameters, and
        //            out-of-section photons use the same global response.
        //            If global prefit is disabled/low-stat, falls back to nominal
        //            (a=0,b=1,c=0,sigma=0,sigma_pos=0).
        //   sweep N: each section fits in-section parameters while out-of-section
        //            photons use completed section parameters from sweep N-1.
        //
        // Each section fit still uses Sobol -> keep best N -> MIGRAD -> HESSE/profile.
        const bool ENABLE_GLOBAL_PREFIT_SEED = true;
        const int ITERATIVE_SECTION_SWEEPS = 10;
        const bool ENABLE_COUPLED_SWEEP_CONVERGENCE_STOP = true;
        const double COUPLED_CONV_MU = 5e-4;
        const double COUPLED_CONV_SIGMA = 5e-4;
        const double COUPLED_CONV_SIGMA_POS = 5e-4;
        const double COUPLED_ACCEPT_REL_TOL = 1e-8;
        const double COUPLED_ACCEPT_ABS_TOL = 1e-6;
        const double COUPLED_REPEAT_NORM_TOL = 1e-7;
        const bool ENABLE_COUPLED_REJECTED_REPEAT_STOP = true;
        const int COUPLED_REJECTED_REPEAT_PATIENCE = 1;
        const int COUPLED_REJECTED_CYCLE_PATIENCE = 2;
        
        // ========================================================================
        // OUTPUT FILE SETTINGS
        // ========================================================================
        
        // Default output directory (prepended to relative paths)
        // const string OUTPUT_DIR = "../output/plots/x60_4b/";
        const string OUTPUT_DIR = "../output/plots/x60_4b/production_wfpi0/";
        
        // Output file naming
        const string CSV_FILENAME = "section_map.csv";
        const string OPTIMIZER_SUMMARY_CSV_FILENAME = "smearing_optimizer_summary.csv";
        const string OPTIMIZER_SEEDS_CSV_FILENAME = "smearing_optimizer_seeds.csv";
        const string OPTIMIZER_PROFILE_CSV_FILENAME = "smearing_optimizer_profiles.csv";
        const string CLOSURE_SUMMARY_CSV_FILENAME = "smearing_closure_summary.csv";
        const string SWEEP_HISTORY_CSV_FILENAME = "smearing_sweep_history.csv";
        const string OBJECTIVE_BREAKDOWN_CSV_FILENAME = "smearing_objective_breakdown.csv";
        const string CACHE_FINGERPRINT_FILENAME = "smearing_config_fingerprint.txt";
        const string CHI2_PDF_FILENAME = "chi2_scans.pdf";
        const string INTERPOLATED_SUFFIX = "_interpolated";
        
        // ========================================================================
        // CONFIGURATION SUMMARY
        // ========================================================================
        // This section helps verify your configuration intent.
        // 
        // FOR THE CURRENT COUPLED-SWEEP PROFILE, verify these settings:
        //   - ITERATIVE_SECTION_SWEEPS set intentionally
        //   - Sobol multistart + MIGRAD/HESSE/profile
        //   - ENERGY_SMEARING_HISTOGRAM = HIST_BOTH
        //   - ENABLE_POSITION_SMEARING = false
        //   - ENABLE_ELECTRON_MOMENTUM_SCALING = false
        //   - W_MPI0, W_MMISS, and W_MPGG2 set intentionally for your strategy
        //
        // CURRENT DEFAULTS IN THIS FILE:
        //   - iterative coupled section sweeps
        //   - ENERGY_SMEARING_HISTOGRAM = HIST_BOTH
        //   - ENABLE_POSITION_SMEARING = false
        //   - W_MPI0 = 2.0, W_MMISS = 1.0, W_MPGG2 = 1.0
        //
        // CAUTION:
        //   ✗ ENERGY_SMEARING_HISTOGRAM = HIST_MMISS_ONLY (M_miss only)
        //   ✗ very wide a/c bounds can create unphysical low-energy response
        //   See MISSING_MASS_FITTING_ISSUE.txt for detailed physics explanation
        // ========================================================================

        inline const char* histogram_mode_label() {
            if (ENERGY_SMEARING_HISTOGRAM == HIST_MPI0_ONLY) return "M_gg only";
            if (ENERGY_SMEARING_HISTOGRAM == HIST_MMISS_ONLY) return "M_miss only";
            if (ENERGY_SMEARING_HISTOGRAM == HIST_MPGG2_ONLY) return "(p_target + #gamma#gamma)^2 only";
            return "M_gg + M_miss + (p_target + #gamma#gamma)^2";
        }

        inline bool fit_objective_uses_mmiss() {
            return (ENERGY_SMEARING_HISTOGRAM == HIST_MMISS_ONLY ||
                    ENERGY_SMEARING_HISTOGRAM == HIST_BOTH);
        }

        inline string sweep_acceptance_strategy() {
            const char *env = std::getenv("NPS_SMEARING_SWEEP_ACCEPTANCE");
            string value = (env && env[0]) ? string(env) : "rollback";
            return "jacobi_global_accept_rollback";
        }

        inline void print_configuration_summary() {
            cout << "\n==== Active Configuration ====\n";
            cout << "Optimization mode: Sobol multistart + MIGRAD/HESSE/profile\n";
            cout << "Section fit staging: energy response + energy resolution, then sigma_pos\n";
            cout << "Global seeds: " << GLOBAL_MULTISTART_SEEDS
                << "  keep-best: " << GLOBAL_MULTISTART_KEEP_BEST
                << "  HESSE=" << (ENABLE_HESSE_DIAGNOSTICS ? "on" : "off")
                << "  profile=" << (ENABLE_PROFILE_SCANS ? "on" : "off") << "\n";
            cout << "Parallel section fits: "
                << (ENABLE_PARALLEL_SECTION_FITS ? "enabled" : "disabled") << "\n";
            cout << "Section orchestration: iterative coupled sweeps\n"
                 << "  sweeps=" << ITERATIVE_SECTION_SWEEPS
                 << "  convergence_stop=" << (ENABLE_COUPLED_SWEEP_CONVERGENCE_STOP ? "on" : "off") << "\n"
                 << "  global prefit seed: " << (ENABLE_GLOBAL_PREFIT_SEED ? "enabled" : "disabled") << "\n"
                 << "  sweep 1 external photons: global prefit response (nominal fallback)\n"
                 << "  sweep acceptance strategy: " << sweep_acceptance_strategy() << "\n"
                 << "  sweep N external photons: previous accepted sweep results\n";
            cout << "Fit observable: " << histogram_mode_label() << "\n";
            cout << "Optimizer acceleration: Nsmear <= " << OPTIMIZATION_NSMEAR
                 << ", max sim events/section=" << OPT_MAX_SIM_EVENTS_PER_SECTION
                 << ", max global prefit sim events=" << OPT_MAX_SIM_EVENTS_GLOBAL_PREFIT << "\n";
            cout << "Position smearing fit: "
                << (ENABLE_POSITION_SMEARING ? "enabled" : "disabled") << "\n";
            if (ENABLE_POSITION_SMEARING) {
                cout << "Position smearing energy dependence: "
                    << (ENABLE_ENERGY_DEPENDENT_SIGMA_POS ? "enabled" : "disabled") << "\n";
                if (ENABLE_ENERGY_DEPENDENT_SIGMA_POS) {
                    cout << "  sigma_pos(E)=sigma_pos_ref*sqrt(E0/E) with E0="
                        << SIGMA_POS_ENERGY_E0_GEV << " GeV\n";
                }
            }
            cout << "Electron momentum scaling: "
                << (ENABLE_ELECTRON_MOMENTUM_SCALING ? "enabled" : "disabled") << "\n";
            if (ENABLE_ELECTRON_MOMENTUM_SCALING && fit_objective_uses_mmiss()) {
                if (ENABLE_FINAL_GLOBAL_PE_SCALE_STAGE4) {
                    cout << "p_e_scale mode: per-section in coupled sweeps + final global refinement\n";
                } else {
                    cout << "p_e_scale mode: per-section in coupled sweeps only (final global refinement disabled)\n";
                }
            }
            cout << "Weights (W_MPI0, W_MMISS, W_MPGG2): " << W_MPI0 << ", " << W_MMISS << ", " << W_MPGG2 << "\n";
            cout << "SIMC cross-section de-modeling: enabled; sim event weight = full_weight/"
                 << SIM_MODEL_XSEC_BRANCH << "\n";
            cout << "Exclusive gating in weights: "
                << (APPLY_IS_EXCLUSIVE_SELECTION ? "enabled" : "disabled (all events)") << "\n";
            if (APPLY_IS_EXCLUSIVE_SELECTION) {
                cout << "Data exclusive branch: " << DATA_EXCLUSIVITY_BRANCH << "\n";
                cout << "Sim exclusive branch:  " << SIM_EXCLUSIVITY_BRANCH << "\n";
            }
            cout << "mpgg2 energy scaling (w_mpgg2_energy): " << W_MPGG2_ENERGY << "\n";
            cout << "Energy-smearing PDF: "
                << (ENERGY_SMEAR_SHAPE == SMEAR_SHAPE_LANDAU ? "Landau (sigma=FWHM)" : "Gaussian (sigma=std-dev)")
                << "\n";
            cout << "m_gg linear pair-energy correction: "
                << (ENABLE_MGG_LINEAR_ENERGY_CORRECTION ? "enabled" : "disabled") << "\n";
            if (ENABLE_MGG_LINEAR_ENERGY_CORRECTION) {
                cout << "  slope=" << MGG_LINEAR_SLOPE
                    << " GeV^-1, pivot=" << MGG_LINEAR_PIVOT_GEV
                    << " GeV, mode=" << (MGG_LINEAR_USE_INVERSE ? "inverse" : "direct")
                    << ", clamp=[" << MGG_LINEAR_FACTOR_MIN << ", " << MGG_LINEAR_FACTOR_MAX << "]\n";
            }
            if (ENERGY_SMEAR_SHAPE == SMEAR_SHAPE_LANDAU) {
                cout << "  Landau conversion: scale = FWHM / " << LANDAU_FWHM_TO_SCALE << "\n";
                cout << "  Landau truncation: "
                    << (ENABLE_TRUNCATED_LANDAU ? "enabled" : "disabled")
                    << ", max redraws=" << LANDAU_MAX_REDRAWS
                    << ", upper tail cap=" << LANDAU_MAX_FWHM_ABOVE_MPV << " x FWHM\n";
            }
            cout << "Energy-dependent mu(E): "
                << (ENABLE_ENERGY_DEPENDENT_MU ? "enabled" : "disabled") << "\n";
            if (ENABLE_ENERGY_DEPENDENT_MU) {
                cout << "  mu_eff(E) = a + b*E + c*ln(E)  (all three fitted per section)\n"
                     << "  initial: a=" << MU_ENERGY_A_INIT
                     << "  b=" << MU_ENERGY_B_INIT
                     << "  c=" << MU_ENERGY_C_INIT << "\n"
                     << "  energy floor: E_floor=" << MU_ENERGY_MIN_GEV << " GeV\n";
                cout << "  a/c MIGRAD multistart: "
                     << (ENABLE_MU_AC_MULTISTART ? "enabled" : "disabled")
                     << "  seed spans: da=" << MU_A_MULTISTART_SPAN
                     << " GeV, dc=" << MU_C_MULTISTART_SPAN << " GeV\n";
            } else {
                cout << "  mu_eff(E) = b*E  (a=0, c=0 fixed; b is the fitted scalar mu)\n";
            }
            cout << "Coupled sweep convergence thresholds: "
                 << "mu=" << COUPLED_CONV_MU
                 << " sigma=" << COUPLED_CONV_SIGMA
                 << " sigma_pos=" << COUPLED_CONV_SIGMA_POS << "\n";
            cout << "Rejected repeated-candidate stop: "
                 << (ENABLE_COUPLED_REJECTED_REPEAT_STOP ? "enabled" : "disabled")
                 << " repeat_patience=" << COUPLED_REJECTED_REPEAT_PATIENCE
                 << " cycle_patience=" << COUPLED_REJECTED_CYCLE_PATIENCE << "\n";

            if (ENABLE_ELECTRON_MOMENTUM_SCALING && !fit_objective_uses_mmiss()) {
                cout << "Note: p_e_scale is disabled when fit objective does not use M_miss\n";
            }
            if (!ENABLE_POSITION_SMEARING && SIGMA_POS_MAX > 0.0) {
                cout << "Note: sigma_pos scan range is configured but unused because ENABLE_POSITION_SMEARING=false\n";
            }
            cout << "==============================\n";
        }
        
    } // namespace Config

    // ============================================================================
    //                      END OF USER CONFIGURATION
    // ============================================================================

    static string currentTimestampTag() {
        std::time_t now = std::time(nullptr);
        std::tm tm_now;
    #if defined(_WIN32)
        localtime_s(&tm_now, &now);
    #else
        localtime_r(&now, &tm_now);
    #endif
        char buf[32];
        std::strftime(buf, sizeof(buf), "%Y%m%d_%H%M%S", &tm_now);
        return string(buf);
    }

    static string sanitizeFileTag(const string &tag) {
        string out;
        out.reserve(tag.size());
        for (char ch : tag) {
            unsigned char uch = static_cast<unsigned char>(ch);
            if (std::isalnum(uch) || ch == '_' || ch == '-' || ch == '.') {
                out.push_back(ch);
            } else {
                out.push_back('_');
            }
        }
        return out.empty() ? currentTimestampTag() : out;
    }

    static string getRunTag() {
        const char *env_tag = std::getenv("RUN_TAG");
        if (env_tag && env_tag[0]) return sanitizeFileTag(env_tag);
        return currentTimestampTag();
    }

    static string insertTagBeforeExtension(const string &path, const string &tag) {
        if (tag.empty()) return path;
        size_t slash = path.find_last_of('/');
        size_t dot = path.find_last_of('.');
        size_t basename_start = (slash == string::npos) ? 0 : slash + 1;
        string suffix = "_" + tag;
        if (dot == string::npos || dot < basename_start) {
            if (path.size() >= suffix.size() &&
                path.compare(path.size() - suffix.size(), suffix.size(), suffix) == 0) {
                return path;
            }
            return path + suffix;
        }
        string stem = path.substr(0, dot);
        if (stem.size() >= suffix.size() &&
            stem.compare(stem.size() - suffix.size(), suffix.size(), suffix) == 0) {
            return path;
        }
        return stem + suffix + path.substr(dot);
    }

    static string interpolatedPathForOutput(const string &out_file) {
        size_t dot = out_file.find_last_of('.');
        size_t slash = out_file.find_last_of('/');
        size_t basename_start = (slash == string::npos) ? 0 : slash + 1;
        if (dot != string::npos && dot >= basename_start) {
            return out_file.substr(0, dot) + Config::INTERPOLATED_SUFFIX + out_file.substr(dot);
        }
        return out_file + Config::INTERPOLATED_SUFFIX + ".root";
    }

    static void copyFileIfDifferent(const string &src, const string &dst, const char *label) {
        if (src.empty() || dst.empty() || src == dst) return;
        int rc = gSystem->CopyFile(src.c_str(), dst.c_str(), true);
        if (rc == 0) {
            cout << "Saved timestamped " << label << ": " << dst << "\n";
        } else {
            cerr << "WARNING: Could not copy " << label << " from " << src
                 << " to " << dst << " (rc=" << rc << ")\n";
        }
    }

    static void writeNamedString(TDirectory *dir, const string &name, const string &value) {
        if (!dir) return;
        TDirectory *save_dir = gDirectory;
        dir->cd();
        TNamed named(name.c_str(), value.c_str());
        named.Write(name.c_str(), TObject::kOverwrite);
        if (save_dir) save_dir->cd();
    }

    static void writeCanvasToDir(TDirectory *dir, TCanvas *canvas,
                                 const string &name, const string &title = "") {
        if (!dir || !canvas) return;
        TDirectory *save_dir = gDirectory;
        dir->cd();
        canvas->SetName(name.c_str());
        if (!title.empty()) canvas->SetTitle(title.c_str());
        canvas->Write(name.c_str(), TObject::kOverwrite);
        if (save_dir) save_dir->cd();
    }

    static void writeHistToDir(TDirectory *dir, TH1 *hist, const string &name = "") {
        if (!dir || !hist) return;
        TDirectory *save_dir = gDirectory;
        dir->cd();
        string write_name = name.empty() ? hist->GetName() : name;
        hist->Write(write_name.c_str(), TObject::kOverwrite);
        if (save_dir) save_dir->cd();
    }

    static void writeSmearingManifest(const string &filename,
                                      const string &run_tag,
                                      const string &created_at,
                                      const string &data_file,
                                      const string &sim_file,
                                      const string &out_file,
                                      const string &timestamped_out_file,
                                      const string &interp_file,
                                      const string &timestamped_interp_file,
                                      const string &pdf_file,
                                      const string &canonical_pdf_file,
                                      const string &progress_pdf_dir,
                                      const string &current_cache_fingerprint) {
        ofstream meta(filename.c_str());
        if (!meta.is_open()) {
            cerr << "WARNING: Could not write smearing metadata manifest: " << filename << "\n";
            return;
        }
        meta << "{\n";
        meta << "  \"run_tag\": \"" << run_tag << "\",\n";
        meta << "  \"created_at_local\": \"" << created_at << "\",\n";
        meta << "  \"data_file\": \"" << data_file << "\",\n";
        meta << "  \"sim_file\": \"" << sim_file << "\",\n";
        meta << "  \"out_file\": \"" << out_file << "\",\n";
        meta << "  \"timestamped_out_file\": \"" << timestamped_out_file << "\",\n";
        meta << "  \"interpolated_file\": \"" << interp_file << "\",\n";
        meta << "  \"timestamped_interpolated_file\": \"" << timestamped_interp_file << "\",\n";
        meta << "  \"chi2_pdf\": \"" << pdf_file << "\",\n";
        meta << "  \"canonical_chi2_pdf\": \"" << canonical_pdf_file << "\",\n";
        meta << "  \"chi2_progress_dir\": \"" << progress_pdf_dir << "\",\n";
        meta << "  \"cache_fingerprint\": ";
        meta << "\"";
        for (char ch : current_cache_fingerprint) {
            if (ch == '\\') meta << "\\\\";
            else if (ch == '"') meta << "\\\"";
            else if (ch == '\n') meta << "\\n";
            else meta << ch;
        }
        meta << "\"\n";
        meta << "}\n";
        cout << "Smearing metadata manifest saved to " << filename << "\n";
    }

    inline double computeEnergyResolution(double E_scaled, double sigma,
                                          double res_A, double res_B, double res_C);

    struct ModelParameter {
        string name;
        double value;
        double min_value;
        double max_value;
        string unit;
        string description;

        bool isValid() const {
            return std::isfinite(value) && value >= min_value && value <= max_value;
        }
    };

    struct SmearingModel1D {
        enum class Type {
            Constant,
            Linear,
            Polynomial,
            APlusBEPlusCLnE
        };

        Type type = Type::APlusBEPlusCLnE;
        string name = "a_plus_bE_plus_clnE";
        vector<ModelParameter> parameters;
        double x_floor = Config::MU_ENERGY_MIN_GEV;

        double evaluate(double x) const {
            const double x_safe = std::max(x, x_floor);
            if (type == Type::Constant) {
                return parameters.empty() ? 0.0 : parameters[0].value;
            }
            if (type == Type::Linear) {
                const double a = parameters.size() > 0 ? parameters[0].value : 0.0;
                const double b = parameters.size() > 1 ? parameters[1].value : 1.0;
                return a + b * x_safe;
            }
            if (type == Type::Polynomial) {
                double y = 0.0;
                double xp = 1.0;
                for (const auto &p : parameters) {
                    y += p.value * xp;
                    xp *= x_safe;
                }
                return y;
            }
            const double a = parameters.size() > 0 ? parameters[0].value : 0.0;
            const double b = parameters.size() > 1 ? parameters[1].value : 1.0;
            const double c = parameters.size() > 2 ? parameters[2].value : 0.0;
            return a + b * x_safe + c * std::log(x_safe);
        }

        bool isValid() const {
            if (!(x_floor > 0.0) || !std::isfinite(x_floor)) return false;
            for (const auto &p : parameters) {
                if (!p.isValid()) return false;
            }
            return true;
        }
    };

    struct PhotonSmearingParameters {
        double energy_mean_a = 0.0;
        double energy_mean_b = 1.0;
        double energy_mean_c = 0.0;
        double energy_sigma = 0.0;
        double position_sigma = 0.0;
        double res_A = Config::RESOLUTION_A_DEFAULT;
        double res_B = Config::RESOLUTION_B_DEFAULT;
        double res_C = Config::RESOLUTION_C_DEFAULT;
    };

    inline SmearingModel1D makeEnergyMeanModel(double a, double b, double c) {
        SmearingModel1D model;
        model.type = SmearingModel1D::Type::APlusBEPlusCLnE;
        model.name = "a_plus_bE_plus_clnE";
        model.x_floor = Config::MU_ENERGY_MIN_GEV;
        model.parameters = {
            {"a", a, Config::MU_A_MIN, Config::MU_A_MAX, "GeV", "constant reconstructed-energy offset"},
            {"b", b, Config::MU_MIN, Config::MU_MAX, "dimensionless", "linear energy-response coefficient"},
            {"c", c, Config::MU_C_MIN, Config::MU_C_MAX, "GeV", "logarithmic energy-response coefficient"}
        };
        return model;
    }

    struct PhotonSmearingModel {
        PhotonSmearingParameters params;
        SmearingModel1D energy_mean_model;

        explicit PhotonSmearingModel(const PhotonSmearingParameters &p)
            : params(p), energy_mean_model(makeEnergyMeanModel(p.energy_mean_a,
                                                               p.energy_mean_b,
                                                               p.energy_mean_c)) {}

        double meanEnergy(double E) const {
            double E_mean = energy_mean_model.evaluate(E);
            if (!std::isfinite(E_mean) || E_mean <= 0.0) return Config::NONPOSITIVE_CLAMP;
            return E_mean;
        }

        double sigmaEnergy(double E_mean) const {
            if (!(params.energy_sigma > 0.0) || !std::isfinite(params.energy_sigma)) return 0.0;
            return computeEnergyResolution(E_mean, params.energy_sigma,
                                           params.res_A, params.res_B, params.res_C);
        }

        double sigmaPosition(double E_mean) const {
            if (!(params.position_sigma > 0.0) || !std::isfinite(params.position_sigma)) return 0.0;
            if (!Config::ENABLE_ENERGY_DEPENDENT_SIGMA_POS) return params.position_sigma;
            const double E_for_pos = std::max(E_mean, Config::SIGMA_POS_ENERGY_MIN_GEV);
            return params.position_sigma * sqrt(Config::SIGMA_POS_ENERGY_E0_GEV / E_for_pos);
        }
    };

    // ============================================================================
	    // HMS ELECTRON CUTS (for simulation - consistent with simc_pi0_analysis.C)
    // NOTE: These cuts are applied in simc_pi0_analysis.C, not in this script.
    //       This function is kept for reference/documentation purposes only.
    // ============================================================================
    // Same cuts as in simc_pi0_analysis.C - simulation doesn't have npesum/etotnorm
    inline bool hms_electron_cuts_simulation(double h_delta, double h_gtr_th,
                                            double h_gtr_ph, double h_react_z) noexcept
    {
        // Simulation doesn't have npesum and etotnorm, so we skip those cuts
        if (h_react_z < -8.0 || h_react_z > 8.0) return false;
        if (h_delta < -15.0 || h_delta > 15.0) return false;
        if (h_gtr_th < -0.1 || h_gtr_th > 0.1) return false;
        if (h_gtr_ph < -0.04 || h_gtr_ph > 0.04) return false;
        return true;
    }


    // Helper structures
    struct Section {
        int ix, iy;               // grid indices
        double xlo, xhi, ylo, yhi; // bounds for this section
        string name() const {
            ostringstream os; os << "sec_" << ix << "_" << iy; return os.str();
        }
    };

    // Check whether a point (x,y) is inside the section (inclusive)
    bool inSection(const Section &s, double x, double y) {
        return (x >= s.xlo && x <= s.xhi && y >= s.ylo && y <= s.yhi);
    }

    // Event structure to hold cluster information
    struct ClusterPair {
        double e1, x1, y1;
        double e2, x2, y2;
        double weight;
        // Electron kinematics for missing mass calculation
        double Ee, px_e, py_e, pz_e;
        // Section-membership flags.
        // photon1_in_section = true  → photon 1 hit position is inside the current section.
        // photon2_in_section = false → photon 2 belongs to a different section.
        bool photon1_in_section = false;
        bool photon2_in_section = false;
        // External calibration for out-of-section photons.
        // Sweep 0/global prefit seeds these values; later sweeps refresh them
        // from the previous accepted/completed section map before optimization.
        // Defaults: mu_eff(E)=b*E with b=1 (no scale change), sigma=0 (no resolution smearing).
        // mu_eff(E) = mu_a_ext + mu_ext*E + mu_c_ext*ln(E)
        double mu_a1_ext = 0.0, mu1_ext = 1.0, mu_c1_ext = 0.0;
        double sigma1_ext = 0.0, sigma_pos1_ext = 0.0;
        double mu_a2_ext = 0.0, mu2_ext = 1.0, mu_c2_ext = 0.0;
        double sigma2_ext = 0.0, sigma_pos2_ext = 0.0;
    };

    // Lightweight data-event record for global all-sections diagnostics
    struct DataGlobalEvent {
        double x1, y1;
        double x2, y2;
        double weight;
    };

    // ============================================================================
    // HELPER FUNCTIONS FOR RESOLUTION CALCULATIONS
    // ============================================================================

    // Calculate energy-dependent resolution for a photon
    // Returns σ_E based on selected model (simple stochastic or 3-term)
    inline double computeEnergyResolution(double E_scaled, double sigma,
                                        double res_A, double res_B, double res_C) {
        if (!(E_scaled > 0.0) || !std::isfinite(E_scaled) ||
            !(sigma >= 0.0) || !std::isfinite(sigma)) {
            return 0.0;
        }
        if (Config::USE_SIMPLE_STOCHASTIC_MODEL) {
            // Simple stochastic: σ_E = σ × √E
            return sigma * sqrt(E_scaled);
        } else {
            // Full 3-term: σ_E = σ × E × √((A/√E)² + (B/E)² + C²)
            double A_sq = res_A * res_A;
            double B_sq = res_B * res_B;
            double C_sq = res_C * res_C;
            double sigma_rel_sq = A_sq / E_scaled + B_sq / (E_scaled * E_scaled) + C_sq;
            if (!(sigma_rel_sq >= 0.0) || !std::isfinite(sigma_rel_sq)) return 0.0;
            return sigma * E_scaled * sqrt(sigma_rel_sq);
        }
    }

    // ============================================================================
    // PER-PHOTON SMEARING HELPER
    // ============================================================================
    // Applies energy scale + optional Gaussian resolution smearing to one photon.
    // Uses the provided calibration parameters (mu, sigma, sigma_pos):
    //   mu        — energy scale factor
    //   sigma     — stochastic resolution coefficient (0 → scale only, no noise)
    //   sigma_pos — position resolution in cm (0 → no position shift)
    // Called every Nsmear iteration; parameters are pre-resolved (in-section vs ext).
    // smearPhoton: apply energy-scale model + optional resolution smearing to one photon.
    //
    // mu_a, mu_b, mu_c — coefficients of  mu_eff(E) = mu_a + mu_b*E + mu_c*ln(E).
    //   Disabled mode: mu_a=0, mu_c=0 -> E_sc = mu_b * E.
    //   Enabled  mode: all three non-zero → E_sc = mu_a + mu_b*E + mu_c*ln(E_safe).
    //
    // An energy floor (Config::MU_ENERGY_MIN_GEV) is applied before ln(E) to prevent
    // log singularities at very low energies.
    inline void smearPhoton(double E, double x, double y,
                            double mu_a, double mu_b, double mu_c,
                            double sigma, double sigma_pos,
                            double res_A, double res_B, double res_C,
                            TRandom3 &rng,
                            double &E_out, double &x_out, double &y_out) {
	        PhotonSmearingParameters photon_params;
	        photon_params.energy_mean_a = mu_a;
	        photon_params.energy_mean_b = mu_b;
	        photon_params.energy_mean_c = mu_c;
	        photon_params.energy_sigma = sigma;
	        photon_params.position_sigma = sigma_pos;
	        photon_params.res_A = res_A;
	        photon_params.res_B = res_B;
	        photon_params.res_C = res_C;
	        PhotonSmearingModel photon_model(photon_params);

	        // Single active energy-response convention:
	        // E_mean(E) = a + b*E_safe + c*ln(E_safe), in GeV.
	        double E_sc = photon_model.meanEnergy(E);

        // --- Deterministic-draw block ------------------------------------------------
        // Always consume a FIXED number of RNG values per photon regardless of the
        // current parameter values (sigma, sigma_pos).  This guarantees that the RNG
        // state after processing each event is independent of the parameters being
        // scanned, producing a smooth, deterministic chi2 surface when the RNG seed
        // is reset before each evaluation.
        //   Gaussian mode: 3 draws per photon  (1 energy + 2 position)
        //   Landau  mode:  LANDAU_MAX_REDRAWS + 2 draws per photon
        // --------------------------------------------------------------------------

        if (Config::ENERGY_SMEAR_SHAPE == Config::SMEAR_SHAPE_LANDAU) {
            // Always consume exactly LANDAU_MAX_REDRAWS Landau draws for energy.
            // Only the first accepted draw is used; the rest keep the RNG state
            // deterministic across parameter evaluations.
	            double fwhmE = photon_model.sigmaEnergy(E_sc);
            double landau_scale = fwhmE / max(Config::LANDAU_FWHM_TO_SCALE, 1e-9);
            double e_max_allowed = E_sc + Config::LANDAU_MAX_FWHM_ABOVE_MPV * fwhmE;
            e_max_allowed = max(e_max_allowed, Config::NONPOSITIVE_CLAMP);

            bool accepted = false;
            double draw = E_sc;
            int n_try = max(1, Config::LANDAU_MAX_REDRAWS);
            for (int itry = 0; itry < n_try; ++itry) {
                double this_draw = rng.Landau(E_sc, landau_scale);
                if (!accepted && sigma > 0.0) {
                    if (this_draw > Config::NONPOSITIVE_CLAMP &&
                        (!Config::ENABLE_TRUNCATED_LANDAU || this_draw <= e_max_allowed)) {
                        accepted = true;
                        draw = this_draw;
                    }
                }
                // Continue drawing even after acceptance to consume all RNG values
            }
            E_out = (sigma > 0.0 && accepted) ? draw : E_sc;
            if (E_out <= 0.0) E_out = Config::NONPOSITIVE_CLAMP;
        } else {
            // Gaussian mode — always draw exactly 1 unit-Gaussian pull, even if
            // sigma == 0 (in which case the pull is unused).
            double pull_e = rng.Gaus(0.0, 1.0);
            if (sigma > 0.0) {
	                double sigE = photon_model.sigmaEnergy(E_sc);
	                E_out = E_sc + sigE * pull_e;
                if (E_out <= 0.0) E_out = Config::NONPOSITIVE_CLAMP;
            } else {
                E_out = E_sc;
            }
        }

        // Position smearing — always draw 2 unit-Gaussian pulls to keep the RNG
        // state independent of sigma_pos.
        double pull_x = rng.Gaus(0.0, 1.0);
        double pull_y = rng.Gaus(0.0, 1.0);
        if (sigma_pos > 0.0) {
	            double sigma_pos_eff = photon_model.sigmaPosition(E_sc);
	            x_out = x + sigma_pos_eff * pull_x;
            y_out = y + sigma_pos_eff * pull_y;
        } else {
            x_out = x;
            y_out = y;
        }
    }

    // Calculate photon direction vector and momentum components
    // Given energy and position, returns (px, py, pz) components
    struct PhotonMomentum {
        double px, py, pz;
    };

    struct ElectronKinematics {
        double Ee;
        double px;
        double py;
        double pz;
    };

    inline PhotonMomentum computePhotonMomentum(double E, double x, double y, double z) {
        double r_sq = x*x + y*y + z*z;
        double r_inv = 1.0 / sqrt(r_sq);
        return {E * x * r_inv, E * y * r_inv, E * z * r_inv};
    }

    inline ElectronKinematics scaleElectronKinematicsFromMomentum(double px_e, double py_e, double pz_e,
                                                                double p_e_scale) {
        constexpr double kElectronMassGeV = 0.00051099895;
        double px_scaled = p_e_scale * px_e;
        double py_scaled = p_e_scale * py_e;
        double pz_scaled = p_e_scale * pz_e;
        double p2_scaled = px_scaled*px_scaled + py_scaled*py_scaled + pz_scaled*pz_scaled;
        double Ee_scaled = sqrt(p2_scaled + kElectronMassGeV * kElectronMassGeV);
        return {Ee_scaled, px_scaled, py_scaled, pz_scaled};
    }

    // Calculate invariant mass from two photons (direct calculation, avoids TLorentzVector overhead)
    inline double computeInvariantMass(double E1, double px1, double py1, double pz1,
                                    double E2, double px2, double py2, double pz2) {
        double E_tot = E1 + E2;
        double px_tot = px1 + px2;
        double py_tot = py1 + py2;
        double pz_tot = pz1 + pz2;
        double m_sq = E_tot*E_tot - px_tot*px_tot - py_tot*py_tot - pz_tot*pz_tot;
        return (m_sq > 0) ? sqrt(m_sq) : 0.0;
    }

    inline double computeTargetPlusDiphotonMass2(double E1, const PhotonMomentum &p1,
                                                double E2, const PhotonMomentum &p2) {
        constexpr double kProtonMassGeV = 0.9382720813;
        double denom = E1 * E2;
        if (denom <= 0.0) return kProtonMassGeV * kProtonMassGeV;

        double dot = p1.px * p2.px + p1.py * p2.py + p1.pz * p2.pz;
        double cos_theta = dot / denom;
        cos_theta = max(-1.0, min(1.0, cos_theta));

        double diphoton_mass2 = 2.0 * E1 * E2 * (1.0 - cos_theta);
        return kProtonMassGeV * kProtonMassGeV + diphoton_mass2 + 2.0 * kProtonMassGeV * (E1 + E2);
    }

    inline void apply_mgg_linear_pair_energy_correction(double &E1,
                                                        double &E2,
                                                        double x1,
                                                        double y1,
                                                        double x2,
                                                        double y2) {
        if (!Config::ENABLE_MGG_LINEAR_ENERGY_CORRECTION) return;
        if (E1 <= 0.0 || E2 <= 0.0) return;

        const double mgg = nps::invariant_mass_pi0(E1, E2, x1, x2, y1, y2, nps::kDefaultZ_NPS_cm);
        if (!std::isfinite(mgg)) return;

        double factor = 1.0 + Config::MGG_LINEAR_SLOPE * (mgg - Config::MGG_LINEAR_PIVOT_GEV);
        if (Config::MGG_LINEAR_USE_INVERSE) {
            if (std::abs(factor) < 1e-9) return;
            factor = 1.0 / factor;
        }

        factor = max(Config::MGG_LINEAR_FACTOR_MIN, min(Config::MGG_LINEAR_FACTOR_MAX, factor));
        E1 = max(Config::NONPOSITIVE_CLAMP, E1 * factor);
        E2 = max(Config::NONPOSITIVE_CLAMP, E2 * factor);
    }

    inline vector<ClusterPair> makeOptimizationSubset(const vector<ClusterPair> &events,
                                                      int max_events,
                                                      int n_bins) {
        if (max_events <= 0 || (int)events.size() <= max_events) return events;

        n_bins = max(1, n_bins);
        vector<vector<int>> bins(n_bins);
        const double lo = Config::OPT_SUBSET_MGG_MIN;
        const double hi = Config::OPT_SUBSET_MGG_MAX;
        const double width = max(1e-12, hi - lo);

        for (int i = 0; i < (int)events.size(); ++i) {
            const auto &ev = events[i];
            double mgg = nps::invariant_mass_pi0(ev.e1, ev.e2, ev.x1, ev.x2, ev.y1, ev.y2,
                                                 nps::kDefaultZ_NPS_cm);
            int ibin = 0;
            if (std::isfinite(mgg)) {
                ibin = (int)floor((mgg - lo) / width * n_bins);
            }
            ibin = max(0, min(n_bins - 1, ibin));
            bins[ibin].push_back(i);
        }

        vector<int> take(n_bins, 0);
        int total_take = 0;
        for (int ibin = 0; ibin < n_bins; ++ibin) {
            const int n = (int)bins[ibin].size();
            if (n == 0) continue;
            int t = (int)floor((double)max_events * n / (double)events.size());
            t = max(1, min(n, t));
            take[ibin] = t;
            total_take += t;
        }

        while (total_take > max_events) {
            int drop_bin = -1;
            int largest_take = 1;
            for (int ibin = 0; ibin < n_bins; ++ibin) {
                if (take[ibin] > largest_take) {
                    largest_take = take[ibin];
                    drop_bin = ibin;
                }
            }
            if (drop_bin < 0) break;
            --take[drop_bin];
            --total_take;
        }
        while (total_take < max_events) {
            int add_bin = -1;
            int largest_room = 0;
            for (int ibin = 0; ibin < n_bins; ++ibin) {
                int room = (int)bins[ibin].size() - take[ibin];
                if (room > largest_room) {
                    largest_room = room;
                    add_bin = ibin;
                }
            }
            if (add_bin < 0) break;
            ++take[add_bin];
            ++total_take;
        }

        vector<ClusterPair> subset;
        subset.reserve(total_take);
        for (int ibin = 0; ibin < n_bins; ++ibin) {
            const int n = (int)bins[ibin].size();
            const int t = take[ibin];
            if (n == 0 || t == 0) continue;

            vector<int> chosen;
            chosen.reserve(t);
            for (int j = 0; j < t; ++j) {
                int local = (int)floor(((double)j + 0.5) * n / t);
                local = max(0, min(n - 1, local));
                chosen.push_back(bins[ibin][local]);
            }

            double total_w = 0.0;
            double chosen_w = 0.0;
            for (int idx : bins[ibin]) total_w += events[idx].weight;
            for (int idx : chosen) chosen_w += events[idx].weight;
            double scale = (std::isfinite(total_w) && std::isfinite(chosen_w) && chosen_w > 0.0)
                         ? total_w / chosen_w
                         : (double)n / (double)t;

            for (int idx : chosen) {
                ClusterPair ev = events[idx];
                ev.weight *= scale;
                subset.push_back(ev);
            }
        }
        return subset;
    }

    inline void fillSmearedHistogramsAtParams(const vector<ClusterPair> &simEvents,
                                            TH1D &hmpi0,
                                            TH1D &hmmiss,
                                            TH1D &hmpgg2,
                                            double mu_a,
                                            double mu_b,
                                            double mu_c,
                                            double sigma,
                                            double sigma_pos,
                                            double p_e_scale,
                                            TRandom3 &rng,
                                            int Nsmear,
                                            double res_A,
                                            double res_B,
                                            double res_C);

    // ============================================================================
    // CHI-SQUARED CALCULATION
    // ============================================================================

    struct HistObjectiveMetrics {
        double chi2 = 1e300;
        int informative_bins = 0;
        int empty_sim_data_positive_bins = 0;
        double data_integral = 0.0;
        double sim_integral = 0.0;
    };

    struct ObjectiveBreakdown {
        HistObjectiveMetrics mpi0;
        HistObjectiveMetrics mmiss;
        HistObjectiveMetrics mpgg2;

        double total(double w_mpi0, double w_mmiss, double w_mpgg2) const {
            return w_mpi0 * mpi0.chi2 + w_mmiss * mmiss.chi2 + w_mpgg2 * mpgg2.chi2;
        }
    };

    // compute chi2 using per-bin errors: sum ( (s-d)^2 / (sigma_s^2 + sigma_d^2) )
    // Requires Sumw2() to be called on both histograms before filling so that
    // GetBinError() returns sqrt(sum_of_weights^2) rather than sqrt(N).
    // This correctly handles weighted data and weighted simulation histograms.
    // NOTE: If both errors are zero in a populated bin, use a Poisson-like fallback
    // variance max(1, s+d) to avoid both divide-by-zero and overly weak penalties.
    HistObjectiveMetrics computeChi2MetricsFromHist(const TH1D &hsim, const TH1D &hdata) {
        HistObjectiveMetrics metrics;
        metrics.data_integral = hdata.Integral();
        metrics.sim_integral = hsim.Integral();
        if (hsim.GetNbinsX() != hdata.GetNbinsX() ||
            fabs(hsim.GetXaxis()->GetXmin() - hdata.GetXaxis()->GetXmin()) > 1e-9 ||
            fabs(hsim.GetXaxis()->GetXmax() - hdata.GetXaxis()->GetXmax()) > 1e-9) {
            cerr << "Histogram binning mismatch in computeChi2FromHist" << endl;
            return metrics;
        }
        double chi2 = 0.;
        int nb = hsim.GetNbinsX();
        const double empty_sim_floor = std::max(Config::BAKER_COUSINS_EMPTY_SIM_FLOOR_ABS,
                                                Config::BAKER_COUSINS_EMPTY_SIM_FLOOR_FRAC *
                                                std::max(1.0, metrics.data_integral));

        if (Config::USE_BAKER_COUSINS_CHI2) {
            // Baker-Cousins log-likelihood ratio:
            //   chi2 = 2 * sum_i [ s_i - d_i + d_i * ln(d_i / s_i) ]
            // Handles low-count bins correctly; unbiased; follows a chi-squared
            // distribution more closely than Pearson chi2.  Does not depend on
            // simulation error bars, eliminating stochastic smearing noise.
            for (int i = 1; i <= nb; ++i) {
                double s = hsim.GetBinContent(i);
                double d = hdata.GetBinContent(i);
                if (d > 0.0) {
                    ++metrics.informative_bins;
                    double s_eval = s;
                    if (!(s_eval > 0.0)) {
                        s_eval = empty_sim_floor;
                        ++metrics.empty_sim_data_positive_bins;
                    }
                    chi2 += 2.0 * (s_eval - d + d * log(d / s_eval));
                } else if (s > 0.0) {
                    chi2 += 2.0 * s;
                }
                // Both zero: no contribution (correct).  Data-positive,
                // sim-empty bins are finite but strongly penalized above.
            }
        } else {
            // Pearson chi2 with combined sim+data errors:
            //   chi2 = sum (s-d)^2 / (sigma_s^2 + sigma_d^2)
            for (int i = 1; i <= nb; ++i) {
                double s    = hsim.GetBinContent(i);
                double d    = hdata.GetBinContent(i);
                if (d > 0.0) ++metrics.informative_bins;
                if (d > 0.0 && !(s > 0.0)) ++metrics.empty_sim_data_positive_bins;
                double es   = hsim.GetBinError(i);
                double ed   = hdata.GetBinError(i);
                double denom = es * es + ed * ed;
                if (denom <= 0.0 && (s > 0.0 || d > 0.0)) {
                    denom = max(1.0, s + d);
                }
                if (denom <= 0) continue;
                double diff = s - d;
                chi2 += (diff * diff) / denom;
            }
        }
        metrics.chi2 = chi2;
        return metrics;
    }

    double computeChi2FromHist(const TH1D &hsim, const TH1D &hdata) {
        return computeChi2MetricsFromHist(hsim, hdata).chi2;
    }

    inline int countInformativeBinsData(const TH1D &h) {
        int n = 0;
        for (int i = 1; i <= h.GetNbinsX(); ++i) {
            if (h.GetBinContent(i) > 0.0) ++n;
        }
        return max(1, n);
    }

    inline void fillUnsmearedHistogramsForNormalization(const vector<ClusterPair> &simEvents,
                                                        TH1D &hmpi0,
                                                        TH1D &hmmiss,
                                                        TH1D &hmpgg2,
                                                        double p_e_scale) {
        const double z_nps = nps::kDefaultZ_NPS_cm;
        for (const auto &ev : simEvents) {
            if (ev.e1 <= 0 || ev.e2 <= 0) continue;

            PhotonMomentum p1 = computePhotonMomentum(ev.e1, ev.x1, ev.y1, z_nps);
            PhotonMomentum p2 = computePhotonMomentum(ev.e2, ev.x2, ev.y2, z_nps);
            hmpi0.Fill(computeInvariantMass(ev.e1, p1.px, p1.py, p1.pz,
                                            ev.e2, p2.px, p2.py, p2.pz), ev.weight);

            ElectronKinematics e_kin = scaleElectronKinematicsFromMomentum(
                ev.px_e, ev.py_e, ev.pz_e, p_e_scale);
            double mmiss = nps::missing_mass_proton_pi0(
                Config::BEAM_ENERGY,
                e_kin.Ee, e_kin.px, e_kin.py, e_kin.pz,
                ev.e1, ev.e2, ev.x1, ev.y1, ev.x2, ev.y2);
            hmmiss.Fill(mmiss, ev.weight);

            double E1_mpgg2 = Config::W_MPGG2_ENERGY * ev.e1;
            double E2_mpgg2 = Config::W_MPGG2_ENERGY * ev.e2;
            PhotonMomentum p1_mpgg2 = computePhotonMomentum(E1_mpgg2, ev.x1, ev.y1, z_nps);
            PhotonMomentum p2_mpgg2 = computePhotonMomentum(E2_mpgg2, ev.x2, ev.y2, z_nps);
            hmpgg2.Fill(computeTargetPlusDiphotonMass2(E1_mpgg2, p1_mpgg2,
                                                       E2_mpgg2, p2_mpgg2), ev.weight);
        }
    }

    inline bool scaleSmearedFromUnsmearedIntegral(TH1D &hsmeared,
                                                  const TH1D &hunsmeared,
                                                  const TH1D &hdata) {
        const double integral_unsmeared = hunsmeared.Integral();
        const double integral_data = hdata.Integral();
        if (integral_unsmeared <= 0.0 || integral_data <= 0.0) return false;
        hsmeared.Scale(integral_data / integral_unsmeared);
        return true;
    }

    // ============================================================================
    // COMBINED CHI2 EVALUATION
    // ============================================================================
    ObjectiveBreakdown eval_objective_breakdown(double mu_a, double mu_b, double mu_c,
                            double sigma, double sigma_pos, double p_e_scale,
                            const vector<ClusterPair> &simEvents,
                            const TH1D &hdata_mpi0,
                            const TH1D &hdata_mmiss,
                            const TH1D &hdata_mpgg2,
                            TRandom3 &rng,
                            int Nsmear,
                            double res_A, double res_B, double res_C) {

        ObjectiveBreakdown out;
        TH1D hsim_mpi0("hsim_mpi0_breakdown_tmp","smeared sim invariant mass (tmp)",
                    hdata_mpi0.GetNbinsX(), hdata_mpi0.GetXaxis()->GetXmin(), hdata_mpi0.GetXaxis()->GetXmax());
        TH1D hsim_mmiss("hsim_mmiss_breakdown_tmp","smeared sim missing mass (tmp)",
                        hdata_mmiss.GetNbinsX(), hdata_mmiss.GetXaxis()->GetXmin(), hdata_mmiss.GetXaxis()->GetXmax());
        TH1D hsim_mpgg2("hsim_mpgg2_breakdown_tmp","smeared sim (p_{target}+#gamma#gamma)^{2} (tmp)",
                        hdata_mpgg2.GetNbinsX(), hdata_mpgg2.GetXaxis()->GetXmin(), hdata_mpgg2.GetXaxis()->GetXmax());
        TH1D hnorm_mpi0("hnorm_mpi0_breakdown_tmp","unsmeared sim invariant mass norm (tmp)",
                    hdata_mpi0.GetNbinsX(), hdata_mpi0.GetXaxis()->GetXmin(), hdata_mpi0.GetXaxis()->GetXmax());
        TH1D hnorm_mmiss("hnorm_mmiss_breakdown_tmp","unsmeared sim missing mass norm (tmp)",
                        hdata_mmiss.GetNbinsX(), hdata_mmiss.GetXaxis()->GetXmin(), hdata_mmiss.GetXaxis()->GetXmax());
        TH1D hnorm_mpgg2("hnorm_mpgg2_breakdown_tmp","unsmeared sim (p_{target}+#gamma#gamma)^{2} norm (tmp)",
                        hdata_mpgg2.GetNbinsX(), hdata_mpgg2.GetXaxis()->GetXmin(), hdata_mpgg2.GetXaxis()->GetXmax());
        hsim_mpi0.Sumw2();
        hsim_mmiss.Sumw2();
        hsim_mpgg2.Sumw2();
        hnorm_mpi0.Sumw2();
        hnorm_mmiss.Sumw2();
        hnorm_mpgg2.Sumw2();

        fillUnsmearedHistogramsForNormalization(simEvents, hnorm_mpi0, hnorm_mmiss, hnorm_mpgg2, p_e_scale);

        fillSmearedHistogramsAtParams(simEvents,
                                    hsim_mpi0, hsim_mmiss, hsim_mpgg2,
                                    mu_a, mu_b, mu_c,
                                    sigma, sigma_pos,
                                    p_e_scale,
                                    rng, Nsmear,
                                    res_A, res_B, res_C);

        if (!scaleSmearedFromUnsmearedIntegral(hsim_mpi0, hnorm_mpi0, hdata_mpi0)) return out;
        if (!scaleSmearedFromUnsmearedIntegral(hsim_mmiss, hnorm_mmiss, hdata_mmiss)) return out;
        if (!scaleSmearedFromUnsmearedIntegral(hsim_mpgg2, hnorm_mpgg2, hdata_mpgg2)) return out;

        out.mpi0 = computeChi2MetricsFromHist(hsim_mpi0, hdata_mpi0);
        out.mmiss = computeChi2MetricsFromHist(hsim_mmiss, hdata_mmiss);
        out.mpgg2 = computeChi2MetricsFromHist(hsim_mpgg2, hdata_mpgg2);
        return out;
    }

    // COMBINED CHI2: Simultaneously evaluate both M_γγ and M_miss observables
    // This accounts for correlations between observables (both depend on photon energies)
    // Returns: w_mpi0 * χ²_M_γγ + w_mmiss * χ²_M_miss
    double eval_chi2_combined(double mu_a, double mu_b, double mu_c,
                            double sigma, double sigma_pos, double p_e_scale,
                            const vector<ClusterPair> &simEvents,
                            const TH1D &hdata_mpi0,   // invariant mass data histogram
                            const TH1D &hdata_mmiss,  // missing mass data histogram
                            const TH1D &hdata_mpgg2,  // (p_target + gamma+gamma)^2 data histogram
                            TRandom3 &rng,
                            int Nsmear,
                            double res_A, double res_B, double res_C,
                            double w_mpi0 = 1.0,      // weight for invariant mass chi2
                            double w_mmiss = 1.0,     // weight for missing mass chi2
                            double w_mpgg2 = 1.0) {   // weight for (p_target + gamma+gamma)^2 chi2
        
        ObjectiveBreakdown breakdown = eval_objective_breakdown(mu_a, mu_b, mu_c,
                                    sigma, sigma_pos, p_e_scale,
                                    simEvents,
                                    hdata_mpi0, hdata_mmiss, hdata_mpgg2,
                                    rng, Nsmear,
                                    res_A, res_B, res_C);
        return breakdown.total(w_mpi0, w_mmiss, w_mpgg2);
    }

    // ============================================================================
    // INVARIANT MASS CHI2 EVALUATION
    // ============================================================================
    // Evaluate chi2 using ONLY invariant mass M_γγ
    // Used by selected-observable diagnostic evaluators.
    double eval_chi2_mpi0_only(double mu_a, double mu_b, double mu_c,
                            double sigma, double sigma_pos,
                            const vector<ClusterPair> &simEvents,
                            const TH1D &hdata_mpi0,
                            TRandom3 &rng,
                            int Nsmear,
                            double res_A, double res_B, double res_C) {

        TH1D hsim_mpi0("hsim_mpi0_stage1_tmp","smeared sim invariant mass (tmp)",
                    hdata_mpi0.GetNbinsX(), hdata_mpi0.GetXaxis()->GetXmin(), hdata_mpi0.GetXaxis()->GetXmax());
        TH1D hnorm_mpi0("hnorm_mpi0_stage1_tmp","unsmeared sim invariant mass norm (tmp)",
                    hdata_mpi0.GetNbinsX(), hdata_mpi0.GetXaxis()->GetXmin(), hdata_mpi0.GetXaxis()->GetXmax());
        hsim_mpi0.Sumw2();
        hnorm_mpi0.Sumw2();

        rng.SetSeed(Config::SMEAR_DETERMINISTIC_SEED);
        const double z_nps = nps::kDefaultZ_NPS_cm;

        for (const auto &ev : simEvents) {
            if (ev.e1 <= 0 || ev.e2 <= 0) continue;
            PhotonMomentum p1u = computePhotonMomentum(ev.e1, ev.x1, ev.y1, z_nps);
            PhotonMomentum p2u = computePhotonMomentum(ev.e2, ev.x2, ev.y2, z_nps);
            double mass_u = computeInvariantMass(ev.e1, p1u.px, p1u.py, p1u.pz,
                                                 ev.e2, p2u.px, p2u.py, p2u.pz);
            hnorm_mpi0.Fill(mass_u, ev.weight);

            double ma1   = ev.photon1_in_section ? mu_a      : ev.mu_a1_ext;
            double mb1   = ev.photon1_in_section ? mu_b      : ev.mu1_ext;
            double mc1   = ev.photon1_in_section ? mu_c      : ev.mu_c1_ext;
            double sig1  = ev.photon1_in_section ? sigma     : ev.sigma1_ext;
            double spos1 = ev.photon1_in_section ? sigma_pos : ev.sigma_pos1_ext;
            double ma2   = ev.photon2_in_section ? mu_a      : ev.mu_a2_ext;
            double mb2   = ev.photon2_in_section ? mu_b      : ev.mu2_ext;
            double mc2   = ev.photon2_in_section ? mu_c      : ev.mu_c2_ext;
            double sig2  = ev.photon2_in_section ? sigma     : ev.sigma2_ext;
            double spos2 = ev.photon2_in_section ? sigma_pos : ev.sigma_pos2_ext;

            double event_weight = ev.weight / Nsmear;
            for (int k = 0; k < Nsmear; ++k) {
                double E1_sm, x1_sm, y1_sm, E2_sm, x2_sm, y2_sm;
                smearPhoton(ev.e1, ev.x1, ev.y1, ma1, mb1, mc1, sig1, spos1,
                            res_A, res_B, res_C, rng, E1_sm, x1_sm, y1_sm);
                smearPhoton(ev.e2, ev.x2, ev.y2, ma2, mb2, mc2, sig2, spos2,
                            res_A, res_B, res_C, rng, E2_sm, x2_sm, y2_sm);
                apply_mgg_linear_pair_energy_correction(E1_sm, E2_sm, x1_sm, y1_sm, x2_sm, y2_sm);
                PhotonMomentum p1 = computePhotonMomentum(E1_sm, x1_sm, y1_sm, z_nps);
                PhotonMomentum p2 = computePhotonMomentum(E2_sm, x2_sm, y2_sm, z_nps);
                double mass = computeInvariantMass(E1_sm, p1.px, p1.py, p1.pz,
                                                E2_sm, p2.px, p2.py, p2.pz);
                hsim_mpi0.Fill(mass, event_weight);
            }
        }

        double integral_sim = hnorm_mpi0.Integral();
        double integral_data = hdata_mpi0.Integral();
        if (integral_sim <= 0 || integral_data <= 0) return 1e300;
        hsim_mpi0.Scale(integral_data / integral_sim);

        return computeChi2FromHist(hsim_mpi0, hdata_mpi0);
    }

    // ============================================================================
    // MISSING MASS CHI2 EVALUATION
    // ============================================================================
    // Evaluate chi2 using ONLY missing mass M_miss
    // Used by the selected-observable dispatcher and optional p_e_scale refinement.
    double eval_chi2_mmiss_only(double mu_a, double mu_b, double mu_c,
                                double sigma, double sigma_pos, double p_e_scale,
                                const vector<ClusterPair> &simEvents,
                                const TH1D &hdata_mmiss,
                                TRandom3 &rng,
                                int Nsmear,
                                double res_A, double res_B, double res_C) {

        TH1D hsim_mmiss("hsim_mmiss_stage2_tmp","smeared sim missing mass (tmp)",
                        hdata_mmiss.GetNbinsX(), hdata_mmiss.GetXaxis()->GetXmin(), hdata_mmiss.GetXaxis()->GetXmax());
        TH1D hnorm_mmiss("hnorm_mmiss_stage2_tmp","unsmeared sim missing mass norm (tmp)",
                        hdata_mmiss.GetNbinsX(), hdata_mmiss.GetXaxis()->GetXmin(), hdata_mmiss.GetXaxis()->GetXmax());
        hsim_mmiss.Sumw2();
        hnorm_mmiss.Sumw2();

        rng.SetSeed(Config::SMEAR_DETERMINISTIC_SEED);

        for (const auto &ev : simEvents) {
            if (ev.e1 <= 0 || ev.e2 <= 0) continue;
            ElectronKinematics e_kin_u = scaleElectronKinematicsFromMomentum(
                ev.px_e, ev.py_e, ev.pz_e, p_e_scale);
            double mmiss_u = nps::missing_mass_proton_pi0(
                Config::BEAM_ENERGY,
                e_kin_u.Ee, e_kin_u.px, e_kin_u.py, e_kin_u.pz,
                ev.e1, ev.e2, ev.x1, ev.y1, ev.x2, ev.y2);
            hnorm_mmiss.Fill(mmiss_u, ev.weight);

            double ma1   = ev.photon1_in_section ? mu_a      : ev.mu_a1_ext;
            double mb1   = ev.photon1_in_section ? mu_b      : ev.mu1_ext;
            double mc1   = ev.photon1_in_section ? mu_c      : ev.mu_c1_ext;
            double sig1  = ev.photon1_in_section ? sigma     : ev.sigma1_ext;
            double spos1 = ev.photon1_in_section ? sigma_pos : ev.sigma_pos1_ext;
            double ma2   = ev.photon2_in_section ? mu_a      : ev.mu_a2_ext;
            double mb2   = ev.photon2_in_section ? mu_b      : ev.mu2_ext;
            double mc2   = ev.photon2_in_section ? mu_c      : ev.mu_c2_ext;
            double sig2  = ev.photon2_in_section ? sigma     : ev.sigma2_ext;
            double spos2 = ev.photon2_in_section ? sigma_pos : ev.sigma_pos2_ext;

            double event_weight = ev.weight / Nsmear;
            for (int k = 0; k < Nsmear; ++k) {
                double E1_sm, x1_sm, y1_sm, E2_sm, x2_sm, y2_sm;
                smearPhoton(ev.e1, ev.x1, ev.y1, ma1, mb1, mc1, sig1, spos1,
                            res_A, res_B, res_C, rng, E1_sm, x1_sm, y1_sm);
                smearPhoton(ev.e2, ev.x2, ev.y2, ma2, mb2, mc2, sig2, spos2,
                            res_A, res_B, res_C, rng, E2_sm, x2_sm, y2_sm);
                apply_mgg_linear_pair_energy_correction(E1_sm, E2_sm, x1_sm, y1_sm, x2_sm, y2_sm);
                ElectronKinematics e_kin = scaleElectronKinematicsFromMomentum(
                    ev.px_e, ev.py_e, ev.pz_e, p_e_scale);
                double mmiss = nps::missing_mass_proton_pi0(
                    Config::BEAM_ENERGY,
                    e_kin.Ee, e_kin.px, e_kin.py, e_kin.pz,
                    E1_sm, E2_sm, x1_sm, y1_sm, x2_sm, y2_sm
                );
                hsim_mmiss.Fill(mmiss, event_weight);
            }
        }

        double integral_sim = hnorm_mmiss.Integral();
        double integral_data = hdata_mmiss.Integral();
        if (integral_sim <= 0 || integral_data <= 0) return 1e300;
        hsim_mmiss.Scale(integral_data / integral_sim);

        return computeChi2FromHist(hsim_mmiss, hdata_mmiss);
    }

    // ============================================================================
    // TARGET+2γ MASS-SQUARED CHI2 EVALUATION
    // ============================================================================
    double eval_chi2_mpgg2_only(double mu_a, double mu_b, double mu_c,
                                double sigma, double sigma_pos,
                                const vector<ClusterPair> &simEvents,
                                const TH1D &hdata_mpgg2,
                                TRandom3 &rng,
                                int Nsmear,
                                double res_A, double res_B, double res_C) {

        TH1D hsim_mpgg2("hsim_mpgg2_stage2_tmp","smeared sim (p_{target}+#gamma#gamma)^{2} (tmp)",
                        hdata_mpgg2.GetNbinsX(), hdata_mpgg2.GetXaxis()->GetXmin(), hdata_mpgg2.GetXaxis()->GetXmax());
        TH1D hnorm_mpgg2("hnorm_mpgg2_stage2_tmp","unsmeared sim (p_{target}+#gamma#gamma)^{2} norm (tmp)",
                        hdata_mpgg2.GetNbinsX(), hdata_mpgg2.GetXaxis()->GetXmin(), hdata_mpgg2.GetXaxis()->GetXmax());
        hsim_mpgg2.Sumw2();
        hnorm_mpgg2.Sumw2();

        rng.SetSeed(Config::SMEAR_DETERMINISTIC_SEED);
        const double z_nps = nps::kDefaultZ_NPS_cm;

        for (const auto &ev : simEvents) {
            if (ev.e1 <= 0 || ev.e2 <= 0) continue;
            double E1u_mpgg2 = Config::W_MPGG2_ENERGY * ev.e1;
            double E2u_mpgg2 = Config::W_MPGG2_ENERGY * ev.e2;
            PhotonMomentum p1u_mpgg2 = computePhotonMomentum(E1u_mpgg2, ev.x1, ev.y1, z_nps);
            PhotonMomentum p2u_mpgg2 = computePhotonMomentum(E2u_mpgg2, ev.x2, ev.y2, z_nps);
            double mpgg2_u = computeTargetPlusDiphotonMass2(E1u_mpgg2, p1u_mpgg2,
                                                            E2u_mpgg2, p2u_mpgg2);
            hnorm_mpgg2.Fill(mpgg2_u, ev.weight);

            double ma1   = ev.photon1_in_section ? mu_a      : ev.mu_a1_ext;
            double mb1   = ev.photon1_in_section ? mu_b      : ev.mu1_ext;
            double mc1   = ev.photon1_in_section ? mu_c      : ev.mu_c1_ext;
            double sig1  = ev.photon1_in_section ? sigma     : ev.sigma1_ext;
            double spos1 = ev.photon1_in_section ? sigma_pos : ev.sigma_pos1_ext;
            double ma2   = ev.photon2_in_section ? mu_a      : ev.mu_a2_ext;
            double mb2   = ev.photon2_in_section ? mu_b      : ev.mu2_ext;
            double mc2   = ev.photon2_in_section ? mu_c      : ev.mu_c2_ext;
            double sig2  = ev.photon2_in_section ? sigma     : ev.sigma2_ext;
            double spos2 = ev.photon2_in_section ? sigma_pos : ev.sigma_pos2_ext;

            double event_weight = ev.weight / Nsmear;
            for (int k = 0; k < Nsmear; ++k) {
                double E1_sm, x1_sm, y1_sm, E2_sm, x2_sm, y2_sm;
                smearPhoton(ev.e1, ev.x1, ev.y1, ma1, mb1, mc1, sig1, spos1,
                            res_A, res_B, res_C, rng, E1_sm, x1_sm, y1_sm);
                smearPhoton(ev.e2, ev.x2, ev.y2, ma2, mb2, mc2, sig2, spos2,
                            res_A, res_B, res_C, rng, E2_sm, x2_sm, y2_sm);
                apply_mgg_linear_pair_energy_correction(E1_sm, E2_sm, x1_sm, y1_sm, x2_sm, y2_sm);
                double E1_mpgg2 = Config::W_MPGG2_ENERGY * E1_sm;
                double E2_mpgg2 = Config::W_MPGG2_ENERGY * E2_sm;
                PhotonMomentum p1_mpgg2 = computePhotonMomentum(E1_mpgg2, x1_sm, y1_sm, z_nps);
                PhotonMomentum p2_mpgg2 = computePhotonMomentum(E2_mpgg2, x2_sm, y2_sm, z_nps);
                double mpgg2 = computeTargetPlusDiphotonMass2(E1_mpgg2, p1_mpgg2, E2_mpgg2, p2_mpgg2);
                hsim_mpgg2.Fill(mpgg2, event_weight);
            }
        }

        double integral_sim = hnorm_mpgg2.Integral();
        double integral_data = hdata_mpgg2.Integral();
        if (integral_sim <= 0 || integral_data <= 0) return 1e300;
        hsim_mpgg2.Scale(integral_data / integral_sim);

        return computeChi2FromHist(hsim_mpgg2, hdata_mpgg2);
    }

    // ============================================================================
    // CHI2 EVALUATION DISPATCHER
    // ============================================================================
    // Wrapper function that calls the appropriate chi2 evaluation function
    // based on Config::ENERGY_SMEARING_HISTOGRAM setting:
    //   0 = combined histograms (M_γγ + M_miss + (p_target+γγ)^2)
    //   1 = invariant mass only (M_γγ)
    //   2 = missing mass only (M_miss)
    //   3 = (p_target+γγ)^2 only
    double eval_chi2_selected(double mu_a, double mu_b, double mu_c,
                            double sigma, double sigma_pos, double p_e_scale,
                            const vector<ClusterPair> &simEvents,
                            const TH1D &hdata_mpi0,
                            const TH1D &hdata_mmiss,
                            const TH1D &hdata_mpgg2,
                            TRandom3 &rng,
                            int Nsmear,
                            double res_A, double res_B, double res_C,
                            double w_mpi0 = 1.0,
                            double w_mmiss = 1.0,
                            double w_mpgg2 = 1.0) {

        if (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_MPI0_ONLY) {
            return eval_chi2_mpi0_only(mu_a, mu_b, mu_c, sigma, sigma_pos, simEvents,
                                       hdata_mpi0, rng, Nsmear, res_A, res_B, res_C);
        } else if (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_MMISS_ONLY) {
            return eval_chi2_mmiss_only(mu_a, mu_b, mu_c, sigma, sigma_pos, p_e_scale, simEvents,
                                        hdata_mmiss, rng, Nsmear, res_A, res_B, res_C);
        } else if (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_MPGG2_ONLY) {
            return eval_chi2_mpgg2_only(mu_a, mu_b, mu_c, sigma, sigma_pos, simEvents,
                                        hdata_mpgg2, rng, Nsmear, res_A, res_B, res_C);
        } else {
            return eval_chi2_combined(mu_a, mu_b, mu_c, sigma, sigma_pos, p_e_scale, simEvents,
                                      hdata_mpi0, hdata_mmiss, hdata_mpgg2, rng, Nsmear,
                                      res_A, res_B, res_C, w_mpi0, w_mmiss, w_mpgg2);
        }
    }

    // ============================================================================
    // SHARED SMEARING ROUTINE (for diagnostics + consistency checks)
    // ============================================================================
    // Fills M_γγ, M_miss, and M_pgg2 histograms using the exact same smearing
    // path as eval_chi2_combined: per-photon in/out-of-section params, smearPhoton,
    // weighted fills with ev.weight/Nsmear, and p_e_scale in missing mass.
    inline void fillSmearedHistogramsAtParams(const vector<ClusterPair> &simEvents,
                                            TH1D &hsim_mpi0,
                                            TH1D &hsim_mmiss,
                                            TH1D &hsim_mpgg2,
                                            double mu_a, double mu_b, double mu_c,
                                            double sigma, double sigma_pos,
                                            double p_e_scale,
                                            TRandom3 &rng,
                                            int Nsmear,
                                            double res_A, double res_B, double res_C) {
        rng.SetSeed(Config::SMEAR_DETERMINISTIC_SEED);
        const double z_nps = nps::kDefaultZ_NPS_cm;
        for (const auto &ev : simEvents) {
            if (ev.e1 <= 0 || ev.e2 <= 0) continue;

            double ma1   = ev.photon1_in_section ? mu_a      : ev.mu_a1_ext;
            double mb1   = ev.photon1_in_section ? mu_b      : ev.mu1_ext;
            double mc1   = ev.photon1_in_section ? mu_c      : ev.mu_c1_ext;
            double sig1  = ev.photon1_in_section ? sigma     : ev.sigma1_ext;
            double spos1 = ev.photon1_in_section ? sigma_pos : ev.sigma_pos1_ext;
            double ma2   = ev.photon2_in_section ? mu_a      : ev.mu_a2_ext;
            double mb2   = ev.photon2_in_section ? mu_b      : ev.mu2_ext;
            double mc2   = ev.photon2_in_section ? mu_c      : ev.mu_c2_ext;
            double sig2  = ev.photon2_in_section ? sigma     : ev.sigma2_ext;
            double spos2 = ev.photon2_in_section ? sigma_pos : ev.sigma_pos2_ext;

            double event_weight = ev.weight / Nsmear;
            for (int k = 0; k < Nsmear; ++k) {
                double E1_sm, x1_sm, y1_sm, E2_sm, x2_sm, y2_sm;
                smearPhoton(ev.e1, ev.x1, ev.y1, ma1, mb1, mc1, sig1, spos1,
                            res_A, res_B, res_C, rng, E1_sm, x1_sm, y1_sm);
                smearPhoton(ev.e2, ev.x2, ev.y2, ma2, mb2, mc2, sig2, spos2,
                            res_A, res_B, res_C, rng, E2_sm, x2_sm, y2_sm);
                apply_mgg_linear_pair_energy_correction(E1_sm, E2_sm, x1_sm, y1_sm, x2_sm, y2_sm);

                PhotonMomentum p1 = computePhotonMomentum(E1_sm, x1_sm, y1_sm, z_nps);
                PhotonMomentum p2 = computePhotonMomentum(E2_sm, x2_sm, y2_sm, z_nps);
                double mass = computeInvariantMass(E1_sm, p1.px, p1.py, p1.pz,
                                                E2_sm, p2.px, p2.py, p2.pz);
                hsim_mpi0.Fill(mass, event_weight);

                ElectronKinematics e_kin = scaleElectronKinematicsFromMomentum(
                    ev.px_e, ev.py_e, ev.pz_e, p_e_scale);
                double mmiss = nps::missing_mass_proton_pi0(
                    Config::BEAM_ENERGY,
                    e_kin.Ee, e_kin.px, e_kin.py, e_kin.pz,
                    E1_sm, E2_sm, x1_sm, y1_sm, x2_sm, y2_sm
                );
                hsim_mmiss.Fill(mmiss, event_weight);

                double E1_mpgg2 = Config::W_MPGG2_ENERGY * E1_sm;
                double E2_mpgg2 = Config::W_MPGG2_ENERGY * E2_sm;
                PhotonMomentum p1_mpgg2 = computePhotonMomentum(E1_mpgg2, x1_sm, y1_sm, z_nps);
                PhotonMomentum p2_mpgg2 = computePhotonMomentum(E2_mpgg2, x2_sm, y2_sm, z_nps);
                double mpgg2 = computeTargetPlusDiphotonMass2(E1_mpgg2, p1_mpgg2, E2_mpgg2, p2_mpgg2);
                hsim_mpgg2.Fill(mpgg2, event_weight);
            }
        }
    }

    // ============================================================================
    // CalibrationMap: Provides bilinear-interpolated smearing parameters across calorimeter
    // ============================================================================
    class CalibrationMap {
    private:
        int nx_, ny_;
        double x_min_, x_max_, y_min_, y_max_;
        double base_wx_, base_wy_;
        
        // Storage for section centers and fitted parameters
        vector<double> centers_x_;  // section center x coordinates
        vector<double> centers_y_;  // section center y coordinates
        // mu_eff(E) = mu_a + mu_b*E + mu_c*ln(E);  disabled: mu_a=0, mu_c=0
        vector<double> mu_a_values_;
        vector<double> mu_values_;   // = mu_b
        vector<double> mu_c_values_;
        vector<double> sigma_values_;
        vector<double> sigma_pos_values_;

        // GP/Kriging-style interpolation cache
        mutable bool gp_model_ready_;
        mutable bool gp_model_valid_;
        mutable vector<double> gp_alpha_mu_a_;
        mutable vector<double> gp_alpha_mu_;   // = mu_b
        mutable vector<double> gp_alpha_mu_c_;
        mutable vector<double> gp_alpha_sigma_;
        mutable vector<double> gp_alpha_sigma_pos_;
        double gp_lx_;
        double gp_ly_;
        double gp_nugget_;
        
        // Get linear index from (ix, iy)
        int getIndex(int ix, int iy) const {
            return iy * nx_ + ix;
        }
        
        // Clamp value to range [min, max]
        double clamp(double val, double min, double max) const {
            return std::max(min, std::min(max, val));
        }

        // Squared-exponential kernel with anisotropic length scales
        double kernel(double x1, double y1, double x2, double y2) const {
            double dx = x1 - x2;
            double dy = y1 - y2;
            double r2 = (dx * dx) / (gp_lx_ * gp_lx_) + (dy * dy) / (gp_ly_ * gp_ly_);
            return std::exp(-0.5 * r2);
        }

        static bool solve_linear_system(vector<double> A, vector<double> b, vector<double> &x, int n) {
            if ((int)b.size() != n || (int)A.size() != n * n) return false;

            for (int col = 0; col < n; ++col) {
                int pivot = col;
                double pivot_abs = std::fabs(A[col * n + col]);
                for (int row = col + 1; row < n; ++row) {
                    double cand = std::fabs(A[row * n + col]);
                    if (cand > pivot_abs) {
                        pivot_abs = cand;
                        pivot = row;
                    }
                }

                if (pivot_abs < 1e-14) return false;

                if (pivot != col) {
                    for (int j = col; j < n; ++j) {
                        std::swap(A[col * n + j], A[pivot * n + j]);
                    }
                    std::swap(b[col], b[pivot]);
                }

                double diag = A[col * n + col];
                for (int row = col + 1; row < n; ++row) {
                    double factor = A[row * n + col] / diag;
                    if (factor == 0.0) continue;
                    for (int j = col; j < n; ++j) {
                        A[row * n + j] -= factor * A[col * n + j];
                    }
                    b[row] -= factor * b[col];
                }
            }

            x.assign(n, 0.0);
            for (int row = n - 1; row >= 0; --row) {
                double sum = b[row];
                for (int j = row + 1; j < n; ++j) {
                    sum -= A[row * n + j] * x[j];
                }
                double diag = A[row * n + row];
                if (std::fabs(diag) < 1e-14) return false;
                x[row] = sum / diag;
            }
            return true;
        }

        void compute_bilinear_fallback(double x, double y,
                                       double &mu_a, double &mu, double &mu_c,
                                       double &sigma, double &sigma_pos) const {
            x = clamp(x, x_min_, x_max_);
            y = clamp(y, y_min_, y_max_);

            double fx = (x - x_min_) / base_wx_;
            double fy = (y - y_min_) / base_wy_;

            int ix0 = (int)floor(fx);
            int iy0 = (int)floor(fy);
            ix0 = max(0, min(nx_ - 1, ix0));
            iy0 = max(0, min(ny_ - 1, iy0));

            int ix1 = min(nx_ - 1, ix0 + 1);
            int iy1 = min(ny_ - 1, iy0 + 1);

            double tx = clamp(fx - ix0, 0.0, 1.0);
            double ty = clamp(fy - iy0, 0.0, 1.0);

            int idx00 = getIndex(ix0, iy0);
            int idx10 = getIndex(ix1, iy0);
            int idx01 = getIndex(ix0, iy1);
            int idx11 = getIndex(ix1, iy1);

            auto blerp = [&](const vector<double> &v) {
                return (1-tx)*(1-ty)*v[idx00] + tx*(1-ty)*v[idx10]
                     + (1-tx)*ty   *v[idx01] + tx*ty   *v[idx11];
            };
            mu_a = blerp(mu_a_values_);
            mu   = blerp(mu_values_);
            mu_c = blerp(mu_c_values_);
            sigma     = blerp(sigma_values_);
            sigma_pos = blerp(sigma_pos_values_);
        }

        void ensure_gp_model() const {
            if (gp_model_ready_) return;

            gp_model_ready_ = true;
            gp_model_valid_ = false;

            int nsec = nx_ * ny_;
            if (nsec <= 0) return;

            vector<double> K(nsec * nsec, 0.0);
            for (int i = 0; i < nsec; ++i) {
                for (int j = 0; j < nsec; ++j) {
                    K[i * nsec + j] = kernel(centers_x_[i], centers_y_[i], centers_x_[j], centers_y_[j]);
                }
                K[i * nsec + i] += gp_nugget_;
            }

            vector<double> alpha_mu_a, alpha_mu, alpha_mu_c, alpha_sigma, alpha_sigma_pos;
            bool ok_mu_a = solve_linear_system(K, mu_a_values_, alpha_mu_a, nsec);
            bool ok_mu   = solve_linear_system(K, mu_values_,   alpha_mu,   nsec);
            bool ok_mu_c = solve_linear_system(K, mu_c_values_, alpha_mu_c, nsec);
            bool ok_sigma     = solve_linear_system(K, sigma_values_,     alpha_sigma,     nsec);
            bool ok_sigma_pos = solve_linear_system(K, sigma_pos_values_, alpha_sigma_pos, nsec);

            if (!(ok_mu_a && ok_mu && ok_mu_c && ok_sigma && ok_sigma_pos)) return;

            gp_alpha_mu_a_ = alpha_mu_a;
            gp_alpha_mu_   = alpha_mu;
            gp_alpha_mu_c_ = alpha_mu_c;
            gp_alpha_sigma_     = alpha_sigma;
            gp_alpha_sigma_pos_ = alpha_sigma_pos;
            gp_model_valid_ = true;
        }

    public:
        CalibrationMap(int nx, int ny, double x_min, double x_max, double y_min, double y_max)
            : nx_(nx), ny_(ny), x_min_(x_min), x_max_(x_max), y_min_(y_min), y_max_(y_max),
            gp_model_ready_(false), gp_model_valid_(false), gp_lx_(0.0), gp_ly_(0.0), gp_nugget_(1e-6)
        {
            double total_wx = x_max - x_min;
            double total_wy = y_max - y_min;
            base_wx_ = total_wx / nx;
            base_wy_ = total_wy / ny;
            gp_lx_ = 1.5 * base_wx_;
            gp_ly_ = 1.5 * base_wy_;

            int nsec = nx * ny;
            centers_x_.resize(nsec);
            centers_y_.resize(nsec);
            mu_a_values_.resize(nsec, 0.0);    // default: mu_eff(E)=1*E (no shift)
            mu_values_.resize(nsec, 1.0);      // = mu_b; default: b=1
            mu_c_values_.resize(nsec, 0.0);    // default: no log term
            sigma_values_.resize(nsec, 0.05);
            sigma_pos_values_.resize(nsec, 0.0);
            
            // Compute section centers
            for (int iy = 0; iy < ny; ++iy) {
                for (int ix = 0; ix < nx; ++ix) {
                    int idx = getIndex(ix, iy);
                    centers_x_[idx] = x_min + (ix + 0.5) * base_wx_;
                    centers_y_[idx] = y_min + (iy + 0.5) * base_wy_;
                }
            }
        }
        
        // Set fitted parameters for a specific section (full mu_a/mu_b/mu_c form)
        void setParams(int ix, int iy, double mu_a, double mu_b, double mu_c,
                       double sigma, double sigma_pos = 0.0) {
            int idx = getIndex(ix, iy);
            mu_a_values_[idx] = mu_a;
            mu_values_[idx]   = mu_b;
            mu_c_values_[idx] = mu_c;
            sigma_values_[idx]     = sigma;
            sigma_pos_values_[idx] = sigma_pos;
            gp_model_ready_ = false;
            gp_model_valid_ = false;
        }
        // Convenience overload: disabled mode — a=0, c=0
        void setParams(int ix, int iy, double mu_b, double sigma, double sigma_pos = 0.0) {
            setParams(ix, iy, 0.0, mu_b, 0.0, sigma, sigma_pos);
        }
        
        // Get interpolated parameters at position (x, y) — full (a, b, c) form
        void getInterpolatedParams(double x, double y,
                                   double &mu_a, double &mu_b, double &mu_c,
                                   double &sigma, double &sigma_pos) const {
            compute_bilinear_fallback(x, y, mu_a, mu_b, mu_c, sigma, sigma_pos);
            sigma     = std::max(1e-6, sigma);
            sigma_pos = std::max(0.0, sigma_pos);
        }
        // Convenience overload: returns only (mu_b, sigma) — disabled-mode callers
        void getInterpolatedParams(double x, double y, double &mu, double &sigma) const {
            double mu_a_dummy = 0.0, mu_c_dummy = 0.0, sigma_pos_dummy = 0.0;
            getInterpolatedParams(x, y, mu_a_dummy, mu, mu_c_dummy, sigma, sigma_pos_dummy);
        }
        // Convenience overload: (mu_b, sigma, sigma_pos)
        void getInterpolatedParams(double x, double y, double &mu, double &sigma, double &sigma_pos) const {
            double mu_a_dummy = 0.0, mu_c_dummy = 0.0;
            getInterpolatedParams(x, y, mu_a_dummy, mu, mu_c_dummy, sigma, sigma_pos);
        }
        
        // Load calibration from CSV file
        bool loadFromCSV(const string &filename) {
            ifstream file(filename);
            if (!file.is_open()) {
                cerr << "Error: Cannot open calibration file " << filename << endl;
                return false;
            }
            
            string line;
            getline(file, line); // skip header
            
            int count = 0;
            while (getline(file, line)) {
                if (line.empty()) continue;
                istringstream ss(line);
                string token;
                vector<string> tokens;
                while (getline(ss, token, ',')) {
                    tokens.push_back(token);
                }
                
                if (tokens.size() < 14) continue;

                int ix = atoi(tokens[0].c_str());
                int iy = atoi(tokens[1].c_str());

                // CSV columns: ix,iy,xlo,xhi,ylo,yhi,n_data,n_sim,
                //              best_mu_a(8),best_mu_b(9),best_mu_c(10),best_sigma(11),best_sigma_pos(12),...
                if (tokens[8].empty() || tokens[9].empty() || tokens[10].empty() ||
                    tokens[11].empty() || tokens[12].empty()) {
                    continue;
                }

                double mu_a     = atof(tokens[8].c_str());
                double mu_b     = atof(tokens[9].c_str());
                double mu_c     = atof(tokens[10].c_str());
                double sigma    = atof(tokens[11].c_str());
                double sigma_pos = atof(tokens[12].c_str());

                setParams(ix, iy, mu_a, mu_b, mu_c, sigma, sigma_pos);
                count++;
            }
            
            file.close();
            cout << "Loaded calibration for " << count << " sections from " << filename << endl;
            return true;
        }
        
        // Save interpolated map to 2D histogram (for visualization)
        void saveAsHistogram(const string &filename,
                             const string &run_tag = "",
                             const string &created_at_local = "") const {
            TFile fout(filename.c_str(), "RECREATE");
            fout.cd();
            if (!run_tag.empty()) TNamed("run_tag", run_tag.c_str()).Write();
            if (!created_at_local.empty()) TNamed("created_at_local", created_at_local.c_str()).Write();
            TNamed("interpolated_output_file", filename.c_str()).Write();

            int nbinsx = 100, nbinsy = 100;

            TH2D *h_mu_a = new TH2D("h_mu_a_interp",
                "Interpolated #mu_{a} map (const offset);x [cm];y [cm]",
                nbinsx, x_min_, x_max_, nbinsy, y_min_, y_max_);
	            TH2D *h_mu = new TH2D("h_mu_interp",
	                "Interpolated #mu_{b} map (linear coeff);x [cm];y [cm]",
	                nbinsx, x_min_, x_max_, nbinsy, y_min_, y_max_);
	            TH2D *h_mu_b = new TH2D("h_mu_b_interp",
	                "Interpolated #mu_{b} map (linear coeff);x [cm];y [cm]",
	                nbinsx, x_min_, x_max_, nbinsy, y_min_, y_max_);
            TH2D *h_mu_c = new TH2D("h_mu_c_interp",
                "Interpolated #mu_{c} map (log coeff);x [cm];y [cm]",
                nbinsx, x_min_, x_max_, nbinsy, y_min_, y_max_);
            TH2D *h_sigma = new TH2D("h_sigma_interp",
                "Interpolated #sigma map;x [cm];y [cm]",
                nbinsx, x_min_, x_max_, nbinsy, y_min_, y_max_);
            TH2D *h_sigma_pos = new TH2D("h_sigma_pos_interp",
                "Interpolated #sigma_{pos} map;x [cm];y [cm]",
                nbinsx, x_min_, x_max_, nbinsy, y_min_, y_max_);
            vector<TH2D*> response_ratio_maps;
            const double E_fixed[] = {1.0, 2.0, 3.0, 4.0, 5.0};
            const int N_E_fixed = 5;
            for (int ie = 0; ie < N_E_fixed; ++ie) {
                double E = E_fixed[ie];
                response_ratio_maps.push_back(new TH2D(
                    Form("h_response_ratio_interp_E%.0fGeV", E),
                    Form("Interpolated energy response #mu_{eff}/E at E=%.1f GeV;x [cm];y [cm];#mu_{eff}/E", E),
                    nbinsx, x_min_, x_max_, nbinsy, y_min_, y_max_));
            }

            for (int ix = 1; ix <= nbinsx; ++ix) {
                for (int iy = 1; iy <= nbinsy; ++iy) {
                    double x = h_mu->GetXaxis()->GetBinCenter(ix);
                    double y = h_mu->GetYaxis()->GetBinCenter(iy);

                    double ma, mb, mc, sigma, sigma_pos;
                    getInterpolatedParams(x, y, ma, mb, mc, sigma, sigma_pos);

	                    h_mu_a->SetBinContent(ix, iy, ma);
	                    h_mu->SetBinContent(ix, iy, mb);
	                    h_mu_b->SetBinContent(ix, iy, mb);
	                    h_mu_c->SetBinContent(ix, iy, mc);
	                    h_sigma->SetBinContent(ix, iy, sigma);
	                    h_sigma_pos->SetBinContent(ix, iy, sigma_pos);
                    for (int ie = 0; ie < N_E_fixed; ++ie) {
                        double E = E_fixed[ie];
                        double E_s = std::max(E, Config::MU_ENERGY_MIN_GEV);
                        double mu_eff = ma + mb * E_s + mc * std::log(E_s);
                        response_ratio_maps[ie]->SetBinContent(ix, iy, (E_s > 0.0) ? mu_eff / E_s : 0.0);
                    }
	                }
	            }

	            TNamed("smearing_model_energy_mean", "a_plus_bE_plus_clnE").Write();
	            TNamed("energy_mean_convention", "reconstructed_energy_GeV").Write();
	            TNamed("energy_mean_formula", "E_mean = a + b*E_safe + c*ln(E_safe/1 GeV)").Write();
	            TNamed("energy_log_floor_GeV", Form("%.8g", Config::MU_ENERGY_MIN_GEV)).Write();
	            h_mu_a->Write(); h_mu->Write(); h_mu_b->Write(); h_mu_c->Write();
	            h_sigma->Write(); h_sigma_pos->Write();
            for (TH2D *h : response_ratio_maps) h->Write();

            TDirectory *canvas_dir = fout.mkdir("interpolated_canvases");
            TCanvas *c_params = new TCanvas("c_interpolated_parameter_maps", "Interpolated parameter maps", 1800, 900);
            c_params->Divide(3, 2);
            vector<TH2D*> param_maps = {h_mu_a, h_mu_b, h_mu_c, h_sigma, h_sigma_pos};
            for (int i = 0; i < (int)param_maps.size(); ++i) {
                c_params->cd(i + 1);
                gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.10); gPad->SetBottomMargin(0.12);
                param_maps[i]->SetStats(0);
                param_maps[i]->SetMarkerSize(0.6);
                param_maps[i]->Draw("COLZ");
            }
            writeCanvasToDir(canvas_dir, c_params, "c_interpolated_parameter_maps");

            TCanvas *c_resp = new TCanvas("c_interpolated_response_ratio_maps",
                                          "Interpolated energy response ratio maps", 2200, 800);
            c_resp->Divide(N_E_fixed, 1);
            for (int ie = 0; ie < N_E_fixed; ++ie) {
                c_resp->cd(ie + 1);
                gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.10); gPad->SetBottomMargin(0.14);
                response_ratio_maps[ie]->SetStats(0);
                response_ratio_maps[ie]->SetMarkerSize(0.6);
                response_ratio_maps[ie]->GetZaxis()->SetTitle("#mu_{eff}/E");
                response_ratio_maps[ie]->Draw("COLZ");
            }
            writeCanvasToDir(canvas_dir, c_resp, "c_interpolated_response_ratio_maps");

	            fout.Close();

            cout << "Saved interpolated calibration maps to " << filename << endl;
            delete c_params;
            delete c_resp;
            for (TH2D *h : response_ratio_maps) delete h;
        }
        
        // Print parameters at a specific position
        void printParamsAt(double x, double y) const {
            double mu, sigma, sigma_pos;
            getInterpolatedParams(x, y, mu, sigma, sigma_pos);
            cout << "Calibration at (" << x << ", " << y << "): "
                << "mu=" << mu << ", sigma=" << sigma << ", sigma_pos=" << sigma_pos << endl;
        }
    };

    // Structure to store chi2 scan data for plotting
    struct Chi2ScanData {
        vector<double> mu_values;
        vector<double> chi2_vs_mu;
        vector<double> sigma_values;
        vector<double> chi2_vs_sigma;
        vector<double> mu_2d;
        vector<double> sigma_2d;
        vector<double> chi2_2d;
        vector<double> sigma_pos_values;
        vector<double> chi2_vs_sigma_pos;
    };

    struct SeedDiagnostic {
        int rank = -1;
        int seed_index = -1;
        double seed_chi2 = 1e300;
        double migrad_chi2 = 1e300;
        double mu_a = 0.0;
        double mu = 1.0;
        double mu_c = 0.0;
        double sigma = 0.05;
        double sigma_pos = 0.0;
        double p_e_scale = 1.0;
        bool minimized = false;
        bool hesse_ok = false;
        double max_abs_corr = 0.0;
        string max_corr_pair = "";
    };

    struct ProfileDiagnostic {
        string parameter;
        double fixed_value = 0.0;
        double chi2 = 1e300;
        bool minimized = false;
    };

    // Structure to store fit results
    struct FitResult {
        // mu_eff(E) = mu_a + mu * E + mu_c * ln(E_safe)
        // Disabled mode: mu_a = 0 (fixed), mu_c = 0 (fixed), only mu (=b) is fitted.
        // Enabled  mode: all three are per-section fitted parameters.
        double mu_a = 0.0;  // constant offset  (GeV; forced 0 in disabled mode)
        double mu   = 1.0;  // linear coeff     (dimensionless; this is the scalar "b")
        double mu_c = 0.0;  // logarithmic coeff (GeV; forced 0 in disabled mode)
        double sigma = 0.05, sigma_pos = 0.0, p_e_scale = 1.0, chi2 = 1e300;
        // 3-term model parameters (only used if USE_SIMPLE_STOCHASTIC_MODEL = false)
        double res_A = Config::RESOLUTION_A_DEFAULT;
        double res_B = Config::RESOLUTION_B_DEFAULT;
        double res_C = Config::RESOLUTION_C_DEFAULT;
        int ndf = 1;
        Chi2ScanData scan_data;
        vector<SeedDiagnostic> seed_diagnostics;
        vector<ProfileDiagnostic> profile_diagnostics;
        string optimizer_mode = "unset";
        int n_seed_trials = 0;
        int n_migrad_trials = 0;
        bool hesse_ok = false;
        double max_abs_corr = 0.0;
        string max_corr_pair = "";
        double chi2_per_ndf() const { return (ndf > 0) ? chi2 / ndf : 0.0; }
        // Evaluate effective scale at energy E: mu_eff(E) = a + b*E + c*ln(E_safe)
        double muEff(double E) const {
            double E_safe = std::max(E, Config::MU_ENERGY_MIN_GEV);
            return mu_a + mu * E_safe + mu_c * std::log(E_safe);
        }
    };

    inline void generate_visualization_scan_data(const vector<ClusterPair> &simEvents,
                                                const TH1D &hdata_mpi0,
                                                const TH1D &hdata_mmiss,
                                                const TH1D &hdata_mpgg2,
                                                TRandom3 &rng,
                                                int Nsmear,
                                                FitResult &res) {
        res.scan_data = Chi2ScanData();

        const bool do_position_scan = Config::ENABLE_POSITION_SMEARING;
        const double mu_a0      = res.mu_a;
        const double mu0        = res.mu;   // mu_b
        const double mu_c0      = res.mu_c;
        const double sigma0     = res.sigma;
        const double sigma_pos0 = res.sigma_pos;
        const double p_e_scale0 = res.p_e_scale;

        const int mu_scan_points = max(2, Config::VIS_MU_SCAN_POINTS);
        const int sigma_scan_points = max(2, Config::VIS_SIGMA_SCAN_POINTS);
        const int slice_points = max(2, Config::VIS_SLICE_POINTS);
        const int sigma_pos_scan_points = max(2, Config::VIS_SIGMA_POS_POINTS);

        const double mu_scan_step = (Config::MU_MAX - Config::MU_MIN) / (mu_scan_points - 1);
        const double sigma_scan_step = (Config::SIGMA_MAX - Config::SIGMA_MIN) / (sigma_scan_points - 1);

        // 2D chi2 landscape: scan mu_b (x-axis) vs sigma (y-axis), holding a,c fixed
        for (int i_mu = 0; i_mu < mu_scan_points; ++i_mu) {
            double mu = Config::MU_MIN + i_mu * mu_scan_step;
            for (int i_sig = 0; i_sig < sigma_scan_points; ++i_sig) {
                double sig = Config::SIGMA_MIN + i_sig * sigma_scan_step;
                double chi2 = eval_chi2_selected(mu_a0, mu, mu_c0, sig, sigma_pos0, p_e_scale0,
                                                simEvents,
                                                hdata_mpi0, hdata_mmiss, hdata_mpgg2, rng, Nsmear,
                                                Config::RESOLUTION_A_DEFAULT,
                                                Config::RESOLUTION_B_DEFAULT,
                                                Config::RESOLUTION_C_DEFAULT,
                                                Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);
                res.scan_data.mu_2d.push_back(mu);
                res.scan_data.sigma_2d.push_back(sig);
                res.scan_data.chi2_2d.push_back(chi2);
            }
        }

        // 1D chi2 vs mu_b slice
        const double mu_slice_step = (Config::MU_MAX - Config::MU_MIN) / (slice_points - 1);
        for (int i = 0; i < slice_points; ++i) {
            double mu = Config::MU_MIN + i * mu_slice_step;
            double chi2 = eval_chi2_selected(mu_a0, mu, mu_c0, sigma0, sigma_pos0, p_e_scale0,
                                            simEvents,
                                            hdata_mpi0, hdata_mmiss, hdata_mpgg2, rng, Nsmear,
                                            Config::RESOLUTION_A_DEFAULT,
                                            Config::RESOLUTION_B_DEFAULT,
                                            Config::RESOLUTION_C_DEFAULT,
                                            Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);
            res.scan_data.mu_values.push_back(mu);
            res.scan_data.chi2_vs_mu.push_back(chi2);
        }

        // 1D chi2 vs sigma slice
        const double sigma_slice_step = (Config::SIGMA_MAX - Config::SIGMA_MIN) / (slice_points - 1);
        for (int i = 0; i < slice_points; ++i) {
            double sigma = Config::SIGMA_MIN + i * sigma_slice_step;
            double chi2 = eval_chi2_selected(mu_a0, mu0, mu_c0, sigma, sigma_pos0, p_e_scale0,
                                            simEvents,
                                            hdata_mpi0, hdata_mmiss, hdata_mpgg2, rng, Nsmear,
                                            Config::RESOLUTION_A_DEFAULT,
                                            Config::RESOLUTION_B_DEFAULT,
                                            Config::RESOLUTION_C_DEFAULT,
                                            Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);
            res.scan_data.sigma_values.push_back(sigma);
            res.scan_data.chi2_vs_sigma.push_back(chi2);
        }

        if (do_position_scan) {
            const double sigma_pos_scan_step = (Config::SIGMA_POS_MAX - Config::SIGMA_POS_MIN) / (sigma_pos_scan_points - 1);
            for (int i = 0; i < sigma_pos_scan_points; ++i) {
                double s_pos = Config::SIGMA_POS_MIN + i * sigma_pos_scan_step;
                double chi2 = eval_chi2_selected(mu_a0, mu0, mu_c0, sigma0, s_pos, p_e_scale0,
                                                simEvents,
                                                hdata_mpi0, hdata_mmiss, hdata_mpgg2, rng, Nsmear,
                                                Config::RESOLUTION_A_DEFAULT,
                                                Config::RESOLUTION_B_DEFAULT,
                                                Config::RESOLUTION_C_DEFAULT,
                                                Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);
                res.scan_data.sigma_pos_values.push_back(s_pos);
                res.scan_data.chi2_vs_sigma_pos.push_back(chi2);
            }
        }
    }

    struct MigradRefineResult {
        double mu_a, mu, mu_c;  // mu_eff(E) = mu_a + mu*E + mu_c*ln(E)
        double sigma, sigma_pos, p_e_scale;
        double chi2;
        bool minimized;
        int n_starts = 1;
        bool used_nonzero_ac_seed = false;
        bool hesse_ok = false;
        double max_abs_corr = 0.0;
        string max_corr_pair = "";
    };

    // chi2_fn signature: (mu_a, mu_b, mu_c, sigma, sigma_pos, p_e_scale) -> double
    // In disabled mode: pass fit_mu_a=false, fit_mu_c=false so those are held fixed.
    template <typename Chi2Fn>
    MigradRefineResult run_migrad_refinement(Chi2Fn chi2_fn,
                                            double mu_a_seed,
                                            double mu_seed,    // = mu_b
                                            double mu_c_seed,
                                            double sigma_seed,
                                            double sigma_pos_seed,
                                            double p_e_scale_seed,
                                            bool fit_mu_a,
                                            bool fit_mu,
                                            bool fit_mu_c,
                                            bool fit_sigma,
                                            bool fit_sigma_pos,
                                            bool fit_p_e_scale) {
        MigradRefineResult out{mu_a_seed, mu_seed, mu_c_seed,
                               sigma_seed, sigma_pos_seed, p_e_scale_seed,
                               chi2_fn(mu_a_seed, mu_seed, mu_c_seed,
                                       sigma_seed, sigma_pos_seed, p_e_scale_seed),
                               false};

        if (!Config::USE_MIGRAD_REFINEMENT) return out;
        if (!(fit_mu_a || fit_mu || fit_mu_c || fit_sigma || fit_sigma_pos || fit_p_e_scale)) return out;

        unique_ptr<ROOT::Math::Minimizer> minimizer(
            ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
        if (!minimizer) return out;

        minimizer->SetMaxFunctionCalls(Config::MIGRAD_MAX_FUNCTION_CALLS);
        minimizer->SetMaxIterations(Config::MIGRAD_MAX_ITERATIONS);
        minimizer->SetTolerance(Config::MIGRAD_TOLERANCE);
        minimizer->SetStrategy(Config::MIGRAD_STRATEGY);
        minimizer->SetPrintLevel(0);

        // 6-parameter functor: [mu_a, mu_b, mu_c, sigma, sigma_pos, p_e_scale]
        auto fcn = [&](const double *x) {
            return chi2_fn(x[0], x[1], x[2], x[3], x[4], x[5]);
        };
        ROOT::Math::Functor functor(fcn, 6);
        minimizer->SetFunction(functor);

        if (fit_mu_a) {
            minimizer->SetLimitedVariable(0, "mu_a", mu_a_seed,
                                          Config::MIGRAD_STEP_MU_A,
                                          Config::MU_A_MIN, Config::MU_A_MAX);
        } else {
            minimizer->SetFixedVariable(0, "mu_a", mu_a_seed);
        }

        if (fit_mu) {
            minimizer->SetLimitedVariable(1, "mu_b", mu_seed,
                                          Config::MIGRAD_STEP_MU,
                                          Config::MU_MIN, Config::MU_MAX);
        } else {
            minimizer->SetFixedVariable(1, "mu_b", mu_seed);
        }

        if (fit_mu_c) {
            minimizer->SetLimitedVariable(2, "mu_c", mu_c_seed,
                                          Config::MIGRAD_STEP_MU_C,
                                          Config::MU_C_MIN, Config::MU_C_MAX);
        } else {
            minimizer->SetFixedVariable(2, "mu_c", mu_c_seed);
        }

        double sigma_lower = max(1e-6, Config::SIGMA_MIN);
        if (fit_sigma) {
            minimizer->SetLimitedVariable(3, "sigma", sigma_seed,
                                          Config::MIGRAD_STEP_SIGMA,
                                          sigma_lower, Config::SIGMA_MAX);
        } else {
            minimizer->SetFixedVariable(3, "sigma", sigma_seed);
        }

        if (fit_sigma_pos) {
            minimizer->SetLimitedVariable(4, "sigma_pos", sigma_pos_seed,
                                          Config::MIGRAD_STEP_SIGMA_POS,
                                          Config::SIGMA_POS_MIN, Config::SIGMA_POS_MAX);
        } else {
            minimizer->SetFixedVariable(4, "sigma_pos", sigma_pos_seed);
        }

        if (fit_p_e_scale) {
            minimizer->SetLimitedVariable(5, "p_e_scale", p_e_scale_seed,
                                          Config::MIGRAD_STEP_PE_SCALE,
                                          Config::PE_SCALE_MIN, Config::PE_SCALE_MAX);
        } else {
            minimizer->SetFixedVariable(5, "p_e_scale", p_e_scale_seed);
        }

        bool ok = minimizer->Minimize();
        const double *xbest = minimizer->X();
        if (!xbest) return out;

        double chi2_best = chi2_fn(xbest[0], xbest[1], xbest[2],
                                   xbest[3], xbest[4], xbest[5]);
        // Keep a finite lower-chi2 point even if Minuit reports a non-success
        // status; failures often mean covariance/convergence trouble, not that
        // the returned coordinates are worse than the seed.
        if (std::isfinite(chi2_best) && chi2_best < out.chi2) {
            out.mu_a = xbest[0]; out.mu = xbest[1]; out.mu_c = xbest[2];
            out.sigma = xbest[3]; out.sigma_pos = xbest[4]; out.p_e_scale = xbest[5];
            out.chi2 = chi2_best;
            out.minimized = ok;
            if (ok && Config::ENABLE_HESSE_DIAGNOSTICS) {
                out.hesse_ok = minimizer->Hesse();
                const bool active[6] = {fit_mu_a, fit_mu, fit_mu_c, fit_sigma, fit_sigma_pos, fit_p_e_scale};
                const char* names[6] = {"mu_a", "mu_b", "mu_c", "sigma", "sigma_pos", "p_e_scale"};
                for (int i = 0; i < 6; ++i) {
                    if (!active[i]) continue;
                    for (int j = i + 1; j < 6; ++j) {
                        if (!active[j]) continue;
                        double corr = minimizer->Correlation(i, j);
                        if (std::isfinite(corr) && std::fabs(corr) > out.max_abs_corr) {
                            out.max_abs_corr = std::fabs(corr);
                            out.max_corr_pair = string(names[i]) + ":" + names[j];
                        }
                    }
                }
            }
        }
        return out;
    }

    template <typename Chi2Fn>
    MigradRefineResult run_migrad_ac_multistart_refinement(Chi2Fn chi2_fn,
                                                           double mu_a_seed,
                                                           double mu_seed,
                                                           double mu_c_seed,
                                                           double sigma_seed,
                                                           double sigma_pos_seed,
                                                           double p_e_scale_seed,
                                                           bool fit_mu_a,
                                                           bool fit_mu,
                                                           bool fit_mu_c,
                                                           bool fit_sigma,
                                                           bool fit_sigma_pos,
                                                           bool fit_p_e_scale) {
        MigradRefineResult best = run_migrad_refinement(chi2_fn,
                                                        mu_a_seed, mu_seed, mu_c_seed,
                                                        sigma_seed, sigma_pos_seed, p_e_scale_seed,
                                                        fit_mu_a, fit_mu, fit_mu_c,
                                                        fit_sigma, fit_sigma_pos, fit_p_e_scale);
        int n_starts = 1;

        if (!Config::ENABLE_MU_AC_MULTISTART ||
            !(fit_mu_a || fit_mu_c) ||
            !Config::ENABLE_ENERGY_DEPENDENT_MU) {
            best.n_starts = n_starts;
            return best;
        }

        auto clamp_to = [](double v, double lo, double hi) {
            return std::max(lo, std::min(hi, v));
        };

        const double da = Config::MU_A_MULTISTART_SPAN;
        const double dc = Config::MU_C_MULTISTART_SPAN;
        const vector<pair<double, double>> ac_offsets = {
            { da,  0.0}, {-da,  0.0},
            {0.0,   dc}, {0.0,  -dc},
            { da,   dc}, {-da, -dc},
            { da,  -dc}, {-da,  dc}
        };

        for (const auto &off : ac_offsets) {
            const double a_start = fit_mu_a
                ? clamp_to(mu_a_seed + off.first, Config::MU_A_MIN, Config::MU_A_MAX)
                : mu_a_seed;
            const double c_start = fit_mu_c
                ? clamp_to(mu_c_seed + off.second, Config::MU_C_MIN, Config::MU_C_MAX)
                : mu_c_seed;

            if (std::fabs(a_start - mu_a_seed) < 1e-15 &&
                std::fabs(c_start - mu_c_seed) < 1e-15) {
                continue;
            }

            MigradRefineResult cand = run_migrad_refinement(chi2_fn,
                                                            a_start, mu_seed, c_start,
                                                            sigma_seed, sigma_pos_seed, p_e_scale_seed,
                                                            fit_mu_a, fit_mu, fit_mu_c,
                                                            fit_sigma, fit_sigma_pos, fit_p_e_scale);
            ++n_starts;
            if (cand.chi2 < best.chi2) {
                best = cand;
                best.used_nonzero_ac_seed = true;
            }
        }

        best.n_starts = n_starts;
        return best;
    }

    inline double lerp_bound(double lo, double hi, double u) {
        return lo + (hi - lo) * std::min(1.0, std::max(0.0, u));
    }

    inline double sobol_unit_safe(double u) {
        const double eps = 1e-12;
        if (!std::isfinite(u)) return 0.5;
        return std::min(1.0 - eps, std::max(eps, u));
    }

    struct FitFlags {
        bool mu_a = false;
        bool mu = true;
        bool mu_c = false;
        bool sigma = true;
        bool sigma_pos = false;
        bool p_e_scale = false;
    };

    struct ParamPoint {
        double mu_a = Config::MU_ENERGY_A_INIT;
        double mu = Config::MU_ENERGY_B_INIT;
        double mu_c = Config::MU_ENERGY_C_INIT;
        double sigma = 0.05;
        double sigma_pos = 0.0;
        double p_e_scale = 1.0;
    };

    inline void print_progress_bar(const string &label, int done, int total) {
        if (label.empty() || total <= 0) return;
        const int width = 36;
        done = max(0, min(done, total));
        double frac = (double)done / (double)total;
        int filled = (int)floor(frac * width);
        cout << "\r" << label << " [";
        for (int i = 0; i < width; ++i) cout << (i < filled ? '=' : ' ');
        cout << "] " << done << "/" << total << flush;
        if (done >= total) cout << endl;
    }

    template <typename Chi2Fn>
    MigradRefineResult run_global_multistart_refinement(Chi2Fn chi2_fn,
                                                        const ParamPoint &center,
                                                        const FitFlags &flags,
                                                        vector<SeedDiagnostic> *seed_debug,
                                                        const string &progress_label = "") {
        int ndim = 0;
        if (flags.mu_a) ++ndim;
        if (flags.mu) ++ndim;
        if (flags.mu_c) ++ndim;
        if (flags.sigma) ++ndim;
        if (flags.sigma_pos) ++ndim;
        if (flags.p_e_scale) ++ndim;
        ndim = std::max(1, ndim);

        const int n_seeds = std::max(1, Config::GLOBAL_MULTISTART_SEEDS);
        vector<array<double, 6>> sobol_points(n_seeds + 1);
        for (auto &pt : sobol_points) pt.fill(0.5);

        ROOT::Math::QuasiRandomSobol sobol((unsigned int)ndim);
        for (int iseed = 1; iseed <= n_seeds; ++iseed) {
            double u[6] = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
            if (!sobol.Next(u)) {
                std::cerr << "[WARN] Sobol generator failed at seed " << iseed
                          << "; using centered fallback for this seed." << std::endl;
            }
            for (int idim = 0; idim < ndim && idim < 6; ++idim) {
                sobol_points[iseed][idim] = sobol_unit_safe(u[idim]);
            }
        }

        vector<SeedDiagnostic> candidates;
        candidates.reserve(Config::GLOBAL_MULTISTART_SEEDS + 1);

        auto make_point = [&](int seed_index) {
            ParamPoint p = center;
            int dim = 0;
            if (flags.mu_a) {
                p.mu_a = (seed_index == 0) ? center.mu_a :
                    lerp_bound(Config::MU_A_MIN, Config::MU_A_MAX, sobol_points[seed_index][dim]);
                ++dim;
            }
            if (flags.mu) {
                p.mu = (seed_index == 0) ? center.mu :
                    lerp_bound(Config::MU_MIN, Config::MU_MAX, sobol_points[seed_index][dim]);
                ++dim;
            }
            if (flags.mu_c) {
                p.mu_c = (seed_index == 0) ? center.mu_c :
                    lerp_bound(Config::MU_C_MIN, Config::MU_C_MAX, sobol_points[seed_index][dim]);
                ++dim;
            }
            if (flags.sigma) {
                p.sigma = (seed_index == 0) ? center.sigma :
                    lerp_bound(Config::SIGMA_MIN, Config::SIGMA_MAX, sobol_points[seed_index][dim]);
                ++dim;
            }
            if (flags.sigma_pos) {
                p.sigma_pos = (seed_index == 0) ? center.sigma_pos :
                    lerp_bound(Config::SIGMA_POS_MIN, Config::SIGMA_POS_MAX, sobol_points[seed_index][dim]);
                ++dim;
            }
            if (flags.p_e_scale) {
                p.p_e_scale = (seed_index == 0) ? center.p_e_scale :
                    lerp_bound(Config::PE_SCALE_MIN, Config::PE_SCALE_MAX, sobol_points[seed_index][dim]);
                ++dim;
            }
            p.mu_a = std::min(std::max(p.mu_a, Config::MU_A_MIN), Config::MU_A_MAX);
            p.mu = std::min(std::max(p.mu, Config::MU_MIN), Config::MU_MAX);
            p.mu_c = std::min(std::max(p.mu_c, Config::MU_C_MIN), Config::MU_C_MAX);
            p.sigma = std::min(std::max(p.sigma, Config::SIGMA_MIN), Config::SIGMA_MAX);
            p.sigma_pos = std::min(std::max(p.sigma_pos, Config::SIGMA_POS_MIN), Config::SIGMA_POS_MAX);
            p.p_e_scale = std::min(std::max(p.p_e_scale, Config::PE_SCALE_MIN), Config::PE_SCALE_MAX);
            return p;
        };

        for (int iseed = 0; iseed <= n_seeds; ++iseed) {
            ParamPoint p = make_point(iseed);
            SeedDiagnostic d;
            d.seed_index = iseed;
            d.mu_a = p.mu_a;
            d.mu = p.mu;
            d.mu_c = p.mu_c;
            d.sigma = p.sigma;
            d.sigma_pos = p.sigma_pos;
            d.p_e_scale = p.p_e_scale;
            d.seed_chi2 = chi2_fn(p.mu_a, p.mu, p.mu_c, p.sigma, p.sigma_pos, p.p_e_scale);
            d.migrad_chi2 = d.seed_chi2;
            candidates.push_back(d);
            if (!progress_label.empty()) {
                print_progress_bar(progress_label + " Sobol seed scan", iseed + 1, n_seeds + 1);
            }
        }

        std::sort(candidates.begin(), candidates.end(),
                  [](const SeedDiagnostic &a, const SeedDiagnostic &b) {
                      return a.seed_chi2 < b.seed_chi2;
                  });

        const int keep = std::min((int)candidates.size(),
                                  std::max(1, Config::GLOBAL_MULTISTART_KEEP_BEST));
        candidates.resize(keep);

        MigradRefineResult best{candidates[0].mu_a, candidates[0].mu, candidates[0].mu_c,
                                candidates[0].sigma, candidates[0].sigma_pos,
                                candidates[0].p_e_scale, candidates[0].seed_chi2, false};
        best.n_starts = keep;

        for (int i = 0; i < keep; ++i) {
            if (!progress_label.empty()) {
                print_progress_bar(progress_label + " MIGRAD retained seeds", i, keep);
            }
            SeedDiagnostic &d = candidates[i];
            d.rank = i + 1;
            MigradRefineResult mr = run_migrad_refinement(chi2_fn,
                                                          d.mu_a, d.mu, d.mu_c,
                                                          d.sigma, d.sigma_pos, d.p_e_scale,
                                                          flags.mu_a, flags.mu, flags.mu_c,
                                                          flags.sigma, flags.sigma_pos, flags.p_e_scale);
            d.migrad_chi2 = mr.chi2;
            d.mu_a = mr.mu_a;
            d.mu = mr.mu;
            d.mu_c = mr.mu_c;
            d.sigma = mr.sigma;
            d.sigma_pos = mr.sigma_pos;
            d.p_e_scale = mr.p_e_scale;
            d.minimized = mr.minimized;
            d.hesse_ok = mr.hesse_ok;
            d.max_abs_corr = mr.max_abs_corr;
            d.max_corr_pair = mr.max_corr_pair;
            if (mr.chi2 < best.chi2) best = mr;
            if (!progress_label.empty()) {
                print_progress_bar(progress_label + " MIGRAD retained seeds", i + 1, keep);
            }
        }

        if (seed_debug) *seed_debug = candidates;
        return best;
    }

    template <typename Chi2Fn>
    vector<ProfileDiagnostic> build_profile_diagnostics(Chi2Fn chi2_fn,
                                                        const MigradRefineResult &best,
                                                        const FitFlags &flags) {
        vector<ProfileDiagnostic> out;
        if (!Config::ENABLE_PROFILE_SCANS ||
            !best.hesse_ok ||
            best.max_abs_corr < Config::PROFILE_CORR_THRESHOLD) {
            return out;
        }

        auto add_profile = [&](const string &name, int ipar, double center,
                               double lo, double hi, const FitFlags &base_flags) {
            int npoints = std::max(3, Config::PROFILE_SCAN_POINTS);
            double span = Config::PROFILE_RANGE_FRACTION * (hi - lo);
            double plo = std::max(lo, center - span);
            double phi = std::min(hi, center + span);
            for (int i = 0; i < npoints; ++i) {
                double v = (npoints == 1) ? center : plo + (phi - plo) * i / (npoints - 1);
                double a = best.mu_a, b = best.mu, c = best.mu_c;
                double sig = best.sigma, spos = best.sigma_pos, pe = best.p_e_scale;
                FitFlags pf = base_flags;
                if (ipar == 0) { a = v; pf.mu_a = false; }
                if (ipar == 1) { b = v; pf.mu = false; }
                if (ipar == 2) { c = v; pf.mu_c = false; }
                if (ipar == 3) { sig = v; pf.sigma = false; }
                if (ipar == 4) { spos = v; pf.sigma_pos = false; }
                if (ipar == 5) { pe = v; pf.p_e_scale = false; }
                MigradRefineResult pr = run_migrad_refinement(chi2_fn, a, b, c, sig, spos, pe,
                                                              pf.mu_a, pf.mu, pf.mu_c,
                                                              pf.sigma, pf.sigma_pos, pf.p_e_scale);
                ProfileDiagnostic d;
                d.parameter = name;
                d.fixed_value = v;
                d.chi2 = pr.chi2;
                d.minimized = pr.minimized;
                out.push_back(d);
            }
        };

        auto pair_has = [&](const string &name) {
            size_t sep = best.max_corr_pair.find(':');
            if (sep == string::npos) return best.max_corr_pair == name;
            return best.max_corr_pair.substr(0, sep) == name ||
                   best.max_corr_pair.substr(sep + 1) == name;
        };

        int nprofiled = 0;
        if (flags.mu_a && nprofiled < Config::PROFILE_MAX_PARAMETERS &&
            pair_has("mu_a")) {
            add_profile("mu_a", 0, best.mu_a, Config::MU_A_MIN, Config::MU_A_MAX, flags);
            ++nprofiled;
        }
        if (flags.mu && nprofiled < Config::PROFILE_MAX_PARAMETERS &&
            pair_has("mu_b")) {
            add_profile("mu_b", 1, best.mu, Config::MU_MIN, Config::MU_MAX, flags);
            ++nprofiled;
        }
        if (flags.mu_c && nprofiled < Config::PROFILE_MAX_PARAMETERS &&
            pair_has("mu_c")) {
            add_profile("mu_c", 2, best.mu_c, Config::MU_C_MIN, Config::MU_C_MAX, flags);
            ++nprofiled;
        }
        if (flags.sigma && nprofiled < Config::PROFILE_MAX_PARAMETERS &&
            pair_has("sigma")) {
            add_profile("sigma", 3, best.sigma, Config::SIGMA_MIN, Config::SIGMA_MAX, flags);
            ++nprofiled;
        }
        if (flags.sigma_pos && nprofiled < Config::PROFILE_MAX_PARAMETERS &&
            pair_has("sigma_pos")) {
            add_profile("sigma_pos", 4, best.sigma_pos, Config::SIGMA_POS_MIN, Config::SIGMA_POS_MAX, flags);
            ++nprofiled;
        }
        if (flags.p_e_scale && nprofiled < Config::PROFILE_MAX_PARAMETERS &&
            pair_has("p_e_scale")) {
            add_profile("p_e_scale", 5, best.p_e_scale, Config::PE_SCALE_MIN, Config::PE_SCALE_MAX, flags);
        }
        return out;
    }

    // rotate vector v around axis a_unit by angle theta (radians)
    TVector3 rotateAroundAxis(const TVector3 &v, const TVector3 &a_unit, double theta) {
        TVector3 term1 = v * cos(theta);
        TVector3 term2 = a_unit.Cross(v) * sin(theta);
        TVector3 term3 = a_unit * (a_unit.Dot(v)) * (1.0 - cos(theta));
        return term1 + term2 + term3;
    }

    // Plot chi2 scans for a section
    // Helper: draw a comparison pad (data vs unsmeared vs smeared) with legend
    static void drawComparisonPad(const TH1D& hdata, const TH1D& hsim, const TH1D& hunsm,
                                const char* xtitle, const char* label) {
        TH1D* hd = (TH1D*)hdata.Clone(Form("_diag_d_%s", label));
        TH1D* hs = (TH1D*)hsim.Clone(Form("_diag_s_%s", label));
        TH1D* hu = (TH1D*)hunsm.Clone(Form("_diag_u_%s", label));
        hd->SetLineColor(kBlack); hd->SetMarkerStyle(20); hd->SetMarkerSize(0.5); hd->SetLineWidth(2);
        hu->SetLineColor(kBlue);  hu->SetLineWidth(2); hu->SetLineStyle(2);
        hs->SetLineColor(kRed);   hs->SetLineWidth(2);
        hd->GetXaxis()->SetTitle(xtitle);
        hd->GetYaxis()->SetTitle("Counts");
        hd->GetYaxis()->SetTitleOffset(1.3);
        hd->GetXaxis()->SetTitleSize(0.05); hd->GetYaxis()->SetTitleSize(0.05);
        hd->GetXaxis()->SetLabelSize(0.04); hd->GetYaxis()->SetLabelSize(0.04);
        double mx = max({hd->GetMaximum(), hs->GetMaximum(), hu->GetMaximum()});
        hd->SetMaximum(mx * 1.25);
        hd->SetMinimum(0);
        hd->Draw("E");
        hu->Draw("HIST SAME");
        hs->Draw("HIST SAME");
        TLegend *lg = new TLegend(0.48, 0.67, 0.88, 0.87);
        lg->SetBorderSize(1); lg->SetFillColor(0); lg->SetTextSize(0.033);
        lg->AddEntry(hd, "Data", "lep");
        lg->AddEntry(hu, "Unsmeared Sim", "l");
        lg->AddEntry(hs, "Smeared Sim", "l");
        lg->Draw();
    }

    // Helper: draw 1D chi2 scan with minimum marker and Δχ²=1 reference line
    static void drawChi2Scan1D(const vector<double>& xvals, const vector<double>& chi2vals,
                                double best_x, double best_chi2,
                                const char* xtitle, const char* title) {
        if (xvals.empty()) return;
        TGraph *g = new TGraph(xvals.size(), &xvals[0], &chi2vals[0]);
        g->SetMarkerStyle(20); g->SetMarkerSize(0.7);
        g->SetTitle(Form("%s;%s;#chi^{2}", title, xtitle));
        g->GetXaxis()->SetTitleOffset(1.2);
        g->GetYaxis()->SetTitleOffset(1.4);
        g->GetXaxis()->SetTitleSize(0.05); g->GetYaxis()->SetTitleSize(0.05);
        g->Draw("APL");
        // Δχ² = 1 line (1σ parameter error)
        double xlo = *min_element(xvals.begin(), xvals.end());
        double xhi = *max_element(xvals.begin(), xvals.end());
        TLine *l1 = new TLine(xlo, best_chi2 + 1.0, xhi, best_chi2 + 1.0);
        l1->SetLineColor(kOrange+1); l1->SetLineWidth(2); l1->SetLineStyle(2); l1->Draw();
        TLine *l4 = new TLine(xlo, best_chi2 + 4.0, xhi, best_chi2 + 4.0);
        l4->SetLineColor(kMagenta);  l4->SetLineWidth(1); l4->SetLineStyle(3); l4->Draw();
        // Mark minimum
        TMarker *mk = new TMarker(best_x, best_chi2, 29);
        mk->SetMarkerColor(kRed); mk->SetMarkerSize(2.0); mk->Draw();
        TLatex *tx = new TLatex(); tx->SetNDC(); tx->SetTextSize(0.038);
        tx->DrawLatex(0.15, 0.87, Form("Best = %.5f", best_x));
        tx->DrawLatex(0.15, 0.82, Form("#chi^{2}_{min} = %.2f", best_chi2));
        // Legend for reference lines
        TLegend *ll = new TLegend(0.50, 0.72, 0.88, 0.82);
        ll->SetBorderSize(1); ll->SetFillColor(0); ll->SetTextSize(0.030);
        ll->AddEntry(l1, "#Delta#chi^{2}=1  (1#sigma)", "l");
        ll->AddEntry(l4, "#Delta#chi^{2}=4  (2#sigma)", "l");
        ll->Draw();
    }

    static void drawPullPad(const TH1D& hdata, const TH1D& hsim,
                            const char* xtitle, const char* title, const char* label) {
        int nb = hdata.GetNbinsX();
        TH1D *hpull = new TH1D(Form("_diag_pull_%s", label),
                               Form("%s;%s;Pull", title, xtitle),
                               nb, hdata.GetXaxis()->GetXmin(), hdata.GetXaxis()->GetXmax());
        hpull->GetXaxis()->SetTitle(xtitle);
        hpull->GetYaxis()->SetTitle("(Data-Sm.)/#sigma");
        for (int ib = 1; ib <= nb; ++ib) {
            double d = hdata.GetBinContent(ib);
            double s = hsim.GetBinContent(ib);
            double ed = hdata.GetBinError(ib);
            double es = hsim.GetBinError(ib);
            double denom = sqrt(ed*ed + es*es);
            if (denom > 0.0) hpull->SetBinContent(ib, (d - s) / denom);
        }
        hpull->SetLineColor(kBlack);
        hpull->SetFillColor(kCyan-9);
        double pmax = max(fabs(hpull->GetMinimum()), fabs(hpull->GetMaximum()));
        pmax = max(pmax * 1.20, 3.0);
        hpull->SetMaximum(pmax);
        hpull->SetMinimum(-pmax);
        hpull->Draw("HIST");
        TLine *lz = new TLine(hdata.GetXaxis()->GetXmin(), 0.0, hdata.GetXaxis()->GetXmax(), 0.0);
        lz->SetLineColor(kRed); lz->SetLineWidth(2); lz->Draw();
        TLine *lp2 = new TLine(hdata.GetXaxis()->GetXmin(), 2.0, hdata.GetXaxis()->GetXmax(), 2.0);
        TLine *lm2 = new TLine(hdata.GetXaxis()->GetXmin(), -2.0, hdata.GetXaxis()->GetXmax(), -2.0);
        lp2->SetLineColor(kOrange+1); lp2->SetLineStyle(2); lp2->Draw();
        lm2->SetLineColor(kOrange+1); lm2->SetLineStyle(2); lm2->Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.045);
        tx.DrawLatex(0.15, 0.88, title);
    }

    static void drawFractionalResidualPad(const TH1D& hdata, const TH1D& hsim,
                                          const char* xtitle, const char* title,
                                          const char* label) {
        int nb = hdata.GetNbinsX();
        TH1D *hres = new TH1D(Form("_diag_fracres_%s", label),
                              Form("%s;%s;(Sm.-Data)/Data", title, xtitle),
                              nb, hdata.GetXaxis()->GetXmin(), hdata.GetXaxis()->GetXmax());
        for (int ib = 1; ib <= nb; ++ib) {
            double d = hdata.GetBinContent(ib);
            double s = hsim.GetBinContent(ib);
            if (d > 0.0) hres->SetBinContent(ib, (s - d) / d);
        }
        hres->SetLineColor(kBlack);
        hres->SetFillColor(kOrange-9);
        double rmax = max(fabs(hres->GetMinimum()), fabs(hres->GetMaximum()));
        rmax = max(rmax * 1.20, 0.25);
        hres->SetMaximum(rmax);
        hres->SetMinimum(-rmax);
        hres->Draw("HIST");
        TLine *lz = new TLine(hdata.GetXaxis()->GetXmin(), 0.0, hdata.GetXaxis()->GetXmax(), 0.0);
        lz->SetLineColor(kRed); lz->SetLineWidth(2); lz->Draw();
        TLatex tx; tx.SetNDC(); tx.SetTextSize(0.045);
        tx.DrawLatex(0.15, 0.88, title);
    }

    static void drawNoDataText(const char *msg) {
        TLatex tx;
        tx.SetNDC();
        tx.SetTextSize(0.045);
        tx.DrawLatex(0.15, 0.55, msg);
    }

    static double normalizedPosition(double v, double lo, double hi) {
        if (!(hi > lo) || !std::isfinite(v)) return 0.0;
        return min(1.0, max(0.0, (v - lo) / (hi - lo)));
    }

    void plot_section_diagnostics(const Section& sec, const FitResult& fres, TCanvas* c,
                        const string &pdf_file,
                        int &page_count,
                        TDirectory *canvas_dir,
                        const TH1D& hdata,          const TH1D& hsim_final,       const TH1D& hsim_unsmeared,
                        const TH1D& hdata_mmiss,    const TH1D& hsim_mmiss,       const TH1D& hunsmeared_mmiss,
                        const TH1D& hdata_mpgg2,    const TH1D& hsim_mpgg2,       const TH1D& hunsmeared_mpgg2,
                        int n_data, int n_data_selected, int n_sim,
                        int nsmear,
                        double global_p_e_scale,
                        bool use_global_p_e_scale) {
        HistObjectiveMetrics mgg_metrics = computeChi2MetricsFromHist(hsim_final, hdata);
        HistObjectiveMetrics mmiss_metrics = computeChi2MetricsFromHist(hsim_mmiss, hdata_mmiss);
        HistObjectiveMetrics mpgg2_metrics = computeChi2MetricsFromHist(hsim_mpgg2, hdata_mpgg2);
        double selected_cost = Config::W_MPI0 * mgg_metrics.chi2
                             + Config::W_MMISS * mmiss_metrics.chi2
                             + Config::W_MPGG2 * mpgg2_metrics.chi2;

        c->Clear();
        c->Divide(3, 3);
        c->cd(1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        drawComparisonPad(hdata, hsim_final, hsim_unsmeared,
                        "M_{#gamma#gamma} [GeV/c^{2}]", Form("%s_mgg_cmp", sec.name().c_str()));
        TLatex tx_title; tx_title.SetNDC(); tx_title.SetTextSize(0.045);
        tx_title.DrawLatex(0.15, 0.93, Form("%s  x=[%.1f, %.1f] y=[%.1f, %.1f] cm",
                                            sec.name().c_str(), sec.xlo, sec.xhi, sec.ylo, sec.yhi));

        c->cd(2); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        drawComparisonPad(hdata_mmiss, hsim_mmiss, hunsmeared_mmiss,
                        "M_{miss} [GeV/c^{2}]", Form("%s_mmiss_cmp", sec.name().c_str()));

        c->cd(3); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        drawComparisonPad(hdata_mpgg2, hsim_mpgg2, hunsmeared_mpgg2,
                        "(p_{target}+#gamma#gamma)^{2} [GeV^{2}]", Form("%s_mpgg2_cmp", sec.name().c_str()));

        c->cd(4); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        drawPullPad(hdata, hsim_final, "M_{#gamma#gamma} [GeV/c^{2}]",
                    "M_{#gamma#gamma} pull", Form("%s_mgg_pull", sec.name().c_str()));
        c->cd(5); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        drawPullPad(hdata_mmiss, hsim_mmiss, "M_{miss} [GeV/c^{2}]",
                    "M_{miss} pull", Form("%s_mmiss_pull", sec.name().c_str()));
        c->cd(6); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        drawPullPad(hdata_mpgg2, hsim_mpgg2, "(p_{target}+#gamma#gamma)^{2} [GeV^{2}]",
                    "M_{p#gamma#gamma}^{2} pull", Form("%s_mpgg2_pull", sec.name().c_str()));

        c->cd(7); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        drawFractionalResidualPad(hdata, hsim_final, "M_{#gamma#gamma} [GeV/c^{2}]",
                                  "M_{#gamma#gamma} fractional residual", Form("%s_mgg_frac", sec.name().c_str()));
        c->cd(8); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        drawFractionalResidualPad(hdata_mmiss, hsim_mmiss, "M_{miss} [GeV/c^{2}]",
                                  "M_{miss} fractional residual", Form("%s_mmiss_frac", sec.name().c_str()));
        c->cd(9); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        drawFractionalResidualPad(hdata_mpgg2, hsim_mpgg2, "(p_{target}+#gamma#gamma)^{2} [GeV^{2}]",
                                  "M_{p#gamma#gamma}^{2} fractional residual", Form("%s_mpgg2_frac", sec.name().c_str()));
        c->Print(pdf_file.c_str());
        writeCanvasToDir(canvas_dir, c,
                         Form("c_%s_observable_comparison", sec.name().c_str()),
                         Form("%s observable comparison", sec.name().c_str()));
        ++page_count;

        c->Clear();
        c->Divide(3, 2);
        c->cd(1);
        gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.13); gPad->SetGridx(); gPad->SetGridy();
        {
            const int NE = 120;
            const double E_lo = Config::RESPONSE_CURVE_E_MIN_GEV;
            const double E_hi = Config::RESPONSE_CURVE_E_MAX_GEV;
            TGraph *g_mean = new TGraph(NE);
            TGraph *g_identity = new TGraph(NE);
            for (int i = 0; i < NE; ++i) {
                double E = E_lo + (E_hi - E_lo) * i / (NE - 1);
                g_mean->SetPoint(i, E, fres.muEff(E));
                g_identity->SetPoint(i, E, E);
            }
            g_mean->SetTitle("#mu_{eff}(E): reconstructed mean energy;E_{true} [GeV];E_{mean} [GeV]");
            g_mean->SetLineColor(kRed); g_mean->SetLineWidth(3);
            g_identity->SetLineColor(kBlack); g_identity->SetLineStyle(2); g_identity->SetLineWidth(2);
            g_mean->Draw("AL");
            g_identity->Draw("L SAME");
            TLegend *lg = new TLegend(0.16, 0.70, 0.55, 0.88);
            lg->SetBorderSize(1); lg->SetFillColor(0); lg->SetTextSize(0.033);
            lg->AddEntry(g_mean, "#mu_{eff}(E)", "l");
            lg->AddEntry(g_identity, "identity", "l");
            lg->Draw();
        }

        c->cd(2);
        gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.13); gPad->SetGridx(); gPad->SetGridy();
        {
            const int NE = 120;
            const double E_lo = Config::RESPONSE_CURVE_E_MIN_GEV;
            const double E_hi = Config::RESPONSE_CURVE_E_MAX_GEV;
            TGraph *g_ratio = new TGraph(NE);
            for (int i = 0; i < NE; ++i) {
                double E = E_lo + (E_hi - E_lo) * i / (NE - 1);
                g_ratio->SetPoint(i, E, (E > 0.0) ? fres.muEff(E) / E : 0.0);
            }
            g_ratio->SetTitle("Energy response ratio;E_{true} [GeV];#mu_{eff}(E)/E");
            g_ratio->SetLineColor(kBlue+1); g_ratio->SetLineWidth(3);
            g_ratio->Draw("AL");
            TLine *l1 = new TLine(E_lo, 1.0, E_hi, 1.0);
            l1->SetLineColor(kBlack); l1->SetLineStyle(2); l1->Draw();
        }

        c->cd(3);
        gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.13); gPad->SetGridx(); gPad->SetGridy();
        {
            const int NE = 120;
            const double E_lo = Config::RESPONSE_CURVE_E_MIN_GEV;
            const double E_hi = Config::RESPONSE_CURVE_E_MAX_GEV;
            TGraph *g_rel = new TGraph(NE);
            TGraph *g_pos = new TGraph(NE);
            for (int i = 0; i < NE; ++i) {
                double E = E_lo + (E_hi - E_lo) * i / (NE - 1);
                double Emean = max(Config::NONPOSITIVE_CLAMP, fres.muEff(E));
                double sigE = computeEnergyResolution(Emean, fres.sigma, fres.res_A, fres.res_B, fres.res_C);
                double rel = (Emean > 0.0) ? sigE / Emean : 0.0;
                double spos = fres.sigma_pos;
                if (Config::ENABLE_ENERGY_DEPENDENT_SIGMA_POS) {
                    double E_for_pos = max(Emean, Config::SIGMA_POS_ENERGY_MIN_GEV);
                    spos *= sqrt(Config::SIGMA_POS_ENERGY_E0_GEV / E_for_pos);
                }
                g_rel->SetPoint(i, E, rel);
                g_pos->SetPoint(i, E, spos);
            }
            TMultiGraph *mg = new TMultiGraph("mg_section_res", "Resolution model;E_{true} [GeV];value");
            g_rel->SetLineColor(kRed); g_rel->SetLineWidth(3);
            g_pos->SetLineColor(kGreen+2); g_pos->SetLineWidth(3);
            mg->Add(g_rel, "L");
            mg->Add(g_pos, "L");
            mg->Draw("A");
            TLegend *lg = new TLegend(0.15, 0.70, 0.65, 0.88);
            lg->SetBorderSize(1); lg->SetFillColor(0); lg->SetTextSize(0.030);
            lg->AddEntry(g_rel, "#sigma_{E}/E_{mean}", "l");
            lg->AddEntry(g_pos, "#sigma_{pos,eff} [cm]", "l");
            lg->Draw();
        }

        c->cd(4);
        gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.15);
        {
            TH1D *hcost = new TH1D(Form("hcost_%s", sec.name().c_str()),
                                   "Objective contribution after area normalization;observable;weighted #chi^{2}",
                                   3, 0.5, 3.5);
            hcost->GetXaxis()->SetBinLabel(1, "Mgg");
            hcost->GetXaxis()->SetBinLabel(2, "Mmiss");
            hcost->GetXaxis()->SetBinLabel(3, "Mpgg2");
            hcost->SetBinContent(1, Config::W_MPI0 * mgg_metrics.chi2);
            hcost->SetBinContent(2, Config::W_MMISS * mmiss_metrics.chi2);
            hcost->SetBinContent(3, Config::W_MPGG2 * mpgg2_metrics.chi2);
            hcost->SetFillColor(kAzure-9);
            hcost->SetStats(0);
            hcost->Draw("BAR2 TEXT");
            TLatex tx; tx.SetNDC(); tx.SetTextSize(0.040);
            tx.DrawLatex(0.14, 0.87, Form("Selected cost = %.3g", selected_cost));
            tx.DrawLatex(0.14, 0.81, Form("raw: Mgg %.3g, Mmiss %.3g, Mpgg2 %.3g",
                                          mgg_metrics.chi2, mmiss_metrics.chi2, mpgg2_metrics.chi2));
        }

        c->cd(5);
        gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.20);
        {
            vector<string> names;
            vector<double> vals;
            names.push_back("a"); vals.push_back(normalizedPosition(fres.mu_a, Config::MU_A_MIN, Config::MU_A_MAX));
            names.push_back("b"); vals.push_back(normalizedPosition(fres.mu, Config::MU_MIN, Config::MU_MAX));
            names.push_back("c"); vals.push_back(normalizedPosition(fres.mu_c, Config::MU_C_MIN, Config::MU_C_MAX));
            names.push_back("sig"); vals.push_back(normalizedPosition(fres.sigma, Config::SIGMA_MIN, Config::SIGMA_MAX));
            names.push_back("pos"); vals.push_back(normalizedPosition(fres.sigma_pos, Config::SIGMA_POS_MIN, Config::SIGMA_POS_MAX));
            if (Config::ENABLE_ELECTRON_MOMENTUM_SCALING && Config::fit_objective_uses_mmiss()) {
                names.push_back("pe"); vals.push_back(normalizedPosition(fres.p_e_scale, Config::PE_SCALE_MIN, Config::PE_SCALE_MAX));
            }
            TH1D *hbound = new TH1D(Form("hbound_%s", sec.name().c_str()),
                                    "Parameter location inside fit bounds;parameter;(value-min)/(max-min)",
                                    names.size(), 0.5, names.size() + 0.5);
            for (int i = 0; i < (int)names.size(); ++i) {
                hbound->GetXaxis()->SetBinLabel(i + 1, names[i].c_str());
                hbound->SetBinContent(i + 1, vals[i]);
            }
            hbound->SetMinimum(0.0);
            hbound->SetMaximum(1.0);
            hbound->SetFillColor(kOrange-9);
            hbound->SetStats(0);
            hbound->Draw("BAR2 TEXT");
            TLine *llo = new TLine(0.5, 0.05, names.size() + 0.5, 0.05);
            TLine *lhi = new TLine(0.5, 0.95, names.size() + 0.5, 0.95);
            llo->SetLineColor(kRed); llo->SetLineStyle(2); llo->Draw();
            lhi->SetLineColor(kRed); lhi->SetLineStyle(2); lhi->Draw();
        }

        c->cd(6);
        gPad->SetLeftMargin(0.05); gPad->SetRightMargin(0.05);
        {
            double chi2_ndf_val = fres.chi2_per_ndf();
            bool good = (chi2_ndf_val <= Config::MAX_CHI2_PER_NDF);
            TPaveText *info = new TPaveText(0.03, 0.03, 0.97, 0.97, "brNDC");
            info->SetBorderSize(1); info->SetFillColor(good ? kGreen-9 : kRed-9);
            info->SetTextAlign(12); info->SetTextSize(0.030);
            info->AddText(Form("Section %s", sec.name().c_str()));
            info->AddText(Form("mu_eff(E)=a+bE+c ln(E): a=%.5f b=%.5f c=%.5f", fres.mu_a, fres.mu, fres.mu_c));
            info->AddText(Form("sigma=%.5f, sigma_pos=%.5f cm", fres.sigma, fres.sigma_pos));
            if (Config::ENABLE_ELECTRON_MOMENTUM_SCALING && Config::fit_objective_uses_mmiss())
                info->AddText(Form("p_e_scale=%.6f", fres.p_e_scale));
            info->AddText(Form("chi2=%.3g, ndf=%d, chi2/ndf=%.4f", fres.chi2, fres.ndf, chi2_ndf_val));
            info->AddText(Form("HESSE=%s, max |corr|=%.3f %s", fres.hesse_ok ? "ok" : "not ok",
                               fres.max_abs_corr, fres.max_corr_pair.c_str()));
            info->AddText(Form("seeds=%d, MIGRAD trials=%d, Nsmear=%d", fres.n_seed_trials, fres.n_migrad_trials, nsmear));
            info->AddText(Form("Ndata=%d, Nselected=%d, Nsim=%d", n_data, n_data_selected, n_sim));
            info->AddText(Form("Objective: %s, weights %.2f %.2f %.2f",
                               Config::histogram_mode_label(),
                               Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2));
            if (use_global_p_e_scale) info->AddText(Form("Global p_e_scale=%.6f", global_p_e_scale));
            info->Draw();
        }
        c->Print(pdf_file.c_str());
        writeCanvasToDir(canvas_dir, c,
                         Form("c_%s_response_summary", sec.name().c_str()),
                         Form("%s response and fit summary", sec.name().c_str()));
        ++page_count;

        c->Clear();
        c->Divide(3, 2);
        c->cd(1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        drawChi2Scan1D(fres.scan_data.mu_values, fres.scan_data.chi2_vs_mu,
                       fres.mu, fres.chi2, "#mu_{b}", "#chi^{2} vs #mu_{b}");
        c->cd(2); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        drawChi2Scan1D(fres.scan_data.sigma_values, fres.scan_data.chi2_vs_sigma,
                       fres.sigma, fres.chi2, "#sigma", "#chi^{2} vs #sigma");
        c->cd(3); gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.15); gPad->SetBottomMargin(0.12);
        if (!fres.scan_data.mu_2d.empty()) {
            TGraph2D *g2d = new TGraph2D(fres.scan_data.mu_2d.size(),
                                        const_cast<double*>(&fres.scan_data.mu_2d[0]),
                                        const_cast<double*>(&fres.scan_data.sigma_2d[0]),
                                        const_cast<double*>(&fres.scan_data.chi2_2d[0]));
            g2d->SetTitle("#chi^{2} landscape;#mu_{b};#sigma;#chi^{2}");
            g2d->Draw("COLZ");
            TMarker *mk2d = new TMarker(fres.mu, fres.sigma, 29);
            mk2d->SetMarkerColor(kRed); mk2d->SetMarkerSize(2.0); mk2d->Draw();
        } else {
            drawNoDataText("No 2D scan stored");
        }
        c->cd(4); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        if (Config::ENABLE_POSITION_SMEARING && !fres.scan_data.sigma_pos_values.empty()) {
            drawChi2Scan1D(fres.scan_data.sigma_pos_values, fres.scan_data.chi2_vs_sigma_pos,
                           fres.sigma_pos, fres.chi2, "#sigma_{pos} [cm]", "#chi^{2} vs #sigma_{pos}");
        } else {
            drawNoDataText("Position scan not available");
        }
        c->cd(5); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        if (!fres.profile_diagnostics.empty()) {
            TMultiGraph *mg = new TMultiGraph("mg_profile_diag", "Profile scans;fixed parameter value;#Delta#chi^{2}");
            TLegend *lg = new TLegend(0.60, 0.62, 0.90, 0.88);
            lg->SetBorderSize(1); lg->SetFillColor(0); lg->SetTextSize(0.030);
            const char* pname[6] = {"mu_a", "mu_b", "mu_c", "sigma", "sigma_pos", "p_e_scale"};
            int colors[6] = {kRed, kBlue, kGreen+2, kMagenta, kOrange+7, kCyan+2};
            for (int ip = 0; ip < 6; ++ip) {
                vector<double> x, y;
                for (const auto &pd : fres.profile_diagnostics) {
                    if (pd.parameter == pname[ip] && std::isfinite(pd.chi2)) {
                        x.push_back(pd.fixed_value);
                        y.push_back(pd.chi2 - fres.chi2);
                    }
                }
                if (x.empty()) continue;
                TGraph *g = new TGraph(x.size(), &x[0], &y[0]);
                g->SetLineColor(colors[ip]); g->SetMarkerColor(colors[ip]);
                g->SetMarkerStyle(20); g->SetLineWidth(2);
                mg->Add(g, "LP");
                lg->AddEntry(g, pname[ip], "lp");
            }
            mg->Draw("A");
            lg->Draw();
        } else {
            drawNoDataText("No profile scans stored");
        }
        c->cd(6); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        if (!fres.seed_diagnostics.empty()) {
            vector<double> x, y_seed, y_migrad;
            int n = min(40, (int)fres.seed_diagnostics.size());
            for (int i = 0; i < n; ++i) {
                const auto &sd = fres.seed_diagnostics[i];
                x.push_back(i + 1);
                y_seed.push_back(sd.seed_chi2);
                y_migrad.push_back(sd.migrad_chi2);
            }
            TGraph *g_seed = new TGraph(x.size(), &x[0], &y_seed[0]);
            TGraph *g_migrad = new TGraph(x.size(), &x[0], &y_migrad[0]);
            TMultiGraph *mg = new TMultiGraph("mg_seed_diag", "Retained optimizer seeds;rank;#chi^{2}");
            g_seed->SetMarkerStyle(24); g_seed->SetMarkerColor(kBlue); g_seed->SetLineColor(kBlue);
            g_migrad->SetMarkerStyle(20); g_migrad->SetMarkerColor(kRed); g_migrad->SetLineColor(kRed);
            mg->Add(g_seed, "LP");
            mg->Add(g_migrad, "LP");
            mg->Draw("A");
            TLegend *lg = new TLegend(0.55, 0.70, 0.88, 0.88);
            lg->SetBorderSize(1); lg->SetFillColor(0); lg->SetTextSize(0.030);
            lg->AddEntry(g_seed, "seed #chi^{2}", "lp");
            lg->AddEntry(g_migrad, "post-MIGRAD #chi^{2}", "lp");
            lg->Draw();
        } else {
            drawNoDataText("No seed diagnostics stored");
        }
        c->Print(pdf_file.c_str());
        writeCanvasToDir(canvas_dir, c,
                         Form("c_%s_optimizer_diagnostics", sec.name().c_str()),
                         Form("%s optimizer diagnostics", sec.name().c_str()));
        ++page_count;
    }

    // ============================================================================
    // MAIN FUNCTION
    // ============================================================================
    int main(int argc, char** argv) {
        if (argc < 13) {
            cout << "Usage: " << argv[0] << " data.root dataTree sim.root simTree out.root nx ny x_min x_max y_min y_max overlap_frac [Nsmear]\n";
            cout << "Example: " << argv[0] << " data.root dataTree sim.root simTree out.root 8 8 -20 20 -20 20 0.2 20\n";
            return 1;
        }
        string data_file = argv[1];
        string data_tree_name = argv[2];
        string sim_file = argv[3];
        string sim_tree_name = argv[4];
        string out_file = argv[5];
        int nx = atoi(argv[6]);
        int ny = atoi(argv[7]);
        double x_min = atof(argv[8]);
        double x_max = atof(argv[9]);
        double y_min = atof(argv[10]);
        double y_max = atof(argv[11]);
        double overlap_frac = atof(argv[12]);
        int Nsmear = (argc >= 14) ? atoi(argv[13]) : 20;

        Config::print_configuration_summary();

        // Set output directory and prepend to relative paths
        if (out_file[0] != '/') {  // if relative path
            out_file = Config::OUTPUT_DIR + out_file;
        }
        const string run_tag = getRunTag();
        const string created_at_local = currentTimestampTag();
        const string timestamped_out_file = insertTagBeforeExtension(out_file, run_tag);
        cout << "Output file: " << out_file << endl;
        cout << "Diagnostic run tag: " << run_tag << endl;
        if (timestamped_out_file != out_file) {
            cout << "Timestamped output copy: " << timestamped_out_file << endl;
        }

        // open data and sim files
        TFile fdata(data_file.c_str(), "READ");
        if (fdata.IsZombie()) { cerr << "Cannot open data file " << data_file << endl; return 2; }
        TTree *tdata = dynamic_cast<TTree*>(fdata.Get(data_tree_name.c_str()));
        if (!tdata) { cerr << "Cannot find data tree '" << data_tree_name << "' in " << data_file << endl; return 3; }

        TFile fsim(sim_file.c_str(), "READ");
        if (fsim.IsZombie()) { cerr << "Cannot open sim file " << sim_file << endl; return 4; }
        TTree *tsim = dynamic_cast<TTree*>(fsim.Get(sim_tree_name.c_str()));
        if (!tsim) { cerr << "Cannot find sim tree '" << sim_tree_name << "' in " << sim_file << endl; return 5; }

        // ===== DATA TREE BRANCHES =====
        Double_t d_cluster_x_1, d_cluster_y_1, d_cluster_e_1;
        Double_t d_cluster_x_2, d_cluster_y_2, d_cluster_e_2;
        Double_t d_pi0_weight;
        Float_t d_scale;
        Int_t d_exclusive_flag = 1;
        Double_t d_mmiss_all = 0;  // Pre-calculated missing mass from data tree
        
        tdata->SetBranchAddress("cluster_x_1", &d_cluster_x_1);
        tdata->SetBranchAddress("cluster_y_1", &d_cluster_y_1);
        tdata->SetBranchAddress("cluster_e_1", &d_cluster_e_1);
        tdata->SetBranchAddress("cluster_x_2", &d_cluster_x_2);
        tdata->SetBranchAddress("cluster_y_2", &d_cluster_y_2);
        tdata->SetBranchAddress("cluster_e_2", &d_cluster_e_2);
        tdata->SetBranchAddress("pi0_weight", &d_pi0_weight);
        tdata->SetBranchAddress("scale", &d_scale);
        bool has_data_exclusive_branch = false;
        if (tdata->GetBranch(Config::DATA_EXCLUSIVITY_BRANCH)) {
            tdata->SetBranchAddress(Config::DATA_EXCLUSIVITY_BRANCH, &d_exclusive_flag);
            has_data_exclusive_branch = true;
        }
        
        // Read pre-calculated missing mass from data tree
        if (tdata->GetBranch("mmiss_all")) {
            tdata->SetBranchAddress("mmiss_all", &d_mmiss_all);
            cout << "Reading missing mass from branch: mmiss_all" << endl;
        } else {
            cerr << "WARNING: Missing mass branch 'mmiss_all' not found in data tree." << endl;
            cerr << "         Missing mass histograms will be empty." << endl;
        }

        // ===== SIMULATION TREE BRANCHES =====
        Float_t s_clust_X[10];  // Allow up to 10 clusters
        Float_t s_clust_Y[10];
        Float_t s_clust_E[10];
        Int_t s_nclust = 0;  // Number of clusters
        Float_t s_full_weight = 0.0f;
        Float_t s_sigcm = 0.0f;
        Int_t s_exclusive_flag = 1;
        
        // SIMC electron kinematics (for missing mass calculation)
        Float_t sc_e_E = 0, sc_e_Px = 0, sc_e_Py = 0, sc_e_Pz = 0;  // SIMC frame (MeV)
        Float_t s_e_E = 0, s_e_Px = 0, s_e_Py = 0, s_e_Pz = 0;  // Transformed (GeV)
        
        tsim->SetBranchAddress("clust_X", s_clust_X);
        tsim->SetBranchAddress("clust_Y", s_clust_Y);
        tsim->SetBranchAddress("clust_E", s_clust_E);
        tsim->SetBranchAddress("nclust", &s_nclust);
        if (!tsim->GetBranch("full_weight")) {
            cerr << "ERROR: SIMC branch 'full_weight' not found. "
                 << "Model-independent smearing requires full_weight/"
                 << Config::SIM_MODEL_XSEC_BRANCH << "." << endl;
            return 6;
        }
        if (!tsim->GetBranch(Config::SIM_MODEL_XSEC_BRANCH)) {
            cerr << "ERROR: SIMC model cross-section branch '"
                 << Config::SIM_MODEL_XSEC_BRANCH << "' not found. "
                 << "Model-independent smearing requires full_weight/"
                 << Config::SIM_MODEL_XSEC_BRANCH << "." << endl;
            return 6;
        }
        tsim->SetBranchAddress("full_weight", &s_full_weight);
        tsim->SetBranchAddress(Config::SIM_MODEL_XSEC_BRANCH, &s_sigcm);
        cout << "SIMC cross-section de-modeling: enabled\n"
             << "  sim weight = full_weight/" << Config::SIM_MODEL_XSEC_BRANCH << "\n"
             << "  invalid if |" << Config::SIM_MODEL_XSEC_BRANCH << "| < "
             << Config::SIM_MODEL_XSEC_MIN_ABS << endl;

        bool has_sim_exclusive_branch = false;
        if (tsim->GetBranch(Config::SIM_EXCLUSIVITY_BRANCH)) {
            tsim->SetBranchAddress(Config::SIM_EXCLUSIVITY_BRANCH, &s_exclusive_flag);
            has_sim_exclusive_branch = true;
        }

        if (Config::APPLY_IS_EXCLUSIVE_SELECTION) {
            cout << "Exclusivity weighting mode: enabled\n";
            cout << "  data branch: '" << Config::DATA_EXCLUSIVITY_BRANCH << "'\n";
            cout << "  sim branch:  '" << Config::SIM_EXCLUSIVITY_BRANCH << "'" << endl;
            if (!has_data_exclusive_branch) {
                cerr << "ERROR: Data branch '" << Config::DATA_EXCLUSIVITY_BRANCH
                     << "' not found. Re-run analysis/combine with 2D mass-cut branches." << endl;
                return 6;
            }
            if (!has_sim_exclusive_branch) {
                cerr << "ERROR: Simulation branch '" << Config::SIM_EXCLUSIVITY_BRANCH
                     << "' not found. Re-run scripts/simc_pi0_analysis.C before smearing." << endl;
                return 6;
            }
        } else {
            cout << "Exclusivity weighting mode: disabled (using all events)" << endl;
        }
        
        // Read SIMC electron kinematics
        bool has_sim_electron_branches = false;
        if (tsim->GetBranch("sc_e_E") && tsim->GetBranch("sc_e_Px") && 
            tsim->GetBranch("sc_e_Py") && tsim->GetBranch("sc_e_Pz")) {
            tsim->SetBranchAddress("sc_e_E", &sc_e_E);
            tsim->SetBranchAddress("sc_e_Px", &sc_e_Px);
            tsim->SetBranchAddress("sc_e_Py", &sc_e_Py);
            tsim->SetBranchAddress("sc_e_Pz", &sc_e_Pz);
            has_sim_electron_branches = true;
            cout << "Reading SIMC electron kinematics from branches: sc_e_E, sc_e_Px, sc_e_Py, sc_e_Pz" << endl;
        } else {
            cerr << "WARNING: SIMC electron kinematics branches (sc_e_E, sc_e_Px, sc_e_Py, sc_e_Pz) not found." << endl;
            cerr << "         Missing mass calculation will use zero momentum." << endl;
        }

        // NOTE: HMS electron cuts are applied in simc_pi0_analysis.C before this script
        // Input simulation data is pre-filtered, so no additional HMS cuts needed here

        // build sections with overlap
        vector<Section> sections;
        double total_wx = x_max - x_min;
        double total_wy = y_max - y_min;
        double base_wx = total_wx / nx;
        double base_wy = total_wy / ny;
        double wx = base_wx * (1.0 + overlap_frac);
        double wy = base_wy * (1.0 + overlap_frac);
        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                // center of cell without overlap
                double cx = x_min + (ix + 0.5) * base_wx;
                double cy = y_min + (iy + 0.5) * base_wy;
                Section s; s.ix = ix; s.iy = iy;
                s.xlo = cx - wx/2.0; s.xhi = cx + wx/2.0;
                s.ylo = cy - wy/2.0; s.yhi = cy + wy/2.0;
                // clamp to calorimeter bounds
                s.xlo = max(s.xlo, x_min); s.xhi = min(s.xhi, x_max);
                s.ylo = max(s.ylo, y_min); s.yhi = min(s.yhi, y_max);
                sections.push_back(s);
            }
        }
        int nsec = sections.size();
        cout << "Built " << nsec << " overlapping sections (nx="<<nx<<" ny="<<ny<<" overlap="<<overlap_frac<<")"<<endl;

        auto in_geometry = [&](double x, double y) {
            return (x >= x_min && x <= x_max && y >= y_min && y <= y_max);
        };
        auto event_touches_geometry = [&](double x1, double y1, double x2, double y2) {
            return in_geometry(x1, y1) || in_geometry(x2, y2);
        };

        // Global histograms/buffer for optional global p_e_scale fit (touches-geometry)
        TH1D hdata_global_mpi0("h_data_global_mpi0", "global data invariant mass",
                            Config::MGGAMMA_NBINS, Config::MGGAMMA_MIN, Config::MGGAMMA_MAX);
        TH1D hdata_global_mmiss("h_data_global_mmiss", "global data missing mass",
                                Config::MMISS_NBINS, Config::MMISS_MIN, Config::MMISS_MAX);
        TH1D hdata_global_mpgg2("h_data_global_mpgg2", "global data (p_{target}+#gamma#gamma)^{2}",
                                Config::MPGG2_NBINS, Config::MPGG2_MIN, Config::MPGG2_MAX);

        // Dedicated all-sections summary histograms/buffer (strict BOTH photons in geometry)
        TH1D hdata_summary_mpi0("h_data_summary_mpi0", "summary data invariant mass",
                                Config::MGGAMMA_NBINS, Config::MGGAMMA_MIN, Config::MGGAMMA_MAX);
        TH1D hdata_summary_mmiss("h_data_summary_mmiss", "summary data missing mass",
                                Config::MMISS_NBINS, Config::MMISS_MIN, Config::MMISS_MAX);
        TH1D hdata_summary_mpgg2("h_data_summary_mpgg2", "summary data (p_{target}+#gamma#gamma)^{2}",
                                Config::MPGG2_NBINS, Config::MPGG2_MIN, Config::MPGG2_MAX);

        // Track sum-of-weights-squared for correct weighted chi2
        hdata_global_mpi0.Sumw2();
        hdata_global_mmiss.Sumw2();
        hdata_global_mpgg2.Sumw2();
        hdata_summary_mpi0.Sumw2();
        hdata_summary_mmiss.Sumw2();
        hdata_summary_mpgg2.Sumw2();
        vector<ClusterPair> sim_events_global;
        vector<ClusterPair> sim_events_summary;
        vector<DataGlobalEvent> data_events_summary;

        // Prepare per-section data histograms and sim event buffers
        cout << "Scanning data tree and building per-section data histograms..." << endl;
        vector<vector<double>> data_mass_per_section(nsec);          // invariant mass M_γγ
        vector<vector<double>> data_mmiss_per_section(nsec);         // missing mass M_miss
        vector<vector<double>> data_mpgg2_per_section(nsec);         // (p_target + gamma+gamma)^2
        vector<vector<double>> data_weight_per_section(nsec);
        vector<int> data_selected_count_per_section(nsec, 0); // count with positive selected weight

        Long64_t ndata = tdata->GetEntries();
        for (Long64_t i=0;i<ndata;++i) {
            tdata->GetEntry(i);
            
            // Calculate invariant mass using nps_helper functions
            double m = nps::invariant_mass_pi0(d_cluster_e_1, d_cluster_e_2, 
                                            d_cluster_x_1, d_cluster_x_2, 
                                            d_cluster_y_1, d_cluster_y_2, 
                                            nps::kDefaultZ_NPS_cm);
            
            // Use pre-calculated missing mass from data tree
            double mmiss = d_mmiss_all;

            double d_cluster_e_1_mpgg2 = Config::W_MPGG2_ENERGY * d_cluster_e_1;
            double d_cluster_e_2_mpgg2 = Config::W_MPGG2_ENERGY * d_cluster_e_2;
            PhotonMomentum dp1 = computePhotonMomentum(d_cluster_e_1_mpgg2, d_cluster_x_1, d_cluster_y_1, nps::kDefaultZ_NPS_cm);
            PhotonMomentum dp2 = computePhotonMomentum(d_cluster_e_2_mpgg2, d_cluster_x_2, d_cluster_y_2, nps::kDefaultZ_NPS_cm);
            double mpgg2 = computeTargetPlusDiphotonMass2(d_cluster_e_1_mpgg2, dp1, d_cluster_e_2_mpgg2, dp2);
            
            // Calculate data event weight with optional exclusivity factor.
            double data_exclusive_factor = 1.0;
            if (Config::APPLY_IS_EXCLUSIVE_SELECTION && has_data_exclusive_branch) {
                data_exclusive_factor = static_cast<double>(d_exclusive_flag);
            }
            double weight = d_pi0_weight * d_scale * data_exclusive_factor;

            bool touches_geometry = event_touches_geometry(d_cluster_x_1, d_cluster_y_1,
                                    d_cluster_x_2, d_cluster_y_2);
            bool both_in_geometry = in_geometry(d_cluster_x_1, d_cluster_y_1) &&
                        in_geometry(d_cluster_x_2, d_cluster_y_2);

            // Fill global histograms only for events inside configured calorimeter geometry
            if (touches_geometry) {
                hdata_global_mpi0.Fill(m, weight);
                hdata_global_mmiss.Fill(mmiss, weight);
                hdata_global_mpgg2.Fill(mpgg2, weight);
            }

            // Fill dedicated all-sections summary pool only when BOTH photons are inside geometry
            if (both_in_geometry) {
                hdata_summary_mpi0.Fill(m, weight);
                hdata_summary_mmiss.Fill(mmiss, weight);
                hdata_summary_mpgg2.Fill(mpgg2, weight);
                data_events_summary.push_back({d_cluster_x_1, d_cluster_y_1,
                                            d_cluster_x_2, d_cluster_y_2,
                                            weight});
            }
            
            // assign to all sections where either photon cluster lies inside
            // only for events with BOTH photons inside the configured geometry
            for (int is=0; is<nsec; ++is) {
                if (both_in_geometry &&
                    (inSection(sections[is], d_cluster_x_1, d_cluster_y_1) || 
                    inSection(sections[is], d_cluster_x_2, d_cluster_y_2))) {
                    data_mass_per_section[is].push_back(m);
                    data_mmiss_per_section[is].push_back(mmiss);
                    data_mpgg2_per_section[is].push_back(mpgg2);
                    data_weight_per_section[is].push_back(weight);
                    if (weight > 0.0) ++data_selected_count_per_section[is];
                }
            }
        }

        // Check that we have data
        bool has_data = false;
        for (int is=0; is<nsec; ++is) {
            if (!data_mass_per_section[is].empty()) {
                has_data = true;
                break;
            }
        }
        if (!has_data) { 
            cerr << " No data events found in any section. Check cluster branches and calorimeter bounds." << endl; 
            return 6; 
        }
        
        cout << "Global invariant mass window for histograms: ["<<Config::MGGAMMA_MIN<<","<<Config::MGGAMMA_MAX<<"] nbins="<<Config::MGGAMMA_NBINS<<"\n";
        cout << "Missing mass window for histograms: ["<<Config::MMISS_MIN<<","<<Config::MMISS_MAX<<"] nbins="<<Config::MMISS_NBINS<<"\n";
        cout << "(p_target + #gamma#gamma)^2 window for histograms: ["<<Config::MPGG2_MIN<<","<<Config::MPGG2_MAX<<"] nbins="<<Config::MPGG2_NBINS<<"\n";

        // create TH1D per section for data (both invariant mass and missing mass) and fill with weights
        vector<TH1D*> hdata_sec(nsec,nullptr);        // Invariant mass M_γγ
        vector<TH1D*> hdata_mmiss_sec(nsec,nullptr);  // Missing mass M_miss
        vector<TH1D*> hdata_mpgg2_sec(nsec,nullptr);  // (p_target + gamma+gamma)^2
        
        for (int is=0; is<nsec; ++is) {
            // Invariant mass histogram
            string hname = string("h_data_mpi0_") + sections[is].name();
            hdata_sec[is] = new TH1D(hname.c_str(), hname.c_str(), Config::MGGAMMA_NBINS, Config::MGGAMMA_MIN, Config::MGGAMMA_MAX);
            hdata_sec[is]->Sumw2();
            for (size_t j=0; j<data_mass_per_section[is].size(); ++j) {
                hdata_sec[is]->Fill(data_mass_per_section[is][j], data_weight_per_section[is][j]);
            }
            
            // Missing mass histogram
            string hname_mmiss = string("h_data_mmiss_") + sections[is].name();
            hdata_mmiss_sec[is] = new TH1D(hname_mmiss.c_str(), hname_mmiss.c_str(), Config::MMISS_NBINS, Config::MMISS_MIN, Config::MMISS_MAX);
            hdata_mmiss_sec[is]->Sumw2();
            for (size_t j=0; j<data_mmiss_per_section[is].size(); ++j) {
                hdata_mmiss_sec[is]->Fill(data_mmiss_per_section[is][j], data_weight_per_section[is][j]);
            }

            string hname_mpgg2 = string("h_data_mpgg2_") + sections[is].name();
            hdata_mpgg2_sec[is] = new TH1D(hname_mpgg2.c_str(), hname_mpgg2.c_str(), Config::MPGG2_NBINS, Config::MPGG2_MIN, Config::MPGG2_MAX);
            hdata_mpgg2_sec[is]->Sumw2();
            for (size_t j=0; j<data_mpgg2_per_section[is].size(); ++j) {
                hdata_mpgg2_sec[is]->Fill(data_mpgg2_per_section[is][j], data_weight_per_section[is][j]);
            }
            
            cout << "Section "<<sections[is].name()<<" data entries="<<data_mass_per_section[is].size()
                << "  (selected=" << data_selected_count_per_section[is] << ")\n";
        }

        // Now load sim tree and build per-section sim event lists
        cout << "Scanning sim tree and building per-section sim event buffers..." << endl;
        vector<vector<ClusterPair>> sim_events_per_section(nsec);
        vector<int> sim_selected_count_per_section(nsec, 0);
        long long sim_skipped_invalid_model_weight = 0;
        Long64_t nsim = tsim->GetEntries();
        for (Long64_t i=0;i<nsim;++i) {
            tsim->GetEntry(i);
            
            // Transform SIMC electron kinematics to Hall C frame
            if (has_sim_electron_branches) {
                s_e_E = sc_e_E / 1000.0;        // MeV → GeV
                s_e_Px = sc_e_Py / 1000.0;      // Axis transformation: px_e = sc_e_Py
                s_e_Py = -sc_e_Px / 1000.0;     // Axis transformation: py_e = -sc_e_Px
                s_e_Pz = sc_e_Pz / 1000.0;      // MeV → GeV
            }
            
            // Check that we have at least 2 clusters
            if (s_nclust < 2) continue;
            
            // Get first two clusters
            ClusterPair pair;
            pair.e1 = s_clust_E[0];
            pair.x1 = s_clust_X[0];
            pair.y1 = s_clust_Y[0];
            pair.e2 = s_clust_E[1];
            pair.x2 = s_clust_X[1];
            pair.y2 = s_clust_Y[1];
            double sim_exclusive_factor = 1.0;
            if (Config::APPLY_IS_EXCLUSIVE_SELECTION && has_sim_exclusive_branch) {
                sim_exclusive_factor = static_cast<double>(s_exclusive_flag);
            }
            if (!std::isfinite(s_full_weight) ||
                !std::isfinite(s_sigcm) ||
                std::fabs(static_cast<double>(s_sigcm)) < Config::SIM_MODEL_XSEC_MIN_ABS) {
                ++sim_skipped_invalid_model_weight;
                continue;
            }
            double sim_base_weight = static_cast<double>(s_full_weight) / static_cast<double>(s_sigcm);
            if (!std::isfinite(sim_base_weight)) {
                ++sim_skipped_invalid_model_weight;
                continue;
            }
            pair.weight = sim_base_weight * sim_exclusive_factor;
            
            // Store transformed electron kinematics
            pair.Ee = s_e_E;
            pair.px_e = s_e_Px;
            pair.py_e = s_e_Py;
            pair.pz_e = s_e_Pz;

            bool touches_geometry = event_touches_geometry(pair.x1, pair.y1, pair.x2, pair.y2);
            bool both_in_geometry = in_geometry(pair.x1, pair.y1) && in_geometry(pair.x2, pair.y2);

            // Keep one global simulation buffer (no section duplication), geometry-gated
            if (touches_geometry) {
                sim_events_global.push_back(pair);
            }
            if (both_in_geometry) {
                sim_events_summary.push_back(pair);
            }
            
            for (int is=0; is<nsec; ++is) {
                const bool in1 = inSection(sections[is], pair.x1, pair.y1);
                const bool in2 = inSection(sections[is], pair.x2, pair.y2);
                if (both_in_geometry && (in1 || in2)) {
                    // Tag which photon(s) belong to this section so per-photon
                    // smearing is applied correctly in the chi2 evaluation functions.
                    ClusterPair sec_pair = pair;
                    sec_pair.photon1_in_section = in1;
                    sec_pair.photon2_in_section = in2;
                    sim_events_per_section[is].push_back(sec_pair);
                    if (sec_pair.weight > 0.0) ++sim_selected_count_per_section[is];
                }
            }
        }
        for (int is=0; is<nsec; ++is)
            cout << "Section "<<sections[is].name()<<" sim entries="<<sim_events_per_section[is].size()
                << "  (selected=" << sim_selected_count_per_section[is] << ")\n";
        cout << "SIMC cross-section de-modeling summary: weight=full_weight/"
             << Config::SIM_MODEL_XSEC_BRANCH
             << ", skipped invalid model-weight events="
             << sim_skipped_invalid_model_weight << "\n";

        const int optimizer_nsmear = max(1, min(Nsmear, Config::OPTIMIZATION_NSMEAR));
        vector<vector<ClusterPair>> sim_events_opt_per_section(nsec);
        long long opt_full_events_total = 0;
        long long opt_subset_events_total = 0;
        for (int is = 0; is < nsec; ++is) {
            sim_events_opt_per_section[is] = makeOptimizationSubset(
                sim_events_per_section[is],
                Config::OPT_MAX_SIM_EVENTS_PER_SECTION,
                Config::OPT_SUBSET_MGG_BINS);
            opt_full_events_total += (long long)sim_events_per_section[is].size();
            opt_subset_events_total += (long long)sim_events_opt_per_section[is].size();
        }
        cout << "\n==== Optimizer event thinning ====\n"
             << "Optimizer Nsmear=" << optimizer_nsmear
             << " (final Nsmear=" << Nsmear << ")\n"
             << "Section optimizer sim events: " << opt_subset_events_total
             << " / " << opt_full_events_total
             << " after deterministic M_gg-stratified, weight-compensated thinning.\n"
             << "Final chi2 recompute and output histograms still use full sim events.\n";

        // RNG - will be created per-thread in parallel region
        TRandom3 rng(0);  // Main thread RNG

        // Determine global electron momentum scale once (shared across sections)
        bool use_global_p_e_scale = Config::ENABLE_ELECTRON_MOMENTUM_SCALING &&
                                    Config::fit_objective_uses_mmiss();
        double global_p_e_scale = Config::GLOBAL_PE_SCALE_DEFAULT;

        if (use_global_p_e_scale) {
            cout << "\n==== Global electron momentum scale mode enabled ====" << endl;
            if (Config::ENABLE_FINAL_GLOBAL_PE_SCALE_STAGE4) {
                cout << "Coupled sweeps will jointly refine per-section p_e_scale with (mu, sigma, sigma_pos), then final global p_e_scale refinement will run." << endl;
            } else {
                cout << "Coupled sweeps will jointly refine per-section p_e_scale with (mu, sigma, sigma_pos); final global p_e_scale refinement is disabled." << endl;
            }
        }

        const string sweep_acceptance_strategy = Config::sweep_acceptance_strategy();
        const bool use_accepted_state_rollback =
            (sweep_acceptance_strategy == "jacobi_global_accept_rollback");

	        // Output file
	        TFile fout(out_file.c_str(), "RECREATE");
	        fout.cd();
	        TNamed("run_tag", run_tag.c_str()).Write();
	        TNamed("created_at_local", created_at_local.c_str()).Write();
	        TNamed("timestamped_output_file", timestamped_out_file.c_str()).Write();
	        TNamed("smearing_model_energy_mean", "a_plus_bE_plus_clnE").Write();
	        TNamed("energy_mean_convention", "reconstructed_energy_GeV").Write();
	        TNamed("energy_mean_formula", "E_mean = a + b*E_safe + c*ln(E_safe/1 GeV)").Write();
	        TNamed("energy_log_floor_GeV", Form("%.8g", Config::MU_ENERGY_MIN_GEV)).Write();
	        TNamed("energy_sigma_model", Config::USE_SIMPLE_STOCHASTIC_MODEL ? "sigma_sqrt_Emean" : "scaled_three_term").Write();
	        TNamed("position_sigma_model", Config::ENABLE_ENERGY_DEPENDENT_SIGMA_POS ? "sigma_pos_sqrt_E0_over_Emean" : "constant").Write();
	        TNamed("optimizer_model", "sobol_multistart_global_joint_section_energy_sigma_then_sigma_pos").Write();
	        TNamed("optimizer_global_seed_count", Form("%d", Config::GLOBAL_MULTISTART_SEEDS)).Write();
	        TNamed("optimizer_keep_best_seeds", Form("%d", Config::GLOBAL_MULTISTART_KEEP_BEST)).Write();
	        TNamed("optimizer_nsmear", Form("%d", optimizer_nsmear)).Write();
	        TNamed("final_nsmear", Form("%d", Nsmear)).Write();
	        TNamed("optimizer_max_sim_events_per_section", Form("%d", Config::OPT_MAX_SIM_EVENTS_PER_SECTION)).Write();
	        TNamed("optimizer_max_sim_events_global_prefit", Form("%d", Config::OPT_MAX_SIM_EVENTS_GLOBAL_PREFIT)).Write();
	        TNamed("coupled_sweep_acceptance_strategy", sweep_acceptance_strategy.c_str()).Write();
	        TNamed("sim_weight_mode", "full_weight_over_sigcm").Write();
	        TNamed("sim_model_xsec_branch", Config::SIM_MODEL_XSEC_BRANCH).Write();
	        TNamed("sim_model_xsec_min_abs", Form("%.17g", Config::SIM_MODEL_XSEC_MIN_ABS)).Write();
	        TNamed("sim_skipped_invalid_model_weight_events", Form("%lld", sim_skipped_invalid_model_weight)).Write();
            TDirectory *diagnostic_canvas_dir = fout.mkdir("diagnostic_canvases");
            TDirectory *diagnostic_map_dir = fout.mkdir("diagnostic_maps");
            fout.cd();

	        // CSV summary - place in same directory as ROOT output
        string csv_file = out_file;
        size_t last_slash = csv_file.find_last_of('/');
        if (last_slash != string::npos) {
            csv_file = csv_file.substr(0, last_slash + 1) + Config::CSV_FILENAME;
        } else {
            csv_file = Config::CSV_FILENAME;
        }
        auto sibling_output_path = [&](const string &filename) {
            if (last_slash != string::npos) {
                return out_file.substr(0, last_slash + 1) + filename;
            }
            return filename;
        };
        string optimizer_summary_csv_file = sibling_output_path(Config::OPTIMIZER_SUMMARY_CSV_FILENAME);
        string optimizer_seeds_csv_file = sibling_output_path(Config::OPTIMIZER_SEEDS_CSV_FILENAME);
        string optimizer_profile_csv_file = sibling_output_path(Config::OPTIMIZER_PROFILE_CSV_FILENAME);
        string closure_summary_csv_file = sibling_output_path(Config::CLOSURE_SUMMARY_CSV_FILENAME);
        string sweep_history_csv_file = sibling_output_path(Config::SWEEP_HISTORY_CSV_FILENAME);
        string objective_breakdown_csv_file = sibling_output_path(Config::OBJECTIVE_BREAKDOWN_CSV_FILENAME);
        string cache_fingerprint_file = sibling_output_path(Config::CACHE_FINGERPRINT_FILENAME);
        auto fnv1a_update = [](unsigned long long h, const char *data, size_t n) {
            const unsigned long long prime = 1099511628211ull;
            for (size_t i = 0; i < n; ++i) {
                h ^= (unsigned char)data[i];
                h *= prime;
            }
            return h;
        };
        auto hash_file = [&](const string &path) {
            unsigned long long h = 1469598103934665603ull;
            ifstream in(path.c_str(), std::ios::binary);
            if (!in.good()) return 0ull;
            char buf[8192];
            while (in.good()) {
                in.read(buf, sizeof(buf));
                std::streamsize n = in.gcount();
                if (n > 0) h = fnv1a_update(h, buf, (size_t)n);
            }
            return h;
        };
        auto file_stat_signature = [&](const string &path) {
            struct stat st;
            ostringstream os;
            if (stat(path.c_str(), &st) == 0) {
                os << path << "|size=" << (long long)st.st_size
                   << "|mtime=" << (long long)st.st_mtime;
            } else {
                os << path << "|missing";
            }
            return os.str();
        };
        string source_path_for_hash = "scripts/nps_sim_smearing_new_try.C";
        unsigned long long source_hash = hash_file(source_path_for_hash);
        if (source_hash == 0ull) {
            source_path_for_hash = __FILE__;
            source_hash = hash_file(source_path_for_hash);
        }
        ostringstream cache_fp;
        cache_fp << "cache_version=section_sweep_mass_plots_v1\n"
                 << "source_path=" << source_path_for_hash << "\n"
                 << "source_hash_fnv1a=" << source_hash << "\n"
                 << "data=" << file_stat_signature(data_file) << "\n"
                 << "sim=" << file_stat_signature(sim_file) << "\n"
                 << "data_tree=" << data_tree_name << "\n"
                 << "sim_tree=" << sim_tree_name << "\n"
                 << "nx=" << nx << "\n"
                 << "ny=" << ny << "\n"
                 << "x_min=" << std::setprecision(17) << x_min << "\n"
                 << "x_max=" << std::setprecision(17) << x_max << "\n"
                 << "y_min=" << std::setprecision(17) << y_min << "\n"
                 << "y_max=" << std::setprecision(17) << y_max << "\n"
                 << "overlap_frac=" << std::setprecision(17) << overlap_frac << "\n"
                 << "Nsmear=" << Nsmear << "\n";
        const string current_cache_fingerprint = cache_fp.str();
        fout.cd();
        TNamed("cache_fingerprint_current", current_cache_fingerprint.c_str()).Write();
        // Store fit results for interpolation
        vector<FitResult> fit_results(nsec);
        vector<bool> fit_success(nsec, false);
        vector<string> fit_status(nsec, "not_fitted");

        auto optimize_global_p_e_scale_for_fit =
            [&](const vector<FitResult> &src_results,
                const vector<bool> &src_success,
                double &io_global_p_e_scale,
                const string &tag,
                double *out_best_chi2 = nullptr) -> bool {
                auto section_center_x = [&](int is) {
                    return x_min + (sections[is].ix + 0.5) * base_wx;
                };
                auto section_center_y = [&](int is) {
                    return y_min + (sections[is].iy + 0.5) * base_wy;
                };
                auto best_section_for_point = [&](double x, double y) {
                    int best_is = -1;
                    double best_d2 = std::numeric_limits<double>::max();
                    for (int is = 0; is < nsec; ++is) {
                        if (!src_success[is]) continue;
                        if (!inSection(sections[is], x, y)) continue;
                        double dx = x - section_center_x(is);
                        double dy = y - section_center_y(is);
                        double d2 = dx * dx + dy * dy;
                        if (d2 < best_d2) {
                            best_d2 = d2;
                            best_is = is;
                        }
                    }
                    return best_is;
                };

                vector<ClusterPair> sim_events_global_stage4;
                sim_events_global_stage4.reserve(sim_events_global.size());
                for (const auto &ev : sim_events_global) {
                    ClusterPair gev = ev;
                    gev.photon1_in_section = false;
                    gev.photon2_in_section = false;
                    gev.mu_a1_ext = 0.0; gev.mu1_ext = 1.0; gev.mu_c1_ext = 0.0;
                    gev.sigma1_ext = 0.0; gev.sigma_pos1_ext = 0.0;
                    gev.mu_a2_ext = 0.0; gev.mu2_ext = 1.0; gev.mu_c2_ext = 0.0;
                    gev.sigma2_ext = 0.0; gev.sigma_pos2_ext = 0.0;

                    int is1 = in_geometry(ev.x1, ev.y1) ? best_section_for_point(ev.x1, ev.y1) : -1;
                    if (is1 >= 0) {
                        gev.mu_a1_ext = src_results[is1].mu_a;
                        gev.mu1_ext   = src_results[is1].mu;
                        gev.mu_c1_ext = src_results[is1].mu_c;
                        gev.sigma1_ext = src_results[is1].sigma;
                        gev.sigma_pos1_ext = src_results[is1].sigma_pos;
                    }
                    int is2 = in_geometry(ev.x2, ev.y2) ? best_section_for_point(ev.x2, ev.y2) : -1;
                    if (is2 >= 0) {
                        gev.mu_a2_ext = src_results[is2].mu_a;
                        gev.mu2_ext   = src_results[is2].mu;
                        gev.mu_c2_ext = src_results[is2].mu_c;
                        gev.sigma2_ext = src_results[is2].sigma;
                        gev.sigma_pos2_ext = src_results[is2].sigma_pos;
                    }
                    sim_events_global_stage4.push_back(gev);
                }

                if (sim_events_global_stage4.empty() || hdata_global_mmiss.Integral() <= 0.0) {
                    cout << tag << ": skipped global p_e_scale optimization (insufficient global M_miss statistics)." << endl;
                    return false;
                }

                double best_chi2_global = 1e300;
                TRandom3 rng_stage4(12345);
                double best_scale = io_global_p_e_scale;

                int PE_COARSE = max(9, Config::PE_SCALE_NSTEPS / Config::COARSE_GRID_DIVISOR);
                double pe_step = (Config::PE_SCALE_MAX - Config::PE_SCALE_MIN) / (PE_COARSE - 1);
                for (int i = 0; i < PE_COARSE; ++i) {
                    double pe = Config::PE_SCALE_MIN + i * pe_step;
                    double chi2 = eval_chi2_mmiss_only(0.0, 1.0, 0.0, 0.0, 0.0, pe,
                                                    sim_events_global_stage4,
                                                    hdata_global_mmiss,
                                                    rng_stage4,
                                                    Nsmear,
                                                    Config::RESOLUTION_A_DEFAULT,
                                                    Config::RESOLUTION_B_DEFAULT,
                                                    Config::RESOLUTION_C_DEFAULT);
                    if (chi2 < best_chi2_global) {
                        best_chi2_global = chi2;
                        best_scale = pe;
                    }
                }

                double refine_step = (Config::PE_SCALE_MAX - Config::PE_SCALE_MIN) / 20.0;
                int refine_iter = 0;
                while (refine_step >= 1e-4 && refine_iter < Config::MAX_REFINEMENT_ITERATIONS) {
                    double grid_best_scale = best_scale;
                    double grid_best_chi2 = best_chi2_global;
                    for (double pe = max(Config::PE_SCALE_MIN, best_scale - 2 * refine_step);
                        pe <= min(Config::PE_SCALE_MAX, best_scale + 2 * refine_step) + 1e-15;
                        pe += refine_step) {
                        double chi2 = eval_chi2_mmiss_only(0.0, 1.0, 0.0, 0.0, 0.0, pe,
                                                        sim_events_global_stage4,
                                                        hdata_global_mmiss,
                                                        rng_stage4,
                                                        Nsmear,
                                                        Config::RESOLUTION_A_DEFAULT,
                                                        Config::RESOLUTION_B_DEFAULT,
                                                        Config::RESOLUTION_C_DEFAULT);
                        if (chi2 < grid_best_chi2) {
                            grid_best_chi2 = chi2;
                            grid_best_scale = pe;
                        }
                    }
                    best_scale = grid_best_scale;
                    best_chi2_global = grid_best_chi2;
                    refine_step /= 2.0;
                    refine_iter++;
                }

                {
                    auto chi2_eval = [&](double /*mu_a*/, double /*mu_b*/, double /*mu_c*/,
                                        double /*sigma*/, double /*sigma_pos*/, double p_e_scale) {
                        TRandom3 rng_eval(456789);
                        return eval_chi2_mmiss_only(0.0, 1.0, 0.0, 0.0, 0.0, p_e_scale,
                                                    sim_events_global_stage4,
                                                    hdata_global_mmiss,
                                                    rng_eval,
                                                    Nsmear,
                                                    Config::RESOLUTION_A_DEFAULT,
                                                    Config::RESOLUTION_B_DEFAULT,
                                                    Config::RESOLUTION_C_DEFAULT);
                    };
                    MigradRefineResult mr = run_migrad_refinement(chi2_eval,
                                                                0.0, 1.0, 0.0, 0.0, 0.0, best_scale,
                                                                false, false, false, false, false, true);
                    if (mr.minimized) {
                        best_scale = mr.p_e_scale;
                        best_chi2_global = mr.chi2;
                    }
                }

                if (refine_iter >= Config::MAX_REFINEMENT_ITERATIONS) {
                    cout << "[WARN] " << tag << ": global p_e_scale refinement reached MAX_REFINEMENT_ITERATIONS="
                        << Config::MAX_REFINEMENT_ITERATIONS << endl;
                }

                io_global_p_e_scale = best_scale;
                if (out_best_chi2) *out_best_chi2 = best_chi2_global;
                return true;
            };

        // Create canvas and PDF for chi^2 plots
        string canonical_pdf_file = sibling_output_path(Config::CHI2_PDF_FILENAME);
        string pdf_file = insertTagBeforeExtension(canonical_pdf_file, run_tag);
        TCanvas *c_chi2 = new TCanvas("c_chi2", "Section Smearing Diagnostics", 1600, 1100);
        string pdf_open = pdf_file + "[";
        c_chi2->Print(pdf_open.c_str());
        int page_count = 0;
        string progress_pdf_dir = sibling_output_path("chi2_scans_progress_" + run_tag);
        gSystem->mkdir(progress_pdf_dir.c_str(), true);
        string metadata_manifest_file = sibling_output_path("smearing_artifacts_" + run_tag + ".json");
        string interp_file = interpolatedPathForOutput(out_file);
        string timestamped_interp_file = interpolatedPathForOutput(timestamped_out_file);
        fout.cd();
        TNamed("chi2_pdf_file", pdf_file.c_str()).Write();
        TNamed("canonical_chi2_pdf_file", canonical_pdf_file.c_str()).Write();
        TNamed("chi2_progress_dir", progress_pdf_dir.c_str()).Write();
        TNamed("smearing_artifact_manifest", metadata_manifest_file.c_str()).Write();
        TNamed("interpolated_output_file", interp_file.c_str()).Write();
        TNamed("timestamped_interpolated_output_file", timestamped_interp_file.c_str()).Write();

        // Loop over sections and fit each where we have enough stats
        // OpenMP parallelization: each section is fitted independently
        
        // ROOT is not thread-safe by default. Disable automatic histogram registration
        Bool_t original_adddir_status = TH1::AddDirectoryStatus();
        TH1::AddDirectory(kFALSE);
        
        // Report number of threads
        int num_threads = 1;
        #pragma omp parallel if(Config::ENABLE_PARALLEL_SECTION_FITS)
        {
            #pragma omp single
            {
                num_threads = omp_get_num_threads();
                cout << "\n==== Starting parallel fitting with " << num_threads << " threads ====" << endl;
            }
        }

        cout << "\n==== Iterative coupled section fitting ====\n";
        cout << "Global all-calorimeter prefit seeds coupled sweeps." << endl;

        // =========================================================================
        // === Cross-boundary ext parameter assignment
        // =========================================================================
        // Shared helper for sweep N -> sweep N+1 handoff.  It assigns
        // out-of-section photon coefficients from completed section results.
        auto refresh_cross_boundary_ext = [&](const vector<FitResult> &src_results,
                                            const vector<bool> &src_success) {
            CalibrationMap coupled_map(nx, ny, x_min, x_max, y_min, y_max);
            for (int js = 0; js < nsec; ++js) {
                if (!src_success[js]) continue;
                coupled_map.setParams(sections[js].ix, sections[js].iy,
                                    src_results[js].mu_a,
                                    src_results[js].mu,
                                    src_results[js].mu_c,
                                    src_results[js].sigma,
                                    src_results[js].sigma_pos);
            }

            auto best_fitted_overlap_section = [&](double x, double y) {
                int best_is = -1;
                double best_d2 = std::numeric_limits<double>::max();
                for (int js = 0; js < nsec; ++js) {
                    if (!src_success[js]) continue;
                    if (!inSection(sections[js], x, y)) continue;
                    double cx = x_min + (sections[js].ix + 0.5) * base_wx;
                    double cy = y_min + (sections[js].iy + 0.5) * base_wy;
                    double dx = x - cx;
                    double dy = y - cy;
                    double d2 = dx * dx + dy * dy;
                    if (d2 < best_d2) {
                        best_d2 = d2;
                        best_is = js;
                    }
                }
                return best_is;
            };

            auto base_grid_section_index = [&](double x, double y) {
                if (!in_geometry(x, y)) return -1;
                int ix = (int)floor((x - x_min) / base_wx);
                int iy = (int)floor((y - y_min) / base_wy);
                ix = max(0, min(nx - 1, ix));
                iy = max(0, min(ny - 1, iy));
                return iy * nx + ix;
            };

            int assigned = 0;
            int interpolated = 0;
            auto refresh_event = [&](ClusterPair &ev, bool count_this) {
                if (!ev.photon1_in_section) {
                    int js = best_fitted_overlap_section(ev.x1, ev.y1);
                    if (js < 0) {
                        int jb = base_grid_section_index(ev.x1, ev.y1);
                        if (jb >= 0 && src_success[jb]) js = jb;
                    }
                    if (js >= 0) {
                        ev.mu_a1_ext      = src_results[js].mu_a;
                        ev.mu1_ext        = src_results[js].mu;
                        ev.mu_c1_ext      = src_results[js].mu_c;
                        ev.sigma1_ext     = src_results[js].sigma;
                        ev.sigma_pos1_ext = src_results[js].sigma_pos;
                    } else if (in_geometry(ev.x1, ev.y1)) {
                        coupled_map.getInterpolatedParams(ev.x1, ev.y1,
                                                        ev.mu_a1_ext, ev.mu1_ext, ev.mu_c1_ext,
                                                        ev.sigma1_ext, ev.sigma_pos1_ext);
                        if (count_this) ++interpolated;
                    }
                    if (count_this) ++assigned;
                }
                if (!ev.photon2_in_section) {
                    int js = best_fitted_overlap_section(ev.x2, ev.y2);
                    if (js < 0) {
                        int jb = base_grid_section_index(ev.x2, ev.y2);
                        if (jb >= 0 && src_success[jb]) js = jb;
                    }
                    if (js >= 0) {
                        ev.mu_a2_ext      = src_results[js].mu_a;
                        ev.mu2_ext        = src_results[js].mu;
                        ev.mu_c2_ext      = src_results[js].mu_c;
                        ev.sigma2_ext     = src_results[js].sigma;
                        ev.sigma_pos2_ext = src_results[js].sigma_pos;
                    } else if (in_geometry(ev.x2, ev.y2)) {
                        coupled_map.getInterpolatedParams(ev.x2, ev.y2,
                                                        ev.mu_a2_ext, ev.mu2_ext, ev.mu_c2_ext,
                                                        ev.sigma2_ext, ev.sigma_pos2_ext);
                        if (count_this) ++interpolated;
                    }
                    if (count_this) ++assigned;
                }
            };
            for (int is = 0; is < nsec; ++is) {
                for (auto &ev : sim_events_per_section[is]) refresh_event(ev, true);
                for (auto &ev : sim_events_opt_per_section[is]) refresh_event(ev, false);
            }
            return std::make_pair(assigned, interpolated);
        };

        auto reset_cross_boundary_ext_to_seed = [&](const FitResult &seed) {
            int assigned = 0;
            auto reset_event = [&](ClusterPair &ev, bool count_this) {
                auto reset_photon = [&](double &ma, double &mb, double &mc,
                                        double &sig, double &spos) {
                    ma = seed.mu_a;
                    mb = seed.mu;
                    mc = seed.mu_c;
                    sig = seed.sigma;
                    spos = seed.sigma_pos;
                    if (count_this) ++assigned;
                };
                if (!ev.photon1_in_section) {
                    reset_photon(ev.mu_a1_ext, ev.mu1_ext, ev.mu_c1_ext,
                                 ev.sigma1_ext, ev.sigma_pos1_ext);
                }
                if (!ev.photon2_in_section) {
                    reset_photon(ev.mu_a2_ext, ev.mu2_ext, ev.mu_c2_ext,
                                 ev.sigma2_ext, ev.sigma_pos2_ext);
                }
            };
            for (int is = 0; is < nsec; ++is) {
                for (auto &ev : sim_events_per_section[is]) reset_event(ev, true);
                for (auto &ev : sim_events_opt_per_section[is]) reset_event(ev, false);
            }
            return assigned;
        };

        auto build_global_events_with_fit =
            [&](const vector<FitResult> &src_results,
                const vector<bool> &src_success,
                int *out_interpolated_count = nullptr) {
            auto section_center_x_local = [&](int is) {
                return x_min + (sections[is].ix + 0.5) * base_wx;
            };
            auto section_center_y_local = [&](int is) {
                return y_min + (sections[is].iy + 0.5) * base_wy;
            };
            auto base_grid_section_index_local = [&](double x, double y) {
                if (!in_geometry(x, y)) return -1;
                int ix = (int)floor((x - x_min) / base_wx);
                int iy = (int)floor((y - y_min) / base_wy);
                ix = max(0, min(nx - 1, ix));
                iy = max(0, min(ny - 1, iy));
                return iy * nx + ix;
            };
            auto section_index_for_photon = [&](double x, double y) {
                if (!in_geometry(x, y)) return -1;
                int best_is = -1;
                double best_d2 = std::numeric_limits<double>::max();
                for (int is = 0; is < nsec; ++is) {
                    if (!src_success[is]) continue;
                    if (!inSection(sections[is], x, y)) continue;
                    double dx = x - section_center_x_local(is);
                    double dy = y - section_center_y_local(is);
                    double d2 = dx * dx + dy * dy;
                    if (d2 < best_d2) {
                        best_d2 = d2;
                        best_is = is;
                    }
                }
                if (best_is >= 0) return best_is;
                int base_is = base_grid_section_index_local(x, y);
                if (base_is >= 0 && src_success[base_is]) return base_is;
                return -1;
            };

            CalibrationMap cal_map(nx, ny, x_min, x_max, y_min, y_max);
            for (int is = 0; is < nsec; ++is) {
                if (!src_success[is]) continue;
                cal_map.setParams(sections[is].ix, sections[is].iy,
                                  src_results[is].mu_a,
                                  src_results[is].mu,
                                  src_results[is].mu_c,
                                  src_results[is].sigma,
                                  src_results[is].sigma_pos);
            }

            int interpolated_count = 0;
            vector<ClusterPair> out;
            out.reserve(sim_events_summary.size());
            for (const auto &ev : sim_events_summary) {
                ClusterPair gev = ev;
                gev.photon1_in_section = false;
                gev.photon2_in_section = false;
                gev.mu_a1_ext = 0.0; gev.mu1_ext = 1.0; gev.mu_c1_ext = 0.0;
                gev.sigma1_ext = 0.0; gev.sigma_pos1_ext = 0.0;
                gev.mu_a2_ext = 0.0; gev.mu2_ext = 1.0; gev.mu_c2_ext = 0.0;
                gev.sigma2_ext = 0.0; gev.sigma_pos2_ext = 0.0;

                int is1 = section_index_for_photon(ev.x1, ev.y1);
                if (is1 >= 0) {
                    gev.mu_a1_ext = src_results[is1].mu_a;
                    gev.mu1_ext = src_results[is1].mu;
                    gev.mu_c1_ext = src_results[is1].mu_c;
                    gev.sigma1_ext = src_results[is1].sigma;
                    gev.sigma_pos1_ext = src_results[is1].sigma_pos;
                } else if (in_geometry(ev.x1, ev.y1)) {
                    cal_map.getInterpolatedParams(ev.x1, ev.y1,
                                                  gev.mu_a1_ext, gev.mu1_ext, gev.mu_c1_ext,
                                                  gev.sigma1_ext, gev.sigma_pos1_ext);
                    ++interpolated_count;
                }

                int is2 = section_index_for_photon(ev.x2, ev.y2);
                if (is2 >= 0) {
                    gev.mu_a2_ext = src_results[is2].mu_a;
                    gev.mu2_ext = src_results[is2].mu;
                    gev.mu_c2_ext = src_results[is2].mu_c;
                    gev.sigma2_ext = src_results[is2].sigma;
                    gev.sigma_pos2_ext = src_results[is2].sigma_pos;
                } else if (in_geometry(ev.x2, ev.y2)) {
                    cal_map.getInterpolatedParams(ev.x2, ev.y2,
                                                  gev.mu_a2_ext, gev.mu2_ext, gev.mu_c2_ext,
                                                  gev.sigma2_ext, gev.sigma_pos2_ext);
                    ++interpolated_count;
                }
                out.push_back(gev);
            }
            if (out_interpolated_count) *out_interpolated_count = interpolated_count;
            return out;
        };

        auto evaluate_global_objective =
            [&](const vector<FitResult> &src_results,
                const vector<bool> &src_success,
                int nsmear_eval,
                ObjectiveBreakdown *out_breakdown = nullptr,
                int *out_interpolated_count = nullptr) {
            vector<ClusterPair> global_events = build_global_events_with_fit(src_results,
                                                                             src_success,
                                                                             out_interpolated_count);
            TRandom3 rng_global_obj(910247);
            ObjectiveBreakdown breakdown = eval_objective_breakdown(
                0.0, 1.0, 0.0,
                0.0, 0.0,
                global_p_e_scale,
                global_events,
                hdata_summary_mpi0, hdata_summary_mmiss, hdata_summary_mpgg2,
                rng_global_obj, nsmear_eval,
                Config::RESOLUTION_A_DEFAULT,
                Config::RESOLUTION_B_DEFAULT,
                Config::RESOLUTION_C_DEFAULT);
            if (out_breakdown) *out_breakdown = breakdown;
            return breakdown.total(Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);
        };

        auto plot_iteration_candidate_histograms =
            [&](const string &iteration_label,
                const string &file_tag,
                const vector<FitResult> &src_results,
                const vector<bool> &src_success,
                const ObjectiveBreakdown &breakdown,
                double objective,
                bool candidate_accepted,
                bool accepted_best,
                const string &sweep_note,
                int interpolated_count) {
            vector<ClusterPair> global_events = build_global_events_with_fit(src_results,
                                                                             src_success,
                                                                             nullptr);
            TH1D hu_mgg(Form("hu_%s_mgg", file_tag.c_str()), ";M_{#gamma#gamma} [GeV/c^{2}];Counts",
                        Config::MGGAMMA_NBINS, Config::MGGAMMA_MIN, Config::MGGAMMA_MAX);
            TH1D hu_mmiss(Form("hu_%s_mmiss", file_tag.c_str()), ";M_{miss} [GeV/c^{2}];Counts",
                          Config::MMISS_NBINS, Config::MMISS_MIN, Config::MMISS_MAX);
            TH1D hu_mpgg2(Form("hu_%s_mpgg2", file_tag.c_str()), ";(p_{target}+#gamma#gamma)^{2} [GeV^{2}];Counts",
                          Config::MPGG2_NBINS, Config::MPGG2_MIN, Config::MPGG2_MAX);
            TH1D hs_mgg(Form("hs_%s_mgg", file_tag.c_str()), ";M_{#gamma#gamma} [GeV/c^{2}];Counts",
                        Config::MGGAMMA_NBINS, Config::MGGAMMA_MIN, Config::MGGAMMA_MAX);
            TH1D hs_mmiss(Form("hs_%s_mmiss", file_tag.c_str()), ";M_{miss} [GeV/c^{2}];Counts",
                          Config::MMISS_NBINS, Config::MMISS_MIN, Config::MMISS_MAX);
            TH1D hs_mpgg2(Form("hs_%s_mpgg2", file_tag.c_str()), ";(p_{target}+#gamma#gamma)^{2} [GeV^{2}];Counts",
                          Config::MPGG2_NBINS, Config::MPGG2_MIN, Config::MPGG2_MAX);
            hu_mgg.Sumw2(); hu_mmiss.Sumw2(); hu_mpgg2.Sumw2();
            hs_mgg.Sumw2(); hs_mmiss.Sumw2(); hs_mpgg2.Sumw2();

            fillUnsmearedHistogramsForNormalization(global_events, hu_mgg, hu_mmiss, hu_mpgg2,
                                                    global_p_e_scale);
            unsigned int plot_seed = 700000;
            for (char ch : file_tag) plot_seed = plot_seed * 131u + (unsigned int)(unsigned char)ch;
            TRandom3 rng_sweep_plot(plot_seed);
            fillSmearedHistogramsAtParams(global_events,
                                          hs_mgg, hs_mmiss, hs_mpgg2,
                                          0.0, 1.0, 0.0,
                                          0.0, 0.0,
                                          global_p_e_scale,
                                          rng_sweep_plot,
                                          Nsmear,
                                          Config::RESOLUTION_A_DEFAULT,
                                          Config::RESOLUTION_B_DEFAULT,
                                          Config::RESOLUTION_C_DEFAULT);

            auto scale_pair = [](TH1D &hu, TH1D &hs, const TH1D &hd) {
                if (hu.Integral() <= 0.0 || hd.Integral() <= 0.0) return;
                double scale = hd.Integral() / hu.Integral();
                hu.Scale(scale);
                hs.Scale(scale);
            };
            scale_pair(hu_mgg, hs_mgg, hdata_summary_mpi0);
            scale_pair(hu_mmiss, hs_mmiss, hdata_summary_mmiss);
            scale_pair(hu_mpgg2, hs_mpgg2, hdata_summary_mpgg2);

            c_chi2->Clear();
            c_chi2->Divide(3, 2);
            c_chi2->cd(1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
            drawComparisonPad(hdata_summary_mpi0, hs_mgg, hu_mgg,
                              "M_{#gamma#gamma} [GeV/c^{2}]", Form("%s_mgg", file_tag.c_str()));
            TLatex title; title.SetNDC(); title.SetTextSize(0.045);
            title.DrawLatex(0.15, 0.93, Form("%s: M_{#gamma#gamma}", iteration_label.c_str()));

            c_chi2->cd(2); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
            drawComparisonPad(hdata_summary_mmiss, hs_mmiss, hu_mmiss,
                              "M_{miss} [GeV/c^{2}]", Form("%s_mmiss", file_tag.c_str()));
            title.DrawLatex(0.15, 0.93, Form("%s: M_{miss}", iteration_label.c_str()));

            c_chi2->cd(3); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
            drawComparisonPad(hdata_summary_mpgg2, hs_mpgg2, hu_mpgg2,
                              "(p_{target}+#gamma#gamma)^{2} [GeV^{2}]", Form("%s_mpgg2", file_tag.c_str()));
            title.DrawLatex(0.15, 0.93, Form("%s: M_{p#gamma#gamma}^{2}", iteration_label.c_str()));

            c_chi2->cd(4); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
            drawPullPad(hdata_summary_mpi0, hs_mgg, "M_{#gamma#gamma} [GeV/c^{2}]",
                        "M_{#gamma#gamma} pull", Form("%s_mgg_pull", file_tag.c_str()));

            c_chi2->cd(5); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
            drawPullPad(hdata_summary_mmiss, hs_mmiss, "M_{miss} [GeV/c^{2}]",
                        "M_{miss} pull", Form("%s_mmiss_pull", file_tag.c_str()));

            c_chi2->cd(6);
            gPad->SetLeftMargin(0.05); gPad->SetRightMargin(0.05);
            TPaveText *txt = new TPaveText(0.03, 0.03, 0.97, 0.97, "NDC");
            txt->SetBorderSize(1);
            txt->SetFillColor(candidate_accepted ? kGreen-9 : kRed-9);
            txt->SetTextAlign(12);
            txt->SetTextSize(0.030);
            txt->AddText(iteration_label.c_str());
            txt->AddText(Form("accepted=%s  best=%s", candidate_accepted ? "yes" : "no",
                              accepted_best ? "yes" : "no"));
            txt->AddText(Form("objective=%.6g", objective));
            txt->AddText(Form("chi2 Mgg=%.6g  Mmiss=%.6g  Mpgg2=%.6g",
                              breakdown.mpi0.chi2, breakdown.mmiss.chi2, breakdown.mpgg2.chi2));
            txt->AddText(Form("weights %.2f %.2f %.2f",
                              Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2));
            txt->AddText(Form("interp photon lookups=%d", interpolated_count));
            txt->AddText(Form("note: %s", sweep_note.c_str()));
            txt->AddText("Sim norm: unsmeared->data, same scale applied to smeared");
            txt->Draw();

            c_chi2->Print(pdf_file.c_str());
            string progress_file = progress_pdf_dir + "/" + file_tag + ".pdf";
            c_chi2->Print(progress_file.c_str());
            ++page_count;

            auto section_for_point = [&](double x, double y) {
                if (!in_geometry(x, y)) return -1;
                int best_is = -1;
                double best_d2 = std::numeric_limits<double>::max();
                for (int js = 0; js < nsec; ++js) {
                    if (!src_success[js]) continue;
                    if (!inSection(sections[js], x, y)) continue;
                    double cx = x_min + (sections[js].ix + 0.5) * base_wx;
                    double cy = y_min + (sections[js].iy + 0.5) * base_wy;
                    double dx = x - cx;
                    double dy = y - cy;
                    double d2 = dx * dx + dy * dy;
                    if (d2 < best_d2) {
                        best_d2 = d2;
                        best_is = js;
                    }
                }
                if (best_is >= 0) return best_is;
                int ix = (int)floor((x - x_min) / base_wx);
                int iy = (int)floor((y - y_min) / base_wy);
                ix = max(0, min(nx - 1, ix));
                iy = max(0, min(ny - 1, iy));
                int base_is = iy * nx + ix;
                return (base_is >= 0 && base_is < nsec && src_success[base_is]) ? base_is : -1;
            };

            auto attach_ext_params = [&](ClusterPair &ev) {
                auto set_ext = [&](int js, double &ma, double &mb, double &mc,
                                   double &sig, double &spos) {
                    if (js < 0 || js >= nsec || !src_success[js]) return;
                    ma = src_results[js].mu_a;
                    mb = src_results[js].mu;
                    mc = src_results[js].mu_c;
                    sig = src_results[js].sigma;
                    spos = src_results[js].sigma_pos;
                };
                if (!ev.photon1_in_section) {
                    set_ext(section_for_point(ev.x1, ev.y1),
                            ev.mu_a1_ext, ev.mu1_ext, ev.mu_c1_ext,
                            ev.sigma1_ext, ev.sigma_pos1_ext);
                }
                if (!ev.photon2_in_section) {
                    set_ext(section_for_point(ev.x2, ev.y2),
                            ev.mu_a2_ext, ev.mu2_ext, ev.mu_c2_ext,
                            ev.sigma2_ext, ev.sigma_pos2_ext);
                }
            };

            string section_progress_file = progress_pdf_dir + "/" + file_tag + "_sections.pdf";
            bool opened_section_progress = false;
            for (int is = 0; is < nsec; ++is) {
                if (!src_success[is]) continue;
                if (hdata_sec[is]->Integral() <= 0.0 ||
                    hdata_mmiss_sec[is]->Integral() <= 0.0 ||
                    hdata_mpgg2_sec[is]->Integral() <= 0.0) continue;

                vector<ClusterPair> section_events = sim_events_opt_per_section[is];
                for (auto &ev : section_events) attach_ext_params(ev);

                TH1D hu_sec_mgg(Form("hu_%s_%s_mgg", file_tag.c_str(), sections[is].name().c_str()),
                                ";M_{#gamma#gamma} [GeV/c^{2}];Counts",
                                Config::MGGAMMA_NBINS, Config::MGGAMMA_MIN, Config::MGGAMMA_MAX);
                TH1D hu_sec_mmiss(Form("hu_%s_%s_mmiss", file_tag.c_str(), sections[is].name().c_str()),
                                  ";M_{miss} [GeV/c^{2}];Counts",
                                  Config::MMISS_NBINS, Config::MMISS_MIN, Config::MMISS_MAX);
                TH1D hu_sec_mpgg2(Form("hu_%s_%s_mpgg2", file_tag.c_str(), sections[is].name().c_str()),
                                  ";(p_{target}+#gamma#gamma)^{2} [GeV^{2}];Counts",
                                  Config::MPGG2_NBINS, Config::MPGG2_MIN, Config::MPGG2_MAX);
                TH1D hs_sec_mgg(Form("hs_%s_%s_mgg", file_tag.c_str(), sections[is].name().c_str()),
                                ";M_{#gamma#gamma} [GeV/c^{2}];Counts",
                                Config::MGGAMMA_NBINS, Config::MGGAMMA_MIN, Config::MGGAMMA_MAX);
                TH1D hs_sec_mmiss(Form("hs_%s_%s_mmiss", file_tag.c_str(), sections[is].name().c_str()),
                                  ";M_{miss} [GeV/c^{2}];Counts",
                                  Config::MMISS_NBINS, Config::MMISS_MIN, Config::MMISS_MAX);
                TH1D hs_sec_mpgg2(Form("hs_%s_%s_mpgg2", file_tag.c_str(), sections[is].name().c_str()),
                                  ";(p_{target}+#gamma#gamma)^{2} [GeV^{2}];Counts",
                                  Config::MPGG2_NBINS, Config::MPGG2_MIN, Config::MPGG2_MAX);
                hu_sec_mgg.Sumw2(); hu_sec_mmiss.Sumw2(); hu_sec_mpgg2.Sumw2();
                hs_sec_mgg.Sumw2(); hs_sec_mmiss.Sumw2(); hs_sec_mpgg2.Sumw2();

                fillUnsmearedHistogramsForNormalization(section_events, hu_sec_mgg, hu_sec_mmiss, hu_sec_mpgg2,
                                                        src_results[is].p_e_scale);
                TRandom3 rng_section_plot(plot_seed + (unsigned int)is + 17u);
                fillSmearedHistogramsAtParams(section_events,
                                              hs_sec_mgg, hs_sec_mmiss, hs_sec_mpgg2,
                                              src_results[is].mu_a,
                                              src_results[is].mu,
                                              src_results[is].mu_c,
                                              src_results[is].sigma,
                                              src_results[is].sigma_pos,
                                              src_results[is].p_e_scale,
                                              rng_section_plot,
                                              optimizer_nsmear,
                                              src_results[is].res_A,
                                              src_results[is].res_B,
                                              src_results[is].res_C);
                scale_pair(hu_sec_mgg, hs_sec_mgg, *hdata_sec[is]);
                scale_pair(hu_sec_mmiss, hs_sec_mmiss, *hdata_mmiss_sec[is]);
                scale_pair(hu_sec_mpgg2, hs_sec_mpgg2, *hdata_mpgg2_sec[is]);

                c_chi2->Clear();
                c_chi2->Divide(3, 2);
                c_chi2->cd(1); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
                drawComparisonPad(*hdata_sec[is], hs_sec_mgg, hu_sec_mgg,
                                  "M_{#gamma#gamma} [GeV/c^{2}]",
                                  Form("%s_%s_mgg", file_tag.c_str(), sections[is].name().c_str()));
                TLatex st; st.SetNDC(); st.SetTextSize(0.043);
                st.DrawLatex(0.15, 0.93, Form("%s %s: M_{#gamma#gamma}",
                                               iteration_label.c_str(), sections[is].name().c_str()));

                c_chi2->cd(2); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
                drawComparisonPad(*hdata_mmiss_sec[is], hs_sec_mmiss, hu_sec_mmiss,
                                  "M_{miss} [GeV/c^{2}]",
                                  Form("%s_%s_mmiss", file_tag.c_str(), sections[is].name().c_str()));
                st.DrawLatex(0.15, 0.93, Form("%s %s: M_{miss}",
                                               iteration_label.c_str(), sections[is].name().c_str()));

                c_chi2->cd(3); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
                drawComparisonPad(*hdata_mpgg2_sec[is], hs_sec_mpgg2, hu_sec_mpgg2,
                                  "(p_{target}+#gamma#gamma)^{2} [GeV^{2}]",
                                  Form("%s_%s_mpgg2", file_tag.c_str(), sections[is].name().c_str()));
                st.DrawLatex(0.15, 0.93, Form("%s %s: M_{p#gamma#gamma}^{2}",
                                               iteration_label.c_str(), sections[is].name().c_str()));

                c_chi2->cd(4); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
                drawPullPad(*hdata_sec[is], hs_sec_mgg, "M_{#gamma#gamma} [GeV/c^{2}]",
                            "M_{#gamma#gamma} pull",
                            Form("%s_%s_mgg_pull", file_tag.c_str(), sections[is].name().c_str()));

                c_chi2->cd(5); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
                drawPullPad(*hdata_mmiss_sec[is], hs_sec_mmiss, "M_{miss} [GeV/c^{2}]",
                            "M_{miss} pull",
                            Form("%s_%s_mmiss_pull", file_tag.c_str(), sections[is].name().c_str()));

                c_chi2->cd(6);
                TPaveText *info = new TPaveText(0.03, 0.03, 0.97, 0.97, "NDC");
                info->SetBorderSize(1);
                info->SetFillColor(candidate_accepted ? kGreen-9 : kRed-9);
                info->SetTextAlign(12);
                info->SetTextSize(0.030);
                info->AddText(Form("%s section %s", iteration_label.c_str(), sections[is].name().c_str()));
                info->AddText(Form("a=%.5f b=%.5f c=%.5f",
                                   src_results[is].mu_a, src_results[is].mu, src_results[is].mu_c));
                info->AddText(Form("sigma=%.5f sigma_pos=%.5f",
                                   src_results[is].sigma, src_results[is].sigma_pos));
                info->AddText(Form("section chi2=%.6g chi2/ndf=%.6g",
                                   src_results[is].chi2, src_results[is].chi2_per_ndf()));
                info->AddText(Form("global objective=%.6g", objective));
                info->AddText(Form("accepted=%s best=%s", candidate_accepted ? "yes" : "no",
                                   accepted_best ? "yes" : "no"));
                info->AddText("Section plots use optimizer subset + optimizer Nsmear");
                info->Draw();

                c_chi2->Print(pdf_file.c_str());
                if (!opened_section_progress) {
                    c_chi2->Print((section_progress_file + "[").c_str());
                    opened_section_progress = true;
                }
                c_chi2->Print(section_progress_file.c_str());
                ++page_count;
            }
            if (opened_section_progress) {
                c_chi2->Print((section_progress_file + "]").c_str());
            }
        };

        // =========================================================================
        // === Per-section Sobol/MIGRAD fit used by every sweep
        // =========================================================================
        auto fit_section_fast_refine = [&](const vector<ClusterPair> &simEvents,
                                        const TH1D &hdata_mpi0,
                                        const TH1D &hdata_mmiss,
                                        const TH1D &hdata_mpgg2,
                                        TRandom3 &rng_local,
                                        int nsmear_local,
                                        double p_e_scale_local,
                                        const FitResult &seed,
                                        const string &progress_label,
                                        bool separate_sigma_pos_stage) {
            FitResult out = seed;
            out.scan_data = Chi2ScanData();
            const bool fit_p_e_local = Config::ENABLE_ELECTRON_MOMENTUM_SCALING && Config::fit_objective_uses_mmiss();
            double pe_seed = fit_p_e_local ? seed.p_e_scale : p_e_scale_local;
            if (!std::isfinite(pe_seed)) pe_seed = p_e_scale_local;
            pe_seed = min(max(pe_seed, Config::PE_SCALE_MIN), Config::PE_SCALE_MAX);
            out.p_e_scale = pe_seed;
            out.res_A = Config::RESOLUTION_A_DEFAULT;
            out.res_B = Config::RESOLUTION_B_DEFAULT;
            out.res_C = Config::RESOLUTION_C_DEFAULT;

            int n_params = 2 + (Config::ENABLE_POSITION_SMEARING ? 1 : 0) + (fit_p_e_local ? 1 : 0)
                            + (Config::ENABLE_ENERGY_DEPENDENT_MU ? 2 : 0);
            int n_bins_chi2 = 0;
            if (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_MPI0_ONLY) {
                n_bins_chi2 = countInformativeBinsData(hdata_mpi0);
            } else if (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_MMISS_ONLY) {
                n_bins_chi2 = countInformativeBinsData(hdata_mmiss);
            } else if (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_MPGG2_ONLY) {
                n_bins_chi2 = countInformativeBinsData(hdata_mpgg2);
            } else {
                n_bins_chi2 = countInformativeBinsData(hdata_mpi0)
                            + countInformativeBinsData(hdata_mmiss)
                            + countInformativeBinsData(hdata_mpgg2);
            }
            int n_norm_hists = (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_BOTH) ? 3 : 1;
            out.ndf = n_bins_chi2 - n_params - n_norm_hists;

            double mu_a0 = seed.mu_a;
            double mu0 = min(max(seed.mu, Config::MU_MIN), Config::MU_MAX);
            double mu_c0 = seed.mu_c;
            double sigma0 = min(max(seed.sigma, max(1e-6, Config::SIGMA_MIN)), Config::SIGMA_MAX);
            double sigma_pos0 = Config::ENABLE_POSITION_SMEARING
                                    ? min(max(seed.sigma_pos, Config::SIGMA_POS_MIN), Config::SIGMA_POS_MAX)
                                    : 0.0;
            double p_e0 = pe_seed;

            auto chi2_eval = [&](double mu_a, double mu_b, double mu_c, double sigma, double sigma_pos, double p_e_scale) {
                TRandom3 rng_eval(97531);
                return eval_chi2_selected(mu_a, mu_b, mu_c, sigma, sigma_pos, p_e_scale,
                                        simEvents,
                                        hdata_mpi0, hdata_mmiss, hdata_mpgg2,
                                        rng_eval, nsmear_local,
                                        Config::RESOLUTION_A_DEFAULT,
                                        Config::RESOLUTION_B_DEFAULT,
                                        Config::RESOLUTION_C_DEFAULT,
                                        Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);
            };

            ParamPoint center;
            center.mu_a = mu_a0;
            center.mu = mu0;
            center.mu_c = mu_c0;
            center.sigma = sigma0;
            center.sigma_pos = sigma_pos0;
            center.p_e_scale = p_e0;

            auto append_seed_debug = [&](const vector<SeedDiagnostic> &src) {
                out.seed_diagnostics.insert(out.seed_diagnostics.end(), src.begin(), src.end());
            };

            MigradRefineResult final_mr;
            FitFlags final_flags;
            if (separate_sigma_pos_stage && Config::ENABLE_POSITION_SMEARING) {
                FitFlags energy_flags;
                energy_flags.mu_a = Config::ENABLE_ENERGY_DEPENDENT_MU;
                energy_flags.mu = true;
                energy_flags.mu_c = Config::ENABLE_ENERGY_DEPENDENT_MU;
                energy_flags.sigma = true;
                energy_flags.sigma_pos = false;
                energy_flags.p_e_scale = fit_p_e_local;

                vector<SeedDiagnostic> energy_debug;
                MigradRefineResult energy_mr = run_global_multistart_refinement(
                    chi2_eval, center, energy_flags, &energy_debug,
                    progress_label.empty() ? "" : progress_label + " energy+sigma");
                append_seed_debug(energy_debug);

                ParamPoint pos_center;
                pos_center.mu_a = energy_mr.mu_a;
                pos_center.mu = energy_mr.mu;
                pos_center.mu_c = energy_mr.mu_c;
                pos_center.sigma = energy_mr.sigma;
                pos_center.sigma_pos = sigma_pos0;
                pos_center.p_e_scale = energy_mr.p_e_scale;

                FitFlags pos_flags;
                pos_flags.sigma_pos = true;

                vector<SeedDiagnostic> pos_debug;
                MigradRefineResult pos_mr = run_global_multistart_refinement(
                    chi2_eval, pos_center, pos_flags, &pos_debug,
                    progress_label.empty() ? "" : progress_label + " sigma_pos");
                append_seed_debug(pos_debug);

                final_mr = pos_mr;
                final_flags = pos_flags;
                out.optimizer_mode = "sobol_multistart_energy_sigma_then_sigma_pos";
                out.hesse_ok = energy_mr.hesse_ok && pos_mr.hesse_ok;
                out.max_abs_corr = std::max(energy_mr.max_abs_corr, pos_mr.max_abs_corr);
                out.max_corr_pair = (energy_mr.max_abs_corr >= pos_mr.max_abs_corr)
                    ? energy_mr.max_corr_pair
                    : pos_mr.max_corr_pair;
                out.profile_diagnostics = build_profile_diagnostics(chi2_eval, energy_mr, energy_flags);
            } else {
                FitFlags joint_flags;
                joint_flags.mu_a = Config::ENABLE_ENERGY_DEPENDENT_MU;
                joint_flags.mu = true;
                joint_flags.mu_c = Config::ENABLE_ENERGY_DEPENDENT_MU;
                joint_flags.sigma = true;
                joint_flags.sigma_pos = Config::ENABLE_POSITION_SMEARING;
                joint_flags.p_e_scale = fit_p_e_local;

                vector<SeedDiagnostic> joint_debug;
                final_mr = run_global_multistart_refinement(chi2_eval,
                                                            center,
                                                            joint_flags,
                                                            &joint_debug,
                                                            progress_label);
                append_seed_debug(joint_debug);
                final_flags = joint_flags;
                out.optimizer_mode = "sobol_multistart";
                out.hesse_ok = final_mr.hesse_ok;
                out.max_abs_corr = final_mr.max_abs_corr;
                out.max_corr_pair = final_mr.max_corr_pair;
                out.profile_diagnostics = build_profile_diagnostics(chi2_eval, final_mr, final_flags);
            }

            out.mu_a      = final_mr.mu_a;
            out.mu        = final_mr.mu;
            out.mu_c      = final_mr.mu_c;
            out.sigma     = final_mr.sigma;
            out.sigma_pos = Config::ENABLE_POSITION_SMEARING ? final_mr.sigma_pos : 0.0;
            out.p_e_scale = final_mr.p_e_scale;
            out.chi2      = final_mr.chi2;
            out.n_seed_trials = (Config::GLOBAL_MULTISTART_SEEDS + 1) *
                                ((separate_sigma_pos_stage && Config::ENABLE_POSITION_SMEARING) ? 2 : 1);
            out.n_migrad_trials = (int)out.seed_diagnostics.size();
            return out;
        };

        vector<bool> section_active(nsec, false);
        for (int is = 0; is < nsec; ++is) {
            size_t ndata_sec = (size_t)data_selected_count_per_section[is];
            size_t nsim_sec = (size_t)sim_selected_count_per_section[is];
            if (ndata_sec >= Config::MIN_EVENTS_PER_SECTION && nsim_sec >= Config::MIN_EVENTS_PER_SECTION) {
                section_active[is] = true;
            } else {
                fit_status[is] = "low_stats";
            }
        }

        auto make_nominal_fit_seed = [&]() {
            FitResult seed{};
            seed.mu_a = Config::ENABLE_ENERGY_DEPENDENT_MU ? Config::MU_ENERGY_A_INIT : 0.0;
            seed.mu = Config::MU_ENERGY_B_INIT;
            seed.mu_c = Config::ENABLE_ENERGY_DEPENDENT_MU ? Config::MU_ENERGY_C_INIT : 0.0;
            seed.sigma = 0.05;
            seed.sigma_pos = 0.0;
            seed.p_e_scale = global_p_e_scale;
            seed.chi2 = 1e300;
            seed.res_A = Config::RESOLUTION_A_DEFAULT;
            seed.res_B = Config::RESOLUTION_B_DEFAULT;
            seed.res_C = Config::RESOLUTION_C_DEFAULT;
            seed.ndf = 1;
            seed.scan_data = Chi2ScanData();
            return seed;
        };

        FitResult global_prefit_result = make_nominal_fit_seed();
        bool global_prefit_success = false;
        string global_prefit_status = Config::ENABLE_GLOBAL_PREFIT_SEED ? "not_run" : "disabled";

        auto make_default_fit_seed = [&]() {
            FitResult seed = global_prefit_success ? global_prefit_result : make_nominal_fit_seed();
            seed.chi2 = 1e300;
            seed.scan_data = Chi2ScanData();
            seed.seed_diagnostics.clear();
            seed.profile_diagnostics.clear();
            seed.optimizer_mode = global_prefit_success ? "global_prefit_seed" : "nominal_seed";
            seed.n_seed_trials = 0;
            seed.n_migrad_trials = 0;
            seed.hesse_ok = false;
            seed.max_abs_corr = 0.0;
            seed.max_corr_pair = "";
            return seed;
        };

        auto update_final_fit_status = [&]() {
            for (int is = 0; is < nsec; ++is) {
                if (!section_active[is]) continue;
                if (!fit_success[is]) {
                    fit_status[is] = "fit_failed";
                    continue;
                }
                double chi2_ndf = fit_results[is].chi2_per_ndf();
                bool good_fit = (chi2_ndf <= Config::MAX_CHI2_PER_NDF);
                if (Config::SKIP_BAD_FITS && !good_fit) {
                    fit_status[is] = "poor_fit_skipped";
                    fit_success[is] = false;
                } else {
                    fit_status[is] = good_fit ? "fit_ok" : "poor_fit";
                }
            }
        };

        auto file_exists = [](const string &path) {
            ifstream in(path.c_str());
            return in.good();
        };
        auto read_text_file = [](const string &path) {
            ifstream in(path.c_str(), std::ios::binary);
            ostringstream ss;
            ss << in.rdbuf();
            return ss.str();
        };
        auto split_csv_line = [](const string &line) {
            vector<string> out;
            string field;
            bool in_quotes = false;
            for (char ch : line) {
                if (ch == '"') {
                    in_quotes = !in_quotes;
                } else if (ch == ',' && !in_quotes) {
                    out.push_back(field);
                    field.clear();
                } else {
                    field.push_back(ch);
                }
            }
            out.push_back(field);
            return out;
        };
        auto trim_cr = [](string s) {
            while (!s.empty() && (s.back() == '\r' || s.back() == '\n')) s.pop_back();
            return s;
        };
        auto csv_header_map = [&](const vector<string> &header) {
            std::map<string, int> idx;
            for (int i = 0; i < (int)header.size(); ++i) idx[trim_cr(header[i])] = i;
            return idx;
        };
        auto csv_value = [&](const vector<string> &row, const std::map<string, int> &idx, const string &key) {
            auto it = idx.find(key);
            if (it == idx.end() || it->second < 0 || it->second >= (int)row.size()) return string("");
            return trim_cr(row[it->second]);
        };
        auto to_double_csv = [](const string &s, double fallback = 0.0) {
            if (s.empty()) return fallback;
            char *end = nullptr;
            double v = std::strtod(s.c_str(), &end);
            return (end != s.c_str() && std::isfinite(v)) ? v : fallback;
        };
        auto to_int_csv = [](const string &s, int fallback = 0) {
            if (s.empty()) return fallback;
            char *end = nullptr;
            long v = std::strtol(s.c_str(), &end, 10);
            return (end != s.c_str()) ? (int)v : fallback;
        };
        auto fill_fit_result_from_csv = [&](FitResult &fr,
                                            const vector<string> &row,
                                            const std::map<string, int> &idx) {
            fr.mu_a = to_double_csv(csv_value(row, idx, "mu_a"), fr.mu_a);
            fr.mu = to_double_csv(csv_value(row, idx, "mu_b"), fr.mu);
            fr.mu_c = to_double_csv(csv_value(row, idx, "mu_c"), fr.mu_c);
            fr.sigma = to_double_csv(csv_value(row, idx, "sigma"), fr.sigma);
            fr.sigma_pos = to_double_csv(csv_value(row, idx, "sigma_pos"), fr.sigma_pos);
            fr.p_e_scale = to_double_csv(csv_value(row, idx, "p_e_scale"), fr.p_e_scale);
            fr.chi2 = to_double_csv(csv_value(row, idx, "chi2"), fr.chi2);
            double chi2_ndf = to_double_csv(csv_value(row, idx, "chi2_ndf"), 0.0);
            if (chi2_ndf > 0.0 && std::isfinite(fr.chi2) && fr.chi2 > 0.0) {
                fr.ndf = max(1, (int)std::llround(fr.chi2 / chi2_ndf));
            } else {
                fr.ndf = max(1, fr.ndf);
            }
            fr.res_A = Config::RESOLUTION_A_DEFAULT;
            fr.res_B = Config::RESOLUTION_B_DEFAULT;
            fr.res_C = Config::RESOLUTION_C_DEFAULT;
        };

        bool cache_forced_off = false;
        const char *force_reopt_env = std::getenv("NPS_SMEARING_FORCE_REOPT");
        if (force_reopt_env && string(force_reopt_env) == "1") cache_forced_off = true;
        bool cache_fingerprint_ok = false;
        if (!cache_forced_off && file_exists(cache_fingerprint_file)) {
            cache_fingerprint_ok = (read_text_file(cache_fingerprint_file) == current_cache_fingerprint);
        }

        if (cache_fingerprint_ok &&
            file_exists(sweep_history_csv_file) &&
            file_exists(optimizer_summary_csv_file)) {
            cout << "\n==== Cached smearing replay: config/source fingerprint unchanged ====\n"
                 << "Using previous optimizer CSVs to rebuild global + sweep diagnostic plots.\n"
                 << "Set NPS_SMEARING_FORCE_REOPT=1 to force full Sobol/MIGRAD rerun.\n";

            bool cache_loaded_any = false;
            {
                ifstream opt_in(optimizer_summary_csv_file.c_str());
                string header_line;
                if (std::getline(opt_in, header_line)) {
                    vector<string> header = split_csv_line(header_line);
                    auto idx = csv_header_map(header);
                    string line;
                    while (std::getline(opt_in, line)) {
                        if (line.empty()) continue;
                        vector<string> row = split_csv_line(line);
                        if (csv_value(row, idx, "section") != "global_prefit") continue;
                        global_prefit_result = make_nominal_fit_seed();
                        fill_fit_result_from_csv(global_prefit_result, row, idx);
                        global_prefit_result.optimizer_mode = csv_value(row, idx, "optimizer_mode");
                        global_prefit_result.n_seed_trials = to_int_csv(csv_value(row, idx, "n_seed_trials"), 0);
                        global_prefit_result.n_migrad_trials = to_int_csv(csv_value(row, idx, "n_migrad_trials"), 0);
                        global_prefit_result.hesse_ok = (to_int_csv(csv_value(row, idx, "hesse_ok"), 0) != 0);
                        global_prefit_result.max_abs_corr = to_double_csv(csv_value(row, idx, "max_abs_corr"), 0.0);
                        global_prefit_result.max_corr_pair = csv_value(row, idx, "max_corr_pair");
                        global_prefit_status = csv_value(row, idx, "fit_status");
                        global_prefit_success = (global_prefit_status == "fit_ok" &&
                                                 std::isfinite(global_prefit_result.chi2) &&
                                                 global_prefit_result.chi2 < 1e299);
                        if (global_prefit_success) {
                            global_p_e_scale = global_prefit_result.p_e_scale;
                            vector<FitResult> global_prefit_map(nsec, global_prefit_result);
                            vector<bool> global_prefit_map_success(nsec, true);
                            ObjectiveBreakdown global_prefit_breakdown;
                            int global_prefit_interp = 0;
                            double global_prefit_objective = evaluate_global_objective(global_prefit_map,
                                                                                       global_prefit_map_success,
                                                                                       optimizer_nsmear,
                                                                                       &global_prefit_breakdown,
                                                                                       &global_prefit_interp);
                            plot_iteration_candidate_histograms("Cached global prefit",
                                                                "cached_global_prefit",
                                                                global_prefit_map,
                                                                global_prefit_map_success,
                                                                global_prefit_breakdown,
                                                                global_prefit_objective,
                                                                true,
                                                                true,
                                                                "cached_global_prefit",
                                                                global_prefit_interp);
                            cache_loaded_any = true;
                        }
                        break;
                    }
                }
            }

            struct CachedSweepRow {
                int sweep = 0;
                vector<FitResult> results;
                vector<bool> success;
                ObjectiveBreakdown breakdown;
                double objective = 1e300;
                bool candidate_accepted = false;
                bool accepted_best = false;
                string sweep_note = "cached";
            };
            std::map<int, CachedSweepRow> cached_sweeps;
            {
                ifstream sweep_in(sweep_history_csv_file.c_str());
                string header_line;
                if (std::getline(sweep_in, header_line)) {
                    vector<string> header = split_csv_line(header_line);
                    auto idx = csv_header_map(header);
                    string line;
                    while (std::getline(sweep_in, line)) {
                        if (line.empty()) continue;
                        vector<string> row = split_csv_line(line);
                        int sweep_id = to_int_csv(csv_value(row, idx, "sweep"), 0);
                        int ix = to_int_csv(csv_value(row, idx, "ix"), -999);
                        int iy = to_int_csv(csv_value(row, idx, "iy"), -999);
                        if (sweep_id <= 0 || ix < 0 || ix >= nx || iy < 0 || iy >= ny) continue;
                        int is = iy * nx + ix;
                        CachedSweepRow &cs = cached_sweeps[sweep_id];
                        if (cs.results.empty()) {
                            cs.sweep = sweep_id;
                            cs.results.assign(nsec, make_nominal_fit_seed());
                            cs.success.assign(nsec, false);
                            cs.objective = to_double_csv(csv_value(row, idx, "candidate_objective"), 1e300);
                            cs.candidate_accepted = (to_int_csv(csv_value(row, idx, "candidate_accepted"), 0) != 0);
                            cs.accepted_best = (to_int_csv(csv_value(row, idx, "accepted_best"), 0) != 0);
                            cs.sweep_note = csv_value(row, idx, "sweep_note");
                            cs.breakdown.mpi0.chi2 = to_double_csv(csv_value(row, idx, "global_chi2_mpi0"), 1e300);
                            cs.breakdown.mmiss.chi2 = to_double_csv(csv_value(row, idx, "global_chi2_mmiss"), 1e300);
                            cs.breakdown.mpgg2.chi2 = to_double_csv(csv_value(row, idx, "global_chi2_mpgg2"), 1e300);
                            cs.breakdown.mpi0.empty_sim_data_positive_bins = to_int_csv(csv_value(row, idx, "global_empty_sim_mpi0"), 0);
                            cs.breakdown.mmiss.empty_sim_data_positive_bins = to_int_csv(csv_value(row, idx, "global_empty_sim_mmiss"), 0);
                            cs.breakdown.mpgg2.empty_sim_data_positive_bins = to_int_csv(csv_value(row, idx, "global_empty_sim_mpgg2"), 0);
                        }
                        bool row_success = (to_int_csv(csv_value(row, idx, "fit_success"), 0) != 0);
                        if (!row_success) continue;
                        fill_fit_result_from_csv(cs.results[is], row, idx);
                        cs.results[is].optimizer_mode = "cached_sweep_replay";
                        cs.success[is] = true;
                    }
                }
            }

            for (const auto &kv : cached_sweeps) {
                const CachedSweepRow &cs = kv.second;
                plot_iteration_candidate_histograms(Form("Cached coupled sweep %d", cs.sweep),
                                                    Form("cached_sweep_%03d_candidate", cs.sweep),
                                                    cs.results,
                                                    cs.success,
                                                    cs.breakdown,
                                                    cs.objective,
                                                    cs.candidate_accepted,
                                                    cs.accepted_best,
                                                    cs.sweep_note.empty() ? "cached" : cs.sweep_note,
                                                    0);
                cache_loaded_any = true;
            }

            if (cache_loaded_any) {
                string pdf_close = pdf_file + "]";
                c_chi2->Print(pdf_close.c_str());
                delete c_chi2;
                fout.Close();
                copyFileIfDifferent(pdf_file, canonical_pdf_file, "latest chi2 PDF");
                copyFileIfDifferent(out_file, timestamped_out_file, "fitter ROOT output");
                writeSmearingManifest(metadata_manifest_file,
                                      run_tag,
                                      created_at_local,
                                      data_file,
                                      sim_file,
                                      out_file,
                                      timestamped_out_file,
                                      interp_file,
                                      timestamped_interp_file,
                                      pdf_file,
                                      canonical_pdf_file,
                                      progress_pdf_dir,
                                      current_cache_fingerprint);
                cout << "Cached diagnostic plots saved to " << pdf_file << "\n"
                     << "Progress PDFs in " << progress_pdf_dir << "\n";
                return 0;
            }

            cout << "Cache fingerprint matched, but no replayable sweep/global rows found. Running optimizer afresh.\n";
        } else if (cache_forced_off) {
            cout << "\n==== Cache replay disabled by NPS_SMEARING_FORCE_REOPT=1; running optimizer afresh ====\n";
        } else if (file_exists(cache_fingerprint_file)) {
            cout << "\n==== Cache fingerprint changed; running optimizer afresh ====\n";
        } else {
            cout << "\n==== No smearing cache fingerprint found; running optimizer afresh ====\n";
        }

        ofstream sweep_history(sweep_history_csv_file.c_str());
        sweep_history << "sweep,ix,iy,section,active,fit_success,candidate_accepted,accepted_best,"
                      << "candidate_objective,previous_accepted_objective,current_accepted_objective,best_ever_objective,"
                      << "global_chi2_mpi0,global_chi2_mmiss,global_chi2_mpgg2,"
                      << "global_empty_sim_mpi0,global_empty_sim_mmiss,global_empty_sim_mpgg2,"
                      << "sweep_note,accept_reason,repeated_candidate,"
                      << "norm_vs_previous_candidate,norm_vs_accepted,norm_vs_best,"
                      << "consecutive_rejected,consecutive_repeated_rejected,stop_after_sweep,stop_reason,"
                      << "runtime_sec,strategy,"
                      << "chi2,chi2_ndf,delta_mu_b,delta_sigma,delta_sigma_pos,"
                      << "mu_a,mu_b,mu_c,sigma,sigma_pos,p_e_scale\n";

        ObjectiveBreakdown global_prefit_history_breakdown;
        int global_prefit_history_interp = 0;
        double global_prefit_history_objective = 1e300;

        if (Config::ENABLE_GLOBAL_PREFIT_SEED) {
            vector<ClusterPair> sim_events_global_prefit = sim_events_summary;
            int sim_global_prefit_selected = 0;
            for (auto &ev : sim_events_global_prefit) {
                ev.photon1_in_section = true;
                ev.photon2_in_section = true;
                ev.mu_a1_ext = 0.0; ev.mu1_ext = 1.0; ev.mu_c1_ext = 0.0;
                ev.sigma1_ext = 0.0; ev.sigma_pos1_ext = 0.0;
                ev.mu_a2_ext = 0.0; ev.mu2_ext = 1.0; ev.mu_c2_ext = 0.0;
                ev.sigma2_ext = 0.0; ev.sigma_pos2_ext = 0.0;
                if (ev.weight > 0.0) ++sim_global_prefit_selected;
            }

            int data_global_prefit_selected = 0;
            for (const auto &ev : data_events_summary) {
                if (ev.weight > 0.0) ++data_global_prefit_selected;
            }

            if (data_global_prefit_selected >= Config::MIN_EVENTS_PER_SECTION &&
                sim_global_prefit_selected >= Config::MIN_EVENTS_PER_SECTION &&
                hdata_summary_mpi0.Integral() > 0.0 &&
                !sim_events_global_prefit.empty()) {
                vector<ClusterPair> sim_events_global_prefit_opt = makeOptimizationSubset(
                    sim_events_global_prefit,
                    Config::OPT_MAX_SIM_EVENTS_GLOBAL_PREFIT,
                    Config::OPT_SUBSET_MGG_BINS);
                cout << "\n==== GLOBAL PREFIT: all-calorimeter response seed ====\n"
                     << "Using " << data_global_prefit_selected << " selected data events and "
                     << sim_global_prefit_selected << " selected sim events"
                     << " (optimizer subset=" << sim_events_global_prefit_opt.size()
                     << ", optimizer Nsmear=" << optimizer_nsmear
                     << "). Path: Sobol -> keep best N -> MIGRAD -> HESSE/profile.\n";

                TRandom3 rng_global_prefit(24681357);
                FitResult nominal_seed = make_nominal_fit_seed();
                global_prefit_result = fit_section_fast_refine(sim_events_global_prefit_opt,
                                                               hdata_summary_mpi0,
                                                               hdata_summary_mmiss,
                                                               hdata_summary_mpgg2,
                                                               rng_global_prefit,
                                                               optimizer_nsmear,
                                                               global_p_e_scale,
                                                               nominal_seed,
                                                               "Global prefit",
                                                               false);
                global_prefit_result.optimizer_mode = "global_prefit_sobol_multistart";
                global_prefit_success = std::isfinite(global_prefit_result.chi2) &&
                                        global_prefit_result.chi2 < 1e299;
                global_prefit_status = global_prefit_success ? "fit_ok" : "fit_failed";
                if (global_prefit_success) {
                    global_p_e_scale = global_prefit_result.p_e_scale;
                    cout << "Global prefit -> a=" << global_prefit_result.mu_a
                         << " b=" << global_prefit_result.mu
                         << " c=" << global_prefit_result.mu_c
                         << " sigma=" << global_prefit_result.sigma;
                    if (Config::ENABLE_POSITION_SMEARING) {
                        cout << " sigma_pos=" << global_prefit_result.sigma_pos;
                    }
                    if (Config::ENABLE_ELECTRON_MOMENTUM_SCALING && Config::fit_objective_uses_mmiss()) {
                        cout << " p_e_scale=" << global_prefit_result.p_e_scale;
                    }
                    cout << " chi2/ndf=" << global_prefit_result.chi2_per_ndf() << "\n";
                    fout.cd();
                    TNamed("global_prefit_status", global_prefit_status.c_str()).Write();
                    TNamed("global_prefit_mu_a", Form("%.17g", global_prefit_result.mu_a)).Write();
                    TNamed("global_prefit_mu_b", Form("%.17g", global_prefit_result.mu)).Write();
                    TNamed("global_prefit_mu_c", Form("%.17g", global_prefit_result.mu_c)).Write();
                    TNamed("global_prefit_sigma", Form("%.17g", global_prefit_result.sigma)).Write();
                    TNamed("global_prefit_sigma_pos", Form("%.17g", global_prefit_result.sigma_pos)).Write();
                    TNamed("global_prefit_p_e_scale", Form("%.17g", global_prefit_result.p_e_scale)).Write();
                    TNamed("global_prefit_chi2_ndf", Form("%.17g", global_prefit_result.chi2_per_ndf())).Write();

                    vector<FitResult> global_prefit_map(nsec, global_prefit_result);
                    vector<bool> global_prefit_map_success(nsec, true);
                    global_prefit_history_objective = evaluate_global_objective(global_prefit_map,
                                                                               global_prefit_map_success,
                                                                               optimizer_nsmear,
                                                                               &global_prefit_history_breakdown,
                                                                               &global_prefit_history_interp);
                    plot_iteration_candidate_histograms("Global prefit result",
                                                        "global_prefit",
                                                        global_prefit_map,
                                                        global_prefit_map_success,
                                                        global_prefit_history_breakdown,
                                                        global_prefit_history_objective,
                                                        true,
                                                        true,
                                                        "global_prefit",
                                                        global_prefit_history_interp);
                }
            } else {
                global_prefit_status = "low_stats";
                cout << "\n==== GLOBAL PREFIT: skipped (low statistics) ====\n"
                     << "selected data=" << data_global_prefit_selected
                     << " selected sim=" << sim_global_prefit_selected
                     << " minimum=" << Config::MIN_EVENTS_PER_SECTION << "\n";
            }
        }

        {
            const FitResult &fr = global_prefit_result;
            sweep_history << 0 << ","
                          << -1 << "," << -1 << ",global_prefit,"
                          << (Config::ENABLE_GLOBAL_PREFIT_SEED ? 1 : 0) << ","
                          << (global_prefit_success ? 1 : 0) << ","
                          << (global_prefit_success ? 1 : 0) << ","
                          << (global_prefit_success ? 1 : 0) << ","
                          << global_prefit_history_objective << ","
                          << 1e300 << ","
                          << global_prefit_history_objective << ","
                          << global_prefit_history_objective << ","
                          << global_prefit_history_breakdown.mpi0.chi2 << ","
                          << global_prefit_history_breakdown.mmiss.chi2 << ","
                          << global_prefit_history_breakdown.mpgg2.chi2 << ","
                          << global_prefit_history_breakdown.mpi0.empty_sim_data_positive_bins << ","
                          << global_prefit_history_breakdown.mmiss.empty_sim_data_positive_bins << ","
                          << global_prefit_history_breakdown.mpgg2.empty_sim_data_positive_bins << ","
                          << global_prefit_status << ","
                          << "global_prefit" << ","
                          << 0 << ","
                          << 0.0 << "," << 0.0 << "," << 0.0 << ","
                          << 0 << "," << 0 << "," << 0 << ",,"
                          << 0.0 << ",global_prefit,"
                          << fr.chi2 << "," << fr.chi2_per_ndf() << ","
                          << 0.0 << "," << 0.0 << "," << 0.0 << ","
                          << fr.mu_a << "," << fr.mu << "," << fr.mu_c << ","
                          << fr.sigma << "," << fr.sigma_pos << "," << fr.p_e_scale << "\n";
        }

        FitResult sweep0_ext_seed = global_prefit_success ? global_prefit_result : make_nominal_fit_seed();
        int seed_ext_assignments = reset_cross_boundary_ext_to_seed(sweep0_ext_seed);
        cout << "\n==== Coupled sweep 0: "
             << (global_prefit_success ? "global prefit external photon response" : "nominal external photon response")
             << " ====\n"
             << "Set " << seed_ext_assignments
             << " out-of-section photon coefficient assignments to a=" << sweep0_ext_seed.mu_a
             << ",b=" << sweep0_ext_seed.mu
             << ",c=" << sweep0_ext_seed.mu_c
             << ",sigma=" << sweep0_ext_seed.sigma
             << ",sigma_pos=" << sweep0_ext_seed.sigma_pos << ".\n";

        const int max_sweeps = max(1, Config::ITERATIVE_SECTION_SWEEPS);
        vector<FitResult> prev_results(nsec);
        vector<bool> prev_success(nsec, false);
        vector<FitResult> accepted_results(nsec);
        vector<bool> accepted_success(nsec, false);
        double accepted_global_objective = 1e300;
        int accepted_sweep_index = -1;
        vector<FitResult> best_sweep_results(nsec);
        vector<bool> best_sweep_success(nsec, false);
        double best_sweep_global_objective = 1e300;
        int best_sweep_index = -1;
        vector<double> global_objective_history;
        vector<FitResult> previous_candidate_results(nsec);
        vector<bool> previous_candidate_success(nsec, false);
        bool have_previous_candidate = false;
        int consecutive_rejected_sweeps = 0;
        int consecutive_repeated_rejected_sweeps = 0;
        int consecutive_cycle_rejected_sweeps = 0;

        auto sweep_improves = [&](double candidate, double reference) {
            if (!std::isfinite(candidate)) return false;
            if (!std::isfinite(reference) || reference >= 1e299) return true;
            const double tol = std::max(Config::COUPLED_ACCEPT_ABS_TOL,
                                        Config::COUPLED_ACCEPT_REL_TOL * std::max(std::fabs(reference), 1.0));
            return candidate < reference - tol;
        };

        auto parameter_norm_between = [&](const vector<FitResult> &a, const vector<bool> &a_ok,
                                          const vector<FitResult> &b, const vector<bool> &b_ok) {
            double sum2 = 0.0;
            int n = 0;
            for (int is = 0; is < nsec; ++is) {
                if (!(section_active[is] && a_ok[is] && b_ok[is])) continue;
                const double vals_a[6] = {a[is].mu_a, a[is].mu, a[is].mu_c,
                                          a[is].sigma, a[is].sigma_pos, a[is].p_e_scale};
                const double vals_b[6] = {b[is].mu_a, b[is].mu, b[is].mu_c,
                                          b[is].sigma, b[is].sigma_pos, b[is].p_e_scale};
                for (int ip = 0; ip < 6; ++ip) {
                    if (!std::isfinite(vals_a[ip]) || !std::isfinite(vals_b[ip])) continue;
                    double d = vals_a[ip] - vals_b[ip];
                    sum2 += d * d;
                    ++n;
                }
            }
            return (n > 0) ? std::sqrt(sum2 / n) : 0.0;
        };

        cout << "\n==== Iterative coupled section sweeps (parallel within sweep, barrier between sweeps) ====\n";
        for (int sweep = 0; sweep < max_sweeps; ++sweep) {
            cout << "\n-- Coupled sweep " << (sweep + 1) << " / " << max_sweeps << " --" << endl;
            auto sweep_t0 = std::chrono::steady_clock::now();

            vector<FitResult> sweep_results = fit_results;
            vector<bool> sweep_success = fit_success;

            #pragma omp parallel if(Config::ENABLE_PARALLEL_SECTION_FITS)
            {
                int thread_id = omp_get_thread_num();
                TRandom3 thread_rng((sweep + 1) * 1000003 + thread_id * 123456 + time(NULL));

                #pragma omp for schedule(dynamic)
                for (int is = 0; is < nsec; ++is) {
                    if (!section_active[is]) continue;

                    FitResult seed = fit_success[is] ? fit_results[is] : make_default_fit_seed();
                    FitResult fres = fit_section_fast_refine(sim_events_opt_per_section[is],
                                                            *hdata_sec[is], *hdata_mmiss_sec[is], *hdata_mpgg2_sec[is],
                                                            thread_rng, optimizer_nsmear,
                                                            seed.p_e_scale,
                                                            seed,
                                                            "",
                                                            true);

                    double chi2_ndf = fres.chi2_per_ndf();
                    #pragma omp critical(console)
                    {
                        cout << "[sweep " << (sweep + 1) << "] " << sections[is].name()
                            << " -> a=" << fres.mu_a
                            << " b=" << fres.mu
                            << " c=" << fres.mu_c
                            << " sigma=" << fres.sigma;
                        if (Config::ENABLE_POSITION_SMEARING) {
                            cout << " sigma_pos=" << std::fixed << std::setprecision(5) << fres.sigma_pos << " cm";
                        }
                        if (Config::ENABLE_ELECTRON_MOMENTUM_SCALING && Config::fit_objective_uses_mmiss()) {
                            cout << " p_e_scale=" << fres.p_e_scale;
                        }
                        cout << " chi2/ndf=" << chi2_ndf;
                        if (chi2_ndf > Config::MAX_CHI2_PER_NDF) cout << " [POOR]";
                        cout << "\n";
                    }

                    sweep_results[is] = fres;
                    sweep_success[is] = true;
                }
            }

            double max_dmu = 0.0, max_dsig = 0.0, max_dpos = 0.0;
            int ncomp = 0;
            for (int is = 0; is < nsec; ++is) {
                if (!(section_active[is] && sweep_success[is] && accepted_success[is])) continue;
                max_dmu = max(max_dmu, fabs(sweep_results[is].mu - accepted_results[is].mu));
                max_dsig = max(max_dsig, fabs(sweep_results[is].sigma - accepted_results[is].sigma));
                max_dpos = max(max_dpos, fabs(sweep_results[is].sigma_pos - accepted_results[is].sigma_pos));
                ++ncomp;
            }

            auto refresh_counts = refresh_cross_boundary_ext(sweep_results, sweep_success);

            ObjectiveBreakdown sweep_breakdown;
            int sweep_interp_count = 0;
            double sweep_global_objective = evaluate_global_objective(sweep_results, sweep_success,
                                                                      optimizer_nsmear,
                                                                      &sweep_breakdown,
                                                                      &sweep_interp_count);
            const double previous_sweep_objective = global_objective_history.empty()
                                                ? 1e300
                                                : global_objective_history.back();
            const double previous_accepted_objective = accepted_global_objective;
            vector<FitResult> previous_accepted_results = accepted_results;
            vector<bool> previous_accepted_success = accepted_success;
            bool candidate_accepted = use_accepted_state_rollback
                ? sweep_improves(sweep_global_objective, accepted_global_objective)
                : true;
            bool accepted_best = sweep_improves(sweep_global_objective, best_sweep_global_objective);
            if (accepted_best) {
                best_sweep_global_objective = sweep_global_objective;
                best_sweep_index = sweep;
                best_sweep_results = sweep_results;
                best_sweep_success = sweep_success;
            }
            string sweep_note = accepted_best ? "new_best" : "not_best";
            string accept_reason = use_accepted_state_rollback
                ? (candidate_accepted ? "improved_accepted_objective" : "rejected_no_accepted_improvement")
                : "auto_accept_completed_sweep";
            if (std::isfinite(previous_sweep_objective) && sweep_global_objective > previous_sweep_objective) {
                sweep_note += "|worse_than_previous";
            }
            bool possible_two_sweep_cycle = false;
            if (global_objective_history.size() >= 2) {
                double two_back = global_objective_history[global_objective_history.size() - 2];
                double cycle_tol = std::max(1e-9, 1e-4 * std::max(std::fabs(two_back), 1.0));
                if (std::fabs(sweep_global_objective - two_back) <= cycle_tol &&
                    sweep_global_objective > best_sweep_global_objective + cycle_tol) {
                    possible_two_sweep_cycle = true;
                    sweep_note += "|possible_two_sweep_cycle";
                }
            }
            double norm_vs_previous_candidate = have_previous_candidate
                ? parameter_norm_between(sweep_results, sweep_success,
                                         previous_candidate_results, previous_candidate_success)
                : 0.0;
            double norm_vs_accepted = parameter_norm_between(sweep_results, sweep_success,
                                                             accepted_results, accepted_success);
            double norm_vs_best = parameter_norm_between(sweep_results, sweep_success,
                                                         best_sweep_results, best_sweep_success);
            bool repeated_candidate = have_previous_candidate &&
                                      norm_vs_previous_candidate <= Config::COUPLED_REPEAT_NORM_TOL;
            if (repeated_candidate) sweep_note += "|repeated_candidate";
            global_objective_history.push_back(sweep_global_objective);

            if (candidate_accepted) {
                consecutive_rejected_sweeps = 0;
                consecutive_repeated_rejected_sweeps = 0;
                consecutive_cycle_rejected_sweeps = 0;
            } else {
                ++consecutive_rejected_sweeps;
                if (repeated_candidate) {
                    ++consecutive_repeated_rejected_sweeps;
                } else {
                    consecutive_repeated_rejected_sweeps = 0;
                }
                if (possible_two_sweep_cycle) {
                    ++consecutive_cycle_rejected_sweeps;
                } else {
                    consecutive_cycle_rejected_sweeps = 0;
                }
            }

            bool stop_after_sweep = false;
            string stop_reason = "";
            if (use_accepted_state_rollback && Config::ENABLE_COUPLED_REJECTED_REPEAT_STOP) {
                if (!candidate_accepted &&
                    repeated_candidate &&
                    consecutive_repeated_rejected_sweeps >= Config::COUPLED_REJECTED_REPEAT_PATIENCE) {
                    stop_after_sweep = true;
                    stop_reason = "repeated_rejected_candidate";
                } else if (!candidate_accepted &&
                           possible_two_sweep_cycle &&
                           consecutive_cycle_rejected_sweeps >= Config::COUPLED_REJECTED_CYCLE_PATIENCE) {
                    stop_after_sweep = true;
                    stop_reason = "rejected_two_sweep_cycle";
                }
            }
            if (stop_after_sweep) {
                sweep_note += "|early_stop";
            }

            plot_iteration_candidate_histograms(Form("Coupled sweep %d candidate", sweep + 1),
                                                Form("sweep_%03d_candidate", sweep + 1),
                                                sweep_results,
                                                sweep_success,
                                                sweep_breakdown,
                                                sweep_global_objective,
                                                candidate_accepted,
                                                accepted_best,
                                                sweep_note,
                                                sweep_interp_count);

            if (candidate_accepted) {
                accepted_results = sweep_results;
                accepted_success = sweep_success;
                accepted_global_objective = sweep_global_objective;
                accepted_sweep_index = sweep;
                fit_results = accepted_results;
                fit_success = accepted_success;
                refresh_counts = refresh_cross_boundary_ext(fit_results, fit_success);
            } else {
                fit_results = accepted_results;
                fit_success = accepted_success;
                if (accepted_sweep_index >= 0) {
                    refresh_counts = refresh_cross_boundary_ext(fit_results, fit_success);
                } else {
                    refresh_counts = refresh_cross_boundary_ext(sweep_results, sweep_success);
                }
            }

            previous_candidate_results = sweep_results;
            previous_candidate_success = sweep_success;
            have_previous_candidate = true;
            auto sweep_t1 = std::chrono::steady_clock::now();
            double runtime_sec = std::chrono::duration<double>(sweep_t1 - sweep_t0).count();

            for (int is = 0; is < nsec; ++is) {
                double dmu = 0.0, dsig = 0.0, dpos = 0.0;
                if (previous_accepted_success[is] && sweep_success[is]) {
                    dmu = sweep_results[is].mu - previous_accepted_results[is].mu;
                    dsig = sweep_results[is].sigma - previous_accepted_results[is].sigma;
                    dpos = sweep_results[is].sigma_pos - previous_accepted_results[is].sigma_pos;
                }
                const FitResult &fr = sweep_results[is];
                sweep_history << (sweep + 1) << ","
                              << sections[is].ix << "," << sections[is].iy << ","
                              << sections[is].name() << ","
                              << (section_active[is] ? 1 : 0) << ","
                              << (sweep_success[is] ? 1 : 0) << ","
                              << (candidate_accepted ? 1 : 0) << ","
                              << (accepted_best ? 1 : 0) << ","
                              << sweep_global_objective << ","
                              << previous_accepted_objective << ","
                              << accepted_global_objective << ","
                              << best_sweep_global_objective << ","
                              << sweep_breakdown.mpi0.chi2 << ","
                              << sweep_breakdown.mmiss.chi2 << ","
                              << sweep_breakdown.mpgg2.chi2 << ","
                              << sweep_breakdown.mpi0.empty_sim_data_positive_bins << ","
                              << sweep_breakdown.mmiss.empty_sim_data_positive_bins << ","
                              << sweep_breakdown.mpgg2.empty_sim_data_positive_bins << ","
                              << sweep_note << ","
                              << accept_reason << ","
                              << (repeated_candidate ? 1 : 0) << ","
                              << norm_vs_previous_candidate << ","
                              << norm_vs_accepted << ","
                              << norm_vs_best << ","
                              << consecutive_rejected_sweeps << ","
                              << consecutive_repeated_rejected_sweeps << ","
                              << (stop_after_sweep ? 1 : 0) << ","
                              << stop_reason << ","
                              << runtime_sec << ","
                              << sweep_acceptance_strategy << ","
                              << fr.chi2 << "," << fr.chi2_per_ndf() << ","
                              << dmu << "," << dsig << "," << dpos << ","
                              << fr.mu_a << "," << fr.mu << "," << fr.mu_c << ","
                              << fr.sigma << "," << fr.sigma_pos << "," << fr.p_e_scale << "\n";
            }

            cout << "Sweep " << (sweep + 1) << " complete barrier: refreshed "
                 << refresh_counts.first << " out-of-section photon assignments";
            if (refresh_counts.second > 0) cout << " (interpolation fallback=" << refresh_counts.second << ")";
            cout << "." << endl;
            cout << "Sweep " << (sweep + 1) << " global objective=" << sweep_global_objective
                 << " [Mgg=" << sweep_breakdown.mpi0.chi2
                 << ", Mmiss=" << sweep_breakdown.mmiss.chi2
                 << ", Mpgg2=" << sweep_breakdown.mpgg2.chi2
                 << ", interp=" << sweep_interp_count
                 << "] " << sweep_note
                 << " accepted=" << (candidate_accepted ? "yes" : "no")
                 << " accepted_obj=" << accepted_global_objective
                 << " best_obj=" << best_sweep_global_objective
                 << " repeat=" << (repeated_candidate ? "yes" : "no")
                 << " stop=" << (stop_after_sweep ? stop_reason : "no")
                 << " runtime_s=" << runtime_sec << endl;

            if (ncomp > 0) {
                cout << "Sweep " << (sweep + 1)
                     << " candidate-vs-accepted deltas: dmu=" << max_dmu
                     << " dsigma=" << max_dsig
                     << " dsigma_pos=" << max_dpos << endl;
            }

            prev_results = accepted_results;
            prev_success = accepted_success;

            bool converged = (candidate_accepted && ncomp > 0 &&
                              max_dmu <= Config::COUPLED_CONV_MU &&
                              max_dsig <= Config::COUPLED_CONV_SIGMA &&
                              max_dpos <= Config::COUPLED_CONV_SIGMA_POS);
            if (Config::ENABLE_COUPLED_SWEEP_CONVERGENCE_STOP && converged) {
                cout << "Coupled sweeps converged after sweep " << (sweep + 1) << "." << endl;
                break;
            }
            if (stop_after_sweep) {
                cout << "Coupled sweeps stopped after sweep " << (sweep + 1)
                     << " because " << stop_reason
                     << ". Keeping accepted objective=" << accepted_global_objective
                     << " and best objective=" << best_sweep_global_objective << "." << endl;
                break;
            }
        }

        sweep_history.close();
        cout << "  wrote " << sweep_history_csv_file << "\n";

        if (best_sweep_index >= 0) {
            if (best_sweep_index + 1 != (int)global_objective_history.size()) {
                cout << "Restoring best coupled sweep " << (best_sweep_index + 1)
                     << " by global objective=" << best_sweep_global_objective
                     << " instead of last completed sweep." << endl;
            } else {
                cout << "Last completed sweep is also best by global objective." << endl;
            }
            fit_results = best_sweep_results;
            fit_success = best_sweep_success;
            auto refresh_counts = refresh_cross_boundary_ext(fit_results, fit_success);
            cout << "Best-sweep restore barrier: refreshed "
                 << refresh_counts.first << " out-of-section photon assignments";
            if (refresh_counts.second > 0) cout << " (interpolation fallback=" << refresh_counts.second << ")";
            cout << "." << endl;
            fout.cd();
            TNamed("best_coupled_sweep_index", Form("%d", best_sweep_index + 1)).Write();
            TNamed("best_coupled_sweep_global_objective", Form("%.17g", best_sweep_global_objective)).Write();
        }

        update_final_fit_status();

        if (use_global_p_e_scale && Config::ENABLE_FINAL_GLOBAL_PE_SCALE_STAGE4) {
            cout << "\n==== FINAL GLOBAL P_E_SCALE: refinement with final section smearing ====" << endl;
            double best_chi2_global = 0.0;
            bool pe_ok = optimize_global_p_e_scale_for_fit(fit_results, fit_success,
                                                            global_p_e_scale,
                                                            "Final global p_e_scale",
                                                            &best_chi2_global);
            if (pe_ok) {
                for (int is = 0; is < nsec; ++is) {
                    if (fit_success[is]) fit_results[is].p_e_scale = global_p_e_scale;
                }
                cout << "Final global p_e_scale result: global p_e_scale=" << global_p_e_scale
                    << " (chi2_mmiss=" << best_chi2_global << ")" << endl;
            }
        } else if (use_global_p_e_scale) {
            cout << "\n==== FINAL GLOBAL P_E_SCALE: skipped by configuration (ENABLE_FINAL_GLOBAL_PE_SCALE_STAGE4=false) ====" << endl;
        }

        // Recompute final chi2 with finalized coefficients and finalized p_e_scale
        // so reported quality matches the exact model used for final output plots.
        cout << "\n==== Recomputing final section chi2 with finalized model state ====\n";
        vector<double> chi2_before_final_recompute(nsec, 1e300);
        vector<double> chi2_ndf_before_final_recompute(nsec, 1e300);
        for (int is = 0; is < nsec; ++is) {
            if (!fit_success[is]) continue;
            chi2_before_final_recompute[is] = fit_results[is].chi2;
            chi2_ndf_before_final_recompute[is] = fit_results[is].chi2_per_ndf();
        }
        auto selected_objective_from_breakdown = [&](const ObjectiveBreakdown &bd) {
            if (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_MPI0_ONLY) return bd.mpi0.chi2;
            if (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_MMISS_ONLY) return bd.mmiss.chi2;
            if (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_MPGG2_ONLY) return bd.mpgg2.chi2;
            return bd.total(Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);
        };
        ofstream objective_breakdown(objective_breakdown_csv_file.c_str());
        objective_breakdown << "scope,ix,iy,section,fit_status,selected_objective,"
                            << "chi2_mpi0,chi2_mmiss,chi2_mpgg2,"
                            << "informative_bins_mpi0,informative_bins_mmiss,informative_bins_mpgg2,"
                            << "empty_sim_data_bins_mpi0,empty_sim_data_bins_mmiss,empty_sim_data_bins_mpgg2,"
                            << "data_integral_mpi0,data_integral_mmiss,data_integral_mpgg2,"
                            << "sim_integral_mpi0,sim_integral_mmiss,sim_integral_mpgg2,"
                            << "nsmear\n";
        for (int is = 0; is < nsec; ++is) {
            if (!fit_success[is]) continue;
            TRandom3 rng_final_chi2(8000 + is);
            ObjectiveBreakdown bd = eval_objective_breakdown(fit_results[is].mu_a,
                                                fit_results[is].mu,
                                                fit_results[is].mu_c,
                                                fit_results[is].sigma,
                                                fit_results[is].sigma_pos,
                                                fit_results[is].p_e_scale,
                                                sim_events_per_section[is],
                                                *hdata_sec[is], *hdata_mmiss_sec[is], *hdata_mpgg2_sec[is],
                                                rng_final_chi2, Nsmear,
                                                Config::RESOLUTION_A_DEFAULT,
                                                Config::RESOLUTION_B_DEFAULT,
                                                Config::RESOLUTION_C_DEFAULT);
            fit_results[is].chi2 = selected_objective_from_breakdown(bd);

            bool good_fit = (fit_results[is].chi2_per_ndf() <= Config::MAX_CHI2_PER_NDF);
            if (Config::SKIP_BAD_FITS && !good_fit) {
                fit_status[is] = "poor_fit_skipped";
                fit_success[is] = false;
            } else {
                fit_status[is] = good_fit ? "fit_ok" : "poor_fit";
            }
            objective_breakdown << "section,"
                                << sections[is].ix << "," << sections[is].iy << ","
                                << sections[is].name() << ","
                                << fit_status[is] << ","
                                << fit_results[is].chi2 << ","
                                << bd.mpi0.chi2 << "," << bd.mmiss.chi2 << "," << bd.mpgg2.chi2 << ","
                                << bd.mpi0.informative_bins << "," << bd.mmiss.informative_bins << "," << bd.mpgg2.informative_bins << ","
                                << bd.mpi0.empty_sim_data_positive_bins << "," << bd.mmiss.empty_sim_data_positive_bins << "," << bd.mpgg2.empty_sim_data_positive_bins << ","
                                << bd.mpi0.data_integral << "," << bd.mmiss.data_integral << "," << bd.mpgg2.data_integral << ","
                                << bd.mpi0.sim_integral << "," << bd.mmiss.sim_integral << "," << bd.mpgg2.sim_integral << ","
                                << Nsmear << "\n";
        }
        {
            ObjectiveBreakdown global_bd;
            int global_interp = 0;
            double global_total = evaluate_global_objective(fit_results, fit_success,
                                                            Nsmear,
                                                            &global_bd,
                                                            &global_interp);
            double global_selected = selected_objective_from_breakdown(global_bd);
            objective_breakdown << "global,-1,-1,all_sections,final,"
                                << global_selected << ","
                                << global_bd.mpi0.chi2 << "," << global_bd.mmiss.chi2 << "," << global_bd.mpgg2.chi2 << ","
                                << global_bd.mpi0.informative_bins << "," << global_bd.mmiss.informative_bins << "," << global_bd.mpgg2.informative_bins << ","
                                << global_bd.mpi0.empty_sim_data_positive_bins << "," << global_bd.mmiss.empty_sim_data_positive_bins << "," << global_bd.mpgg2.empty_sim_data_positive_bins << ","
                                << global_bd.mpi0.data_integral << "," << global_bd.mmiss.data_integral << "," << global_bd.mpgg2.data_integral << ","
                                << global_bd.mpi0.sim_integral << "," << global_bd.mmiss.sim_integral << "," << global_bd.mpgg2.sim_integral << ","
                                << Nsmear << "\n";
            cout << "Final global objective=" << global_total
                 << " [Mgg=" << global_bd.mpi0.chi2
                 << ", Mmiss=" << global_bd.mmiss.chi2
                 << ", Mpgg2=" << global_bd.mpgg2.chi2
                 << ", interp=" << global_interp << "]\n";
        }
        objective_breakdown.close();
        cout << "  wrote " << objective_breakdown_csv_file << "\n";

        {
            ofstream closure_summary(closure_summary_csv_file.c_str());
            closure_summary << "ix,iy,section,fit_status,optimizer_mode,optimizer_chi2_subset,"
                            << "full_final_chi2,delta_full_minus_optimizer,optimizer_chi2_ndf_subset,"
                            << "full_final_chi2_ndf,delta_full_minus_optimizer_chi2_ndf,"
                            << "optimizer_nsmear,final_nsmear,n_sim_optimizer,n_sim_full\n";
            for (int is = 0; is < nsec; ++is) {
                const FitResult &fr = fit_results[is];
                const double chi2_before = chi2_before_final_recompute[is];
                const double chi2_ndf_before = chi2_ndf_before_final_recompute[is];
                const bool have_before = std::isfinite(chi2_before) && chi2_before < 1e299;
                const bool have_ndf_before = std::isfinite(chi2_ndf_before) && chi2_ndf_before < 1e299;
                closure_summary << sections[is].ix << "," << sections[is].iy << ","
                                << sections[is].name() << ","
                                << fit_status[is] << ","
                                << fr.optimizer_mode << ","
                                << chi2_before << ","
                                << fr.chi2 << ","
                                << (have_before ? fr.chi2 - chi2_before : 0.0) << ","
                                << chi2_ndf_before << ","
                                << fr.chi2_per_ndf() << ","
                                << (have_ndf_before ? fr.chi2_per_ndf() - chi2_ndf_before : 0.0)
                                << "," << optimizer_nsmear
                                << "," << Nsmear
                                << "," << sim_events_opt_per_section[is].size()
                                << "," << sim_events_per_section[is].size()
                                << "\n";
            }
            cout << "  wrote " << closure_summary_csv_file << " (optimizer subset vs full final chi2)\n";
        }

        cout << "\n==== Writing optimizer diagnostics ====\n";
        {
            ofstream opt_summary(optimizer_summary_csv_file.c_str());
            opt_summary << "ix,iy,section,fit_status,optimizer_mode,n_seed_trials,n_migrad_trials,"
                        << "hesse_ok,max_abs_corr,max_corr_pair,chi2,chi2_ndf,"
                        << "mu_a,mu_b,mu_c,sigma,sigma_pos,p_e_scale\n";
            if (Config::ENABLE_GLOBAL_PREFIT_SEED) {
                const FitResult &fr = global_prefit_result;
                opt_summary << "-1,-1,global_prefit," << global_prefit_status << ","
                            << fr.optimizer_mode << ","
                            << fr.n_seed_trials << "," << fr.n_migrad_trials << ","
                            << (fr.hesse_ok ? 1 : 0) << ","
                            << fr.max_abs_corr << "," << fr.max_corr_pair << ","
                            << fr.chi2 << "," << fr.chi2_per_ndf() << ","
                            << fr.mu_a << "," << fr.mu << "," << fr.mu_c << ","
                            << fr.sigma << "," << fr.sigma_pos << "," << fr.p_e_scale << "\n";
            }
            for (int is = 0; is < nsec; ++is) {
                const FitResult &fr = fit_results[is];
                opt_summary << sections[is].ix << "," << sections[is].iy << ","
                            << sections[is].name() << "," << fit_status[is] << ","
                            << fr.optimizer_mode << ","
                            << fr.n_seed_trials << "," << fr.n_migrad_trials << ","
                            << (fr.hesse_ok ? 1 : 0) << ","
                            << fr.max_abs_corr << "," << fr.max_corr_pair << ","
                            << fr.chi2 << "," << fr.chi2_per_ndf() << ","
                            << fr.mu_a << "," << fr.mu << "," << fr.mu_c << ","
                            << fr.sigma << "," << fr.sigma_pos << "," << fr.p_e_scale << "\n";
            }
            cout << "  wrote " << optimizer_summary_csv_file << "\n";
        }

        {
            ofstream opt_seeds(optimizer_seeds_csv_file.c_str());
            opt_seeds << "ix,iy,section,rank,seed_index,seed_chi2,migrad_chi2,minimized,hesse_ok,"
                      << "max_abs_corr,max_corr_pair,mu_a,mu_b,mu_c,sigma,sigma_pos,p_e_scale\n";
            if (Config::ENABLE_GLOBAL_PREFIT_SEED) {
                const FitResult &fr = global_prefit_result;
                for (const auto &sd : fr.seed_diagnostics) {
                    opt_seeds << "-1,-1,global_prefit,"
                              << sd.rank << "," << sd.seed_index << ","
                              << sd.seed_chi2 << "," << sd.migrad_chi2 << ","
                              << (sd.minimized ? 1 : 0) << "," << (sd.hesse_ok ? 1 : 0) << ","
                              << sd.max_abs_corr << "," << sd.max_corr_pair << ","
                              << sd.mu_a << "," << sd.mu << "," << sd.mu_c << ","
                              << sd.sigma << "," << sd.sigma_pos << "," << sd.p_e_scale << "\n";
                }
            }
            for (int is = 0; is < nsec; ++is) {
                const FitResult &fr = fit_results[is];
                for (const auto &sd : fr.seed_diagnostics) {
                    opt_seeds << sections[is].ix << "," << sections[is].iy << ","
                              << sections[is].name() << ","
                              << sd.rank << "," << sd.seed_index << ","
                              << sd.seed_chi2 << "," << sd.migrad_chi2 << ","
                              << (sd.minimized ? 1 : 0) << "," << (sd.hesse_ok ? 1 : 0) << ","
                              << sd.max_abs_corr << "," << sd.max_corr_pair << ","
                              << sd.mu_a << "," << sd.mu << "," << sd.mu_c << ","
                              << sd.sigma << "," << sd.sigma_pos << "," << sd.p_e_scale << "\n";
                }
            }
            cout << "  wrote " << optimizer_seeds_csv_file << "\n";
        }

        {
            ofstream opt_profile(optimizer_profile_csv_file.c_str());
            opt_profile << "ix,iy,section,parameter,fixed_value,chi2,minimized\n";
            if (Config::ENABLE_GLOBAL_PREFIT_SEED) {
                const FitResult &fr = global_prefit_result;
                for (const auto &pd : fr.profile_diagnostics) {
                    opt_profile << "-1,-1,global_prefit,"
                                << pd.parameter << ","
                                << pd.fixed_value << "," << pd.chi2 << ","
                                << (pd.minimized ? 1 : 0) << "\n";
                }
            }
            for (int is = 0; is < nsec; ++is) {
                const FitResult &fr = fit_results[is];
                for (const auto &pd : fr.profile_diagnostics) {
                    opt_profile << sections[is].ix << "," << sections[is].iy << ","
                                << sections[is].name() << ","
                                << pd.parameter << ","
                                << pd.fixed_value << "," << pd.chi2 << ","
                                << (pd.minimized ? 1 : 0) << "\n";
                }
            }
            cout << "  wrote " << optimizer_profile_csv_file << "\n";
        }

        // Regenerate visualization scan data from FINAL model state only
        // (latest coupled sweep + final p_e_scale/status state).
        cout << "\n==== Regenerating visualization scans from final sweep values ====\n";
        #pragma omp parallel if(Config::ENABLE_PARALLEL_SECTION_FITS)
        {
            int thread_id_scan = omp_get_thread_num();
            TRandom3 rng_final_scan(8500 + thread_id_scan * 104729 + time(NULL));

            #pragma omp for schedule(dynamic)
            for (int is = 0; is < nsec; ++is) {
                if (!fit_success[is]) continue;
                generate_visualization_scan_data(sim_events_opt_per_section[is],
                                                *hdata_sec[is], *hdata_mmiss_sec[is], *hdata_mpgg2_sec[is],
                                                rng_final_scan, optimizer_nsmear,
                                                fit_results[is]);
            }
        }

        // =========================================================================
        // === Final per-section outputs with full-stat model state
        // =========================================================================
        cout << "\n==== Writing final per-section histograms and PDF pages ====" << endl;
        ofstream csv(csv_file.c_str());
        csv << "ix,iy,xlo,xhi,ylo,yhi,n_data,n_sim,best_mu_a,best_mu_b,best_mu_c,best_sigma,best_sigma_pos_cm,best_p_e_scale,res_A,res_B,res_C,best_chi2,ndf,chi2_ndf,fit_status\n";
        for (int is = 0; is < nsec; ++is) {
            size_t ndata_sec = (size_t)data_selected_count_per_section[is];
            size_t nsim_sec = (size_t)sim_selected_count_per_section[is];

            if (ndata_sec < Config::MIN_EVENTS_PER_SECTION || nsim_sec < Config::MIN_EVENTS_PER_SECTION || !fit_success[is]) {
                TDirectory *dir = fout.mkdir(sections[is].name().c_str());
                dir->cd();
                hdata_sec[is]->Write(hdata_sec[is]->GetName());
                hdata_mmiss_sec[is]->Write(hdata_mmiss_sec[is]->GetName());
                hdata_mpgg2_sec[is]->Write(hdata_mpgg2_sec[is]->GetName());
                fout.cd();
                if (fit_status[is] == "not_fitted") {
                    fit_status[is] = (ndata_sec < Config::MIN_EVENTS_PER_SECTION || nsim_sec < Config::MIN_EVENTS_PER_SECTION)
                                    ? "low_stats" : "fit_failed";
                }
                csv << sections[is].ix << "," << sections[is].iy << "," << sections[is].xlo << "," << sections[is].xhi << ","
                    << sections[is].ylo << "," << sections[is].yhi << "," << ndata_sec << "," << nsim_sec << ",,,,,,,,,,,,,"
                    << fit_status[is] << "\n";
                continue;
            }

            const FitResult &fres = fit_results[is];
            TRandom3 rng_out(9000 + is);

            TH1D hunsmeared((string("h_unsmeared_") + sections[is].name()).c_str(), "h_unsmeared", Config::MGGAMMA_NBINS, Config::MGGAMMA_MIN, Config::MGGAMMA_MAX);
            TH1D hfinal((string("h_smeared_best_") + sections[is].name()).c_str(), "h_smeared_best", Config::MGGAMMA_NBINS, Config::MGGAMMA_MIN, Config::MGGAMMA_MAX);
            TH1D hfinal_mmiss((string("h_smeared_best_mmiss_") + sections[is].name()).c_str(), "h_smeared_best_mmiss", Config::MMISS_NBINS, Config::MMISS_MIN, Config::MMISS_MAX);
            TH1D hfinal_mpgg2((string("h_smeared_best_mpgg2_") + sections[is].name()).c_str(), "h_smeared_best_mpgg2", Config::MPGG2_NBINS, Config::MPGG2_MIN, Config::MPGG2_MAX);
            hunsmeared.Sumw2(); hfinal.Sumw2(); hfinal_mmiss.Sumw2(); hfinal_mpgg2.Sumw2();

            for (const auto &ev : sim_events_per_section[is]) {
                if (ev.e1 <= 0 || ev.e2 <= 0) continue;
                PhotonMomentum p1u = computePhotonMomentum(ev.e1, ev.x1, ev.y1, nps::kDefaultZ_NPS_cm);
                PhotonMomentum p2u = computePhotonMomentum(ev.e2, ev.x2, ev.y2, nps::kDefaultZ_NPS_cm);
                double mass_u = computeInvariantMass(ev.e1, p1u.px, p1u.py, p1u.pz,
                                                    ev.e2, p2u.px, p2u.py, p2u.pz);
                hunsmeared.Fill(mass_u, ev.weight);
            }

            fillSmearedHistogramsAtParams(sim_events_per_section[is],
                                        hfinal, hfinal_mmiss, hfinal_mpgg2,
                                        fres.mu_a, fres.mu, fres.mu_c,
                                        fres.sigma, fres.sigma_pos,
                                        fres.p_e_scale,
                                        rng_out, Nsmear,
                                        fres.res_A, fres.res_B, fres.res_C);

            const double z_nps_diag = nps::kDefaultZ_NPS_cm;
            string sec_tag = sections[is].name();
            TH1D hunsmeared_mmiss(("h_unsmeared_mmiss_" + sec_tag).c_str(), "h_unsmeared_mmiss",
                                Config::MMISS_NBINS, Config::MMISS_MIN, Config::MMISS_MAX);
            TH1D hunsmeared_mpgg2(("h_unsmeared_mpgg2_" + sec_tag).c_str(), "h_unsmeared_mpgg2",
                                Config::MPGG2_NBINS, Config::MPGG2_MIN, Config::MPGG2_MAX);
            hunsmeared_mmiss.Sumw2();
            hunsmeared_mpgg2.Sumw2();
            for (const auto &ev : sim_events_per_section[is]) {
                if (ev.e1 <= 0 || ev.e2 <= 0) continue;
                double mmiss_u = nps::missing_mass_proton_pi0(
                    Config::BEAM_ENERGY, ev.Ee, ev.px_e, ev.py_e, ev.pz_e,
                    ev.e1, ev.e2, ev.x1, ev.y1, ev.x2, ev.y2);
                hunsmeared_mmiss.Fill(mmiss_u, ev.weight);
                double E1u = Config::W_MPGG2_ENERGY * ev.e1;
                double E2u = Config::W_MPGG2_ENERGY * ev.e2;
                PhotonMomentum p1u = computePhotonMomentum(E1u, ev.x1, ev.y1, z_nps_diag);
                PhotonMomentum p2u = computePhotonMomentum(E2u, ev.x2, ev.y2, z_nps_diag);
                hunsmeared_mpgg2.Fill(computeTargetPlusDiphotonMass2(E1u, p1u, E2u, p2u), ev.weight);
            }
            auto scaleUnsmearedAndSmeared = [](TH1D &hunsm, TH1D &hsm, const TH1D &hdata_ref) {
                const double integral_unsm = hunsm.Integral();
                const double integral_data = hdata_ref.Integral();
                if (integral_unsm <= 0.0 || integral_data <= 0.0) return;
                const double scale = integral_data / integral_unsm;
                hunsm.Scale(scale);
                hsm.Scale(scale);
            };
            scaleUnsmearedAndSmeared(hunsmeared, hfinal, *hdata_sec[is]);
            scaleUnsmearedAndSmeared(hunsmeared_mmiss, hfinal_mmiss, *hdata_mmiss_sec[is]);
            scaleUnsmearedAndSmeared(hunsmeared_mpgg2, hfinal_mpgg2, *hdata_mpgg2_sec[is]);

            size_t ndata_sec_diag = data_mass_per_section[is].size();
            int nsel_sec_diag = data_selected_count_per_section[is];
            int nsim_sec_diag = sim_selected_count_per_section[is];
            plot_section_diagnostics(sections[is], fres, c_chi2,
                                     pdf_file, page_count,
                                     diagnostic_canvas_dir,
                                     *hdata_sec[is], hfinal, hunsmeared,
                                     *hdata_mmiss_sec[is], hfinal_mmiss, hunsmeared_mmiss,
                                     *hdata_mpgg2_sec[is], hfinal_mpgg2, hunsmeared_mpgg2,
                                     (int)ndata_sec_diag, nsel_sec_diag, nsim_sec_diag,
                                     Nsmear, global_p_e_scale, use_global_p_e_scale);

            TDirectory *dir = fout.mkdir(sections[is].name().c_str());
            dir->cd();
            hdata_sec[is]->Write(hdata_sec[is]->GetName());
            hdata_mmiss_sec[is]->Write(hdata_mmiss_sec[is]->GetName());
            hdata_mpgg2_sec[is]->Write(hdata_mpgg2_sec[is]->GetName());
            hunsmeared.Write(hunsmeared.GetName());
            hfinal.Write(hfinal.GetName());
            hfinal_mmiss.Write(hfinal_mmiss.GetName());
            hfinal_mpgg2.Write(hfinal_mpgg2.GetName());
            hunsmeared_mmiss.Write(hunsmeared_mmiss.GetName());
            hunsmeared_mpgg2.Write(hunsmeared_mpgg2.GetName());

            TPaveText *pt = new TPaveText(0.1, 0.7, 0.9, 0.9, "brNDC");
            if (Config::ENABLE_ENERGY_DEPENDENT_MU) {
                pt->AddText(Form("mu_eff(E) = a + b*E + c*ln(E):  a=%.5f  b=%.5f  c=%.5f",
                                fres.mu_a, fres.mu, fres.mu_c));
            } else {
                pt->AddText(Form("Best mu (b) = %.5f  (a=0, c=0 forced)", fres.mu));
            }
            pt->AddText(Form("Best sigma = %.5f", fres.sigma));
            pt->AddText(Form("Best sigma_pos = %.5f cm", fres.sigma_pos));
            if (Config::ENABLE_ELECTRON_MOMENTUM_SCALING && Config::fit_objective_uses_mmiss()) {
                pt->AddText(Form("Best p_e_scale = %.6f", fres.p_e_scale));
            }
            pt->AddText(Form("Chi2 = %.2f, ndf = %d, chi2/ndf = %.2f", fres.chi2, fres.ndf, fres.chi2_per_ndf()));
            pt->Write("fit_params");
            delete pt;

            fout.cd();
            csv << sections[is].ix << "," << sections[is].iy << ","
                << sections[is].xlo << "," << sections[is].xhi << "," << sections[is].ylo << "," << sections[is].yhi << ","
                << ndata_sec << "," << nsim_sec << ","
                << fres.mu_a << "," << fres.mu << "," << fres.mu_c << ","
                << fres.sigma << "," << fres.sigma_pos << "," << fres.p_e_scale << ","
                << fres.res_A << "," << fres.res_B << "," << fres.res_C << ","
                << fres.chi2 << "," << fres.ndf << "," << fres.chi2_per_ndf() << "," << fit_status[is] << "\n";
        }
        
        // Restore ROOT's histogram directory status
        TH1::AddDirectory(original_adddir_status);
        
        // ========================================================================
        // FIT QUALITY SUMMARY
        // ========================================================================
        cout << "\n==== Fit Quality Summary ====" << endl;
        int n_fitted = 0;
        int n_good_fits = 0;
        int n_poor_fits = 0;
        double sum_chi2_ndf = 0.0;
        double max_chi2_ndf = 0.0;
        double min_chi2_ndf = 1e10;
        
        for (int is = 0; is < nsec; ++is) {
            if (fit_success[is]) {
                n_fitted++;
                double chi2_ndf = fit_results[is].chi2_per_ndf();
                sum_chi2_ndf += chi2_ndf;
                max_chi2_ndf = max(max_chi2_ndf, chi2_ndf);
                min_chi2_ndf = min(min_chi2_ndf, chi2_ndf);
                
                if (chi2_ndf <= Config::MAX_CHI2_PER_NDF) {
                    n_good_fits++;
                } else {
                    n_poor_fits++;
                    cout << "  Section " << sections[is].name() << ": chi2/ndf = " << chi2_ndf 
                        << " [EXCEEDS THRESHOLD " << Config::MAX_CHI2_PER_NDF << "]" << endl;
                }
            }
        }
        
        cout << "Total sections: " << nsec << endl;
        cout << "Fitted sections: " << n_fitted << " (" << (100.0 * n_fitted / nsec) << "%)" << endl;
        cout << "Good fits (chi2/ndf <= " << Config::MAX_CHI2_PER_NDF << "): " << n_good_fits 
            << " (" << (n_fitted > 0 ? 100.0 * n_good_fits / n_fitted : 0.0) << "%)" << endl;
        cout << "Poor fits (chi2/ndf > " << Config::MAX_CHI2_PER_NDF << "): " << n_poor_fits 
            << " (" << (n_fitted > 0 ? 100.0 * n_poor_fits / n_fitted : 0.0) << "%)" << endl;
        
        if (n_fitted > 0) {
            cout << "Chi2/ndf range: [" << min_chi2_ndf << ", " << max_chi2_ndf << "]" << endl;
            cout << "Average chi2/ndf: " << (sum_chi2_ndf / n_fitted) << endl;
        }
        
        if (n_poor_fits > 0 && !Config::SKIP_BAD_FITS) {
            cout << "\n*** WARNING: " << n_poor_fits << " section(s) have poor fit quality ***" << endl;
            cout << "*** Consider: ***" << endl;
            cout << "***   1. Increasing statistics (check MIN_EVENTS_PER_SECTION)" << endl;
            cout << "***   2. Adjusting fit ranges (MU_MIN/MAX, SIGMA_MIN/MAX)" << endl;
            cout << "***   3. Checking data quality in these sections" << endl;
            cout << "***   4. Setting SKIP_BAD_FITS=true to exclude poor fits" << endl;
        }
        cout << "========================================" << endl;

        // ========================================================================
        // SUMMARY PAGES — Combined comparisons (all sections merged)
        // ========================================================================
        cout << "\n==== Creating combined summary pages ====" << endl;

        // Allocate combined histograms for all three observables
        TH1D *hd_c_mgg  = new TH1D("hd_c_mgg",  ";M_{#gamma#gamma} [GeV/c^{2}];Counts",
                                    Config::MGGAMMA_NBINS, Config::MGGAMMA_MIN, Config::MGGAMMA_MAX);
        TH1D *hd_c_mmiss = new TH1D("hd_c_mmiss", ";M_{miss} [GeV/c^{2}];Counts",
                                    Config::MMISS_NBINS, Config::MMISS_MIN, Config::MMISS_MAX);
        TH1D *hd_c_mpgg2 = new TH1D("hd_c_mpgg2", ";(p_{target}+#gamma#gamma)^{2} [GeV^{2}];Counts",
                                    Config::MPGG2_NBINS, Config::MPGG2_MIN, Config::MPGG2_MAX);

        TH1D *hu_c_mgg  = new TH1D("hu_c_mgg",  ";M_{#gamma#gamma} [GeV/c^{2}];Counts",
                                    Config::MGGAMMA_NBINS, Config::MGGAMMA_MIN, Config::MGGAMMA_MAX);
        TH1D *hu_c_mmiss = new TH1D("hu_c_mmiss", ";M_{miss} [GeV/c^{2}];Counts",
                                    Config::MMISS_NBINS, Config::MMISS_MIN, Config::MMISS_MAX);
        TH1D *hu_c_mpgg2 = new TH1D("hu_c_mpgg2", ";(p_{target}+#gamma#gamma)^{2} [GeV^{2}];Counts",
                                    Config::MPGG2_NBINS, Config::MPGG2_MIN, Config::MPGG2_MAX);
        hu_c_mgg->Sumw2(); hu_c_mmiss->Sumw2(); hu_c_mpgg2->Sumw2();

        TH1D *hs_c_mgg  = new TH1D("hs_c_mgg",  ";M_{#gamma#gamma} [GeV/c^{2}];Counts",
                                    Config::MGGAMMA_NBINS, Config::MGGAMMA_MIN, Config::MGGAMMA_MAX);
        TH1D *hs_c_mmiss = new TH1D("hs_c_mmiss", ";M_{miss} [GeV/c^{2}];Counts",
                                    Config::MMISS_NBINS, Config::MMISS_MIN, Config::MMISS_MAX);
        TH1D *hs_c_mpgg2 = new TH1D("hs_c_mpgg2", ";(p_{target}+#gamma#gamma)^{2} [GeV^{2}];Counts",
                                    Config::MPGG2_NBINS, Config::MPGG2_MIN, Config::MPGG2_MAX);
        hs_c_mgg->Sumw2(); hs_c_mmiss->Sumw2(); hs_c_mpgg2->Sumw2();

        // Fill combined data from unique summary histograms (strict BOTH photons in geometry)
        hd_c_mgg->Sumw2();
        hd_c_mmiss->Sumw2();
        hd_c_mpgg2->Sumw2();
        hd_c_mgg->Add(&hdata_summary_mpi0);
        hd_c_mmiss->Add(&hdata_summary_mmiss);
        hd_c_mpgg2->Add(&hdata_summary_mpgg2);

        // Build globally unique sim buffer (strict BOTH photons in geometry) with
        // per-photon section coefficients attached deterministically from section grid.
        auto section_center_x = [&](int is) {
            return x_min + (sections[is].ix + 0.5) * base_wx;
        };
        auto section_center_y = [&](int is) {
            return y_min + (sections[is].iy + 0.5) * base_wy;
        };
        auto section_index_from_base_grid = [&](double x, double y) {
            if (!in_geometry(x, y)) return -1;
            int ix = (int)floor((x - x_min) / base_wx);
            int iy = (int)floor((y - y_min) / base_wy);
            ix = max(0, min(nx - 1, ix));
            iy = max(0, min(ny - 1, iy));
            return iy * nx + ix;
        };

        // Per-photon section ownership for summary smearing:
        //  1) prefer fitted sections that contain the hit in overlap geometry
        //     (nearest section center wins if multiple overlaps)
        //  2) fallback to fitted base-grid cell index
        //  3) otherwise use interpolation fallback
        auto section_index_for_photon_coeff = [&](double x, double y) {
            if (!in_geometry(x, y)) return -1;

            int best_is = -1;
            double best_d2 = std::numeric_limits<double>::max();
            for (int is = 0; is < nsec; ++is) {
                if (!fit_success[is]) continue;
                if (!inSection(sections[is], x, y)) continue;
                double dx = x - section_center_x(is);
                double dy = y - section_center_y(is);
                double d2 = dx * dx + dy * dy;
                if (d2 < best_d2) {
                    best_d2 = d2;
                    best_is = is;
                }
            }
            if (best_is >= 0) return best_is;

            int base_is = section_index_from_base_grid(x, y);
            if (base_is >= 0 && fit_success[base_is]) return base_is;
            return -1;
        };

        CalibrationMap all_sections_cal_map(nx, ny, x_min, x_max, y_min, y_max);
        for (int is = 0; is < nsec; ++is) {
            if (!fit_success[is]) continue;
            all_sections_cal_map.setParams(sections[is].ix, sections[is].iy,
                                        fit_results[is].mu_a,
                                        fit_results[is].mu,
                                        fit_results[is].mu_c,
                                        fit_results[is].sigma,
                                        fit_results[is].sigma_pos);
        }

        vector<ClusterPair> sim_events_global_with_fit;
        sim_events_global_with_fit.reserve(sim_events_summary.size());
        int n_summary_missing_section_coeff = 0;
        for (const auto &ev : sim_events_summary) {
            ClusterPair gev = ev;
            gev.photon1_in_section = false;
            gev.photon2_in_section = false;
            gev.mu_a1_ext = 0.0; gev.mu1_ext = 1.0; gev.mu_c1_ext = 0.0;
            gev.sigma1_ext = 0.0; gev.sigma_pos1_ext = 0.0;
            gev.mu_a2_ext = 0.0; gev.mu2_ext = 1.0; gev.mu_c2_ext = 0.0;
            gev.sigma2_ext = 0.0; gev.sigma_pos2_ext = 0.0;

            int is1 = section_index_for_photon_coeff(ev.x1, ev.y1);
            if (is1 >= 0) {
                gev.mu_a1_ext = fit_results[is1].mu_a;
                gev.mu1_ext   = fit_results[is1].mu;
                gev.mu_c1_ext = fit_results[is1].mu_c;
                gev.sigma1_ext = fit_results[is1].sigma;
                gev.sigma_pos1_ext = fit_results[is1].sigma_pos;
            } else if (in_geometry(ev.x1, ev.y1)) {
                all_sections_cal_map.getInterpolatedParams(ev.x1, ev.y1,
                                                        gev.mu_a1_ext, gev.mu1_ext, gev.mu_c1_ext,
                                                        gev.sigma1_ext, gev.sigma_pos1_ext);
                ++n_summary_missing_section_coeff;
            }
            int is2 = section_index_for_photon_coeff(ev.x2, ev.y2);
            if (is2 >= 0) {
                gev.mu_a2_ext = fit_results[is2].mu_a;
                gev.mu2_ext   = fit_results[is2].mu;
                gev.mu_c2_ext = fit_results[is2].mu_c;
                gev.sigma2_ext = fit_results[is2].sigma;
                gev.sigma_pos2_ext = fit_results[is2].sigma_pos;
            } else if (in_geometry(ev.x2, ev.y2)) {
                all_sections_cal_map.getInterpolatedParams(ev.x2, ev.y2,
                                                        gev.mu_a2_ext, gev.mu2_ext, gev.mu_c2_ext,
                                                        gev.sigma2_ext, gev.sigma_pos2_ext);
                ++n_summary_missing_section_coeff;
            }
            sim_events_global_with_fit.push_back(gev);
        }
        if (n_summary_missing_section_coeff > 0) {
            cout << "All-sections summary: " << n_summary_missing_section_coeff
                << " photon coefficient lookups used interpolated fallback due to unfitted sections." << endl;
        }

        // Fill unsmeared and smeared combined histograms in one pass (global unique events)
        TRandom3 rng_combined(42);
        const double z_nps_comb = nps::kDefaultZ_NPS_cm;
        for (const auto &ev : sim_events_summary) {
            if (ev.e1 <= 0 || ev.e2 <= 0) continue;

            PhotonMomentum p1u = computePhotonMomentum(ev.e1, ev.x1, ev.y1, z_nps_comb);
            PhotonMomentum p2u = computePhotonMomentum(ev.e2, ev.x2, ev.y2, z_nps_comb);
            hu_c_mgg->Fill(computeInvariantMass(ev.e1, p1u.px, p1u.py, p1u.pz,
                                            ev.e2, p2u.px, p2u.py, p2u.pz), ev.weight);
            double mmiss_u = nps::missing_mass_proton_pi0(
                Config::BEAM_ENERGY, ev.Ee, ev.px_e, ev.py_e, ev.pz_e,
                ev.e1, ev.e2, ev.x1, ev.y1, ev.x2, ev.y2);
            hu_c_mmiss->Fill(mmiss_u, ev.weight);
            double E1um = Config::W_MPGG2_ENERGY * ev.e1, E2um = Config::W_MPGG2_ENERGY * ev.e2;
            PhotonMomentum p1um = computePhotonMomentum(E1um, ev.x1, ev.y1, z_nps_comb);
            PhotonMomentum p2um = computePhotonMomentum(E2um, ev.x2, ev.y2, z_nps_comb);
            hu_c_mpgg2->Fill(computeTargetPlusDiphotonMass2(E1um, p1um, E2um, p2um), ev.weight);
        }

        fillSmearedHistogramsAtParams(sim_events_global_with_fit,
                                    *hs_c_mgg, *hs_c_mmiss, *hs_c_mpgg2,
                                    0.0, 1.0, 0.0,   // mu_a, mu_b, mu_c (identity; ext coeffs used for all photons)
                                    0.0, 0.0,          // sigma, sigma_pos
                                    global_p_e_scale,
                                    rng_combined,
                                    Nsmear,
                                    Config::RESOLUTION_A_DEFAULT,
                                    Config::RESOLUTION_B_DEFAULT,
                                    Config::RESOLUTION_C_DEFAULT);

        // Normalize the unsmeared simulation to data first, then apply the same
        // scale factor to the smeared histogram. This keeps smearing-induced
        // migration into/out of the plotted window visible in the diagnostics.
        auto scaleUnsmearedAndSmearedToData = [](TH1D* hu, TH1D* hs, const TH1D* hd) {
            if (hu->Integral() <= 0.0 || hd->Integral() <= 0.0) return;
            const double scale = hd->Integral() / hu->Integral();
            hu->Scale(scale);
            hs->Scale(scale);
        };
        scaleUnsmearedAndSmearedToData(hu_c_mgg,   hs_c_mgg,   hd_c_mgg);
        scaleUnsmearedAndSmearedToData(hu_c_mmiss, hs_c_mmiss, hd_c_mmiss);
        scaleUnsmearedAndSmearedToData(hu_c_mpgg2, hs_c_mpgg2, hd_c_mpgg2);

        // ---- Summary diagnostics: attribute all-sections mismatch by section ----
        vector<double> mismatch_score(nsec, 0.0);
        vector<double> sim_owner_weight(nsec, 0.0);
        vector<double> data_owner_weight(nsec, 0.0);

        auto owner_section_for_centroid = [&](double x, double y) {
            int is = section_index_for_photon_coeff(x, y);
            if (is < 0) return -1;
            return is;
        };

        // Ownership by diphoton centroid (for attribution only)
        vector<vector<ClusterPair>> sim_events_by_owner(nsec);
        for (const auto &ev : sim_events_global_with_fit) {
            double cx = 0.5 * (ev.x1 + ev.x2);
            double cy = 0.5 * (ev.y1 + ev.y2);
            int owner = owner_section_for_centroid(cx, cy);
            if (owner >= 0) {
                sim_events_by_owner[owner].push_back(ev);
                sim_owner_weight[owner] += ev.weight;
            }
        }
        for (const auto &ev : data_events_summary) {
            double cx = 0.5 * (ev.x1 + ev.x2);
            double cy = 0.5 * (ev.y1 + ev.y2);
            int owner = owner_section_for_centroid(cx, cy);
            if (owner >= 0) data_owner_weight[owner] += ev.weight;
        }

        vector<TH1D*> hs_owner_mgg(nsec, nullptr);
        vector<TH1D*> hs_owner_mmiss(nsec, nullptr);
        vector<TH1D*> hs_owner_mpgg2(nsec, nullptr);

        TRandom3 rng_owner(424242);
        for (int is = 0; is < nsec; ++is) {
            if (!fit_success[is]) continue;
            hs_owner_mgg[is] = new TH1D(Form("hs_owner_mgg_%d", is), "", Config::MGGAMMA_NBINS, Config::MGGAMMA_MIN, Config::MGGAMMA_MAX);
            hs_owner_mmiss[is] = new TH1D(Form("hs_owner_mmiss_%d", is), "", Config::MMISS_NBINS, Config::MMISS_MIN, Config::MMISS_MAX);
            hs_owner_mpgg2[is] = new TH1D(Form("hs_owner_mpgg2_%d", is), "", Config::MPGG2_NBINS, Config::MPGG2_MIN, Config::MPGG2_MAX);
            hs_owner_mgg[is]->Sumw2();
            hs_owner_mmiss[is]->Sumw2();
            hs_owner_mpgg2[is]->Sumw2();

            fillSmearedHistogramsAtParams(sim_events_by_owner[is],
                                        *hs_owner_mgg[is], *hs_owner_mmiss[is], *hs_owner_mpgg2[is],
                                        0.0, 1.0, 0.0,   // mu_a, mu_b, mu_c (identity; ext coeffs used)
                                        0.0, 0.0,          // sigma, sigma_pos
                                        global_p_e_scale,
                                        rng_owner,
                                        Nsmear,
                                        Config::RESOLUTION_A_DEFAULT,
                                        Config::RESOLUTION_B_DEFAULT,
                                        Config::RESOLUTION_C_DEFAULT);
        }

        auto accumulate_mismatch_scores = [&](const TH1D* hd_global,
                                            const TH1D* hs_global,
                                            const vector<TH1D*> &hs_owner) {
            int nb = hd_global->GetNbinsX();
            for (int ib = 1; ib <= nb; ++ib) {
                double residual = fabs(hs_global->GetBinContent(ib) - hd_global->GetBinContent(ib));
                if (residual <= 0.0) continue;
                double owner_sum = 0.0;
                for (int is = 0; is < nsec; ++is) {
                    if (!hs_owner[is]) continue;
                    owner_sum += hs_owner[is]->GetBinContent(ib);
                }
                if (owner_sum <= 0.0) continue;
                for (int is = 0; is < nsec; ++is) {
                    if (!hs_owner[is]) continue;
                    double frac = hs_owner[is]->GetBinContent(ib) / owner_sum;
                    mismatch_score[is] += frac * residual;
                }
            }
        };

        accumulate_mismatch_scores(hd_c_mgg, hs_c_mgg, hs_owner_mgg);
        accumulate_mismatch_scores(hd_c_mmiss, hs_c_mmiss, hs_owner_mmiss);
        accumulate_mismatch_scores(hd_c_mpgg2, hs_c_mpgg2, hs_owner_mpgg2);

        vector<int> score_order;
        for (int is = 0; is < nsec; ++is) {
            if (fit_success[is]) score_order.push_back(is);
        }
        sort(score_order.begin(), score_order.end(), [&](int a, int b) {
            return mismatch_score[a] > mismatch_score[b];
        });

        double total_data_owner_w = 0.0;
        double total_sim_owner_w = 0.0;
        for (int is = 0; is < nsec; ++is) {
            total_data_owner_w += data_owner_weight[is];
            total_sim_owner_w += sim_owner_weight[is];
        }

        cout << "\n==== All-sections mismatch attribution (top contributors) ====\n";
        int top_n_print = min(12, (int)score_order.size());
        for (int ir = 0; ir < top_n_print; ++ir) {
            int is = score_order[ir];
            double f_data = (total_data_owner_w > 0.0) ? data_owner_weight[is] / total_data_owner_w : 0.0;
            double f_sim  = (total_sim_owner_w > 0.0) ? sim_owner_weight[is] / total_sim_owner_w : 0.0;
            cout << "  #" << (ir + 1)
                << " " << sections[is].name()
                << " score=" << mismatch_score[is]
                << " data_frac=" << f_data
                << " sim_frac=" << f_sim
                << " delta(sim-data)=" << (f_sim - f_data)
                << " chi2/ndf=" << fit_results[is].chi2_per_ndf()
                << "\n";
        }

        // ---- Summary page: 3 observables side by side ----
        {
            TCanvas *c_obs = new TCanvas("c_obs_summary", "All Observables — All Sections", 1600, 600);
            c_obs->Divide(3, 1);
            auto drawObs = [&](int pad, TH1D* hd, TH1D* hu, TH1D* hs, const char* title) {
                c_obs->cd(pad);
                gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.14);
                drawComparisonPad(*hd, *hs, *hu, hd->GetXaxis()->GetTitle(), Form("cobs%d", pad));
                TLatex tx; tx.SetNDC(); tx.SetTextSize(0.050);
                tx.DrawLatex(0.15, 0.93, title);
            };
            drawObs(1, hd_c_mgg,   hu_c_mgg,   hs_c_mgg,   "M_{#gamma#gamma} — All Sections");
            drawObs(2, hd_c_mmiss, hu_c_mmiss, hs_c_mmiss, "M_{miss} — All Sections");
            drawObs(3, hd_c_mpgg2, hu_c_mpgg2, hs_c_mpgg2, "(p_{target}+#gamma#gamma)^{2} — All Sections");
            c_obs->Print(pdf_file.c_str());
            writeCanvasToDir(diagnostic_canvas_dir, c_obs, "c_obs_summary");
            delete c_obs;
        }

        // ---- Summary page: Section-level attribution of all-sections mismatch ----
        {
            TCanvas *c_attr = new TCanvas("c_mismatch_attr", "All-sections mismatch attribution", 1500, 900);
            c_attr->Divide(2, 1);

            c_attr->cd(1);
            gPad->SetRightMargin(0.15); gPad->SetLeftMargin(0.10); gPad->SetBottomMargin(0.12);
            TH2D *h_attr_map = new TH2D("h_attr_map_diag",
                                        "All-sections mismatch contribution score;x [cm];y [cm]",
                                        nx, x_min, x_max, ny, y_min, y_max);
            for (int is = 0; is < nsec; ++is) {
                if (!fit_success[is]) continue;
                double cx = x_min + (sections[is].ix + 0.5) * base_wx;
                double cy = y_min + (sections[is].iy + 0.5) * base_wy;
                int bx = h_attr_map->GetXaxis()->FindBin(cx);
                int by = h_attr_map->GetYaxis()->FindBin(cy);
                h_attr_map->SetBinContent(bx, by, mismatch_score[is]);
            }
            h_attr_map->SetStats(0);
            h_attr_map->Draw("COLZ TEXT");

            c_attr->cd(2);
            gPad->SetLeftMargin(0.03); gPad->SetRightMargin(0.03);
            TPaveText *pt_attr = new TPaveText(0.02, 0.02, 0.98, 0.98, "NDC");
            pt_attr->SetBorderSize(1);
            pt_attr->SetFillColor(0);
            pt_attr->SetTextAlign(12);
            pt_attr->SetTextSize(0.026);
            pt_attr->AddText("Top section contributors to all-sections mismatch");
            pt_attr->AddText("(score = share-weighted |Smeared_global - Data_global| across all 3 observables)");
            pt_attr->AddText(" ");
            int top_n = min(10, (int)score_order.size());
            for (int ir = 0; ir < top_n; ++ir) {
                int is = score_order[ir];
                double f_data = (total_data_owner_w > 0.0) ? data_owner_weight[is] / total_data_owner_w : 0.0;
                double f_sim  = (total_sim_owner_w > 0.0) ? sim_owner_weight[is] / total_sim_owner_w : 0.0;
                pt_attr->AddText(Form("#%d  %s   score=%.4g   data_frac=%.4f   sim_frac=%.4f   #Delta=%.4f   #chi^{2}/ndf=%.3f",
                                    ir + 1,
                                    sections[is].name().c_str(),
                                    mismatch_score[is],
                                    f_data,
                                    f_sim,
                                    f_sim - f_data,
                                    fit_results[is].chi2_per_ndf()));
            }
            pt_attr->Draw();

            c_attr->Print(pdf_file.c_str());
            writeCanvasToDir(diagnostic_canvas_dir, c_attr, "c_mismatch_attr");
            writeHistToDir(diagnostic_map_dir, h_attr_map);
            delete pt_attr;
            delete h_attr_map;
            delete c_attr;
        }

        // ---- Summary page: M_γγ pull distribution (all sections combined) ----
        {
            TCanvas *c_pull = new TCanvas("c_pull_summary", "M_gg Pull — All Sections", 1200, 550);
            c_pull->Divide(2, 1);
            TH1D *hpull_c = nullptr;
            TH1D *hpull_dist = nullptr;

            c_pull->cd(1);
            gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.14);
            {
                int nb = hd_c_mgg->GetNbinsX();
                hpull_c = new TH1D("hpull_comb", "(Data#minusSmeared)/#sigma  (All Sections);M_{#gamma#gamma} [GeV/c^{2}];Pull",
                                   nb, hd_c_mgg->GetXaxis()->GetXmin(), hd_c_mgg->GetXaxis()->GetXmax());
                for (int ib = 1; ib <= nb; ++ib) {
                    double d = hd_c_mgg->GetBinContent(ib), s = hs_c_mgg->GetBinContent(ib);
                    double ed = hd_c_mgg->GetBinError(ib), es = hs_c_mgg->GetBinError(ib);
                    double denom = sqrt(ed*ed + es*es);
                    if (denom > 0) hpull_c->SetBinContent(ib, (d - s) / denom);
                }
                hpull_c->SetFillColor(kCyan-9); hpull_c->SetLineColor(kBlack);
                double pm = max(fabs(hpull_c->GetMinimum()), fabs(hpull_c->GetMaximum()));
                hpull_c->SetMaximum( max(pm*1.2, 3.0)); hpull_c->SetMinimum(-max(pm*1.2, 3.0));
                hpull_c->Draw("HIST");
                TLine *lz = new TLine(hd_c_mgg->GetXaxis()->GetXmin(), 0, hd_c_mgg->GetXaxis()->GetXmax(), 0);
                lz->SetLineColor(kRed); lz->SetLineWidth(2); lz->Draw();
            }
            c_pull->cd(2);
            gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.14);
            {
                // Pull distribution histogram (Gaussian check)
                int nb = hd_c_mgg->GetNbinsX();
                hpull_dist = new TH1D("hpull_dist", "Pull distribution;Pull value;Bins",
                                      30, -5.0, 5.0);
                for (int ib = 1; ib <= nb; ++ib) {
                    double d = hd_c_mgg->GetBinContent(ib), s = hs_c_mgg->GetBinContent(ib);
                    double ed = hd_c_mgg->GetBinError(ib), es = hs_c_mgg->GetBinError(ib);
                    double denom = sqrt(ed*ed + es*es);
                    if (denom > 0) hpull_dist->Fill((d - s) / denom);
                }
                hpull_dist->SetFillColor(kOrange-9); hpull_dist->SetLineColor(kBlack);
                hpull_dist->Draw("HIST");
                TLatex tx; tx.SetNDC(); tx.SetTextSize(0.045);
                tx.DrawLatex(0.15, 0.87, "Pull value distribution");
                tx.DrawLatex(0.15, 0.81, "(expect: Gaussian, #mu=0, #sigma=1)");
            }
            c_pull->Print(pdf_file.c_str());
            writeCanvasToDir(diagnostic_canvas_dir, c_pull, "c_pull_summary");
            writeHistToDir(diagnostic_map_dir, hpull_c);
            writeHistToDir(diagnostic_map_dir, hpull_dist);
            delete c_pull;
            delete hpull_c;
            delete hpull_dist;
        }

        // ---- Summary page: 2D parameter maps (mu_a, mu_b, mu_c, sigma, sigma_pos, p_e_scale, chi2/ndf) ----
        {
            TCanvas *c_maps = new TCanvas("c_param_maps", "Fit Parameter Maps", 2100, 1000);
            c_maps->Divide(4, 2);

            double base_wx_map = (x_max - x_min) / nx;
            double base_wy_map = (y_max - y_min) / ny;

            TH2D *h_mu_a_map  = new TH2D("h_mu_a_map_diag",  "#mu_{a} (const offset) per section;x [cm];y [cm]",
                                        nx, x_min, x_max, ny, y_min, y_max);
            TH2D *h_mu_map    = new TH2D("h_mu_map_diag",    "#mu_{b} (linear coeff) per section;x [cm];y [cm]",
                                        nx, x_min, x_max, ny, y_min, y_max);
            TH2D *h_mu_c_map  = new TH2D("h_mu_c_map_diag",  "#mu_{c} (log coeff) per section;x [cm];y [cm]",
                                        nx, x_min, x_max, ny, y_min, y_max);
            TH2D *h_sig_map   = new TH2D("h_sig_map_diag",   "#sigma per section;x [cm];y [cm]",
                                        nx, x_min, x_max, ny, y_min, y_max);
            TH2D *h_sigpos_map= new TH2D("h_sigpos_map_diag","#sigma_{pos} [cm] per section;x [cm];y [cm]",
                                        nx, x_min, x_max, ny, y_min, y_max);
            TH2D *h_pe_map    = new TH2D("h_pe_map_diag",    "p_{e} scale per section;x [cm];y [cm]",
                                        nx, x_min, x_max, ny, y_min, y_max);
            TH2D *h_chi2_map  = new TH2D("h_chi2_map_diag",  "#chi^{2}/ndf per section;x [cm];y [cm]",
                                        nx, x_min, x_max, ny, y_min, y_max);

            for (int is = 0; is < nsec; ++is) {
                if (!fit_success[is]) continue;
                double cx = x_min + (sections[is].ix + 0.5) * base_wx_map;
                double cy = y_min + (sections[is].iy + 0.5) * base_wy_map;
                int bx = h_mu_map->GetXaxis()->FindBin(cx);
                int by = h_mu_map->GetYaxis()->FindBin(cy);
                h_mu_a_map->SetBinContent(bx, by, fit_results[is].mu_a);
                h_mu_map->SetBinContent(bx, by, fit_results[is].mu);
                h_mu_c_map->SetBinContent(bx, by, fit_results[is].mu_c);
                h_sig_map->SetBinContent(bx, by, fit_results[is].sigma);
                h_sigpos_map->SetBinContent(bx, by, fit_results[is].sigma_pos);
                h_pe_map->SetBinContent(bx, by, fit_results[is].p_e_scale);
                h_chi2_map->SetBinContent(bx, by, fit_results[is].chi2_per_ndf());
            }

            auto drawMap = [](TH2D* h, int pad, TCanvas* c) {
                c->cd(pad);
                gPad->SetRightMargin(0.15); gPad->SetLeftMargin(0.10); gPad->SetBottomMargin(0.12);
                h->SetStats(0);
                h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
                h->GetZaxis()->SetLabelSize(0.035);
                h->Draw("COLZ TEXT");
            };
            drawMap(h_mu_a_map,   1, c_maps);
            drawMap(h_mu_map,     2, c_maps);
            drawMap(h_mu_c_map,   3, c_maps);
            drawMap(h_sig_map,    4, c_maps);
            drawMap(h_sigpos_map, 5, c_maps);
            drawMap(h_pe_map,     6, c_maps);
            drawMap(h_chi2_map,   7, c_maps);

            c_maps->Print(pdf_file.c_str());
            writeCanvasToDir(diagnostic_canvas_dir, c_maps, "c_param_maps");
            writeHistToDir(diagnostic_map_dir, h_mu_a_map);
            writeHistToDir(diagnostic_map_dir, h_mu_map);
            writeHistToDir(diagnostic_map_dir, h_mu_c_map);
            writeHistToDir(diagnostic_map_dir, h_sig_map);
            writeHistToDir(diagnostic_map_dir, h_sigpos_map);
            writeHistToDir(diagnostic_map_dir, h_pe_map);
            writeHistToDir(diagnostic_map_dir, h_chi2_map);
            delete c_maps;
            delete h_mu_a_map; delete h_mu_map; delete h_mu_c_map;
            delete h_sig_map; delete h_sigpos_map; delete h_pe_map; delete h_chi2_map;
        }

        // ---- Visualization: energy response ratio curves for a selection of sections ----
        {
            // Select up to 8 fitted sections evenly spaced in section index
            vector<int> sel_sections;
            {
                int step = max(1, (int)score_order.size() > 0 ? nsec / min(8, nsec) : 1);
                for (int is = 0; is < nsec && (int)sel_sections.size() < 8; is += step) {
                    if (fit_success[is]) sel_sections.push_back(is);
                }
                if (sel_sections.empty()) {
                    for (int is = 0; is < nsec && (int)sel_sections.size() < 8; ++is)
                        if (fit_success[is]) sel_sections.push_back(is);
                }
            }
            if (!sel_sections.empty()) {
                TCanvas *c_resp = new TCanvas("c_response_ratio_curves", "Energy Response Ratio Curves per Section", 1400, 700);
                c_resp->Divide(1, 1);
                c_resp->cd(1);
                gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.14);
                gPad->SetGridx(); gPad->SetGridy();

                const int NE = 80;
                const double E_lo = Config::RESPONSE_CURVE_E_MIN_GEV;
                const double E_hi = Config::RESPONSE_CURVE_E_MAX_GEV;
                int colors[] = {kBlue, kRed, kGreen+2, kMagenta, kOrange+7, kCyan+2, kViolet+5, kGray+2};
                TMultiGraph *mg = new TMultiGraph("mg_response_ratio", "Energy response ratio R(E)=#mu_{eff}(E)/E per section;E [GeV];#mu_{eff}(E)/E");
                TLegend *lg_me = new TLegend(0.65, 0.15, 0.95, 0.55);
                lg_me->SetBorderSize(1); lg_me->SetFillColor(0); lg_me->SetTextSize(0.026);

                for (int ii = 0; ii < (int)sel_sections.size(); ++ii) {
                    int is = sel_sections[ii];
                    const FitResult &fr = fit_results[is];
                    TGraph *g = new TGraph(NE);
                    for (int ie = 0; ie < NE; ++ie) {
                        double E = E_lo + (E_hi - E_lo) * ie / (NE - 1);
                        double E_s = std::max(E, Config::MU_ENERGY_MIN_GEV);
                        double mu_eff = fr.mu_a + fr.mu * E_s + fr.mu_c * std::log(E_s);
                        g->SetPoint(ie, E, (E_s > 0.0) ? mu_eff / E_s : 0.0);
                    }
                    g->SetLineColor(colors[ii % 8]); g->SetLineWidth(2);
                    mg->Add(g, "L");
                    double E_ref = 2.0;
                    double E_ref_s = std::max(E_ref, Config::MU_ENERGY_MIN_GEV);
                    double r_ref = (fr.mu_a + fr.mu * E_ref_s + fr.mu_c * std::log(E_ref_s)) / E_ref_s;
                    lg_me->AddEntry(g, Form("Sec %s  R(2 GeV)=%.4f",
                                            sections[is].name().c_str(), r_ref), "l");
                }
                mg->Draw("A");
                mg->GetXaxis()->SetLimits(E_lo, E_hi);
                TLine *l1 = new TLine(E_lo, 1.0, E_hi, 1.0);
                l1->SetLineColor(kBlack); l1->SetLineStyle(2); l1->SetLineWidth(2); l1->Draw();
                lg_me->Draw();
                c_resp->Print(pdf_file.c_str());
                writeCanvasToDir(diagnostic_canvas_dir, c_resp, "c_response_ratio_curves");
                delete c_resp;
            }
        }

        // ---- Visualization: 2D energy response ratio maps at fixed energies ----
        if (Config::ENABLE_ENERGY_DEPENDENT_MU) {
            const double E_fixed[] = {1.0, 2.0, 3.0, 4.0, 5.0};
            const int N_E_fixed = 5;
            TCanvas *c_resp2d = new TCanvas("c_response_ratio_maps", "Energy response ratio maps at fixed energies", 2200, 800);
            c_resp2d->Divide(N_E_fixed, 1);
            vector<TH2D*> response_maps;

            double base_wx_me = (x_max - x_min) / nx;
            double base_wy_me = (y_max - y_min) / ny;

            for (int ie = 0; ie < N_E_fixed; ++ie) {
                double E = E_fixed[ie];
                double E_s = std::max(E, Config::MU_ENERGY_MIN_GEV);
                TH2D *h = new TH2D(Form("h_response_ratio2d_E%.0fGeV", E),
                                   Form("Energy response ratio #mu_{eff}/E at E=%.1f GeV;x [cm];y [cm];#mu_{eff}/E", E),
                                   nx, x_min, x_max, ny, y_min, y_max);
                for (int is = 0; is < nsec; ++is) {
                    if (!fit_success[is]) continue;
                    const FitResult &fr = fit_results[is];
                    double cx = x_min + (sections[is].ix + 0.5) * base_wx_me;
                    double cy = y_min + (sections[is].iy + 0.5) * base_wy_me;
                    int bx = h->GetXaxis()->FindBin(cx);
                    int by = h->GetYaxis()->FindBin(cy);
                    double mu_eff = fr.mu_a + fr.mu * E_s + fr.mu_c * std::log(E_s);
                    h->SetBinContent(bx, by, (E_s > 0.0) ? mu_eff / E_s : 0.0);
                }
                response_maps.push_back(h);
                c_resp2d->cd(ie + 1);
                gPad->SetRightMargin(0.15); gPad->SetLeftMargin(0.10); gPad->SetBottomMargin(0.14);
                h->SetStats(0);
                h->SetMarkerSize(1.0);
                h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
                h->GetZaxis()->SetTitle("#mu_{eff}/E");
                h->GetZaxis()->SetLabelSize(0.035);
                h->Draw("COLZ TEXT");
            }
            c_resp2d->Print(pdf_file.c_str());
            writeCanvasToDir(diagnostic_canvas_dir, c_resp2d, "c_response_ratio_maps");
            for (TH2D *h : response_maps) {
                writeHistToDir(diagnostic_map_dir, h);
                delete h;
            }
            delete c_resp2d;
        }

        // Totals shared by summary pages (stats map page + final text report)
        long long total_ndata = 0;
        long long total_nsel = 0;
        long long total_nsim = 0;
        long long total_events_all_sections = 0;

        // ---- Summary page: Section statistics (N_data all, N_data selected, N_sim selected) ----
        {
            TCanvas *c_stats = new TCanvas("c_stats_maps", "Section Statistics", 2100, 550);
            c_stats->Divide(3, 1);

            double base_wx_st = (x_max - x_min) / nx;
            double base_wy_st = (y_max - y_min) / ny;

            TH2D *h_ndata_map = new TH2D("h_ndata_map_diag",
                                        "N_{data} (all) per section;x [cm];y [cm]",
                                        nx, x_min, x_max, ny, y_min, y_max);
            TH2D *h_nsel_map = new TH2D("h_nsel_map_diag",
                                        "N_{data} (selected by weight mode) per section;x [cm];y [cm]",
                                        nx, x_min, x_max, ny, y_min, y_max);
            TH2D *h_nsim_map  = new TH2D("h_nsim_map_diag",  "N_{sim} (selected by weight mode) per section;x [cm];y [cm]",
                                        nx, x_min, x_max, ny, y_min, y_max);

            total_ndata = 0;
            total_nsel = 0;
            total_nsim = 0;
            for (int is = 0; is < nsec; ++is) {
                double cx = x_min + (sections[is].ix + 0.5) * base_wx_st;
                double cy = y_min + (sections[is].iy + 0.5) * base_wy_st;
                int bx = h_ndata_map->GetXaxis()->FindBin(cx);
                int by = h_ndata_map->GetYaxis()->FindBin(cy);
                long long nd  = (long long)data_mass_per_section[is].size();
                long long ne  = (long long)data_selected_count_per_section[is];
                long long ns  = (long long)sim_selected_count_per_section[is];
                h_ndata_map->SetBinContent(bx, by, (double)nd);
                h_nsel_map->SetBinContent(bx, by, (double)ne);
                h_nsim_map->SetBinContent(bx, by,  (double)ns);
                total_ndata += nd;
                total_nsel += ne;
                total_nsim  += ns;
            }
            total_events_all_sections = total_ndata + total_nsim;

            auto drawStatsPad = [](TH2D* h, const char* total_label) {
                gPad->SetRightMargin(0.15); gPad->SetLeftMargin(0.10); gPad->SetBottomMargin(0.12);
                h->SetStats(0);
                h->GetXaxis()->SetTitleSize(0.05); h->GetYaxis()->SetTitleSize(0.05);
                h->Draw("COLZ TEXT");
                TLatex *ltot = new TLatex();
                ltot->SetNDC(); ltot->SetTextSize(0.045); ltot->SetTextFont(62);
                ltot->DrawLatex(0.12, 0.94, total_label);
            };

            c_stats->cd(1);
            drawStatsPad(h_ndata_map, Form("Total N_{{data}} (all) = %lld", total_ndata));

            c_stats->cd(2);
            drawStatsPad(h_nsel_map, Form("Total N_{{data}} (selected) = %lld", total_nsel));

            c_stats->cd(3);
            drawStatsPad(h_nsim_map, Form("Total N_{{sim}} (selected) = %lld", total_nsim));

            c_stats->Print(pdf_file.c_str());
            writeCanvasToDir(diagnostic_canvas_dir, c_stats, "c_stats_maps");
            writeHistToDir(diagnostic_map_dir, h_ndata_map);
            writeHistToDir(diagnostic_map_dir, h_nsel_map);
            writeHistToDir(diagnostic_map_dir, h_nsim_map);
            delete c_stats;
            delete h_ndata_map; delete h_nsel_map; delete h_nsim_map;
        }

        // ---- Final summary page: text-only fit/smearing report (no overlap with plots) ----
        {
            TCanvas *c_final = new TCanvas("c_final_summary", "Fit and Smearing Summary", 1600, 1100);
            c_final->cd();
            gPad->SetLeftMargin(0.04);
            gPad->SetRightMargin(0.04);
            gPad->SetTopMargin(0.04);
            gPad->SetBottomMargin(0.04);

            TPaveText *txt = new TPaveText(0.04, 0.04, 0.96, 0.96, "NDC");
            txt->SetFillColor(0);
            txt->SetBorderSize(1);
            txt->SetTextAlign(12);
            txt->SetTextSize(0.024);

            txt->AddText("NPS #pi^{0} Smearing Fit Summary");
            txt->AddText(" ");
            txt->AddText("Optimization mode: Sobol multistart + staged sigma_{pos}");
            txt->AddText(Form("Observable mode: %s", Config::histogram_mode_label()));
            txt->AddText(Form("Weights: w_{M_{#gamma#gamma}}=%.3f, w_{M_{miss}}=%.3f, w_{(p+#gamma#gamma)^{2}}=%.3f",
                            Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2));
            txt->AddText(Form("SIMC de-modeling: sim weight = full_weight/%s",
                            Config::SIM_MODEL_XSEC_BRANCH));
            txt->AddText(Form("Invalid SIMC model-weight events skipped: %lld",
                            sim_skipped_invalid_model_weight));
            txt->AddText(Form("Energy model: %s", Config::USE_SIMPLE_STOCHASTIC_MODEL ? "#sigma_{E}=#sigma#sqrt{E}" : "3-term resolution model"));
            txt->AddText(Form("Position smearing: %s", Config::ENABLE_POSITION_SMEARING ? "enabled" : "disabled"));
            txt->AddText(Form("N_{smear} per event: %d", Nsmear));
            txt->AddText(Form("Calorimeter grid: nx=%d, ny=%d, x=[%.1f, %.1f] cm, y=[%.1f, %.1f] cm",
                            nx, ny, x_min, x_max, y_min, y_max));
            txt->AddText(" ");
            txt->AddText(Form("Electron momentum scaling mode: %s",
                            use_global_p_e_scale
                                ? (Config::ENABLE_FINAL_GLOBAL_PE_SCALE_STAGE4
                                        ? "per-section coupled sweeps + final global p_e_scale"
                                        : "per-section coupled sweeps only (final global p_e_scale disabled)")
                                : "disabled"));
            if (use_global_p_e_scale) {
                if (Config::ENABLE_FINAL_GLOBAL_PE_SCALE_STAGE4) {
                    txt->AddText(Form("Global fitted p_{e} scale (final): %.6f", global_p_e_scale));
                } else {
                    double pe_min = std::numeric_limits<double>::max();
                    double pe_max = -std::numeric_limits<double>::max();
                    double pe_sum = 0.0;
                    int pe_n = 0;
                    for (int is = 0; is < nsec; ++is) {
                        if (!fit_success[is]) continue;
                        pe_min = min(pe_min, fit_results[is].p_e_scale);
                        pe_max = max(pe_max, fit_results[is].p_e_scale);
                        pe_sum += fit_results[is].p_e_scale;
                        ++pe_n;
                    }
                    if (pe_n > 0) {
                        txt->AddText(Form("Per-section p_{e} scale range: [%.6f, %.6f], mean=%.6f",
                                        pe_min, pe_max, pe_sum / pe_n));
                    }
                }
            }
            txt->AddText(" ");
            txt->AddText(Form("Sections total=%d, fitted=%d, good=%d, poor=%d",
                    nsec, n_fitted, n_good_fits, n_poor_fits));
            txt->AddText(Form("Fit-quality threshold: #chi^{2}/ndf <= %.2f", Config::MAX_CHI2_PER_NDF));
            if (n_fitted > 0) {
                txt->AddText(Form("Fit-quality range: min #chi^{2}/ndf=%.4f, max #chi^{2}/ndf=%.4f, average=%.4f",
                                min_chi2_ndf, max_chi2_ndf, sum_chi2_ndf / n_fitted));
            } else {
                txt->AddText("Fit-quality range: no successful section fits");
            }
            txt->AddText(Form("Section-aggregated totals: N_{data}(all)=%lld, N_{data}(selected)=%lld, N_{sim}(selected)=%lld",
                    total_ndata, total_nsel, total_nsim));
            txt->AddText(Form("All-sections total events: N_{data}(all) + N_{sim}(selected) = %lld", total_events_all_sections));
            txt->AddText(Form("PDF content: %d section diagnostic pages + summary pages (this page is the final text report)",
                            page_count));
            txt->AddText(Form("Run tag: %s", run_tag.c_str()));
            txt->AddText(Form("Diagnostics PDF: %s", pdf_file.c_str()));
            txt->AddText(Form("Progress PDFs: %s", progress_pdf_dir.c_str()));
            txt->AddText(Form("Artifact manifest: %s", metadata_manifest_file.c_str()));
            txt->AddText(" ");
            txt->AddText("Note: all-sections histograms use unique global event buffers with BOTH photons in geometry.");
            txt->AddText("All-sections smearing assigns coefficients per photon from its own base-grid section.");

            txt->Draw();
            c_final->Print(pdf_file.c_str());
            writeCanvasToDir(diagnostic_canvas_dir, c_final, "c_final_summary");

            delete txt;
            delete c_final;
        }

        // Clean up combined histograms
        delete hd_c_mgg;  delete hd_c_mmiss;  delete hd_c_mpgg2;
        delete hu_c_mgg;  delete hu_c_mmiss;  delete hu_c_mpgg2;
        delete hs_c_mgg;  delete hs_c_mmiss;  delete hs_c_mpgg2;
        for (int is = 0; is < nsec; ++is) {
            delete hs_owner_mgg[is];
            delete hs_owner_mmiss[is];
            delete hs_owner_mpgg2[is];
        }

        // Close the PDF file
        string pdf_close = pdf_file + "]";
        c_chi2->Print(pdf_close.c_str());
        delete c_chi2;
        cout << "\nChi^2 scan plots saved to " << pdf_file << endl;
        copyFileIfDifferent(pdf_file, canonical_pdf_file, "latest chi2 PDF");

        {
            ofstream cache_out(cache_fingerprint_file.c_str(), std::ios::binary);
            cache_out << current_cache_fingerprint;
            cout << "Cache fingerprint saved to " << cache_fingerprint_file << "\n";
        }

        fout.Close();
        copyFileIfDifferent(out_file, timestamped_out_file, "fitter ROOT output");
        csv.close();
        cout << "All done. Results written to "<<out_file<<" and "<<csv_file<<endl;
        
        // ============================================================================
        // Create interpolated calibration map for smooth parameter variation
        // ============================================================================
        cout << "\n==== Building interpolated calibration map ====\n";
        CalibrationMap calMap(nx, ny, x_min, x_max, y_min, y_max);
        
        // Populate with fitted parameters
        for (int is=0; is<nsec; ++is) {
            if (fit_success[is]) {
                calMap.setParams(sections[is].ix, sections[is].iy,
                            fit_results[is].mu_a,
                            fit_results[is].mu,
                            fit_results[is].mu_c,
                            fit_results[is].sigma,
                            fit_results[is].sigma_pos);
            }
        }
        
        // Save 2D interpolated maps for visualization
        calMap.saveAsHistogram(interp_file, run_tag, created_at_local);
        copyFileIfDifferent(interp_file, timestamped_interp_file, "interpolated ROOT output");
        writeSmearingManifest(metadata_manifest_file,
                              run_tag,
                              created_at_local,
                              data_file,
                              sim_file,
                              out_file,
                              timestamped_out_file,
                              interp_file,
                              timestamped_interp_file,
                              pdf_file,
                              canonical_pdf_file,
                              progress_pdf_dir,
                              current_cache_fingerprint);
        
        // Demonstrate usage: print interpolated values at a few test points
        cout << "\n==== Example interpolated values ====\n";
        calMap.printParamsAt(0.0, 0.0);   // center
        calMap.printParamsAt(-15.0, -15.0); // lower-left quadrant
        calMap.printParamsAt(15.0, 15.0);   // upper-right quadrant
        calMap.printParamsAt(5.5, -8.3);    // arbitrary point
        
        cout << "\n==== Usage for event-by-event correction ====\n";
        cout << "// Example: Apply smearing to a photon at position (x, y):\n";
        cout << "double mu_a, mu_b, mu_c, sigma, sigma_pos;\n";
        cout << "calMap.getInterpolatedParams(cluster_x, cluster_y, mu_a, mu_b, mu_c, sigma, sigma_pos);\n";
        cout << "// Energy-scale model: mu_eff(E) = mu_a + mu_b*E + mu_c*ln(E)\n";
        cout << "double E_safe = std::max(original_energy, " << Config::MU_ENERGY_MIN_GEV << ");\n";
        if (Config::ENABLE_ENERGY_DEPENDENT_MU) {
            cout << "double E_sc = mu_a + mu_b * E_safe + mu_c * std::log(E_safe);\n";
        } else {
            cout << "double E_sc = mu_b * E_safe;  // disabled mode: a=0, c=0\n";
        }
        if (Config::USE_SIMPLE_STOCHASTIC_MODEL) {
            cout << "double smeared_energy = rng.Gaus(E_sc, sigma * sqrt(E_sc));\n";
        } else {
            cout << "double sigma_E = computeEnergyResolution(E_sc, sigma, "
                << "Config::RESOLUTION_A_DEFAULT, Config::RESOLUTION_B_DEFAULT, Config::RESOLUTION_C_DEFAULT);\n";
            cout << "double smeared_energy = rng.Gaus(E_sc, sigma_E);\n";
        }
        if (Config::ENABLE_ENERGY_DEPENDENT_SIGMA_POS) {
            cout << "double smeared_x = cluster_x + rng.Gaus(0.0, sigma_pos_eff);\n";
            cout << "double smeared_y = cluster_y + rng.Gaus(0.0, sigma_pos_eff);\n";
        } else {
            cout << "double smeared_x = cluster_x + rng.Gaus(0.0, sigma_pos);\n";
            cout << "double smeared_y = cluster_y + rng.Gaus(0.0, sigma_pos);\n";
        }
        
        return 0;
    }
