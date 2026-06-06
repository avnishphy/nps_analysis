    // ============================================================================
    // File: nps_sim_smearing_new.C
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
    //   - Optional global p_e_scale fit from M_miss-only Stage 4 (single value over calorimeter)
    //   - Bilinear interpolation for smooth parameter variation across calorimeter
    //   - Outputs: discrete section parameters (CSV) + interpolated 2D maps (ROOT)
    //
    // Compile:
    //   g++ nps_sim_smearing_new.C `root-config --cflags --libs` -O2 -std=c++17 -fopenmp -I../src -o nps_sim_smearing_new
    //
    // Usage example: 
    //   ./nps_sim_smearing_new /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/combined_branches_LH2.root physics /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/simc_pi0_analysis_output.root simulation out_smear.root 11 16 -24 28 -34 34 0.2 50
    //
    //   ./nps_sim_smearing_new /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/production_wfpi0/combined_branches_LH2_wfpi0.root physics /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/production_wfpi0/simc_pi0_analysis_output.root simulation out_smear.root 11 16 -24 28 -34 34 0.2 50
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
    #include <memory>
    #include <omp.h>

    // Bring standard library types into global namespace for legacy headers
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
    #include "TDirectory.h"
    #include "TFitResult.h"
    #include "TGraph.h"
    #include "TGraph2D.h"
    #include "TMarker.h"
    #include "Math/Factory.h"
    #include "Math/Functor.h"
    #include "Math/Minimizer.h"

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
    //   - USE_THREE_STAGE_OPTIMIZATION = true   (Stage 1: energy, Stage 2: position, Stage 3: joint)
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
            // Control whether to use the simulation full_weight branch for event weights
            const bool USE_SIM_FULL_WEIGHT = false; // Set false to use unweighted simulation
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
        // Choose between staged sequential or simultaneous optimization
        //
        // STAGED (TRUE) [RECOMMENDED]:
        //   Stage 1: Fit mu, sigma using selected observable (sigma_pos=0)
        //   Stage 2: Fit sigma_pos using M_γγ only (mu, sigma fixed from Stage 1)
        //   Stage 3: Joint refinement of all (mu, sigma, sigma_pos) simultaneously
        //   Benefits: Physically motivated ordering; energy calibrated first (largest
        //             effect), then angular resolution, then joint cleanup.
        //
        // SIMULTANEOUS (FALSE):
        //   Fit all parameters together using weighted combined chi2
        //   Benefits: Can find global minimum in parameter space
        //   Drawback: Couples angular and energy systematics, slower convergence
        //
        const bool USE_THREE_STAGE_OPTIMIZATION = true;  // RECOMMENDED: true for best physics separation
        
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
        const double W_MPI0 = 3.0;   // Weight for invariant mass chi2 (M_γγ)
        const double W_MMISS = 0.0;  // Weight for missing mass chi2 (M_miss)
        const double W_MPGG2 = 1.0;  // Weight for (p_target + γγ)^2 chi2
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
        // NOTE: In three-stage mode, this controls Stage 1 and Stage 3 objectives.
        //       Stage 2 always fits sigma_pos using M_γγ only.
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

        // Optional energy-dependent nonlinearity in photon energy scale:
        //   mu_eff(E) = mu0 * ( (E0 / (E + E_shift))^b + 1 )
        // With E0=2 and b=4, low-energy photons are scaled more strongly and
        // the correction falls off for high-energy photons.
        const bool ENABLE_ENERGY_DEPENDENT_MU = false;
        const double MU_ENERGY_E0_GEV = 1.2;
        const double MU_ENERGY_MIN_GEV = 0.2;
        const double MU_ENERGY_SHIFT_GEV = 3.0;
        const double MU_ENERGY_B_EXPONENT = 4.0;

        const double MU_ENERGY_A = 1.0;   // initial guess
        const double MU_ENERGY_B = 0.0;   // GeV^-1
        const double MU_ENERGY_C = 0.0;   // dimensionless

        
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
        const bool ENABLE_POSITION_SMEARING = false;    // RECOMMENDED: true
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
        // Set ENABLE_ELECTRON_MOMENTUM_SCALING = false to keep legacy behavior.
        const bool ENABLE_ELECTRON_MOMENTUM_SCALING = false;  // RECOMMENDED: false unless you have specific reasons to suspect electron momentum scale issues
        const bool ENABLE_FINAL_GLOBAL_PE_SCALE_STAGE4 = false;  // If false, keep per-section p_e_scale from Pass-2 and skip final global Stage-4 fit
        const double GLOBAL_PE_SCALE_DEFAULT = 0.995;              // initial seed for global p_e_scale fit

        // Initial seeds used before global two-step fit:
        //   Step 1: fit (mu, sigma) with selected combined objective
        //   Step 2: Stage 4 fit p_e_scale with M_miss only, keeping (mu, sigma, sigma_pos) fixed
        const double GLOBAL_PE_FIT_MU = 1.0;
        const double GLOBAL_PE_FIT_SIGMA = 0.02;
        const double GLOBAL_PE_FIT_SIGMA_POS = 0.0;

        const double PE_SCALE_MIN = 0.98;
        const double PE_SCALE_MAX = 1.02;
        const int PE_SCALE_NSTEPS = 100;

        
        // Coarse grid search density (speeds up initial parameter space exploration)
        // Higher divisor = coarser initial search = faster but might miss narrow minima
        // RECOMMENDED: 3-4 for good balance of speed vs accuracy
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
        const double MGGAMMA_MIN = 0.128;    // GeV/c²
        const double MGGAMMA_MAX = 0.142;    // GeV/c²  
        const int MGGAMMA_NBINS = 50;      // ~0.4 MeV/bin for fine shape comparison
        
        // Missing mass (M_miss) - focused on proton mass region (938 MeV/c²)
        const double MMISS_MIN = 0.8;       // GeV/c²
        const double MMISS_MAX = 1.0;       // GeV/c²
        const int MMISS_NBINS = 120;        // ~4 MeV/bin for detailed peak comparison

        // (p_target + γγ)^2 in GeV²
        const double MPGG2_MIN = 7.8;
        const double MPGG2_MAX = 12.2;
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

        // Optional exclusivity gating in event weights:
        //   false (default): use all events
        //   true: multiply weights by is_exclusive when branch is available
        const bool APPLY_IS_EXCLUSIVE_SELECTION = true;
        
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

        // Iterative coupled section fitting (Pass-2)
        // Refit sections in sweeps; after each sweep, refresh out-of-section photon
        // coefficients from latest section fits to propagate cross-boundary coupling.
        const bool ENABLE_ITERATIVE_COUPLED_PASS2 = true;
        const int COUPLED_PASS2_MAX_SWEEPS = 5;
        const double COUPLED_CONV_MU = 5e-4;
        const double COUPLED_CONV_SIGMA = 5e-4;
        const double COUPLED_CONV_SIGMA_POS = 5e-4;
        
        // ========================================================================
        // OUTPUT FILE SETTINGS
        // ========================================================================
        
        // Default output directory (prepended to relative paths)
        // const string OUTPUT_DIR = "../output/plots/x60_4b/";
        const string OUTPUT_DIR = "../output/plots/x60_4b/production_wfpi0/";
        
        // Output file naming
        const string CSV_FILENAME = "section_map.csv";
        const string CHI2_PDF_FILENAME = "chi2_scans.pdf";
        const string INTERPOLATED_SUFFIX = "_interpolated";
        
        // ========================================================================
        // CONFIGURATION SUMMARY
        // ========================================================================
        // This section helps verify your configuration intent.
        // 
        // FOR THE RECOMMENDED PROFILE, verify these settings:
        //   ✓ USE_THREE_STAGE_OPTIMIZATION = true  (separates angular/energy systematics)
        //   ✓ ENERGY_SMEARING_HISTOGRAM = HIST_BOTH      (uses both observables)
        //   ✓ ENABLE_POSITION_SMEARING = true    (calibrates angular resolution)
        //   ✓ W_MPI0 and W_MMISS set intentionally for your strategy
        //
        // CURRENT DEFAULTS IN THIS FILE:
        //   - ENERGY_SMEARING_HISTOGRAM = HIST_BOTH
        //   - ENABLE_POSITION_SMEARING = true
        //   - W_MPI0 = 1.5, W_MMISS = 0.0, W_MPGG2 = 1.0
        //
        // AVOID THIS CONFIGURATION (poor M_γγ agreement):
        //   ✗ ENERGY_SMEARING_HISTOGRAM = HIST_MMISS_ONLY (M_miss only)
        //   ✗ ENABLE_POSITION_SMEARING = false
        //   See MISSING_MASS_FITTING_ISSUE.txt for detailed physics explanation
        // ========================================================================

        inline const char* histogram_mode_label() {
            if (ENERGY_SMEARING_HISTOGRAM == HIST_MPI0_ONLY) return "M_gg only";
            if (ENERGY_SMEARING_HISTOGRAM == HIST_MMISS_ONLY) return "M_miss only";
            if (ENERGY_SMEARING_HISTOGRAM == HIST_MPGG2_ONLY) return "(p_target + #gamma#gamma)^2 only";
            return "M_gg + M_miss + (p_target + #gamma#gamma)^2";
        }

        inline bool stage2_uses_mmiss() {
            return (ENERGY_SMEARING_HISTOGRAM == HIST_MMISS_ONLY ||
                    ENERGY_SMEARING_HISTOGRAM == HIST_BOTH);
        }

        inline void print_configuration_summary() {
            cout << "\n==== Active Configuration ====\n";
            cout << "Optimization mode: "
                << (USE_THREE_STAGE_OPTIMIZATION ? "staged (1-3 + optional 4)" : "simultaneous") << "\n";
            cout << "Stage-2 observable: " << histogram_mode_label() << "\n";
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
            if (ENABLE_ELECTRON_MOMENTUM_SCALING && stage2_uses_mmiss()) {
                if (ENABLE_FINAL_GLOBAL_PE_SCALE_STAGE4) {
                    cout << "p_e_scale mode: per-section in Pass-2 sweeps + final global Stage-4 refinement\n";
                } else {
                    cout << "p_e_scale mode: per-section in Pass-2 sweeps only (final global Stage-4 disabled)\n";
                }
            }
            cout << "Weights (W_MPI0, W_MMISS, W_MPGG2): " << W_MPI0 << ", " << W_MMISS << ", " << W_MPGG2 << "\n";
            cout << "Exclusive gating in weights: "
                << (APPLY_IS_EXCLUSIVE_SELECTION ? "enabled" : "disabled (all events)") << "\n";
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
                    cout << "  mu(E)=mu0*((E0/(E+Eshift))^b+1) with E0=" << MU_ENERGY_E0_GEV
                        << " GeV, Eshift=" << MU_ENERGY_SHIFT_GEV
                        << " GeV, b=" << MU_ENERGY_B_EXPONENT << "\n";
            }
            cout << "Pass-2 mode: "
                << (ENABLE_ITERATIVE_COUPLED_PASS2 ? "iterative coupled" : "single sweep") << "\n";
            if (ENABLE_ITERATIVE_COUPLED_PASS2) {
                cout << "Pass-2 max sweeps: " << COUPLED_PASS2_MAX_SWEEPS
                    << "  conv(mu,sigma,sigma_pos)=(" << COUPLED_CONV_MU
                    << ", " << COUPLED_CONV_SIGMA
                    << ", " << COUPLED_CONV_SIGMA_POS << ")\n";
                cout << "Pass-2 sweep mode: iterative seeded refine (uses previous sweep parameters)\n";
            }

            if (ENABLE_ELECTRON_MOMENTUM_SCALING && !stage2_uses_mmiss()) {
                cout << "Note: p_e_scale is disabled when Stage-2 mode does not use M_miss\n";
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
        if (h_react_z < -4.0 || h_react_z > 4.0) return false;
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
        // Pass-1 (first-pass) calibration for out-of-section photons, used in pass 2.
        // Set from the pass-1 fit result of whichever center-cell section owns each photon.
        // Defaults: mu=1 (no scale change), sigma=0 (no resolution smearing).
        double mu1_ext = 1.0, sigma1_ext = 0.0, sigma_pos1_ext = 0.0;
        double mu2_ext = 1.0, sigma2_ext = 0.0, sigma_pos2_ext = 0.0;
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
        if (Config::USE_SIMPLE_STOCHASTIC_MODEL) {
            // Simple stochastic: σ_E = σ × √E
            return sigma * sqrt(E_scaled);
        } else {
            // Full 3-term: σ_E = σ × E × √((A/√E)² + (B/E)² + C²)
            double A_sq = res_A * res_A;
            double B_sq = res_B * res_B;
            double C_sq = res_C * res_C;
            double sigma_rel_sq = A_sq / E_scaled + B_sq / (E_scaled * E_scaled) + C_sq;
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
    inline void smearPhoton(double E, double x, double y,
                            double mu, double sigma, double sigma_pos,
                            double res_A, double res_B, double res_C,
                            TRandom3 &rng,
                            double &E_out, double &x_out, double &y_out) {
        double mu_eff = mu;
        // if (Config::ENABLE_ENERGY_DEPENDENT_MU) {
        //     double E_safe = max(E, Config::MU_ENERGY_MIN_GEV);
        //     double mu_shape = pow(Config::MU_ENERGY_E0_GEV /
        //                   (E_safe + Config::MU_ENERGY_SHIFT_GEV),
        //                   Config::MU_ENERGY_B_EXPONENT) + 1.0;
        //     mu_eff = mu * mu_shape;
        // }

        if (Config::ENABLE_ENERGY_DEPENDENT_MU) {
            double E_safe = std::max(E, Config::MU_ENERGY_MIN_GEV);
            mu_eff = Config::MU_ENERGY_A
                + Config::MU_ENERGY_B * E_safe
                + Config::MU_ENERGY_C * std::log(E_safe);
        }
        double E_sc = mu_eff * E;

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
            double fwhmE = computeEnergyResolution(E_sc, sigma, res_A, res_B, res_C);
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
                double sigE = computeEnergyResolution(E_sc, sigma, res_A, res_B, res_C);
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
            double sigma_pos_eff = sigma_pos;
            if (Config::ENABLE_ENERGY_DEPENDENT_SIGMA_POS) {
                double E_for_pos = max(E_sc, Config::SIGMA_POS_ENERGY_MIN_GEV);
                sigma_pos_eff = sigma_pos * sqrt(Config::SIGMA_POS_ENERGY_E0_GEV / E_for_pos);
            }
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

    inline void fillSmearedHistogramsAtParams(const vector<ClusterPair> &simEvents,
                                            TH1D &hmpi0,
                                            TH1D &hmmiss,
                                            TH1D &hmpgg2,
                                            double mu,
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

    // compute chi2 using per-bin errors: sum ( (s-d)^2 / (sigma_s^2 + sigma_d^2) )
    // Requires Sumw2() to be called on both histograms before filling so that
    // GetBinError() returns sqrt(sum_of_weights^2) rather than sqrt(N).
    // This correctly handles weighted data and weighted simulation histograms.
    // NOTE: If both errors are zero in a populated bin, use a Poisson-like fallback
    // variance max(1, s+d) to avoid both divide-by-zero and overly weak penalties.
    double computeChi2FromHist(const TH1D &hsim, const TH1D &hdata) {
        if (hsim.GetNbinsX() != hdata.GetNbinsX() ||
            fabs(hsim.GetXaxis()->GetXmin() - hdata.GetXaxis()->GetXmin()) > 1e-9 ||
            fabs(hsim.GetXaxis()->GetXmax() - hdata.GetXaxis()->GetXmax()) > 1e-9) {
            cerr << "Histogram binning mismatch in computeChi2FromHist" << endl;
            return 1e300;
        }
        double chi2 = 0.;
        int nb = hsim.GetNbinsX();

        if (Config::USE_BAKER_COUSINS_CHI2) {
            // Baker-Cousins log-likelihood ratio:
            //   chi2 = 2 * sum_i [ s_i - d_i + d_i * ln(d_i / s_i) ]
            // Handles low-count bins correctly; unbiased; follows a chi-squared
            // distribution more closely than Pearson chi2.  Does not depend on
            // simulation error bars, eliminating stochastic smearing noise.
            for (int i = 1; i <= nb; ++i) {
                double s = hsim.GetBinContent(i);
                double d = hdata.GetBinContent(i);
                if (d > 0.0 && s > 0.0) {
                    chi2 += 2.0 * (s - d + d * log(d / s));
                } else if (s > 0.0) {
                    chi2 += 2.0 * s;
                }
                // d > 0, s <= 0: pathological after area normalisation — skip to
                // avoid log(0).  Both zero: no contribution (correct).
            }
        } else {
            // Pearson chi2 with combined sim+data errors:
            //   chi2 = sum (s-d)^2 / (sigma_s^2 + sigma_d^2)
            for (int i = 1; i <= nb; ++i) {
                double s    = hsim.GetBinContent(i);
                double d    = hdata.GetBinContent(i);
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
        return chi2;
    }

    inline int countInformativeBinsData(const TH1D &h) {
        int n = 0;
        for (int i = 1; i <= h.GetNbinsX(); ++i) {
            if (h.GetBinContent(i) > 0.0) ++n;
        }
        return max(1, n);
    }

    // ============================================================================
    // COMBINED CHI2 EVALUATION
    // ============================================================================
    // COMBINED CHI2: Simultaneously evaluate both M_γγ and M_miss observables
    // This accounts for correlations between observables (both depend on photon energies)
    // Returns: w_mpi0 * χ²_M_γγ + w_mmiss * χ²_M_miss
    double eval_chi2_combined(double mu, double sigma, double sigma_pos, double p_e_scale,
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
        
        // Create temporary histograms for both observables
        TH1D hsim_mpi0("hsim_mpi0_comb_tmp","smeared sim invariant mass (tmp)", 
                    hdata_mpi0.GetNbinsX(), hdata_mpi0.GetXaxis()->GetXmin(), hdata_mpi0.GetXaxis()->GetXmax());
        TH1D hsim_mmiss("hsim_mmiss_comb_tmp","smeared sim missing mass (tmp)", 
                        hdata_mmiss.GetNbinsX(), hdata_mmiss.GetXaxis()->GetXmin(), hdata_mmiss.GetXaxis()->GetXmax());
        TH1D hsim_mpgg2("hsim_mpgg2_comb_tmp","smeared sim (p_{target}+#gamma#gamma)^{2} (tmp)",
                        hdata_mpgg2.GetNbinsX(), hdata_mpgg2.GetXaxis()->GetXmin(), hdata_mpgg2.GetXaxis()->GetXmax());
        hsim_mpi0.Sumw2();
        hsim_mmiss.Sumw2();
        hsim_mpgg2.Sumw2();

        fillSmearedHistogramsAtParams(simEvents,
                                    hsim_mpi0, hsim_mmiss, hsim_mpgg2,
                                    mu, sigma, sigma_pos,
                                    p_e_scale,
                                    rng, Nsmear,
                                    res_A, res_B, res_C);

        // Scale both histograms to match data
        double integral_sim_mpi0 = hsim_mpi0.Integral();
        double integral_data_mpi0 = hdata_mpi0.Integral();
        if (integral_sim_mpi0 <= 0 || integral_data_mpi0 <= 0) return 1e300;
        hsim_mpi0.Scale(integral_data_mpi0 / integral_sim_mpi0);
        
        double integral_sim_mmiss = hsim_mmiss.Integral();
        double integral_data_mmiss = hdata_mmiss.Integral();
        if (integral_sim_mmiss <= 0 || integral_data_mmiss <= 0) return 1e300;
        hsim_mmiss.Scale(integral_data_mmiss / integral_sim_mmiss);

        double integral_sim_mpgg2 = hsim_mpgg2.Integral();
        double integral_data_mpgg2 = hdata_mpgg2.Integral();
        if (integral_sim_mpgg2 <= 0 || integral_data_mpgg2 <= 0) return 1e300;
        hsim_mpgg2.Scale(integral_data_mpgg2 / integral_sim_mpgg2);
        
        // Calculate chi2 for both observables
        double chi2_mpi0 = computeChi2FromHist(hsim_mpi0, hdata_mpi0);
        double chi2_mmiss = computeChi2FromHist(hsim_mmiss, hdata_mmiss);
        double chi2_mpgg2 = computeChi2FromHist(hsim_mpgg2, hdata_mpgg2);
        
        // Combined chi2 with weights
        return w_mpi0 * chi2_mpi0 + w_mmiss * chi2_mmiss + w_mpgg2 * chi2_mpgg2;
    }

    // ============================================================================
    // INVARIANT MASS CHI2 EVALUATION (Stage 1 - Position Smearing)
    // ============================================================================
    // Evaluate chi2 using ONLY invariant mass M_γγ
    // This is used in Stage 1 to fit position smearing parameters
    double eval_chi2_mpi0_only(double mu, double sigma, double sigma_pos,
                            const vector<ClusterPair> &simEvents,
                            const TH1D &hdata_mpi0,
                            TRandom3 &rng,
                            int Nsmear,
                            double res_A, double res_B, double res_C) {
        
        TH1D hsim_mpi0("hsim_mpi0_stage1_tmp","smeared sim invariant mass (tmp)",
                    hdata_mpi0.GetNbinsX(), hdata_mpi0.GetXaxis()->GetXmin(), hdata_mpi0.GetXaxis()->GetXmax());
        hsim_mpi0.Sumw2();

        rng.SetSeed(Config::SMEAR_DETERMINISTIC_SEED);
        const double z_nps = nps::kDefaultZ_NPS_cm;

        for (const auto &ev : simEvents) {
            if (ev.e1 <= 0 || ev.e2 <= 0) continue;
            double mu1   = ev.photon1_in_section ? mu        : ev.mu1_ext;
            double sig1  = ev.photon1_in_section ? sigma     : ev.sigma1_ext;
            double spos1 = ev.photon1_in_section ? sigma_pos : ev.sigma_pos1_ext;
            double mu2   = ev.photon2_in_section ? mu        : ev.mu2_ext;
            double sig2  = ev.photon2_in_section ? sigma     : ev.sigma2_ext;
            double spos2 = ev.photon2_in_section ? sigma_pos : ev.sigma_pos2_ext;

            double event_weight = ev.weight / Nsmear;
            for (int k = 0; k < Nsmear; ++k) {
                double E1_sm, x1_sm, y1_sm, E2_sm, x2_sm, y2_sm;
                smearPhoton(ev.e1, ev.x1, ev.y1, mu1, sig1, spos1,
                            res_A, res_B, res_C, rng, E1_sm, x1_sm, y1_sm);
                smearPhoton(ev.e2, ev.x2, ev.y2, mu2, sig2, spos2,
                            res_A, res_B, res_C, rng, E2_sm, x2_sm, y2_sm);
                apply_mgg_linear_pair_energy_correction(E1_sm, E2_sm, x1_sm, y1_sm, x2_sm, y2_sm);
                PhotonMomentum p1 = computePhotonMomentum(E1_sm, x1_sm, y1_sm, z_nps);
                PhotonMomentum p2 = computePhotonMomentum(E2_sm, x2_sm, y2_sm, z_nps);
                double mass = computeInvariantMass(E1_sm, p1.px, p1.py, p1.pz,
                                                E2_sm, p2.px, p2.py, p2.pz);
                hsim_mpi0.Fill(mass, event_weight);
            }
        }

        double integral_sim = hsim_mpi0.Integral();
        double integral_data = hdata_mpi0.Integral();
        if (integral_sim <= 0 || integral_data <= 0) return 1e300;
        hsim_mpi0.Scale(integral_data / integral_sim);

        return computeChi2FromHist(hsim_mpi0, hdata_mpi0);
    }

    // ============================================================================
    // MISSING MASS CHI2 EVALUATION (Stage 2 - Energy Smearing)
    // ============================================================================
    // Evaluate chi2 using ONLY missing mass M_miss
    // This is used in Stage 2 to fit energy smearing parameters with fixed position smearing
    double eval_chi2_mmiss_only(double mu, double sigma, double sigma_pos, double p_e_scale,
                                const vector<ClusterPair> &simEvents,
                                const TH1D &hdata_mmiss,
                                TRandom3 &rng,
                                int Nsmear,
                                double res_A, double res_B, double res_C) {
        
        TH1D hsim_mmiss("hsim_mmiss_stage2_tmp","smeared sim missing mass (tmp)",
                        hdata_mmiss.GetNbinsX(), hdata_mmiss.GetXaxis()->GetXmin(), hdata_mmiss.GetXaxis()->GetXmax());
        hsim_mmiss.Sumw2();

        rng.SetSeed(Config::SMEAR_DETERMINISTIC_SEED);

        for (const auto &ev : simEvents) {
            if (ev.e1 <= 0 || ev.e2 <= 0) continue;
            double mu1   = ev.photon1_in_section ? mu        : ev.mu1_ext;
            double sig1  = ev.photon1_in_section ? sigma     : ev.sigma1_ext;
            double spos1 = ev.photon1_in_section ? sigma_pos : ev.sigma_pos1_ext;
            double mu2   = ev.photon2_in_section ? mu        : ev.mu2_ext;
            double sig2  = ev.photon2_in_section ? sigma     : ev.sigma2_ext;
            double spos2 = ev.photon2_in_section ? sigma_pos : ev.sigma_pos2_ext;

            double event_weight = ev.weight / Nsmear;
            for (int k = 0; k < Nsmear; ++k) {
                double E1_sm, x1_sm, y1_sm, E2_sm, x2_sm, y2_sm;
                smearPhoton(ev.e1, ev.x1, ev.y1, mu1, sig1, spos1,
                            res_A, res_B, res_C, rng, E1_sm, x1_sm, y1_sm);
                smearPhoton(ev.e2, ev.x2, ev.y2, mu2, sig2, spos2,
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

        double integral_sim = hsim_mmiss.Integral();
        double integral_data = hdata_mmiss.Integral();
        if (integral_sim <= 0 || integral_data <= 0) return 1e300;
        hsim_mmiss.Scale(integral_data / integral_sim);

        return computeChi2FromHist(hsim_mmiss, hdata_mmiss);
    }

    // ============================================================================
    // TARGET+2γ MASS-SQUARED CHI2 EVALUATION (Stage 2 - Energy Smearing)
    // ============================================================================
    double eval_chi2_mpgg2_only(double mu, double sigma, double sigma_pos,
                                const vector<ClusterPair> &simEvents,
                                const TH1D &hdata_mpgg2,
                                TRandom3 &rng,
                                int Nsmear,
                                double res_A, double res_B, double res_C) {

        TH1D hsim_mpgg2("hsim_mpgg2_stage2_tmp","smeared sim (p_{target}+#gamma#gamma)^{2} (tmp)",
                        hdata_mpgg2.GetNbinsX(), hdata_mpgg2.GetXaxis()->GetXmin(), hdata_mpgg2.GetXaxis()->GetXmax());
        hsim_mpgg2.Sumw2();

        rng.SetSeed(Config::SMEAR_DETERMINISTIC_SEED);
        const double z_nps = nps::kDefaultZ_NPS_cm;

        for (const auto &ev : simEvents) {
            if (ev.e1 <= 0 || ev.e2 <= 0) continue;
            double mu1   = ev.photon1_in_section ? mu        : ev.mu1_ext;
            double sig1  = ev.photon1_in_section ? sigma     : ev.sigma1_ext;
            double spos1 = ev.photon1_in_section ? sigma_pos : ev.sigma_pos1_ext;
            double mu2   = ev.photon2_in_section ? mu        : ev.mu2_ext;
            double sig2  = ev.photon2_in_section ? sigma     : ev.sigma2_ext;
            double spos2 = ev.photon2_in_section ? sigma_pos : ev.sigma_pos2_ext;

            double event_weight = ev.weight / Nsmear;
            for (int k = 0; k < Nsmear; ++k) {
                double E1_sm, x1_sm, y1_sm, E2_sm, x2_sm, y2_sm;
                smearPhoton(ev.e1, ev.x1, ev.y1, mu1, sig1, spos1,
                            res_A, res_B, res_C, rng, E1_sm, x1_sm, y1_sm);
                smearPhoton(ev.e2, ev.x2, ev.y2, mu2, sig2, spos2,
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

        double integral_sim = hsim_mpgg2.Integral();
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
    double eval_chi2_selected(double mu, double sigma, double sigma_pos, double p_e_scale,
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
            // Use invariant mass only
            return eval_chi2_mpi0_only(mu, sigma, sigma_pos, simEvents, 
                                    hdata_mpi0, rng, Nsmear, res_A, res_B, res_C);
        } else if (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_MMISS_ONLY) {
            // Use missing mass only
            return eval_chi2_mmiss_only(mu, sigma, sigma_pos, p_e_scale, simEvents, 
                                        hdata_mmiss, rng, Nsmear, res_A, res_B, res_C);
        } else if (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_MPGG2_ONLY) {
            // Use (p_target + gamma+gamma)^2 only
            return eval_chi2_mpgg2_only(mu, sigma, sigma_pos, simEvents,
                                        hdata_mpgg2, rng, Nsmear, res_A, res_B, res_C);
        } else {
            // Use both histograms (combined chi2)
            return eval_chi2_combined(mu, sigma, sigma_pos, p_e_scale, simEvents, 
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
                                            double mu, double sigma, double sigma_pos,
                                            double p_e_scale,
                                            TRandom3 &rng,
                                            int Nsmear,
                                            double res_A, double res_B, double res_C) {
        rng.SetSeed(Config::SMEAR_DETERMINISTIC_SEED);
        const double z_nps = nps::kDefaultZ_NPS_cm;
        for (const auto &ev : simEvents) {
            if (ev.e1 <= 0 || ev.e2 <= 0) continue;

            double mu1   = ev.photon1_in_section ? mu        : ev.mu1_ext;
            double sig1  = ev.photon1_in_section ? sigma     : ev.sigma1_ext;
            double spos1 = ev.photon1_in_section ? sigma_pos : ev.sigma_pos1_ext;
            double mu2   = ev.photon2_in_section ? mu        : ev.mu2_ext;
            double sig2  = ev.photon2_in_section ? sigma     : ev.sigma2_ext;
            double spos2 = ev.photon2_in_section ? sigma_pos : ev.sigma_pos2_ext;

            double event_weight = ev.weight / Nsmear;
            for (int k = 0; k < Nsmear; ++k) {
                double E1_sm, x1_sm, y1_sm, E2_sm, x2_sm, y2_sm;
                smearPhoton(ev.e1, ev.x1, ev.y1, mu1, sig1, spos1,
                            res_A, res_B, res_C, rng, E1_sm, x1_sm, y1_sm);
                smearPhoton(ev.e2, ev.x2, ev.y2, mu2, sig2, spos2,
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
        vector<double> mu_values_;
        vector<double> sigma_values_;
        vector<double> sigma_pos_values_;

        // GP/Kriging-style interpolation cache
        mutable bool gp_model_ready_;
        mutable bool gp_model_valid_;
        mutable vector<double> gp_alpha_mu_;
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

        void compute_bilinear_fallback(double x, double y, double &mu, double &sigma, double &sigma_pos) const {
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

            double mu00 = mu_values_[idx00], mu10 = mu_values_[idx10], mu01 = mu_values_[idx01], mu11 = mu_values_[idx11];
            mu = (1 - tx) * (1 - ty) * mu00 + tx * (1 - ty) * mu10 + (1 - tx) * ty * mu01 + tx * ty * mu11;

            double s00 = sigma_values_[idx00], s10 = sigma_values_[idx10], s01 = sigma_values_[idx01], s11 = sigma_values_[idx11];
            sigma = (1 - tx) * (1 - ty) * s00 + tx * (1 - ty) * s10 + (1 - tx) * ty * s01 + tx * ty * s11;

            double p00 = sigma_pos_values_[idx00], p10 = sigma_pos_values_[idx10], p01 = sigma_pos_values_[idx01], p11 = sigma_pos_values_[idx11];
            sigma_pos = (1 - tx) * (1 - ty) * p00 + tx * (1 - ty) * p10 + (1 - tx) * ty * p01 + tx * ty * p11;
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

            vector<double> alpha_mu, alpha_sigma, alpha_sigma_pos;
            bool ok_mu = solve_linear_system(K, mu_values_, alpha_mu, nsec);
            bool ok_sigma = solve_linear_system(K, sigma_values_, alpha_sigma, nsec);
            bool ok_sigma_pos = solve_linear_system(K, sigma_pos_values_, alpha_sigma_pos, nsec);

            if (!(ok_mu && ok_sigma && ok_sigma_pos)) return;

            gp_alpha_mu_ = alpha_mu;
            gp_alpha_sigma_ = alpha_sigma;
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
            mu_values_.resize(nsec, 1.0);      // default values
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
        
        // Set fitted parameters for a specific section
        void setParams(int ix, int iy, double mu, double sigma, double sigma_pos = 0.0) {
            int idx = getIndex(ix, iy);
            mu_values_[idx] = mu;
            sigma_values_[idx] = sigma;
            sigma_pos_values_[idx] = sigma_pos;
            gp_model_ready_ = false;
            gp_model_valid_ = false;
        }
        
        // Get interpolated parameters at position (x, y) using bilinear interpolation
        void getInterpolatedParams(double x, double y, double &mu, double &sigma) const {
            double sigma_pos_dummy = 0.0;
            getInterpolatedParams(x, y, mu, sigma, sigma_pos_dummy);
        }

        // Get interpolated parameters at position (x, y) using bilinear interpolation
        void getInterpolatedParams(double x, double y, double &mu, double &sigma, double &sigma_pos) const {
            compute_bilinear_fallback(x, y, mu, sigma, sigma_pos);
            sigma = std::max(1e-6, sigma);
            sigma_pos = std::max(0.0, sigma_pos);
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
                
                if (tokens.size() < 12) continue;
                
                int ix = atoi(tokens[0].c_str());
                int iy = atoi(tokens[1].c_str());
                
                // Check if fit succeeded (has valid parameters)
                if (tokens[8].empty() || tokens[9].empty() || tokens[10].empty()) {
                    // No valid fit for this section - keep defaults
                    continue;
                }
                
                double mu = atof(tokens[8].c_str());
                double sigma = atof(tokens[9].c_str());
                double sigma_pos = atof(tokens[10].c_str());
                
                setParams(ix, iy, mu, sigma, sigma_pos);
                count++;
            }
            
            file.close();
            cout << "Loaded calibration for " << count << " sections from " << filename << endl;
            return true;
        }
        
        // Save interpolated map to 2D histogram (for visualization)
        void saveAsHistogram(const string &filename) const {
            TFile fout(filename.c_str(), "RECREATE");
            
            int nbinsx = 100;
            int nbinsy = 100;
            
            TH2D *h_mu = new TH2D("h_mu_interp", "Interpolated #mu map;x [cm];y [cm]", 
                                nbinsx, x_min_, x_max_, nbinsy, y_min_, y_max_);
            TH2D *h_sigma = new TH2D("h_sigma_interp", "Interpolated #sigma map;x [cm];y [cm]", 
                                    nbinsx, x_min_, x_max_, nbinsy, y_min_, y_max_);
            TH2D *h_sigma_pos = new TH2D("h_sigma_pos_interp", "Interpolated #sigma_{pos} map;x [cm];y [cm]", 
                            nbinsx, x_min_, x_max_, nbinsy, y_min_, y_max_);
            
            for (int ix = 1; ix <= nbinsx; ++ix) {
                for (int iy = 1; iy <= nbinsy; ++iy) {
                    double x = h_mu->GetXaxis()->GetBinCenter(ix);
                    double y = h_mu->GetYaxis()->GetBinCenter(iy);
                    
                    double mu, sigma, sigma_pos;
                    getInterpolatedParams(x, y, mu, sigma, sigma_pos);
                    
                    h_mu->SetBinContent(ix, iy, mu);
                    h_sigma->SetBinContent(ix, iy, sigma);
                    h_sigma_pos->SetBinContent(ix, iy, sigma_pos);
                }
            }
            
            h_mu->Write();
            h_sigma->Write();
            h_sigma_pos->Write();
            fout.Close();
            
            cout << "Saved interpolated calibration maps to " << filename << endl;
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

    // Structure to store fit results
    struct FitResult { 
        double mu, sigma, sigma_pos, p_e_scale, chi2;
        // 3-term model parameters (only used if USE_SIMPLE_STOCHASTIC_MODEL = false)
        double res_A, res_B, res_C;
        int ndf;
        Chi2ScanData scan_data;
        double chi2_per_ndf() const { return (ndf > 0) ? chi2 / ndf : 0.0; }
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
        const double mu0 = res.mu;
        const double sigma0 = res.sigma;
        const double sigma_pos0 = res.sigma_pos;
        const double p_e_scale0 = res.p_e_scale;

        const int mu_scan_points = max(2, Config::VIS_MU_SCAN_POINTS);
        const int sigma_scan_points = max(2, Config::VIS_SIGMA_SCAN_POINTS);
        const int slice_points = max(2, Config::VIS_SLICE_POINTS);
        const int sigma_pos_scan_points = max(2, Config::VIS_SIGMA_POS_POINTS);

        const double mu_scan_step = (Config::MU_MAX - Config::MU_MIN) / (mu_scan_points - 1);
        const double sigma_scan_step = (Config::SIGMA_MAX - Config::SIGMA_MIN) / (sigma_scan_points - 1);

        for (int i_mu = 0; i_mu < mu_scan_points; ++i_mu) {
            double mu = Config::MU_MIN + i_mu * mu_scan_step;
            for (int i_sig = 0; i_sig < sigma_scan_points; ++i_sig) {
                double sig = Config::SIGMA_MIN + i_sig * sigma_scan_step;
                double chi2 = eval_chi2_selected(mu, sig, sigma_pos0, p_e_scale0, simEvents,
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

        const double mu_slice_step = (Config::MU_MAX - Config::MU_MIN) / (slice_points - 1);
        for (int i = 0; i < slice_points; ++i) {
            double mu = Config::MU_MIN + i * mu_slice_step;
            double chi2 = eval_chi2_selected(mu, sigma0, sigma_pos0, p_e_scale0, simEvents,
                                            hdata_mpi0, hdata_mmiss, hdata_mpgg2, rng, Nsmear,
                                            Config::RESOLUTION_A_DEFAULT,
                                            Config::RESOLUTION_B_DEFAULT,
                                            Config::RESOLUTION_C_DEFAULT,
                                            Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);
            res.scan_data.mu_values.push_back(mu);
            res.scan_data.chi2_vs_mu.push_back(chi2);
        }

        const double sigma_slice_step = (Config::SIGMA_MAX - Config::SIGMA_MIN) / (slice_points - 1);
        for (int i = 0; i < slice_points; ++i) {
            double sigma = Config::SIGMA_MIN + i * sigma_slice_step;
            double chi2 = eval_chi2_selected(mu0, sigma, sigma_pos0, p_e_scale0, simEvents,
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
                double chi2 = eval_chi2_selected(mu0, sigma0, s_pos, p_e_scale0, simEvents,
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
        double mu;
        double sigma;
        double sigma_pos;
        double p_e_scale;
        double chi2;
        bool minimized;
    };

    template <typename Chi2Fn>
    MigradRefineResult run_migrad_refinement(Chi2Fn chi2_fn,
                                            double mu_seed,
                                            double sigma_seed,
                                            double sigma_pos_seed,
                                            double p_e_scale_seed,
                                            bool fit_mu,
                                            bool fit_sigma,
                                            bool fit_sigma_pos,
                                            bool fit_p_e_scale) {
        MigradRefineResult out{mu_seed, sigma_seed, sigma_pos_seed, p_e_scale_seed,
                            chi2_fn(mu_seed, sigma_seed, sigma_pos_seed, p_e_scale_seed),
                            false};

        if (!Config::USE_MIGRAD_REFINEMENT) return out;
        if (!(fit_mu || fit_sigma || fit_sigma_pos || fit_p_e_scale)) return out;

        unique_ptr<ROOT::Math::Minimizer> minimizer(
            ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
        if (!minimizer) return out;

        minimizer->SetMaxFunctionCalls(Config::MIGRAD_MAX_FUNCTION_CALLS);
        minimizer->SetMaxIterations(Config::MIGRAD_MAX_ITERATIONS);
        minimizer->SetTolerance(Config::MIGRAD_TOLERANCE);
        minimizer->SetStrategy(Config::MIGRAD_STRATEGY);
        minimizer->SetPrintLevel(0);

        auto fcn = [&](const double *x) {
            const double mu = x[0];
            const double sigma = x[1];
            const double sigma_pos = x[2];
            const double p_e_scale = x[3];
            return chi2_fn(mu, sigma, sigma_pos, p_e_scale);
        };
        ROOT::Math::Functor functor(fcn, 4);
        minimizer->SetFunction(functor);

        // Always define all 4 coordinates; fixed dimensions are held constant.
        if (fit_mu) {
            minimizer->SetLimitedVariable(0, "mu", mu_seed,
                                        Config::MIGRAD_STEP_MU,
                                        Config::MU_MIN, Config::MU_MAX);
        } else {
            minimizer->SetFixedVariable(0, "mu", mu_seed);
        }

        double sigma_lower = max(1e-6, Config::SIGMA_MIN);
        if (fit_sigma) {
            minimizer->SetLimitedVariable(1, "sigma", sigma_seed,
                                        Config::MIGRAD_STEP_SIGMA,
                                        sigma_lower, Config::SIGMA_MAX);
        } else {
            minimizer->SetFixedVariable(1, "sigma", sigma_seed);
        }

        if (fit_sigma_pos) {
            minimizer->SetLimitedVariable(2, "sigma_pos", sigma_pos_seed,
                                        Config::MIGRAD_STEP_SIGMA_POS,
                                        Config::SIGMA_POS_MIN, Config::SIGMA_POS_MAX);
        } else {
            minimizer->SetFixedVariable(2, "sigma_pos", sigma_pos_seed);
        }

        if (fit_p_e_scale) {
            minimizer->SetLimitedVariable(3, "p_e_scale", p_e_scale_seed,
                                        Config::MIGRAD_STEP_PE_SCALE,
                                        Config::PE_SCALE_MIN, Config::PE_SCALE_MAX);
        } else {
            minimizer->SetFixedVariable(3, "p_e_scale", p_e_scale_seed);
        }

        bool ok = minimizer->Minimize();
        const double *xbest = minimizer->X();
        if (!xbest) return out;

        double mu_best = xbest[0];
        double sigma_best = xbest[1];
        double sigma_pos_best = xbest[2];
        double p_e_scale_best = xbest[3];
        double chi2_best = chi2_fn(mu_best, sigma_best, sigma_pos_best, p_e_scale_best);

        if (ok && chi2_best < out.chi2) {
            out.mu = mu_best;
            out.sigma = sigma_best;
            out.sigma_pos = sigma_pos_best;
            out.p_e_scale = p_e_scale_best;
            out.chi2 = chi2_best;
            out.minimized = true;
        }
        return out;
    }

    // ============================================================================
    // FITTING FUNCTION
    // ============================================================================
    FitResult fit_section(const vector<ClusterPair> &simEvents,
                        const TH1D &hdata_mpi0,     // Invariant mass M_γγ
                        const TH1D &hdata_mmiss,    // Missing mass M_miss
                        const TH1D &hdata_mpgg2,    // (p_target + gamma+gamma)^2
                        TRandom3 &rng,
                        int Nsmear,
                        double global_p_e_scale = 1.0,
                        bool build_visualization_scans = true) {
        FitResult res; 
        res.mu = 1.0; 
        res.sigma = 0.05;
        res.sigma_pos = 0.0;
        res.p_e_scale = global_p_e_scale;
        res.res_A = Config::RESOLUTION_A_DEFAULT;
        res.res_B = Config::RESOLUTION_B_DEFAULT;
        res.res_C = Config::RESOLUTION_C_DEFAULT;
        res.chi2 = 1e300;
        
        // ndf depends on fitted free parameters per section.
        // A/B/C are fixed configuration constants in this macro (not fitted here).
        int n_params = 2 + (Config::ENABLE_POSITION_SMEARING ? 1 : 0);  // (mu, sigma[, sigma_pos])
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
        // Subtract one degree of freedom per independently area-normalized histogram,
        // since scaling sim to match data integral removes one constraint per histogram.
        int n_norm_hists = (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_BOTH) ? 3 : 1;
        res.ndf = n_bins_chi2 - n_params - n_norm_hists;
        
        // Determine if we're actually doing position smearing (used by both approaches)
        const bool do_position_scan = Config::ENABLE_POSITION_SMEARING;

        // Deterministic evaluators for MIGRAD. Re-seeding the RNG per function call
        // stabilizes the objective seen by the derivative-based minimizer.
        auto chi2_selected_eval = [&](double mu, double sigma, double sigma_pos, double p_e_scale) {
            TRandom3 rng_eval(1234567);
            return eval_chi2_selected(mu, sigma, sigma_pos, p_e_scale, simEvents,
                                    hdata_mpi0, hdata_mmiss, hdata_mpgg2, rng_eval, Nsmear,
                                    Config::RESOLUTION_A_DEFAULT,
                                    Config::RESOLUTION_B_DEFAULT,
                                    Config::RESOLUTION_C_DEFAULT,
                                    Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);
        };
        auto chi2_mpi0_eval = [&](double mu, double sigma, double sigma_pos, double /*p_e_scale*/) {
            TRandom3 rng_eval(223344);
            return eval_chi2_mpi0_only(mu, sigma, sigma_pos, simEvents,
                                    hdata_mpi0, rng_eval, Nsmear,
                                    Config::RESOLUTION_A_DEFAULT,
                                    Config::RESOLUTION_B_DEFAULT,
                                    Config::RESOLUTION_C_DEFAULT);
        };
        
        // ========================================================================
        // OPTIMIZATION APPROACH SELECTION
        // ========================================================================
        
        if (Config::USE_THREE_STAGE_OPTIMIZATION) {
            // THREE-STAGE APPROACH
            // Stage 1 : Fit energy (mu, sigma)        — sigma_pos fixed at 0
            // Stage 2 : Fit position (sigma_pos)      — mu, sigma fixed from Stage 1, M_γγ only
            // Stage 3 : Fit all (mu, sigma, sigma_pos) simultaneously — joint refinement
            cout << "\n==== THREE-STAGE OPTIMIZATION ====" << endl;
            cout << "Stage 1: Energy smearing (mu, sigma) with sigma_pos=0" << endl;
            cout << "Stage 2: Position smearing (sigma_pos) using M_γγ, fixed (mu, sigma) from Stage 1" << endl;
            cout << "Stage 3: Joint refinement (mu, sigma, sigma_pos) using selected observable(s)" << endl;

            string sel_obs = (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_MPI0_ONLY) ? "M_γγ only" :
                            (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_MMISS_ONLY) ? "M_miss only" :
                            (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_MPGG2_ONLY) ? "(p_target + γγ)^2 only" :
                            "M_γγ + M_miss + (p_target + γγ)^2";

            // ====================================================================
            // STAGE 1: FIT ENERGY (mu, sigma) WITH sigma_pos = 0
            // ====================================================================
            cout << "\n--- STAGE 1: Fitting energy (mu, sigma) using " << sel_obs << ", sigma_pos=0 ---" << endl;

            double mu0 = 1.0, sigma0 = 0.05, sigma_pos0 = 0.0, p_e_scale0 = global_p_e_scale;
            double best_chi2_s1 = 1e300;

            int MU_COARSE    = max(15, Config::MU_NSTEPS    / Config::COARSE_GRID_DIVISOR);
            int SIGMA_COARSE = max(15, Config::SIGMA_NSTEPS  / Config::COARSE_GRID_DIVISOR);

            double mu_cs    = (Config::MU_MAX    - Config::MU_MIN)    / (MU_COARSE    - 1);
            double sigma_cs = (Config::SIGMA_MAX - Config::SIGMA_MIN)  / (SIGMA_COARSE - 1);

            for (int i_mu = 0; i_mu < MU_COARSE; ++i_mu) {
                double mu = Config::MU_MIN + i_mu * mu_cs;
                for (int i_sig = 0; i_sig < SIGMA_COARSE; ++i_sig) {
                    double sig = Config::SIGMA_MIN + i_sig * sigma_cs;
                    double chi2 = eval_chi2_selected(mu, sig, 0.0, p_e_scale0, simEvents,
                                                    hdata_mpi0, hdata_mmiss, hdata_mpgg2, rng, Nsmear,
                                                    Config::RESOLUTION_A_DEFAULT,
                                                    Config::RESOLUTION_B_DEFAULT,
                                                    Config::RESOLUTION_C_DEFAULT,
                                                    Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);
                    if (chi2 < best_chi2_s1) { best_chi2_s1 = chi2; mu0 = mu; sigma0 = sig; }
                }
            }
            cout << "Stage 1 coarse min: mu=" << mu0 << " sigma=" << sigma0;
            cout << " chi2=" << best_chi2_s1 << endl;

            // Fine refinement
            {
                double step_mu = 0.005, step_sig = 0.01;
                int iref = 0;
                while ((step_mu >= 0.0001 || step_sig >= 0.0001) && iref < Config::MAX_REFINEMENT_ITERATIONS) {
                    double bmu = mu0, bsig = sigma0, bchi = best_chi2_s1;
                    for (double mu = mu0 - 2*step_mu; mu <= mu0 + 2*step_mu + 1e-15; mu += step_mu)
                    for (double sig = max(1e-6, sigma0 - 2*step_sig); sig <= sigma0 + 2*step_sig + 1e-15; sig += step_sig) {
                        double chi2 = eval_chi2_selected(mu, sig, 0.0, p_e_scale0, simEvents,
                                                        hdata_mpi0, hdata_mmiss, hdata_mpgg2, rng, Nsmear,
                                                        Config::RESOLUTION_A_DEFAULT,
                                                        Config::RESOLUTION_B_DEFAULT,
                                                        Config::RESOLUTION_C_DEFAULT,
                                                        Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);
                        if (chi2 < bchi) { bchi = chi2; bmu = mu; bsig = sig; }
                    }
                    mu0 = bmu; sigma0 = bsig; best_chi2_s1 = bchi;
                    step_mu /= 2.0; step_sig /= 2.0; ++iref;
                }
                if (iref >= Config::MAX_REFINEMENT_ITERATIONS) {
                    cout << "[WARN] Stage 1 refinement reached MAX_REFINEMENT_ITERATIONS="
                        << Config::MAX_REFINEMENT_ITERATIONS << endl;
                }
            }
            cout << "Stage 1 final: mu=" << mu0 << " sigma=" << sigma0;
            cout << " chi2=" << best_chi2_s1 << endl;

            {
                MigradRefineResult mr = run_migrad_refinement(chi2_selected_eval,
                                                            mu0, sigma0, 0.0, p_e_scale0,
                                                            true, true, false, false);
                if (mr.minimized) {
                    mu0 = mr.mu;
                    sigma0 = mr.sigma;
                    best_chi2_s1 = mr.chi2;
                    cout << "Stage 1 MIGRAD refine: mu=" << mu0
                        << " sigma=" << sigma0
                        << " chi2=" << best_chi2_s1 << endl;
                }
            }

            // ====================================================================
            // STAGE 2: FIT POSITION (sigma_pos) USING M_γγ — energy fixed
            // ====================================================================
            cout << "\n--- STAGE 2: Fitting sigma_pos using M_γγ, fixed mu=" << mu0 << " sigma=" << sigma0 << " ---" << endl;

            double best_chi2_s2 = 1e300;

            if (do_position_scan) {
                int POS_COARSE = max(8, Config::SIGMA_POS_NSTEPS / Config::COARSE_GRID_DIVISOR);
                double pos_cs = (Config::SIGMA_POS_MAX - Config::SIGMA_POS_MIN) / (POS_COARSE - 1);
                for (int i_pos = 0; i_pos < POS_COARSE; ++i_pos) {
                    double s_pos = Config::SIGMA_POS_MIN + i_pos * pos_cs;
                    double chi2 = eval_chi2_mpi0_only(mu0, sigma0, s_pos, simEvents,
                                                    hdata_mpi0, rng, Nsmear,
                                                    Config::RESOLUTION_A_DEFAULT,
                                                    Config::RESOLUTION_B_DEFAULT,
                                                    Config::RESOLUTION_C_DEFAULT);
                    if (chi2 < best_chi2_s2) { best_chi2_s2 = chi2; sigma_pos0 = s_pos; }
                }
                cout << "Stage 2 coarse min: sigma_pos=" << sigma_pos0 << " chi2=" << best_chi2_s2 << endl;

                // Fine refinement
                double step_pos = 0.01;
                int iref = 0;
                while (step_pos >= 0.0001 && iref < Config::MAX_REFINEMENT_ITERATIONS) {
                    double bpos = sigma_pos0, bchi = best_chi2_s2;
                    for (double s_pos = max(0.0, sigma_pos0 - 2*step_pos);
                        s_pos <= sigma_pos0 + 2*step_pos + 1e-15; s_pos += step_pos) {
                        double chi2 = eval_chi2_mpi0_only(mu0, sigma0, s_pos, simEvents,
                                                        hdata_mpi0, rng, Nsmear,
                                                        Config::RESOLUTION_A_DEFAULT,
                                                        Config::RESOLUTION_B_DEFAULT,
                                                        Config::RESOLUTION_C_DEFAULT);
                        if (chi2 < bchi) { bchi = chi2; bpos = s_pos; }
                    }
                    sigma_pos0 = bpos; best_chi2_s2 = bchi; step_pos /= 2.0; ++iref;
                }
                if (iref >= Config::MAX_REFINEMENT_ITERATIONS) {
                    cout << "[WARN] Stage 2 refinement reached MAX_REFINEMENT_ITERATIONS="
                        << Config::MAX_REFINEMENT_ITERATIONS << endl;
                }
                cout << "Stage 2 final: sigma_pos=" << sigma_pos0 << " chi2=" << best_chi2_s2 << endl;

                {
                    MigradRefineResult mr = run_migrad_refinement(chi2_mpi0_eval,
                                                                mu0, sigma0, sigma_pos0, p_e_scale0,
                                                                false, false, true, false);
                    if (mr.minimized) {
                        sigma_pos0 = mr.sigma_pos;
                        best_chi2_s2 = mr.chi2;
                        cout << "Stage 2 MIGRAD refine: sigma_pos=" << sigma_pos0
                            << " chi2=" << best_chi2_s2 << endl;
                    }
                }
            } else {
                sigma_pos0 = 0.0;
                cout << "Stage 2: Skipped (position smearing disabled)" << endl;
            }

            // ====================================================================
            // STAGE 3: JOINT REFINEMENT (mu, sigma, sigma_pos) — all free
            // ====================================================================
            cout << "\n--- STAGE 3: Joint refinement (mu, sigma";
            if (do_position_scan) cout << ", sigma_pos";
            cout << ") using " << sel_obs << " ---" << endl;
            cout << "  Starting from Stage 1+2 values: mu=" << mu0 << " sigma=" << sigma0
                << " sigma_pos=" << sigma_pos0 << endl;

            double best_chi2_s3 = 1e300;
            // Evaluate chi2 at the Stage 1+2 starting point
            best_chi2_s3 = eval_chi2_selected(mu0, sigma0, sigma_pos0, p_e_scale0, simEvents,
                                            hdata_mpi0, hdata_mmiss, hdata_mpgg2, rng, Nsmear,
                                            Config::RESOLUTION_A_DEFAULT,
                                            Config::RESOLUTION_B_DEFAULT,
                                            Config::RESOLUTION_C_DEFAULT,
                                            Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);

            {
                double step_mu  = 0.005, step_sig = 0.01;
                double step_pos = do_position_scan ? 0.01 : 0.0;
                int iref = 0;
                while ((step_mu >= 0.0001 || step_sig >= 0.0001 || (do_position_scan && step_pos >= 0.0001)) &&
                    iref < Config::MAX_REFINEMENT_ITERATIONS) {
                    double bmu = mu0, bsig = sigma0, bpos = sigma_pos0, bchi = best_chi2_s3;
                    for (double mu  = mu0  - 2*step_mu;  mu  <= mu0  + 2*step_mu  + 1e-15; mu  += step_mu)
                    for (double sig = max(1e-6, sigma0 - 2*step_sig); sig <= sigma0 + 2*step_sig + 1e-15; sig += step_sig) {
                        double ps0 = do_position_scan ? max(0.0, sigma_pos0 - 2*step_pos) : sigma_pos0;
                        double ps1 = do_position_scan ? sigma_pos0 + 2*step_pos + 1e-15    : sigma_pos0;
                        double pss = do_position_scan ? step_pos : 1.0;
                        for (double s_pos = ps0; s_pos <= ps1; s_pos += pss) {
                            double chi2 = eval_chi2_selected(mu, sig, s_pos, p_e_scale0, simEvents,
                                                            hdata_mpi0, hdata_mmiss, hdata_mpgg2, rng, Nsmear,
                                                            Config::RESOLUTION_A_DEFAULT,
                                                            Config::RESOLUTION_B_DEFAULT,
                                                            Config::RESOLUTION_C_DEFAULT,
                                                            Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);
                            if (chi2 < bchi) { bchi = chi2; bmu = mu; bsig = sig; bpos = s_pos; }
                        }
                    }
                    mu0 = bmu; sigma0 = bsig; sigma_pos0 = bpos; best_chi2_s3 = bchi;
                    step_mu /= 2.0; step_sig /= 2.0;
                    if (do_position_scan) step_pos /= 2.0;
                    ++iref;
                }
                if (iref >= Config::MAX_REFINEMENT_ITERATIONS) {
                    cout << "[WARN] Stage 3 refinement reached MAX_REFINEMENT_ITERATIONS="
                        << Config::MAX_REFINEMENT_ITERATIONS << endl;
                }
            }

            cout << "Stage 3 final: mu=" << mu0 << " sigma=" << sigma0
                << " sigma_pos=" << sigma_pos0 << " cm";
            cout << " chi2=" << best_chi2_s3 << endl;

            {
                MigradRefineResult mr = run_migrad_refinement(chi2_selected_eval,
                                                            mu0, sigma0, sigma_pos0, p_e_scale0,
                                                            true, true, do_position_scan, false);
                if (mr.minimized) {
                    mu0 = mr.mu;
                    sigma0 = mr.sigma;
                    sigma_pos0 = mr.sigma_pos;
                    best_chi2_s3 = mr.chi2;
                    cout << "Stage 3 MIGRAD refine: mu=" << mu0
                        << " sigma=" << sigma0
                        << " sigma_pos=" << sigma_pos0
                        << " chi2=" << best_chi2_s3 << endl;
                }
            }
            if (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_BOTH)
                cout << "  (w_M_γγ=" << Config::W_MPI0
                    << ", w_M_miss=" << Config::W_MMISS
                    << ", w_(p+γγ)^2=" << Config::W_MPGG2 << ")" << endl;

            // Store final results
            res.mu        = mu0;
            res.sigma     = sigma0;
            res.sigma_pos = sigma_pos0;
            res.p_e_scale = p_e_scale0;
            res.chi2      = best_chi2_s3;
            
        } else {
            // ====================================================================
            // SIMULTANEOUS OPTIMIZATION (ORIGINAL APPROACH)
            // ====================================================================
            // Fit all parameters together using selected histogram(s)
            
            string sim_observable = (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_MPI0_ONLY) ? "M_γγ only" :
                    (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_MMISS_ONLY) ? "M_miss only" :
                    (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_MPGG2_ONLY) ? "(p_target + γγ)^2 only" :
                    "M_γγ + M_miss + (p_target + γγ)^2";
            cout << "\n==== SIMULTANEOUS OPTIMIZATION ====" << endl;
            cout << "Fitting all parameters together using " << sim_observable << "..." << endl;
            if (Config::ENERGY_SMEARING_HISTOGRAM == Config::HIST_BOTH) {
                cout << "  Observable weights: w_M_γγ = " << Config::W_MPI0
                    << ", w_M_miss = " << Config::W_MMISS
                    << ", w_(p+γγ)^2 = " << Config::W_MPGG2 << endl;
            }
            
            double best_chi2 = 1e300;
            double mu0 = 1.0, sigma0 = 0.05, sigma_pos0 = 0.0, p_e_scale0 = global_p_e_scale;
            
            // Coarse 3D grid search
            int MU_COARSE = max(15, Config::MU_NSTEPS / Config::COARSE_GRID_DIVISOR);
            int SIGMA_COARSE = max(15, Config::SIGMA_NSTEPS / Config::COARSE_GRID_DIVISOR);
            int SIGMA_POS_COARSE = do_position_scan ? max(5, Config::SIGMA_POS_NSTEPS / Config::COARSE_GRID_DIVISOR) : 1;
            
            double mu_coarse_step = (Config::MU_MAX - Config::MU_MIN) / (MU_COARSE - 1);
            double sigma_coarse_step = (Config::SIGMA_MAX - Config::SIGMA_MIN) / (SIGMA_COARSE - 1);
            double sigma_pos_coarse_step = do_position_scan ? (Config::SIGMA_POS_MAX - Config::SIGMA_POS_MIN) / (SIGMA_POS_COARSE - 1) : 0.0;
            
            cout << "Performing coarse grid search (mu, sigma" << (do_position_scan ? ", sigma_pos" : "")
                << ") with fixed global p_e_scale=" << p_e_scale0 << "..." << endl;
            
            for (int i_mu = 0; i_mu < MU_COARSE; ++i_mu) {
                double mu = Config::MU_MIN + i_mu * mu_coarse_step;
                for (int i_sig = 0; i_sig < SIGMA_COARSE; ++i_sig) {
                    double sig = Config::SIGMA_MIN + i_sig * sigma_coarse_step;
                    for (int i_pos = 0; i_pos < SIGMA_POS_COARSE; ++i_pos) {
                        double s_pos = Config::SIGMA_POS_MIN + i_pos * sigma_pos_coarse_step;
                        
                        double chi2 = eval_chi2_selected(mu, sig, s_pos, p_e_scale0, simEvents,
                                                        hdata_mpi0, hdata_mmiss, hdata_mpgg2, rng, Nsmear,
                                                        Config::RESOLUTION_A_DEFAULT,
                                                        Config::RESOLUTION_B_DEFAULT,
                                                        Config::RESOLUTION_C_DEFAULT,
                                                        Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);

                        if (chi2 < best_chi2) {
                            best_chi2 = chi2;
                            mu0 = mu;
                            sigma0 = sig;
                            sigma_pos0 = s_pos;
                        }
                    }
                }
            }
            cout << "Coarse minimum at mu=" << mu0 << ", sigma=" << sigma0;
            if (do_position_scan) cout << ", sigma_pos=" << sigma_pos0 << " cm";
            cout << ", chi2=" << best_chi2 << endl;
            
            // Fine refinement
            cout << "Refining minimum..." << endl;
            double step_mu = 0.005, step_sigma = 0.01, step_sigma_pos = do_position_scan ? 0.01 : 0.0;
            int refinement_iterations = 0;
            
            while ((step_mu >= 0.0001 || step_sigma >= 0.0001 || (do_position_scan && step_sigma_pos >= 0.0001)) && 
                refinement_iterations < Config::MAX_REFINEMENT_ITERATIONS) {
                double grid_best_mu = mu0, grid_best_sigma = sigma0, grid_best_sigma_pos = sigma_pos0;
                double grid_best_chi2 = best_chi2;
                
                for (double mu = mu0 - 2*step_mu; mu <= mu0 + 2*step_mu + 1e-15; mu += step_mu) {
                    for (double sig = max(1e-6, sigma0 - 2*step_sigma); sig <= sigma0 + 2*step_sigma + 1e-15; sig += step_sigma) {
                        double pos_start = do_position_scan ? max(0.0, sigma_pos0 - 2*step_sigma_pos) : 0.0;
                        double pos_end = do_position_scan ? sigma_pos0 + 2*step_sigma_pos + 1e-15 : 0.0;
                        double pos_step = do_position_scan ? step_sigma_pos : 1.0;
                        
                        for (double s_pos = pos_start; s_pos <= pos_end; s_pos += pos_step) {
                            if (sig <= 0) continue;
                            
                            double chi2 = eval_chi2_selected(mu, sig, s_pos, p_e_scale0, simEvents,
                                                            hdata_mpi0, hdata_mmiss, hdata_mpgg2, rng, Nsmear,
                                                            Config::RESOLUTION_A_DEFAULT,
                                                            Config::RESOLUTION_B_DEFAULT,
                                                            Config::RESOLUTION_C_DEFAULT,
                                                            Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);
                            if (chi2 < grid_best_chi2) {
                                grid_best_chi2 = chi2;
                                grid_best_mu = mu;
                                grid_best_sigma = sig;
                                grid_best_sigma_pos = s_pos;
                            }
                        }
                    }
                }
                mu0 = grid_best_mu;
                sigma0 = grid_best_sigma;
                sigma_pos0 = grid_best_sigma_pos;
                best_chi2 = grid_best_chi2;
                step_mu /= 2.0;
                step_sigma /= 2.0;
                if (do_position_scan) step_sigma_pos /= 2.0;
                refinement_iterations++;
            }
            if (refinement_iterations >= Config::MAX_REFINEMENT_ITERATIONS) {
                cout << "[WARN] Simultaneous refinement reached MAX_REFINEMENT_ITERATIONS="
                    << Config::MAX_REFINEMENT_ITERATIONS << endl;
            }
            cout << "Refined minimum at mu=" << mu0 << ", sigma=" << sigma0;
            if (do_position_scan) cout << ", sigma_pos=" << sigma_pos0 << " cm";
            cout << ", chi2=" << best_chi2 << endl;

            {
                MigradRefineResult mr = run_migrad_refinement(chi2_selected_eval,
                                                            mu0, sigma0, sigma_pos0, p_e_scale0,
                                                            true, true, do_position_scan, false);
                if (mr.minimized) {
                    mu0 = mr.mu;
                    sigma0 = mr.sigma;
                    sigma_pos0 = mr.sigma_pos;
                    best_chi2 = mr.chi2;
                    cout << "Simultaneous MIGRAD refine: mu=" << mu0
                        << ", sigma=" << sigma0;
                    if (do_position_scan) cout << ", sigma_pos=" << sigma_pos0;
                    cout << ", chi2=" << best_chi2 << endl;
                }
            }
            
            // Store final results from simultaneous optimization
            res.mu = mu0;
            res.sigma = sigma0;
            res.sigma_pos = sigma_pos0;
            res.p_e_scale = p_e_scale0;
            res.chi2 = best_chi2;
        }
        
        // ========================================================================
        // GENERATE VISUALIZATION DATA (from final fitted parameters)
        // ========================================================================
        if (build_visualization_scans) {
            cout << "\n--- Generating visualization data ---" << endl;
            generate_visualization_scan_data(simEvents,
                                            hdata_mpi0, hdata_mmiss, hdata_mpgg2,
                                            rng, Nsmear,
                                            res);
        }
        
        return res;
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

    void plot_chi2_scans(const Section& sec, const FitResult& fres, TCanvas* c,
                        const TH1D& hdata,          const TH1D& hsim_final,       const TH1D& hsim_unsmeared,
                        const TH1D& hdata_mmiss,    const TH1D& hsim_mmiss,       const TH1D& hunsmeared_mmiss,
                        const TH1D& hdata_mpgg2,    const TH1D& hsim_mpgg2,       const TH1D& hunsmeared_mpgg2,
                        int n_data, int n_data_selected, int n_sim,
                        int nsmear,
                        double global_p_e_scale,
                        bool use_global_p_e_scale) {
        c->Clear();
        c->Divide(3, 3);  // 3×3 grid: 9 diagnostic panels

        // -------------------------------------------------------------------------
        // Row 1: Observable comparisons
        // -------------------------------------------------------------------------

        // Pad 1 — M_γγ: Data vs Unsmeared vs Smeared
        c->cd(1);
        gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        drawComparisonPad(hdata, hsim_final, hsim_unsmeared,
                        "M_{#gamma#gamma} [GeV/c^{2}]", "mgg");
        TLatex *t_sec = new TLatex(); t_sec->SetNDC(); t_sec->SetTextSize(0.046);
        t_sec->DrawLatex(0.15, 0.93, Form("Section %s  [x: %.1f#rightarrow%.1f  y: %.1f#rightarrow%.1f cm]",
                        sec.name().c_str(), sec.xlo, sec.xhi, sec.ylo, sec.yhi));

        // Pad 2 — M_miss: Data vs Unsmeared vs Smeared
        c->cd(2);
        gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        drawComparisonPad(hdata_mmiss, hsim_mmiss, hunsmeared_mmiss,
                        "M_{miss} [GeV/c^{2}]", "mmiss");
        TLatex *t_mm = new TLatex(); t_mm->SetNDC(); t_mm->SetTextSize(0.046);
        t_mm->DrawLatex(0.15, 0.93, "Missing Mass M_{miss}");

        // Pad 3 — M_pgg2: Data vs Unsmeared vs Smeared
        c->cd(3);
        gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        drawComparisonPad(hdata_mpgg2, hsim_mpgg2, hunsmeared_mpgg2,
                        "(p_{target}+#gamma#gamma)^{2} [GeV^{2}]", "mpgg2");
        TLatex *t_mp = new TLatex(); t_mp->SetNDC(); t_mp->SetTextSize(0.046);
        t_mp->DrawLatex(0.15, 0.93, "(p_{target}+#gamma#gamma)^{2}");

        // -------------------------------------------------------------------------
        // Row 2: χ² scans
        // -------------------------------------------------------------------------

        // Pad 4 — χ² vs μ
        c->cd(4);
        gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        drawChi2Scan1D(fres.scan_data.mu_values, fres.scan_data.chi2_vs_mu,
                    fres.mu, fres.chi2, "#mu", "#chi^{2} vs #mu");

        // Pad 5 — χ² vs σ
        c->cd(5);
        gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        drawChi2Scan1D(fres.scan_data.sigma_values, fres.scan_data.chi2_vs_sigma,
                    fres.sigma, fres.chi2, "#sigma", "#chi^{2} vs #sigma");

        // Pad 6 — 2D χ² map (μ × σ)
        c->cd(6);
        gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.15); gPad->SetBottomMargin(0.12);
        if (!fres.scan_data.mu_2d.empty()) {
            TGraph2D *g2d = new TGraph2D(fres.scan_data.mu_2d.size(),
                                        const_cast<double*>(&fres.scan_data.mu_2d[0]),
                                        const_cast<double*>(&fres.scan_data.sigma_2d[0]),
                                        const_cast<double*>(&fres.scan_data.chi2_2d[0]));
            g2d->SetTitle("#chi^{2} landscape;#mu;#sigma;#chi^{2}");
            g2d->Draw("COLZ");
            TMarker *mk2d = new TMarker(fres.mu, fres.sigma, 29);
            mk2d->SetMarkerColor(kRed); mk2d->SetMarkerSize(2.0); mk2d->Draw();
        }

        // -------------------------------------------------------------------------
        // Row 3: Position scan, residuals, info
        // -------------------------------------------------------------------------

        // Pad 7 — χ² vs σ_pos (if enabled), else note
        c->cd(7);
        gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        if (Config::ENABLE_POSITION_SMEARING && !fres.scan_data.sigma_pos_values.empty()) {
            drawChi2Scan1D(fres.scan_data.sigma_pos_values, fres.scan_data.chi2_vs_sigma_pos,
                        fres.sigma_pos, fres.chi2, "#sigma_{pos} [cm]", "#chi^{2} vs #sigma_{pos}");
        } else {
            TLatex *tdis = new TLatex(); tdis->SetNDC(); tdis->SetTextSize(0.04);
            tdis->DrawLatex(0.15, 0.55, "Position smearing disabled");
            tdis->DrawLatex(0.15, 0.45, Form("#sigma_{pos} = 0 (fixed)"));
        }

        // Pad 8 — M_γγ pull: (Data − SmSim) / σ_bin
        c->cd(8);
        gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.12);
        {
            int nb = hdata.GetNbinsX();
            TH1D* hpull = new TH1D(Form("_diag_pull_%s", sec.name().c_str()),
                                    "(Data#minusSmeared)/#sigma;M_{#gamma#gamma} [GeV/c^{2}];Pull", nb,
                                    hdata.GetXaxis()->GetXmin(), hdata.GetXaxis()->GetXmax());
            hpull->GetXaxis()->SetTitleSize(0.05); hpull->GetYaxis()->SetTitleSize(0.05);
            hpull->GetXaxis()->SetLabelSize(0.04); hpull->GetYaxis()->SetLabelSize(0.04);
            for (int ib = 1; ib <= nb; ++ib) {
                double d = hdata.GetBinContent(ib);
                double s = hsim_final.GetBinContent(ib);
                double ed = hdata.GetBinError(ib);
                double es = hsim_final.GetBinError(ib);
                double denom = sqrt(ed*ed + es*es);
                if (denom > 0) hpull->SetBinContent(ib, (d - s) / denom);
            }
            hpull->SetLineColor(kBlack); hpull->SetFillColor(kCyan-9);
            double pmax = max(fabs(hpull->GetMinimum()), fabs(hpull->GetMaximum()));
            hpull->SetMaximum( max(pmax * 1.2, 3.0));
            hpull->SetMinimum(-max(pmax * 1.2, 3.0));
            hpull->Draw("HIST");
            TLine *lz = new TLine(hdata.GetXaxis()->GetXmin(), 0, hdata.GetXaxis()->GetXmax(), 0);
            lz->SetLineColor(kRed); lz->SetLineWidth(2); lz->Draw();
            TLine *lp2 = new TLine(hdata.GetXaxis()->GetXmin(),  2, hdata.GetXaxis()->GetXmax(),  2);
            TLine *lm2 = new TLine(hdata.GetXaxis()->GetXmin(), -2, hdata.GetXaxis()->GetXmax(), -2);
            lp2->SetLineColor(kOrange+1); lp2->SetLineStyle(2); lp2->Draw();
            lm2->SetLineColor(kOrange+1); lm2->SetLineStyle(2); lm2->Draw();
            TLatex *tpl = new TLatex(); tpl->SetNDC(); tpl->SetTextSize(0.040);
            tpl->DrawLatex(0.15, 0.87, "M_{#gamma#gamma} pull: (Data#minusSm.)/#sigma");
        }

        // Pad 9 — Info summary text
        c->cd(9);
        gPad->SetLeftMargin(0.05); gPad->SetRightMargin(0.05);
        {
            double chi2_ndf_val = fres.chi2_per_ndf();
            bool good = (chi2_ndf_val <= Config::MAX_CHI2_PER_NDF);
            TPaveText *info = new TPaveText(0.03, 0.03, 0.97, 0.97, "brNDC");
            info->SetBorderSize(1); info->SetFillColor(good ? kGreen-9 : kRed-9);
            info->SetTextAlign(12); info->SetTextSize(0.042);
            info->AddText(Form("Section  %s", sec.name().c_str()));
            info->AddText(Form("x: [%.1f, %.1f] cm", sec.xlo, sec.xhi));
            info->AddText(Form("y: [%.1f, %.1f] cm", sec.ylo, sec.yhi));
            info->AddText("─────────────────");
            info->AddText(Form("#mu        = %.5f", fres.mu));
            if (Config::ENABLE_ENERGY_DEPENDENT_MU)
                info->AddText(Form("#mu(E)=#mu_{0}((E_{0}/(E+E_{shift}))^{b}+1), E_{0}=%.3f GeV, E_{shift}=%.3f GeV, b=%.1f",
                                Config::MU_ENERGY_E0_GEV,
                                Config::MU_ENERGY_SHIFT_GEV,
                                Config::MU_ENERGY_B_EXPONENT));
            info->AddText(Form("#sigma     = %.5f", fres.sigma));
            if (Config::ENABLE_POSITION_SMEARING)
                info->AddText(Form("#sigma_{pos} = %.4f cm", fres.sigma_pos));
            if (Config::ENABLE_ELECTRON_MOMENTUM_SCALING && Config::stage2_uses_mmiss())
                info->AddText(Form("p_{e} scale = %.6f", fres.p_e_scale));
            info->AddText("─────────────────");
            info->AddText(Form("#chi^{2}     = %.2f", fres.chi2));
            info->AddText(Form("ndf       = %d", fres.ndf));
            info->AddText(Form("#chi^{2}/ndf = %.3f  %s", chi2_ndf_val, good ? "(OK)" : "[POOR]"));
            info->AddText("─────────────────");
            info->AddText(Form("N_{data}  = %d (all)", n_data));
            info->AddText(Form("N_{sel}   = %d (selected)", n_data_selected));
            info->AddText(Form("N_{sim}   = %d", n_sim));
            info->AddText(Form("N_{smear} = %d", nsmear));
            info->AddText(Form("Opt. mode: %s", Config::USE_THREE_STAGE_OPTIMIZATION ? "3-stage" : "simult."));
            info->AddText(Form("Obs: %s", Config::histogram_mode_label()));
            info->AddText(Form("Weights: w_{M_{#gamma#gamma}}=%.2f, w_{M_{miss}}=%.2f, w_{(p+#gamma#gamma)^{2}}=%.2f",
                            Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2));
            if (use_global_p_e_scale) {
                info->AddText(Form("Global p_{e} scale = %.6f", global_p_e_scale));
            }
            info->Draw();
        }
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
        cout << "Output file: " << out_file << endl;

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
        Int_t d_is_exclusive = 1;
        Double_t d_mmiss_all = 0;  // Pre-calculated missing mass from data tree
        
        tdata->SetBranchAddress("cluster_x_1", &d_cluster_x_1);
        tdata->SetBranchAddress("cluster_y_1", &d_cluster_y_1);
        tdata->SetBranchAddress("cluster_e_1", &d_cluster_e_1);
        tdata->SetBranchAddress("cluster_x_2", &d_cluster_x_2);
        tdata->SetBranchAddress("cluster_y_2", &d_cluster_y_2);
        tdata->SetBranchAddress("cluster_e_2", &d_cluster_e_2);
        tdata->SetBranchAddress("pi0_weight", &d_pi0_weight);
        tdata->SetBranchAddress("scale", &d_scale);
        bool has_data_is_exclusive = false;
        if (tdata->GetBranch("is_exclusive")) {
            tdata->SetBranchAddress("is_exclusive", &d_is_exclusive);
            has_data_is_exclusive = true;
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
        Float_t s_full_weight;
        Int_t s_is_exclusive = 1;
        
        // SIMC electron kinematics (for missing mass calculation)
        Float_t sc_e_E = 0, sc_e_Px = 0, sc_e_Py = 0, sc_e_Pz = 0;  // SIMC frame (MeV)
        Float_t s_e_E = 0, s_e_Px = 0, s_e_Py = 0, s_e_Pz = 0;  // Transformed (GeV)
        
        tsim->SetBranchAddress("clust_X", s_clust_X);
        tsim->SetBranchAddress("clust_Y", s_clust_Y);
        tsim->SetBranchAddress("clust_E", s_clust_E);
        tsim->SetBranchAddress("nclust", &s_nclust);
        tsim->SetBranchAddress("full_weight", &s_full_weight);
        // Optionally ignore sim full_weight if configured

        bool has_sim_is_exclusive = false;
        if (tsim->GetBranch("is_exclusive")) {
            tsim->SetBranchAddress("is_exclusive", &s_is_exclusive);
            has_sim_is_exclusive = true;
        }

        if (Config::APPLY_IS_EXCLUSIVE_SELECTION) {
            cout << "Exclusivity weighting mode: enabled" << endl;
            if (!has_data_is_exclusive) {
                cerr << "WARNING: Data branch 'is_exclusive' not found; using factor=1 for data." << endl;
            }
            if (!has_sim_is_exclusive) {
                cerr << "WARNING: Simulation branch 'is_exclusive' not found; using factor=1 for simulation." << endl;
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
            if (Config::APPLY_IS_EXCLUSIVE_SELECTION && has_data_is_exclusive) {
                data_exclusive_factor = static_cast<double>(d_is_exclusive);
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
            if (Config::APPLY_IS_EXCLUSIVE_SELECTION && has_sim_is_exclusive) {
                sim_exclusive_factor = static_cast<double>(s_is_exclusive);
            }
            pair.weight = (Config::USE_SIM_FULL_WEIGHT ? s_full_weight : 1.0) * sim_exclusive_factor;
            
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

        // RNG - will be created per-thread in parallel region
        TRandom3 rng(0);  // Main thread RNG

        // Determine global electron momentum scale once (shared across sections)
        bool use_global_p_e_scale = Config::ENABLE_ELECTRON_MOMENTUM_SCALING &&
                                    Config::stage2_uses_mmiss();
        double global_p_e_scale = Config::GLOBAL_PE_SCALE_DEFAULT;

        if (use_global_p_e_scale) {
            cout << "\n==== Global electron momentum scale mode enabled ====" << endl;
            if (Config::ENABLE_FINAL_GLOBAL_PE_SCALE_STAGE4) {
                cout << "Pass-2 will jointly refine per-section p_e_scale with (mu, sigma, sigma_pos), then Stage 4 will do final global p_e_scale refinement." << endl;
            } else {
                cout << "Pass-2 will jointly refine per-section p_e_scale with (mu, sigma, sigma_pos); final global Stage-4 p_e_scale refinement is disabled." << endl;
            }
        }

        // Output file
        TFile fout(out_file.c_str(), "RECREATE");

        // CSV summary - place in same directory as ROOT output
        string csv_file = out_file;
        size_t last_slash = csv_file.find_last_of('/');
        if (last_slash != string::npos) {
            csv_file = csv_file.substr(0, last_slash + 1) + Config::CSV_FILENAME;
        } else {
            csv_file = Config::CSV_FILENAME;
        }
        ofstream csv(csv_file.c_str());
        csv << "ix,iy,xlo,xhi,ylo,yhi,n_data,n_sim,best_mu,best_sigma,best_sigma_pos_cm,best_p_e_scale,res_A,res_B,res_C,best_chi2,ndf,chi2_ndf,fit_status\n";

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
                    gev.mu1_ext = 1.0; gev.sigma1_ext = 0.0; gev.sigma_pos1_ext = 0.0;
                    gev.mu2_ext = 1.0; gev.sigma2_ext = 0.0; gev.sigma_pos2_ext = 0.0;

                    int is1 = in_geometry(ev.x1, ev.y1) ? best_section_for_point(ev.x1, ev.y1) : -1;
                    if (is1 >= 0) {
                        gev.mu1_ext = src_results[is1].mu;
                        gev.sigma1_ext = src_results[is1].sigma;
                        gev.sigma_pos1_ext = src_results[is1].sigma_pos;
                    }
                    int is2 = in_geometry(ev.x2, ev.y2) ? best_section_for_point(ev.x2, ev.y2) : -1;
                    if (is2 >= 0) {
                        gev.mu2_ext = src_results[is2].mu;
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
                    double chi2 = eval_chi2_mmiss_only(1.0,
                                                    0.0,
                                                    0.0,
                                                    pe,
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
                        double chi2 = eval_chi2_mmiss_only(1.0,
                                                        0.0,
                                                        0.0,
                                                        pe,
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
                    auto chi2_eval = [&](double /*mu*/, double /*sigma*/, double /*sigma_pos*/, double p_e_scale) {
                        TRandom3 rng_eval(456789);
                        return eval_chi2_mmiss_only(1.0,
                                                    0.0,
                                                    0.0,
                                                    p_e_scale,
                                                    sim_events_global_stage4,
                                                    hdata_global_mmiss,
                                                    rng_eval,
                                                    Nsmear,
                                                    Config::RESOLUTION_A_DEFAULT,
                                                    Config::RESOLUTION_B_DEFAULT,
                                                    Config::RESOLUTION_C_DEFAULT);
                    };
                    MigradRefineResult mr = run_migrad_refinement(chi2_eval,
                                                                1.0, 0.0, 0.0, best_scale,
                                                                false, false, false, true);
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
        string pdf_file = out_file;
        size_t pdf_slash = pdf_file.find_last_of('/');
        if (pdf_slash != string::npos) {
            pdf_file = pdf_file.substr(0, pdf_slash + 1) + Config::CHI2_PDF_FILENAME;
        } else {
            pdf_file = Config::CHI2_PDF_FILENAME;
        }
        TCanvas *c_chi2 = new TCanvas("c_chi2", "Chi2 Scans", 1400, 1000);
        string pdf_open = pdf_file + "[";
        c_chi2->Print(pdf_open.c_str());
        int page_count = 0;

        // Loop over sections and fit each where we have enough stats
        // OpenMP parallelization: each section is fitted independently
        
        // ROOT is not thread-safe by default. Disable automatic histogram registration
        Bool_t original_adddir_status = TH1::AddDirectoryStatus();
        TH1::AddDirectory(kFALSE);
        
        // Report number of threads
        int num_threads = 1;
        #pragma omp parallel
        {
            #pragma omp single
            {
                num_threads = omp_get_num_threads();
                cout << "\n==== Starting parallel fitting with " << num_threads << " threads ====" << endl;
            }
        }

        // =========================================================================
        // === PASS 1: Fit each section using only "both-photons-in-section" pairs
        // =========================================================================
        // Build per-section subsets containing only pairs where BOTH photons hit
        // the same (current) section.  These pairs are the cleanest calibration
        // sample — no cross-calibration ambiguity.
        vector<vector<ClusterPair>> sim_events_both_in_sec(nsec);
        vector<int> sim_both_selected_count(nsec, 0);
        for (int is = 0; is < nsec; ++is) {
            for (const auto &ev : sim_events_per_section[is]) {
                if (ev.photon1_in_section && ev.photon2_in_section) {
                    sim_events_both_in_sec[is].push_back(ev);
                    if (ev.weight > 0.0) ++sim_both_selected_count[is];
                }
            }
            cout << "Section " << sections[is].name()
                << "  both-in-sec sim=" << sim_events_both_in_sec[is].size()
                << " (selected=" << sim_both_selected_count[is] << ")"
                << "  all-sec sim="     << sim_events_per_section[is].size() << "\n";
        }

        vector<FitResult> fit_results_pass1(nsec);
        vector<bool>      fit_success_pass1(nsec, false);
        vector<int>       pass1_seed_mode(nsec, 0); // 0=skip, 1=both-in-sec, 2=all-sec fallback

        cout << "\n==== PASS 1: fitting sections with both-photon-in-section pairs ====\n";
        #pragma omp parallel
        {
            int thread_id_p1 = omp_get_thread_num();
            TRandom3 rng_p1(thread_id_p1 * 234567 + time(NULL));

            #pragma omp for schedule(dynamic)
            for (int is = 0; is < nsec; ++is) {
                size_t nsim_both = (size_t)sim_both_selected_count[is];
                size_t nsim_all = (size_t)sim_selected_count_per_section[is];
                size_t ndata_sec = (size_t)data_selected_count_per_section[is];

                const vector<ClusterPair>* pass1_events = &sim_events_both_in_sec[is];
                size_t nsim_seed = nsim_both;
                int seed_mode = 1;
                if (nsim_seed < (size_t)Config::MIN_EVENTS_PER_SECTION &&
                    nsim_all >= (size_t)Config::MIN_EVENTS_PER_SECTION) {
                    pass1_events = &sim_events_per_section[is];
                    nsim_seed = nsim_all;
                    seed_mode = 2;
                }

                if (ndata_sec < Config::MIN_EVENTS_PER_SECTION || nsim_seed < (size_t)Config::MIN_EVENTS_PER_SECTION) {
                    #pragma omp critical(console)
                    {
                        cout << "[Pass1] Skipping section " << sections[is].name()
                            << "  data_sel=" << ndata_sec
                            << "  both-in-sec sim=" << nsim_both
                            << "  all-sec sim=" << nsim_all << "\n";
                    }
                    continue;
                }

                FitResult fres1 = fit_section(*pass1_events,
                                            *hdata_sec[is], *hdata_mmiss_sec[is], *hdata_mpgg2_sec[is],
                                            rng_p1, Nsmear, global_p_e_scale, false);
                #pragma omp critical(console)
                {
                    cout << "[Pass1] Section " << sections[is].name()
                        << "  seed=" << (seed_mode == 1 ? "both-in-sec" : "all-sec-fallback")
                        << "  mu=" << fres1.mu << "  sigma=" << fres1.sigma
                        << "  chi2/ndf=" << fres1.chi2_per_ndf() << "\n";
                }
                fit_results_pass1[is] = fres1;
                fit_success_pass1[is] = true;
                pass1_seed_mode[is] = seed_mode;
            }
        }

        {
            int n_pass1_both = 0;
            int n_pass1_fallback = 0;
            int n_pass1_skipped = 0;
            for (int is = 0; is < nsec; ++is) {
                if (!fit_success_pass1[is]) { ++n_pass1_skipped; continue; }
                if (pass1_seed_mode[is] == 1) ++n_pass1_both;
                else if (pass1_seed_mode[is] == 2) ++n_pass1_fallback;
            }
            cout << "Pass-1 seed summary: both-in-sec=" << n_pass1_both
                << ", all-sec-fallback=" << n_pass1_fallback
                << ", skipped=" << n_pass1_skipped << "\n";
        }

        // =========================================================================
        // === Cross-boundary ext parameter assignment
        // =========================================================================
        // Shared helper used both for the initial Pass-1→Pass-2 handoff and for
        // all subsequent Pass-2 coupled sweeps, ensuring a single consistent
        // algorithm: nearest-center match → base-grid fallback → interpolation.
        auto refresh_cross_boundary_ext = [&](const vector<FitResult> &src_results,
                                            const vector<bool> &src_success) {
            CalibrationMap coupled_map(nx, ny, x_min, x_max, y_min, y_max);
            for (int js = 0; js < nsec; ++js) {
                if (!src_success[js]) continue;
                coupled_map.setParams(sections[js].ix, sections[js].iy,
                                    src_results[js].mu,
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
            for (int is = 0; is < nsec; ++is) {
                for (auto &ev : sim_events_per_section[is]) {
                    if (!ev.photon1_in_section) {
                        int js = best_fitted_overlap_section(ev.x1, ev.y1);
                        if (js < 0) {
                            int jb = base_grid_section_index(ev.x1, ev.y1);
                            if (jb >= 0 && src_success[jb]) js = jb;
                        }
                        if (js >= 0) {
                            ev.mu1_ext        = src_results[js].mu;
                            ev.sigma1_ext     = src_results[js].sigma;
                            ev.sigma_pos1_ext = src_results[js].sigma_pos;
                        } else if (in_geometry(ev.x1, ev.y1)) {
                            coupled_map.getInterpolatedParams(ev.x1, ev.y1,
                                                            ev.mu1_ext, ev.sigma1_ext, ev.sigma_pos1_ext);
                            ++interpolated;
                        }
                        ++assigned;
                    }
                    if (!ev.photon2_in_section) {
                        int js = best_fitted_overlap_section(ev.x2, ev.y2);
                        if (js < 0) {
                            int jb = base_grid_section_index(ev.x2, ev.y2);
                            if (jb >= 0 && src_success[jb]) js = jb;
                        }
                        if (js >= 0) {
                            ev.mu2_ext        = src_results[js].mu;
                            ev.sigma2_ext     = src_results[js].sigma;
                            ev.sigma_pos2_ext = src_results[js].sigma_pos;
                        } else if (in_geometry(ev.x2, ev.y2)) {
                            coupled_map.getInterpolatedParams(ev.x2, ev.y2,
                                                            ev.mu2_ext, ev.sigma2_ext, ev.sigma_pos2_ext);
                            ++interpolated;
                        }
                        ++assigned;
                    }
                }
            }
            return std::make_pair(assigned, interpolated);
        };

        // =========================================================================
        // === Between passes: assign ext calibration to cross-boundary photons
        // =========================================================================
        // Use the same nearest-center + base-grid + interpolation algorithm that
        // Pass-2 sweeps use, so the initial ext assignment is fully consistent.
        cout << "\n==== Between passes: assigning ext calibration for cross-boundary photons ====\n";
        {
            auto init_counts = refresh_cross_boundary_ext(fit_results_pass1, fit_success_pass1);
            cout << "Initial ext assignment: " << init_counts.first
                << " cross-boundary photon coefficient assignments";
            if (init_counts.second > 0) cout << " (interpolation fallback=" << init_counts.second << ")";
            cout << "." << endl;
        }

        // =========================================================================
        // === PASS 2: Iterative coupled fit using all section events
        // =========================================================================
        auto fit_section_fast_refine = [&](const vector<ClusterPair> &simEvents,
                                        const TH1D &hdata_mpi0,
                                        const TH1D &hdata_mmiss,
                                        const TH1D &hdata_mpgg2,
                                        TRandom3 &rng_local,
                                        int nsmear_local,
                                        double p_e_scale_local,
                                        const FitResult &seed) {
            FitResult out = seed;
            out.scan_data = Chi2ScanData();
            const bool fit_p_e_local = Config::ENABLE_ELECTRON_MOMENTUM_SCALING && Config::stage2_uses_mmiss();
            double pe_seed = fit_p_e_local ? seed.p_e_scale : p_e_scale_local;
            if (!std::isfinite(pe_seed)) pe_seed = p_e_scale_local;
            pe_seed = min(max(pe_seed, Config::PE_SCALE_MIN), Config::PE_SCALE_MAX);
            out.p_e_scale = pe_seed;
            out.res_A = Config::RESOLUTION_A_DEFAULT;
            out.res_B = Config::RESOLUTION_B_DEFAULT;
            out.res_C = Config::RESOLUTION_C_DEFAULT;

            int n_params = 2 + (Config::ENABLE_POSITION_SMEARING ? 1 : 0) + (fit_p_e_local ? 1 : 0);
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

            double mu0 = min(max(seed.mu, Config::MU_MIN), Config::MU_MAX);
            double sigma0 = min(max(seed.sigma, max(1e-6, Config::SIGMA_MIN)), Config::SIGMA_MAX);
            double sigma_pos0 = Config::ENABLE_POSITION_SMEARING
                                    ? min(max(seed.sigma_pos, Config::SIGMA_POS_MIN), Config::SIGMA_POS_MAX)
                                    : 0.0;
            double p_e0 = pe_seed;

            double best_chi2 = eval_chi2_selected(mu0, sigma0, sigma_pos0, p_e0,
                                                simEvents,
                                                hdata_mpi0, hdata_mmiss, hdata_mpgg2,
                                                rng_local, nsmear_local,
                                                Config::RESOLUTION_A_DEFAULT,
                                                Config::RESOLUTION_B_DEFAULT,
                                                Config::RESOLUTION_C_DEFAULT,
                                                Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);

            double step_mu = 0.0025;
            double step_sig = 0.005;
            double step_pos = Config::ENABLE_POSITION_SMEARING ? 0.005 : 0.0;
            double step_pe = fit_p_e_local ? (Config::PE_SCALE_MAX - Config::PE_SCALE_MIN) / 20.0 : 0.0;
            int iref = 0;
            while ((step_mu >= 0.0001 || step_sig >= 0.0001 || (Config::ENABLE_POSITION_SMEARING && step_pos >= 0.0001) || (fit_p_e_local && step_pe >= 1e-4)) &&
                iref < Config::MAX_REFINEMENT_ITERATIONS) {
                double bmu = mu0, bsig = sigma0, bpos = sigma_pos0, bpe = p_e0, bchi = best_chi2;
                for (double mu = max(Config::MU_MIN, mu0 - 2*step_mu);
                    mu <= min(Config::MU_MAX, mu0 + 2*step_mu) + 1e-15; mu += step_mu) {
                    for (double sig = max(max(1e-6, Config::SIGMA_MIN), sigma0 - 2*step_sig);
                        sig <= min(Config::SIGMA_MAX, sigma0 + 2*step_sig) + 1e-15; sig += step_sig) {
                        double ps0 = Config::ENABLE_POSITION_SMEARING ? max(Config::SIGMA_POS_MIN, sigma_pos0 - 2*step_pos) : sigma_pos0;
                        double ps1 = Config::ENABLE_POSITION_SMEARING ? min(Config::SIGMA_POS_MAX, sigma_pos0 + 2*step_pos) : sigma_pos0;
                        double pss = Config::ENABLE_POSITION_SMEARING ? step_pos : 1.0;
                        for (double spos = ps0; spos <= ps1 + 1e-15; spos += pss) {
                            double pe0 = fit_p_e_local ? max(Config::PE_SCALE_MIN, p_e0 - 2 * step_pe) : p_e0;
                            double pe1 = fit_p_e_local ? min(Config::PE_SCALE_MAX, p_e0 + 2 * step_pe) : p_e0;
                            double pes = fit_p_e_local ? step_pe : 1.0;
                            for (double pe = pe0; pe <= pe1 + 1e-15; pe += pes) {
                                double chi2 = eval_chi2_selected(mu, sig, spos, pe,
                                                                simEvents,
                                                                hdata_mpi0, hdata_mmiss, hdata_mpgg2,
                                                                rng_local, nsmear_local,
                                                                Config::RESOLUTION_A_DEFAULT,
                                                                Config::RESOLUTION_B_DEFAULT,
                                                                Config::RESOLUTION_C_DEFAULT,
                                                                Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);
                                if (chi2 < bchi) { bchi = chi2; bmu = mu; bsig = sig; bpos = spos; bpe = pe; }
                            }
                        }
                    }
                }
                mu0 = bmu; sigma0 = bsig; sigma_pos0 = bpos; p_e0 = bpe; best_chi2 = bchi;
                step_mu /= 2.0;
                step_sig /= 2.0;
                if (Config::ENABLE_POSITION_SMEARING) step_pos /= 2.0;
                if (fit_p_e_local) step_pe /= 2.0;
                ++iref;
            }

            out.mu = mu0;
            out.sigma = sigma0;
            out.sigma_pos = sigma_pos0;
            out.p_e_scale = p_e0;
            out.chi2 = best_chi2;

            {
                auto chi2_eval = [&](double mu, double sigma, double sigma_pos, double p_e_scale) {
                    TRandom3 rng_eval(97531);
                    return eval_chi2_selected(mu, sigma, sigma_pos, p_e_scale,
                                            simEvents,
                                            hdata_mpi0, hdata_mmiss, hdata_mpgg2,
                                            rng_eval, nsmear_local,
                                            Config::RESOLUTION_A_DEFAULT,
                                            Config::RESOLUTION_B_DEFAULT,
                                            Config::RESOLUTION_C_DEFAULT,
                                            Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);
                };
                MigradRefineResult mr = run_migrad_refinement(chi2_eval,
                                                            mu0, sigma0, sigma_pos0, p_e0,
                                                            true, true,
                                                            Config::ENABLE_POSITION_SMEARING,
                                                            fit_p_e_local);
                if (mr.minimized) {
                    out.mu = mr.mu;
                    out.sigma = mr.sigma;
                    out.sigma_pos = mr.sigma_pos;
                    out.p_e_scale = mr.p_e_scale;
                    out.chi2 = mr.chi2;
                }
            }
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

        // Initialize Pass-2 seeds from Pass-1 values:
        //   - direct: section has a successful Pass-1 fit
        //   - derived: nearest successful Pass-1 section (center distance)
        //   - fallback: deterministic defaults only if no successful Pass-1 seeds exist
        int n_seeded_from_pass1 = 0;
        int n_seeded_from_nearest_pass1 = 0;
        int n_seeded_default = 0;
        auto nearest_successful_pass1_section = [&](int is_ref) {
            int best_js = -1;
            double best_d2 = std::numeric_limits<double>::max();
            double x_ref = x_min + (sections[is_ref].ix + 0.5) * base_wx;
            double y_ref = y_min + (sections[is_ref].iy + 0.5) * base_wy;
            for (int js = 0; js < nsec; ++js) {
                if (!fit_success_pass1[js]) continue;
                double x_js = x_min + (sections[js].ix + 0.5) * base_wx;
                double y_js = y_min + (sections[js].iy + 0.5) * base_wy;
                double dx = x_ref - x_js;
                double dy = y_ref - y_js;
                double d2 = dx * dx + dy * dy;
                if (d2 < best_d2) {
                    best_d2 = d2;
                    best_js = js;
                }
            }
            return best_js;
        };

        for (int is = 0; is < nsec; ++is) {
            if (!section_active[is]) continue;

            if (fit_success_pass1[is]) {
                fit_results[is] = fit_results_pass1[is];
                fit_results[is].p_e_scale = global_p_e_scale;
                fit_results[is].scan_data = Chi2ScanData();
                fit_success[is] = true;
                ++n_seeded_from_pass1;
                continue;
            }

            int js = nearest_successful_pass1_section(is);
            if (js >= 0) {
                fit_results[is] = fit_results_pass1[js];
                fit_results[is].p_e_scale = global_p_e_scale;
                fit_results[is].scan_data = Chi2ScanData();
                fit_success[is] = false;
                ++n_seeded_from_nearest_pass1;
            } else {
                FitResult seed{};
                seed.mu = 1.0;
                seed.sigma = 0.05;
                seed.sigma_pos = 0.0;
                seed.p_e_scale = global_p_e_scale;
                seed.chi2 = 1e300;
                seed.res_A = Config::RESOLUTION_A_DEFAULT;
                seed.res_B = Config::RESOLUTION_B_DEFAULT;
                seed.res_C = Config::RESOLUTION_C_DEFAULT;
                seed.ndf = 1;
                seed.scan_data = Chi2ScanData();
                fit_results[is] = seed;
                fit_success[is] = false;
                ++n_seeded_default;
            }
        }
        cout << "Pass-2 sweep-1 seeding: from Pass-1=" << n_seeded_from_pass1
            << ", from nearest Pass-1=" << n_seeded_from_nearest_pass1
            << ", default-fallback=" << n_seeded_default << "\n";

        // Start coupled iteration from pass-1 ext assignment already performed above.
        const int max_sweeps = Config::ENABLE_ITERATIVE_COUPLED_PASS2 ? max(1, Config::COUPLED_PASS2_MAX_SWEEPS) : 1;
        vector<FitResult> prev_results(nsec);
        vector<bool> prev_success(nsec, false);

        cout << "\n==== PASS 2: coupled section fitting (max sweeps=" << max_sweeps << ") ====\n";
        for (int sweep = 0; sweep < max_sweeps; ++sweep) {
            cout << "\n-- Pass-2 sweep " << (sweep + 1) << " / " << max_sweeps << " --" << endl;

            vector<FitResult> sweep_results = fit_results;
            vector<bool> sweep_success = fit_success;

            #pragma omp parallel
            {
                int thread_id = omp_get_thread_num();
                TRandom3 thread_rng((sweep + 1) * 1000003 + thread_id * 123456 + time(NULL));

                #pragma omp for schedule(dynamic)
                for (int is = 0; is < nsec; ++is) {
                    if (!section_active[is]) continue;

                    FitResult fres = fit_section_fast_refine(sim_events_per_section[is],
                                                            *hdata_sec[is], *hdata_mmiss_sec[is], *hdata_mpgg2_sec[is],
                                                            thread_rng, Nsmear,
                                                            fit_results[is].p_e_scale,
                                                            fit_results[is]);

                    double chi2_ndf = fres.chi2_per_ndf();
                    #pragma omp critical(console)
                    {
                        cout << "[sweep " << (sweep + 1) << "] " << sections[is].name()
                            << " -> mu=" << fres.mu
                            << " sigma=" << fres.sigma
                            << " sigma_pos=" << std::fixed << std::setprecision(5) << fres.sigma_pos << " cm";
                        if (Config::ENABLE_ELECTRON_MOMENTUM_SCALING && Config::stage2_uses_mmiss()) {
                            cout << " p_e_scale=" << fres.p_e_scale;
                        }
                        cout << " chi2/ndf=" << chi2_ndf;
                        if (chi2_ndf > Config::MAX_CHI2_PER_NDF) {
                            cout << " [POOR]";
                        }
                        cout << "\n";
                    }

                    sweep_results[is] = fres;
                    sweep_success[is] = true;
                }
            }

            double max_dmu = 0.0, max_dsig = 0.0, max_dpos = 0.0;
            int ncomp = 0;
            for (int is = 0; is < nsec; ++is) {
                if (!(section_active[is] && sweep_success[is] && prev_success[is])) continue;
                max_dmu = max(max_dmu, fabs(sweep_results[is].mu - prev_results[is].mu));
                max_dsig = max(max_dsig, fabs(sweep_results[is].sigma - prev_results[is].sigma));
                max_dpos = max(max_dpos, fabs(sweep_results[is].sigma_pos - prev_results[is].sigma_pos));
                ++ncomp;
            }

            fit_results = sweep_results;
            fit_success = sweep_success;
            auto refresh_counts = refresh_cross_boundary_ext(fit_results, fit_success);
            int n_assigned = refresh_counts.first;
            int n_interp = refresh_counts.second;
            cout << "Sweep " << (sweep + 1) << ": refreshed " << n_assigned
                << " cross-boundary photon coefficient assignments";
            if (n_interp > 0) cout << " (interpolation fallback=" << n_interp << ")";
            cout << "." << endl;

            bool converged = (ncomp > 0 &&
                            max_dmu <= Config::COUPLED_CONV_MU &&
                            max_dsig <= Config::COUPLED_CONV_SIGMA &&
                            max_dpos <= Config::COUPLED_CONV_SIGMA_POS);
            if (ncomp > 0) {
                cout << "Sweep " << (sweep + 1)
                    << " deltas: dmu=" << max_dmu
                    << " dsigma=" << max_dsig
                    << " dsigma_pos=" << max_dpos << endl;
            }

            prev_results = fit_results;
            prev_success = fit_success;

            if (Config::ENABLE_ITERATIVE_COUPLED_PASS2 && converged) {
                cout << "Pass-2 coupled fit converged after sweep " << (sweep + 1) << "." << endl;
                break;
            }
        }

        // Final fit status labeling and optional poor-fit skipping
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

        // Final coherence refresh: enforce cross-boundary ext coefficients to be
        // consistent with the final accepted section mask (after optional skipping).
        {
            auto final_refresh_counts = refresh_cross_boundary_ext(fit_results, fit_success);
            cout << "Final Pass-2 refresh: " << final_refresh_counts.first
                << " cross-boundary photon assignments";
            if (final_refresh_counts.second > 0) {
                cout << " (interpolation fallback=" << final_refresh_counts.second << ")";
            }
            cout << "." << endl;
        }

        // =========================================================================
        // === STAGE 4: Global p_e_scale using final Stage-3 section parameters
        // =========================================================================
        if (use_global_p_e_scale && Config::ENABLE_FINAL_GLOBAL_PE_SCALE_STAGE4) {
            cout << "\n==== STAGE 4 (GLOBAL): final p_e_scale refinement with final section smearing ====" << endl;
            double best_chi2_global = 0.0;
            bool pe_ok = optimize_global_p_e_scale_for_fit(fit_results, fit_success,
                                                            global_p_e_scale,
                                                            "Stage 4",
                                                            &best_chi2_global);
            if (pe_ok) {
                for (int is = 0; is < nsec; ++is) {
                    if (fit_success[is]) fit_results[is].p_e_scale = global_p_e_scale;
                }

                cout << "Stage 4 result: global p_e_scale=" << global_p_e_scale
                    << " (chi2_mmiss=" << best_chi2_global << ")" << endl;
            }
        } else if (use_global_p_e_scale) {
            cout << "\n==== STAGE 4 (GLOBAL): skipped by configuration (ENABLE_FINAL_GLOBAL_PE_SCALE_STAGE4=false) ====" << endl;
        }

        // Recompute final chi2 with finalized coefficients and finalized p_e_scale
        // so reported quality matches the exact model used for final output plots.
        cout << "\n==== Recomputing final section chi2 with finalized model state ====\n";
        for (int is = 0; is < nsec; ++is) {
            if (!fit_success[is]) continue;
            TRandom3 rng_final_chi2(8000 + is);
            double chi2_final = eval_chi2_selected(fit_results[is].mu,
                                                fit_results[is].sigma,
                                                fit_results[is].sigma_pos,
                                                fit_results[is].p_e_scale,
                                                sim_events_per_section[is],
                                                *hdata_sec[is], *hdata_mmiss_sec[is], *hdata_mpgg2_sec[is],
                                                rng_final_chi2, Nsmear,
                                                Config::RESOLUTION_A_DEFAULT,
                                                Config::RESOLUTION_B_DEFAULT,
                                                Config::RESOLUTION_C_DEFAULT,
                                                Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2);
            fit_results[is].chi2 = chi2_final;

            bool good_fit = (fit_results[is].chi2_per_ndf() <= Config::MAX_CHI2_PER_NDF);
            if (Config::SKIP_BAD_FITS && !good_fit) {
                fit_status[is] = "poor_fit_skipped";
                fit_success[is] = false;
            } else {
                fit_status[is] = good_fit ? "fit_ok" : "poor_fit";
            }
        }

        // Regenerate visualization scan data from FINAL model state only
        // (latest coupled sweep + final p_e_scale/status state).
        cout << "\n==== Regenerating visualization scans from final sweep values ====\n";
        #pragma omp parallel
        {
            int thread_id_scan = omp_get_thread_num();
            TRandom3 rng_final_scan(8500 + thread_id_scan * 104729 + time(NULL));

            #pragma omp for schedule(dynamic)
            for (int is = 0; is < nsec; ++is) {
                if (!fit_success[is]) continue;
                generate_visualization_scan_data(sim_events_per_section[is],
                                                *hdata_sec[is], *hdata_mmiss_sec[is], *hdata_mpgg2_sec[is],
                                                rng_final_scan, Nsmear,
                                                fit_results[is]);
            }
        }

        // =========================================================================
        // === Final per-section outputs (using final Stage-3 + Stage-4 parameters)
        // =========================================================================
        cout << "\n==== Writing final per-section histograms and PDF pages ====" << endl;
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
                    << sections[is].ylo << "," << sections[is].yhi << "," << ndata_sec << "," << nsim_sec << ",,,,,,,,,,"
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
                                        fres.mu, fres.sigma, fres.sigma_pos,
                                        fres.p_e_scale,
                                        rng_out, Nsmear,
                                        fres.res_A, fres.res_B, fres.res_C);

            if (hunsmeared.Integral() > 0 && hdata_sec[is]->Integral() > 0)
                hunsmeared.Scale(hdata_sec[is]->Integral() / hunsmeared.Integral());
            if (hfinal.Integral() > 0 && hdata_sec[is]->Integral() > 0)
                hfinal.Scale(hdata_sec[is]->Integral() / hfinal.Integral());
            if (hfinal_mmiss.Integral() > 0 && hdata_mmiss_sec[is]->Integral() > 0)
                hfinal_mmiss.Scale(hdata_mmiss_sec[is]->Integral() / hfinal_mmiss.Integral());
            if (hfinal_mpgg2.Integral() > 0 && hdata_mpgg2_sec[is]->Integral() > 0)
                hfinal_mpgg2.Scale(hdata_mpgg2_sec[is]->Integral() / hfinal_mpgg2.Integral());

            const double z_nps_diag = nps::kDefaultZ_NPS_cm;
            string sec_tag = sections[is].name();
            TH1D hunsmeared_mmiss(("h_unsmeared_mmiss_" + sec_tag).c_str(), "h_unsmeared_mmiss",
                                Config::MMISS_NBINS, Config::MMISS_MIN, Config::MMISS_MAX);
            TH1D hunsmeared_mpgg2(("h_unsmeared_mpgg2_" + sec_tag).c_str(), "h_unsmeared_mpgg2",
                                Config::MPGG2_NBINS, Config::MPGG2_MIN, Config::MPGG2_MAX);
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
            if (hunsmeared_mmiss.Integral() > 0 && hdata_mmiss_sec[is]->Integral() > 0)
                hunsmeared_mmiss.Scale(hdata_mmiss_sec[is]->Integral() / hunsmeared_mmiss.Integral());
            if (hunsmeared_mpgg2.Integral() > 0 && hdata_mpgg2_sec[is]->Integral() > 0)
                hunsmeared_mpgg2.Scale(hdata_mpgg2_sec[is]->Integral() / hunsmeared_mpgg2.Integral());

            size_t ndata_sec_diag = data_mass_per_section[is].size();
            int nsel_sec_diag = data_selected_count_per_section[is];
            int nsim_sec_diag = sim_selected_count_per_section[is];
            plot_chi2_scans(sections[is], fres, c_chi2,
                            *hdata_sec[is], hfinal, hunsmeared,
                            *hdata_mmiss_sec[is], hfinal_mmiss, hunsmeared_mmiss,
                            *hdata_mpgg2_sec[is], hfinal_mpgg2, hunsmeared_mpgg2,
                    (int)ndata_sec_diag, nsel_sec_diag, nsim_sec_diag,
                            Nsmear, global_p_e_scale, use_global_p_e_scale);
            ++page_count;
            c_chi2->Print(pdf_file.c_str());

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
            pt->AddText(Form("Best mu = %.5f", fres.mu));
            if (Config::ENABLE_ENERGY_DEPENDENT_MU) {
                pt->AddText(Form("mu(E)=mu0*((E0/(E+Eshift))^b+1), E0=%.3f GeV, Eshift=%.3f GeV, b=%.1f",
                                Config::MU_ENERGY_E0_GEV,
                                Config::MU_ENERGY_SHIFT_GEV,
                                Config::MU_ENERGY_B_EXPONENT));
            }
            pt->AddText(Form("Best sigma = %.5f", fres.sigma));
            pt->AddText(Form("Best sigma_pos = %.5f cm", fres.sigma_pos));
            if (Config::ENABLE_ELECTRON_MOMENTUM_SCALING && Config::stage2_uses_mmiss()) {
                pt->AddText(Form("Best p_e_scale = %.6f", fres.p_e_scale));
            }
            pt->AddText(Form("Chi2 = %.2f, ndf = %d, chi2/ndf = %.2f", fres.chi2, fres.ndf, fres.chi2_per_ndf()));
            pt->Write("fit_params");
            delete pt;

            fout.cd();
            csv << sections[is].ix << "," << sections[is].iy << ","
                << sections[is].xlo << "," << sections[is].xhi << "," << sections[is].ylo << "," << sections[is].yhi << ","
                << ndata_sec << "," << nsim_sec << ","
                << fres.mu << "," << fres.sigma << "," << fres.sigma_pos << "," << fres.p_e_scale << ","
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
                                        fit_results[is].mu,
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
            gev.mu1_ext = 1.0; gev.sigma1_ext = 0.0; gev.sigma_pos1_ext = 0.0;
            gev.mu2_ext = 1.0; gev.sigma2_ext = 0.0; gev.sigma_pos2_ext = 0.0;

            int is1 = section_index_for_photon_coeff(ev.x1, ev.y1);
            if (is1 >= 0) {
                gev.mu1_ext = fit_results[is1].mu;
                gev.sigma1_ext = fit_results[is1].sigma;
                gev.sigma_pos1_ext = fit_results[is1].sigma_pos;
            } else if (in_geometry(ev.x1, ev.y1)) {
                all_sections_cal_map.getInterpolatedParams(ev.x1, ev.y1,
                                                        gev.mu1_ext, gev.sigma1_ext, gev.sigma_pos1_ext);
                ++n_summary_missing_section_coeff;
            }
            int is2 = section_index_for_photon_coeff(ev.x2, ev.y2);
            if (is2 >= 0) {
                gev.mu2_ext = fit_results[is2].mu;
                gev.sigma2_ext = fit_results[is2].sigma;
                gev.sigma_pos2_ext = fit_results[is2].sigma_pos;
            } else if (in_geometry(ev.x2, ev.y2)) {
                all_sections_cal_map.getInterpolatedParams(ev.x2, ev.y2,
                                                        gev.mu2_ext, gev.sigma2_ext, gev.sigma_pos2_ext);
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
                                    1.0, 0.0, 0.0,
                                    global_p_e_scale,
                                    rng_combined,
                                    Nsmear,
                                    Config::RESOLUTION_A_DEFAULT,
                                    Config::RESOLUTION_B_DEFAULT,
                                    Config::RESOLUTION_C_DEFAULT);

        // Scale to data integrals
        auto scaleToData = [](TH1D* hs, const TH1D* hd) {
            if (hs->Integral() > 0 && hd->Integral() > 0)
                hs->Scale(hd->Integral() / hs->Integral());
        };
        scaleToData(hu_c_mgg,   hd_c_mgg);
        scaleToData(hu_c_mmiss, hd_c_mmiss);
        scaleToData(hu_c_mpgg2, hd_c_mpgg2);
        scaleToData(hs_c_mgg,   hd_c_mgg);
        scaleToData(hs_c_mmiss, hd_c_mmiss);
        scaleToData(hs_c_mpgg2, hd_c_mpgg2);

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
                                        1.0, 0.0, 0.0,
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
            delete pt_attr;
            delete h_attr_map;
            delete c_attr;
        }

        // ---- Summary page: M_γγ pull distribution (all sections combined) ----
        {
            TCanvas *c_pull = new TCanvas("c_pull_summary", "M_gg Pull — All Sections", 1200, 550);
            c_pull->Divide(2, 1);

            c_pull->cd(1);
            gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.14);
            {
                int nb = hd_c_mgg->GetNbinsX();
                TH1D hpull_c("hpull_comb", "(Data#minusSmeared)/#sigma  (All Sections);M_{#gamma#gamma} [GeV/c^{2}];Pull",
                            nb, hd_c_mgg->GetXaxis()->GetXmin(), hd_c_mgg->GetXaxis()->GetXmax());
                for (int ib = 1; ib <= nb; ++ib) {
                    double d = hd_c_mgg->GetBinContent(ib), s = hs_c_mgg->GetBinContent(ib);
                    double ed = hd_c_mgg->GetBinError(ib), es = hs_c_mgg->GetBinError(ib);
                    double denom = sqrt(ed*ed + es*es);
                    if (denom > 0) hpull_c.SetBinContent(ib, (d - s) / denom);
                }
                hpull_c.SetFillColor(kCyan-9); hpull_c.SetLineColor(kBlack);
                double pm = max(fabs(hpull_c.GetMinimum()), fabs(hpull_c.GetMaximum()));
                hpull_c.SetMaximum( max(pm*1.2, 3.0)); hpull_c.SetMinimum(-max(pm*1.2, 3.0));
                hpull_c.Draw("HIST");
                TLine lz(hd_c_mgg->GetXaxis()->GetXmin(), 0, hd_c_mgg->GetXaxis()->GetXmax(), 0);
                lz.SetLineColor(kRed); lz.SetLineWidth(2); lz.Draw();
            }
            c_pull->cd(2);
            gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.14);
            {
                // Pull distribution histogram (Gaussian check)
                int nb = hd_c_mgg->GetNbinsX();
                TH1D hpull_dist("hpull_dist", "Pull distribution;Pull value;Bins",
                                30, -5.0, 5.0);
                for (int ib = 1; ib <= nb; ++ib) {
                    double d = hd_c_mgg->GetBinContent(ib), s = hs_c_mgg->GetBinContent(ib);
                    double ed = hd_c_mgg->GetBinError(ib), es = hs_c_mgg->GetBinError(ib);
                    double denom = sqrt(ed*ed + es*es);
                    if (denom > 0) hpull_dist.Fill((d - s) / denom);
                }
                hpull_dist.SetFillColor(kOrange-9); hpull_dist.SetLineColor(kBlack);
                hpull_dist.Draw("HIST");
                TLatex tx; tx.SetNDC(); tx.SetTextSize(0.045);
                tx.DrawLatex(0.15, 0.87, "Pull value distribution");
                tx.DrawLatex(0.15, 0.81, "(expect: Gaussian, #mu=0, #sigma=1)");
            }
            c_pull->Print(pdf_file.c_str());
            delete c_pull;
        }

        // ---- Summary page: 2D parameter maps (mu, sigma, sigma_pos, p_e_scale, chi2/ndf) ----
        {
            TCanvas *c_maps = new TCanvas("c_param_maps", "Fit Parameter Maps", 1800, 1000);
            c_maps->Divide(3, 2);

            double base_wx_map = (x_max - x_min) / nx;
            double base_wy_map = (y_max - y_min) / ny;

            TH2D *h_mu_map    = new TH2D("h_mu_map_diag",    "#mu per section;x [cm];y [cm]",
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
                h_mu_map->SetBinContent(bx, by, fit_results[is].mu);
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
            drawMap(h_mu_map,     1, c_maps);
            drawMap(h_sig_map,    2, c_maps);
            drawMap(h_sigpos_map, 3, c_maps);
            drawMap(h_pe_map,     4, c_maps);
            drawMap(h_chi2_map,   5, c_maps);

            c_maps->Print(pdf_file.c_str());
            delete c_maps;
            delete h_mu_map; delete h_sig_map; delete h_sigpos_map; delete h_pe_map; delete h_chi2_map;
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
            txt->AddText(Form("Optimization mode: %s", Config::USE_THREE_STAGE_OPTIMIZATION ? "three-stage" : "simultaneous"));
            txt->AddText(Form("Observable mode: %s", Config::histogram_mode_label()));
            txt->AddText(Form("Weights: w_{M_{#gamma#gamma}}=%.3f, w_{M_{miss}}=%.3f, w_{(p+#gamma#gamma)^{2}}=%.3f",
                            Config::W_MPI0, Config::W_MMISS, Config::W_MPGG2));
            txt->AddText(Form("Energy model: %s", Config::USE_SIMPLE_STOCHASTIC_MODEL ? "#sigma_{E}=#sigma#sqrt{E}" : "3-term resolution model"));
            txt->AddText(Form("Position smearing: %s", Config::ENABLE_POSITION_SMEARING ? "enabled" : "disabled"));
            txt->AddText(Form("N_{smear} per event: %d", Nsmear));
            txt->AddText(Form("Calorimeter grid: nx=%d, ny=%d, x=[%.1f, %.1f] cm, y=[%.1f, %.1f] cm",
                            nx, ny, x_min, x_max, y_min, y_max));
            txt->AddText(" ");
            txt->AddText(Form("Electron momentum scaling mode: %s",
                            use_global_p_e_scale
                                ? (Config::ENABLE_FINAL_GLOBAL_PE_SCALE_STAGE4
                                        ? "per-section Pass-2 + final global Stage-4"
                                        : "per-section Pass-2 only (global Stage-4 disabled)")
                                : "disabled"));
            if (use_global_p_e_scale) {
                if (Config::ENABLE_FINAL_GLOBAL_PE_SCALE_STAGE4) {
                    txt->AddText(Form("Global fitted p_{e} scale (Stage-4): %.6f", global_p_e_scale));
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
            txt->AddText(Form("PDF content: %d section pages + 6 summary pages (this page is the final text report)",
                            page_count));
            txt->AddText(" ");
            txt->AddText("Note: all-sections histograms use unique global event buffers with BOTH photons in geometry.");
            txt->AddText("All-sections smearing assigns coefficients per photon from its own base-grid section.");

            txt->Draw();
            c_final->Print(pdf_file.c_str());

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

        fout.Close();
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
                            fit_results[is].mu,
                            fit_results[is].sigma,
                            fit_results[is].sigma_pos);
            }
        }
        
        // Save 2D interpolated maps for visualization
        string interp_file = out_file;
        size_t dot_pos = interp_file.find_last_of('.');
        if (dot_pos != string::npos) {
            interp_file.insert(dot_pos, Config::INTERPOLATED_SUFFIX);
        } else {
            interp_file += Config::INTERPOLATED_SUFFIX + ".root";
        }
        calMap.saveAsHistogram(interp_file);
        
        // Demonstrate usage: print interpolated values at a few test points
        cout << "\n==== Example interpolated values ====\n";
        calMap.printParamsAt(0.0, 0.0);   // center
        calMap.printParamsAt(-15.0, -15.0); // lower-left quadrant
        calMap.printParamsAt(15.0, 15.0);   // upper-right quadrant
        calMap.printParamsAt(5.5, -8.3);    // arbitrary point
        
        cout << "\n==== Usage for event-by-event correction ====\n";
        cout << "// Example: Apply smearing to a photon at position (x, y):\n";
        cout << "double mu, sigma, sigma_pos;\n";
        cout << "calMap.getInterpolatedParams(cluster_x, cluster_y, mu, sigma, sigma_pos);\n";
        if (Config::ENABLE_ENERGY_DEPENDENT_MU) {
            cout << "double E_safe = std::max(original_energy, " << Config::MU_ENERGY_MIN_GEV << ");\n";
            cout << "double mu_eff = mu * (pow(" << Config::MU_ENERGY_E0_GEV
                << "/(E_safe+" << Config::MU_ENERGY_SHIFT_GEV << "),"
                << Config::MU_ENERGY_B_EXPONENT << ") + 1.0);\n";
        }
        if (Config::USE_SIMPLE_STOCHASTIC_MODEL) {
            if (Config::ENABLE_ENERGY_DEPENDENT_MU) {
                cout << "double E_sc = mu_eff * original_energy;\n";
                cout << "double smeared_energy = rng.Gaus(E_sc, sigma * sqrt(E_sc));\n";
            } else {
                cout << "double smeared_energy = rng.Gaus(mu * original_energy, sigma * sqrt(mu * original_energy));\n";
            }
        } else {
            if (Config::ENABLE_ENERGY_DEPENDENT_MU) {
                cout << "double E_sc = mu_eff * original_energy;\n";
            } else {
                cout << "double E_sc = mu * original_energy;\n";
            }
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
