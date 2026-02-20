// ============================================================================
// File: nps_sim_smearing.C
// Purpose: NPS π0 simulation smearing analysis with per-section optimization
//
// This code performs chi2 minimization to find optimal smearing parameters
// (mu, sigma) for each section of the NPS calorimeter by comparing
// data and simulation invariant mass distributions.
//
// Uses physics calculations exactly as in nps_analysis.C
//
// Features:
//   - Per-section optimization of energy scale (mu) and resolution (sigma)
//   - Bilinear interpolation for smooth parameter variation across calorimeter
//   - Outputs: discrete section parameters (CSV) + interpolated 2D maps (ROOT)
//
// Compile:
//   g++ nps_sim_smearing.C `root-config --cflags --libs` -O2 -std=c++17 -fopenmp -I../src -o nps_sim_smearing
//
// Usage example: 
//   ./nps_sim_smearing /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/combined_branches_LH2.root physics /w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/simc_pi0_analysis_output.root simulation out_smear.root 8 8 -30 30 -36 36 0.2 20
//
// Output files:
//   - out.root: Per-section histograms and fit parameters
//   - out_interpolated.root: 2D maps (h_mu_interp, h_sigma_interp)
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
#include "TNamed.h"
#include "TDirectory.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMarker.h"

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

// ----- GLOBAL SETTINGS -----
const double NONPOSITIVE_CLAMP = 1e-6; // clamp for non-positive smearing draws

// ============================================================================
// ENERGY RESOLUTION MODEL SELECTION
// ============================================================================
// Choose between two resolution models:
//
// MODEL 1 (Simple Stochastic): σ_E = σ × √E
//   - Data-driven, 2 parameters per section (μ, σ)
//   - σ typically 0.02-0.04 GeV^(1/2) for calorimeters
//   - Fast, good for empirical calibration
//
// MODEL 2 (Full 3-term): σ_E = σ × E × √((A/√E)² + (B/E)² + C²)
//   - Physics-motivated: stochastic + noise + constant terms
//   - Requires A,B,C from external measurements or fit
//   - σ ≈ 1.0 if A,B,C are correct for your detector
//
const bool USE_SIMPLE_STOCHASTIC_MODEL = true;  // true = Model 1, false = Model 2

// Model 2 constants (only used if USE_SIMPLE_STOCHASTIC_MODEL = false)
// Reference values from Wasim's NIM paper draft (E in GeV, ⊕ means quadrature sum):
//   RESOLUTION_A = 0.0097  (Stochastic term, 0.97%)
//   RESOLUTION_B = 0.011   (Noise term, 1.1%)
//   RESOLUTION_C = 0.0114  (Constant term, 1.14%)
// When using 3-term model, these are FITTED parameters, not fixed constants.
const double RESOLUTION_A_DEFAULT = 0.0097;  // Stochastic term starting value
const double RESOLUTION_B_DEFAULT = 0.011;   // Noise term starting value
const double RESOLUTION_C_DEFAULT = 0.0114;  // Constant term starting value
// ============================================================================

// rotate vector v around axis a_unit by angle theta (radians)
TVector3 rotateAroundAxis(const TVector3 &v, const TVector3 &a_unit, double theta) {
    TVector3 term1 = v * cos(theta);
    TVector3 term2 = a_unit.Cross(v) * sin(theta);
    TVector3 term3 = a_unit * (a_unit.Dot(v)) * (1.0 - cos(theta));
    return term1 + term2 + term3;
}

// compute chi2 using formula: sum ( (s-d)^2 / (s+d) ) skipping bins with denom<=0
double computeChi2FromHist(const TH1D &hsim, const TH1D &hdata) {
    if (hsim.GetNbinsX() != hdata.GetNbinsX() ||
        fabs(hsim.GetXaxis()->GetXmin() - hdata.GetXaxis()->GetXmin()) > 1e-9 ||
        fabs(hsim.GetXaxis()->GetXmax() - hdata.GetXaxis()->GetXmax()) > 1e-9) {
        cerr << "Histogram binning mismatch in computeChi2FromHist" << endl;
        return 1e300;
    }
    double chi2 = 0.;
    int nb = hsim.GetNbinsX();
    for (int i = 1; i <= nb; ++i) {
        double s = hsim.GetBinContent(i);
        double d = hdata.GetBinContent(i);
        double denom = s + d;
        if (denom <= 0) continue;
        double diff = s - d;
        chi2 += (diff * diff) / denom;
    }
    return chi2;
}

// Section structure
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
};

// ============================================================================
// CalibrationMap: Provides interpolated smearing parameters across calorimeter
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
    
    // Get linear index from (ix, iy)
    int getIndex(int ix, int iy) const {
        return iy * nx_ + ix;
    }
    
    // Clamp value to range [min, max]
    double clamp(double val, double min, double max) const {
        return std::max(min, std::min(max, val));
    }

public:
    CalibrationMap(int nx, int ny, double x_min, double x_max, double y_min, double y_max)
        : nx_(nx), ny_(ny), x_min_(x_min), x_max_(x_max), y_min_(y_min), y_max_(y_max)
    {
        double total_wx = x_max - x_min;
        double total_wy = y_max - y_min;
        base_wx_ = total_wx / nx;
        base_wy_ = total_wy / ny;
        
        int nsec = nx * ny;
        centers_x_.resize(nsec);
        centers_y_.resize(nsec);
        mu_values_.resize(nsec, 1.0);      // default values
        sigma_values_.resize(nsec, 0.05);
        
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
    void setParams(int ix, int iy, double mu, double sigma) {
        int idx = getIndex(ix, iy);
        mu_values_[idx] = mu;
        sigma_values_[idx] = sigma;
    }
    
    // Get interpolated parameters at position (x, y) using bilinear interpolation
    void getInterpolatedParams(double x, double y, double &mu, double &sigma) const {
        // Clamp position to calorimeter bounds
        x = clamp(x, x_min_, x_max_);
        y = clamp(y, y_min_, y_max_);
        
        // Find fractional position in grid
        double fx = (x - x_min_) / base_wx_;  // fractional x index
        double fy = (y - y_min_) / base_wy_;  // fractional y index
        
        // Get integer indices of surrounding grid points
        int ix0 = (int)floor(fx);
        int iy0 = (int)floor(fy);
        
        // Clamp to valid section indices
        ix0 = max(0, min(nx_ - 1, ix0));
        iy0 = max(0, min(ny_ - 1, iy0));
        
        int ix1 = min(nx_ - 1, ix0 + 1);
        int iy1 = min(ny_ - 1, iy0 + 1);
        
        // Get fractional parts for interpolation
        double tx = fx - ix0;  // 0 to 1 within cell
        double ty = fy - iy0;
        
        // Clamp interpolation weights
        tx = clamp(tx, 0.0, 1.0);
        ty = clamp(ty, 0.0, 1.0);
        
        // Get values at four corners
        int idx00 = getIndex(ix0, iy0);
        int idx10 = getIndex(ix1, iy0);
        int idx01 = getIndex(ix0, iy1);
        int idx11 = getIndex(ix1, iy1);
        
        // Bilinear interpolation for mu
        double mu00 = mu_values_[idx00];
        double mu10 = mu_values_[idx10];
        double mu01 = mu_values_[idx01];
        double mu11 = mu_values_[idx11];
        mu = (1-tx)*(1-ty)*mu00 + tx*(1-ty)*mu10 + (1-tx)*ty*mu01 + tx*ty*mu11;
        
        // Bilinear interpolation for sigma
        double sig00 = sigma_values_[idx00];
        double sig10 = sigma_values_[idx10];
        double sig01 = sigma_values_[idx01];
        double sig11 = sigma_values_[idx11];
        sigma = (1-tx)*(1-ty)*sig00 + tx*(1-ty)*sig10 + (1-tx)*ty*sig01 + tx*ty*sig11;
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
            // tokens[2-7] are bounds and counts
            
            // Check if fit succeeded (has valid parameters)
            if (tokens[8].empty() || tokens[9].empty() || tokens[10].empty()) {
                // No valid fit for this section - keep defaults
                continue;
            }
            
            double mu = atof(tokens[8].c_str());
            double sigma = atof(tokens[9].c_str());
            
            setParams(ix, iy, mu, sigma);
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
        
        for (int ix = 1; ix <= nbinsx; ++ix) {
            for (int iy = 1; iy <= nbinsy; ++iy) {
                double x = h_mu->GetXaxis()->GetBinCenter(ix);
                double y = h_mu->GetYaxis()->GetBinCenter(iy);
                
                double mu, sigma;
                getInterpolatedParams(x, y, mu, sigma);
                
                h_mu->SetBinContent(ix, iy, mu);
                h_sigma->SetBinContent(ix, iy, sigma);
            }
        }
        
        h_mu->Write();
        h_sigma->Write();
        fout.Close();
        
        cout << "Saved interpolated calibration maps to " << filename << endl;
    }
    
    // Print parameters at a specific position
    void printParamsAt(double x, double y) const {
        double mu, sigma;
        getInterpolatedParams(x, y, mu, sigma);
        cout << "Calibration at (" << x << ", " << y << "): "
             << "mu=" << mu << ", sigma=" << sigma << endl;
    }
};

// Evaluate chi2 for (mu, sigma) for a given set of sim events that belong to a section
// Evaluate chi2 with energy smearing only (optimized version)
double eval_chi2_section(double mu, double sigma,
                          const vector<ClusterPair> &simEvents,
                          const TH1D &hdata,
                          TRandom3 &rng,
                          int Nsmear,
                          double res_A, double res_B, double res_C) {
    TH1D hsim("hsim_tmp","smeared sim (tmp)", hdata.GetNbinsX(), hdata.GetXaxis()->GetXmin(), hdata.GetXaxis()->GetXmax());

    // Pre-allocate random number arrays for better cache performance
    vector<double> random_draws(Nsmear * 2);
    
    // Pre-compute constants for 3-term model if needed
    const double A_sq = USE_SIMPLE_STOCHASTIC_MODEL ? 0.0 : (res_A * res_A);
    const double B_sq = USE_SIMPLE_STOCHASTIC_MODEL ? 0.0 : (res_B * res_B);
    const double C_sq = USE_SIMPLE_STOCHASTIC_MODEL ? 0.0 : (res_C * res_C);
    
    for (const auto &ev : simEvents) {
        double E1 = ev.e1;
        double E2 = ev.e2;
        if (E1 <= 0 || E2 <= 0) continue;
        
        // Apply energy scale calibration
        double E1_scaled = mu * E1;
        double E2_scaled = mu * E2;
        
        // Compute energy-dependent resolution based on selected model
        double sigma_E1, sigma_E2;
        if (USE_SIMPLE_STOCHASTIC_MODEL) {
            // Simple stochastic: σ_E = σ × √E
            sigma_E1 = sigma * sqrt(E1);
            sigma_E2 = sigma * sqrt(E2);
        } else {
            // Full 3-term: σ_E = σ × E × √((A/√E)² + (B/E)² + C²)
            double sigma_E1_rel_sq = A_sq/E1_scaled + B_sq/(E1_scaled*E1_scaled) + C_sq;
            double sigma_E2_rel_sq = A_sq/E2_scaled + B_sq/(E2_scaled*E2_scaled) + C_sq;
            sigma_E1 = sigma * E1_scaled * sqrt(sigma_E1_rel_sq);
            sigma_E2 = sigma * E2_scaled * sqrt(sigma_E2_rel_sq);
        }
        
        // Pre-compute photon direction vectors (reused for all smears)
        double r1_sq = ev.x1*ev.x1 + ev.y1*ev.y1 + nps::kDefaultZ_NPS_cm*nps::kDefaultZ_NPS_cm;
        double r1_inv = 1.0 / sqrt(r1_sq);
        double u1x = ev.x1 * r1_inv;
        double u1y = ev.y1 * r1_inv;
        double u1z = nps::kDefaultZ_NPS_cm * r1_inv;
        
        double r2_sq = ev.x2*ev.x2 + ev.y2*ev.y2 + nps::kDefaultZ_NPS_cm*nps::kDefaultZ_NPS_cm;
        double r2_inv = 1.0 / sqrt(r2_sq);
        double u2x = ev.x2 * r2_inv;
        double u2y = ev.y2 * r2_inv;
        double u2z = nps::kDefaultZ_NPS_cm * r2_inv;
        
        // Batch generate all random numbers for this event with energy-dependent resolution
        for (int k = 0; k < Nsmear * 2; k += 2) {
            random_draws[k] = rng.Gaus(E1_scaled, sigma_E1);
            if (random_draws[k] <= 0) random_draws[k] = NONPOSITIVE_CLAMP;
            
            random_draws[k + 1] = rng.Gaus(E2_scaled, sigma_E2);
            if (random_draws[k + 1] <= 0) random_draws[k + 1] = NONPOSITIVE_CLAMP;
        }
        
        double event_weight = 1.0;
        for (int k = 0; k < Nsmear; ++k) {
            double E1p = random_draws[k*2];
            double E2p = random_draws[k*2 + 1];
            
            // Direct invariant mass calculation (avoid TLorentzVector overhead)
            // p1 = (E1p, E1p*u1x, E1p*u1y, E1p*u1z)
            // p2 = (E2p, E2p*u2x, E2p*u2y, E2p*u2z)
            // m^2 = (E1+E2)^2 - (px1+px2)^2 - (py1+py2)^2 - (pz1+pz2)^2
            double E_tot = E1p + E2p;
            double px_tot = E1p*u1x + E2p*u2x;
            double py_tot = E1p*u1y + E2p*u2y;
            double pz_tot = E1p*u1z + E2p*u2z;
            
            double m_sq = E_tot*E_tot - px_tot*px_tot - py_tot*py_tot - pz_tot*pz_tot;
            double mass = (m_sq > 0) ? sqrt(m_sq) : 0.0;
            
            hsim.Fill(mass, event_weight);
        }
    }

    double integral_sim = hsim.Integral();
    double integral_data = hdata.Integral();
    if (integral_sim <= 0) return 1e300;
    double scale = (integral_data > 0) ? (integral_data / integral_sim) : 1.0;
    hsim.Scale(scale);
    double chi2 = computeChi2FromHist(hsim, hdata);
    return chi2;
}

// Evaluate chi2 with both energy and position smearing (optimized version)
double eval_chi2_section_with_pos(double mu, double sigma, double sigma_pos,
                                    const vector<ClusterPair> &simEvents,
                                    const TH1D &hdata,
                                    TRandom3 &rng,
                                    int Nsmear,
                                    double res_A, double res_B, double res_C) {
    TH1D hsim("hsim_tmp_pos","smeared sim with position (tmp)", hdata.GetNbinsX(), hdata.GetXaxis()->GetXmin(), hdata.GetXaxis()->GetXmax());

    // Pre-allocate random number arrays (2 energy + 4 position per smear)
    vector<double> random_energy(Nsmear * 2);
    vector<double> random_pos(Nsmear * 4);
    
    const double z_nps = nps::kDefaultZ_NPS_cm;
    const double z_sq = z_nps * z_nps;
    
    // Pre-compute constants for 3-term model if needed
    const double A_sq = USE_SIMPLE_STOCHASTIC_MODEL ? 0.0 : (res_A * res_A);
    const double B_sq = USE_SIMPLE_STOCHASTIC_MODEL ? 0.0 : (res_B * res_B);
    const double C_sq = USE_SIMPLE_STOCHASTIC_MODEL ? 0.0 : (res_C * res_C);
    
    for (const auto &ev : simEvents) {
        double E1 = ev.e1;
        double E2 = ev.e2;
        if (E1 <= 0 || E2 <= 0) continue;
        
        // Apply energy scale calibration first
        double E1_scaled = mu * E1;
        double E2_scaled = mu * E2;
        
        // Compute energy-dependent resolution based on selected model
        double sigma_E1, sigma_E2;
        if (USE_SIMPLE_STOCHASTIC_MODEL) {
            // Simple stochastic: σ_E = σ × √E
            sigma_E1 = sigma * sqrt(E1);
            sigma_E2 = sigma * sqrt(E2);
        } else {
            // Full 3-term: σ_E = σ × E × √((A/√E)² + (B/E)² + C²)
            double sigma_E1_rel_sq = A_sq/E1_scaled + B_sq/(E1_scaled*E1_scaled) + C_sq;
            double sigma_E2_rel_sq = A_sq/E2_scaled + B_sq/(E2_scaled*E2_scaled) + C_sq;
            sigma_E1 = sigma * E1_scaled * sqrt(sigma_E1_rel_sq);
            sigma_E2 = sigma * E2_scaled * sqrt(sigma_E2_rel_sq);
        }
        
        // Batch generate all random numbers for this event
        // Batch generate all random numbers for this event
        for (int k = 0; k < Nsmear * 2; k += 2) {
            random_energy[k] = rng.Gaus(E1_scaled, sigma_E1);
            if (random_energy[k] <= 0) random_energy[k] = NONPOSITIVE_CLAMP;
            
            random_energy[k + 1] = rng.Gaus(E2_scaled, sigma_E2);
            if (random_energy[k + 1] <= 0) random_energy[k + 1] = NONPOSITIVE_CLAMP;
        }
        
        if (sigma_pos > 0) {
            for (int k = 0; k < Nsmear * 4; ++k) {
                random_pos[k] = rng.Gaus(0, sigma_pos);
            }
        }
        
        double event_weight = 1.0;
        for (int k = 0; k < Nsmear; ++k) {
            // 1. Energy smearing (already applied in random_energy generation)
            double E1_smeared = random_energy[k*2];
            double E2_smeared = random_energy[k*2 + 1];
            
            // 2. Position smearing
            double x1_smeared, y1_smeared, x2_smeared, y2_smeared;
            if (sigma_pos > 0) {
                x1_smeared = ev.x1 + random_pos[k*4];
                y1_smeared = ev.y1 + random_pos[k*4 + 1];
                x2_smeared = ev.x2 + random_pos[k*4 + 2];
                y2_smeared = ev.y2 + random_pos[k*4 + 3];
            } else {
                x1_smeared = ev.x1;
                y1_smeared = ev.y1;
                x2_smeared = ev.x2;
                y2_smeared = ev.y2;
            }
            
            // 3. Compute photon direction vectors and invariant mass directly
            double r1_sq = x1_smeared*x1_smeared + y1_smeared*y1_smeared + z_sq;
            double r1_inv = 1.0 / sqrt(r1_sq);
            double px1 = E1_smeared * x1_smeared * r1_inv;
            double py1 = E1_smeared * y1_smeared * r1_inv;
            double pz1 = E1_smeared * z_nps * r1_inv;
            
            double r2_sq = x2_smeared*x2_smeared + y2_smeared*y2_smeared + z_sq;
            double r2_inv = 1.0 / sqrt(r2_sq);
            double px2 = E2_smeared * x2_smeared * r2_inv;
            double py2 = E2_smeared * y2_smeared * r2_inv;
            double pz2 = E2_smeared * z_nps * r2_inv;
            
            // Direct invariant mass calculation
            double E_tot = E1_smeared + E2_smeared;
            double px_tot = px1 + px2;
            double py_tot = py1 + py2;
            double pz_tot = pz1 + pz2;
            
            double m_sq = E_tot*E_tot - px_tot*px_tot - py_tot*py_tot - pz_tot*pz_tot;
            double mass = (m_sq > 0) ? sqrt(m_sq) : 0.0;
            
            hsim.Fill(mass, event_weight);
        }
    }

    double integral_sim = hsim.Integral();
    double integral_data = hdata.Integral();
    if (integral_sim <= 0) return 1e300;
    double scale = (integral_data > 0) ? (integral_data / integral_sim) : 1.0;
    hsim.Scale(scale);
    double chi2 = computeChi2FromHist(hsim, hdata);
    return chi2;
}

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

// Minimization strategy per section -> returns best (mu,sigma,sigma_pos,chi2)
struct FitResult { 
    double mu, sigma, sigma_pos, chi2;
    // 3-term model parameters (only used if USE_SIMPLE_STOCHASTIC_MODEL = false)
    double res_A, res_B, res_C;
    int ndf;
    Chi2ScanData scan_data;
    double chi2_per_ndf() const { return (ndf > 0) ? chi2 / ndf : 0.0; }
};

FitResult fit_section(const vector<ClusterPair> &simEvents,
                      const TH1D &hdata,
                      TRandom3 &rng,
                      int Nsmear) {
    FitResult res; 
    res.mu = 1.0; 
    res.sigma = 0.05;
    res.sigma_pos = 0.0;
    res.res_A = RESOLUTION_A_DEFAULT;
    res.res_B = RESOLUTION_B_DEFAULT;
    res.res_C = RESOLUTION_C_DEFAULT;
    res.chi2 = 1e300;
    // ndf depends on model: simple=2 params (μ,σ), 3-term=5 params (μ,σ,A,B,C)
    int n_params = USE_SIMPLE_STOCHASTIC_MODEL ? 2 : 5;
    res.ndf = hdata.GetNbinsX() - n_params;
    
    // ========================================================================
    // USER CONFIGURATION: Scan ranges and resolution
    // ========================================================================
    
    // Energy-dependent resolution model: σ_E/E = A/√E ⊕ B/E ⊕ C
    // A, B, C constants defined globally (see top of file)
    
    // Energy scale (mu) scan range and resolution
    const double MU_MIN = 1.015;
    const double MU_MAX = 1.04;
    const int MU_NSTEPS = 60;  // Number of points to scan
    
    // Resolution parameter (sigma) scan range
    // Units depend on model (see USE_SIMPLE_STOCHASTIC_MODEL at top of file):
    //   Simple model: σ in GeV^(1/2), typically 0.02-0.04
    //   3-term model: σ dimensionless, typically 0.5-1.5
    const double SIGMA_MIN = USE_SIMPLE_STOCHASTIC_MODEL ? 0.01 : 0.5;
    const double SIGMA_MAX = USE_SIMPLE_STOCHASTIC_MODEL ? 0.09 : 1.5;
    const int SIGMA_NSTEPS = 60;  // Number of points to scan
    
    // 3-term model coefficient scan ranges (only used if USE_SIMPLE_STOCHASTIC_MODEL = false)
    // These are fractional values (A=0.01 means 1%)
    const double RES_A_MIN = 0.005;  // Stochastic: 0.5% to 2%
    const double RES_A_MAX = 0.020;
    const int RES_A_NSTEPS = 15;
    
    const double RES_B_MIN = 0.005;  // Noise: 0.5% to 2%
    const double RES_B_MAX = 0.020;
    const int RES_B_NSTEPS = 15;
    
    const double RES_C_MIN = 0.005;  // Constant: 0.5% to 2.5%
    const double RES_C_MAX = 0.025;
    const int RES_C_NSTEPS = 15;
    
    // Position smearing toggle and parameters
    const bool ENABLE_POSITION_SMEARING = false;  // Set to false to disable position smearing
    const double SIGMA_POS_MIN = 0.0;  // cm
    const double SIGMA_POS_MAX = 0.1;  // cm
    const int SIGMA_POS_NSTEPS = 60;  // Number of points to scan
    // ========================================================================
    
    // Calculate step sizes from number of steps
    double mu_2d_step = (MU_MAX - MU_MIN) / (MU_NSTEPS - 1);
    double sigma_2d_step = (SIGMA_MAX - SIGMA_MIN) / (SIGMA_NSTEPS - 1);
    double sigma_pos_step = (SIGMA_POS_MAX - SIGMA_POS_MIN) / (SIGMA_POS_NSTEPS - 1);
    
    // 1) Coarse 2D grid search to find global minimum region
    // NOTE: Using eval_chi2_section_with_pos with sigma_pos=0 for consistency
    double best_chi2 = 1e300;
    double mu0 = 1.0, sigma0 = 0.05;
    
    cout << "Performing coarse 2D grid search..." << endl;
    for (double mu = MU_MIN; mu <= MU_MAX + 1e-12; mu += mu_2d_step) {
        for (double sig = SIGMA_MIN; sig <= SIGMA_MAX + 1e-12; sig += sigma_2d_step) {
            double chi2 = eval_chi2_section_with_pos(mu, sig, 0.0, simEvents, hdata, rng, Nsmear, 
                                                     RESOLUTION_A_DEFAULT, RESOLUTION_B_DEFAULT, RESOLUTION_C_DEFAULT);
            // Store for visualization
            res.scan_data.mu_2d.push_back(mu);
            res.scan_data.sigma_2d.push_back(sig);
            res.scan_data.chi2_2d.push_back(chi2);
            
            if (chi2 < best_chi2) {
                best_chi2 = chi2;
                mu0 = mu;
                sigma0 = sig;
            }
        }
    }
    cout << "Coarse 2D minimum at mu=" << mu0 << ", sigma=" << sigma0 << ", chi2=" << best_chi2 << endl;

    // 2) Fine 2D refinement around global minimum from coarse search
    double twoD_step_mu = 0.005, twoD_step_sigma = 0.001;
    while (twoD_step_mu >= 0.0001 || twoD_step_sigma >= 0.0001) {
        double grid_best_mu = mu0, grid_best_sigma = sigma0;
        double grid_best_chi2 = best_chi2;
        for (double mu = mu0 - 2*twoD_step_mu; mu <= mu0 + 2*twoD_step_mu + 1e-15; mu += twoD_step_mu) {
            for (double sig = max(1e-6, sigma0 - 2*twoD_step_sigma); sig <= sigma0 + 2*twoD_step_sigma + 1e-15; sig += twoD_step_sigma) {
                if (sig <= 0) continue;
                double chi2 = eval_chi2_section_with_pos(mu, sig, 0.0, simEvents, hdata, rng, Nsmear, 
                                                         RESOLUTION_A_DEFAULT, RESOLUTION_B_DEFAULT, RESOLUTION_C_DEFAULT);
                if (chi2 < grid_best_chi2) {
                    grid_best_chi2 = chi2; grid_best_mu = mu; grid_best_sigma = sig;
                }
            }
        }
        mu0 = grid_best_mu; sigma0 = grid_best_sigma; best_chi2 = grid_best_chi2;
        twoD_step_mu /= 2.0; twoD_step_sigma /= 2.0;
        if (twoD_step_mu < 1e-6 && twoD_step_sigma < 1e-6) break;
    }
    cout << "Refined 2D minimum at mu=" << mu0 << ", sigma=" << sigma0 << ", chi2=" << best_chi2 << endl;

    // 3) Generate 1D slices through optimum for visualization
    cout << "Generating 1D slices for visualization..." << endl;
    int mu_vis_step = max(1, MU_NSTEPS / 20);  // ~20 points for visualization
    for (int i = 0; i < MU_NSTEPS; i += mu_vis_step) {
        double mu = MU_MIN + i * mu_2d_step;
        double chi2 = eval_chi2_section_with_pos(mu, sigma0, 0.0, simEvents, hdata, rng, Nsmear, 
                                                 RESOLUTION_A_DEFAULT, RESOLUTION_B_DEFAULT, RESOLUTION_C_DEFAULT);
        res.scan_data.mu_values.push_back(mu);
        res.scan_data.chi2_vs_mu.push_back(chi2);
    }
    
    int sigma_vis_step = max(1, SIGMA_NSTEPS / 20);  // ~20 points for visualization
    for (int i = 0; i < SIGMA_NSTEPS; i += sigma_vis_step) {
        double sigma = SIGMA_MIN + i * sigma_2d_step;
        double chi2 = eval_chi2_section_with_pos(mu0, sigma, 0.0, simEvents, hdata, rng, Nsmear, 
                                                 RESOLUTION_A_DEFAULT, RESOLUTION_B_DEFAULT, RESOLUTION_C_DEFAULT);
        res.scan_data.sigma_values.push_back(sigma);
        res.scan_data.chi2_vs_sigma.push_back(chi2);
    }

    // 3.5) Optimize A, B, C parameters (only for 3-term model)
    double best_res_A = RESOLUTION_A_DEFAULT;
    double best_res_B = RESOLUTION_B_DEFAULT;
    double best_res_C = RESOLUTION_C_DEFAULT;
    
    if (!USE_SIMPLE_STOCHASTIC_MODEL) {
        cout << "Optimizing 3-term resolution parameters (A, B, C)..." << endl;
        double res_A_step = (RES_A_MAX - RES_A_MIN) / (RES_A_NSTEPS - 1);
        double res_B_step = (RES_B_MAX - RES_B_MIN) / (RES_B_NSTEPS - 1);
        double res_C_step = (RES_C_MAX - RES_C_MIN) / (RES_C_NSTEPS - 1);
        
        // Optimize one parameter at a time to avoid 5D search
        // Iteration 1: Optimize A with fixed B, C
        double best_chi2_A = best_chi2;
        for (int i = 0; i < RES_A_NSTEPS; ++i) {
            double res_A = RES_A_MIN + i * res_A_step;
            double chi2 = eval_chi2_section_with_pos(mu0, sigma0, 0.0, simEvents, hdata, rng, Nsmear, 
                                                     res_A, best_res_B, best_res_C);
            if (chi2 < best_chi2_A) {
                best_chi2_A = chi2;
                best_res_A = res_A;
            }
        }
        
        // Iteration 2: Optimize B with updated A, fixed C
        double best_chi2_B = best_chi2_A;
        for (int i = 0; i < RES_B_NSTEPS; ++i) {
            double res_B = RES_B_MIN + i * res_B_step;
            double chi2 = eval_chi2_section_with_pos(mu0, sigma0, 0.0, simEvents, hdata, rng, Nsmear,
                                                     best_res_A, res_B, best_res_C);
            if (chi2 < best_chi2_B) {
                best_chi2_B = chi2;
                best_res_B = res_B;
            }
        }
        
        // Iteration 3: Optimize C with updated A, B
        double best_chi2_C = best_chi2_B;
        for (int i = 0; i < RES_C_NSTEPS; ++i) {
            double res_C = RES_C_MIN + i * res_C_step;
            double chi2 = eval_chi2_section_with_pos(mu0, sigma0, 0.0, simEvents, hdata, rng, Nsmear,
                                                     best_res_A, best_res_B, res_C);
            if (chi2 < best_chi2_C) {
                best_chi2_C = chi2;
                best_res_C = res_C;
            }
        }
        
        cout << "Optimized resolution parameters: A=" << best_res_A*100 << "%, B=" << best_res_B*100 
             << "%, C=" << best_res_C*100 << "%, chi2=" << best_chi2_C << endl;
        
        // Re-optimize (mu, sigma) with new A, B, C
        cout << "Re-optimizing (mu, sigma) with updated A,B,C..." << endl;
        double reopt_best_chi2 = best_chi2_C;
        double reopt_mu = mu0, reopt_sigma = sigma0;
        double reopt_step_mu = 0.002, reopt_step_sigma = 0.01;
        
        for (double mu = mu0 - 0.01; mu <= mu0 + 0.01 + 1e-15; mu += reopt_step_mu) {
            for (double sig = max(1e-6, sigma0 - 0.2); sig <= sigma0 + 0.2 + 1e-15; sig += reopt_step_sigma) {
                if (sig <= 0) continue;
                double chi2 = eval_chi2_section_with_pos(mu, sig, 0.0, simEvents, hdata, rng, Nsmear,
                                                         best_res_A, best_res_B, best_res_C);
                if (chi2 < reopt_best_chi2) {
                    reopt_best_chi2 = chi2;
                    reopt_mu = mu;
                    reopt_sigma = sig;
                }
            }
        }
        
        mu0 = reopt_mu;
        sigma0 = reopt_sigma;
        best_chi2 = reopt_best_chi2;
        cout << "Re-optimized: mu=" << mu0 << ", sigma=" << sigma0 << ", chi2=" << best_chi2 << endl;
    }
    
    // Store fitted resolution parameters
    res.res_A = best_res_A;
    res.res_B = best_res_B;
    res.res_C = best_res_C;

    // 4) Optimize position smearing parameter with current (mu, sigma)
    //    (Only if ENABLE_POSITION_SMEARING is true)
    double best_sigma_pos = 0.0;
    double best_chi2_pos = best_chi2;
    double mu_refine = mu0, sigma_refine = sigma0;
    double chi2_refine = best_chi2;
    
    if (ENABLE_POSITION_SMEARING) {
        // Use the chi2 from stage 2 as baseline (already at sigma_pos=0)
        cout << "Optimizing position smearing parameter..." << endl;
        cout << "Baseline chi2 at optimal (mu,sigma) with sigma_pos=0: " << best_chi2 << endl;
        
        for (double s_pos = SIGMA_POS_MIN; s_pos <= SIGMA_POS_MAX + 1e-12; s_pos += sigma_pos_step) {
            double chi2 = eval_chi2_section_with_pos(mu0, sigma0, s_pos, simEvents, hdata, rng, Nsmear, 
                                                     best_res_A, best_res_B, best_res_C);
            res.scan_data.sigma_pos_values.push_back(s_pos);
            res.scan_data.chi2_vs_sigma_pos.push_back(chi2);
            
            if (chi2 < best_chi2_pos) {
                best_chi2_pos = chi2;
                best_sigma_pos = s_pos;
            }
        }
        
        double chi2_improvement = best_chi2 - best_chi2_pos;
        cout << "Optimal sigma_pos=" << std::fixed << std::setprecision(5) << best_sigma_pos << " cm, chi2=" << best_chi2_pos;
        if (chi2_improvement > 0) {
            cout << " (improved by " << std::setprecision(2) << chi2_improvement << " from baseline)";
        } else {
            cout << " (no improvement, keeping sigma_pos=0)";
            best_sigma_pos = 0.0;
            best_chi2_pos = best_chi2;
        }
        cout << endl;

        // 5) Re-optimize (mu, sigma) with the new sigma_pos to handle parameter coupling
        // Only proceed if position smearing actually improved the fit
        chi2_refine = best_chi2_pos;
        
        if (best_sigma_pos > 0 && chi2_improvement > 0.5) {  // Require at least 0.5 chi2 improvement
            cout << "Re-optimizing (mu, sigma) with sigma_pos=" << best_sigma_pos << " cm..." << endl;
            
            // Fine 2D refinement with position smearing included (use consistent Nsmear)
            double refine_step_mu = 0.005, refine_step_sigma = 0.001;
            int max_refine_iter = 3;
            
            for (int iter = 0; iter < max_refine_iter; ++iter) {
                double iter_best_mu = mu_refine, iter_best_sigma = sigma_refine;
                double iter_best_chi2 = chi2_refine;
                
                for (double mu = mu_refine - 2*refine_step_mu; mu <= mu_refine + 2*refine_step_mu + 1e-15; mu += refine_step_mu) {
                    for (double sig = max(1e-6, sigma_refine - 2*refine_step_sigma); sig <= sigma_refine + 2*refine_step_sigma + 1e-15; sig += refine_step_sigma) {
                        if (sig <= 0) continue;
                        double chi2 = eval_chi2_section_with_pos(mu, sig, best_sigma_pos, simEvents, hdata, rng, Nsmear, 
                                                                 best_res_A, best_res_B, best_res_C);
                        if (chi2 < iter_best_chi2) {
                            iter_best_chi2 = chi2;
                            iter_best_mu = mu;
                            iter_best_sigma = sig;
                        }
                    }
                }
                
                // Check convergence
                if (fabs(iter_best_mu - mu_refine) < 1e-6 && fabs(iter_best_sigma - sigma_refine) < 1e-6) {
                    mu_refine = iter_best_mu;
                    mu_refine = iter_best_mu;
                    sigma_refine = iter_best_sigma;
                    chi2_refine = iter_best_chi2;
                    break;
                }
                
                mu_refine = iter_best_mu;
                sigma_refine = iter_best_sigma;
                chi2_refine = iter_best_chi2;
                refine_step_mu /= 2.0;
                refine_step_sigma /= 2.0;
            }
            
            cout << "Refined: mu=" << mu_refine << ", sigma=" << sigma_refine << ", chi2=" << chi2_refine << endl;
            
            // 6) One final sigma_pos optimization pass with refined (mu, sigma) if chi2 improved significantly
            if (chi2_refine < best_chi2_pos * 0.99) {  // Only if >1% improvement
                cout << "Final sigma_pos refinement..." << endl;
                double final_best_sigma_pos = best_sigma_pos;
                double final_best_chi2 = chi2_refine;
                
                // Refine around ±20% of current sigma_pos value using same step size
                double refine_range = max(0.2, best_sigma_pos * 0.2);  // At least 0.2 cm range
                for (double s_pos = max(0.0, best_sigma_pos - refine_range); s_pos <= best_sigma_pos + refine_range + 1e-12; s_pos += sigma_pos_step) {
                    double chi2 = eval_chi2_section_with_pos(mu_refine, sigma_refine, s_pos, simEvents, hdata, rng, Nsmear, 
                                                             best_res_A, best_res_B, best_res_C);
                    if (chi2 < final_best_chi2) {
                        final_best_chi2 = chi2;
                        final_best_sigma_pos = s_pos;
                    }
                }
                
                if (final_best_chi2 < chi2_refine) {
                    best_sigma_pos = final_best_sigma_pos;
                    chi2_refine = final_best_chi2;
                    cout << "Final optimal: sigma_pos=" << std::fixed << std::setprecision(5) << best_sigma_pos << " cm, chi2=" << chi2_refine << endl;
                }
            }
        } else {
            cout << "Skipping re-optimization (position smearing did not improve fit)" << endl;
        }
    
    } else {
        // Position smearing disabled - use stage 2 results directly
        cout << "Position smearing disabled (ENABLE_POSITION_SMEARING=false)" << endl;
        cout << "Using energy-only optimization results: mu=" << mu0 << ", sigma=" << sigma0 << ", chi2=" << best_chi2 << endl;
        // mu_refine, sigma_refine, chi2_refine already initialized to mu0, sigma0, best_chi2
        // best_sigma_pos already initialized to 0.0
    }

    res.mu = mu_refine; res.sigma = sigma_refine; res.sigma_pos = best_sigma_pos; res.chi2 = chi2_refine;
    // ndf = nbins - nparams; position smearing counts as parameter only if non-zero
    int nparams = (best_sigma_pos > 0) ? 3 : 2;  // (mu, sigma) or (mu, sigma, sigma_pos)
    res.ndf = hdata.GetNbinsX() - nparams;
    return res;
}

// Function to create chi^2 plots for a section
void plot_chi2_scans(const Section &sec, const FitResult &fres, TCanvas *c, int page_num, 
                     const TH1D &hdata, const TH1D &hsim_smeared, const TH1D &hsim_unsmeared) {
    c->Clear();
    c->Divide(3, 2);
    
    // Plot 1: chi^2 vs mu
    c->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    
    int n_mu = fres.scan_data.mu_values.size();
    if (n_mu > 0) {
        TGraph *gr_mu = new TGraph(n_mu);
        gr_mu->SetName(Form("gr_mu_%d_%d", sec.ix, sec.iy));
        for (int i = 0; i < n_mu; ++i) {
            gr_mu->SetPoint(i, fres.scan_data.mu_values[i], fres.scan_data.chi2_vs_mu[i]);
        }
        gr_mu->SetTitle(Form("Section %d,%d: #chi^{2} vs #mu (at optimal #sigma=%.4f)", sec.ix, sec.iy, fres.sigma));
        gr_mu->GetXaxis()->SetTitle("#mu (energy scale)");
        gr_mu->GetYaxis()->SetTitle("#chi^{2}");
        gr_mu->SetMarkerStyle(20);
        gr_mu->SetMarkerSize(0.8);
        gr_mu->SetMarkerColor(kBlue);
        gr_mu->SetLineColor(kBlue);
        gr_mu->Draw("APL");
        
        // Mark the best point
        TMarker *m_mu = new TMarker(fres.mu, fres.chi2, 29);
        m_mu->SetMarkerColor(kRed);
        m_mu->SetMarkerSize(2.0);
        m_mu->Draw();
        
        TLatex *t_mu = new TLatex();
        t_mu->SetNDC();
        t_mu->SetTextSize(0.04);
        t_mu->DrawLatex(0.15, 0.85, Form("Best #mu = %.4f", fres.mu));
        t_mu->DrawLatex(0.15, 0.80, Form("#chi^{2}/ndf = %.2f/%d", fres.chi2, fres.ndf));
    }
    
    // Plot 2: chi^2 vs sigma
    c->cd(2);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    
    int n_sigma = fres.scan_data.sigma_values.size();
    if (n_sigma > 0) {
        TGraph *gr_sigma = new TGraph(n_sigma);
        gr_sigma->SetName(Form("gr_sigma_%d_%d", sec.ix, sec.iy));
        for (int i = 0; i < n_sigma; ++i) {
            gr_sigma->SetPoint(i, fres.scan_data.sigma_values[i], fres.scan_data.chi2_vs_sigma[i]);
        }
        gr_sigma->SetTitle(Form("Section %d,%d: #chi^{2} vs #sigma (at optimal #mu=%.4f)", sec.ix, sec.iy, fres.mu));
        gr_sigma->GetXaxis()->SetTitle("#sigma (energy resolution)");
        gr_sigma->GetYaxis()->SetTitle("#chi^{2}");
        gr_sigma->SetMarkerStyle(20);
        gr_sigma->SetMarkerSize(0.8);
        gr_sigma->SetMarkerColor(kGreen+2);
        gr_sigma->SetLineColor(kGreen+2);
        gr_sigma->Draw("APL");
        
        // Mark the best point
        TMarker *m_sigma = new TMarker(fres.sigma, fres.chi2, 29);
        m_sigma->SetMarkerColor(kRed);
        m_sigma->SetMarkerSize(2.0);
        m_sigma->Draw();
        
        TLatex *t_sigma = new TLatex();
        t_sigma->SetNDC();
        t_sigma->SetTextSize(0.04);
        t_sigma->DrawLatex(0.15, 0.85, Form("Best #sigma = %.4f", fres.sigma));
        t_sigma->DrawLatex(0.15, 0.80, Form("#chi^{2}/ndf = %.2f/%d", fres.chi2, fres.ndf));
    }
    
    // Plot 3: 2D chi^2 map
    c->cd(3);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.15);
    
    int n_2d = fres.scan_data.mu_2d.size();
    if (n_2d > 0) {
        // Find ranges and determine bin structure from grid
        double mu_min_2d = *min_element(fres.scan_data.mu_2d.begin(), fres.scan_data.mu_2d.end());
        double mu_max_2d = *max_element(fres.scan_data.mu_2d.begin(), fres.scan_data.mu_2d.end());
        double sig_min_2d = *min_element(fres.scan_data.sigma_2d.begin(), fres.scan_data.sigma_2d.end());
        double sig_max_2d = *max_element(fres.scan_data.sigma_2d.begin(), fres.scan_data.sigma_2d.end());
        
        // Count unique values to determine number of bins
        std::set<double> unique_mu(fres.scan_data.mu_2d.begin(), fres.scan_data.mu_2d.end());
        std::set<double> unique_sigma(fres.scan_data.sigma_2d.begin(), fres.scan_data.sigma_2d.end());
        int nbins_mu = unique_mu.size();
        int nbins_sigma = unique_sigma.size();
        
        // Create TH2D for proper grid visualization
        TH2D *h_2d = new TH2D(Form("h_chi2_2d_%d_%d", sec.ix, sec.iy),
                              Form("Section %d,%d: #chi^{2}(#mu, #sigma)", sec.ix, sec.iy),
                              nbins_mu, mu_min_2d - 0.005, mu_max_2d + 0.005,
                              nbins_sigma, sig_min_2d - 0.001, sig_max_2d + 0.001);
        
        // Fill histogram
        for (int i = 0; i < n_2d; ++i) {
            h_2d->Fill(fres.scan_data.mu_2d[i], fres.scan_data.sigma_2d[i], fres.scan_data.chi2_2d[i]);
        }
        
        h_2d->GetXaxis()->SetTitle("#mu");
        h_2d->GetYaxis()->SetTitle("#sigma");
        h_2d->GetZaxis()->SetTitle("#chi^{2}");
        h_2d->Draw("COLZ");
        
        // Mark the best point
        TMarker *m_2d = new TMarker(fres.mu, fres.sigma, 29);
        m_2d->SetMarkerColor(kRed);
        m_2d->SetMarkerSize(2.5);
        m_2d->Draw();
    }
    
    // Plot 4: Invariant mass comparison (data vs unsmeared vs smeared sim)
    c->cd(4);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    gPad->SetRightMargin(0.05);
    
    // Clone histograms for drawing
    TH1D *h_data_draw = (TH1D*)hdata.Clone("h_data_tmp");
    TH1D *h_sim_smeared_draw = (TH1D*)hsim_smeared.Clone("h_sim_smeared_tmp");
    TH1D *h_sim_unsmeared_draw = (TH1D*)hsim_unsmeared.Clone("h_sim_unsmeared_tmp");
    
    h_data_draw->SetTitle(Form("Section %d,%d: Invariant Mass Comparison", sec.ix, sec.iy));
    h_data_draw->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    h_data_draw->GetYaxis()->SetTitle("Counts");
    // h_data_draw->GetXaxis()->SetRangeUser(0.09, 0.15);  // fixed range for clearer visualization
    h_data_draw->SetLineColor(kBlack);
    h_data_draw->SetLineWidth(2);
    h_data_draw->SetMarkerStyle(20);
    h_data_draw->SetMarkerSize(0.6);
    
    h_sim_unsmeared_draw->SetLineColor(kBlue);
    h_sim_unsmeared_draw->SetLineWidth(2);
    h_sim_unsmeared_draw->SetLineStyle(2);  // dashed
    
    h_sim_smeared_draw->SetLineColor(kRed);
    h_sim_smeared_draw->SetLineWidth(2);
    h_sim_smeared_draw->SetLineStyle(1);  // solid
    
    double max_val = max({h_data_draw->GetMaximum(), h_sim_smeared_draw->GetMaximum(), h_sim_unsmeared_draw->GetMaximum()});
    h_data_draw->SetMaximum(max_val * 1.2);
    
    h_data_draw->Draw("E");
    h_sim_unsmeared_draw->Draw("HIST SAME");
    h_sim_smeared_draw->Draw("HIST SAME");
    
    TLegend *leg = new TLegend(0.60, 0.65, 0.88, 0.88);
    leg->SetBorderSize(1);
    leg->SetFillColor(0);
    leg->SetTextSize(0.035);
    leg->AddEntry(h_data_draw, "Data", "lep");
    leg->AddEntry(h_sim_unsmeared_draw, "Unsmeared Sim", "l");
    leg->AddEntry(h_sim_smeared_draw, "Smeared Sim (E+Pos)", "l");
    leg->Draw();
    
    TLatex *t_mass = new TLatex();
    t_mass->SetNDC();
    t_mass->SetTextSize(0.032);
    t_mass->DrawLatex(0.15, 0.88, Form("#mu=%.4f, #sigma=%.4f", fres.mu, fres.sigma));
    t_mass->DrawLatex(0.15, 0.84, Form("#sigma_{pos}=%.5f cm", fres.sigma_pos));
    t_mass->DrawLatex(0.15, 0.80, Form("#chi^{2}/ndf=%.2f/%d", fres.chi2, fres.ndf));
    
    // Plot 5: Summary text
    c->cd(5);
    gPad->SetLeftMargin(0.05);
    gPad->SetRightMargin(0.05);
    
    TPaveText *pt = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
    pt->SetTextAlign(12);
    pt->SetTextFont(42);
    pt->SetTextSize(0.05);
    pt->SetFillColor(0);
    pt->SetBorderSize(1);
    
    pt->AddText(Form("Section %d, %d", sec.ix, sec.iy));
    pt->AddText(" ");
    pt->AddText(Form("x: [%.2f, %.2f] cm", sec.xlo, sec.xhi));
    pt->AddText(Form("y: [%.2f, %.2f] cm", sec.ylo, sec.yhi));
    pt->AddText(" ");
    pt->AddText("Fit Results:");
    pt->AddText(Form("#mu = %.5f", fres.mu));
    pt->AddText(Form("#sigma = %.5f", fres.sigma));
    pt->AddText(Form("#sigma_{pos} = %.5f cm", fres.sigma_pos));
    pt->AddText(" ");
    pt->AddText(Form("#chi^{2} = %.2f", fres.chi2));
    pt->AddText(Form("ndf = %d", fres.ndf));
    pt->AddText(Form("#chi^{2}/ndf = %.3f", fres.chi2_per_ndf()));
    
    pt->Draw();
    
    // Plot 6: chi^2 vs sigma_pos
    c->cd(6);
    gPad->SetLeftMargin(0.12);
    gPad->SetBottomMargin(0.12);
    
    int n_sigma_pos = fres.scan_data.sigma_pos_values.size();
    if (n_sigma_pos > 0) {
        TGraph *gr_sigma_pos = new TGraph(n_sigma_pos);
        gr_sigma_pos->SetName(Form("gr_sigma_pos_%d_%d", sec.ix, sec.iy));
        for (int i = 0; i < n_sigma_pos; ++i) {
            gr_sigma_pos->SetPoint(i, fres.scan_data.sigma_pos_values[i], fres.scan_data.chi2_vs_sigma_pos[i]);
        }
        gr_sigma_pos->SetTitle(Form("Section %d,%d: #chi^{2} vs #sigma_{pos}", sec.ix, sec.iy));
        gr_sigma_pos->GetXaxis()->SetTitle("#sigma_{pos} (position smearing) [cm]");
        gr_sigma_pos->GetYaxis()->SetTitle("#chi^{2}");
        gr_sigma_pos->SetMarkerStyle(20);
        gr_sigma_pos->SetMarkerSize(0.8);
        gr_sigma_pos->SetMarkerColor(kMagenta);
        gr_sigma_pos->SetLineColor(kMagenta);
        gr_sigma_pos->Draw("APL");
        
        // Mark the best point
        TMarker *m_sigma_pos = new TMarker(fres.sigma_pos, fres.chi2, 29);
        m_sigma_pos->SetMarkerColor(kRed);
        m_sigma_pos->SetMarkerSize(2.0);
        m_sigma_pos->Draw();
        
        TLatex *t_sigma_pos = new TLatex();
        t_sigma_pos->SetNDC();
        t_sigma_pos->SetTextSize(0.04);
        t_sigma_pos->DrawLatex(0.15, 0.85, Form("Best #sigma_{pos} = %.5f cm", fres.sigma_pos));
        t_sigma_pos->DrawLatex(0.15, 0.80, Form("#chi^{2}/ndf = %.2f/%d", fres.chi2, fres.ndf));
    }
}

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

    // Set output directory and prepend to relative paths
    string output_dir = "../output/plots/x60_4b/";
    if (out_file[0] != '/') {  // if relative path
        out_file = output_dir + out_file;
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
    // As specified: cluster_x_1, cluster_y_1, cluster_e_1, cluster_x_2, cluster_y_2, cluster_e_2, pi0_weight, scale, is_exclusive
    Double_t d_cluster_x_1, d_cluster_y_1, d_cluster_e_1;
    Double_t d_cluster_x_2, d_cluster_y_2, d_cluster_e_2;
    Double_t d_pi0_weight;
    Float_t d_scale;
    Int_t d_is_exclusive;
    
    tdata->SetBranchAddress("cluster_x_1", &d_cluster_x_1);
    tdata->SetBranchAddress("cluster_y_1", &d_cluster_y_1);
    tdata->SetBranchAddress("cluster_e_1", &d_cluster_e_1);
    tdata->SetBranchAddress("cluster_x_2", &d_cluster_x_2);
    tdata->SetBranchAddress("cluster_y_2", &d_cluster_y_2);
    tdata->SetBranchAddress("cluster_e_2", &d_cluster_e_2);
    tdata->SetBranchAddress("pi0_weight", &d_pi0_weight);
    tdata->SetBranchAddress("scale", &d_scale);
    tdata->SetBranchAddress("is_exclusive", &d_is_exclusive);

    // ===== SIMULATION TREE BRANCHES =====
    // Note: branches are Float_t arrays, not vectors
    Float_t s_clust_X[10];  // Allow up to 10 clusters
    Float_t s_clust_Y[10];
    Float_t s_clust_E[10];
    Int_t s_nclust = 0;  // Number of clusters
    Float_t s_full_weight;
    
    tsim->SetBranchAddress("clust_X", s_clust_X);
    tsim->SetBranchAddress("clust_Y", s_clust_Y);
    tsim->SetBranchAddress("clust_E", s_clust_E);
    tsim->SetBranchAddress("nclust", &s_nclust);
    tsim->SetBranchAddress("full_weight", &s_full_weight);

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

    // Prepare per-section data histograms and sim event buffers
    cout << "Scanning data tree and building per-section data histograms..." << endl;
    vector<vector<double>> data_mass_per_section(nsec);
    vector<vector<double>> data_weight_per_section(nsec);

    Long64_t ndata = tdata->GetEntries();
    for (Long64_t i=0;i<ndata;++i) {
        tdata->GetEntry(i);
        
        // Calculate invariant mass using nps_helper functions (exactly as in nps_analysis.C)
        double m = nps::invariant_mass_pi0(d_cluster_e_1, d_cluster_e_2, 
                                          d_cluster_x_1, d_cluster_x_2, 
                                          d_cluster_y_1, d_cluster_y_2, 
                                          nps::kDefaultZ_NPS_cm);
        
        // Calculate weight (pi0_weight * scale * is_exclusive)
        double weight = d_pi0_weight * d_scale * d_is_exclusive;
        
        // assign to all sections where either photon cluster lies inside
        for (int is=0; is<nsec; ++is) {
            if (inSection(sections[is], d_cluster_x_1, d_cluster_y_1) || 
                inSection(sections[is], d_cluster_x_2, d_cluster_y_2)) {
                data_mass_per_section[is].push_back(m);
                data_weight_per_section[is].push_back(weight);
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
        cerr << "No data events found in any section. Check cluster branches and calorimeter bounds." << endl; 
        return 6; 
    }
    
    // Hardcoded mass window optimized for pi0 peak region
    double mlo = 0.11;
    double mhi = 0.15;
    int nbins = 200;
    cout << "Global mass window for histograms: ["<<mlo<<","<<mhi<<"] nbins="<<nbins<<" (hardcoded)\n";

    // create TH1D per section for data and fill with weights
    vector<TH1D*> hdata_sec(nsec,nullptr);
    for (int is=0; is<nsec; ++is) {
        string hname = string("h_data_") + sections[is].name();
        hdata_sec[is] = new TH1D(hname.c_str(), hname.c_str(), nbins, mlo, mhi);
        for (size_t j=0; j<data_mass_per_section[is].size(); ++j) {
            hdata_sec[is]->Fill(data_mass_per_section[is][j], data_weight_per_section[is][j]);
        }
        cout << "Section "<<sections[is].name()<<" data entries="<<data_mass_per_section[is].size()<<"\n";
    }

    // Now load sim tree and build per-section sim event lists
    cout << "Scanning sim tree and building per-section sim event buffers..." << endl;
    vector<vector<ClusterPair>> sim_events_per_section(nsec);
    Long64_t nsim = tsim->GetEntries();
    for (Long64_t i=0;i<nsim;++i) {
        tsim->GetEntry(i);
        
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
        pair.weight = s_full_weight;
        
        for (int is=0; is<nsec; ++is) {
            if (inSection(sections[is], pair.x1, pair.y1) || 
                inSection(sections[is], pair.x2, pair.y2)) {
                sim_events_per_section[is].push_back(pair);
            }
        }
    }
    for (int is=0; is<nsec; ++is) 
        cout << "Section "<<sections[is].name()<<" sim entries="<<sim_events_per_section[is].size()<<"\n";

    // RNG - will be created per-thread in parallel region
    TRandom3 rng(0);  // Main thread RNG

    // Output file
    TFile fout(out_file.c_str(), "RECREATE");

    // CSV summary - place in same directory as ROOT output
    string csv_file = out_file;
    size_t last_slash = csv_file.find_last_of('/');
    if (last_slash != string::npos) {
        csv_file = csv_file.substr(0, last_slash + 1) + "section_map.csv";
    } else {
        csv_file = "section_map.csv";
    }
    ofstream csv(csv_file.c_str());
    csv << "ix,iy,xlo,xhi,ylo,yhi,n_data,n_sim,best_mu,best_sigma,best_sigma_pos_cm,res_A,res_B,res_C,best_chi2,ndf,chi2_ndf\n";

    // Store fit results for interpolation
    vector<FitResult> fit_results(nsec);
    vector<bool> fit_success(nsec, false);

    // Create canvas and PDF for chi^2 plots
    string pdf_file = out_file;
    size_t pdf_slash = pdf_file.find_last_of('/');
    if (pdf_slash != string::npos) {
        pdf_file = pdf_file.substr(0, pdf_slash + 1) + "chi2_scans.pdf";
    } else {
        pdf_file = "chi2_scans.pdf";
    }
    TCanvas *c_chi2 = new TCanvas("c_chi2", "Chi2 Scans", 1400, 1000);
    string pdf_open = pdf_file + "[";
    c_chi2->Print(pdf_open.c_str());
    int page_count = 0;

    // Loop over sections and fit each where we have enough stats
    // OpenMP parallelization: each section is fitted independently
    
    // ROOT is not thread-safe by default. Disable automatic histogram registration
    // to prevent race conditions when multiple threads create histograms
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
    
    #pragma omp parallel
    {
        // Each thread gets its own random number generator with unique seed
        int thread_id = omp_get_thread_num();
        TRandom3 thread_rng(thread_id * 123456 + time(NULL));
        
        #pragma omp for schedule(dynamic)
        for (int is=0; is<nsec; ++is) {
            #pragma omp critical(console)
            {
                cout << "\n==== [Thread " << thread_id << "] Fitting section "<<sections[is].name()<<" ===="<<endl;
            }
            size_t ndata_sec = data_mass_per_section[is].size();
            size_t nsim_sec = sim_events_per_section[is].size();
            if (ndata_sec < 50 || nsim_sec < 50) {
                #pragma omp critical(console)
                {
                    cout << "Skipping section "<<sections[is].name()<<" due to insufficient statistics (data="<<ndata_sec<<" sim="<<nsim_sec<<")"<<endl;
                }
                #pragma omp critical(file_io)
                {
                    // still write empty histograms for bookkeeping
                    TDirectory *dir = fout.mkdir(sections[is].name().c_str()); dir->cd();
                    hdata_sec[is]->Write(hdata_sec[is]->GetName());
                    fout.cd();
                    csv << sections[is].ix<<","<<sections[is].iy<<","<<sections[is].xlo<<","<<sections[is].xhi<<","<<sections[is].ylo<<","<<sections[is].yhi<<","<<ndata_sec<<","<<nsim_sec<<",,,,,,,,"<<"\n";
                }
                continue;
            }

            // Fit using thread-local RNG
            FitResult fres = fit_section(sim_events_per_section[is], *hdata_sec[is], thread_rng, Nsmear);
            #pragma omp critical(console)
            {
                cout << "Section "<<sections[is].name()<<" -> mu="<<fres.mu<<" sigma="<<fres.sigma<<" sigma_pos="<<std::fixed<<std::setprecision(5)<<fres.sigma_pos<<" cm chi2="<<fres.chi2<<" ndf="<<fres.ndf<<" chi2/ndf="<<fres.chi2_per_ndf()<<"\n";
            }
            
            // Store fit result
            fit_results[is] = fres;
            fit_success[is] = true;

            // Build unsmeared sim histogram for comparison
            TH1D hunsmeared((string("h_unsmeared_")+sections[is].name()).c_str(), "h_unsmeared", nbins, mlo, mhi);
            for (const auto &ev : sim_events_per_section[is]) {
                Vec4 ph1 = nps::photon4vector(ev.e1, ev.x1, ev.y1, nps::kDefaultZ_NPS_cm);
                Vec4 ph2 = nps::photon4vector(ev.e2, ev.x2, ev.y2, nps::kDefaultZ_NPS_cm);
                
                TLorentzVector g1(ph1[1], ph1[2], ph1[3], ph1[0]);
                TLorentzVector g2(ph2[1], ph2[2], ph2[3], ph2[0]);
                TLorentzVector pi0 = g1 + g2;
                
                // hunsmeared.Fill(pi0.M(), ev.weight);  // COMMENTED OUT: simulation weighting disabled
                hunsmeared.Fill(pi0.M(), 1.0);
            }
            if (hunsmeared.Integral()>0 && hdata_sec[is]->Integral()>0) 
                hunsmeared.Scale(hdata_sec[is]->Integral() / hunsmeared.Integral());

            // Build final smeared sim histogram at best-fit for this section (with energy and position smearing)
            TH1D hfinal((string("h_smeared_best_")+sections[is].name()).c_str(), "h_smeared_best", nbins, mlo, mhi);
            
            // Pre-allocate random arrays outside event loop
            vector<double> rand_energy(Nsmear * 2);
            vector<double> rand_pos(Nsmear * 4);
            const double z_nps = nps::kDefaultZ_NPS_cm;
            const double z_sq = z_nps * z_nps;
            
            for (const auto &ev : sim_events_per_section[is]) {
                double E1 = ev.e1, E2 = ev.e2; 
                if (E1<=0||E2<=0) continue;
            
                // Pre-scale energies by mu
                double E1_scaled = fres.mu * E1;
                double E2_scaled = fres.mu * E2;
                
                // Compute energy-dependent resolutions based on selected model
                double sigma_E1, sigma_E2;
                if (USE_SIMPLE_STOCHASTIC_MODEL) {
                    // Simple stochastic: σ_E = σ × √E
                    sigma_E1 = fres.sigma * sqrt(E1);
                    sigma_E2 = fres.sigma * sqrt(E2);
                } else {
                    // Full 3-term: σ_E = σ × E × √((A/√E)² + (B/E)² + C²)
                    double A_sq = fres.res_A * fres.res_A;
                    double B_sq = fres.res_B * fres.res_B;
                    double C_sq = fres.res_C * fres.res_C;
                    sigma_E1 = fres.sigma * E1 * sqrt(A_sq / E1 + B_sq / (E1*E1) + C_sq);
                    sigma_E2 = fres.sigma * E2 * sqrt(A_sq / E2 + B_sq / (E2*E2) + C_sq);
                }
                
                // Batch generate random energies using thread-local RNG
                for (int k=0; k<Nsmear; ++k) {
                    rand_energy[k*2] = thread_rng.Gaus(E1_scaled, sigma_E1);
                    rand_energy[k*2+1] = thread_rng.Gaus(E2_scaled, sigma_E2);
                    if (rand_energy[k*2]<=0) rand_energy[k*2] = NONPOSITIVE_CLAMP;
                    if (rand_energy[k*2+1]<=0) rand_energy[k*2+1] = NONPOSITIVE_CLAMP;
                }
                if (fres.sigma_pos > 0) {
                    for (int k=0; k<Nsmear*4; ++k) {
                        rand_pos[k] = thread_rng.Gaus(0, fres.sigma_pos);
                    }
                }
                
                for (int k=0;k<Nsmear;++k) {
                    double E1_smeared = rand_energy[k*2]; 
                    double E2_smeared = rand_energy[k*2+1];
                    
                    double x1_smeared, y1_smeared, x2_smeared, y2_smeared;
                    if (fres.sigma_pos > 0) {
                        x1_smeared = ev.x1 + rand_pos[k*4];
                        y1_smeared = ev.y1 + rand_pos[k*4+1];
                        x2_smeared = ev.x2 + rand_pos[k*4+2];
                        y2_smeared = ev.y2 + rand_pos[k*4+3];
                    } else {
                        x1_smeared = ev.x1; y1_smeared = ev.y1;
                        x2_smeared = ev.x2; y2_smeared = ev.y2;
                    }
                    
                    // Direct mass calculation
                    double r1_sq = x1_smeared*x1_smeared + y1_smeared*y1_smeared + z_sq;
                    double r1_inv = 1.0 / sqrt(r1_sq);
                    double px1 = E1_smeared * x1_smeared * r1_inv;
                    double py1 = E1_smeared * y1_smeared * r1_inv;
                    double pz1 = E1_smeared * z_nps * r1_inv;
                    
                    double r2_sq = x2_smeared*x2_smeared + y2_smeared*y2_smeared + z_sq;
                    double r2_inv = 1.0 / sqrt(r2_sq);
                    double px2 = E2_smeared * x2_smeared * r2_inv;
                    double py2 = E2_smeared * y2_smeared * r2_inv;
                    double pz2 = E2_smeared * z_nps * r2_inv;
                    
                    double E_tot = E1_smeared + E2_smeared;
                    double px_tot = px1 + px2;
                    double py_tot = py1 + py2;
                    double pz_tot = pz1 + pz2;
                    double m_sq = E_tot*E_tot - px_tot*px_tot - py_tot*py_tot - pz_tot*pz_tot;
                    double mass = (m_sq > 0) ? sqrt(m_sq) : 0.0;
                    
                    hfinal.Fill(mass, 1.0);
                }
            }
            if (hfinal.Integral()>0 && hdata_sec[is]->Integral()>0) 
                hfinal.Scale(hdata_sec[is]->Integral() / hfinal.Integral());

            // Thread-safe file I/O: only one thread writes at a time
            #pragma omp critical(file_io)
            {
                // Create chi^2 plots for this section
                plot_chi2_scans(sections[is], fres, c_chi2, page_count++, *hdata_sec[is], hfinal, hunsmeared);
                c_chi2->Print(pdf_file.c_str());

                // write per-section outputs
                TDirectory *dir = fout.mkdir(sections[is].name().c_str()); dir->cd();
                hdata_sec[is]->Write(hdata_sec[is]->GetName());
                hunsmeared.Write(hunsmeared.GetName());
                hfinal.Write(hfinal.GetName());
                TNamed meta_mu((string("best_mu_")+sections[is].name()).c_str(), to_string(fres.mu).c_str()); 
                meta_mu.Write();
                TNamed meta_sigma((string("best_sigma_")+sections[is].name()).c_str(), to_string(fres.sigma).c_str()); 
                meta_sigma.Write();
                TNamed meta_sigma_pos((string("best_sigma_pos_")+sections[is].name()).c_str(), to_string(fres.sigma_pos).c_str()); 
                meta_sigma_pos.Write();
                TNamed meta_chi((string("best_chi2_")+sections[is].name()).c_str(), to_string(fres.chi2).c_str()); 
                meta_chi.Write();
                TNamed meta_ndf((string("ndf_")+sections[is].name()).c_str(), to_string(fres.ndf).c_str()); 
                meta_ndf.Write();
                TNamed meta_chi2ndf((string("chi2_ndf_")+sections[is].name()).c_str(), to_string(fres.chi2_per_ndf()).c_str()); 
                meta_chi2ndf.Write();
                fout.cd();

                csv << sections[is].ix<<","<<sections[is].iy<<","<<sections[is].xlo<<","<<sections[is].xhi<<","<<sections[is].ylo<<","<<sections[is].yhi<<","<<ndata_sec<<","<<nsim_sec<<","<<fres.mu<<","<<fres.sigma<<","<<fres.sigma_pos<<","<<fres.res_A<<","<<fres.res_B<<","<<fres.res_C<<","<<fres.chi2<<","<<fres.ndf<<","<<fres.chi2_per_ndf()<<"\n";
            }
        }  // end parallel for
    }  // end parallel region
    
    // Restore original ROOT histogram directory behavior
    TH1::AddDirectory(original_adddir_status);

    // ============================================================================
    // Additional summary plots
    // ============================================================================
    
    // Plot 1: Grid of pi0 invariant mass for each section with chi2 values
    cout << "\nCreating grid plot of section invariant masses..." << endl;
    
    // Find global maximum for uniform y-axis scaling
    double global_max = 0.0;
    for (int is = 0; is < nsec; ++is) {
        if (hdata_sec[is]) {
            double hmax = hdata_sec[is]->GetMaximum();
            if (hmax > global_max) global_max = hmax;
        }
    }
    global_max *= 1.15;  // Add 15% headroom for text
    
    TCanvas *c_grid = new TCanvas("c_grid", "Section Grid", 1600, 1400);
    c_grid->Divide(nx, ny, 0.001, 0.001);
    
    for (int iy = ny-1; iy >= 0; --iy) {  // top to bottom
        for (int ix = 0; ix < nx; ++ix) {   // left to right
            int sec_idx = iy * nx + ix;
            int pad_idx = (ny-1-iy) * nx + ix + 1;  // ROOT pad numbering: top-left is 1
            
            c_grid->cd(pad_idx);
            gPad->SetLeftMargin(0.12);
            gPad->SetRightMargin(0.05);
            gPad->SetTopMargin(0.12);
            gPad->SetBottomMargin(0.12);
            
            if (sec_idx < nsec && hdata_sec[sec_idx]) {
                TH1D *h_copy = (TH1D*)hdata_sec[sec_idx]->Clone();
                h_copy->SetTitle(Form("Sec(%d,%d)", sections[sec_idx].ix, sections[sec_idx].iy));
                h_copy->GetXaxis()->SetTitle("M_{#gamma#gamma}");
                h_copy->GetYaxis()->SetTitle("Counts");
                h_copy->GetXaxis()->SetTitleSize(0.06);
                h_copy->GetYaxis()->SetTitleSize(0.06);
                h_copy->GetXaxis()->SetLabelSize(0.05);
                h_copy->GetYaxis()->SetLabelSize(0.05);
                h_copy->SetLineColor(kBlack);
                h_copy->SetLineWidth(2);
                h_copy->SetMaximum(global_max);  // Set uniform y-axis
                h_copy->Draw("HIST");
                
                // Add chi2 text if available
                if (fit_success[sec_idx]) {
                    TLatex *t = new TLatex();
                    t->SetNDC();
                    t->SetTextSize(0.06);
                    t->SetTextColor(kRed);
                    t->DrawLatex(0.15, 0.85, Form("#chi^{2}/ndf=%.1f", fit_results[sec_idx].chi2_per_ndf()));
                }
            }
        }
    }
    
    c_grid->Print(pdf_file.c_str());
    delete c_grid;
    
    // Plot 2: Combined invariant mass for entire calorimeter
    cout << "Creating combined calorimeter plot..." << endl;
    TCanvas *c_combined = new TCanvas("c_combined", "Combined Calorimeter", 1200, 800);
    c_combined->Divide(1, 1);
    c_combined->cd(1);
    gPad->SetLeftMargin(0.12);
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    gPad->SetBottomMargin(0.12);
    
    // Combine all data
    TH1D *h_data_combined = new TH1D("h_data_combined", "Combined Calorimeter: Data vs Smeared Sim", nbins, mlo, mhi);
    TH1D *h_sim_combined = new TH1D("h_sim_combined", "h_sim_combined", nbins, mlo, mhi);
    TH1D *h_unsmeared_combined = new TH1D("h_unsmeared_combined", "h_unsmeared_combined", nbins, mlo, mhi);
    
    // Fill combined histograms - use only successfully fitted sections
    for (int is = 0; is < nsec; ++is) {
        if (!fit_success[is]) continue;
        
        // Add data histogram
        h_data_combined->Add(hdata_sec[is]);
        
        // Build and add unsmeared sim
        for (const auto &ev : sim_events_per_section[is]) {
            Vec4 ph1 = nps::photon4vector(ev.e1, ev.x1, ev.y1, nps::kDefaultZ_NPS_cm);
            Vec4 ph2 = nps::photon4vector(ev.e2, ev.x2, ev.y2, nps::kDefaultZ_NPS_cm);
            TLorentzVector g1(ph1[1], ph1[2], ph1[3], ph1[0]);
            TLorentzVector g2(ph2[1], ph2[2], ph2[3], ph2[0]);
            TLorentzVector pi0 = g1 + g2;
            // h_unsmeared_combined->Fill(pi0.M(), ev.weight);  // COMMENTED OUT: simulation weighting disabled
            h_unsmeared_combined->Fill(pi0.M(), 1.0);
        }
        
        // Build and add smeared sim (with energy and position smearing)
        // Pre-allocate random arrays outside event loop  
        vector<double> rand_energy(Nsmear * 2);
        vector<double> rand_pos(Nsmear * 4);
        const double z_nps = nps::kDefaultZ_NPS_cm;
        const double z_sq = z_nps * z_nps;
        
        for (const auto &ev : sim_events_per_section[is]) {
            double E1 = ev.e1, E2 = ev.e2;
            if (E1 <= 0 || E2 <= 0) continue;
            
            // Batch generate random numbers
            for (int k=0; k<Nsmear*2; ++k) {
                rand_energy[k] = rng.Gaus(fit_results[is].mu, fit_results[is].sigma);
                if (rand_energy[k]<=0) rand_energy[k] = NONPOSITIVE_CLAMP;
            }
            if (fit_results[is].sigma_pos > 0) {
                for (int k=0; k<Nsmear*4; ++k) {
                    rand_pos[k] = rng.Gaus(0, fit_results[is].sigma_pos);
                }
            }
            
            for (int k = 0; k < Nsmear; ++k) {
                double E1_smeared = E1 * rand_energy[k*2];
                double E2_smeared = E2 * rand_energy[k*2+1];
                
                double x1_smeared, y1_smeared, x2_smeared, y2_smeared;
                if (fit_results[is].sigma_pos > 0) {
                    x1_smeared = ev.x1 + rand_pos[k*4];
                    y1_smeared = ev.y1 + rand_pos[k*4+1];
                    x2_smeared = ev.x2 + rand_pos[k*4+2];
                    y2_smeared = ev.y2 + rand_pos[k*4+3];
                } else {
                    x1_smeared = ev.x1; y1_smeared = ev.y1;
                    x2_smeared = ev.x2; y2_smeared = ev.y2;
                }
                
                // Direct mass calculation
                double r1_sq = x1_smeared*x1_smeared + y1_smeared*y1_smeared + z_sq;
                double r1_inv = 1.0 / sqrt(r1_sq);
                double px1 = E1_smeared * x1_smeared * r1_inv;
                double py1 = E1_smeared * y1_smeared * r1_inv;
                double pz1 = E1_smeared * z_nps * r1_inv;
                
                double r2_sq = x2_smeared*x2_smeared + y2_smeared*y2_smeared + z_sq;
                double r2_inv = 1.0 / sqrt(r2_sq);
                double px2 = E2_smeared * x2_smeared * r2_inv;
                double py2 = E2_smeared * y2_smeared * r2_inv;
                double pz2 = E2_smeared * z_nps * r2_inv;
                
                double E_tot = E1_smeared + E2_smeared;
                double px_tot = px1 + px2;
                double py_tot = py1 + py2;
                double pz_tot = pz1 + pz2;
                double m_sq = E_tot*E_tot - px_tot*px_tot - py_tot*py_tot - pz_tot*pz_tot;
                double mass = (m_sq > 0) ? sqrt(m_sq) : 0.0;
                
                h_sim_combined->Fill(mass, 1.0);
            }
        }
    }
    
    // Scale simulation to match data
    if (h_unsmeared_combined->Integral() > 0 && h_data_combined->Integral() > 0)
        h_unsmeared_combined->Scale(h_data_combined->Integral() / h_unsmeared_combined->Integral());
    if (h_sim_combined->Integral() > 0 && h_data_combined->Integral() > 0)
        h_sim_combined->Scale(h_data_combined->Integral() / h_sim_combined->Integral());
    
    // Style histograms
    h_data_combined->SetLineColor(kBlack);
    h_data_combined->SetLineWidth(2);
    h_data_combined->SetMarkerStyle(20);
    h_data_combined->SetMarkerSize(0.6);
    h_data_combined->GetXaxis()->SetTitle("M_{#gamma#gamma} [GeV/c^{2}]");
    h_data_combined->GetYaxis()->SetTitle("Counts");
    
    h_unsmeared_combined->SetLineColor(kBlue);
    h_unsmeared_combined->SetLineWidth(2);
    h_unsmeared_combined->SetLineStyle(2);
    
    h_sim_combined->SetLineColor(kRed);
    h_sim_combined->SetLineWidth(2);
    h_sim_combined->SetLineStyle(1);
    
    double max_combined = max({h_data_combined->GetMaximum(), 
                               h_sim_combined->GetMaximum(), 
                               h_unsmeared_combined->GetMaximum()});
    h_data_combined->SetMaximum(max_combined * 1.2);
    
    h_data_combined->Draw("E");
    h_unsmeared_combined->Draw("HIST SAME");
    h_sim_combined->Draw("HIST SAME");
    
    TLegend *leg_combined = new TLegend(0.60, 0.65, 0.88, 0.88);
    leg_combined->SetBorderSize(1);
    leg_combined->SetFillColor(0);
    leg_combined->SetTextSize(0.035);
    leg_combined->AddEntry(h_data_combined, "Data (All Sections)", "lep");
    leg_combined->AddEntry(h_unsmeared_combined, "Unsmeared Sim", "l");
    leg_combined->AddEntry(h_sim_combined, "Smeared Sim (Optimized)", "l");
    leg_combined->Draw();
    
    c_combined->Print(pdf_file.c_str());
    delete c_combined;
    delete h_data_combined;
    delete h_sim_combined;
    delete h_unsmeared_combined;

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
                           fit_results[is].sigma);
        }
    }
    
    // Save 2D interpolated maps for visualization
    string interp_file = out_file;
    size_t dot_pos = interp_file.find_last_of('.');
    if (dot_pos != string::npos) {
        interp_file.insert(dot_pos, "_interpolated");
    } else {
        interp_file += "_interpolated.root";
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
    cout << "double mu, sigma;\n";
    cout << "calMap.getInterpolatedParams(cluster_x, cluster_y, mu, sigma);\n";
    cout << "double smeared_energy = original_energy * rng.Gaus(mu, sigma);\n";
    
    return 0;
}
