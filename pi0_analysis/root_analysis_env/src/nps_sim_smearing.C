// File: minimize_chi2_per_section.cpp
// Purpose: Perform per-section (overlapping grid) smearing and chi2 minimization
//          for calorimeter MC tuning (mu, sigma, theta_s) exactly as described.
//
// Compile:
//   g++ minimize_chi2_per_section.cpp `root-config --cflags --libs` -O2 -std=c++17 -o minimize_chi2_per_section
//
// Usage example:
//   ./minimize_chi2_per_section data.root dataTree sim.root simTree out.root 8 8 -20 20 -20 20 0.2 20
//
// Arguments:
// 1) data_file      - ROOT file with experimental events (TTree)
// 2) data_tree      - name of TTree inside data_file (events must contain per-photon cluster X,Y and 4-vector branches)
// 3) sim_file       - ROOT file with simulated events (TTree)
// 4) sim_tree       - name of TTree inside sim_file
// 5) out_file       - output ROOT file to write per-section histograms and best-fit metadata
// 6) nx             - number of sections in x (columns)
// 7) ny             - number of sections in y (rows)
// 8) x_min          - calorimeter x lower bound (same units as cluster positions)
// 9) x_max          - calorimeter x upper bound
//10) y_min          - calorimeter y lower bound
//11) y_max          - calorimeter y upper bound
//12) overlap_frac   - fraction of section width used for overlap (e.g. 0.2 means 20% overlap)
//13) Nsmear         - number of smearing replicas per sim event (default 20)
//
// Required branches (edit the names below if your tree uses different names):
// For data tree: (per-event)
//    phot1_px, phot1_py, phot1_pz, phot1_E,
//    phot2_px, phot2_py, phot2_pz, phot2_E,
//    phot1_clust_x, phot1_clust_y, phot2_clust_x, phot2_clust_y
// For sim tree: same branch names but from simulation (we need cluster positions too)
//
// Note: The program assigns an event to any section for which AT LEAST ONE photon cluster
// lies inside that section's (overlapping) area. In other words, events can contribute to
// multiple overlapping sections — this matches the overlapping mapping strategy.
//
// Output: out_file contains per-section histograms named:
//    h_data_sec_<ix>_<iy>
//    h_sim_sec_<ix>_<iy>
//    h_smeared_best_sec_<ix>_<iy>
// And TNamed metadata for each section: best_mu_sec_<ix>_<iy>, best_sigma_sec_<ix>_<iy>, best_theta_deg_sec_<ix>_<iy>, best_chi2_sec_<ix>_<iy>
// Also a CSV named "section_map.csv" is saved to disk (in running directory) with a table of coefficients.

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include <limits>
#include <algorithm>
#include <sstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom3.h"
#include "TLorentzVector.h"
#include "TNamed.h"
#include "TDirectory.h"

using namespace std;

// ----- CONFIGURABLE BRANCH NAMES (edit if your trees use different names) -----
const string BR_ph1_px = "phot1_px";
const string BR_ph1_py = "phot1_py";
const string BR_ph1_pz = "phot1_pz";
const string BR_ph1_E  = "phot1_E";
const string BR_ph2_px = "phot2_px";
const string BR_ph2_py = "phot2_py";
const string BR_ph2_pz = "phot2_pz";
const string BR_ph2_E  = "phot2_E";
const string BR_ph1_cx = "phot1_clust_x"; // cluster X position (calorimeter coords)
const string BR_ph1_cy = "phot1_clust_y";
const string BR_ph2_cx = "phot2_clust_x";
const string BR_ph2_cy = "phot2_clust_y";

// ----- GLOBAL SETTINGS -----
const double NONPOSITIVE_CLAMP = 1e-6; // clamp for non-positive smearing draws

// rotate vector v around axis a_unit by angle theta (radians)
TVector3 rotateAroundAxis(const TVector3 &v, const TVector3 &a_unit, double theta) {
    TVector3 term1 = v * cos(theta);
    TVector3 term2 = a_unit.Cross(v) * sin(theta);
    TVector3 term3 = a_unit * (a_unit.Dot(v)) * (1.0 - cos(theta));
    return term1 + term2 + term3;
}

// compute chi2 using your formula: sum ( (s-d)^2 / (s+d) ) skipping bins with denom<=0
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

// Evaluate chi2 for (mu, sigma, theta_s) for a given set of sim events that belong to a section
// - simEvents: vector of pairs (g1,g2) already filtered to events that belong to this section (but may include events that also belong elsewhere)
// - hdata: data histogram for that section (must already be filled)
// - theta_s is in radians
// - Nsmear: replicas per event

double eval_chi2_section(double mu, double sigma, double theta_s,
                          const vector<pair<TLorentzVector,TLorentzVector>> &simEvents,
                          const TH1D &hdata,
                          TRandom3 &rng,
                          int Nsmear) {
    TH1D hsim("hsim_tmp","smeared sim (tmp)", hdata.GetNbinsX(), hdata.GetXaxis()->GetXmin(), hdata.GetXaxis()->GetXmax());

    for (const auto &ev : simEvents) {
        TLorentzVector g1 = ev.first;
        TLorentzVector g2 = ev.second;
        TLorentzVector pi0 = g1 + g2;
        TVector3 beta = pi0.BoostVector();
        // boost to rest
        TLorentzVector g1_rest = g1; g1_rest.Boost(-beta);
        TLorentzVector g2_rest = g2; g2_rest.Boost(-beta);
        TVector3 v1 = g1_rest.Vect();
        TVector3 v2 = g2_rest.Vect();
        TVector3 axis = v1.Cross(v2);
        if (axis.Mag() >= 1e-12) {
            TVector3 axis_u = axis.Unit();
            TVector3 v1rot = rotateAroundAxis(v1, axis_u, +0.5*theta_s);
            TVector3 v2rot = rotateAroundAxis(v2, axis_u, -0.5*theta_s);
            g1_rest.SetVectM(v1rot.Unit() * v1.Mag(), 0.0);
            g2_rest.SetVectM(v2rot.Unit() * v2.Mag(), 0.0);
            g1_rest.SetE(g1_rest.E());
            g2_rest.SetE(g2_rest.E());
        }
        TLorentzVector g1_lab_rot = g1_rest; g1_lab_rot.Boost(beta);
        TLorentzVector g2_lab_rot = g2_rest; g2_lab_rot.Boost(beta);
        double E1 = g1_lab_rot.E();
        double E2 = g2_lab_rot.E();
        if (E1 <= 0 || E2 <= 0) continue;
        for (int k=0;k<Nsmear;++k) {
            double draw1 = rng.Gaus(mu, sigma);
            double draw2 = rng.Gaus(mu, sigma);
            if (draw1 <= 0) draw1 = NONPOSITIVE_CLAMP;
            if (draw2 <= 0) draw2 = NONPOSITIVE_CLAMP;
            double E1p = E1 * draw1;
            double E2p = E2 * draw2;
            double s1 = E1p / E1;
            double s2 = E2p / E2;
            TLorentzVector g1p(g1_lab_rot.Px()*s1, g1_lab_rot.Py()*s1, g1_lab_rot.Pz()*s1, E1p);
            TLorentzVector g2p(g2_lab_rot.Px()*s2, g2_lab_rot.Py()*s2, g2_lab_rot.Pz()*s2, E2p);
            TLorentzVector pi0s = g1p + g2p;
            hsim.Fill(pi0s.M());
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

// Minimization strategy per section -> returns best (mu,sigma,theta,chi2)
struct FitResult { double mu, sigma, theta_rad, chi2; };

FitResult fit_section(const vector<pair<TLorentzVector,TLorentzVector>> &simEvents,
                      const TH1D &hdata,
                      TRandom3 &rng,
                      int Nsmear) {
    FitResult res; res.mu = 1.0; res.sigma = 0.05; res.theta_rad = 0.0; res.chi2 = 1e300;
    // 1) mu scan (sigma fixed = 0.05)
    double sigma_fix_for_mu = 0.05;
    double mu_min = 0.8, mu_max = 1.2, mu_step = 0.04;
    double best_mu = mu_min; double best_chi2 = 1e300;
    for (double mu = mu_min; mu <= mu_max + 1e-12; mu += mu_step) {
        double chi2 = eval_chi2_section(mu, sigma_fix_for_mu, 0.0, simEvents, hdata, rng, Nsmear);
        if (chi2 < best_chi2) { best_chi2 = chi2; best_mu = mu; }
    }
    double step = mu_step;
    while (step >= 0.001) {
        double start = max(0.0, best_mu - 2*step);
        double end = best_mu + 2*step;
        double new_best_mu = best_mu; double new_best_chi2 = best_chi2;
        for (double mu = start; mu <= end + 1e-15; mu += step/2.5) {
            double chi2 = eval_chi2_section(mu, sigma_fix_for_mu, 0.0, simEvents, hdata, rng, Nsmear);
            if (chi2 < new_best_chi2) { new_best_chi2 = chi2; new_best_mu = mu; }
        }
        best_mu = new_best_mu; best_chi2 = new_best_chi2; step = step / 2.5;
    }
    double mu_star = best_mu;

    // 2) sigma scan (mu fixed)
    double mu_fix_for_sigma = mu_star;
    double sigma_min = 0.03, sigma_max = 0.08, sigma_step = 0.005;
    double best_sigma = sigma_min; best_chi2 = 1e300;
    for (double sigma = sigma_min; sigma <= sigma_max + 1e-12; sigma += sigma_step) {
        double chi2 = eval_chi2_section(mu_fix_for_sigma, sigma, 0.0, simEvents, hdata, rng, Nsmear);
        if (chi2 < best_chi2) { best_chi2 = chi2; best_sigma = sigma; }
    }
    step = sigma_step;
    while (step >= 0.0001) {
        double start = max(0.0, best_sigma - 2*step);
        double end = best_sigma + 2*step;
        double new_best_sigma = best_sigma; double new_best_chi2 = best_chi2;
        for (double sigma = start; sigma <= end + 1e-15; sigma += step/2.5) {
            double chi2 = eval_chi2_section(mu_fix_for_sigma, sigma, 0.0, simEvents, hdata, rng, Nsmear);
            if (chi2 < new_best_chi2) { new_best_chi2 = chi2; new_best_sigma = sigma; }
        }
        best_sigma = new_best_sigma; best_chi2 = new_best_chi2; step = step / 2.5;
    }
    double sigma_star = best_sigma;

    // 3) theta scan (mu*,sigma* fixed)
    double theta_deg_min = -5.0, theta_deg_max = 5.0; // degrees
    double theta_step_deg = 1.0;
    double best_theta_deg = 0.0; best_chi2 = 1e300;
    for (double thdeg = theta_deg_min; thdeg <= theta_deg_max + 1e-12; thdeg += theta_step_deg) {
        double thrad = thdeg * M_PI / 180.0;
        double chi2 = eval_chi2_section(mu_star, sigma_star, thrad, simEvents, hdata, rng, Nsmear);
        if (chi2 < best_chi2) { best_chi2 = chi2; best_theta_deg = thdeg; }
    }
    double step_deg = theta_step_deg;
    while ((step_deg/180.0*M_PI) >= (0.0005)) {
        double start_deg = best_theta_deg - 2*step_deg;
        double end_deg = best_theta_deg + 2*step_deg;
        double new_best_theta_deg = best_theta_deg; double new_best_chi2 = best_chi2;
        for (double thdeg = start_deg; thdeg <= end_deg + 1e-12; thdeg += step_deg/2.5) {
            double thrad = thdeg * M_PI / 180.0;
            double chi2 = eval_chi2_section(mu_star, sigma_star, thrad, simEvents, hdata, rng, Nsmear);
            if (chi2 < new_best_chi2) { new_best_chi2 = chi2; new_best_theta_deg = thdeg; }
        }
        best_theta_deg = new_best_theta_deg; best_chi2 = new_best_chi2; step_deg = step_deg / 2.5;
    }
    double theta_star = best_theta_deg * M_PI / 180.0;

    // 4) final 3D refine around (mu_star, sigma_star, theta_star)
    double mu0 = mu_star, sigma0 = sigma_star, theta0 = theta_star;
    double twoD_step_mu = 0.0025, twoD_step_sigma = 0.0025, twoD_step_theta = 0.5 * M_PI / 180.0; // 0.5 deg
    best_chi2 = eval_chi2_section(mu0, sigma0, theta0, simEvents, hdata, rng, Nsmear);
    while (twoD_step_mu >= 0.0001 || twoD_step_sigma >= 0.0001 || twoD_step_theta >= (0.0005)) {
        double grid_best_mu = mu0, grid_best_sigma = sigma0, grid_best_theta = theta0;
        double grid_best_chi2 = best_chi2;
        for (double mu = mu0 - 2*twoD_step_mu; mu <= mu0 + 2*twoD_step_mu + 1e-15; mu += twoD_step_mu) {
            for (double sig = max(1e-6, sigma0 - 2*twoD_step_sigma); sig <= sigma0 + 2*twoD_step_sigma + 1e-15; sig += twoD_step_sigma) {
                for (double th = theta0 - 2*twoD_step_theta; th <= theta0 + 2*twoD_step_theta + 1e-15; th += twoD_step_theta) {
                    if (sig <= 0) continue;
                    double chi2 = eval_chi2_section(mu, sig, th, simEvents, hdata, rng, Nsmear);
                    if (chi2 < grid_best_chi2) {
                        grid_best_chi2 = chi2; grid_best_mu = mu; grid_best_sigma = sig; grid_best_theta = th;
                    }
                }
            }
        }
        mu0 = grid_best_mu; sigma0 = grid_best_sigma; theta0 = grid_best_theta; best_chi2 = grid_best_chi2;
        twoD_step_mu /= 2.0; twoD_step_sigma /= 2.0; twoD_step_theta /= 2.0;
        if (twoD_step_mu < 1e-6 && twoD_step_sigma < 1e-6 && twoD_step_theta < 1e-6) break;
    }

    res.mu = mu0; res.sigma = sigma0; res.theta_rad = theta0; res.chi2 = best_chi2;
    return res;
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

    // open data and sim files
    TFile fdata(data_file.c_str(), "READ");
    if (fdata.IsZombie()) { cerr << "Cannot open data file " << data_file << endl; return 2; }
    TTree *tdata = dynamic_cast<TTree*>(fdata.Get(data_tree_name.c_str()));
    if (!tdata) { cerr << "Cannot find data tree '" << data_tree_name << "' in " << data_file << endl; return 3; }

    TFile fsim(sim_file.c_str(), "READ");
    if (fsim.IsZombie()) { cerr << "Cannot open sim file " << sim_file << endl; return 4; }
    TTree *tsim = dynamic_cast<TTree*>(fsim.Get(sim_tree_name.c_str()));
    if (!tsim) { cerr << "Cannot find sim tree '" << sim_tree_name << "' in " << sim_file << endl; return 5; }

    // prepare branch variables for data tree
    float d_ph1_px, d_ph1_py, d_ph1_pz, d_ph1_E;
    float d_ph2_px, d_ph2_py, d_ph2_pz, d_ph2_E;
    float d_ph1_cx, d_ph1_cy, d_ph2_cx, d_ph2_cy;
    tdata->SetBranchAddress(BR_ph1_px.c_str(), &d_ph1_px);
    tdata->SetBranchAddress(BR_ph1_py.c_str(), &d_ph1_py);
    tdata->SetBranchAddress(BR_ph1_pz.c_str(), &d_ph1_pz);
    tdata->SetBranchAddress(BR_ph1_E.c_str(),  &d_ph1_E);
    tdata->SetBranchAddress(BR_ph2_px.c_str(), &d_ph2_px);
    tdata->SetBranchAddress(BR_ph2_py.c_str(), &d_ph2_py);
    tdata->SetBranchAddress(BR_ph2_pz.c_str(), &d_ph2_pz);
    tdata->SetBranchAddress(BR_ph2_E.c_str(),  &d_ph2_E);
    tdata->SetBranchAddress(BR_ph1_cx.c_str(), &d_ph1_cx);
    tdata->SetBranchAddress(BR_ph1_cy.c_str(), &d_ph1_cy);
    tdata->SetBranchAddress(BR_ph2_cx.c_str(), &d_ph2_cx);
    tdata->SetBranchAddress(BR_ph2_cy.c_str(), &d_ph2_cy);

    // prepare branch variables for sim tree
    float s_ph1_px, s_ph1_py, s_ph1_pz, s_ph1_E;
    float s_ph2_px, s_ph2_py, s_ph2_pz, s_ph2_E;
    float s_ph1_cx, s_ph1_cy, s_ph2_cx, s_ph2_cy;
    tsim->SetBranchAddress(BR_ph1_px.c_str(), &s_ph1_px);
    tsim->SetBranchAddress(BR_ph1_py.c_str(), &s_ph1_py);
    tsim->SetBranchAddress(BR_ph1_pz.c_str(), &s_ph1_pz);
    tsim->SetBranchAddress(BR_ph1_E.c_str(),  &s_ph1_E);
    tsim->SetBranchAddress(BR_ph2_px.c_str(), &s_ph2_px);
    tsim->SetBranchAddress(BR_ph2_py.c_str(), &s_ph2_py);
    tsim->SetBranchAddress(BR_ph2_pz.c_str(), &s_ph2_pz);
    tsim->SetBranchAddress(BR_ph2_E.c_str(),  &s_ph2_E);
    tsim->SetBranchAddress(BR_ph1_cx.c_str(), &s_ph1_cx);
    tsim->SetBranchAddress(BR_ph1_cy.c_str(), &s_ph1_cy);
    tsim->SetBranchAddress(BR_ph2_cx.c_str(), &s_ph2_cx);
    tsim->SetBranchAddress(BR_ph2_cy.c_str(), &s_ph2_cy);

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
    // Choose histogram binning centrally: we'll use global mass window from data tree's M distribution
    // First, scan data tree once to determine mass window and fill temporary vector of data events by section
    cout << "Scanning data tree and building per-section data histograms..." << endl;
    vector<vector<double>> data_mass_per_section(nsec); // store masses for filling hist later

    Long64_t ndata = tdata->GetEntries();
    for (Long64_t i=0;i<ndata;++i) {
        tdata->GetEntry(i);
        TLorentzVector g1(d_ph1_px, d_ph1_py, d_ph1_pz, d_ph1_E);
        TLorentzVector g2(d_ph2_px, d_ph2_py, d_ph2_pz, d_ph2_E);
        double m = (g1+g2).M();
        // assign to all sections where either photon cluster lies inside
        for (int is=0; is<nsec; ++is) {
            if (inSection(sections[is], d_ph1_cx, d_ph1_cy) || inSection(sections[is], d_ph2_cx, d_ph2_cy)) {
                data_mass_per_section[is].push_back(m);
            }
        }
    }

    // Determine a sensible global mass window from data mass percentiles (to avoid extreme tails)
    vector<double> all_masses;
    for (int is=0; is<nsec; ++is) {
        all_masses.insert(all_masses.end(), data_mass_per_section[is].begin(), data_mass_per_section[is].end());
    }
    if (all_masses.empty()) { cerr << "No data events found in any section. Check cluster branches and calorimeter bounds." << endl; return 6; }
    sort(all_masses.begin(), all_masses.end());
    double mlo = all_masses[ max(0, (int)(0.001*all_masses.size())) ];
    double mhi = all_masses[ min((int)all_masses.size()-1, (int)(0.999*all_masses.size())) ];
    // expand slightly
    double margin = 0.05*(mhi - mlo);
    mlo = max(0.0, mlo - margin);
    mhi = mhi + margin;
    int nbins = 200;
    cout << "Global mass window for histograms: ["<<mlo<<","<<mhi<<"] nbins="<<nbins<<"\n";

    // create TH1D per section for data and fill
    vector<TH1D*> hdata_sec(nsec,nullptr);
    for (int is=0; is<nsec; ++is) {
        string hname = string("h_data_") + sections[is].name();
        hdata_sec[is] = new TH1D(hname.c_str(), hname.c_str(), nbins, mlo, mhi);
        for (double m : data_mass_per_section[is]) hdata_sec[is]->Fill(m);
        cout << "Section "<<sections[is].name()<<" data entries="<<data_mass_per_section[is].size()<<"\n";
    }

    // Now load sim tree and build per-section sim event lists
    cout << "Scanning sim tree and building per-section sim event buffers..." << endl;
    vector<vector<pair<TLorentzVector,TLorentzVector>>> sim_events_per_section(nsec);
    Long64_t nsim = tsim->GetEntries();
    for (Long64_t i=0;i<nsim;++i) {
        tsim->GetEntry(i);
        TLorentzVector g1(s_ph1_px, s_ph1_py, s_ph1_pz, s_ph1_E);
        TLorentzVector g2(s_ph2_px, s_ph2_py, s_ph2_pz, s_ph2_E);
        for (int is=0; is<nsec; ++is) {
            if (inSection(sections[is], s_ph1_cx, s_ph1_cy) || inSection(sections[is], s_ph2_cx, s_ph2_cy)) {
                sim_events_per_section[is].emplace_back(make_pair(g1,g2));
            }
        }
    }
    for (int is=0; is<nsec; ++is) cout << "Section "<<sections[is].name()<<" sim entries="<<sim_events_per_section[is].size()<<"\n";

    // RNG
    TRandom3 rng(0);

    // Output file
    TFile fout(out_file.c_str(), "RECREATE");

    // CSV summary
    ofstream csv("section_map.csv");
    csv << "ix,iy,xlo,xhi,ylo,yhi,n_data,n_sim,best_mu,best_sigma,best_theta_deg,best_chi2\n";

    // Loop over sections and fit each where we have enough stats
    for (int is=0; is<nsec; ++is) {
        cout << "\n==== Fitting section "<<sections[is].name()<<" ===="<<endl;
        size_t ndata_sec = data_mass_per_section[is].size();
        size_t nsim_sec = sim_events_per_section[is].size();
        if (ndata_sec < 50 || nsim_sec < 50) {
            cout << "Skipping section "<<sections[is].name()<<" due to insufficient statistics (data="<<ndata_sec<<" sim="<<nsim_sec<<")"<<endl;
            // still write empty histograms for bookkeeping
            TDirectory *dir = fout.mkdir(sections[is].name().c_str()); dir->cd();
            hdata_sec[is]->Write(hdata_sec[is]->GetName());
            fout.cd();
            csv << sections[is].ix<<","<<sections[is].iy<<","<<sections[is].xlo<<","<<sections[is].xhi<<","<<sections[is].ylo<<","<<sections[is].yhi<<","<<ndata_sec<<","<<nsim_sec<<",,,,"<<"\n";
            continue;
        }

        // Fit
        FitResult fres = fit_section(sim_events_per_section[is], *hdata_sec[is], rng, Nsmear);
        cout << "Section "<<sections[is].name()<<" -> mu="<<fres.mu<<" sigma="<<fres.sigma<<" theta_deg="<<fres.theta_rad*180.0/M_PI<<" chi2="<<fres.chi2<<"\n";

        // Build final smeared sim histogram at best-fit for this section
        TH1D hfinal((string("h_smeared_best_")+sections[is].name()).c_str(), "h_smeared_best", nbins, mlo, mhi);
        for (const auto &ev : sim_events_per_section[is]) {
            TLorentzVector g1 = ev.first; TLorentzVector g2 = ev.second;
            TLorentzVector pi0 = g1 + g2; TVector3 beta = pi0.BoostVector();
            TLorentzVector g1_rest = g1; g1_rest.Boost(-beta);
            TLorentzVector g2_rest = g2; g2_rest.Boost(-beta);
            TVector3 v1 = g1_rest.Vect(); TVector3 v2 = g2_rest.Vect(); TVector3 axis = v1.Cross(v2);
            if (axis.Mag() >= 1e-12) {
                TVector3 axis_u = axis.Unit();
                TVector3 v1rot = rotateAroundAxis(v1, axis_u, +0.5*fres.theta_rad);
                TVector3 v2rot = rotateAroundAxis(v2, axis_u, -0.5*fres.theta_rad);
                g1_rest.SetVectM(v1rot.Unit()*v1.Mag(), 0.0);
                g2_rest.SetVectM(v2rot.Unit()*v2.Mag(), 0.0);
            }
            TLorentzVector g1_lab_rot = g1_rest; g1_lab_rot.Boost(beta);
            TLorentzVector g2_lab_rot = g2_rest; g2_lab_rot.Boost(beta);
            double E1 = g1_lab_rot.E(), E2 = g2_lab_rot.E(); if (E1<=0||E2<=0) continue;
            for (int k=0;k<Nsmear;++k) {
                double draw1 = rng.Gaus(fres.mu, fres.sigma); if (draw1<=0) draw1 = NONPOSITIVE_CLAMP;
                double draw2 = rng.Gaus(fres.mu, fres.sigma); if (draw2<=0) draw2 = NONPOSITIVE_CLAMP;
                double E1p = E1*draw1; double E2p = E2*draw2; double s1 = E1p/E1; double s2 = E2p/E2;
                TLorentzVector g1p(g1_lab_rot.Px()*s1, g1_lab_rot.Py()*s1, g1_lab_rot.Pz()*s1, E1p);
                TLorentzVector g2p(g2_lab_rot.Px()*s2, g2_lab_rot.Py()*s2, g2_lab_rot.Pz()*s2, E2p);
                hfinal.Fill((g1p+g2p).M());
            }
        }
        if (hfinal.Integral()>0 && hdata_sec[is]->Integral()>0) hfinal.Scale(hdata_sec[is]->Integral() / hfinal.Integral());

        // write per-section outputs
        TDirectory *dir = fout.mkdir(sections[is].name().c_str()); dir->cd();
        hdata_sec[is]->Write(hdata_sec[is]->GetName());
        hfinal.Write(hfinal.GetName());
        TNamed meta_mu((string("best_mu_")+sections[is].name()).c_str(), to_string(fres.mu).c_str()); meta_mu.Write();
        TNamed meta_sigma((string("best_sigma_")+sections[is].name()).c_str(), to_string(fres.sigma).c_str()); meta_sigma.Write();
        TNamed meta_theta((string("best_theta_deg_")+sections[is].name()).c_str(), to_string(fres.theta_rad*180.0/M_PI).c_str()); meta_theta.Write();
        TNamed meta_chi((string("best_chi2_")+sections[is].name()).c_str(), to_string(fres.chi2).c_str()); meta_chi.Write();
        fout.cd();

        csv << sections[is].ix<<","<<sections[is].iy<<","<<sections[is].xlo<<","<<sections[is].xhi<<","<<sections[is].ylo<<","<<sections[is].yhi<<","<<ndata_sec<<","<<nsim_sec<<","<<fres.mu<<","<<fres.sigma<<","<<fres.theta_rad*180.0/M_PI<<","<<fres.chi2<<"\n";
    }

    fout.Close();
    csv.close();
    cout << "All done. Results written to "<<out_file<<" and section_map.csv"<<endl;
    return 0;
}
