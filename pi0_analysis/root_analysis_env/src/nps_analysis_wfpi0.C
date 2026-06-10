// ============================================================================
// File: nps_analysis_wfpi0.C
// Purpose: Full diagnostic pipeline for NPS π0 analysis with robust outputs (Production/WF data)
// Author: Avnish Singh (physics), ChatGPT (refactoring 2026)
//
// Design:
//   - RAII wrappers for all ROOT objects
//   - Per-run memory isolation (histograms deleted after each run)
//   - Exception-safe event processing with catch-and-continue
//   - Seamless integration with header file utilities
//   - Physics calculations preserved exactly (no algorithm changes)
//
// Build: root -l -b nps_analysis_wfpi0.C
// ============================================================================

#include "utils.C"
#include <TF1.h>
#include "nps_helper.h"
#include "nps_time_bg.h"
#include "nps_comb_bg.h"
#include "nps_mmiss_cor.h"
#include "nps_2d_mass_cut.h"
#include "nps_hms_track_eff.h"
#include "nps_yao_database_reader.h"

// Suppress empty-body warning in physics_var.h
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wempty-body"
#include "nps_physics_var.h"
#pragma GCC diagnostic pop

// ROOT headers
#include <TFile.h>
#include <TApplication.h>
#include <TTree.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TParameter.h>
#include <TBox.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <TF1.h>
#include <TPaveText.h>

// C++ standard library
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstring>
#include <memory>
#include <sstream>
#include <iomanip>
#include <numeric>
#include <limits>
#include <stdexcept>
#include <cerrno>
#include <ctime>

using namespace std;

// ============================================================================
// Configuration Constants
// ============================================================================
constexpr int MAX_CLUS = 20;
constexpr double DEFAULT_TIME_THRESH_NS = 13.0;
constexpr double DEFAULT_TIME_WINDOW_WRT_150 = 13.0;
constexpr double EBEAM_DEFAULT = 10.538;
constexpr double HMS_MOM_OFFSET_SCALE = 1.0;
constexpr Long64_t MIN_PRINT_EVERY = 1000;

// ============================================================================
// Utility Functions
// ============================================================================

/// Template for safe deletion of pointers
template <typename T>
inline void safe_delete(T*& ptr)
{
    if (ptr) {
        delete ptr;
        ptr = nullptr;
    }
}

/// ============================================================================
/// RAII File Guard for Automatic TFile Resource Management
/// ============================================================================

class FileGuard {
private:
    TFile* file;
public:
    explicit FileGuard(TFile* f) : file(f) {}
    ~FileGuard() {
        if (file) {
            try {
                file->Close();
            } catch (...) {
                // Ignore errors during close
            }
            delete file;
        }
    }
    TFile* get() const { return file; }
    TFile* operator->() { return file; }
    
    FileGuard(const FileGuard&) = delete;
    FileGuard& operator=(const FileGuard&) = delete;
};

/// Print a progress bar with percentage
inline void print_progress(int run, const char* label, Long64_t current, Long64_t total) 
{
    const int bar_width = 40;
    double progress = (total > 0) ? (double)current / total : 0.0;
    int pos = (int)(bar_width * progress);
    
    cout << " Run " << run << " " << label << " [";
    for (int i = 0; i < bar_width; ++i) {
        if (i < pos) cout << "=";
        else if (i == pos) cout << ">";
        else cout << " ";
    }
    cout << "] " << int(progress * 100.0) << "% (" << current << "/" << total << ")\r";
    cout.flush();
}

/// Write global CSV header (appended to by each run)
inline bool write_global_csv_header(const TString &path) 
{
    std::ofstream f(path.Data(), std::ios::out | std::ios::trunc);
    if (!f.is_open()) {
        std::cerr << "[ERROR] Cannot open CSV for writing: " << path.Data()
                  << " (errno=" << errno << ", " << std::strerror(errno) << ")\n";
        return false;
    }
    f << "run,accumulated_charge(mC),current_mean_uA,CPUT_LT,Beam_Time(s),"
      << "total_entries,pass_hms,pass_hms_nps,total_coin_entries,"
      << "estimated_time_accidentals,chi2_ndf_comb_bg,pi0_mu_MeV,pi0_sigma_MeV,"
      << "pi0_signal_counts,mmiss_p_mean_GeV,mmiss_p_sigma_GeV,"
      << "hms_track_eff,hms_track_eff_err,s1x_peak,s1x_err,s1y_peak,s1y_err,s2x_peak,s2x_err,s2y_peak,s2y_err,run_status\n";
    f.close();
    return true;
}

// ============================================================================
// Main Analysis Macro
// ============================================================================
/// Comprehensive NPS π0 analysis with per-run processing and global summary (Production data)
///
/// @param skimDir_in          Directory containing production ROOT files
/// @param outPlotDir_in       Output directory for plots and histograms  
/// @param runlistFile         File listing runs to process (one per line)
/// @param Ebeam              Beam energy in GeV
void nps_analysis_wfpi0(const TString &skimDir_in = "/lustre24/expphy/volatile/hallc/nps/hhuang/farmFile/Production/DVCS",
                   const TString &outPlotDir_in = "output/plots/x60_4b/production_wfpi0",
                   const TString &runlistFile = "config/runlist_x60_4b.txt",
                   const double Ebeam = EBEAM_DEFAULT)
{
    // Timing and logging
    TStopwatch sw_total; 
    sw_total.Start();
    logmsg(INFO, "=========== NPS π0 FULL diagnostic analysis ===========");
    
    // Optimization: Set ROOT to batch mode for faster canvas operations
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    // Read parameters from environment variables if set (bypass shell quoting issues)
    TString skimDir = skimDir_in;
    TString outPlotDir = outPlotDir_in;
    TString runlist = runlistFile;
    double beam_energy = Ebeam;
    
    const char* env_skim = gSystem->Getenv("NPS_SKIM_DIR");
    const char* env_output = gSystem->Getenv("NPS_OUTPUT_DIR");
    const char* env_runlist = gSystem->Getenv("NPS_RUNLIST");
    const char* env_ebeam = gSystem->Getenv("NPS_EBEAM");
    
    if (env_skim) skimDir = TString(env_skim);
    if (env_output) outPlotDir = TString(env_output);
    if (env_runlist) runlist = TString(env_runlist);
    if (env_ebeam) beam_energy = atof(env_ebeam);
    
    // Normalize directory paths
    skimDir = skimDir.EndsWith("/") ? skimDir : skimDir + "/";
    outPlotDir = outPlotDir.EndsWith("/") ? outPlotDir : outPlotDir + "/";
    
    // Create output directories
    gSystem->mkdir("output", true);
    gSystem->mkdir(outPlotDir, true);

    // Load run list
    std::vector<int> runs = readRunList(runlist.Data());
    if (runs.empty()) {
        logmsg(ERROR, std::string("No runs found in: ") + runlist.Data());
        return;
    }
    logmsg(INFO, Form("Processing %zu runs", runs.size()));

    TString global_csv = outPlotDir + "summary_all_runs.csv";
    if (!write_global_csv_header(global_csv)) {
        TString fallback_csv = Form("/tmp/summary_all_runs_%d.csv", (int)time(nullptr));
        logmsg(WARN, Form("Falling back global CSV path to %s", fallback_csv.Data()));
        if (!write_global_csv_header(fallback_csv)) {
            logmsg(ERROR, "Failed to create global summary CSV in both output directory and /tmp");
            return;
        }
        global_csv = fallback_csv;
    }

    // ========================================================================
    // Per-run loop (process each run independently with memory cleanup)
    // ========================================================================
    for (int run : runs) {
        TStopwatch sw_run;
        sw_run.Start();

        // Construct input filename pattern (production files are split into multiple files per run)
        TString infile_pattern = skimDir + Form("nps_production_%d_*_wf_calib.root", run);

        // Validate at least one file exists for this run
        void* dirp = gSystem->OpenDirectory(skimDir);
        if (!dirp) {
            logmsg(WARN, Form("Cannot open directory %s", skimDir.Data()));
            continue;
        }
        
        bool found_file = false;
        const char* entry;
        while ((entry = gSystem->GetDirEntry(dirp))) {
            TString fname(entry);
            if (fname.Contains(Form("nps_production_%d_", run)) && fname.EndsWith("_wf_calib.root")) {
                found_file = true;
                break;
            }
        }
        gSystem->FreeDirectory(dirp);
        
        if (!found_file) {
            logmsg(WARN, Form("Skipping run %d: no files found matching pattern", run));
            continue;
        }

        logmsg(INFO, Form("Processing run %d", run));
        // ====================================================================
        // Per-run processing with RAII file management and exception safety
        // ====================================================================
        std::string run_status = "OK";  // Track whether run completes without errors
        try {
            // Open files with TChain to handle split files (automatic cleanup via unique_ptr)
            std::unique_ptr<TChain> chain(new TChain("t_prod"));
            
            // Add all split files for this run to the chain
            TString file_pattern = skimDir + Form("nps_production_%d_*_wf_calib.root", run);
            Int_t nfiles_added = chain->Add(file_pattern);
            
            if (nfiles_added == 0) {
                logmsg(ERROR, Form("Failed to add files to chain for run %d", run));
                continue;
            }
            
            logmsg(INFO, Form("Run %d: Added %d files to chain", run, nfiles_added));
            
            TTree *T = chain.get();  // Use the chain as a TTree
            if (!T) {
                logmsg(ERROR, Form("Chain is null for run %d", run));
                continue;
            }

            Long64_t nentries = T->GetEntries();
            logmsg(INFO, Form("Run %d: %lld entries", run, nentries));

            // -------------------------
            // Branch variables (no static arrays)
            // -------------------------
            Double_t HgtrX=0, HgtrY=0, HgtrTh=0, HgtrPh=0, hdelta=0, HgtrP=0, hreactz=0, hcernpeSum=0, hcaletotnorm=0;
            Double_t HgtrPx=0, HgtrPy=0, HgtrPz=0;
            Double_t edtmtdc=0;
            Int_t nclust_dbl = 0;
            Double_t BCM4A_scalerCurrent = 0, BCM4A_scalerCharge = 0, H_1MHz_scaler = 0;
            Double_t h_hodbetanotrack = 0, h_hodgoodscinhit = 0, hdcntrack = 0;
            Double_t s1x_rate = 0, s1y_rate = 0, s2x_rate = 0, s2y_rate = 0;
            
            // HMS focal plane variables
            Double_t hxfp=0, hyfp=0, hxpfp=0, hypfp=0;

            // local fixed arrays for convenience (fresh per run)
            Double_t clusE[MAX_CLUS] = {0}, clusX[MAX_CLUS] = {0}, clusY[MAX_CLUS] = {0}, clusT[MAX_CLUS] = {0};

            // vector bindings (safe if branches are stored as vectors)
            std::vector<double>* clusE_vec = nullptr;
            std::vector<double>* clusX_vec = nullptr;
            std::vector<double>* clusY_vec = nullptr;
            std::vector<double>* clusT_vec = nullptr;

            // set branch status / addresses (enable only needed branches)
            T->SetBranchStatus("*", 0);

            auto enable_scalar = [&](const char* b, void* addr){
                TBranch *br = T->GetBranch(b);
                if (br) { T->SetBranchStatus(b,1); T->SetBranchAddress(b, addr); }
            };

            // set scalar branches
            enable_scalar("H.gtr.x", &HgtrX);
            enable_scalar("H.gtr.y", &HgtrY);
            enable_scalar("H.gtr.p", &HgtrP);
            enable_scalar("H.gtr.px", &HgtrPx);
            enable_scalar("H.gtr.py", &HgtrPy);
            enable_scalar("H.gtr.pz", &HgtrPz);
            enable_scalar("H.gtr.dp", &hdelta);
            enable_scalar("H.gtr.th", &HgtrTh);
            enable_scalar("H.gtr.ph", &HgtrPh);
            enable_scalar("H.react.z", &hreactz);
            enable_scalar("H.cer.npeSum", &hcernpeSum);
            enable_scalar("H.cal.etotnorm", &hcaletotnorm);
            enable_scalar("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);
            enable_scalar("H.BCM4A.scalerCurrent", &BCM4A_scalerCurrent);
            enable_scalar("H.BCM4A.scalerCharge", &BCM4A_scalerCharge);
            enable_scalar("H.1MHz.scaler", &H_1MHz_scaler);
            enable_scalar("NPS.prod.nclust", &nclust_dbl);

            // HMS efficiencies
            enable_scalar("H.hod.betanotrack", &h_hodbetanotrack);
            enable_scalar("H.hod.goodscinhit", &h_hodgoodscinhit);
            enable_scalar("H.dc.ntrack", &hdcntrack);
            enable_scalar("H.S1X.scalerRate", &s1x_rate);
            enable_scalar("H.S1Y.scalerRate", &s1y_rate);
            enable_scalar("H.S2X.scalerRate", &s2x_rate);
            enable_scalar("H.S2Y.scalerRate", &s2y_rate);
            
            // HMS focal plane variables
            enable_scalar("H.dc.x_fp", &hxfp);
            enable_scalar("H.dc.y_fp", &hyfp);
            enable_scalar("H.dc.xp_fp", &hxpfp);
            enable_scalar("H.dc.yp_fp", &hypfp);

            // ----- clusters: try to bind to std::vector<double> if present; otherwise bind to arrays
            auto try_bind_cluster = [&](const char* bname, std::vector<double>** vecptr, Double_t arr[]) {
                TBranch *br = T->GetBranch(bname);
                if (!br) return;
                T->SetBranchStatus(bname,1);
                const char *cls = br->GetClassName();
                if (cls && (strstr(cls, "vector") || strstr(cls, "std::vector"))) {
                    T->SetBranchAddress(bname, vecptr);
                } else {
                    // fallback: attempt to bind to C-array (works if the file stored C-array like Double_t NPS.cal.clusE[MAX_CLUS])
                    T->SetBranchAddress(bname, arr);
                }
            };

            try_bind_cluster("NPS.prod.clusE", &clusE_vec, clusE);
            try_bind_cluster("NPS.prod.clusXcorr", &clusX_vec, clusX);
            try_bind_cluster("NPS.prod.clusYcorr", &clusY_vec, clusY);
            try_bind_cluster("NPS.prod.clusT", &clusT_vec, clusT);

            cout << "Run " << run << " entries: " << nentries << endl;

            // Get beam current mean from database
            double run_current_mean = nps::getAve_BeamCurr_or_nan(run);

            // Read charge & livetimes from DB. If missing or invalid, skip this run
        double accumulated_charge_uC = nps::getChargeTot_or_nan(run);  // expects uC
        if (!std::isfinite(accumulated_charge_uC) || accumulated_charge_uC <= 0.0) {
            logmsg(WARN, Form("Charge not found or invalid for run %d in database - skipping run", run));
            continue;  // FileGuard handles cleanup automatically
        }
        double accumulated_charge_mC = accumulated_charge_uC / 1000.0;

        double cpu_lt = nps::getCPU_LT_or_nan(run);
        if (!std::isfinite(cpu_lt) || cpu_lt <= 0 || cpu_lt > 1.05) {
            logmsg(WARN, Form("Bad CPU_LT for run %d: %.6g - skipping run", run, cpu_lt));
            continue;  // FileGuard handles cleanup automatically
        }

        double beam_time = nps::getBeam_Time_or_nan(run);
        if (!std::isfinite(beam_time) || beam_time <= 0) {
            logmsg(WARN, Form("Bad Beam_Time for run %d: %.6g - skipping run", run, beam_time));
            continue;  // FileGuard handles cleanup automatically
        }

        cout << "Charge = " << accumulated_charge_mC << " mC\n";
        cout << "CPU_LT = " << cpu_lt << "\n";
        cout << "Beam_Time = " << beam_time << " s\n";

        // -------------------------
        // Histograms (unique names per run). SetDirectory(nullptr) to avoid ROOT file ownership issues
        // -------------------------
        auto name = [&](const char* base){ return TString::Format("%s_run%d", base, run); };

        auto make1D = [&](const char* base, const char* title, int nb, double xlo, double xhi)->TH1D* {
        TH1D *h = new TH1D(name(base), title, nb, xlo, xhi);
        h->SetDirectory(nullptr);
        return h;
        };
        auto make2D = [&](const char* base, const char* title, int nbx, double xlo, double xhi, int nby, double ylo, double yhi)->TH2D* {
        TH2D *h = new TH2D(name(base), title, nbx, xlo, xhi, nby, ylo, yhi);
        h->SetDirectory(nullptr);
        return h;
        };

        TH1D *h_nclusters = make1D("h_nclusters","NPS clusters per event;N_{clus};Events",21, -0.5, 20.5);
        TH1D *h_clustE = make1D("h_clustE","Cluster energy;E_{clus} [GeV];Counts",200, 0.0, 4.0);
        TH1D *h_clustT = make1D("h_clustT","Cluster time; t [ns];Counts",200, 120, 180);
        TH2D *h_clustE_vs_T = make2D("h_clustE_vs_T","Cluster E vs t; t [ns];E [GeV]",200,120,180,200,0,8);
        TH2D *h_clustXY = make2D("h_clustXY","Cluster X vs Y; X [cm]; Y [cm]",30,-30,30,36,-36,36);
        TH1D *h_clustE_sum = make1D("h_clustE_sum","Sum cluster E per event;E_{sum} [GeV];Events",200,0,9);
        TH1D *h_opening_angle = make1D("h_opening_angle","Opening angle between two photons;#theta_{open} [rad];Counts",200,0,0.5);
        TH1D *h_photon_Eratio = make1D("h_photon_Eratio","Photon energy ratio E1/E2;E1/E2;Counts",100,0,5);

        TH1D *h_mpi0_all = make1D("h_mpi0_all","Invariant mass (all best-pairs);M_{#gamma#gamma} [GeV];Events",200,0.0,0.4);
        TH1D *h_mpi0_2 = make1D("h_mpi0_2","Invariant mass (2-cluster);M [GeV];Events",200,0.0,0.4);
        TH1D *h_mpi0_3 = make1D("h_mpi0_3","Invariant mass (3-cluster best pair);M [GeV];Events",200,0.0,0.4);
        TH1D *h_mpi0_4 = make1D("h_mpi0_4","Invariant mass (4-cluster best pair);M [GeV];Events",200,0.0,0.4);

        TH1D *h_mmiss_2 = make1D("h_mmiss_2","Missing mass (2-cluster);M_{miss} [GeV];Events",200,0.0,2.0);
        TH1D *h_mmiss_3 = make1D("h_mmiss_3","Missing mass (3-cluster);M_{miss} [GeV];Events",200,0.0,2.0);
        TH1D *h_mmiss_4 = make1D("h_mmiss_4","Missing mass (4-cluster);M_{miss} [GeV];Events",200,0.0,2.0);
        TH1D *h_mmiss_dvcs = make1D("h_mmiss_dvcs","Missing mass (DVCS-like single photon);M_{miss} [GeV];Events",200,0.0,2.0);

        TH1D *h_mmiss_all = make1D("h_mmiss_all","Missing mass;M_{miss} [GeV];Events",200,0.0,2.0);
        TH1D *h_mmiss_all_corr = make1D("h_mmiss_all_corr","Missing mass;M_{miss} [GeV];Events",200,0.0,2.0);

        const double t_min = 140.0, t_max = 160.0;
        const int nbins_t = 200;
        TH2D *h_t1_t2 = make2D("h_t1_t2", "t1 (y) vs t2 (x);t2 [ns];t1 [ns]", nbins_t, t_min, t_max, nbins_t, t_min, t_max);
        TH1D *h_t1_proj = make1D("h_t1_proj", "t1 projection; t1 [ns];Entries", nbins_t, t_min, t_max);
        TH1D *h_t2_proj = make1D("h_t2_proj", "t2 projection; t2 [ns];Entries", nbins_t, t_min, t_max);

        // coin / acc histos and per-box mgg histograms
        TH1D *h_m_pi0_coin = make1D("h_m_pi0_coin","Pi0 mass (Coincidence);M [GeV];Counts",200,0.0,0.4);
        TH1D *h_m_pi0_acc  = make1D("h_m_pi0_acc","Pi0 mass (Accidentals - outside coin);M [GeV];Counts",200,0.0,0.4);
        TH1D *h_coin_bgsub = nullptr; // assigned by helper; caller owns pointer if created

        // Physics variable histograms (unique names)
        TH1D *h_Q2 = make1D("h_Q2","Q^{2};Q^{2}  [GeV^{2}];Counts",200,0,10);
        TH1D *h_W  = make1D("h_W","W;W  [GeV];Counts",200,0.5,4.5);
        TH1D *h_t  = make1D("h_t","t;t  [GeV^{2}];Counts",200,-5.0,0.5);
        TH1D *h_tmin = make1D("h_tmin","t_{min};t_{min}  [GeV^{2}];Counts",200,-5.0,0.0);
        TH1D *h_pt = make1D("h_pt","p_{T};p_{T}  [GeV];Counts",200,0.0,1.0);
        TH1D *h_theta = make1D("h_theta","#theta;#theta  [rad];Counts",180,0.0,0.5);
        TH1D *h_phi = make1D("h_phi","#phi;#phi  [rad];Counts",180,-3.2,3.2);
        TH1D *h_s = make1D("h_s","s;s  [GeV^{2}];Counts",200,5.0,30.0);
        TH1D *h_xB = make1D("h_xB","x_{B};x_{B};Counts",200,0.0,1.0);
        TH1D *h_z = make1D("h_z","z;z;Counts",200,0.0,1.2);

        // -----------------------------------------------
        // π⁰ statistical background subtraction histograms
        // -----------------------------------------------
        // Store tree data in memory indexed by event_id
        // Tree will be created and filled at the end during file write
        // -----------------------------------------------
        struct TreeEntry {
            Double_t mpi0_all, mmiss_all, mmiss_all_corr, mmiss_all_no_mom_offset, mmiss_all_corr_no_mom_offset, pi0_weight;
            Int_t is_exclusive, is_exclusive_ellipse, is_exclusive_mcd, is_weighted;
            Double_t Q2, W, t, tmin, pt, theta, phi, s, xB, z;
            Int_t nclust_selected, event_id;
            Double_t cluster_x_1, cluster_y_1, cluster_e_1, cluster_x_2, cluster_y_2, cluster_e_2;
            // HMS tracking variables
            Double_t delta, xptar, yptar, xtar, ytar, xfp, yfp, xpfp, ypfp;
        };
        std::map<Long64_t, TreeEntry> treeData;  // Indexed by event_id

        // -----------------------------------------------
        // Bin-by-bin weight histogram (weight_i = N_final(i) / N_all(i))
        // Used to statistically subtract π⁰ background on event-by-event basis
        TH1D *h_pi0_weight = make1D("h_pi0_weight","π⁰ bin-by-bin weight;M_{γγ} [GeV];Weight",200,0.0,0.4);

        // Weighted physics distributions (without exclusivity cut - for overlay comparison)
        TH1D *h_Q2_weighted = make1D("h_Q2_weighted","Q^{2} (weighted);Q^{2} [GeV^{2}];Counts",200,0,10);
        TH1D *h_W_weighted = make1D("h_W_weighted","W (weighted);W [GeV];Counts",200,0.5,4.5);
        TH1D *h_t_weighted = make1D("h_t_weighted","t (weighted);t [GeV^{2}];Counts",200,-5.0,0.5);
        TH1D *h_tmin_weighted = make1D("h_tmin_weighted","t_{min} (weighted);t_{min} [GeV^{2}];Counts",200,-5.0,0.0);
        TH1D *h_pt_weighted = make1D("h_pt_weighted","p_{T} (weighted);p_{T} [GeV];Counts",200,0.0,1.0);
        TH1D *h_theta_weighted = make1D("h_theta_weighted","θ (weighted);θ [rad];Counts",180,0.0,0.5);
        TH1D *h_phi_weighted = make1D("h_phi_weighted","φ (weighted);φ [rad];Counts",180,-3.2,3.2);
        TH1D *h_s_weighted = make1D("h_s_weighted","s (weighted);s [GeV^{2}];Counts",200,5.0,30.0);
        TH1D *h_xB_weighted = make1D("h_xB_weighted","x_{B} (weighted);x_{B};Counts",200,0.0,1.0);
        TH1D *h_z_weighted = make1D("h_z_weighted","z (weighted);z;Counts",200,0.0,1.2);

        // Weighted+exclusive physics distributions (final signal sample)
        TH1D *h_Q2_excl_weighted = make1D("h_Q2_excl_weighted","Q^{2} (exclusive, weighted);Q^{2} [GeV^{2}];Counts",200,0,10);
        TH1D *h_W_excl_weighted = make1D("h_W_excl_weighted","W (exclusive, weighted);W [GeV];Counts",200,0.5,4.5);
        TH1D *h_t_excl_weighted = make1D("h_t_excl_weighted","t (exclusive, weighted);t [GeV^{2}];Counts",200,-5.0,0.5);
        TH1D *h_tmin_excl_weighted = make1D("h_tmin_excl_weighted","t_{min} (exclusive, weighted);t_{min} [GeV^{2}];Counts",200,-5.0,0.0);
        TH1D *h_pt_excl_weighted = make1D("h_pt_excl_weighted","p_{T} (exclusive, weighted);p_{T} [GeV];Counts",200,0.0,1.0);
        TH1D *h_theta_excl_weighted = make1D("h_theta_excl_weighted","θ (exclusive, weighted);θ [rad];Counts",180,0.0,0.5);
        TH1D *h_phi_excl_weighted = make1D("h_phi_excl_weighted","φ (exclusive, weighted);φ [rad];Counts",180,-3.2,3.2);
        TH1D *h_s_excl_weighted = make1D("h_s_excl_weighted","s (exclusive, weighted);s [GeV^{2}];Counts",200,5.0,30.0);
        TH1D *h_xB_excl_weighted = make1D("h_xB_excl_weighted","x_{B} (exclusive, weighted);x_{B};Counts",200,0.0,1.0);
        TH1D *h_z_excl_weighted = make1D("h_z_excl_weighted","z (exclusive, weighted);z;Counts",200,0.0,1.2);

        // Missing mass histograms (for overlay comparison)
        TH1D *h_mmiss_all_weighted = make1D("h_mmiss_all_weighted","Missing mass (weighted);M_{miss} [GeV];Counts",200,0.0,2.0);

        // -----------------------------------------------
        // 2D Physics variables on x-y cluster map (detector heatmaps)
        // -----------------------------------------------
        TH2D *h_xB_xy = make2D("h_xB_xy", "x_{B} distribution on detector; Cluster X [cm]; Cluster Y [cm]", 60, -30, 30, 72, -36, 36);
        TH2D *h_nu_xy = make2D("h_nu_xy", "#nu (E_{beam} - E_{scatter}) on detector; Cluster X [cm]; Cluster Y [cm]", 60, -30, 30, 72, -36, 36);
        TH2D *h_mx_xy = make2D("h_mx_xy", "Missing Mass on detector; Cluster X [cm]; Cluster Y [cm]", 60, -30, 30, 72, -36, 36);
        TH2D *h_phi_xy = make2D("h_phi_xy", "#phi distribution on detector; Cluster X [cm]; Cluster Y [cm]", 60, -30, 30, 72, -36, 36);
        TH2D *h_t_xy = make2D("h_t_xy", "-t distribution on detector; Cluster X [cm]; Cluster Y [cm]", 60, -30, 30, 72, -36, 36);

        // -----------------------------------------------
        // 2D correlations of final weighted+exclusive physics variables
        // -----------------------------------------------
        TH2D *h_Q2_vs_W = make2D("h_Q2_vs_W", "Q^{2} vs W (excl weighted); W [GeV]; Q^{2} [GeV^{2}]", 150, 0.5, 4.5, 150, 0, 10);
        TH2D *h_Q2_vs_t = make2D("h_Q2_vs_t", "Q^{2} vs -t (excl weighted); -t [GeV^{2}]; Q^{2} [GeV^{2}]", 150, 0, 5, 150, 0, 10);
        TH2D *h_Q2_vs_xB = make2D("h_Q2_vs_xB", "Q^{2} vs x_{B} (excl weighted); x_{B}; Q^{2} [GeV^{2}]", 150, 0, 1, 150, 0, 10);
        TH2D *h_W_vs_t = make2D("h_W_vs_t", "W vs -t (excl weighted); -t [GeV^{2}]; W [GeV]", 150, 0, 5, 150, 0.5, 4.5);
        TH2D *h_W_vs_xB = make2D("h_W_vs_xB", "W vs x_{B} (excl weighted); x_{B}; W [GeV]", 150, 0, 1, 150, 0.5, 4.5);
        TH2D *h_t_vs_xB = make2D("h_t_vs_xB", "-t vs x_{B} (excl weighted); x_{B}; -t [GeV^{2}]", 150, 0, 1, 150, 0, 5);
        TH2D *h_t_vs_phi = make2D("h_t_vs_phi", "-t vs #phi (excl weighted); #phi [rad]; -t [GeV^{2}]", 180, -3.2, 3.2, 150, 0, 5);
        TH2D *h_xB_vs_nu = make2D("h_xB_vs_nu", "x_{B} vs #nu (excl weighted); #nu [GeV]; x_{B}", 200, 0, 12, 150, 0, 1);
        TH2D *h_Q2_vs_nu = make2D("h_Q2_vs_nu", "Q^{2} vs #nu (excl weighted); #nu [GeV]; Q^{2} [GeV^{2}]", 200, 0, 12, 150, 0, 10);
        TH2D *h_W_vs_nu = make2D("h_W_vs_nu", "W vs #nu (excl weighted); #nu [GeV]; W [GeV]", 200, 0, 12, 150, 0.5, 4.5);

        // Per-window mgg histograms (diag/side/full)
        // Adjusted timing windows for production data (lower windows -2ns, upper windows +2ns)
        vector<pair<double,double>> diag_windows = {
            {139.0, 141.0}, {141.0, 143.0}, {143.0, 145.0},
            {155.0, 157.0}, {157.0, 159.0}, {159.0, 161.0}
        };
        vector<pair<double,double>> side_windows = diag_windows;  // same as diag
        auto coin_win = nps::default_coin_window();  // Keep unchanged: {149.0, 151.0}
        auto full1_t1 = std::make_pair(155.0, 161.0);  // +2ns from default {153.0, 159.0}
        auto full1_t2 = std::make_pair(139.0, 145.0);  // -2ns from default {141.0, 147.0}
        auto full2_t1 = std::make_pair(139.0, 145.0);  // -2ns from default {141.0, 147.0}
        auto full2_t2 = std::make_pair(155.0, 161.0);  // +2ns from default {153.0, 159.0}

        vector<TH1D*> h_mgg_diag;
        for (size_t i=0;i<diag_windows.size();++i) {
        h_mgg_diag.push_back(make1D(TString::Format("h_mgg_diag%d",(int)i).Data(),
                                    TString::Format("m_{#gamma#gamma} diag %d;M [GeV];Counts",(int)i).Data(),
                                    200,0.0,0.4));
        }
        vector<TH1D*> h_mgg_hor; // horizontal sideboxes (t1=coin, t2 in side)
        vector<TH1D*> h_mgg_ver; // vertical sideboxes (t2=coin, t1 in side)
        for (size_t i=0;i<side_windows.size();++i) {
        h_mgg_hor.push_back(make1D(TString::Format("h_mgg_hor%d",(int)i).Data(),
                                   TString::Format("m_{#gamma#gamma} hor %d;M [GeV];Counts",(int)i).Data(),200,0.0,0.4));
        h_mgg_ver.push_back(make1D(TString::Format("h_mgg_ver%d",(int)i).Data(),
                                   TString::Format("m_{#gamma#gamma} ver %d;M [GeV];Counts",(int)i).Data(),200,0.0,0.4));
        }
        TH1D *h_mgg_full1 = make1D("h_mgg_full1","m_{#gamma#gamma} fullbox1;M [GeV];Counts",200,0.0,0.4);
        TH1D *h_mgg_full2 = make1D("h_mgg_full2","m_{#gamma#gamma} fullbox2;M [GeV];Counts",200,0.0,0.4);

        // mm vs mgg scatter
        TH2D *h_mmiss_vs_mgg = make2D("h_mmiss_vs_mgg","M_{miss} vs M_{#gamma#gamma};M_{miss} [GeV];M_{#gamma#gamma} [GeV]",200,0.7,1.8,200,0.10,0.16);
        TH2D *h_mmiss_vs_mgg_corr = make2D("h_mmiss_vs_mgg_corr","M_{miss} vs M_{#gamma#gamma};M_{miss} [GeV];M_{#gamma#gamma} [GeV]",200,0.7,1.8,200,0.10,0.16);

        // mm vs t1/t2
        TH2D *h_mmiss_vs_t1 = make2D("h_mmiss_vs_t1","M_{miss} vs t1; t1 [ns]; M_{miss} [GeV]", nbins_t, t_min, t_max, 200, 0.0, 2.0);
        TH2D *h_mmiss_vs_t2 = make2D("h_mmiss_vs_t2","M_{miss} vs t2; t2 [ns]; M_{miss} [GeV]", nbins_t, t_min, t_max, 200, 0.0, 2.0);

        // bookkeeping counters
        Long64_t n_total = 0, n_pass_hms = 0, n_ge2_hms = 0;
        Long64_t n_mult2_hms = 0, n_mult3_hms = 0, n_mult4_hms = 0;
        Long64_t n_mult2_hms_nps = 0, n_mult3_hms_nps = 0, n_mult4_hms_nps = 0;
        Long64_t n_dvcs_flagged = 0, n_selected_for_analysis = 0;

        // hms efficiency first call
        nps::hms_track_eff::HMSTrackingEffCounter hmstrackEff(run);

        // -------------------------
        // Event loop (exception-safe processing)
        // -------------------------
        const Long64_t print_every = std::max(MIN_PRINT_EVERY, nentries/10);  // Print 10 times instead of 100
        Long64_t n_event_errors = 0;
        
        try {
        for (Long64_t ev=0; ev<nentries; ++ev) {
            try {
                if ((ev % print_every) == 0) {
                    print_progress(run, "Processing", ev, nentries);
                }

                // read entry
                T->GetEntry(ev);
                ++n_total;

                // DEBUG: Print cluster info for first few events
                static int debug_count = 0;
                bool do_debug = (debug_count < 10);

                // If vectors were used for clusters, copy into arrays (safe)
                int nclust = 0;
                if (clusE_vec && clusX_vec && clusY_vec && clusT_vec) {
                    nclust = (int) std::min<size_t>(clusE_vec->size(), MAX_CLUS);
                    if (do_debug && nclust > 0) {
                        cout << "\n[DEBUG] Event " << ev << ": Vector mode, nclust=" << nclust 
                             << " (vec sizes: E=" << clusE_vec->size() 
                             << " X=" << clusX_vec->size() 
                             << " Y=" << clusY_vec->size() 
                             << " T=" << clusT_vec->size() << ")\n";
                        debug_count++;
                    }
                    for (int i=0;i<nclust;++i) {
                        clusE[i] = (*clusE_vec)[i];
                        clusX[i] = (*clusX_vec)[i];
                        clusY[i] = (*clusY_vec)[i];
                        clusT[i] = (*clusT_vec)[i] + 150.0;  // Add 150 ns offset for production data (reference time is 0, not 150)
                    }
                } else {
                    // fallback: use nclust_dbl + arrays (if the tree stored arrays)
                    nclust = nclust_dbl;
                    if (do_debug && ev < 100) {
                        cout << "\n[DEBUG] Event " << ev << ": Array mode, nclust_dbl=" << nclust_dbl 
                             << " (vec pointers: E=" << (void*)clusE_vec 
                             << " X=" << (void*)clusX_vec 
                             << " Y=" << (void*)clusY_vec 
                             << " T=" << (void*)clusT_vec << ")\n";
                        debug_count++;
                    }
                    if (nclust < 0) nclust = 0;
                    if (nclust > MAX_CLUS) nclust = MAX_CLUS;
                }

                // quick validation to catch bad branch bindings early
                if (nclust > MAX_CLUS) {
                    std::cerr << "Warning: nclust > MAX_CLUS at run " << run << " ev " << ev << " (" << nclust << ")\n";
                    nclust = std::min(nclust, MAX_CLUS);
                }

                // HMS electron selection
                if (!nps::hms_electron_cuts(edtmtdc, hdelta, HgtrTh, HgtrPh, hcernpeSum, hcaletotnorm, hreactz)) continue;
                ++n_pass_hms;

                // HMS efficiencies per run
                hmstrackEff.processEvent(
                    h_hodbetanotrack,
                    hcaletotnorm,
                    hcernpeSum,
                    h_hodgoodscinhit,
                    hdcntrack,
                    s1x_rate,
                    s1y_rate,
                    s2x_rate,
                    s2y_rate
                );

                // fill nclusters
                h_nclusters->Fill(nclust);
                if (nclust < 2) continue;

                ++n_ge2_hms;
                if (nclust == 2) ++n_mult2_hms;
                else if (nclust == 3) ++n_mult3_hms;
                else if (nclust == 4) ++n_mult4_hms;

                // per-cluster diagnostics
                double Esum = 0;
                for (int i=0;i<nclust;++i) {
                    if (std::isfinite(clusE[i])) h_clustE->Fill(clusE[i]);
                    if (std::isfinite(clusT[i])) h_clustT->Fill(clusT[i]);
                    h_clustE_vs_T->Fill(clusT[i], clusE[i]);
                    h_clustXY->Fill(clusX[i], clusY[i]);
                    Esum += clusE[i];
                }
                h_clustE_sum->Fill(Esum);

                // pack clusters (your helper) - returns number after packing
                const int n_after = nps::packClusters(clusE, clusX, clusY, clusT, nclust);

                // select good clusters by spatial & energy cuts
                vector<int> good_idx; good_idx.reserve(8);
                for (int i=0;i<n_after;++i) {
                    if (nps::nps_spatial_energy_cuts(clusE[i], clusX[i], clusY[i], clusT[i], DEFAULT_TIME_WINDOW_WRT_150))
                        good_idx.push_back(i);
                }

                // DEBUG: Show why clusters are failing cuts
                static int debug_cuts = 0;
                static int debug_timing = 0;
                if (debug_timing < 10 && n_after >= 2) {
                    cout << "\n[DEBUG TIMING] Event " << ev << ": n_after=" << n_after << "\n";
                    for (int i=0; i<n_after; ++i) {
                        bool passes = nps::nps_spatial_energy_cuts(clusE[i], clusX[i], clusY[i], clusT[i], DEFAULT_TIME_WINDOW_WRT_150);
                        cout << "  Cluster " << i << ": E=" << clusE[i] << " GeV, X=" << clusX[i] 
                             << " cm, Y=" << clusY[i] << " cm, T=" << clusT[i] << " ns"
                             << " [" << (passes ? "PASS" : "FAIL") << "]\n";
                    }
                    debug_timing++;
                }
                if (debug_cuts < 5 && n_after >= 2 && good_idx.size() < 2) {
                    cout << "\n[DEBUG CUT FAILURE] Event " << ev << ": n_after=" << n_after 
                         << " but good_idx.size()=" << good_idx.size() << "\n";
                    for (int i=0; i<n_after; ++i) {
                        cout << "  Cluster " << i << ": E=" << clusE[i] << " GeV, X=" << clusX[i] 
                             << " cm, Y=" << clusY[i] << " cm, T=" << clusT[i] << " ns\n";
                    }
                    debug_cuts++;
                }

                if (good_idx.size() < 2) continue;
                if (good_idx.size() > 4) {
                    sort(good_idx.begin(), good_idx.end(), [&](int a,int b){ return clusE[a] > clusE[b]; });
                    good_idx.resize(4);
                }

                if (good_idx.size() == 2) ++n_mult2_hms_nps;
                else if (good_idx.size() == 3) ++n_mult3_hms_nps;
                else if (good_idx.size() == 4) ++n_mult4_hms_nps;

                const double px_e = HgtrPx, py_e = HgtrPy, pz_e = HgtrPz;
                const double p_e_mom = sqrt(max(0.0, px_e*px_e + py_e*py_e + pz_e*pz_e));
                const double Ee = sqrt(max(0.0, p_e_mom*p_e_mom + nps::kElectronMass_GeV*nps::kElectronMass_GeV));

                // select pair for pi0
                int sel_i=-1, sel_j=-1;
                if (good_idx.size() == 2) { sel_i = good_idx[0]; sel_j = good_idx[1]; }
                else {
                    auto pr = nps::choose_best_pair_closest_pi0(good_idx, clusE, clusX, clusY, clusT, nps::kDefaultZ_NPS_cm, nps::kPi0Mass_GeV, DEFAULT_TIME_THRESH_NS);
                    sel_i = pr.first; sel_j = pr.second;
                    if (sel_i<0 || sel_j<0) continue;
                }

                // fill timing vs timing
                const double t1 = clusT[sel_i], t2 = clusT[sel_j];
                h_t1_t2->Fill(t2, t1);
                h_t1_proj->Fill(t1);
                h_t2_proj->Fill(t2);

                // invariant mass
                const double mgg = nps::invariant_mass_pi0(clusE[sel_i], clusE[sel_j], clusX[sel_i], clusX[sel_j], clusY[sel_i], clusY[sel_j], nps::kDefaultZ_NPS_cm);
                h_mpi0_all->Fill(mgg);
                if (good_idx.size()==2) h_mpi0_2->Fill(mgg);
                else if (good_idx.size()==3) h_mpi0_3->Fill(mgg);
                else if (good_idx.size()==4) h_mpi0_4->Fill(mgg);

                // opening angle & energy ratio
                {
                    const double z_nps = nps::kDefaultZ_NPS_cm;
                    const double r1x = clusX[sel_i], r1y = clusY[sel_i], r1z = z_nps;
                    const double r2x = clusX[sel_j], r2y = clusY[sel_j], r2z = z_nps;
                    const double n1norm = sqrt(r1x*r1x + r1y*r1y + r1z*r1z);
                    const double n2norm = sqrt(r2x*r2x + r2y*r2y + r2z*r2z);
                    double theta_open = 0.0;
                    if (n1norm>0.0 && n2norm>0.0) {
                        double ux1=r1x/n1norm, uy1=r1y/n1norm, uz1=r1z/n1norm;
                        double ux2=r2x/n2norm, uy2=r2y/n2norm, uz2=r2z/n2norm;
                        double dp = ux1*ux2 + uy1*uy2 + uz1*uz2;
                        dp = std::min(1.0, std::max(-1.0, dp));
                        theta_open = acos(dp);
                    }
                    h_opening_angle->Fill(theta_open);
                }

                double eratio = (clusE[sel_j] > 0) ? (clusE[sel_i]/clusE[sel_j]) : 0.0;
                h_photon_Eratio->Fill(eratio);

                // missing mass (proton)
                const double px_e_scaled = px_e * HMS_MOM_OFFSET_SCALE;
                const double py_e_scaled = py_e * HMS_MOM_OFFSET_SCALE;
                const double pz_e_scaled = pz_e * HMS_MOM_OFFSET_SCALE;
                const double p_e_mom_scaled = sqrt(max(0.0, px_e_scaled*px_e_scaled + py_e_scaled*py_e_scaled + pz_e_scaled*pz_e_scaled));
                const double Ee_scaled = sqrt(max(0.0, p_e_mom_scaled*p_e_mom_scaled + nps::kElectronMass_GeV*nps::kElectronMass_GeV));

                const double mm_p_no_offset = nps::missing_mass_proton_pi0(Ebeam, Ee, px_e, py_e, pz_e,
                                                                          clusE[sel_i], clusE[sel_j],
                                                                          clusX[sel_i], clusY[sel_i],
                                                                          clusX[sel_j], clusY[sel_j],
                                                                          nps::kDefaultZ_NPS_cm, -17.51);

                const double mm_p = nps::missing_mass_proton_pi0(Ebeam, Ee_scaled, px_e_scaled, py_e_scaled, pz_e_scaled,
                                                                clusE[sel_i], clusE[sel_j],
                                                                clusX[sel_i], clusY[sel_i],
                                                                clusX[sel_j], clusY[sel_j],
                                                                nps::kDefaultZ_NPS_cm, -17.51);

                const double mm_p_corr = nps::invariant_missing_mass_corrected_avnish_from_detector(mm_p*mm_p, Ebeam, Ee_scaled, px_e_scaled, py_e_scaled, pz_e_scaled,
                                                                clusE[sel_i], clusE[sel_j],
                                                                clusX[sel_i], clusY[sel_i],
                                                                clusX[sel_j], clusY[sel_j],
                                                                nps::kDefaultZ_NPS_cm, -17.51);

                const double mm_p_corr_no_mom_offset = nps::invariant_missing_mass_corrected_avnish_from_detector(mm_p_no_offset*mm_p_no_offset, Ebeam, Ee, px_e, py_e, pz_e,
                                                                clusE[sel_i], clusE[sel_j],
                                                                clusX[sel_i], clusY[sel_i],
                                                                clusX[sel_j], clusY[sel_j],
                                                                nps::kDefaultZ_NPS_cm, -17.51);

                if (good_idx.size()==2) h_mmiss_2->Fill(mm_p);
                else if (good_idx.size()==3) h_mmiss_3->Fill(mm_p);
                else if (good_idx.size()==4) h_mmiss_4->Fill(mm_p);

                h_mmiss_vs_mgg->Fill(mm_p, mgg);
                h_mmiss_vs_mgg_corr->Fill(mm_p_corr, mgg);
                h_mmiss_all->Fill(mm_p);
                h_mmiss_all_corr->Fill(mm_p_corr);

                h_mmiss_vs_t1->Fill(clusT[sel_i], mm_p);
                h_mmiss_vs_t2->Fill(clusT[sel_j], mm_p);
                ++n_selected_for_analysis;

                // coin vs acc classification (both t1 & t2 in coin -> coin)
                bool in_coin = (t1 > coin_win.first && t1 < coin_win.second && t2 > coin_win.first && t2 < coin_win.second);
                if (in_coin) h_m_pi0_coin->Fill(mgg);
                else h_m_pi0_acc->Fill(mgg);

                // side/diag/full windows fill
                for (size_t i=0;i<diag_windows.size();++i) {
                    const auto &w = diag_windows[i];
                    if (t1 > w.first && t1 < w.second && t2 > w.first && t2 < w.second) { h_mgg_diag[i]->Fill(mgg); break; }
                }
                for (size_t i=0;i<side_windows.size();++i) {
                    const auto &w = side_windows[i];
                    if ( (t1 > coin_win.first && t1 < coin_win.second) && (t2 > w.first && t2 < w.second) ) { h_mgg_hor[i]->Fill(mgg); break; }
                }
                for (size_t i=0;i<side_windows.size();++i) {
                    const auto &w = side_windows[i];
                    if ( (t2 > coin_win.first && t2 < coin_win.second) && (t1 > w.first && t1 < w.second) ) { h_mgg_ver[i]->Fill(mgg); break; }
                }
                if ( (t1 > full1_t1.first && t1 < full1_t1.second && t2 > full1_t2.first && t2 < full1_t2.second) ) h_mgg_full1->Fill(mgg);
                if ( (t1 > full2_t1.first && t1 < full2_t1.second && t2 > full2_t2.first && t2 < full2_t2.second) ) h_mgg_full2->Fill(mgg);

                auto phys = nps::compute_physics_vars_from_detector(
                    beam_energy, Ee, px_e, py_e, pz_e,
                    clusE[sel_i], clusX[sel_i], clusY[sel_i],
                    clusE[sel_j], clusX[sel_j], clusY[sel_j],
                    nps::kDefaultZ_NPS_cm, -17.51, false);

                // fill histograms you created
                if (std::isfinite(phys.t)) h_t->Fill(phys.t);
                if (std::isfinite(phys.tmin)) h_tmin->Fill(phys.tmin);
                if (std::isfinite(phys.pt)) h_pt->Fill(phys.pt);
                if (std::isfinite(phys.Q2)) h_Q2->Fill(phys.Q2);
                if (std::isfinite(phys.W)) h_W->Fill(phys.W);
                if (std::isfinite(phys.s)) h_s->Fill(phys.s);
                if (std::isfinite(phys.xB)) h_xB->Fill(phys.xB);
                if (std::isfinite(phys.z)) h_z->Fill(phys.z);
                if (std::isfinite(phys.theta)) h_theta->Fill(phys.theta);
                if (std::isfinite(phys.phi)) h_phi->Fill(phys.phi);
                
                // Store tree data in memory (will be updated in second pass with weights/exclusivity)
                TreeEntry entry;
                entry.mpi0_all = mgg;  // Invariant mass
                entry.mmiss_all = mm_p;
                entry.mmiss_all_corr = mm_p_corr;
                entry.mmiss_all_no_mom_offset = mm_p_no_offset;
                entry.mmiss_all_corr_no_mom_offset = mm_p_corr_no_mom_offset;
                entry.pi0_weight = 0.0;  // Will be filled in second pass
                entry.is_exclusive = 0;  // Will be set in second pass
                entry.is_exclusive_ellipse = 0;
                entry.is_exclusive_mcd = 0;
                entry.is_weighted = 0;   // Will be set in second pass
                entry.Q2 = phys.Q2;
                entry.W = phys.W;
                entry.t = phys.t;
                entry.tmin = phys.tmin;
                entry.pt = phys.pt;
                entry.theta = phys.theta;
                entry.phi = phys.phi;
                entry.s = phys.s;
                entry.xB = phys.xB;
                entry.z = phys.z;
                entry.nclust_selected = good_idx.size();
                entry.cluster_x_1 = clusX[sel_i];
                entry.cluster_y_1 = clusY[sel_i];
                entry.cluster_e_1 = clusE[sel_i];
                entry.cluster_x_2 = clusX[sel_j];
                entry.cluster_y_2 = clusY[sel_j];
                entry.cluster_e_2 = clusE[sel_j];
                // HMS tracking variables
                entry.delta = hdelta;
                entry.xptar = HgtrTh;
                entry.yptar = HgtrPh;
                entry.ytar = HgtrY;
                entry.xtar = HgtrX;
                entry.xfp = hxfp;
                entry.yfp = hyfp;
                entry.xpfp = hxpfp;
                entry.ypfp = hypfp;
                entry.event_id = ev;
                
                treeData[ev] = entry;
            } catch (const std::exception& e) {
                // Event-level error: skip this event and continue
                if (n_event_errors == 0) {  // Log only once per run to avoid spam
                    logmsg(WARN, Form("Event processing error at run %d, event %lld: %s", run, ev, e.what()));
                }
                ++n_event_errors;
                continue;
            } catch (...) {
                // Unknown error: skip and continue
                if (n_event_errors == 0) {
                    logmsg(WARN, Form("Unknown error during event processing at run %d, event %lld", run, ev));
                }
                ++n_event_errors;
                continue;
            }
        } // end for loop
        } catch (...) {
        logmsg(WARN, Form("Critical error in event loop for run %d; processing may be incomplete", run));
        }
        
        if (n_event_errors > 0) {
        logmsg(INFO, Form("Run %d: skipped %lld events due to processing errors", run, n_event_errors));
        }
        
        cout << endl;
        
        // Print detailed cut flow statistics
        cout << "\n===== Run " << run << " Cut Flow =====\n";
        cout << " Total entries:               " << n_total << "\n";
        cout << " Pass HMS cuts:               " << n_pass_hms << " (" << (100.0*n_pass_hms/n_total) << "%)\n";
        cout << " Have >= 2 clusters:          " << n_ge2_hms << " (" << (100.0*n_ge2_hms/n_pass_hms) << "% of HMS)\n";
        cout << " Pass spatial/energy/timing:  " << n_selected_for_analysis << " (" << (100.0*n_selected_for_analysis/n_ge2_hms) << "% of >=2 clusters)\n";
        cout << "=====================================\n\n";


        // -------------------------
        // Summaries & background estimate
        // -------------------------
        // Use custom shifted windows (all except coin shifted by +3ns)
        nps::CoincidenceBGResult bg = nps::estimate_coincidence_background_default(
            h_t1_t2, coin_win, diag_windows, side_windows, 
            full1_t1, full1_t2, full2_t2, full2_t1
        );

        // Data-driven accidental subtraction (returns bg-subtracted histogram)
        h_coin_bgsub = nps::make_and_subtract_accidentals_data_driven(
        h_m_pi0_coin, bg,
        h_mgg_full1, h_mgg_full2,
        h_mgg_hor, h_mgg_ver, h_mgg_diag,
        h_m_pi0_acc,
        h_t1_t2,
        diag_windows, side_windows,
        outPlotDir, run
        );

        // Fit combinatorial BG and subtract (user helper)
        TH1D* h_final = nullptr;
        nps::BGSubtractionResult res;
        res.h_final = nullptr; res.chi2_ndf = 0.0; res.mu_MeV = 0.0; res.sigma_MeV = 0.0; res.signal_counts = 0.0;

        if (h_coin_bgsub) {
        res = nps::FitCombinatorialBGAndSubtract(
            h_coin_bgsub,
            outPlotDir.Data(),
            run,
            4,        // polynomial order
            0.01, 0.11,  // left bg windows
            0.15, 0.40,  // right bg windows
            true      // draw diagnostics inside helper
        );
        h_final = res.h_final; // assigned as in original; caller owns h_final
        }

        // -----------------------------------------------
        // Compute π⁰ bin-by-bin weights for statistical background subtraction
        // -----------------------------------------------
        // weight_i = N_final(i) / N_all(i)
        // This accounts for π⁰ purity in each invariant mass bin
        bool weights_computed = false;
        if (h_mpi0_all && h_final) {
            int nbins = h_mpi0_all->GetNbinsX();
            for (int bin = 1; bin <= nbins; ++bin) {
                double n_all = h_mpi0_all->GetBinContent(bin);
                double n_final = h_final->GetBinContent(bin);
                double weight = (n_all > 0.0) ? (n_final / n_all) : 0.0;
                h_pi0_weight->SetBinContent(bin, weight);
            }
            weights_computed = true;
            logmsg(INFO, Form("Run %d: π⁰ weights computed (weight = signal/all)", run));
        } else {
            logmsg(WARN, Form("Run %d: Could not compute π⁰ weights (h_mpi0_all=%p, h_final=%p)", 
                              run, (void*)h_mpi0_all, (void*)h_final));
        }

        // -----------------------------------------------
        // SECOND PASS: Weighted physics distributions and correlations
        // -----------------------------------------------
        if (weights_computed) {
            logmsg(INFO, Form("Run %d: Starting weighted physics distribution pass", run));
            Long64_t n_weighted_events = 0;
            
            try {
                for (Long64_t ev = 0; ev < nentries; ++ev) {
                    if ((ev % print_every) == 0) {
                        print_progress(run, "Weighted pass", ev, nentries);
                    }

                    T->GetEntry(ev);

                    // Load cluster data (same as first pass, with +150 ns timing offset)
                    int nclust = 0;
                    if (clusE_vec && clusX_vec && clusY_vec && clusT_vec) {
                        nclust = (int) std::min<size_t>(clusE_vec->size(), MAX_CLUS);
                        for (int i=0; i<nclust; ++i) {
                            clusE[i] = (*clusE_vec)[i];
                            clusX[i] = (*clusX_vec)[i];
                            clusY[i] = (*clusY_vec)[i];
                            clusT[i] = (*clusT_vec)[i] + 150.0;  // Add 150 ns offset for production data
                        }
                    } else {
                        nclust = nclust_dbl;
                        if (nclust < 0) nclust = 0;
                        if (nclust > MAX_CLUS) nclust = MAX_CLUS;
                        // If using arrays, timing offset may already be applied or stored directly
                    }

                    if (nclust < 2) continue;
                    if (!nps::hms_electron_cuts(edtmtdc, hdelta, HgtrTh, HgtrPh, hcernpeSum, hcaletotnorm, hreactz)) continue;

                    const int n_after = nps::packClusters(clusE, clusX, clusY, clusT, nclust);
                    vector<int> good_idx; good_idx.reserve(8);
                    for (int i = 0; i < n_after; ++i) {
                        if (nps::nps_spatial_energy_cuts(clusE[i], clusX[i], clusY[i], clusT[i], DEFAULT_TIME_WINDOW_WRT_150))
                            good_idx.push_back(i);
                    }

                    if (good_idx.size() < 2) continue;
                    if (good_idx.size() > 4) {
                        sort(good_idx.begin(), good_idx.end(), [&](int a, int b){ return clusE[a] > clusE[b]; });
                        good_idx.resize(4);
                    }

                    int sel_i = -1, sel_j = -1;
                    if (good_idx.size() == 2) { 
                        sel_i = good_idx[0]; sel_j = good_idx[1]; 
                    } else {
                        auto pr = nps::choose_best_pair_closest_pi0(good_idx, clusE, clusX, clusY, clusT, nps::kDefaultZ_NPS_cm, nps::kPi0Mass_GeV, DEFAULT_TIME_THRESH_NS);
                        sel_i = pr.first; sel_j = pr.second;
                        if (sel_i < 0 || sel_j < 0) continue;
                    }

                    const double px_e = HgtrPx, py_e = HgtrPy, pz_e = HgtrPz;
                    const double p_e_mom = sqrt(max(0.0, px_e*px_e + py_e*py_e + pz_e*pz_e));
                    const double Ee = sqrt(max(0.0, p_e_mom*p_e_mom + nps::kElectronMass_GeV*nps::kElectronMass_GeV));

                    const double mgg = nps::invariant_mass_pi0(clusE[sel_i], clusE[sel_j], clusX[sel_i], clusX[sel_j], clusY[sel_i], clusY[sel_j], nps::kDefaultZ_NPS_cm);

                    const double px_e_scaled = px_e * HMS_MOM_OFFSET_SCALE;
                    const double py_e_scaled = py_e * HMS_MOM_OFFSET_SCALE;
                    const double pz_e_scaled = pz_e * HMS_MOM_OFFSET_SCALE;
                    const double p_e_mom_scaled = sqrt(max(0.0, px_e_scaled*px_e_scaled + py_e_scaled*py_e_scaled + pz_e_scaled*pz_e_scaled));
                    const double Ee_scaled = sqrt(max(0.0, p_e_mom_scaled*p_e_mom_scaled + nps::kElectronMass_GeV*nps::kElectronMass_GeV));

                    const double mm_p_no_offset = nps::missing_mass_proton_pi0(beam_energy, Ee, px_e, py_e, pz_e, clusE[sel_i], clusE[sel_j], clusX[sel_i], clusY[sel_i], clusX[sel_j], clusY[sel_j], nps::kDefaultZ_NPS_cm, -17.51);
                    (void)mm_p_no_offset;

                    const double mm_p = nps::missing_mass_proton_pi0(beam_energy, Ee_scaled, px_e_scaled, py_e_scaled, pz_e_scaled, clusE[sel_i], clusE[sel_j], clusX[sel_i], clusY[sel_i], clusX[sel_j], clusY[sel_j], nps::kDefaultZ_NPS_cm, -17.51);

                    const double mm_p_corr = nps::invariant_missing_mass_corrected_avnish_from_detector(mm_p*mm_p, beam_energy, Ee_scaled, px_e_scaled, py_e_scaled, pz_e_scaled, clusE[sel_i], clusE[sel_j], clusX[sel_i], clusY[sel_i], clusX[sel_j], clusY[sel_j], nps::kDefaultZ_NPS_cm, -17.51);

                    bool isExclusive = (mm_p_corr >= 0.8 && mm_p_corr <= 1.1);
                    int bin = h_pi0_weight->FindBin(mgg);
                    double weight = h_pi0_weight->GetBinContent(bin);

                    if (weight > 0.0 && std::isfinite(weight)) {
                        h_mmiss_all_weighted->Fill(mm_p_corr, weight);
                        
                        // Update tree data with weight and exclusivity information
                        if (treeData.find(ev) != treeData.end()) {
                            treeData[ev].pi0_weight = weight;
                            treeData[ev].is_exclusive = isExclusive ? 1 : 0;
                            treeData[ev].is_weighted = 1;
                        }
                        
                        // Compute physics variables once
                        auto phys = nps::compute_physics_vars_from_detector(beam_energy, Ee, px_e, py_e, pz_e, clusE[sel_i], clusX[sel_i], clusY[sel_i], clusE[sel_j], clusX[sel_j], clusY[sel_j], nps::kDefaultZ_NPS_cm, -17.51, false);
                        
                        // Fill weighted distributions (all events, no exclusivity cut)
                        if (std::isfinite(phys.t)) h_t_weighted->Fill(phys.t, weight);
                        if (std::isfinite(phys.tmin)) h_tmin_weighted->Fill(phys.tmin, weight);
                        if (std::isfinite(phys.pt)) h_pt_weighted->Fill(phys.pt, weight);
                        if (std::isfinite(phys.Q2)) h_Q2_weighted->Fill(phys.Q2, weight);
                        if (std::isfinite(phys.W)) h_W_weighted->Fill(phys.W, weight);
                        if (std::isfinite(phys.s)) h_s_weighted->Fill(phys.s, weight);
                        if (std::isfinite(phys.xB)) h_xB_weighted->Fill(phys.xB, weight);
                        if (std::isfinite(phys.z)) h_z_weighted->Fill(phys.z, weight);
                        if (std::isfinite(phys.theta)) h_theta_weighted->Fill(phys.theta, weight);
                        if (std::isfinite(phys.phi)) h_phi_weighted->Fill(phys.phi, weight);
                        
                        // Fill 2D heatmaps on detector (weighted only)
                        double nu = beam_energy - Ee;
                        double t_value = -phys.t;
                        double mx = mm_p_corr;
                        double clustX = clusX[sel_i];
                        double clustY = clusY[sel_i];
                        
                        if (std::isfinite(phys.xB) && std::isfinite(clustX) && std::isfinite(clustY)) {
                            h_xB_xy->Fill(clustX, clustY, phys.xB * weight);
                        }
                        if (std::isfinite(nu) && std::isfinite(clustX) && std::isfinite(clustY)) {
                            h_nu_xy->Fill(clustX, clustY, nu * weight);
                        }
                        if (std::isfinite(mx) && std::isfinite(clustX) && std::isfinite(clustY)) {
                            h_mx_xy->Fill(clustX, clustY, mx * weight);
                        }
                        if (std::isfinite(phys.phi) && std::isfinite(clustX) && std::isfinite(clustY)) {
                            h_phi_xy->Fill(clustX, clustY, phys.phi * weight);
                        }
                        if (std::isfinite(t_value) && std::isfinite(clustX) && std::isfinite(clustY)) {
                            h_t_xy->Fill(clustX, clustY, t_value * weight);
                        }
                        
                        // Fill exclusive+weighted distributions and correlations (signal sample only)
                        if (isExclusive) {
                            if (std::isfinite(phys.t)) h_t_excl_weighted->Fill(phys.t, weight);
                            if (std::isfinite(phys.tmin)) h_tmin_excl_weighted->Fill(phys.tmin, weight);
                            if (std::isfinite(phys.pt)) h_pt_excl_weighted->Fill(phys.pt, weight);
                            if (std::isfinite(phys.Q2)) h_Q2_excl_weighted->Fill(phys.Q2, weight);
                            if (std::isfinite(phys.W)) h_W_excl_weighted->Fill(phys.W, weight);
                            if (std::isfinite(phys.s)) h_s_excl_weighted->Fill(phys.s, weight);
                            if (std::isfinite(phys.xB)) h_xB_excl_weighted->Fill(phys.xB, weight);
                            if (std::isfinite(phys.z)) h_z_excl_weighted->Fill(phys.z, weight);
                            if (std::isfinite(phys.theta)) h_theta_excl_weighted->Fill(phys.theta, weight);
                            if (std::isfinite(phys.phi)) h_phi_excl_weighted->Fill(phys.phi, weight);
                            
                            // Fill 2D correlations (weighted+exclusive only)
                            if (std::isfinite(phys.Q2) && std::isfinite(phys.W)) {
                                h_Q2_vs_W->Fill(phys.W, phys.Q2, weight);
                            }
                            if (std::isfinite(phys.Q2) && std::isfinite(t_value)) {
                                h_Q2_vs_t->Fill(t_value, phys.Q2, weight);
                            }
                            if (std::isfinite(phys.Q2) && std::isfinite(phys.xB)) {
                                h_Q2_vs_xB->Fill(phys.xB, phys.Q2, weight);
                            }
                            if (std::isfinite(phys.W) && std::isfinite(t_value)) {
                                h_W_vs_t->Fill(t_value, phys.W, weight);
                            }
                            if (std::isfinite(phys.W) && std::isfinite(phys.xB)) {
                                h_W_vs_xB->Fill(phys.xB, phys.W, weight);
                            }
                            if (std::isfinite(t_value) && std::isfinite(phys.xB)) {
                                h_t_vs_xB->Fill(phys.xB, t_value, weight);
                            }
                            if (std::isfinite(t_value) && std::isfinite(phys.phi)) {
                                h_t_vs_phi->Fill(phys.phi, t_value, weight);
                            }
                            if (std::isfinite(phys.xB) && std::isfinite(nu)) {
                                h_xB_vs_nu->Fill(nu, phys.xB, weight);
                            }
                            if (std::isfinite(phys.Q2) && std::isfinite(nu)) {
                                h_Q2_vs_nu->Fill(nu, phys.Q2, weight);
                            }
                            if (std::isfinite(phys.W) && std::isfinite(nu)) {
                                h_W_vs_nu->Fill(nu, phys.W, weight);
                            }
                            ++n_weighted_events;
                        }
                    }
                } // end weighted event loop
            } catch (...) {
                logmsg(WARN, Form("Error during weighted pass for run %d", run));
            }
            cout << endl;
            logmsg(INFO, Form("Run %d: Weighted pass complete (%lld exclusive events)", run, n_weighted_events));
        } else {
            logmsg(WARN, Form("Run %d: Skipping weighted pass (weights not computed)", run));
        }

        // -----------------------------------------------
        // 2D mmiss_all:mpi0_all exclusivity flags.
        // These do not remove events; they add branch-level selectors.
        // -----------------------------------------------
        {
            std::vector<nps2d::Point> mass_cut_points;
            std::vector<Long64_t> mass_cut_event_ids;
            mass_cut_points.reserve(treeData.size());
            mass_cut_event_ids.reserve(treeData.size());
            for (const auto& kv : treeData) {
                const auto& entry = kv.second;
                if (entry.pi0_weight <= 0.0 || !std::isfinite(entry.pi0_weight)) continue;
                nps2d::Point p;
                p.id = kv.first;
                p.mpi0 = entry.mpi0_all;
                p.mmiss = entry.mmiss_all;
                p.weight = entry.pi0_weight;
                mass_cut_points.push_back(p);
                mass_cut_event_ids.push_back(kv.first);
            }

            nps2d::Config mass_cut_cfg;
            mass_cut_cfg.output_dir = std::string(outPlotDir.Data());
            mass_cut_cfg.tag = Form("mass_cut_run%d", run);
            mass_cut_cfg.write_debug = true;
            nps2d::Result mass_cut_result = nps2d::evaluate_mass_cuts(mass_cut_points, mass_cut_cfg);
            if (mass_cut_result.params.valid) {
                Long64_t n_ellipse = 0;
                Long64_t n_mcd = 0;
                for (std::size_t i = 0; i < mass_cut_event_ids.size(); ++i) {
                    auto it = treeData.find(mass_cut_event_ids[i]);
                    if (it == treeData.end()) continue;
                    it->second.is_exclusive_ellipse = mass_cut_result.pass_ellipse[i];
                    it->second.is_exclusive_mcd = mass_cut_result.pass_mcd[i];
                    n_ellipse += mass_cut_result.pass_ellipse[i];
                    n_mcd += mass_cut_result.pass_mcd[i];
                }
                logmsg(INFO, Form("Run %d: 2D mass cuts done (ellipse=%lld, MCD=%lld, peak_fraction=%.3f)",
                                  run, n_ellipse, n_mcd, mass_cut_result.params.peak_fraction));
            } else {
                logmsg(WARN, Form("Run %d: 2D mass cuts failed; ellipse/MCD flags remain 0", run));
            }
        }

        // -----------------------------------------------
        // Create overlay canvases for physics variables (all vs weighted vs exclusive)
        // -----------------------------------------------
        TCanvas *c_Q2_overlay = new TCanvas(Form("c_Q2_overlay_run%d", run), "Q^{2} Overlay", 800, 600);
        h_Q2->SetLineColor(kBlack); h_Q2->SetLineWidth(2);
        h_Q2_weighted->SetLineColor(kBlue); h_Q2_weighted->SetLineWidth(2);
        h_Q2_excl_weighted->SetLineColor(kRed); h_Q2_excl_weighted->SetLineWidth(2);
        h_Q2->Draw("HIST");
        h_Q2_weighted->Draw("HIST SAME");
        h_Q2_excl_weighted->Draw("HIST SAME");
        TLegend *leg_Q2 = new TLegend(0.6, 0.6, 0.88, 0.88);
        leg_Q2->SetBorderSize(0); leg_Q2->SetFillColor(0);
        leg_Q2->AddEntry(h_Q2, "All", "l");
        leg_Q2->AddEntry(h_Q2_weighted, "Weighted", "l");
        leg_Q2->AddEntry(h_Q2_excl_weighted, "Weighted+Exclusive", "l");
        leg_Q2->Draw();
        c_Q2_overlay->SaveAs(Form("%s/Q2_overlay_run%d.png", outPlotDir.Data(), run));

        TCanvas *c_W_overlay = new TCanvas(Form("c_W_overlay_run%d", run), "W Overlay", 800, 600);
        h_W->SetLineColor(kBlack); h_W->SetLineWidth(2);
        h_W_weighted->SetLineColor(kBlue); h_W_weighted->SetLineWidth(2);
        h_W_excl_weighted->SetLineColor(kRed); h_W_excl_weighted->SetLineWidth(2);
        h_W->Draw("HIST");
        h_W_weighted->Draw("HIST SAME");
        h_W_excl_weighted->Draw("HIST SAME");
        TLegend *leg_W = new TLegend(0.6, 0.6, 0.88, 0.88);
        leg_W->SetBorderSize(0); leg_W->SetFillColor(0);
        leg_W->AddEntry(h_W, "All", "l");
        leg_W->AddEntry(h_W_weighted, "Weighted", "l");
        leg_W->AddEntry(h_W_excl_weighted, "Weighted+Exclusive", "l");
        leg_W->Draw();
        c_W_overlay->SaveAs(Form("%s/W_overlay_run%d.png", outPlotDir.Data(), run));

        // Similar for other physics variables (t, tmin, pt, s, xB, z, theta, phi)
        TCanvas *c_t_overlay = new TCanvas(Form("c_t_overlay_run%d", run), "t Overlay", 800, 600);
        h_t->SetLineColor(kBlack); h_t->SetLineWidth(2);
        h_t_weighted->SetLineColor(kBlue); h_t_weighted->SetLineWidth(2);
        h_t_excl_weighted->SetLineColor(kRed); h_t_excl_weighted->SetLineWidth(2);
        h_t->Draw("HIST");
        h_t_weighted->Draw("HIST SAME");
        h_t_excl_weighted->Draw("HIST SAME");
        TLegend *leg_t = new TLegend(0.15, 0.6, 0.43, 0.88);
        leg_t->SetBorderSize(0); leg_t->SetFillColor(0);
        leg_t->AddEntry(h_t, "All", "l");
        leg_t->AddEntry(h_t_weighted, "Weighted", "l");
        leg_t->AddEntry(h_t_excl_weighted, "Weighted+Exclusive", "l");
        leg_t->Draw();
        c_t_overlay->SaveAs(Form("%s/t_overlay_run%d.png", outPlotDir.Data(), run));

        TCanvas *c_tmin_overlay = new TCanvas(Form("c_tmin_overlay_run%d", run), "t_min Overlay", 800, 600);
        h_tmin->SetLineColor(kBlack); h_tmin->SetLineWidth(2);
        h_tmin_weighted->SetLineColor(kBlue); h_tmin_weighted->SetLineWidth(2);
        h_tmin_excl_weighted->SetLineColor(kRed); h_tmin_excl_weighted->SetLineWidth(2);
        h_tmin->Draw("HIST");
        h_tmin_weighted->Draw("HIST SAME");
        h_tmin_excl_weighted->Draw("HIST SAME");
        TLegend *leg_tmin = new TLegend(0.6, 0.6, 0.88, 0.88);
        leg_tmin->SetBorderSize(0); leg_tmin->SetFillColor(0);
        leg_tmin->AddEntry(h_tmin, "All", "l");
        leg_tmin->AddEntry(h_tmin_weighted, "Weighted", "l");
        leg_tmin->AddEntry(h_tmin_excl_weighted, "Weighted+Exclusive", "l");
        leg_tmin->Draw();
        c_tmin_overlay->SaveAs(Form("%s/tmin_overlay_run%d.png", outPlotDir.Data(), run));

        TCanvas *c_pt_overlay = new TCanvas(Form("c_pt_overlay_run%d", run), "p_T Overlay", 800, 600);
        h_pt->SetLineColor(kBlack); h_pt->SetLineWidth(2);
        h_pt_weighted->SetLineColor(kBlue); h_pt_weighted->SetLineWidth(2);
        h_pt_excl_weighted->SetLineColor(kRed); h_pt_excl_weighted->SetLineWidth(2);
        h_pt->Draw("HIST");
        h_pt_weighted->Draw("HIST SAME");
        h_pt_excl_weighted->Draw("HIST SAME");
        TLegend *leg_pt = new TLegend(0.6, 0.6, 0.88, 0.88);
        leg_pt->SetBorderSize(0); leg_pt->SetFillColor(0);
        leg_pt->AddEntry(h_pt, "All", "l");
        leg_pt->AddEntry(h_pt_weighted, "Weighted", "l");
        leg_pt->AddEntry(h_pt_excl_weighted, "Weighted+Exclusive", "l");
        leg_pt->Draw();
        c_pt_overlay->SaveAs(Form("%s/pt_overlay_run%d.png", outPlotDir.Data(), run));

        TCanvas *c_s_overlay = new TCanvas(Form("c_s_overlay_run%d", run), "s Overlay", 800, 600);
        h_s->SetLineColor(kBlack); h_s->SetLineWidth(2);
        h_s_weighted->SetLineColor(kBlue); h_s_weighted->SetLineWidth(2);
        h_s_excl_weighted->SetLineColor(kRed); h_s_excl_weighted->SetLineWidth(2);
        h_s->Draw("HIST");
        h_s_weighted->Draw("HIST SAME");
        h_s_excl_weighted->Draw("HIST SAME");
        TLegend *leg_s = new TLegend(0.6, 0.6, 0.88, 0.88);
        leg_s->SetBorderSize(0); leg_s->SetFillColor(0);
        leg_s->AddEntry(h_s, "All", "l");
        leg_s->AddEntry(h_s_weighted, "Weighted", "l");
        leg_s->AddEntry(h_s_excl_weighted, "Weighted+Exclusive", "l");
        leg_s->Draw();
        c_s_overlay->SaveAs(Form("%s/s_overlay_run%d.png", outPlotDir.Data(), run));

        TCanvas *c_xB_overlay = new TCanvas(Form("c_xB_overlay_run%d", run), "x_B Overlay", 800, 600);
        h_xB->SetLineColor(kBlack); h_xB->SetLineWidth(2);
        h_xB_weighted->SetLineColor(kBlue); h_xB_weighted->SetLineWidth(2);
        h_xB_excl_weighted->SetLineColor(kRed); h_xB_excl_weighted->SetLineWidth(2);
        h_xB->Draw("HIST");
        h_xB_weighted->Draw("HIST SAME");
        h_xB_excl_weighted->Draw("HIST SAME");
        TLegend *leg_xB = new TLegend(0.6, 0.6, 0.88, 0.88);
        leg_xB->SetBorderSize(0); leg_xB->SetFillColor(0);
        leg_xB->AddEntry(h_xB, "All", "l");
        leg_xB->AddEntry(h_xB_weighted, "Weighted", "l");
        leg_xB->AddEntry(h_xB_excl_weighted, "Weighted+Exclusive", "l");
        leg_xB->Draw();
        c_xB_overlay->SaveAs(Form("%s/xB_overlay_run%d.png", outPlotDir.Data(), run));

        TCanvas *c_z_overlay = new TCanvas(Form("c_z_overlay_run%d", run), "z Overlay", 800, 600);
        h_z->SetLineColor(kBlack); h_z->SetLineWidth(2);
        h_z_weighted->SetLineColor(kBlue); h_z_weighted->SetLineWidth(2);
        h_z_excl_weighted->SetLineColor(kRed); h_z_excl_weighted->SetLineWidth(2);
        h_z->Draw("HIST");
        h_z_weighted->Draw("HIST SAME");
        h_z_excl_weighted->Draw("HIST SAME");
        TLegend *leg_z = new TLegend(0.6, 0.6, 0.88, 0.88);
        leg_z->SetBorderSize(0); leg_z->SetFillColor(0);
        leg_z->AddEntry(h_z, "All", "l");
        leg_z->AddEntry(h_z_weighted, "Weighted", "l");
        leg_z->AddEntry(h_z_excl_weighted, "Weighted+Exclusive", "l");
        leg_z->Draw();
        c_z_overlay->SaveAs(Form("%s/z_overlay_run%d.png", outPlotDir.Data(), run));

        TCanvas *c_theta_overlay = new TCanvas(Form("c_theta_overlay_run%d", run), "#theta Overlay", 800, 600);
        h_theta->SetLineColor(kBlack); h_theta->SetLineWidth(2);
        h_theta_weighted->SetLineColor(kBlue); h_theta_weighted->SetLineWidth(2);
        h_theta_excl_weighted->SetLineColor(kRed); h_theta_excl_weighted->SetLineWidth(2);
        h_theta->Draw("HIST");
        h_theta_weighted->Draw("HIST SAME");
        h_theta_excl_weighted->Draw("HIST SAME");
        TLegend *leg_theta = new TLegend(0.6, 0.6, 0.88, 0.88);
        leg_theta->SetBorderSize(0); leg_theta->SetFillColor(0);
        leg_theta->AddEntry(h_theta, "All", "l");
        leg_theta->AddEntry(h_theta_weighted, "Weighted", "l");
        leg_theta->AddEntry(h_theta_excl_weighted, "Weighted+Exclusive", "l");
        leg_theta->Draw();
        c_theta_overlay->SaveAs(Form("%s/theta_overlay_run%d.png", outPlotDir.Data(), run));

        TCanvas *c_phi_overlay = new TCanvas(Form("c_phi_overlay_run%d", run), "#phi Overlay", 800, 600);
        h_phi->SetLineColor(kBlack); h_phi->SetLineWidth(2);
        h_phi_weighted->SetLineColor(kBlue); h_phi_weighted->SetLineWidth(2);
        h_phi_excl_weighted->SetLineColor(kRed); h_phi_excl_weighted->SetLineWidth(2);
        h_phi->Draw("HIST");
        h_phi_weighted->Draw("HIST SAME");
        h_phi_excl_weighted->Draw("HIST SAME");
        TLegend *leg_phi = new TLegend(0.6, 0.6, 0.88, 0.88);
        leg_phi->SetBorderSize(0); leg_phi->SetFillColor(0);
        leg_phi->AddEntry(h_phi, "All", "l");
        leg_phi->AddEntry(h_phi_weighted, "Weighted", "l");
        leg_phi->AddEntry(h_phi_excl_weighted, "Weighted+Exclusive", "l");
        leg_phi->Draw();
        c_phi_overlay->SaveAs(Form("%s/phi_overlay_run%d.png", outPlotDir.Data(), run));

        // Missing mass overlay with weighted (NO FITTING - just overlay display)
        // IMPORTANT: Create this BEFORE fitting h_mmiss_all_weighted so fit curves don't appear
        TCanvas *c_mmiss_overlay = new TCanvas(Form("c_mmiss_overlay_run%d", run), "Missing Mass Overlay", 800, 600);
        c_mmiss_overlay->SetName(Form("c_mmiss_overlay_run%d", run));
        c_mmiss_overlay->SetTitle(Form("Missing Mass Overlay run %d", run));
        h_mmiss_all->SetLineColor(kBlack); h_mmiss_all->SetLineWidth(2);
        h_mmiss_all_corr->SetLineColor(kBlue); h_mmiss_all_corr->SetLineWidth(2);
        h_mmiss_all_weighted->SetLineColor(kGreen+2); h_mmiss_all_weighted->SetLineWidth(2);
        h_mmiss_all->Draw("HIST");
        h_mmiss_all_corr->Draw("HIST SAME");
        h_mmiss_all_weighted->Draw("HISTNOFUNC SAME");  // HISTNOFUNC prevents fitted functions from drawing
        TLegend *leg_mmiss = new TLegend(0.55, 0.6, 0.88, 0.88);
        leg_mmiss->SetBorderSize(0);
        leg_mmiss->SetFillColor(0);
        leg_mmiss->AddEntry(h_mmiss_all, "All candidates", "l");
        leg_mmiss->AddEntry(h_mmiss_all_corr, "Corrected Mx", "l");
        leg_mmiss->AddEntry(h_mmiss_all_weighted, "Weighted (all)", "l");
        leg_mmiss->Draw();
        c_mmiss_overlay->Modified();
        c_mmiss_overlay->Update();
        c_mmiss_overlay->SaveAs(Form("%s/mmiss_overlay_run%d.png", outPlotDir.Data(), run));
        // NOTE: Will redraw this canvas after fitting to ensure no fit functions appear in the ROOT file

        // -------------------------
        // Fit h_mmiss_all_weighted with combined two-Gaussian fit
        // -------------------------
        // Combined fit: fit the sum of two Gaussians in the full range
        // This ensures physical consistency - fit1 and fit2 are components of the total
        TF1 *fit_combined = new TF1(Form("fit_combined_run%d", run), 
            "gaus(0) + gaus(3)", 0.5, 1.8);
        
        // Set initial parameters for the combined fit
        // Component 1: peak around 0.938 GeV
        fit_combined->SetParameter(0, h_mmiss_all_weighted->GetMaximum());     // amplitude
        fit_combined->SetParameter(1, 0.938);                                  // mean
        fit_combined->SetParameter(2, 0.05);                                   // sigma
        
        // Component 2: peak around 1.2 GeV
        fit_combined->SetParameter(3, h_mmiss_all_weighted->GetMaximum()*0.5); // amplitude
        fit_combined->SetParameter(4, 1.2);                                    // mean
        fit_combined->SetParameter(5, 0.08);                                   // sigma
        
        // Fit the combined function
        h_mmiss_all_weighted->Fit(fit_combined, "RQ", "", 0.85, 1.27);
        fit_combined->SetLineColor(kGreen+2);
        fit_combined->SetLineStyle(kSolid);
        fit_combined->SetLineWidth(3);
        
        // Extract the individual components from the combined fit
        // fit1: first Gaussian componentls
        
        TF1 *fit1 = new TF1(Form("fit1_run%d", run), "gaus", 0.5, 1.8);
        fit1->SetParameters(
            fit_combined->GetParameter(0),
            fit_combined->GetParameter(1),
            fit_combined->GetParameter(2)
        );
        fit1->SetLineColor(kBlue);
        fit1->SetLineStyle(kDashed);
        fit1->SetLineWidth(2);
        
        // fit2: second Gaussian component
        TF1 *fit2 = new TF1(Form("fit2_run%d", run), "gaus", 0.5, 1.8);
        fit2->SetParameters(
            fit_combined->GetParameter(3),
            fit_combined->GetParameter(4),
            fit_combined->GetParameter(5)
        );
        fit2->SetLineColor(kMagenta);
        fit2->SetLineStyle(kDashed);
        fit2->SetLineWidth(2);

        // Create canvas with all three fits
        TCanvas *c_mmiss_fit = new TCanvas(Form("c_mmiss_fit_run%d", run), "Missing Mass Weighted Fit", 1000, 700);
        c_mmiss_fit->SetName(Form("c_mmiss_fit_run%d", run));
        c_mmiss_fit->SetTitle(Form("Missing Mass Weighted - Fits run %d", run));
        c_mmiss_fit->SetLeftMargin(0.12);
        c_mmiss_fit->SetRightMargin(0.05);
        c_mmiss_fit->SetBottomMargin(0.12);
        c_mmiss_fit->SetTopMargin(0.08);
        
        h_mmiss_all_weighted->SetLineColor(kBlack);
        h_mmiss_all_weighted->SetLineWidth(2);
        h_mmiss_all_weighted->GetXaxis()->SetLabelSize(0.045);
        h_mmiss_all_weighted->GetYaxis()->SetLabelSize(0.045);
        h_mmiss_all_weighted->GetXaxis()->SetTitleSize(0.050);
        h_mmiss_all_weighted->GetYaxis()->SetTitleSize(0.050);
        h_mmiss_all_weighted->GetXaxis()->SetTitle("Missing Mass (GeV)");
        h_mmiss_all_weighted->GetYaxis()->SetTitle("Counts");
        h_mmiss_all_weighted->Draw("HIST");
        
        fit_combined->Draw("SAME");
        fit1->Draw("SAME");
        fit2->Draw("SAME");
        
        TLegend *leg_fit = new TLegend(0.55, 0.55, 0.92, 0.88);
        leg_fit->SetBorderSize(1);
        leg_fit->SetFillColor(10);
        leg_fit->SetFillStyle(1001);
        leg_fit->SetTextFont(42);
        leg_fit->SetTextSize(0.035);
        leg_fit->AddEntry(h_mmiss_all_weighted, "Data", "l");
        leg_fit->AddEntry(fit_combined, Form("Combined (chi2/ndf = %.2f)", fit_combined->GetChisquare()/fit_combined->GetNDF()), "l");
        leg_fit->AddEntry(fit1, "Gaussian 1 [0.9-1.0]", "l");
        leg_fit->AddEntry(fit2, "Gaussian 2 [1.1-1.3]", "l");
        leg_fit->Draw();

        TLatex *txt_fit = new TLatex();
        txt_fit->SetNDC();
        txt_fit->SetTextFont(42);
        txt_fit->SetTextSize(0.032);
        txt_fit->DrawLatex(0.15, 0.52, Form("Run %d", run));
        txt_fit->DrawLatex(0.15, 0.47, Form("Entries: %.0f", h_mmiss_all_weighted->GetEntries()));
        txt_fit->DrawLatex(0.15, 0.41, Form("G1: #mu=%.4f GeV, #sigma=%.4f GeV", fit1->GetParameter(1), fit1->GetParameter(2)));
        txt_fit->DrawLatex(0.15, 0.35, Form("G2: #mu=%.4f GeV, #sigma=%.4f GeV", fit2->GetParameter(1), fit2->GetParameter(2)));
        
        c_mmiss_fit->Modified();
        c_mmiss_fit->Update();
        c_mmiss_fit->SaveAs(Form("%s/mmiss_weighted_fit_run%d.png", outPlotDir.Data(), run));

        // Keep histogram black for fit canvas in ROOT file (don't change color)
        // The fit canvas will be written to file with the histogram in black
        
        // Clear all fit functions from h_mmiss_all_weighted so they don't appear in overlay
        h_mmiss_all_weighted->GetListOfFunctions()->Delete();
        
        // Restore histogram colors for overlay canvas
        h_mmiss_all->SetLineColor(kBlack); h_mmiss_all->SetLineWidth(2);
        h_mmiss_all_corr->SetLineColor(kBlue); h_mmiss_all_corr->SetLineWidth(2);
        h_mmiss_all_weighted->SetLineColor(kGreen+2); h_mmiss_all_weighted->SetLineWidth(2);
        
        // Redraw the overlay canvas after fitting (so it has no fit curves in the ROOT file)
        c_mmiss_overlay->Clear();
        h_mmiss_all->Draw("HIST");
        h_mmiss_all_corr->Draw("HIST SAME");
        h_mmiss_all_weighted->Draw("HIST SAME");
        TLegend *leg_mmiss_redraw = new TLegend(0.55, 0.6, 0.88, 0.88);
        leg_mmiss_redraw->SetBorderSize(0);
        leg_mmiss_redraw->SetFillColor(0);
        leg_mmiss_redraw->AddEntry(h_mmiss_all, "All candidates", "l");
        leg_mmiss_redraw->AddEntry(h_mmiss_all_corr, "Corrected Mx", "l");
        leg_mmiss_redraw->AddEntry(h_mmiss_all_weighted, "Weighted (all)", "l");
        leg_mmiss_redraw->Draw();
        c_mmiss_overlay->Modified();
        c_mmiss_overlay->Update();

        // -------------------------
        // Open per-run ROOT file (and ensure we write everything there)
        // -------------------------
        TString outf = Form("%s/diagnostics_run%d.root", outPlotDir.Data(), run);
        TFile *fout = TFile::Open(outf, "RECREATE");
        if (!fout || fout->IsZombie()) {
        logmsg(WARN, Form("Could not create ROOT output %s", outf.Data()));
        } else {
        fout->cd();

        // Create and fill event-level physics tree from treeData map
        TTree *treeOut = new TTree("physics", "Event-level physics data with weights and exclusivity flags");
        
        // Declare branch variables fresh for this tree
        Int_t event_id = 0;
        Double_t mpi0_all = 0, mmiss_all = 0, mmiss_all_corr = 0, mmiss_all_no_mom_offset = 0, mmiss_all_corr_no_mom_offset = 0, pi0_weight = 0;
        Int_t is_exclusive = 0, is_exclusive_ellipse = 0, is_exclusive_mcd = 0, is_weighted = 0;
        Double_t Q2 = 0, W = 0, t = 0, tmin = 0, pt = 0;
        Double_t theta = 0, phi = 0, s = 0, xB = 0, z = 0;
        Int_t nclust_selected = 0;
        Double_t cluster_x_1 = 0, cluster_y_1 = 0, cluster_e_1 = 0;
        Double_t cluster_x_2 = 0, cluster_y_2 = 0, cluster_e_2 = 0;
        // HMS tracking variables
        Double_t delta = 0, xptar = 0, yptar = 0, xtar = 0, ytar = 0;
        Double_t xfp = 0, yfp = 0, xpfp = 0, ypfp = 0;
        // Seperating the exclusive peak
        Double_t exclusive_sep=0;
        
        // Create branches
        treeOut->Branch("event_id", &event_id, "event_id/I");
        treeOut->Branch("mpi0_all", &mpi0_all, "mpi0_all/D");
        treeOut->Branch("mmiss_all", &mmiss_all, "mmiss_all/D");
        treeOut->Branch("mmiss_all_corr", &mmiss_all_corr, "mmiss_all_corr/D");
        treeOut->Branch("mmiss_all_no_mom_offset", &mmiss_all_no_mom_offset, "mmiss_all_no_mom_offset/D");
        treeOut->Branch("mmiss_all_corr_no_mom_offset", &mmiss_all_corr_no_mom_offset, "mmiss_all_corr_no_mom_offset/D");
        treeOut->Branch("pi0_weight", &pi0_weight, "pi0_weight/D");
        treeOut->Branch("is_exclusive", &is_exclusive, "is_exclusive/I");
        treeOut->Branch("is_exclusive_ellipse", &is_exclusive_ellipse, "is_exclusive_ellipse/I");
        treeOut->Branch("is_exclusive_mcd", &is_exclusive_mcd, "is_exclusive_mcd/I");
        treeOut->Branch("is_weighted", &is_weighted, "is_weighted/I");
        
        treeOut->Branch("Q2", &Q2, "Q2/D");
        treeOut->Branch("W", &W, "W/D");
        treeOut->Branch("t", &t, "t/D");
        treeOut->Branch("tmin", &tmin, "tmin/D");
        treeOut->Branch("pt", &pt, "pt/D");
        treeOut->Branch("theta", &theta, "theta/D");
        treeOut->Branch("phi", &phi, "phi/D");
        treeOut->Branch("s", &s, "s/D");
        treeOut->Branch("xB", &xB, "xB/D");
        treeOut->Branch("z", &z, "z/D");
        
        treeOut->Branch("nclust_selected", &nclust_selected, "nclust_selected/I");
        treeOut->Branch("cluster_x_1", &cluster_x_1, "cluster_x_1/D");
        treeOut->Branch("cluster_y_1", &cluster_y_1, "cluster_y_1/D");
        treeOut->Branch("cluster_e_1", &cluster_e_1, "cluster_e_1/D");
        treeOut->Branch("cluster_x_2", &cluster_x_2, "cluster_x_2/D");
        treeOut->Branch("cluster_y_2", &cluster_y_2, "cluster_y_2/D");
        treeOut->Branch("cluster_e_2", &cluster_e_2, "cluster_e_2/D");
        
        // HMS tracking branches
        treeOut->Branch("delta", &delta, "delta/D");
        treeOut->Branch("xptar", &xptar, "xptar/D");
        treeOut->Branch("yptar", &yptar, "yptar/D");
        treeOut->Branch("xtar", &xtar, "xtar/D");
        treeOut->Branch("ytar", &ytar, "ytar/D");
        treeOut->Branch("xfp", &xfp, "xfp/D");
        treeOut->Branch("yfp", &yfp, "yfp/D");
        treeOut->Branch("xpfp", &xpfp, "xpfp/D");
        treeOut->Branch("ypfp", &ypfp, "ypfp/D");

        treeOut->Branch("exclusive_sep", &exclusive_sep, "exclusive_sep/D");
        
        // Fill tree from treeData map
        for (const auto& pair : treeData) {
            event_id = pair.first;
            const auto& entry = pair.second;
            mpi0_all = entry.mpi0_all;
            mmiss_all = entry.mmiss_all;
            mmiss_all_corr = entry.mmiss_all_corr;
            mmiss_all_no_mom_offset = entry.mmiss_all_no_mom_offset;
            mmiss_all_corr_no_mom_offset = entry.mmiss_all_corr_no_mom_offset;
            pi0_weight = entry.pi0_weight;
            is_exclusive = entry.is_exclusive;
            is_exclusive_ellipse = entry.is_exclusive_ellipse;
            is_exclusive_mcd = entry.is_exclusive_mcd;
            is_weighted = entry.is_weighted;
            Q2 = entry.Q2;
            W = entry.W;
            t = entry.t;
            tmin = entry.tmin;
            pt = entry.pt;
            theta = entry.theta;
            phi = entry.phi;
            s = entry.s;
            xB = entry.xB;
            z = entry.z;
            nclust_selected = entry.nclust_selected;
            cluster_x_1 = entry.cluster_x_1;
            cluster_y_1 = entry.cluster_y_1;
            cluster_e_1 = entry.cluster_e_1;
            cluster_x_2 = entry.cluster_x_2;
            cluster_y_2 = entry.cluster_y_2;
            cluster_e_2 = entry.cluster_e_2;
            // HMS tracking variables
            delta = entry.delta;
            xptar = entry.xptar;
            yptar = entry.yptar;
            xtar = entry.xtar;
            ytar = entry.ytar;
            xfp = entry.xfp;
            yfp = entry.yfp;
            xpfp = entry.xpfp;
            ypfp = entry.ypfp;
            exclusive_sep = 1.0 / ( pow(std::fabs(entry.mmiss_all_corr - 0.938), 0.5) + 0.01 );  // separation of event from exclusive peak in GeV (for potential future cuts)
            treeOut->Fill();
        }
        
        // Write and cleanup tree
        treeOut->Write();
        safe_delete(treeOut);
        treeData.clear();

        // write histograms (explicit list)
        h_nclusters->Write();
        h_clustE->Write(); h_clustT->Write(); h_clustE_vs_T->Write();
        h_clustXY->Write(); h_clustE_sum->Write(); h_opening_angle->Write(); h_photon_Eratio->Write();
        h_mpi0_all->Write(); h_mpi0_2->Write(); h_mpi0_3->Write(); h_mpi0_4->Write();
        h_mmiss_2->Write(); h_mmiss_3->Write(); h_mmiss_4->Write(); h_mmiss_dvcs->Write();
        h_mmiss_all->Write(); h_mmiss_all_corr->Write();
        // save the canvas under a unique name
        c_mmiss_overlay->Write(Form("c_mmiss_overlay_run%d", run));
        // Set histogram to black before writing fit canvas to file
        h_mmiss_all_weighted->SetLineColor(kBlack);
        c_mmiss_fit->Write(Form("c_mmiss_fit_run%d", run));
        h_t1_t2->Write("h_t1_t2", TObject::kOverwrite); h_t1_proj->Write(); h_t2_proj->Write();
        h_m_pi0_coin->Write(); h_m_pi0_acc->Write();
        if (h_coin_bgsub) h_coin_bgsub->Write("h_pi0_coin_bgsub", TObject::kOverwrite);
        if (h_final) h_final->Write("h_pi0_final", TObject::kOverwrite);

        for (auto *h: h_mgg_diag) if (h) h->Write();
        for (auto *h: h_mgg_hor) if (h) h->Write();
        for (auto *h: h_mgg_ver) if (h) h->Write();
        h_mgg_full1->Write(); h_mgg_full2->Write();
        h_mmiss_vs_mgg->Write();
        h_mmiss_vs_mgg_corr->Write();
        h_mmiss_vs_t1->Write(); h_mmiss_vs_t2->Write();

        h_t->Write();
        h_tmin->Write();
        h_pt->Write();
        h_Q2->Write();
        h_W->Write();
        h_s->Write();
        h_xB->Write();
        h_z->Write();
        h_theta->Write();
        h_phi->Write();

        // Write weighted and exclusive histograms
        h_pi0_weight->Write();
        h_Q2_weighted->Write(); h_W_weighted->Write(); h_t_weighted->Write();
        h_tmin_weighted->Write(); h_pt_weighted->Write(); h_theta_weighted->Write();
        h_phi_weighted->Write(); h_s_weighted->Write(); h_xB_weighted->Write(); h_z_weighted->Write();
        h_Q2_excl_weighted->Write(); h_W_excl_weighted->Write(); h_t_excl_weighted->Write();
        h_tmin_excl_weighted->Write(); h_pt_excl_weighted->Write(); h_theta_excl_weighted->Write();
        h_phi_excl_weighted->Write(); h_s_excl_weighted->Write(); h_xB_excl_weighted->Write(); h_z_excl_weighted->Write();
        h_mmiss_all_weighted->Write();

        // Write 2D heatmaps (physics on detector)
        h_xB_xy->Write(); h_nu_xy->Write(); h_mx_xy->Write(); h_phi_xy->Write(); h_t_xy->Write();

        // Write 2D correlations
        h_Q2_vs_W->Write(); h_Q2_vs_t->Write(); h_Q2_vs_xB->Write();
        h_W_vs_t->Write(); h_W_vs_xB->Write(); h_t_vs_xB->Write();
        h_t_vs_phi->Write(); h_xB_vs_nu->Write(); h_Q2_vs_nu->Write(); h_W_vs_nu->Write();

        // Write overlay canvases
        c_Q2_overlay->Write(Form("c_Q2_overlay_run%d", run));
        c_W_overlay->Write(Form("c_W_overlay_run%d", run));
        c_t_overlay->Write(Form("c_t_overlay_run%d", run));
        c_tmin_overlay->Write(Form("c_tmin_overlay_run%d", run));
        c_pt_overlay->Write(Form("c_pt_overlay_run%d", run));
        c_s_overlay->Write(Form("c_s_overlay_run%d", run));
        c_xB_overlay->Write(Form("c_xB_overlay_run%d", run));
        c_z_overlay->Write(Form("c_z_overlay_run%d", run));
        c_theta_overlay->Write(Form("c_theta_overlay_run%d", run));
        c_phi_overlay->Write(Form("c_phi_overlay_run%d", run));

        // write background and fit scalars (TParameter)
        TParameter<double>("coin_raw", bg.n_coin_raw).Write();
        TParameter<double>("accidental_est", bg.n_accidentals).Write();
        TParameter<double>("accidental_err", bg.n_accidentals_err).Write();
        TParameter<double>("coin_area_ns2", bg.area_coin).Write();

        // write fit summary items returned in res as TParameters (if available)
        TParameter<double>("comb_bg_chi2ndf", res.chi2_ndf).Write();
        TParameter<double>("pi0_mu_MeV", res.mu_MeV).Write();
        TParameter<double>("pi0_sigma_MeV", res.sigma_MeV).Write();
        TParameter<double>("pi0_signal_counts", res.signal_counts).Write();

        // HMS tracking eff finalize/write (keeps your call)
        try {
            auto [hms_track_eff_val, hms_track_eff_unc] = hmstrackEff.finalizeAndWrite(fout, outPlotDir_in);
            std::cout << Form("Run %d HMS tracking eff = %.4f ± %.4f\n", run, hms_track_eff_val, hms_track_eff_unc);
        } catch (...) {
            logmsg(WARN, "hmstrackEff.finalizeAndWrite threw an exception; skipping HMS eff write for this run.");
        }

        fout->Write();
        fout->Close();
        safe_delete(fout);
        logmsg(INFO, Form("Wrote per-run ROOT diagnostics to %s", outf.Data()));
        }

        // -------------------------
        // Draw and save canvases (PNG) for quick look (use unique canvas names)
        // -------------------------
        gStyle->SetOptStat(0);

        // 1) t1 vs t2 canvas with boxes
        TCanvas *c_t12 = new TCanvas(name("c_t12"), "t1 vs t2", 1000,900);
        c_t12->SetName(name("c_t12"));
        h_t1_t2->Draw("COLZ");
        gPad->SetRightMargin(0.15);
        gPad->Update();

        // Draw boxes (coin, diag, side, full) and annotate integrals
        TBox *b_coin = new TBox(coin_win.first, coin_win.first, coin_win.second, coin_win.second);
        b_coin->SetLineColor(kRed); b_coin->SetLineWidth(2); b_coin->SetFillStyle(0); b_coin->Draw("same");

        vector<TBox*> diag_boxes; for (auto &w : diag_windows) { TBox *b=new TBox(w.first,w.first,w.second,w.second); b->SetLineColor(kMagenta); b->SetLineWidth(2); b->SetLineStyle(2); b->SetFillStyle(0); b->Draw("same"); diag_boxes.push_back(b); }
        vector<TBox*> hor_boxes; for (auto &w : side_windows) { TBox *b=new TBox(w.first, coin_win.first, w.second, coin_win.second); b->SetLineColor(kBlue); b->SetLineStyle(3); b->SetFillStyle(0); b->Draw("same"); hor_boxes.push_back(b); }
        vector<TBox*> ver_boxes; for (auto &w : side_windows) { TBox *b=new TBox(coin_win.first, w.first, coin_win.second, w.second); b->SetLineColor(kGreen+2); b->SetLineStyle(3); b->SetFillStyle(0); b->Draw("same"); ver_boxes.push_back(b); }
        TBox *b_full1 = new TBox(full1_t2.first, full1_t1.first, full1_t2.second, full1_t1.second); b_full1->SetLineColor(kOrange+1); b_full1->SetLineWidth(2); b_full1->SetLineStyle(4); b_full1->SetFillStyle(0); b_full1->Draw("same");
        TBox *b_full2 = new TBox(full2_t2.first, full2_t1.first, full2_t2.second, full2_t1.second); b_full2->SetLineColor(kOrange+7); b_full2->SetLineWidth(2); b_full2->SetLineStyle(4); b_full2->SetFillStyle(0); b_full2->Draw("same");

        TLatex *txt = new TLatex();
        txt->SetNDC(); txt->SetTextSize(0.025);
        txt->DrawLatex(0.02, 0.96, Form("Run %d  Coin raw = %.0f  Acc est = %.3f #pm %.3f", run, bg.n_coin_raw, bg.n_accidentals, bg.n_accidentals_err));
        txt->DrawLatex(0.02, 0.92, Form("Coin area (ns^{2}) = %.3f", bg.area_coin));
        c_t12->SaveAs(Form("%s/t1t2_run%d.png", outPlotDir.Data(), run));

        // 2) pi0 mass overlay coin/acc and bgsub (+ final if exists)
        TCanvas *c_pi0 = new TCanvas(name("c_pi0"), "Pi0 inv mass coin vs acc", 900,600);
        c_pi0->SetName(name("c_pi0"));
        h_mpi0_all->SetLineColor(kBlack); h_mpi0_all->SetLineWidth(1);
        h_m_pi0_coin->SetLineColor(kBlue); h_m_pi0_coin->SetLineWidth(2);
        h_mpi0_all->Draw("HIST");
        h_m_pi0_coin->Draw("HIST SAME");
        if (h_coin_bgsub) { h_coin_bgsub->SetLineColor(kGreen+2); h_coin_bgsub->SetLineWidth(2); h_coin_bgsub->Draw("HIST SAME"); }
        if (h_final)      { h_final->SetLineColor(kMagenta+1); h_final->SetLineWidth(2); h_final->Draw("HIST SAME"); }

        TLegend *l2 = new TLegend(0.45,0.60,0.78,0.88); l2->SetBorderSize(0); l2->SetFillColor(0);
        l2->SetHeader("Analysis flow"); l2->SetTextSize(0.030);
        l2->AddEntry(h_mpi0_all,   "1) All #pi^{0} candidates",                     "l");
        l2->AddEntry(h_m_pi0_coin, "2) Within coincidence window (t1 & t2)",      "l");
        if (h_coin_bgsub) l2->AddEntry(h_coin_bgsub, "3) After timing (accidental) subtraction", "l");
        if (h_final)      l2->AddEntry(h_final,      "4) Final after combinatorial (fit) subtraction", "l");
        l2->Draw();
        c_pi0->SaveAs(Form("%s/pi0_mgg_inside_outside_coin_compare_run%d.png", outPlotDir.Data(), run));

        // 3) mm vs mgg
        TCanvas *c_mmiss_mgg = new TCanvas(name("c_mmiss_mgg"), "Mmiss vs Mgg", 900,700);
        c_mmiss_mgg->SetName(name("c_mmiss_mgg"));
        h_mmiss_vs_mgg->Draw("COLZ");
        c_mmiss_mgg->SaveAs(Form("%s/mmiss_vs_mgg_run%d.png", outPlotDir.Data(), run));

        // 4) cluster diagnostics
        TCanvas *c_cluster = new TCanvas(name("c_cluster"), "Cluster diagnostics", 1000,800);
        c_cluster->SetName(name("c_cluster"));
        c_cluster->Divide(2,2);
        c_cluster->cd(1); h_clustE->Draw();
        c_cluster->cd(2); h_clustT->Draw();
        c_cluster->cd(3); h_clustXY->Draw("COLZ");
        c_cluster->cd(4); h_clustE_sum->Draw();
        c_cluster->SaveAs(Form("%s/cluster_E_T_run%d.png", outPlotDir.Data(), run));

        // -------------------------
        // Compute final mmiss proton summary (choose best populated mmiss histogram)
        // -------------------------
        double mmiss_p_mean = 0.0, mmiss_p_sigma = 0.0;
        if (h_mmiss_2->GetEntries() > 50) {
        mmiss_p_mean = h_mmiss_2->GetMean();
        mmiss_p_sigma = h_mmiss_2->GetRMS();
        } else if (h_mmiss_3->GetEntries() > 50) {
        mmiss_p_mean = h_mmiss_3->GetMean();
        mmiss_p_sigma = h_mmiss_3->GetRMS();
        } else if (h_mmiss_4->GetEntries() > 50) {
        mmiss_p_mean = h_mmiss_4->GetMean();
        mmiss_p_sigma = h_mmiss_4->GetRMS();
        } else {
        if (h_mmiss_vs_mgg->GetEntries() > 0) {
            TH1D *hproj = h_mmiss_vs_mgg->ProjectionY(Form("hproj_mmiss_tmp_run%d", run));
            mmiss_p_mean = hproj->GetMean();
            mmiss_p_sigma = hproj->GetRMS();
            safe_delete(hproj);
        } else {
            mmiss_p_mean = mmiss_p_sigma = 0.0;
        }
        }

        // -------------------------
        // Capture HMS tracking efficiency data for integration into main CSV
        // -------------------------
        double hms_eff = 0.0, hms_eff_err = 0.0;
        double s1x_peak = 0.0, s1x_err = 0.0;
        double s1y_peak = 0.0, s1y_err = 0.0;
        double s2x_peak = 0.0, s2x_err = 0.0;
        double s2y_peak = 0.0, s2y_err = 0.0;
        
        try {
            auto last_entry = nps::hms_track_eff::HMSEffAggregator::Instance().GetLastRunEntry();
            if (last_entry.run == run) {
                // Data matches current run, use it
                hms_eff = last_entry.eff;
                hms_eff_err = last_entry.eff_err;
                s1x_peak = last_entry.s1x;
                s1x_err = last_entry.s1x_err;
                s1y_peak = last_entry.s1y;
                s1y_err = last_entry.s1y_err;
                s2x_peak = last_entry.s2x;
                s2x_err = last_entry.s2x_err;
                s2y_peak = last_entry.s2y;
                s2y_err = last_entry.s2y_err;
            }
        } catch (...) {
            // Aggregator access failed, use zeros (will be visible as missing data)
        }

        // -------------------------
        // CRITICAL: Write CSV data to per-run temp file BEFORE summary print
        // Bash will reconstruct global CSV from successful per-run files
        // -------------------------
        {
        auto write_csv_row = [&](std::ostream& out) {
            out << std::fixed;
            out << run << ","
                << std::setprecision(6) << accumulated_charge_mC << ","
                << std::setprecision(6) << run_current_mean << ","
                << std::setprecision(6) << cpu_lt << ","
                << std::setprecision(6) << beam_time << ","
                << std::setprecision(0) << n_total << ","
                << std::setprecision(0) << n_pass_hms << ","
                << std::setprecision(0) << n_selected_for_analysis << ","
                << std::setprecision(6) << bg.n_coin_raw << ","
                << std::setprecision(6) << bg.n_accidentals << ","
                << std::setprecision(6) << res.chi2_ndf << ","
                << std::setprecision(3) << res.mu_MeV << ","
                << std::setprecision(3) << res.sigma_MeV << ","
                << std::setprecision(3) << res.signal_counts << ","
                << std::setprecision(6) << mmiss_p_mean << ","
                << std::setprecision(6) << mmiss_p_sigma << ","
                << std::setprecision(6) << hms_eff << ","
                << std::setprecision(6) << hms_eff_err << ","
                << std::setprecision(6) << s1x_peak << ","
                << std::setprecision(6) << s1x_err << ","
                << std::setprecision(6) << s1y_peak << ","
                << std::setprecision(6) << s1y_err << ","
                << std::setprecision(6) << s2x_peak << ","
                << std::setprecision(6) << s2x_err << ","
                << std::setprecision(6) << s2y_peak << ","
                << std::setprecision(6) << s2y_err << ","
                << run_status << "\n";
        };

        TString per_run_csv = Form("/tmp/nps_csv_run_%d_%d.txt", run, (int)time(nullptr));
        ofstream fg(per_run_csv.Data(), ios::out);
        if (fg.is_open()) {
            write_csv_row(fg);
            fg.flush();
            fg.close();
            cout << "[CSV_WRITTEN] " << per_run_csv << "\n";
        } else {
            logmsg(WARN, Form("Could not write per-run CSV %s", per_run_csv.Data()));
        }

        ofstream fglobal(global_csv.Data(), ios::app);
        if (fglobal.is_open()) {
            write_csv_row(fglobal);
            fglobal.flush();
            fglobal.close();
        } else {
            logmsg(WARN, Form("Could not append global CSV %s (errno=%d, %s)",
                              global_csv.Data(), errno, std::strerror(errno)));
        }
        }

        // -------------------------
        // Print summary to console
        // -------------------------
        cout << "===== Run " << run << " summary =====\n";
        cout << " Total entries: " << n_total << "\n";
        cout << " Pass HMS: " << n_pass_hms << "\n";
        cout << " Pass HMS + NPS selection: " << n_selected_for_analysis << "\n";
        cout << " Coin raw (timing plane): " << bg.n_coin_raw << "\n";
        cout << " Estimated accidentals (time method): " << bg.n_accidentals << " +- " << bg.n_accidentals_err << "\n";
        cout << " Comb. BG fit χ2/ndf: " << std::fixed << std::setprecision(3) << res.chi2_ndf << "\n";
        cout << " Pi0 fit μ (MeV): " << res.mu_MeV << "   σ (MeV): " << res.sigma_MeV << "\n";
        cout << " Pi0 signal counts (final): " << res.signal_counts << "\n";
        cout << " Proton missing mass mean (GeV): " << mmiss_p_mean << "  sigma: " << mmiss_p_sigma << "\n";
        cout << "=====================================\n";

        // -------------------------
        // Save TXT summary (per-run)
        // -------------------------
        {
        TString txtout = Form("%s/summary_run%d.txt", outPlotDir.Data(), run);
        ofstream ftxt(txtout.Data());
        if (ftxt.is_open()) {
            ftxt << "Run " << run << " summary\n";
            ftxt << "Total entries: " << n_total << "\n";
            ftxt << "HMS passed: " << n_pass_hms << "\n";
            ftxt << "Selected for pi0 analysis: " << n_selected_for_analysis << "\n";
            ftxt << Form("Coin raw (timing plane): %.1f\n", bg.n_coin_raw);
            ftxt << Form("Estimated accidentals (time method): %.3f +- %.3f\n", bg.n_accidentals, bg.n_accidentals_err);
            ftxt << Form("Comb. BG fit chi2/ndf: %.3f\n", res.chi2_ndf);
            ftxt << Form("Pi0 fit mean: %.3f MeV\n", res.mu_MeV);
            ftxt << Form("Pi0 fit sigma: %.3f MeV\n", res.sigma_MeV);
            ftxt << Form("Pi0 signal (final): %.1f\n", res.signal_counts);
            ftxt << Form("Proton missing mass mean: %.6f GeV  sigma: %.6f GeV\n", mmiss_p_mean, mmiss_p_sigma);
            ftxt.close();
        } else {
            logmsg(WARN, Form("Could not open TXT summary %s for writing", txtout.Data()));
        }
        }

        // -------------------------
        // Clean up per-run objects (delete created objects to avoid memory leaks)
        // -------------------------
        
        // Batch cleanup: new weighted/exclusive histograms (1D)
        {
        TH1* hist_weighted[] = {
            h_pi0_weight, h_mmiss_all_weighted,
            h_Q2_weighted, h_W_weighted, h_t_weighted, h_tmin_weighted,
            h_pt_weighted, h_theta_weighted, h_phi_weighted, h_s_weighted,
            h_xB_weighted, h_z_weighted,
            h_Q2_excl_weighted, h_W_excl_weighted, h_t_excl_weighted, h_tmin_excl_weighted,
            h_pt_excl_weighted, h_theta_excl_weighted, h_phi_excl_weighted, h_s_excl_weighted,
            h_xB_excl_weighted, h_z_excl_weighted
        };
        for (int i = 0; i < 22; ++i) if (hist_weighted[i]) safe_delete(hist_weighted[i]);
        }

        // Batch cleanup: detector heatmap 2D histograms
        {
        TH1* hist_heatmap[] = {
            h_xB_xy, h_nu_xy, h_mx_xy, h_phi_xy, h_t_xy
        };
        for (int i = 0; i < 5; ++i) if (hist_heatmap[i]) safe_delete(hist_heatmap[i]);
        }

        // Batch cleanup: correlation 2D histograms
        {
        TH1* hist_corr[] = {
            h_Q2_vs_W, h_Q2_vs_t, h_Q2_vs_xB, h_W_vs_t, h_W_vs_xB,
            h_t_vs_xB, h_t_vs_phi, h_xB_vs_nu, h_Q2_vs_nu, h_W_vs_nu
        };
        for (int i = 0; i < 10; ++i) if (hist_corr[i]) safe_delete(hist_corr[i]);
        }

        // Batch cleanup: cluster histograms (mixed TH1D and TH2D)
        {
        TH1* hist_cluster[] = {
            h_nclusters, h_clustE, h_clustT, h_clustE_vs_T, h_clustXY, h_clustE_sum,
            h_opening_angle, h_photon_Eratio
        };
        for (int i = 0; i < 8; ++i) if (hist_cluster[i]) safe_delete(hist_cluster[i]);
        }

        // Batch cleanup: pi0 mass histograms
        {
        TH1* hist_pi0[] = {h_mpi0_all, h_mpi0_2, h_mpi0_3, h_mpi0_4};
        for (int i = 0; i < 4; ++i) if (hist_pi0[i]) safe_delete(hist_pi0[i]);
        }

        // Batch cleanup: missing mass histograms
        {
        TH1* hist_mmiss[] = {
            h_mmiss_2, h_mmiss_3, h_mmiss_4, h_mmiss_dvcs, h_mmiss_all, h_mmiss_all_corr
        };
        for (int i = 0; i < 6; ++i) if (hist_mmiss[i]) safe_delete(hist_mmiss[i]);
        }

        // Batch cleanup: timing & coincidence histograms (mixed TH1D and TH2D)
        {
        TH1* hist_timing[] = {h_t1_t2, h_t1_proj, h_t2_proj, h_m_pi0_coin, h_m_pi0_acc};
        for (int i = 0; i < 5; ++i) if (hist_timing[i]) safe_delete(hist_timing[i]);
        }

        // Clean up conditionally-created backgrounds
        if (h_coin_bgsub) safe_delete(h_coin_bgsub);
        if (h_final) safe_delete(h_final);

        // Batch cleanup: background & correlation histograms (mixed TH1D and TH2D)
        {
        TH1* hist_bg[] = {
            h_mgg_full1, h_mgg_full2, h_mmiss_vs_mgg, h_mmiss_vs_mgg_corr,
            h_mmiss_vs_t1, h_mmiss_vs_t2
        };
        for (int i = 0; i < 6; ++i) if (hist_bg[i]) safe_delete(hist_bg[i]);
        }

        // Batch cleanup: windowed mass histograms (diagonal/horizontal/vertical)
        for (auto *h: h_mgg_diag) if (h) safe_delete(h);
        for (auto *h: h_mgg_hor) if (h) safe_delete(h);
        for (auto *h: h_mgg_ver) if (h) safe_delete(h);

        // Batch cleanup: physics histograms
        {
        TH1* hist_physics[] = {
            h_t, h_tmin, h_pt, h_Q2, h_W, h_s, h_xB, h_z, h_theta, h_phi
        };
        for (int i = 0; i < 10; ++i) if (hist_physics[i]) safe_delete(hist_physics[i]);
        }

        // Clean up canvases (including new overlay canvases)
        {
        TCanvas* canvases[] = {
            c_t12, c_pi0, c_mmiss_mgg, c_cluster, c_mmiss_overlay, c_mmiss_fit,
            c_Q2_overlay, c_W_overlay, c_t_overlay, c_tmin_overlay, c_pt_overlay,
            c_theta_overlay, c_phi_overlay, c_s_overlay, c_xB_overlay, c_z_overlay
        };
        for (int i = 0; i < 16; ++i) if (canvases[i]) safe_delete(canvases[i]);
        }

        // Clean up fit functions
        if (fit1) safe_delete(fit1);
        if (fit2) safe_delete(fit2);
        if (fit_combined) safe_delete(fit_combined);

        // Clean up annotations & boxes
        safe_delete(l2);
        safe_delete(b_coin);
        safe_delete(b_full1);
        safe_delete(b_full2);
        for (auto *b: diag_boxes) if (b) safe_delete(b);
        for (auto *b: hor_boxes) if (b) safe_delete(b);
        for (auto *b: ver_boxes) if (b) safe_delete(b);
        safe_delete(txt);
        safe_delete(txt_fit);
        safe_delete(leg_mmiss);
        safe_delete(leg_mmiss_redraw);
        safe_delete(leg_fit);
        safe_delete(leg_Q2);
        safe_delete(leg_W);
        safe_delete(leg_t);
        safe_delete(leg_tmin);
        safe_delete(leg_pt);
        safe_delete(leg_s);
        safe_delete(leg_xB);
        safe_delete(leg_z);
        safe_delete(leg_theta);
        safe_delete(leg_phi);

        // Chain is automatically cleaned up by unique_ptr destructor
        sw_run.Stop();
        logmsg(INFO, Form("Run %d finished. Runtime: %f s (real)", run, sw_run.RealTime()));
        } catch (const std::exception& e) {
            logmsg(ERROR, Form("Exception during run %d processing: %s", run, e.what()));
            run_status = "ERROR";
            continue;  // unique_ptr cleans up automatically
        } catch (...) {
            logmsg(ERROR, Form("Unknown exception during run %d processing", run));
            run_status = "ERROR";
            continue;  // unique_ptr cleans up automatically
        }
    } // end runs

    // HMS track efficiency global summary (safe-guarded)
    try {
        nps::hms_track_eff::HMSEffAggregator::Instance().WriteGlobalSummary(outPlotDir_in);
    } catch (...) {
        logmsg(WARN, "HMSEffAggregator::WriteGlobalSummary threw an exception or failed - continuing.");
    }

    sw_total.Stop();
    logmsg(INFO, Form("ALL RUNS finished. Total runtime: %f s (real)", sw_total.RealTime()));
    
    logmsg(INFO, Form("Wrote global summary CSV to %s", global_csv.Data()));
    
    // Exit cleanly to avoid ROOT cleanup segfaults
    gApplication->Terminate(0);
}
