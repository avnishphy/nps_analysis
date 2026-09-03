// ============================================================================
// File: nps_analysis_main.C
// Purpose: Unified NPS π0 analysis pipeline for HCANA and waveform ROOT inputs
// Author: Avnish Singh (physics), ChatGPT (refactoring 2026)
//
// Design:
//   - RAII wrappers for all ROOT objects
//   - Per-run memory isolation (histograms deleted after each run)
//   - Exception-safe event processing with catch-and-continue
//   - Seamless integration with header file utilities
//   - Physics calculations preserved exactly (no algorithm changes)
//
// Build: root -l -b src/analysis/nps_analysis_main.C
// ============================================================================
#include "utils.C"
#include <iomanip>
#include <TF1.h>
#include <TFitResult.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLine.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TH2D.h>
#include "nps_helper.h"
#include "nps_time_bg.h"
#include "nps_comb_bg.h"
#include "nps_mmiss_cor.h"
#include "acceptance_cuts.h"
#include "nps_2d_mass_cut.h"
#include "nps_dead_block.h"

// Suppress empty-body warning in physics_var.h
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wempty-body"
#include "nps_physics_var.h"
#pragma GCC diagnostic pop

// ROOT headers
#include <TFile.h>
#include <TApplication.h>
#include <TChain.h>
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
#include <TObjArray.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TMath.h>
#include <TF1.h>
#include <TPaveText.h>

// C++ standard library
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <unordered_set>
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
#include <cctype>
#include <cstdlib>
#include <unistd.h>

using namespace std;

// ============================================================================
// Configuration Constants
// ============================================================================
constexpr int MAX_CLUS = 20;
constexpr int MAX_CLUS_INPUT = 1080;  // NPS block count; ROOT I/O buffer only
constexpr double EBEAM_DEFAULT = 10.538;
constexpr double HMS_MOM_OFFSET_SCALE = 1.0;
constexpr Long64_t MIN_PRINT_EVERY = 1000;
constexpr const char* CANONICAL_OUTPUT_BASE = "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/output";
constexpr const char* MAIN_CONFIG_CSV = "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/config/nps_dvcs_all_kins_main.csv";
constexpr const char* CHARGE_FALLBACK_DB = "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/config/DataBase_production_runs_newBCMOffset.txt";
constexpr const char* ACCEPTANCE_CUTS_CONFIG_DEFAULT = "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/config/acceptance_cuts.conf";
constexpr const char* DEAD_BLOCK_CONFIG_DEFAULT = "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/config/dead_block_per_runs.csv";
constexpr const char* EFFICIENCY_STUFF_DIR = "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/output/efficiency_stuff";
constexpr const char* DEFAULT_WAVEFORM_INPUT_DIR = "/lustre24/expphy/volatile/hallc/nps/hhuang/farmFile/Production/DVCS";
constexpr const char* DEFAULT_HCANA_INPUT_DIR = "/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/updated";

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
    f << "run,accumulated_charge(mC),hel_pos_charge(mC),hel_neg_charge(mC),current_mean_uA,Beam_Time(s),"
      << "total_entries,pass_hms,pass_hms_nps,total_coin_entries,"
      << "estimated_time_accidentals,chi2_ndf_comb_bg,pi0_mu_MeV,pi0_sigma_MeV,"
      << "pi0_signal_counts,mmiss_p_mean_GeV,mmiss_p_sigma_GeV,"
      << "s1x_peak,s1x_err,s1y_peak,s1y_err,s2x_peak,s2x_err,s2y_peak,s2y_err,run_status\n";
    f.close();
    return true;
}

struct RunSummaryRow {
    int run = -1;
    double accumulated_charge_mC = 0.0;
    double hel_pos_charge_mC = 0.0;
    double hel_neg_charge_mC = 0.0;
    double current_mean_uA = 0.0;
    double beam_time_s = 0.0;
    Long64_t total_entries = 0;
    Long64_t pass_hms = 0;
    Long64_t pass_hms_nps = 0;
    double total_coin_entries = 0.0;
    double estimated_time_accidentals = 0.0;
    double chi2_ndf_comb_bg = 0.0;
    double pi0_mu_MeV = 0.0;
    double pi0_sigma_MeV = 0.0;
    double pi0_signal_counts = 0.0;
    double mmiss_p_mean_GeV = 0.0;
    double mmiss_p_sigma_GeV = 0.0;
    double s1x_peak = 0.0;
    double s1x_err = 0.0;
    double s1y_peak = 0.0;
    double s1y_err = 0.0;
    double s2x_peak = 0.0;
    double s2x_err = 0.0;
    double s2y_peak = 0.0;
    double s2y_err = 0.0;
    std::string run_status;
};

struct TreeEntry {
    Double_t mpi0_all = 0.0;
    Double_t mmiss_all = 0.0;
    Double_t mmiss_all_corr = 0.0;
    Double_t mmiss_all_no_mom_offset = 0.0;
    Double_t mmiss_all_corr_no_mom_offset = 0.0;
    Double_t pi0_weight = 0.0;
    Int_t is_exclusive = 0;
    Int_t is_exclusive_ellipse = 0;
    Int_t is_exclusive_mcd = 0;
    Int_t is_weighted = 0;
    Int_t helicity = 0;
    Double_t Q2 = 0.0;
    Double_t W = 0.0;
    Double_t t = 0.0;
    Double_t tmin = 0.0;
    Double_t pt = 0.0;
    Double_t theta = 0.0;
    Double_t phi = 0.0;
    Double_t s = 0.0;
    Double_t xB = 0.0;
    Double_t z = 0.0;
    Int_t nclust_selected = 0;
    Int_t event_id = 0;
    Double_t cluster_x_1 = 0.0;
    Double_t cluster_y_1 = 0.0;
    Double_t cluster_e_1 = 0.0;
    Double_t cluster_x_2 = 0.0;
    Double_t cluster_y_2 = 0.0;
    Double_t cluster_e_2 = 0.0;
    Double_t delta = 0.0;
    Double_t xptar = 0.0;
    Double_t yptar = 0.0;
    Double_t xtar = 0.0;
    Double_t ytar = 0.0;
    Double_t xfp = 0.0;
    Double_t yfp = 0.0;
    Double_t xpfp = 0.0;
    Double_t ypfp = 0.0;
};

inline void write_event_physics_tree(TFile* fout,
                                     const std::map<Long64_t, TreeEntry>& tree_data)
{
    if (!fout) return;
    fout->cd();

    TTree *treeOut = new TTree("physics", "Event-level physics data with weights and exclusivity flags");

    Int_t event_id = 0;
    Double_t mpi0_all = 0, mmiss_all = 0, mmiss_all_corr = 0, mmiss_all_no_mom_offset = 0, mmiss_all_corr_no_mom_offset = 0, pi0_weight = 0;
    Int_t is_exclusive = 0, is_exclusive_ellipse = 0, is_exclusive_mcd = 0, is_weighted = 0, helicity = 0;
    Double_t Q2 = 0, W = 0, t = 0, tmin = 0, pt = 0;
    Double_t theta = 0, phi = 0, s = 0, xB = 0, z = 0;
    Int_t nclust_selected = 0;
    Double_t cluster_x_1 = 0, cluster_y_1 = 0, cluster_e_1 = 0;
    Double_t cluster_x_2 = 0, cluster_y_2 = 0, cluster_e_2 = 0;
    Double_t delta = 0, xptar = 0, yptar = 0, xtar = 0, ytar = 0;
    Double_t xfp = 0, yfp = 0, xpfp = 0, ypfp = 0;
    Double_t exclusive_sep = 0;

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
    treeOut->Branch("helicity", &helicity, "helicity/I");

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

    for (const auto& pair : tree_data) {
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
        helicity = entry.helicity;
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
        delta = entry.delta;
        xptar = entry.xptar;
        yptar = entry.yptar;
        xtar = entry.xtar;
        ytar = entry.ytar;
        xfp = entry.xfp;
        yfp = entry.yfp;
        xpfp = entry.xpfp;
        ypfp = entry.ypfp;
        exclusive_sep = 1.0 / (std::pow(std::fabs(entry.mmiss_all_corr - 0.938), 0.5) + 0.01);
        treeOut->Fill();
    }

    treeOut->Write();
    safe_delete(treeOut);
}

inline void write_summary_csv_row(std::ostream& out, const RunSummaryRow& row)
{
    out << std::fixed;
    out << row.run << ","
        << std::setprecision(6) << row.accumulated_charge_mC << ","
        << std::setprecision(6) << row.hel_pos_charge_mC << ","
        << std::setprecision(6) << row.hel_neg_charge_mC << ","
        << std::setprecision(6) << row.current_mean_uA << ","
        << std::setprecision(6) << row.beam_time_s << ","
        << std::setprecision(0) << row.total_entries << ","
        << std::setprecision(0) << row.pass_hms << ","
        << std::setprecision(0) << row.pass_hms_nps << ","
        << std::setprecision(6) << row.total_coin_entries << ","
        << std::setprecision(6) << row.estimated_time_accidentals << ","
        << std::setprecision(6) << row.chi2_ndf_comb_bg << ","
        << std::setprecision(3) << row.pi0_mu_MeV << ","
        << std::setprecision(3) << row.pi0_sigma_MeV << ","
        << std::setprecision(3) << row.pi0_signal_counts << ","
        << std::setprecision(6) << row.mmiss_p_mean_GeV << ","
        << std::setprecision(6) << row.mmiss_p_sigma_GeV << ","
        << std::setprecision(6) << row.s1x_peak << ","
        << std::setprecision(6) << row.s1x_err << ","
        << std::setprecision(6) << row.s1y_peak << ","
        << std::setprecision(6) << row.s1y_err << ","
        << std::setprecision(6) << row.s2x_peak << ","
        << std::setprecision(6) << row.s2x_err << ","
        << std::setprecision(6) << row.s2y_peak << ","
        << std::setprecision(6) << row.s2y_err << ","
        << row.run_status << "\n";
}

inline void write_summary_txt_file(const TString& out_path,
                                   const RunSummaryRow& row,
                                   double accidental_err)
{
    std::ofstream ftxt(out_path.Data());
    if (!ftxt.is_open()) {
        logmsg(WARN, Form("Could not open TXT summary %s for writing", out_path.Data()));
        return;
    }

    ftxt << "Run " << row.run << " summary\n";
    ftxt << "Total entries: " << row.total_entries << "\n";
    ftxt << Form("Accumulated good-event charge: %.6f mC\n", row.accumulated_charge_mC);
    ftxt << Form("HEL+ good-event charge: %.6f mC\n", row.hel_pos_charge_mC);
    ftxt << Form("HEL- good-event charge: %.6f mC\n", row.hel_neg_charge_mC);
    ftxt << "HMS passed: " << row.pass_hms << "\n";
    ftxt << "Selected for pi0 analysis: " << row.pass_hms_nps << "\n";
    ftxt << Form("Coin raw (timing plane): %.1f\n", row.total_coin_entries);
    ftxt << Form("Estimated accidentals (time method): %.3f +- %.3f\n",
                 row.estimated_time_accidentals,
                 accidental_err);
    ftxt << Form("Comb. BG fit chi2/ndf: %.3f\n", row.chi2_ndf_comb_bg);
    ftxt << Form("Pi0 fit mean: %.3f MeV\n", row.pi0_mu_MeV);
    ftxt << Form("Pi0 fit sigma: %.3f MeV\n", row.pi0_sigma_MeV);
    ftxt << Form("Pi0 signal (final): %.1f\n", row.pi0_signal_counts);
    ftxt << Form("Proton missing mass mean: %.6f GeV  sigma: %.6f GeV\n",
                 row.mmiss_p_mean_GeV,
                 row.mmiss_p_sigma_GeV);
}

enum class AnalysisMode {
    kAuto,
    kHcana,
    kWaveform
};

struct BranchMapping {
    AnalysisMode mode = AnalysisMode::kAuto;
    TString mode_label;
    TString nclust_branch;
    TString clusE_branch;
    TString clusX_branch;
    TString clusY_branch;
    TString clusT_branch;
    double cluster_time_offset_ns = 0.0;
    double time_window_wrt_150_ns = 10.0;
    double time_threshold_ns = 10.0;
    bool use_shifted_timing_windows = false;
};

struct RunMeta {
    int run = -1;
    TString kin;
    TString type;
    double beam_energy = EBEAM_DEFAULT;
    double nps_theta_deg = std::numeric_limits<double>::quiet_NaN();
    double nps_target_distance_cm = nps::kDefaultZ_NPS_cm;
    double charge_mC = std::numeric_limits<double>::quiet_NaN();
    double current_uA = std::numeric_limits<double>::quiet_NaN();
    double beam_time_s = std::numeric_limits<double>::quiet_NaN();
};

struct OutputFolders {
    TString base;
    TString logs;
    TString root;
    TString plots;
    TString tables;
    TString summary;
};

inline std::string nps_trim_copy(const std::string& s)
{
    size_t first = 0;
    while (first < s.size() && std::isspace(static_cast<unsigned char>(s[first]))) {
        ++first;
    }
    size_t last = s.size();
    while (last > first && std::isspace(static_cast<unsigned char>(s[last - 1]))) {
        --last;
    }
    return s.substr(first, last - first);
}

inline std::string nps_to_lower_copy(const std::string& s)
{
    std::string out = s;
    std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c){
        return static_cast<char>(std::tolower(c));
    });
    return out;
}

inline std::vector<std::string> split_csv_quoted(const std::string& line)
{
    std::vector<std::string> fields;
    std::string cur;
    bool in_quotes = false;

    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (in_quotes && (i + 1 < line.size()) && line[i + 1] == '"') {
                cur.push_back('"');
                ++i;
            } else {
                in_quotes = !in_quotes;
            }
        } else if (c == ',' && !in_quotes) {
            fields.push_back(nps_trim_copy(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    fields.push_back(nps_trim_copy(cur));
    return fields;
}

inline bool parse_int_field(const std::vector<std::string>& row, int idx, int& out)
{
    if (idx < 0 || idx >= static_cast<int>(row.size())) return false;
    try {
        out = std::stoi(row[idx]);
        return true;
    } catch (...) {
        return false;
    }
}

inline bool parse_double_field(const std::vector<std::string>& row, int idx, double& out)
{
    if (idx < 0 || idx >= static_cast<int>(row.size())) return false;
    try {
        out = std::stod(row[idx]);
        return true;
    } catch (...) {
        return false;
    }
}

inline std::set<std::string> parse_type_filter(const std::string& types_csv)
{
    std::set<std::string> out;
    std::stringstream ss(types_csv);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        tok = nps_trim_copy(tok);
        if (!tok.empty()) {
            out.insert(nps_to_lower_copy(tok));
        }
    }
    return out;
}

inline AnalysisMode parse_mode_string(const std::string& mode)
{
    const std::string m = nps_to_lower_copy(nps_trim_copy(mode));
    if (m == "hcana" || m == "standard") return AnalysisMode::kHcana;
    if (m == "waveform" || m == "wf" || m == "wfpi0") return AnalysisMode::kWaveform;
    return AnalysisMode::kAuto;
}

inline const char* mode_to_cstr(AnalysisMode mode)
{
    if (mode == AnalysisMode::kHcana) return "hcana";
    if (mode == AnalysisMode::kWaveform) return "waveform";
    return "auto";
}

inline BranchMapping make_branch_map(AnalysisMode mode)
{
    BranchMapping m;
    m.mode = mode;
    if (mode == AnalysisMode::kWaveform) {
        m.mode_label = "waveform";
        m.nclust_branch = "NPS.prod.nclust";
        m.clusE_branch = "NPS.prod.clusE";
        m.clusX_branch = "NPS.prod.clusX";
        m.clusY_branch = "NPS.prod.clusY";
        m.clusT_branch = "NPS.prod.clusT";
        m.cluster_time_offset_ns = 150.0;
        m.time_window_wrt_150_ns = 13.0;
        m.time_threshold_ns = 13.0;
        m.use_shifted_timing_windows = true;
    } else {
        m.mode_label = "hcana";
        m.nclust_branch = "NPS.cal.nclust";
        m.clusE_branch = "NPS.cal.clusE";
        m.clusX_branch = "NPS.cal.clusX";
        m.clusY_branch = "NPS.cal.clusY";
        m.clusT_branch = "NPS.cal.clusT";
        m.cluster_time_offset_ns = 0.0;
        m.time_window_wrt_150_ns = 10.0;
        m.time_threshold_ns = 10.0;
        m.use_shifted_timing_windows = false;
    }
    return m;
}

inline bool detect_branch_mapping(TTree* tree, AnalysisMode requested_mode, BranchMapping& out_map, std::string& error)
{
    if (!tree) {
        error = "null tree pointer";
        return false;
    }

    const bool has_waveform = (tree->GetBranch("NPS.prod.nclust") && tree->GetBranch("NPS.prod.clusE"));
    const bool has_hcana = (tree->GetBranch("NPS.cal.nclust") && tree->GetBranch("NPS.cal.clusE"));

    if (requested_mode == AnalysisMode::kWaveform) {
        if (!has_waveform) {
            error = "requested waveform mode, but waveform branches are not present";
            return false;
        }
        out_map = make_branch_map(AnalysisMode::kWaveform);
        return true;
    }

    if (requested_mode == AnalysisMode::kHcana) {
        if (!has_hcana) {
            error = "requested hcana mode, but HCANA branches are not present";
            return false;
        }
        out_map = make_branch_map(AnalysisMode::kHcana);
        return true;
    }

    if (has_waveform) {
        out_map = make_branch_map(AnalysisMode::kWaveform);
        return true;
    }
    if (has_hcana) {
        out_map = make_branch_map(AnalysisMode::kHcana);
        return true;
    }

    error = "cannot detect input mode from branches (expected NPS.prod.* or NPS.cal.* branches)";
    return false;
}

inline std::string sanitize_token(const std::string& s)
{
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-') {
            out.push_back(c);
        } else {
            out.push_back('_');
        }
    }
    return out;
}

inline bool ensure_output_folders(const TString& output_base, const TString& kin, OutputFolders& out, std::string& error)
{
    if (output_base.IsNull() || kin.IsNull()) {
        error = "empty output base or kinematic setting";
        return false;
    }

    const std::string kin_safe = sanitize_token(kin.Data());
    out.base = TString::Format("%s/%s", output_base.Data(), kin_safe.c_str());
    out.logs = out.base + "/logs";
    out.root = out.base + "/root";
    out.plots = out.base + "/plots";
    out.tables = out.base + "/tables";
    out.summary = out.base + "/summary";

    const TString dirs[] = {out.base, out.logs, out.root, out.plots, out.tables, out.summary};
    for (const auto& d : dirs) {
        if (gSystem->mkdir(d, true) != 0 && gSystem->AccessPathName(d) != 0) {
            error = std::string("failed to create output directory: ") + d.Data();
            return false;
        }
    }
    return true;
}

inline bool load_config_rows(const TString& config_csv,
                             const TString& kin,
                             const std::set<std::string>& allowed_types,
                             int run_filter,
                             std::vector<RunMeta>& rows,
                             std::string& error)
{
    std::ifstream in(config_csv.Data());
    if (!in.is_open()) {
        error = std::string("cannot open config CSV: ") + config_csv.Data();
        return false;
    }

    std::string line;
    if (!std::getline(in, line)) {
        error = "config CSV is empty";
        return false;
    }

    const std::vector<std::string> header = split_csv_quoted(line);
    std::unordered_map<std::string, int> col;
    for (size_t i = 0; i < header.size(); ++i) {
        col[header[i]] = static_cast<int>(i);
    }

    const char* required_cols[] = {
        "Kin_old", "run_number", "Type", "BeamEnergy", "NPS_Thet", "NPS_target_distance"
    };
    for (const char* req : required_cols) {
        if (col.find(req) == col.end()) {
            error = std::string("missing required config column: ") + req;
            return false;
        }
    }

    const int idx_kin = col["Kin_old"];
    const int idx_run = col["run_number"];
    const int idx_type = col["Type"];
    const int idx_beam = col["BeamEnergy"];
    const int idx_nps_theta = col["NPS_Thet"];
    const int idx_nps_dist = col["NPS_target_distance"];
    const int idx_charge = col.count("charge(mC)") ? col["charge(mC)"] : -1;
    const int idx_current = col.count("Ibeam(uA)") ? col["Ibeam(uA)"] : -1;
    const int idx_run_len = col.count("Run_Length(ks)") ? col["Run_Length(ks)"] : -1;

    std::map<int, RunMeta> dedup;
    int malformed = 0;
    while (std::getline(in, line)) {
        if (nps_trim_copy(line).empty()) continue;
        const std::vector<std::string> row = split_csv_quoted(line);
        if (static_cast<int>(row.size()) <= std::max(idx_run, idx_nps_dist)) {
            ++malformed;
            continue;
        }

        RunMeta m;
        m.kin = row[idx_kin].c_str();
        m.type = row[idx_type].c_str();

        if (kin.Length() > 0 && m.kin != kin) continue;
        if (!allowed_types.empty()) {
            const std::string type_l = nps_to_lower_copy(m.type.Data());
            if (allowed_types.find(type_l) == allowed_types.end()) continue;
        }
        if (!parse_int_field(row, idx_run, m.run)) {
            ++malformed;
            continue;
        }
        if (run_filter > 0 && m.run != run_filter) continue;
        if (!parse_double_field(row, idx_beam, m.beam_energy)) {
            ++malformed;
            continue;
        }
        if (!parse_double_field(row, idx_nps_theta, m.nps_theta_deg)) {
            ++malformed;
            continue;
        }
        if (!parse_double_field(row, idx_nps_dist, m.nps_target_distance_cm)) {
            ++malformed;
            continue;
        }
        parse_double_field(row, idx_charge, m.charge_mC);
        parse_double_field(row, idx_current, m.current_uA);
        double run_length_ks = std::numeric_limits<double>::quiet_NaN();
        if (parse_double_field(row, idx_run_len, run_length_ks)) {
            m.beam_time_s = run_length_ks * 1000.0;
        }

        auto it = dedup.find(m.run);
        if (it == dedup.end()) {
            dedup[m.run] = m;
        } else {
            const RunMeta& old = it->second;
            if (old.kin != m.kin || std::fabs(old.beam_energy - m.beam_energy) > 1e-6 ||
                std::fabs(old.nps_target_distance_cm - m.nps_target_distance_cm) > 1e-6) {
                error = Form("inconsistent metadata for run %d in config CSV", m.run);
                return false;
            }
        }
    }

    if (dedup.empty()) {
        error = Form("no runs selected from %s for kin=%s", config_csv.Data(), kin.Data());
        return false;
    }

    rows.clear();
    rows.reserve(dedup.size());
    for (const auto& kv : dedup) {
        rows.push_back(kv.second);
    }

    if (malformed > 0) {
        logmsg(WARN, Form("Skipped %d malformed config rows while parsing %s", malformed, config_csv.Data()));
    }
    return true;
}

inline bool load_charge_fallback_db(const TString& db_path,
                                    std::map<int, double>& charge_mC_by_run,
                                    std::string& error)
{
    std::ifstream in(db_path.Data());
    if (!in.is_open()) {
        error = std::string("cannot open charge fallback DB: ") + db_path.Data();
        return false;
    }

    std::string line;
    std::vector<std::string> header;
    while (std::getline(in, line)) {
        const std::string trimmed = nps_trim_copy(line);
        if (trimmed.empty()) continue;

        std::stringstream ss(trimmed);
        std::vector<std::string> fields;
        std::string tok;
        while (ss >> tok) {
            fields.push_back(tok);
        }
        if (!fields.empty() && fields.front() == "RunNo") {
            header = fields;
            break;
        }
    }

    if (header.empty()) {
        error = std::string("missing header row in charge fallback DB: ") + db_path.Data();
        return false;
    }

    int idx_run = -1;
    int idx_charge_tot = -1;
    for (size_t i = 0; i < header.size(); ++i) {
        if (header[i] == "RunNo") idx_run = static_cast<int>(i);
        if (header[i] == "Charge_tot") idx_charge_tot = static_cast<int>(i);
    }
    if (idx_run < 0 || idx_charge_tot < 0) {
        error = "charge fallback DB must contain RunNo and Charge_tot columns";
        return false;
    }

    int malformed = 0;
    charge_mC_by_run.clear();
    while (std::getline(in, line)) {
        const std::string trimmed = nps_trim_copy(line);
        if (trimmed.empty()) continue;

        std::stringstream ss(trimmed);
        std::vector<std::string> fields;
        std::string tok;
        while (ss >> tok) {
            fields.push_back(tok);
        }

        if (static_cast<int>(fields.size()) <= std::max(idx_run, idx_charge_tot)) {
            ++malformed;
            continue;
        }

        try {
            const int run = std::stoi(fields[idx_run]);
            const double charge_uC = std::stod(fields[idx_charge_tot]);
            if (std::isfinite(charge_uC) && charge_uC > 0.0) {
                charge_mC_by_run[run] = charge_uC / 1000.0;
            } else {
                ++malformed;
            }
        } catch (...) {
            ++malformed;
        }
    }

    if (charge_mC_by_run.empty()) {
        error = Form("no usable Charge_tot values found in %s", db_path.Data());
        return false;
    }

    if (malformed > 0) {
        logmsg(WARN, Form("Skipped %d malformed charge fallback rows while parsing %s", malformed, db_path.Data()));
    }
    return true;
}

inline void fill_missing_charge_from_fallback(std::vector<RunMeta>& rows,
                                              const TString& db_path)
{
    bool need_fallback = false;
    for (const auto& row : rows) {
        if (!std::isfinite(row.charge_mC) || row.charge_mC <= 0.0) {
            need_fallback = true;
            break;
        }
    }
    if (!need_fallback) return;

    std::map<int, double> charge_mC_by_run;
    std::string error;
    if (!load_charge_fallback_db(db_path, charge_mC_by_run, error)) {
        logmsg(WARN, Form("Charge fallback unavailable: %s", error.c_str()));
        return;
    }

    int filled = 0;
    int missing = 0;
    for (auto& row : rows) {
        if (std::isfinite(row.charge_mC) && row.charge_mC > 0.0) continue;

        const auto it = charge_mC_by_run.find(row.run);
        if (it == charge_mC_by_run.end()) {
            ++missing;
            continue;
        }

        row.charge_mC = it->second;
        ++filled;
    }

    if (filled > 0) {
        logmsg(INFO, Form("Filled charge metadata for %d run(s) from %s using Charge_tot/1000", filled, db_path.Data()));
    }
    if (missing > 0) {
        logmsg(WARN, Form("Charge metadata still missing for %d selected run(s); no Charge_tot fallback row found", missing));
    }
}

struct EventRange {
    Long64_t low = 0;
    Long64_t high = -1;
};

struct GoodEventCharge {
    double hel_pos_mC = 0.0;
    double hel_neg_mC = 0.0;

    double total_mC() const { return hel_pos_mC + hel_neg_mC; }
};

inline bool parse_int64_token(const std::string& token, Long64_t& out)
{
    const std::string s = nps_trim_copy(token);
    if (s.empty()) return false;
    try {
        size_t pos = 0;
        const long long val = std::stoll(s, &pos, 10);
        if (pos != s.size()) return false;
        out = static_cast<Long64_t>(val);
        return true;
    } catch (...) {
        return false;
    }
}

inline bool parse_event_range_list(const std::string& field,
                                   std::vector<EventRange>& out,
                                   int& malformed_tokens)
{
    std::stringstream ss(field);
    std::string token;
    bool has_any = false;

    while (std::getline(ss, token, ';')) {
        token = nps_trim_copy(token);
        if (token.empty()) continue;

        Long64_t low = 0;
        Long64_t high = 0;
        const size_t dash = token.find('-');

        if (dash == std::string::npos) {
            if (!parse_int64_token(token, low)) {
                ++malformed_tokens;
                continue;
            }
            high = low;
        } else {
            const std::string left = nps_trim_copy(token.substr(0, dash));
            const std::string right = nps_trim_copy(token.substr(dash + 1));
            if (!parse_int64_token(left, low) || !parse_int64_token(right, high)) {
                ++malformed_tokens;
                continue;
            }
            if (high < low) std::swap(low, high);
        }

        out.push_back(EventRange{low, high});
        has_any = true;
    }

    return has_any;
}

inline void normalize_event_ranges(std::vector<EventRange>& ranges)
{
    if (ranges.empty()) return;

    std::sort(ranges.begin(), ranges.end(), [](const EventRange& a, const EventRange& b) {
        if (a.low != b.low) return a.low < b.low;
        return a.high < b.high;
    });

    std::vector<EventRange> merged;
    merged.reserve(ranges.size());
    merged.push_back(ranges.front());

    for (size_t i = 1; i < ranges.size(); ++i) {
        EventRange& last = merged.back();
        const EventRange& cur = ranges[i];
        // Selection-report event ranges are half-open: [low, high).
        if (cur.low <= last.high) {
            if (cur.high > last.high) last.high = cur.high;
        } else {
            merged.push_back(cur);
        }
    }

    ranges.swap(merged);
}

inline bool load_good_event_ranges_by_run(const TString& report_csv,
                                          const TString& expected_kin,
                                          std::unordered_map<int, std::vector<EventRange>>& ranges_by_run,
                                          std::unordered_map<int, GoodEventCharge>& charge_by_run,
                                          std::string& error)
{
    ranges_by_run.clear();
    charge_by_run.clear();

    std::ifstream in(report_csv.Data());
    if (!in.is_open()) {
        error = std::string("cannot open selection report CSV: ") + report_csv.Data();
        return false;
    }

    std::string line;
    if (!std::getline(in, line)) {
        error = std::string("selection report CSV is empty: ") + report_csv.Data();
        return false;
    }

    const std::vector<std::string> header = split_csv_quoted(line);
    std::unordered_map<std::string, int> col_l;
    for (size_t i = 0; i < header.size(); ++i) {
        col_l[nps_to_lower_copy(nps_trim_copy(header[i]))] = static_cast<int>(i);
    }

    const int idx_run = col_l.count("run_number") ? col_l["run_number"] : -1;
    int idx_kin = -1;
    if (col_l.count("kinematic_setting")) idx_kin = col_l["kinematic_setting"];
    else if (col_l.count("kin_old")) idx_kin = col_l["kin_old"];

    int idx_ranges = -1;
    if (col_l.count("accepted_gevnum_ranges")) {
        idx_ranges = col_l["accepted_gevnum_ranges"];
    }
    const int idx_hel_pos_charge = col_l.count("hel_pos_charge_after_cut_uc") ? col_l["hel_pos_charge_after_cut_uc"] : -1;
    const int idx_hel_neg_charge = col_l.count("hel_neg_charge_after_cut_uc") ? col_l["hel_neg_charge_after_cut_uc"] : -1;

    if (idx_run < 0) {
        error = std::string("selection report is missing run_number column: ") + report_csv.Data();
        return false;
    }
    if (idx_ranges < 0) {
        error = std::string("selection report is missing required accepted_gevnum_ranges column; refusing to use accepted_evnumber_ranges as T->g.evnum ranges: ") + report_csv.Data();
        return false;
    }
    if (idx_hel_pos_charge < 0 || idx_hel_neg_charge < 0) {
        error = std::string("selection report is missing HEL_pos_charge_after_cut_uC or HEL_neg_charge_after_cut_uC column: ") + report_csv.Data();
        return false;
    }

    int malformed_tokens = 0;
    int malformed_charge_rows = 0;
    int rows_for_kin = 0;
    int rows_with_ranges = 0;

    while (std::getline(in, line)) {
        if (nps_trim_copy(line).empty()) continue;
        const std::vector<std::string> row = split_csv_quoted(line);

        if (idx_kin >= 0 && idx_kin < static_cast<int>(row.size()) && expected_kin.Length() > 0) {
            const std::string row_kin = nps_trim_copy(row[idx_kin]);
            if (row_kin != expected_kin.Data()) continue;
        }
        ++rows_for_kin;

        int run = -1;
        if (!parse_int_field(row, idx_run, run)) continue;
        if (idx_ranges < 0 || idx_ranges >= static_cast<int>(row.size())) continue;

        const std::string ranges_field = nps_trim_copy(row[idx_ranges]);
        if (ranges_field.empty()) continue;

        std::vector<EventRange> parsed;
        if (!parse_event_range_list(ranges_field, parsed, malformed_tokens)) continue;
        if (parsed.empty()) continue;

        std::vector<EventRange>& dst = ranges_by_run[run];
        dst.insert(dst.end(), parsed.begin(), parsed.end());

        double hel_pos_charge_uC = 0.0;
        double hel_neg_charge_uC = 0.0;
        if (parse_double_field(row, idx_hel_pos_charge, hel_pos_charge_uC) &&
            parse_double_field(row, idx_hel_neg_charge, hel_neg_charge_uC) &&
            std::isfinite(hel_pos_charge_uC) &&
            std::isfinite(hel_neg_charge_uC) &&
            hel_pos_charge_uC >= 0.0 &&
            hel_neg_charge_uC >= 0.0) {
            GoodEventCharge& charge = charge_by_run[run];
            charge.hel_pos_mC += hel_pos_charge_uC / 1000.0;
            charge.hel_neg_mC += hel_neg_charge_uC / 1000.0;
        } else {
            ++malformed_charge_rows;
        }

        ++rows_with_ranges;
    }

    if (rows_for_kin == 0) {
        error = Form("selection report %s has no rows for kinematic setting %s",
                     report_csv.Data(), expected_kin.Data());
        return false;
    }
    if (rows_with_ranges == 0 || ranges_by_run.empty()) {
        error = Form("selection report %s has no usable T->g.evnum ranges for kinematic setting %s",
                     report_csv.Data(), expected_kin.Data());
        return false;
    }

    for (auto& kv : ranges_by_run) {
        normalize_event_ranges(kv.second);
    }

    if (malformed_tokens > 0) {
        logmsg(WARN, Form("Selection report %s contained %d malformed range token(s) that were ignored",
                          report_csv.Data(), malformed_tokens));
    }
    if (malformed_charge_rows > 0) {
        logmsg(WARN, Form("Selection report %s contained %d row(s) with unusable helicity charge fields",
                          report_csv.Data(), malformed_charge_rows));
    }

    return true;
}

inline bool is_interactive_terminal()
{
    return (::isatty(STDIN_FILENO) != 0) && (::isatty(STDOUT_FILENO) != 0);
}

inline bool choose_use_good_event_cut(const TString& report_csv)
{
    const char* env_choice = gSystem->Getenv("NPS_USE_GEVNUM_CUT");
    if (env_choice && std::strlen(env_choice) > 0) {
        const std::string choice = nps_to_lower_copy(nps_trim_copy(env_choice));
        if (choice == "yes" || choice == "y" || choice == "1" || choice == "true" || choice == "use" || choice == "on") {
            logmsg(INFO, Form("Using optional T->g.evnum cut (NPS_USE_GEVNUM_CUT=%s)", env_choice));
            return true;
        }
        if (choice == "no" || choice == "n" || choice == "0" || choice == "false" || choice == "skip" || choice == "off") {
            logmsg(INFO, Form("Skipping optional T->g.evnum cut (NPS_USE_GEVNUM_CUT=%s)", env_choice));
            return false;
        }
        if (choice != "ask" && choice != "prompt") {
            logmsg(WARN, Form("Unrecognized NPS_USE_GEVNUM_CUT=%s. Falling back to prompt/default behavior.", env_choice));
        }
    }

    if (!is_interactive_terminal()) {
        logmsg(INFO, "Non-interactive mode: defaulting to use optional T->g.evnum cut. Set NPS_USE_GEVNUM_CUT=no to disable.");
        return true;
    }

    std::cout << "\nUse optional T->g.evnum good-event cut from selection report?\n"
              << "Selection report: " << report_csv.Data() << "\n"
              << "Apply cut [Y/n]: " << std::flush;

    std::string reply;
    while (std::getline(std::cin, reply)) {
        reply = nps_to_lower_copy(nps_trim_copy(reply));
        if (reply.empty() || reply == "y" || reply == "yes") {
            logmsg(INFO, "User chose to apply optional T->g.evnum good-event range cut.");
            return true;
        }
        if (reply == "n" || reply == "no") {
            logmsg(INFO, "User chose to skip optional T->g.evnum good-event range cut.");
            return false;
        }
        std::cout << "Please answer y/yes or n/no [Y]: " << std::flush;
    }

    logmsg(WARN, "Input stream closed while confirming T->g.evnum cut usage; defaulting to apply cut.");
    return true;
}

inline bool choose_continue_without_good_event_ranges(const TString& report_csv,
                                                      const std::string& reason)
{
    const char* env_action = gSystem->Getenv("NPS_MISSING_SELECTION_REPORT_ACTION");
    if (env_action && std::strlen(env_action) > 0) {
        const std::string action = nps_to_lower_copy(nps_trim_copy(env_action));
        if (action == "continue" || action == "yes" || action == "y" || action == "1") {
            logmsg(WARN, Form("Good-event range input unavailable (%s): %s. Continuing without event-range cut due to NPS_MISSING_SELECTION_REPORT_ACTION=%s",
                              report_csv.Data(), reason.c_str(), env_action));
            return true;
        }
        if (action == "terminate" || action == "abort" || action == "no" || action == "n" || action == "0") {
            logmsg(ERROR, Form("Good-event range input unavailable (%s): %s. Terminating due to NPS_MISSING_SELECTION_REPORT_ACTION=%s",
                               report_csv.Data(), reason.c_str(), env_action));
            return false;
        }
        logmsg(ERROR, Form("Invalid NPS_MISSING_SELECTION_REPORT_ACTION=%s. Valid values are continue|terminate.", env_action));
        return false;
    }

    if (!is_interactive_terminal()) {
        logmsg(ERROR, Form("Good-event range input unavailable (%s): %s. Non-interactive mode defaults to terminate. Set NPS_MISSING_SELECTION_REPORT_ACTION=continue to proceed without the cut.",
                           report_csv.Data(), reason.c_str()));
        return false;
    }

    std::cout << "\n[ERROR] Good-event range input unavailable: " << report_csv.Data() << "\n"
              << "Reason: " << reason << "\n"
              << "Continue without applying T->g.evnum good-event ranges? [y/N]: " << std::flush;

    std::string reply;
    while (std::getline(std::cin, reply)) {
        reply = nps_to_lower_copy(nps_trim_copy(reply));
        if (reply.empty() || reply == "n" || reply == "no") {
            logmsg(ERROR, "User chose to terminate because good-event range input is unavailable.");
            return false;
        }
        if (reply == "y" || reply == "yes") {
            logmsg(WARN, "User chose to continue without T->g.evnum good-event range cut.");
            return true;
        }
        std::cout << "Please answer y/yes or n/no [N]: " << std::flush;
    }

    logmsg(ERROR, "Input stream closed while prompting for missing good-event range handling. Terminating.");
    return false;
}

inline bool passes_good_event_cut(Double_t event_number,
                                  const std::vector<EventRange>* ranges)
{
    if (!ranges || ranges->empty()) return true;
    if (!std::isfinite(event_number)) return false;

    const Long64_t evnum = static_cast<Long64_t>(std::llround(event_number));
    for (const auto& range : *ranges) {
        if (evnum < range.low) return false;
        if (evnum < range.high) return true;
    }
    return false;
}

struct HelicityInterval {
    Long64_t low = 0;
    Long64_t high = -1;
    int helicity = 0;
};

inline bool build_helicity_intervals_from_file(const TString& file_path,
                                               std::vector<HelicityInterval>& intervals,
                                               std::string& error)
{
    TFile file(file_path.Data(), "READ");
    if (file.IsZombie()) {
        error = Form("could not open ROOT file %s", file_path.Data());
        return false;
    }

    TTree* t = dynamic_cast<TTree*>(file.Get("T"));
    TTree* h = dynamic_cast<TTree*>(file.Get("TSHelH"));
    if (!t || !h) {
        error = Form("missing T or TSHelH tree in %s", file_path.Data());
        return false;
    }
    if (!t->GetBranch("g.evnum") ||
        !h->GetBranch("evcount") ||
        !h->GetBranch("actualHelicity")) {
        error = Form("missing g.evnum, evcount, or actualHelicity branch in %s", file_path.Data());
        return false;
    }

    t->SetBranchStatus("*", 0);
    t->SetBranchStatus("g.evnum", 1);
    double gevnum_raw = 0.0;
    t->SetBranchAddress("g.evnum", &gevnum_raw);

    std::vector<Long64_t> gevnums;
    const Long64_t n_t = t->GetEntries();
    gevnums.reserve(static_cast<size_t>(std::max<Long64_t>(0, n_t)));
    for (Long64_t i = 0; i < n_t; ++i) {
        t->GetEntry(i);
        gevnums.push_back(static_cast<Long64_t>(std::llround(gevnum_raw)));
    }
    if (gevnums.empty()) {
        error = Form("no T entries available for helicity mapping in %s", file_path.Data());
        return false;
    }

    h->SetBranchStatus("*", 0);
    h->SetBranchStatus("evcount", 1);
    h->SetBranchStatus("actualHelicity", 1);
    double evcount_raw = 0.0;
    double helicity_raw = 0.0;
    h->SetBranchAddress("evcount", &evcount_raw);
    h->SetBranchAddress("actualHelicity", &helicity_raw);

    struct ScalerHelRow {
        int evcount = 0;
        int helicity = 0;
    };

    std::vector<ScalerHelRow> rows;
    const Long64_t n_h = h->GetEntries();
    rows.reserve(static_cast<size_t>(std::max<Long64_t>(0, n_h)));
    for (Long64_t i = 0; i < n_h; ++i) {
        h->GetEntry(i);
        const int evcount = static_cast<int>(std::llround(evcount_raw));
        const int helicity = static_cast<int>(std::llround(helicity_raw));
        if (evcount <= 0 || helicity == 0) continue;
        rows.push_back({evcount, helicity > 0 ? 1 : -1});
    }
    if (rows.empty()) {
        error = Form("no nonzero actualHelicity rows found in %s", file_path.Data());
        return false;
    }

    std::sort(rows.begin(), rows.end(), [](const ScalerHelRow& a, const ScalerHelRow& b) {
        return a.evcount < b.evcount;
    });

    const Long64_t max_gevnum = *std::max_element(gevnums.begin(), gevnums.end());
    size_t added = 0;
    for (size_t i = 0; i < rows.size(); ++i) {
        const Long64_t idx = static_cast<Long64_t>(rows[i].evcount) - 1;
        if (idx < 0 || idx >= static_cast<Long64_t>(gevnums.size())) continue;

        const Long64_t low = gevnums[static_cast<size_t>(idx)];
        Long64_t high = max_gevnum;
        for (size_t j = i + 1; j < rows.size(); ++j) {
            const Long64_t next_idx = static_cast<Long64_t>(rows[j].evcount) - 1;
            if (next_idx < 0 || next_idx >= static_cast<Long64_t>(gevnums.size())) continue;
            high = gevnums[static_cast<size_t>(next_idx)] - 1;
            break;
        }
        intervals.push_back({low, std::max(low, high), rows[i].helicity});
        ++added;
    }

    if (added == 0) {
        error = Form("could not map TSHelH evcount rows onto T g.evnum in %s", file_path.Data());
        return false;
    }
    return true;
}

inline bool build_helicity_intervals_from_chain(TChain* chain,
                                                std::vector<HelicityInterval>& intervals,
                                                std::string& error)
{
    intervals.clear();
    if (!chain) {
        error = "null TChain";
        return false;
    }

    TObjArray* files = chain->GetListOfFiles();
    if (!files || files->GetEntries() <= 0) {
        error = "TChain contains no source files";
        return false;
    }

    int ok_files = 0;
    std::string first_error;
    for (int i = 0; i < files->GetEntries(); ++i) {
        TObject* obj = files->At(i);
        if (!obj) continue;
        TString path = obj->GetTitle();
        if (path.IsNull()) path = obj->GetName();
        if (path.IsNull()) continue;

        std::string file_error;
        if (build_helicity_intervals_from_file(path, intervals, file_error)) {
            ++ok_files;
        } else if (first_error.empty()) {
            first_error = file_error;
        }
    }

    if (intervals.empty()) {
        error = first_error.empty() ? "no helicity intervals built from chain" : first_error;
        return false;
    }

    std::sort(intervals.begin(), intervals.end(), [](const HelicityInterval& a, const HelicityInterval& b) {
        if (a.low != b.low) return a.low < b.low;
        return a.high < b.high;
    });
    logmsg(INFO, Form("Built %zu helicity interval(s) from %d source file(s)", intervals.size(), ok_files));
    return true;
}

inline int helicity_for_event(Double_t event_number, const std::vector<HelicityInterval>& intervals)
{
    if (intervals.empty() || !std::isfinite(event_number)) return 0;
    const Long64_t evnum = static_cast<Long64_t>(std::llround(event_number));

    auto it = std::upper_bound(intervals.begin(), intervals.end(), evnum,
        [](Long64_t value, const HelicityInterval& interval) {
            return value < interval.low;
        });
    if (it == intervals.begin()) return 0;
    --it;
    if (evnum >= it->low && evnum <= it->high) return it->helicity;
    return 0;
}

inline bool parse_env_bool(const char* value, bool default_value)
{
    if (!value || std::strlen(value) == 0) return default_value;
    const std::string v = nps_to_lower_copy(nps_trim_copy(value));
    if (v == "1" || v == "true" || v == "yes" || v == "y" || v == "on") return true;
    if (v == "0" || v == "false" || v == "no" || v == "n" || v == "off") return false;
    return default_value;
}

inline bool open_chain_for_run(const TString& input_dir,
                               int run,
                               AnalysisMode requested_mode,
                               std::unique_ptr<TChain>& chain,
                               TString& tree_name_used,
                               int& nfiles_added,
                               std::string& error)
{
    std::vector<TString> patterns;
    if (requested_mode == AnalysisMode::kAuto || requested_mode == AnalysisMode::kWaveform) {
        patterns.push_back(input_dir + Form("nps_production_%d_*_wf_calib.root", run));
    }
    if (requested_mode == AnalysisMode::kAuto || requested_mode == AnalysisMode::kHcana) {
        patterns.push_back(input_dir + Form("skim_run%d.root", run));
        patterns.push_back(input_dir + Form("nps_hms_coin_%d_*_1_-1.root", run));
    }

    std::vector<std::string> tree_candidates;
    if (requested_mode == AnalysisMode::kWaveform) {
        tree_candidates = {"t_prod", "T"};
    } else if (requested_mode == AnalysisMode::kHcana) {
        tree_candidates = {"T", "t_prod"};
    } else {
        tree_candidates = {"t_prod", "T"};
    }

    for (const auto& tree_name : tree_candidates) {
        std::unique_ptr<TChain> trial(new TChain(tree_name.c_str()));
        int added = 0;
        for (const auto& pat : patterns) {
            added += trial->Add(pat);
        }
        if (added <= 0) continue;
        if (trial->GetEntries() <= 0) continue;

        chain = std::move(trial);
        tree_name_used = tree_name.c_str();
        nfiles_added = added;
        return true;
    }

    error = Form("no ROOT files found for run %d under %s", run, input_dir.Data());
    return false;
}

inline bool resolve_run_kinematics(const RunMeta& run_meta,
                                   double beam_energy_override,
                                   double& run_beam_energy,
                                   double& run_nps_theta_deg,
                                   double& run_z_nps_cm,
                                   std::string& error)
{
    run_beam_energy = (beam_energy_override > 0.0) ? beam_energy_override : run_meta.beam_energy;
    run_nps_theta_deg = run_meta.nps_theta_deg;
    run_z_nps_cm = run_meta.nps_target_distance_cm;

    if (!std::isfinite(run_beam_energy) || run_beam_energy <= 0.0) {
        error = Form("invalid beam energy %.6g", run_beam_energy);
        return false;
    }
    if (!std::isfinite(run_z_nps_cm) || run_z_nps_cm <= 0.0) {
        error = Form("invalid NPS_target_distance %.6g cm", run_z_nps_cm);
        return false;
    }
    return true;
}

inline int load_clusters_from_branches(const std::vector<double>* clusE_vec,
                                       const std::vector<double>* clusX_vec,
                                       const std::vector<double>* clusY_vec,
                                       const std::vector<double>* clusT_vec,
                                       Double_t nclust_dbl,
                                       Double_t clusE[],
                                       Double_t clusX[],
                                       Double_t clusY[],
                                       Double_t clusT[],
                                       double cluster_time_offset_ns)
{
    int nclust = 0;
    if (clusE_vec && clusX_vec && clusY_vec && clusT_vec) {
        const size_t nvec = std::min({
            clusE_vec->size(),
            clusX_vec->size(),
            clusY_vec->size(),
            clusT_vec->size(),
            static_cast<size_t>(MAX_CLUS)
        });
        nclust = static_cast<int>(nvec);
        for (int i = 0; i < nclust; ++i) {
            clusE[i] = (*clusE_vec)[i];
            clusX[i] = (*clusX_vec)[i];
            clusY[i] = (*clusY_vec)[i];
            clusT[i] = (*clusT_vec)[i] + cluster_time_offset_ns;
        }
    } else {
        nclust = static_cast<int>(std::lround(nclust_dbl));
        if (nclust < 0) nclust = 0;
        if (nclust > MAX_CLUS) nclust = MAX_CLUS;
    }
    return nclust;
}

inline void collect_good_clusters(const AcceptanceCuts& acceptance_cuts,
                                  const Double_t clusE[],
                                  const Double_t clusX[],
                                  const Double_t clusY[],
                                  const Double_t clusT[],
                                  int n_after,
                                  double cluster_time_halfwidth_ns,
                                  std::vector<int>& good_idx,
                                  const nps_dead_block::Mask* dead_block_mask = nullptr)
{
    good_idx.clear();
    good_idx.reserve(8);
    for (int i = 0; i < n_after; ++i) {
        if (dead_block_mask && dead_block_mask->rejects_xy(clusX[i], clusY[i])) {
            continue;
        }
        if (acceptance_cuts.pass_nps_cluster(clusE[i], clusX[i], clusY[i], clusT[i], cluster_time_halfwidth_ns)) {
            good_idx.push_back(i);
        }
    }
    if (good_idx.size() > 4) {
        std::sort(good_idx.begin(), good_idx.end(), [&](int a, int b) { return clusE[a] > clusE[b]; });
        good_idx.resize(4);
    }
}

inline bool choose_pi0_pair(const std::vector<int>& good_idx,
                            Double_t clusE[],
                            Double_t clusX[],
                            Double_t clusY[],
                            Double_t clusT[],
                            double run_z_nps_cm,
                            double pair_time_diff_max_ns,
                            int& sel_i,
                            int& sel_j)
{
    sel_i = -1;
    sel_j = -1;
    if (good_idx.size() < 2) return false;

    if (good_idx.size() == 2) {
        sel_i = good_idx[0];
        sel_j = good_idx[1];
        return true;
    }

    auto pr = nps::choose_best_pair_closest_pi0(good_idx,
                                                 clusE,
                                                 clusX,
                                                 clusY,
                                                 clusT,
                                                 run_z_nps_cm,
                                                 nps::kPi0Mass_GeV,
                                                 pair_time_diff_max_ns);
    sel_i = pr.first;
    sel_j = pr.second;
    return (sel_i >= 0 && sel_j >= 0);
}

// ============================================================================
// Main Analysis Macro
// ============================================================================
void nps_analysis_main(const TString &kinematic_in = "",
                       const TString &input_dir_in = "",
                       const TString &mode_in = "auto",
                       const TString &config_csv_in = MAIN_CONFIG_CSV,
                       const TString &output_base_in = CANONICAL_OUTPUT_BASE,
                       const TString &types_in = "production,Production",
                       const double Ebeam_override_in = -1.0)
{
    // Timing and logging
    TStopwatch sw_total; 
    sw_total.Start();
    logmsg(INFO, "=========== NPS π0 unified analysis ===========");
    
    // Optimization: Set ROOT to batch mode for faster canvas operations
    gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);

    TString kinematic = kinematic_in;
    TString input_dir = input_dir_in;
    TString mode_str = mode_in;
    TString config_csv = config_csv_in;
    TString acceptance_cuts_config = ACCEPTANCE_CUTS_CONFIG_DEFAULT;
    TString dead_block_config = DEAD_BLOCK_CONFIG_DEFAULT;
    TString output_base = output_base_in;
    TString types = types_in;
    double beam_energy_override = Ebeam_override_in;

    const char* env_kin = gSystem->Getenv("NPS_KIN");
    const char* env_input = gSystem->Getenv("NPS_INPUT_DIR");
    const char* env_input_legacy = gSystem->Getenv("NPS_SKIM_DIR");
    const char* env_mode = gSystem->Getenv("NPS_MODE");
    const char* env_config = gSystem->Getenv("NPS_CONFIG_CSV");
    const char* env_output_base = gSystem->Getenv("NPS_OUTPUT_BASE");
    const char* env_types = gSystem->Getenv("NPS_TYPES");
    const char* env_run = gSystem->Getenv("NPS_RUN");
    const char* env_ebeam = gSystem->Getenv("NPS_EBEAM");
    const char* env_acceptance_cuts = gSystem->Getenv("NPS_ACCEPTANCE_CUTS_CONFIG");
    const char* env_selection_report = gSystem->Getenv("NPS_SELECTION_REPORT_CSV");
    const char* env_dead_block_config = gSystem->Getenv("NPS_DEAD_BLOCK_CONFIG");

    if (env_kin) kinematic = env_kin;
    if (env_input) input_dir = env_input;
    if (input_dir.IsNull() && env_input_legacy) input_dir = env_input_legacy;
    if (env_mode) mode_str = env_mode;
    if (env_config) config_csv = env_config;
    if (env_acceptance_cuts && std::strlen(env_acceptance_cuts) > 0) acceptance_cuts_config = env_acceptance_cuts;
    if (env_dead_block_config && std::strlen(env_dead_block_config) > 0) dead_block_config = env_dead_block_config;
    if (env_output_base) output_base = env_output_base;
    if (env_types) types = env_types;
    if (env_ebeam) beam_energy_override = atof(env_ebeam);

    if (kinematic.IsNull() || TString(kinematic).Strip().Length() == 0) {
        logmsg(ERROR, "NPS_KIN (Kin_old) is required for nps_analysis_main");
        return;
    }

    const AnalysisMode requested_mode = parse_mode_string(mode_str.Data());
    if (input_dir.IsNull() || TString(input_dir).Strip().Length() == 0) {
        input_dir = (requested_mode == AnalysisMode::kHcana) ? DEFAULT_HCANA_INPUT_DIR : DEFAULT_WAVEFORM_INPUT_DIR;
    }

    input_dir = input_dir.EndsWith("/") ? input_dir : input_dir + "/";
    output_base = output_base.EndsWith("/") ? output_base : output_base + "/";

    int run_filter = -1;
    if (env_run && std::strlen(env_run) > 0) {
        run_filter = std::atoi(env_run);
    }

    std::vector<RunMeta> run_meta_rows;
    std::string cfg_error;
    if (!load_config_rows(config_csv, kinematic, parse_type_filter(types.Data()), run_filter, run_meta_rows, cfg_error)) {
        logmsg(ERROR, Form("Config selection failed: %s", cfg_error.c_str()));
        return;
    }
    fill_missing_charge_from_fallback(run_meta_rows, CHARGE_FALLBACK_DB);

    std::unordered_map<int, RunMeta> run_meta_by_run;
    std::vector<int> runs;
    runs.reserve(run_meta_rows.size());
    for (const auto& m : run_meta_rows) {
        run_meta_by_run[m.run] = m;
        runs.push_back(m.run);
    }
    logmsg(INFO, Form("Selected %zu runs for kin=%s from %s", runs.size(), kinematic.Data(), config_csv.Data()));

    OutputFolders out_dirs;
    std::string out_error;
    if (!ensure_output_folders(output_base, kinematic, out_dirs, out_error)) {
        logmsg(ERROR, Form("Output directory setup failed: %s", out_error.c_str()));
        return;
    }

    if (gSystem->AccessPathName(acceptance_cuts_config) != 0) {
        logmsg(ERROR, Form("Acceptance-cuts config not found: %s", acceptance_cuts_config.Data()));
        return;
    }

    std::unique_ptr<AcceptanceCuts> acceptance_cuts;
    try {
        acceptance_cuts.reset(new AcceptanceCuts(acceptance_cuts_config.Data(), kinematic.Data()));
    } catch (const std::exception& ex) {
        logmsg(ERROR, Form("Failed to load acceptance cuts: %s", ex.what()));
        return;
    }

    logmsg(INFO, Form("Using acceptance-cuts config: %s", acceptance_cuts_config.Data()));
    std::cout << acceptance_cuts->summary() << std::endl;

    std::vector<nps_dead_block::RunRange> dead_block_ranges;
    std::string dead_block_error;
    if (!nps_dead_block::load_run_ranges(dead_block_config.Data(), dead_block_ranges, dead_block_error)) {
        logmsg(WARN, Form("Dead-block masking disabled: %s", dead_block_error.c_str()));
    } else {
        logmsg(INFO, Form("Loaded %zu dead-block run range(s) from %s",
                          dead_block_ranges.size(), dead_block_config.Data()));
    }

    TString outPlotDir = out_dirs.plots.EndsWith("/") ? out_dirs.plots : out_dirs.plots + "/";
    TString outRootDir = out_dirs.root.EndsWith("/") ? out_dirs.root : out_dirs.root + "/";
    TString outSummaryDir = out_dirs.summary.EndsWith("/") ? out_dirs.summary : out_dirs.summary + "/";

    TString selection_report_csv;
    if (env_selection_report && std::strlen(env_selection_report) > 0) {
        selection_report_csv = TString(env_selection_report);
    } else {
        const TString direct_name = TString::Format("%s/selection_report_%s.csv", EFFICIENCY_STUFF_DIR, kinematic.Data());
        const TString safe_name = TString::Format("%s/selection_report_%s.csv", EFFICIENCY_STUFF_DIR, sanitize_token(kinematic.Data()).c_str());
        selection_report_csv = (gSystem->AccessPathName(direct_name) == 0) ? direct_name : safe_name;
    }
    std::unordered_map<int, std::vector<EventRange>> good_event_ranges_by_run;
    std::unordered_map<int, GoodEventCharge> good_event_charge_by_run;
    bool apply_good_event_cut = false;
    std::string selection_error;
    const bool user_wants_good_event_cut = choose_use_good_event_cut(selection_report_csv);
    if (user_wants_good_event_cut) {
        if (!load_good_event_ranges_by_run(selection_report_csv, kinematic, good_event_ranges_by_run, good_event_charge_by_run, selection_error)) {
            if (!choose_continue_without_good_event_ranges(selection_report_csv, selection_error)) {
                logmsg(ERROR, "Terminating due to unavailable good-event range input.");
                return;
            }
            logmsg(WARN, "Continuing without applying T->g.evnum good-event range cut.");
        } else {
            apply_good_event_cut = true;
            logmsg(INFO, Form("Loaded T->g.evnum good-event ranges for %zu run(s) from %s",
                              good_event_ranges_by_run.size(), selection_report_csv.Data()));
        }
    } else {
        logmsg(INFO, "Continuing without applying optional T->g.evnum good-event range cut.");
    }

    const bool write_global_csv = parse_env_bool(gSystem->Getenv("NPS_WRITE_GLOBAL_CSV"), false);
    TString global_csv;
    if (write_global_csv) {
        global_csv = outSummaryDir + "summary_all_runs.csv";
        if (!write_global_csv_header(global_csv)) {
            TString fallback_csv = Form("/tmp/summary_all_runs_%d.csv", (int)time(nullptr));
            logmsg(WARN, Form("Falling back global CSV path to %s", fallback_csv.Data()));
            if (!write_global_csv_header(fallback_csv)) {
                logmsg(ERROR, "Failed to create global summary CSV in both output directory and /tmp");
                return;
            }
            global_csv = fallback_csv;
        }
    } else {
        logmsg(INFO, "Global summary CSV append is disabled (NPS_WRITE_GLOBAL_CSV=0). Parallel driver rebuilds summary from per-run rows.");
    }

    // ============================================================================
    // Per-run loop (process each run independently with memory cleanup)
    // ============================================================================
    for (int run : runs) {
        TStopwatch sw_run;
        sw_run.Start();

        auto run_meta_it = run_meta_by_run.find(run);
        if (run_meta_it == run_meta_by_run.end()) {
            logmsg(ERROR, Form("Skipping run %d: metadata row not found after config selection", run));
            continue;
        }
        const RunMeta& run_meta = run_meta_it->second;

        double run_beam_energy = std::numeric_limits<double>::quiet_NaN();
        double run_nps_theta_deg = std::numeric_limits<double>::quiet_NaN();
        double run_z_nps_cm = std::numeric_limits<double>::quiet_NaN();
        std::string run_kinematics_error;
        if (!resolve_run_kinematics(run_meta,
                                    beam_energy_override,
                                    run_beam_energy,
                                    run_nps_theta_deg,
                                    run_z_nps_cm,
                                    run_kinematics_error)) {
            logmsg(ERROR, Form("Skipping run %d: %s", run, run_kinematics_error.c_str()));
            continue;
        }

        logmsg(INFO, Form("Processing run %d [mode request=%s, kin=%s, type=%s, E=%.4f GeV, NPS_Thet=%.3f deg, z=%.2f cm]",
                          run,
                          mode_to_cstr(requested_mode),
                          run_meta.kin.Data(),
                          run_meta.type.Data(),
                          run_beam_energy,
                          run_nps_theta_deg,
                          run_z_nps_cm));

        const std::vector<EventRange>* good_ranges_for_run = nullptr;
        const GoodEventCharge* good_charge_for_run = nullptr;
        if (apply_good_event_cut) {
            auto gr_it = good_event_ranges_by_run.find(run);
            if (gr_it != good_event_ranges_by_run.end() && !gr_it->second.empty()) {
                good_ranges_for_run = &gr_it->second;
                auto gc_it = good_event_charge_by_run.find(run);
                if (gc_it != good_event_charge_by_run.end()) {
                    good_charge_for_run = &gc_it->second;
                }
                logmsg(INFO, Form("Run %d: applying good-event cut with %zu T->g.evnum range(s)",
                                  run, gr_it->second.size()));
            } else {
                logmsg(WARN, Form("Run %d: no T->g.evnum good ranges found in %s for this run; proceeding without event-number cut for this run",
                                  run, selection_report_csv.Data()));
            }
        }

        // ====================================================================
        // Per-run processing with RAII file management and exception safety
        // ====================================================================
        std::string run_status = "OK";  // Track whether run completes without errors
        try {
            std::unique_ptr<TChain> chain;
            TString tree_name_used;
            int nfiles_added = 0;
            std::string open_error;
            if (!open_chain_for_run(input_dir, run, requested_mode, chain, tree_name_used, nfiles_added, open_error)) {
                logmsg(WARN, Form("Skipping run %d: %s", run, open_error.c_str()));
                continue;
            }

            TTree *T = chain.get();  // Use the chain as a TTree
            if (!T) {
                logmsg(ERROR, Form("Chain is null for run %d", run));
                continue;
            }

            BranchMapping branch_map;
            std::string map_error;
            if (!detect_branch_mapping(T, requested_mode, branch_map, map_error)) {
                logmsg(ERROR, Form("Skipping run %d: %s", run, map_error.c_str()));
                continue;
            }

            const double cluster_time_halfwidth_ns =
                acceptance_cuts->resolved_cluster_time_halfwidth_ns(branch_map.time_window_wrt_150_ns);
            const double pair_time_diff_max_ns =
                acceptance_cuts->resolved_pair_time_diff_max_ns(branch_map.time_threshold_ns);
            const bool use_shifted_timing_windows =
                acceptance_cuts->resolved_use_shifted_sidebands(branch_map.use_shifted_timing_windows);

            if (!(cluster_time_halfwidth_ns > 0.0) || !(pair_time_diff_max_ns > 0.0)) {
                logmsg(ERROR, Form("Skipping run %d: resolved timing cuts are invalid (cluster_halfwidth=%.3f, pair_diff_max=%.3f)",
                                  run, cluster_time_halfwidth_ns, pair_time_diff_max_ns));
                continue;
            }

            logmsg(INFO, Form("Run %d: Added %d file(s), tree=%s, detected mode=%s",
                              run, nfiles_added, tree_name_used.Data(), branch_map.mode_label.Data()));
            logmsg(INFO, Form("Run %d: resolved timing cuts (cluster_halfwidth=%.3f ns, pair_diff_max=%.3f ns, shifted_sidebands=%s)",
                              run, cluster_time_halfwidth_ns, pair_time_diff_max_ns,
                              use_shifted_timing_windows ? "true" : "false"));

            Long64_t nentries = T->GetEntries();
            logmsg(INFO, Form("Run %d: %lld entries", run, nentries));

            // -------------------------
            // Branch variables (no static arrays)
            // -------------------------
            Double_t HgtrX=0, HgtrY=0, HgtrTh=0, HgtrPh=0, hdelta=0, HgtrP=0;
            Double_t hreactx=0, hreacty=0, hreactz=0;
            Double_t hfastRasterXA=0, hfastRasterXB=0, hfastRasterYA=0, hfastRasterYB=0;
            Double_t hcernpeSum=0, hcaletotnorm=0;
            Double_t HgtrPx=0, HgtrPy=0, HgtrPz=0;
            Double_t edtmtdc=0;
            Double_t g_evnum=std::numeric_limits<double>::quiet_NaN();
            Double_t t_helicity_value=0;
            Double_t nclust_dbl = 0;
            Double_t BCM4A_scalerCurrent = 0, BCM4A_scalerCharge = 0, H_1MHz_scaler = 0;
            
            // HMS focal plane variables
            Double_t hxfp=0, hyfp=0, hxpfp=0, hypfp=0;

            // ROOT fills variable-length arrays before analysis clips to MAX_CLUS.
            Double_t clusE[MAX_CLUS_INPUT] = {0}, clusX[MAX_CLUS_INPUT] = {0};
            Double_t clusY[MAX_CLUS_INPUT] = {0}, clusT[MAX_CLUS_INPUT] = {0};

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

            auto has_branch = [&](const char* b) {
                return T->GetBranch(b) != nullptr;
            };

            const bool has_react_x = has_branch("H.react.x");
            const bool has_react_y = has_branch("H.react.y");
            const bool has_fast_raster_xa = has_branch("H.rb.raster.fr_xa");
            const bool has_fast_raster_xb = has_branch("H.rb.raster.fr_xb");
            const bool has_fast_raster_ya = has_branch("H.rb.raster.fr_ya");
            const bool has_fast_raster_yb = has_branch("H.rb.raster.fr_yb");
            const bool has_gtr_x = has_branch("H.gtr.x");
            const bool has_gtr_y = has_branch("H.gtr.y");
            const bool has_fp_x = has_branch("H.dc.x_fp");
            const bool has_fp_y = has_branch("H.dc.y_fp");
            const bool has_fp_xp = has_branch("H.dc.xp_fp");
            const bool has_fp_yp = has_branch("H.dc.yp_fp");

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
            enable_scalar("H.react.x", &hreactx);
            enable_scalar("H.react.y", &hreacty);
            enable_scalar("H.react.z", &hreactz);
            enable_scalar("H.rb.raster.fr_xa", &hfastRasterXA);
            enable_scalar("H.rb.raster.fr_xb", &hfastRasterXB);
            enable_scalar("H.rb.raster.fr_ya", &hfastRasterYA);
            enable_scalar("H.rb.raster.fr_yb", &hfastRasterYB);
            enable_scalar("H.cer.npeSum", &hcernpeSum);
            enable_scalar("H.cal.etotnorm", &hcaletotnorm);
            enable_scalar("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);
            enable_scalar("H.BCM4A.scalerCurrent", &BCM4A_scalerCurrent);
            enable_scalar("H.BCM4A.scalerCharge", &BCM4A_scalerCharge);
            enable_scalar("H.1MHz.scaler", &H_1MHz_scaler);
            enable_scalar(branch_map.nclust_branch.Data(), &nclust_dbl);
            const TString helicity_branch =
                (branch_map.mode == AnalysisMode::kWaveform) ? "T.helicity.hel" : "T.helicity.mps";
            const bool has_t_helicity_branch = (T->GetBranch(helicity_branch.Data()) != nullptr);
            if (has_t_helicity_branch) {
                enable_scalar(helicity_branch.Data(), &t_helicity_value);
                logmsg(INFO, Form("Run %d: using %s for output helicity", run, helicity_branch.Data()));
            } else {
                logmsg(WARN, Form("Run %d: %s branch missing; falling back to TSHelH actualHelicity mapping",
                                  run, helicity_branch.Data()));
            }

            bool apply_good_event_cut_for_run = (good_ranges_for_run != nullptr);
            if (apply_good_event_cut_for_run) {
                if (!T->GetBranch("g.evnum")) {
                    const std::string missing_reason = Form("run %d does not contain required branch g.evnum", run);
                    if (choose_continue_without_good_event_ranges(selection_report_csv, missing_reason)) {
                        logmsg(WARN, Form("Run %d: proceeding without event-number cut because g.evnum branch is missing", run));
                        apply_good_event_cut_for_run = false;
                    } else {
                        logmsg(ERROR, Form("Run %d: terminating due to missing g.evnum branch", run));
                        return;
                    }
                } else {
                    enable_scalar("g.evnum", &g_evnum);
                }
            }
            
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

            try_bind_cluster(branch_map.clusE_branch.Data(), &clusE_vec, clusE);
            try_bind_cluster(branch_map.clusX_branch.Data(), &clusX_vec, clusX);
            try_bind_cluster(branch_map.clusY_branch.Data(), &clusY_vec, clusY);
            try_bind_cluster(branch_map.clusT_branch.Data(), &clusT_vec, clusT);

            if (!T->GetBranch(branch_map.nclust_branch.Data()) ||
                !T->GetBranch(branch_map.clusE_branch.Data()) ||
                !T->GetBranch(branch_map.clusX_branch.Data()) ||
                !T->GetBranch(branch_map.clusY_branch.Data()) ||
                !T->GetBranch(branch_map.clusT_branch.Data())) {
                logmsg(ERROR, Form("Skipping run %d: required cluster branches are missing for mode %s",
                                   run, branch_map.mode_label.Data()));
                continue;
            }

            std::vector<HelicityInterval> helicity_intervals;
            if (!has_t_helicity_branch) {
                std::string helicity_error;
                if (!build_helicity_intervals_from_chain(chain.get(), helicity_intervals, helicity_error)) {
                    logmsg(WARN, Form("Run %d: helicity branch will be filled with 0 because mapping failed: %s",
                                      run, helicity_error.c_str()));
                }
            }

            cout << "Run " << run << " entries: " << nentries << endl;

            double run_current_mean = run_meta.current_uA;
            double accumulated_charge_mC = run_meta.charge_mC;
            double hel_pos_charge_mC = 0.0;
            double hel_neg_charge_mC = 0.0;
            double beam_time = run_meta.beam_time_s;

            if (good_ranges_for_run != nullptr) {
                if (good_charge_for_run == nullptr || !(good_charge_for_run->total_mC() > 0.0)) {
                    logmsg(ERROR, Form("Skipping run %d: missing/invalid helicity good-event charge from selection report", run));
                    continue;
                }
                hel_pos_charge_mC = good_charge_for_run->hel_pos_mC;
                hel_neg_charge_mC = good_charge_for_run->hel_neg_mC;
                accumulated_charge_mC = good_charge_for_run->total_mC();
                logmsg(INFO, Form("Run %d: using helicity good-event charge %.6f mC (HEL+=%.6f, HEL-=%.6f)",
                                  run, accumulated_charge_mC, hel_pos_charge_mC, hel_neg_charge_mC));
            }

            if (!std::isfinite(accumulated_charge_mC) || accumulated_charge_mC <= 0.0) {
                logmsg(ERROR, Form("Skipping run %d: missing/invalid charge metadata", run));
                continue;
            }
            if (!std::isfinite(beam_time) || beam_time <= 0.0) {
                logmsg(ERROR, Form("Skipping run %d: missing/invalid beam-time metadata", run));
                continue;
            }
            if (!std::isfinite(run_current_mean) || run_current_mean <= 0.0) {
                logmsg(WARN, Form("Run %d: current metadata unavailable; setting current_mean_uA to 0.0", run));
                run_current_mean = 0.0;
            }

            cout << "Charge = " << accumulated_charge_mC << " mC";
            if (hel_pos_charge_mC > 0.0 || hel_neg_charge_mC > 0.0) {
                cout << " (HEL+=" << hel_pos_charge_mC << " mC, HEL-=" << hel_neg_charge_mC << " mC)";
            }
            cout << "\n";
            cout << "Beam_Time = " << beam_time << " s\n";

            const nps_dead_block::Mask dead_block_mask =
                nps_dead_block::mask_for_run(dead_block_ranges, run);
            if (!dead_block_mask.empty()) {
                logmsg(INFO, Form("Run %d: dead blocks [%s], excluding %zu block(s) in 3x3 mask(s)",
                                  run,
                                  nps_dead_block::join_blocks(dead_block_mask.dead_blocks).c_str(),
                                  dead_block_mask.excluded_blocks.size()));
            }

        // -------------------------
        // Histograms (unique names per run). SetDirectory(nullptr) to avoid ROOT file ownership issues
        // -------------------------
        auto name = [&](const char* base){ return TString::Format("%s_run%d", base, run); };

        auto style1D = [&](TH1D* h) {
        if (!h) return;
        h->SetLineWidth(2);
        h->SetLineColor(kBlue + 1);
        h->SetMarkerStyle(20);
        h->SetMarkerSize(0.70);
        h->SetMarkerColor(kBlue + 1);
        h->GetXaxis()->CenterTitle(true);
        h->GetYaxis()->CenterTitle(true);
        h->GetXaxis()->SetTitleSize(0.045);
        h->GetYaxis()->SetTitleSize(0.045);
        h->GetXaxis()->SetLabelSize(0.040);
        h->GetYaxis()->SetLabelSize(0.040);
        const TString y_title = h->GetYaxis()->GetTitle();
        if (y_title == "Counts" || y_title == "Events" || y_title == "Entries") {
            h->GetYaxis()->SetTitle("Events / bin");
        }
        h->SetStats(0);
        };

        auto style2D = [&](TH2D* h) {
        if (!h) return;
        h->SetStats(0);
        h->GetXaxis()->CenterTitle(true);
        h->GetYaxis()->CenterTitle(true);
        h->GetXaxis()->SetTitleSize(0.045);
        h->GetYaxis()->SetTitleSize(0.045);
        h->GetZaxis()->SetTitleSize(0.040);
        h->GetXaxis()->SetLabelSize(0.040);
        h->GetYaxis()->SetLabelSize(0.040);
        h->GetZaxis()->SetLabelSize(0.035);
        };

        auto style_canvas = [&](TCanvas* c, bool with_colz = false) {
        if (!c) return;
        c->SetLeftMargin(0.12);
        c->SetBottomMargin(0.12);
        c->SetTopMargin(0.08);
        c->SetRightMargin(with_colz ? 0.14 : 0.05);
        c->SetTicks(1, 1);
        };

        auto set_overlay_y_range = [&](const vector<TH1*>& hists,
                                       double x_lo = std::numeric_limits<double>::quiet_NaN(),
                                       double x_hi = std::numeric_limits<double>::quiet_NaN()) {
        const bool use_x_window = std::isfinite(x_lo) && std::isfinite(x_hi) && (x_hi > x_lo);
        double ymin = std::numeric_limits<double>::infinity();
        double ymax = -std::numeric_limits<double>::infinity();
        for (auto *h : hists) {
            if (!h) continue;
            int b_start = 1;
            int b_end = h->GetNbinsX();
            if (use_x_window) {
                b_start = std::max(1, h->FindBin(x_lo));
                b_end = std::min(h->GetNbinsX(), h->FindBin(x_hi));
            }
            for (int b = b_start; b <= b_end; ++b) {
                const double y = h->GetBinContent(b);
                if (!std::isfinite(y)) continue;
                ymin = std::min(ymin, y);
                ymax = std::max(ymax, y);
            }
        }
        if (!std::isfinite(ymin) || !std::isfinite(ymax)) {
            ymin = 0.0;
            ymax = 1.0;
        }

        double span = ymax - ymin;
        if (!(span > 0.0)) {
            span = (std::fabs(ymax) > 0.0) ? std::fabs(ymax) : 1.0;
        }
        const double pad = 0.12 * span;
        const double disp_min = (ymin < 0.0) ? (ymin - pad) : 0.0;
        const double disp_max = ymax + pad;

        for (auto *h : hists) {
            if (!h) continue;
            h->SetMinimum(disp_min);
            h->SetMaximum((disp_max > disp_min) ? disp_max : (disp_min + 1.0));
            h->GetXaxis()->SetTitleSize(0.045);
            h->GetYaxis()->SetTitleSize(0.045);
            h->GetXaxis()->SetLabelSize(0.040);
            h->GetYaxis()->SetLabelSize(0.040);
        }
        };

        auto make1D = [&](const char* base, const char* title, int nb, double xlo, double xhi)->TH1D* {
        TH1D *h = new TH1D(name(base), title, nb, xlo, xhi);
        h->SetDirectory(nullptr);
        style1D(h);
        return h;
        };
        auto make2D = [&](const char* base, const char* title, int nbx, double xlo, double xhi, int nby, double ylo, double yhi)->TH2D* {
        TH2D *h = new TH2D(name(base), title, nbx, xlo, xhi, nby, ylo, yhi);
        h->SetDirectory(nullptr);
        style2D(h);
        return h;
        };

        struct CutDebugPlot {
            TH1D* hist = nullptr;
            TString label;
            bool has_min = false;
            double min_value = 0.0;
            bool has_max = false;
            double max_value = 0.0;
            bool available = true;
        };

        struct Debug2DPlot {
            TH2D* hist = nullptr;
            TString label;
            bool available = true;
        };

        auto makeCut1D = [&](const char* base, const char* title, int nb, double xlo, double xhi)->TH1D* {
        TH1D *h = make1D(base, title, nb, xlo, xhi);
        h->SetLineColor(kBlack);
        h->SetMarkerColor(kBlack);
        h->SetFillColor(kAzure - 9);
        h->SetFillStyle(1001);
        return h;
        };

        auto fill_cut_hist = [](TH1D* h, double value, double weight = 1.0) {
        if (h && std::isfinite(value) && std::isfinite(weight)) h->Fill(value, weight);
        };

        auto fill_debug_hist_2d = [](TH2D* h, double x, double y, double weight = 1.0) {
        if (h && std::isfinite(x) && std::isfinite(y) && std::isfinite(weight)) h->Fill(x, y, weight);
        };

        const HMSDataCuts& hms_debug_cuts = acceptance_cuts->hms_data();
        const NPSClusterCuts& nps_debug_cuts = acceptance_cuts->nps_cluster();
        const WeightedCuts& weighted_debug_cuts = acceptance_cuts->weighted();

        TH1D *h_cut_edtm_tdc = makeCut1D("h_cut_edtm_tdc", "HMS EDTM TDC;EDTM TDC;Events", 220, -5.0, 5.0);
        TH1D *h_cut_react_x = makeCut1D("h_cut_react_x", "HMS reaction x;hreactx [cm];Events", 220, -2.0, 2.0);
        TH1D *h_cut_react_y = makeCut1D("h_cut_react_y", "HMS reaction y;hreacty [cm];Events", 220, -2.0, 2.0);
        TH1D *h_cut_react_z = makeCut1D("h_cut_react_z", "HMS reaction z;hreactz [cm];Events", 220, -12.0, 12.0);
        TH1D *h_cut_fast_raster_xa = makeCut1D("h_cut_fast_raster_xa", "Fast raster xA;fr_{xA} [cm];Events", 220, -0.5, 0.5);
        TH1D *h_cut_fast_raster_xb = makeCut1D("h_cut_fast_raster_xb", "Fast raster xB;fr_{xB} [cm];Events", 220, -0.5, 0.5);
        TH1D *h_cut_fast_raster_ya = makeCut1D("h_cut_fast_raster_ya", "Fast raster yA;fr_{yA} [cm];Events", 220, -0.5, 0.5);
        TH1D *h_cut_fast_raster_yb = makeCut1D("h_cut_fast_raster_yb", "Fast raster yB;fr_{yB} [cm];Events", 220, -0.5, 0.5);
        TH1D *h_cut_fp_x = makeCut1D("h_cut_fp_x", "HMS focal-plane x;x_{fp} [cm];Events", 220, -60.0, 60.0);
        TH1D *h_cut_fp_y = makeCut1D("h_cut_fp_y", "HMS focal-plane y;y_{fp} [cm];Events", 220, -40.0, 40.0);
        TH1D *h_cut_fp_xp = makeCut1D("h_cut_fp_xp", "HMS focal-plane x';x'_{fp} [rad];Events", 220, -0.15, 0.15);
        TH1D *h_cut_fp_yp = makeCut1D("h_cut_fp_yp", "HMS focal-plane y';y'_{fp} [rad];Events", 220, -0.08, 0.08);
        TH1D *h_cut_gtr_x = makeCut1D("h_cut_gtr_x", "HMS target x;h.gtr.x [cm];Events", 220, -5.0, 5.0);
        TH1D *h_cut_gtr_y = makeCut1D("h_cut_gtr_y", "HMS target y;h.gtr.y [cm];Events", 220, -5.0, 5.0);
        TH1D *h_cut_delta = makeCut1D("h_cut_delta", "HMS #delta;#delta [%];Events", 220, -25.0, 25.0);
        TH1D *h_cut_gtr_th = makeCut1D("h_cut_gtr_th", "HMS x' target;h.gtr.th [rad];Events", 220, -0.16, 0.16);
        TH1D *h_cut_gtr_ph = makeCut1D("h_cut_gtr_ph", "HMS y' target;h.gtr.ph [rad];Events", 220, -0.08, 0.08);
        TH1D *h_cut_cer_npe = makeCut1D("h_cut_cer_npe", "HMS Cherenkov npe sum;NPE sum;Events", 220, 0.1, 20.0);
        TH1D *h_cut_cal_etotnorm = makeCut1D("h_cut_cal_etotnorm", "HMS calorimeter E/p;E_{tot}/p;Events", 220, 0.0, 2.0);
        TH1D *h_cut_cluster_e = makeCut1D("h_cut_cluster_e", "NPS cluster energy;E_{clus} [GeV];Clusters", 220, 0.1, 8.0);
        TH1D *h_cut_cluster_x = makeCut1D("h_cut_cluster_x", "NPS cluster x;x [cm];Clusters", 120, -34.0, 34.0);
        TH1D *h_cut_cluster_y = makeCut1D("h_cut_cluster_y", "NPS cluster y;y [cm];Clusters", 120, -40.0, 40.0);
        TH1D *h_cut_cluster_t = makeCut1D("h_cut_cluster_t", "NPS cluster time;t [ns];Clusters", 220, 120.0, 180.0);
        TH1D *h_cut_pair_dt = makeCut1D("h_cut_pair_dt", "Photon pair time difference;t_{1}-t_{2} [ns];Pairs", 220, -35.0, 35.0);
        TH1D *h_cut_mmiss_corr = makeCut1D("h_cut_mmiss_corr", "Corrected missing mass;M_{miss}^{corr} [GeV];Events", 220, 0.0, 2.0);

        TH2D *h_debug_react_xy = make2D("h_debug_react_xy", "HMS reaction vertex;x [cm];y [cm]", 180, -2.0, 2.0, 180, -2.0, 2.0);
        TH2D *h_debug_fast_raster_a_xy = make2D("h_debug_fast_raster_a_xy", "Fast raster A;fr_{xA} [cm];fr_{yA} [cm]", 240, -0.5, 0.5, 240, -0.5, 0.5);
        TH2D *h_debug_fast_raster_b_xy = make2D("h_debug_fast_raster_b_xy", "Fast raster B;fr_{xB} [cm];fr_{yB} [cm]", 240, -0.5, 0.5, 240, -0.5, 0.5);
        TH2D *h_debug_fp_xy = make2D("h_debug_fp_xy", "HMS focal-plane position;x_{fp} [cm];y_{fp} [cm]", 220, -60.0, 60.0, 180, -40.0, 40.0);
        TH2D *h_debug_fp_slopes = make2D("h_debug_fp_slopes", "HMS focal-plane slopes;x'_{fp} [rad];y'_{fp} [rad]", 180, -0.15, 0.15, 180, -0.08, 0.08);
        TH2D *h_debug_target_xy = make2D("h_debug_target_xy", "HMS target position;h.gtr.x [cm];h.gtr.y [cm]", 180, -5.0, 5.0, 180, -5.0, 5.0);
        TH2D *h_debug_target_slopes = make2D("h_debug_target_slopes", "HMS target slopes;h.gtr.th [rad];h.gtr.ph [rad]", 180, -0.16, 0.16, 180, -0.08, 0.08);
        TH2D *h_debug_cluster_xy = make2D("h_debug_cluster_xy", "NPS cluster position;x [cm];y [cm]", 120, -34.0, 34.0, 120, -40.0, 40.0);

        const double cluster_time_center_ns = nps_debug_cuts.time_center_ns;
        std::vector<CutDebugPlot> cut_debug_plots = {
            {h_cut_edtm_tdc, "HMS EDTM veto", false, 0.0, true, hms_debug_cuts.edtm_tdc_max},
            {h_cut_react_x, "HMS reaction x", false, 0.0, false, 0.0, has_react_x},
            {h_cut_react_y, "HMS reaction y", false, 0.0, false, 0.0, has_react_y},
            {h_cut_react_z, "HMS reaction z", true, hms_debug_cuts.react_z_min, true, hms_debug_cuts.react_z_max},
            {h_cut_fast_raster_xa, "Fast raster xA", false, 0.0, false, 0.0, has_fast_raster_xa},
            {h_cut_fast_raster_xb, "Fast raster xB", false, 0.0, false, 0.0, has_fast_raster_xb},
            {h_cut_fast_raster_ya, "Fast raster yA", false, 0.0, false, 0.0, has_fast_raster_ya},
            {h_cut_fast_raster_yb, "Fast raster yB", false, 0.0, false, 0.0, has_fast_raster_yb},
            {h_cut_fp_x, "HMS focal-plane x", false, 0.0, false, 0.0, has_fp_x},
            {h_cut_fp_y, "HMS focal-plane y", false, 0.0, false, 0.0, has_fp_y},
            {h_cut_fp_xp, "HMS focal-plane x'", false, 0.0, false, 0.0, has_fp_xp},
            {h_cut_fp_yp, "HMS focal-plane y'", false, 0.0, false, 0.0, has_fp_yp},
            {h_cut_gtr_x, "HMS target x", false, 0.0, false, 0.0, has_gtr_x},
            {h_cut_gtr_y, "HMS target y", false, 0.0, false, 0.0, has_gtr_y},
            {h_cut_delta, "HMS delta", true, hms_debug_cuts.delta_min, true, hms_debug_cuts.delta_max},
            {h_cut_gtr_th, "HMS x' target", true, hms_debug_cuts.gtr_th_min, true, hms_debug_cuts.gtr_th_max},
            {h_cut_gtr_ph, "HMS y' target", true, hms_debug_cuts.gtr_ph_min, true, hms_debug_cuts.gtr_ph_max},
            {h_cut_cer_npe, "HMS Cherenkov", true, hms_debug_cuts.cer_npe_sum_min, false, 0.0},
            {h_cut_cal_etotnorm, "HMS calorimeter", true, hms_debug_cuts.cal_etotnorm_min, true, hms_debug_cuts.cal_etotnorm_max},
            {h_cut_cluster_e, "NPS cluster energy", true, nps_debug_cuts.energy_min, false, 0.0},
            {h_cut_cluster_x, "NPS cluster x", true, nps_debug_cuts.x_min, true, nps_debug_cuts.x_max},
            {h_cut_cluster_y, "NPS cluster y", true, nps_debug_cuts.y_min, true, nps_debug_cuts.y_max},
            {h_cut_cluster_t, "NPS cluster time", true, cluster_time_center_ns - cluster_time_halfwidth_ns, true, cluster_time_center_ns + cluster_time_halfwidth_ns},
            {h_cut_pair_dt, "Photon pair timing", true, -pair_time_diff_max_ns, true, pair_time_diff_max_ns},
            {h_cut_mmiss_corr, "Weighted exclusivity", true, weighted_debug_cuts.mmiss_corr_min, true, weighted_debug_cuts.mmiss_corr_max}
        };

        std::vector<Debug2DPlot> cut_debug_2d_plots = {
            {h_debug_react_xy, "HMS reaction x vs y", has_react_x && has_react_y},
            {h_debug_fast_raster_a_xy, "Fast raster A: x vs y", has_fast_raster_xa && has_fast_raster_ya},
            {h_debug_fast_raster_b_xy, "Fast raster B: x vs y", has_fast_raster_xb && has_fast_raster_yb},
            {h_debug_fp_xy, "HMS focal-plane x vs y", has_fp_x && has_fp_y},
            {h_debug_fp_slopes, "HMS focal-plane x' vs y'", has_fp_xp && has_fp_yp},
            {h_debug_target_xy, "HMS target x vs y", has_gtr_x && has_gtr_y},
            {h_debug_target_slopes, "HMS target x' vs y'", true},
            {h_debug_cluster_xy, "NPS cluster x vs y", true}
        };

        auto write_cut_debug_pdf = [&](const TString& pdf_path,
                                       const std::vector<CutDebugPlot>& plots,
                                       const std::vector<Debug2DPlot>& plots_2d) {
        if (plots.empty()) return;
        TCanvas *c_cut_debug = new TCanvas(name("c_cut_debug"), Form("Cut Debug Run %d", run), 1600, 1000);
        c_cut_debug->Print(Form("%s[", pdf_path.Data()));
        TLatex latex;
        latex.SetNDC(true);
        latex.SetTextFont(42);
        latex.SetTextSize(0.045);
        for (std::size_t i = 0; i < plots.size(); i += 4) {
            c_cut_debug->Clear();
            c_cut_debug->Divide(2, 2);
            for (std::size_t j = 0; j < 4 && (i + j) < plots.size(); ++j) {
                const CutDebugPlot& plot = plots[i + j];
                c_cut_debug->cd(static_cast<int>(j + 1));
                gPad->SetLeftMargin(0.12);
                gPad->SetBottomMargin(0.12);
                gPad->SetTopMargin(0.10);
                gPad->SetRightMargin(0.05);
                gPad->SetTicks(1, 1);
                if (!plot.hist) continue;

                latex.SetTextSize(0.040);
                if (!plot.available) {
                    plot.hist->Draw("AXIS");
                    latex.DrawLatex(0.15, 0.93, Form("Run %d  %s", run, plot.label.Data()));
                    latex.SetTextColor(kRed + 1);
                    latex.DrawLatex(0.33, 0.50, "Branch unavailable");
                    latex.SetTextColor(kBlack);
                    continue;
                }

                // Crop display to populated structure while retaining configured cut lines.
                // Histogram contents, underflow/overflow, and event selection remain unchanged.
                const int nbins = plot.hist->GetNbinsX();
                const double peak = plot.hist->GetMaximum();
                const double visible_threshold = std::max(1.0, 1.0e-3 * peak);
                int first_bin = 1;
                int last_bin = nbins;
                while (first_bin <= nbins && plot.hist->GetBinContent(first_bin) < visible_threshold) ++first_bin;
                while (last_bin >= 1 && plot.hist->GetBinContent(last_bin) < visible_threshold) --last_bin;
                if (first_bin > last_bin) {
                    first_bin = 1;
                    last_bin = nbins;
                }
                auto include_cut_bin = [&](double x) {
                    if (x < plot.hist->GetXaxis()->GetXmin() || x > plot.hist->GetXaxis()->GetXmax()) return;
                    const int cut_bin = std::max(1, std::min(nbins, plot.hist->GetXaxis()->FindFixBin(x)));
                    first_bin = std::min(first_bin, cut_bin);
                    last_bin = std::max(last_bin, cut_bin);
                };
                if (plot.has_min) include_cut_bin(plot.min_value);
                if (plot.has_max) include_cut_bin(plot.max_value);
                const int pad_bins = std::max(2, (last_bin - first_bin + 1) / 20);
                first_bin = std::max(1, first_bin - pad_bins);
                last_bin = std::min(nbins, last_bin + pad_bins);
                plot.hist->GetXaxis()->SetRange(first_bin, last_bin);

                double visible_ymax = 0.0;
                for (int b = first_bin; b <= last_bin; ++b) {
                    visible_ymax = std::max(visible_ymax, plot.hist->GetBinContent(b));
                }
                plot.hist->SetMaximum(1.20 * std::max(1.0, visible_ymax));
                plot.hist->SetMinimum(0.0);
                plot.hist->Draw("HIST");

                const double yline = plot.hist->GetMaximum();
                auto draw_cut_line = [&](double x, int style) {
                    TLine* line = new TLine(x, 0.0, x, yline);
                    line->SetLineColor(kRed + 1);
                    line->SetLineStyle(style);
                    line->SetLineWidth(2);
                    line->Draw("same");
                };
                if (plot.has_min) draw_cut_line(plot.min_value, 2);
                if (plot.has_max) draw_cut_line(plot.max_value, 2);

                latex.DrawLatex(0.15, 0.93, Form("Run %d  %s", run, plot.label.Data()));
                latex.SetTextSize(0.030);
                TString cut_text;
                if (plot.has_min && plot.has_max) {
                    cut_text = Form("Cut: %.4g < x < %.4g", plot.min_value, plot.max_value);
                } else if (plot.has_min) {
                    cut_text = Form("Cut: x > %.4g", plot.min_value);
                } else if (plot.has_max) {
                    cut_text = Form("Cut: x < %.4g", plot.max_value);
                } else {
                    cut_text = "Cut: none";
                }
                latex.DrawLatex(0.15, 0.87, cut_text.Data());
                latex.DrawLatex(0.56, 0.87,
                                Form("Entries %.0f  UF %.0f  OF %.0f",
                                     plot.hist->GetEntries(),
                                     plot.hist->GetBinContent(0),
                                     plot.hist->GetBinContent(nbins + 1)));
            }
            c_cut_debug->Print(pdf_path.Data());
        }

        for (std::size_t i = 0; i < plots_2d.size(); i += 4) {
            c_cut_debug->Clear();
            c_cut_debug->Divide(2, 2);
            for (std::size_t j = 0; j < 4 && (i + j) < plots_2d.size(); ++j) {
                const Debug2DPlot& plot = plots_2d[i + j];
                c_cut_debug->cd(static_cast<int>(j + 1));
                gPad->SetLeftMargin(0.12);
                gPad->SetBottomMargin(0.12);
                gPad->SetTopMargin(0.10);
                gPad->SetRightMargin(0.15);
                gPad->SetTicks(1, 1);
                if (!plot.hist) continue;

                latex.SetTextSize(0.040);
                if (!plot.available) {
                    plot.hist->Draw("AXIS");
                    latex.DrawLatex(0.15, 0.93, Form("Run %d  %s", run, plot.label.Data()));
                    latex.SetTextColor(kRed + 1);
                    latex.DrawLatex(0.33, 0.50, "Branch unavailable");
                    latex.SetTextColor(kBlack);
                    continue;
                }

                const int nbins_x = plot.hist->GetNbinsX();
                const int nbins_y = plot.hist->GetNbinsY();
                const double threshold = std::max(1.0, 1.0e-3 * plot.hist->GetMaximum());
                int first_x = nbins_x;
                int last_x = 1;
                int first_y = nbins_y;
                int last_y = 1;
                bool found = false;
                for (int bx = 1; bx <= nbins_x; ++bx) {
                    for (int by = 1; by <= nbins_y; ++by) {
                        if (plot.hist->GetBinContent(bx, by) < threshold) continue;
                        found = true;
                        first_x = std::min(first_x, bx);
                        last_x = std::max(last_x, bx);
                        first_y = std::min(first_y, by);
                        last_y = std::max(last_y, by);
                    }
                }
                if (found) {
                    const int pad_x = std::max(2, (last_x - first_x + 1) / 20);
                    const int pad_y = std::max(2, (last_y - first_y + 1) / 20);
                    plot.hist->GetXaxis()->SetRange(std::max(1, first_x - pad_x), std::min(nbins_x, last_x + pad_x));
                    plot.hist->GetYaxis()->SetRange(std::max(1, first_y - pad_y), std::min(nbins_y, last_y + pad_y));
                }
                plot.hist->Draw("COLZ");
                latex.DrawLatex(0.15, 0.93, Form("Run %d  %s", run, plot.label.Data()));
                latex.SetTextSize(0.030);
                latex.DrawLatex(0.70, 0.93, Form("Entries %.0f", plot.hist->GetEntries()));
            }
            c_cut_debug->Print(pdf_path.Data());
        }

        c_cut_debug->Clear();
        c_cut_debug->Divide(1, 1);
        c_cut_debug->cd(1);
        gPad->SetLeftMargin(0.10);
        gPad->SetBottomMargin(0.10);
        gPad->SetTopMargin(0.12);
        gPad->SetRightMargin(0.14);
        gPad->SetTicks(1, 1);

        TH2D* h_dead_grid = new TH2D(name("h_dead_block_grid"),
                                     "NPS dead-block 3x3 mask;Column;Row",
                                     nps_dead_block::Mask::kNCols, -0.5, nps_dead_block::Mask::kNCols - 0.5,
                                     nps_dead_block::Mask::kNRows, -0.5, nps_dead_block::Mask::kNRows - 0.5);
        h_dead_grid->SetDirectory(nullptr);
        h_dead_grid->SetStats(0);
        h_dead_grid->GetZaxis()->SetTitle("0=active, 1=excluded, 2=dead");
        h_dead_grid->SetMinimum(0.0);
        h_dead_grid->SetMaximum(2.0);

        for (int block : dead_block_mask.excluded_blocks) {
            int col = -1;
            int row = -1;
            if (nps_dead_block::Mask::block_to_col_row(block, col, row)) {
                h_dead_grid->SetBinContent(col + 1, row + 1, 1.0);
            }
        }
        for (int block : dead_block_mask.dead_blocks) {
            int col = -1;
            int row = -1;
            if (nps_dead_block::Mask::block_to_col_row(block, col, row)) {
                h_dead_grid->SetBinContent(col + 1, row + 1, 2.0);
            }
        }

        h_dead_grid->Draw("COLZ");

        TLatex grid_latex;
        grid_latex.SetNDC(true);
        grid_latex.SetTextFont(42);
        grid_latex.SetTextSize(0.032);
        grid_latex.DrawLatex(0.12, 0.94, Form("Run %d dead blocks: %s", run, nps_dead_block::join_blocks(dead_block_mask.dead_blocks).c_str()));
        if (dead_block_mask.empty()) {
            grid_latex.DrawLatex(0.12, 0.90, "No dead/missing block IDs configured for this run; all 1080 blocks shown active");
        } else {
            grid_latex.DrawLatex(0.12, 0.90, "Dead block and full 3x3 excluded region shown on 1080-block NPS grid");
        }

        grid_latex.SetTextSize(0.020);
        grid_latex.SetTextAlign(22);
        for (int block : dead_block_mask.dead_blocks) {
            int col = -1;
            int row = -1;
            if (nps_dead_block::Mask::block_to_col_row(block, col, row)) {
                TLatex label;
                label.SetTextFont(42);
                label.SetTextSize(0.018);
                label.SetTextAlign(22);
                label.DrawLatex(col, row, Form("%d", block));
            }
        }

        c_cut_debug->Print(pdf_path.Data());
        safe_delete(h_dead_grid);
        c_cut_debug->Print(Form("%s]", pdf_path.Data()));
        safe_delete(c_cut_debug);
        logmsg(INFO, Form("Wrote cut-debug PDF to %s", pdf_path.Data()));
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
        vector<pair<double,double>> diag_windows;
        vector<pair<double,double>> side_windows;
        auto coin_win = nps::default_coin_window();
        pair<double,double> full1_t1;
        pair<double,double> full1_t2;
        pair<double,double> full2_t1;
        pair<double,double> full2_t2;

        if (use_shifted_timing_windows) {
            diag_windows = {
                {139.0, 141.0}, {141.0, 143.0}, {143.0, 145.0},
                {155.0, 157.0}, {157.0, 159.0}, {159.0, 161.0}
            };
            side_windows = diag_windows;
            full1_t1 = std::make_pair(155.0, 161.0);
            full1_t2 = std::make_pair(139.0, 145.0);
            full2_t1 = std::make_pair(139.0, 145.0);
            full2_t2 = std::make_pair(155.0, 161.0);
        } else {
            diag_windows = nps::default_diag_windows();
            side_windows = nps::default_side_windows();
            full1_t1 = nps::default_full_acc1_t1();
            full1_t2 = nps::default_full_acc1_t2();
            full2_t1 = nps::default_full_acc2_t1();
            full2_t2 = nps::default_full_acc2_t2();
        }

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

                // Apply good-event cut as early as possible to avoid unnecessary processing.
                if (apply_good_event_cut_for_run && !passes_good_event_cut(g_evnum, good_ranges_for_run)) {
                    continue;
                }

                fill_cut_hist(h_cut_edtm_tdc, edtmtdc);
                if (has_react_x) fill_cut_hist(h_cut_react_x, hreactx);
                if (has_react_y) fill_cut_hist(h_cut_react_y, hreacty);
                fill_cut_hist(h_cut_react_z, hreactz);
                if (has_fast_raster_xa) fill_cut_hist(h_cut_fast_raster_xa, hfastRasterXA);
                if (has_fast_raster_xb) fill_cut_hist(h_cut_fast_raster_xb, hfastRasterXB);
                if (has_fast_raster_ya) fill_cut_hist(h_cut_fast_raster_ya, hfastRasterYA);
                if (has_fast_raster_yb) fill_cut_hist(h_cut_fast_raster_yb, hfastRasterYB);
                if (has_fp_x) fill_cut_hist(h_cut_fp_x, hxfp);
                if (has_fp_y) fill_cut_hist(h_cut_fp_y, hyfp);
                if (has_fp_xp) fill_cut_hist(h_cut_fp_xp, hxpfp);
                if (has_fp_yp) fill_cut_hist(h_cut_fp_yp, hypfp);
                if (has_gtr_x) fill_cut_hist(h_cut_gtr_x, HgtrX);
                if (has_gtr_y) fill_cut_hist(h_cut_gtr_y, HgtrY);
                fill_cut_hist(h_cut_delta, hdelta);
                fill_cut_hist(h_cut_gtr_th, HgtrTh);
                fill_cut_hist(h_cut_gtr_ph, HgtrPh);
                fill_cut_hist(h_cut_cer_npe, hcernpeSum);
                fill_cut_hist(h_cut_cal_etotnorm, hcaletotnorm);
                if (has_react_x && has_react_y) fill_debug_hist_2d(h_debug_react_xy, hreactx, hreacty);
                if (has_fast_raster_xa && has_fast_raster_ya) fill_debug_hist_2d(h_debug_fast_raster_a_xy, hfastRasterXA, hfastRasterYA);
                if (has_fast_raster_xb && has_fast_raster_yb) fill_debug_hist_2d(h_debug_fast_raster_b_xy, hfastRasterXB, hfastRasterYB);
                if (has_fp_x && has_fp_y) fill_debug_hist_2d(h_debug_fp_xy, hxfp, hyfp);
                if (has_fp_xp && has_fp_yp) fill_debug_hist_2d(h_debug_fp_slopes, hxpfp, hypfp);
                if (has_gtr_x && has_gtr_y) fill_debug_hist_2d(h_debug_target_xy, HgtrX, HgtrY);
                fill_debug_hist_2d(h_debug_target_slopes, HgtrTh, HgtrPh);

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
                        clusT[i] = (*clusT_vec)[i] + branch_map.cluster_time_offset_ns;
                    }
                } else {
                    // fallback: use nclust_dbl + arrays (if the tree stored arrays)
                    nclust = static_cast<int>(std::lround(nclust_dbl));
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
                if (!acceptance_cuts->pass_hms_data(edtmtdc, hdelta, HgtrTh, HgtrPh, hcernpeSum, hcaletotnorm, hreactz)) continue;
                ++n_pass_hms;

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

                for (int i=0; i<n_after; ++i) {
                    fill_cut_hist(h_cut_cluster_e, clusE[i]);
                    fill_cut_hist(h_cut_cluster_x, clusX[i]);
                    fill_cut_hist(h_cut_cluster_y, clusY[i]);
                    fill_debug_hist_2d(h_debug_cluster_xy, clusX[i], clusY[i]);
                    fill_cut_hist(h_cut_cluster_t, clusT[i]);
                }

                // select good clusters by spatial/energy/timing cuts
                vector<int> good_idx;
                collect_good_clusters(*acceptance_cuts,
                                      clusE,
                                      clusX,
                                      clusY,
                                      clusT,
                                      n_after,
                                      cluster_time_halfwidth_ns,
                                      good_idx,
                                      &dead_block_mask);

                // DEBUG: Show why clusters are failing cuts
                static int debug_cuts = 0;
                static int debug_timing = 0;
                if (debug_timing < 10 && n_after >= 2) {
                    cout << "\n[DEBUG TIMING] Event " << ev << ": n_after=" << n_after << "\n";
                    for (int i=0; i<n_after; ++i) {
                        bool passes = acceptance_cuts->pass_nps_cluster(clusE[i], clusX[i], clusY[i], clusT[i], cluster_time_halfwidth_ns);
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

                if (good_idx.size() == 2) ++n_mult2_hms_nps;
                else if (good_idx.size() == 3) ++n_mult3_hms_nps;
                else if (good_idx.size() == 4) ++n_mult4_hms_nps;

                const double px_e = HgtrPx, py_e = HgtrPy, pz_e = HgtrPz;
                const double p_e_mom = sqrt(max(0.0, px_e*px_e + py_e*py_e + pz_e*pz_e));
                const double Ee = sqrt(max(0.0, p_e_mom*p_e_mom + nps::kElectronMass_GeV*nps::kElectronMass_GeV));

                // select pair for pi0
                int sel_i = -1;
                int sel_j = -1;
                if (!choose_pi0_pair(good_idx,
                                     clusE,
                                     clusX,
                                     clusY,
                                     clusT,
                                     run_z_nps_cm,
                                     pair_time_diff_max_ns,
                                     sel_i,
                                     sel_j)) {
                    continue;
                }

                // fill timing vs timing
                const double t1 = clusT[sel_i], t2 = clusT[sel_j];
                fill_cut_hist(h_cut_pair_dt, t1 - t2);
                h_t1_t2->Fill(t2, t1);
                h_t1_proj->Fill(t1);
                h_t2_proj->Fill(t2);

                // invariant mass
                const double mgg = nps::invariant_mass_pi0(clusE[sel_i], clusE[sel_j], clusX[sel_i], clusX[sel_j], clusY[sel_i], clusY[sel_j], run_z_nps_cm);
                h_mpi0_all->Fill(mgg);
                if (good_idx.size()==2) h_mpi0_2->Fill(mgg);
                else if (good_idx.size()==3) h_mpi0_3->Fill(mgg);
                else if (good_idx.size()==4) h_mpi0_4->Fill(mgg);

                // opening angle & energy ratio
                {
                    const double z_nps = run_z_nps_cm;
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

                const double mm_p_no_offset = nps::missing_mass_proton_pi0(run_beam_energy, Ee, px_e, py_e, pz_e,
                                                                          clusE[sel_i], clusE[sel_j],
                                                                          clusX[sel_i], clusY[sel_i],
                                                                          clusX[sel_j], clusY[sel_j],
                                                                          run_z_nps_cm, -run_nps_theta_deg);

                const double mm_p = nps::missing_mass_proton_pi0(run_beam_energy, Ee_scaled, px_e_scaled, py_e_scaled, pz_e_scaled,
                                                                clusE[sel_i], clusE[sel_j],
                                                                clusX[sel_i], clusY[sel_i],
                                                                clusX[sel_j], clusY[sel_j],
                                                                run_z_nps_cm, -run_nps_theta_deg);

                const double mm_p_corr = nps::invariant_missing_mass_corrected_avnish_from_detector(mm_p*mm_p, run_beam_energy, Ee_scaled, px_e_scaled, py_e_scaled, pz_e_scaled,
                                                                clusE[sel_i], clusE[sel_j],
                                                                clusX[sel_i], clusY[sel_i],
                                                                clusX[sel_j], clusY[sel_j],
                                                                run_z_nps_cm, -run_nps_theta_deg);

                const double mm_p_corr_no_mom_offset = nps::invariant_missing_mass_corrected_avnish_from_detector(mm_p_no_offset*mm_p_no_offset, run_beam_energy, Ee, px_e, py_e, pz_e,
                                                                clusE[sel_i], clusE[sel_j],
                                                                clusX[sel_i], clusY[sel_i],
                                                                clusX[sel_j], clusY[sel_j],
                                                                run_z_nps_cm, -run_nps_theta_deg);

                fill_cut_hist(h_cut_mmiss_corr, mm_p_corr);

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
                    run_beam_energy, Ee, px_e, py_e, pz_e,
                    clusE[sel_i], clusX[sel_i], clusY[sel_i],
                    clusE[sel_j], clusX[sel_j], clusY[sel_j],
                    run_z_nps_cm, -run_nps_theta_deg, false);

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
                entry.is_weighted = 0;   // Will be set in second pass
                const int event_helicity = static_cast<int>(std::llround(t_helicity_value));
                entry.helicity = (has_t_helicity_branch && event_helicity != 0)
                    ? (event_helicity > 0 ? 1 : -1)
                    : helicity_for_event(g_evnum, helicity_intervals);
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

                    if (apply_good_event_cut_for_run && !passes_good_event_cut(g_evnum, good_ranges_for_run)) {
                        continue;
                    }

                    // Load cluster data using the same mapping and timing treatment as first pass.
                    const int nclust = load_clusters_from_branches(
                        clusE_vec,
                        clusX_vec,
                        clusY_vec,
                        clusT_vec,
                        nclust_dbl,
                        clusE,
                        clusX,
                        clusY,
                        clusT,
                        branch_map.cluster_time_offset_ns
                    );

                    if (nclust < 2) continue;
                    if (!acceptance_cuts->pass_hms_data(edtmtdc, hdelta, HgtrTh, HgtrPh, hcernpeSum, hcaletotnorm, hreactz)) continue;

                    const int n_after = nps::packClusters(clusE, clusX, clusY, clusT, nclust);
                    vector<int> good_idx;
                    collect_good_clusters(*acceptance_cuts,
                                          clusE,
                                          clusX,
                                          clusY,
                                          clusT,
                                          n_after,
                                          cluster_time_halfwidth_ns,
                                          good_idx,
                                          &dead_block_mask);

                    if (good_idx.size() < 2) continue;
                    int sel_i = -1;
                    int sel_j = -1;
                    if (!choose_pi0_pair(good_idx,
                                         clusE,
                                         clusX,
                                         clusY,
                                         clusT,
                                         run_z_nps_cm,
                                         pair_time_diff_max_ns,
                                         sel_i,
                                         sel_j)) {
                        continue;
                    }

                    const double px_e = HgtrPx, py_e = HgtrPy, pz_e = HgtrPz;
                    const double p_e_mom = sqrt(max(0.0, px_e*px_e + py_e*py_e + pz_e*pz_e));
                    const double Ee = sqrt(max(0.0, p_e_mom*p_e_mom + nps::kElectronMass_GeV*nps::kElectronMass_GeV));

                    const double mgg = nps::invariant_mass_pi0(clusE[sel_i], clusE[sel_j], clusX[sel_i], clusX[sel_j], clusY[sel_i], clusY[sel_j], run_z_nps_cm);

                    const double px_e_scaled = px_e * HMS_MOM_OFFSET_SCALE;
                    const double py_e_scaled = py_e * HMS_MOM_OFFSET_SCALE;
                    const double pz_e_scaled = pz_e * HMS_MOM_OFFSET_SCALE;
                    const double p_e_mom_scaled = sqrt(max(0.0, px_e_scaled*px_e_scaled + py_e_scaled*py_e_scaled + pz_e_scaled*pz_e_scaled));
                    const double Ee_scaled = sqrt(max(0.0, p_e_mom_scaled*p_e_mom_scaled + nps::kElectronMass_GeV*nps::kElectronMass_GeV));

                    const double mm_p_no_offset = nps::missing_mass_proton_pi0(run_beam_energy, Ee, px_e, py_e, pz_e, clusE[sel_i], clusE[sel_j], clusX[sel_i], clusY[sel_i], clusX[sel_j], clusY[sel_j], run_z_nps_cm, -run_nps_theta_deg);
                    (void)mm_p_no_offset;

                    const double mm_p = nps::missing_mass_proton_pi0(run_beam_energy, Ee_scaled, px_e_scaled, py_e_scaled, pz_e_scaled, clusE[sel_i], clusE[sel_j], clusX[sel_i], clusY[sel_i], clusX[sel_j], clusY[sel_j], run_z_nps_cm, -run_nps_theta_deg);

                    const double mm_p_corr = nps::invariant_missing_mass_corrected_avnish_from_detector(mm_p*mm_p, run_beam_energy, Ee_scaled, px_e_scaled, py_e_scaled, pz_e_scaled, clusE[sel_i], clusE[sel_j], clusX[sel_i], clusY[sel_i], clusX[sel_j], clusY[sel_j], run_z_nps_cm, -run_nps_theta_deg);

                    bool isExclusive = acceptance_cuts->pass_weighted_exclusive(mm_p_corr);
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
                        auto phys = nps::compute_physics_vars_from_detector(run_beam_energy, Ee, px_e, py_e, pz_e, clusE[sel_i], clusX[sel_i], clusY[sel_i], clusE[sel_j], clusX[sel_j], clusY[sel_j], run_z_nps_cm, -run_nps_theta_deg, false);
                        
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
                        double nu = run_beam_energy - Ee;
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
        auto draw_overlay_triplet = [&](const char* canvas_name,
                        const char* canvas_title,
                        TH1D* h_all,
                        TH1D* h_weighted,
                        TH1D* h_excl_weighted,
                        const char* output_stem,
                        double lx1 = 0.6,
                        double ly1 = 0.6,
                        double lx2 = 0.88,
                        double ly2 = 0.88) -> TCanvas* {
            TCanvas *c = new TCanvas(canvas_name, canvas_title, 800, 600);
            style_canvas(c);

            h_all->SetLineColor(kBlack);
            h_all->SetLineWidth(2);
            h_weighted->SetLineColor(kBlue);
            h_weighted->SetLineWidth(2);
            h_excl_weighted->SetLineColor(kRed);
            h_excl_weighted->SetLineWidth(2);

            set_overlay_y_range({h_all, h_weighted, h_excl_weighted});
            h_all->Draw("HIST");
            h_weighted->Draw("HIST SAME");
            h_excl_weighted->Draw("HIST SAME");

            TLegend *leg = new TLegend(lx1, ly1, lx2, ly2);
            leg->SetBorderSize(0);
            leg->SetFillColor(0);
            leg->AddEntry(h_all, "All", "l");
            leg->AddEntry(h_weighted, "Weighted", "l");
            leg->AddEntry(h_excl_weighted, "Weighted+Exclusive", "l");
            leg->Draw();

            c->SaveAs(Form("%s/%s_overlay_run%d.png", outPlotDir.Data(), output_stem, run));
            return c;
        };

        TString cut_debug_pdf = Form("%s/cut_debug_run%d.pdf", outPlotDir.Data(), run);
        write_cut_debug_pdf(cut_debug_pdf, cut_debug_plots, cut_debug_2d_plots);

        TCanvas *c_Q2_overlay = draw_overlay_triplet(Form("c_Q2_overlay_run%d", run), "Q^{2} Overlay",
                                                     h_Q2, h_Q2_weighted, h_Q2_excl_weighted, "Q2");
        TCanvas *c_W_overlay = draw_overlay_triplet(Form("c_W_overlay_run%d", run), "W Overlay",
                                                    h_W, h_W_weighted, h_W_excl_weighted, "W");
        TCanvas *c_t_overlay = draw_overlay_triplet(Form("c_t_overlay_run%d", run), "t Overlay",
                                                    h_t, h_t_weighted, h_t_excl_weighted, "t",
                                                    0.15, 0.6, 0.43, 0.88);
        TCanvas *c_tmin_overlay = draw_overlay_triplet(Form("c_tmin_overlay_run%d", run), "t_min Overlay",
                                                       h_tmin, h_tmin_weighted, h_tmin_excl_weighted, "tmin");
        TCanvas *c_pt_overlay = draw_overlay_triplet(Form("c_pt_overlay_run%d", run), "p_T Overlay",
                                                     h_pt, h_pt_weighted, h_pt_excl_weighted, "pt");
        TCanvas *c_s_overlay = draw_overlay_triplet(Form("c_s_overlay_run%d", run), "s Overlay",
                                                    h_s, h_s_weighted, h_s_excl_weighted, "s");
        TCanvas *c_xB_overlay = draw_overlay_triplet(Form("c_xB_overlay_run%d", run), "x_B Overlay",
                                                     h_xB, h_xB_weighted, h_xB_excl_weighted, "xB");
        TCanvas *c_z_overlay = draw_overlay_triplet(Form("c_z_overlay_run%d", run), "z Overlay",
                                                    h_z, h_z_weighted, h_z_excl_weighted, "z");
        TCanvas *c_theta_overlay = draw_overlay_triplet(Form("c_theta_overlay_run%d", run), "#theta Overlay",
                                                        h_theta, h_theta_weighted, h_theta_excl_weighted, "theta");
        TCanvas *c_phi_overlay = draw_overlay_triplet(Form("c_phi_overlay_run%d", run), "#phi Overlay",
                                                      h_phi, h_phi_weighted, h_phi_excl_weighted, "phi");

        // Missing mass overlay with weighted (NO FITTING - just overlay display)
        // IMPORTANT: Create this BEFORE fitting h_mmiss_all_weighted so fit curves don't appear
        TCanvas *c_mmiss_overlay = new TCanvas(Form("c_mmiss_overlay_run%d", run), "Missing Mass Overlay", 800, 600);
        style_canvas(c_mmiss_overlay);
        c_mmiss_overlay->SetName(Form("c_mmiss_overlay_run%d", run));
        c_mmiss_overlay->SetTitle(Form("Missing Mass Overlay run %d", run));
        h_mmiss_all->SetLineColor(kBlack); h_mmiss_all->SetLineWidth(2);
        h_mmiss_all_corr->SetLineColor(kBlue); h_mmiss_all_corr->SetLineWidth(2);
        h_mmiss_all_weighted->SetLineColor(kGreen+2); h_mmiss_all_weighted->SetLineWidth(2);
        set_overlay_y_range({h_mmiss_all, h_mmiss_all_corr, h_mmiss_all_weighted}, 0.4, 2.0);
        h_mmiss_all->GetXaxis()->SetRangeUser(0.4, 2.0);
        h_mmiss_all_corr->GetXaxis()->SetRangeUser(0.4, 2.0);
        h_mmiss_all_weighted->GetXaxis()->SetRangeUser(0.4, 2.0);
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
        style_canvas(c_mmiss_fit);
        c_mmiss_fit->SetName(Form("c_mmiss_fit_run%d", run));
        c_mmiss_fit->SetTitle(Form("Missing Mass Weighted - Fits run %d", run));
        
        h_mmiss_all_weighted->SetLineColor(kBlack);
        h_mmiss_all_weighted->SetLineWidth(2);
        h_mmiss_all_weighted->GetXaxis()->SetLabelSize(0.045);
        h_mmiss_all_weighted->GetYaxis()->SetLabelSize(0.045);
        h_mmiss_all_weighted->GetXaxis()->SetTitleSize(0.050);
        h_mmiss_all_weighted->GetYaxis()->SetTitleSize(0.050);
        h_mmiss_all_weighted->GetXaxis()->SetTitle("Missing Mass (GeV)");
        h_mmiss_all_weighted->GetYaxis()->SetTitle("Counts");
        {
            // Display scaling should follow the proton-peak region, not distant tails/outliers.
            const double draw_x_lo = 0.4;
            const double draw_x_hi = 2.0;
            const double fit_eval_lo = 0.5;
            const double fit_eval_hi = 1.8;

            const int b_lo = std::max(1, h_mmiss_all_weighted->FindBin(draw_x_lo));
            const int b_hi = std::min(h_mmiss_all_weighted->GetNbinsX(), h_mmiss_all_weighted->FindBin(draw_x_hi));

            double data_peak_max = 0.0;
            double data_peak_min = std::numeric_limits<double>::infinity();
            std::vector<double> peak_contents;
            peak_contents.reserve(std::max(1, b_hi - b_lo + 1));
            for (int b = b_lo; b <= b_hi; ++b) {
                const double y = h_mmiss_all_weighted->GetBinContent(b);
                if (!std::isfinite(y)) continue;
                data_peak_max = std::max(data_peak_max, y);
                data_peak_min = std::min(data_peak_min, y);
                peak_contents.push_back(y);
            }

            if (!std::isfinite(data_peak_max) || data_peak_max <= 0.0) {
                data_peak_max = h_mmiss_all_weighted->GetMaximum();
            }
            if (!std::isfinite(data_peak_max) || data_peak_max <= 0.0) {
                data_peak_max = 1.0;
            }
            if (!std::isfinite(data_peak_min)) {
                data_peak_min = 0.0;
            }

            double robust_peak_max = data_peak_max;
            if (!peak_contents.empty()) {
                std::sort(peak_contents.begin(), peak_contents.end());
                const size_t q_idx = static_cast<size_t>(0.98 * static_cast<double>(peak_contents.size() - 1));
                robust_peak_max = peak_contents[q_idx];
                if (!std::isfinite(robust_peak_max) || robust_peak_max <= 0.0) {
                    robust_peak_max = data_peak_max;
                }
            }

            // If a single outlier bin dominates, scale to the robust peak instead.
            double ref_peak_max = data_peak_max;
            if (robust_peak_max > 0.0 && data_peak_max > 2.5 * robust_peak_max) {
                ref_peak_max = robust_peak_max;
            }
            if (!std::isfinite(ref_peak_max) || ref_peak_max <= 0.0) {
                ref_peak_max = 1.0;
            }

            double fit_peak_max = std::max(fit_combined->GetMaximum(fit_eval_lo, fit_eval_hi),
                                           std::max(fit1->GetMaximum(fit_eval_lo, fit_eval_hi),
                                                    fit2->GetMaximum(fit_eval_lo, fit_eval_hi)));
            if (!std::isfinite(fit_peak_max) || fit_peak_max <= 0.0) {
                fit_peak_max = ref_peak_max;
            }

            // Keep bad fits from exploding the axis and hiding the physical peak.
            fit_peak_max = std::min(fit_peak_max, 1.35 * ref_peak_max);

            const double mmiss_fit_ymax = std::max(ref_peak_max, fit_peak_max);
            const double mmiss_fit_ymin = (data_peak_min < 0.0) ? (1.15 * data_peak_min) : 0.0;

            h_mmiss_all_weighted->GetXaxis()->SetRangeUser(draw_x_lo, draw_x_hi);
            h_mmiss_all_weighted->SetMinimum(mmiss_fit_ymin);
            h_mmiss_all_weighted->SetMaximum(1.20 * mmiss_fit_ymax);
        }
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

        // Restore full x-range for subsequent overlay redraw/output objects.
        h_mmiss_all_weighted->GetXaxis()->UnZoom();

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
        style_canvas(c_mmiss_overlay);
        set_overlay_y_range({h_mmiss_all, h_mmiss_all_corr, h_mmiss_all_weighted}, 0.4, 2.0);
        h_mmiss_all->GetXaxis()->SetRangeUser(0.4, 2.0);
        h_mmiss_all_corr->GetXaxis()->SetRangeUser(0.4, 2.0);
        h_mmiss_all_weighted->GetXaxis()->SetRangeUser(0.4, 2.0);
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
        TString outf = Form("%s/diagnostics_run%d.root", outRootDir.Data(), run);
        TFile *fout = TFile::Open(outf, "RECREATE");
        if (!fout || fout->IsZombie()) {
        logmsg(WARN, Form("Could not create ROOT output %s", outf.Data()));
        } else {
        fout->cd();

        // Create, fill, and write event-level physics tree from treeData map.
        write_event_physics_tree(fout, treeData);
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

        for (const auto& plot : cut_debug_plots) {
            if (plot.hist) plot.hist->Write();
        }
        for (const auto& plot : cut_debug_2d_plots) {
            if (plot.hist) plot.hist->Write();
        }

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
        style_canvas(c_t12, true);
        c_t12->SetName(name("c_t12"));
        h_t1_t2->Draw("COLZ");
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
        style_canvas(c_pi0);
        c_pi0->SetName(name("c_pi0"));
        h_mpi0_all->SetLineColor(kBlack); h_mpi0_all->SetLineWidth(1);
        h_m_pi0_coin->SetLineColor(kBlue); h_m_pi0_coin->SetLineWidth(2);
        {
            vector<TH1*> pi0_overlay_hists = {h_mpi0_all, h_m_pi0_coin};
            if (h_coin_bgsub) pi0_overlay_hists.push_back(h_coin_bgsub);
            if (h_final) pi0_overlay_hists.push_back(h_final);
            set_overlay_y_range(pi0_overlay_hists);
        }
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
        style_canvas(c_mmiss_mgg, true);
        c_mmiss_mgg->SetName(name("c_mmiss_mgg"));
        h_mmiss_vs_mgg->Draw("COLZ");
        c_mmiss_mgg->SaveAs(Form("%s/mmiss_vs_mgg_run%d.png", outPlotDir.Data(), run));

        // 4) cluster diagnostics
        TCanvas *c_cluster = new TCanvas(name("c_cluster"), "Cluster diagnostics", 1000,800);
        c_cluster->SetName(name("c_cluster"));
        c_cluster->Divide(2,2);
        c_cluster->cd(1);
        gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.08); gPad->SetRightMargin(0.05); gPad->SetTicks(1,1);
        set_overlay_y_range({h_clustE});
        h_clustE->Draw();
        c_cluster->cd(2);
        gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.08); gPad->SetRightMargin(0.05); gPad->SetTicks(1,1);
        set_overlay_y_range({h_clustT});
        h_clustT->Draw();
        c_cluster->cd(3);
        gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.08); gPad->SetRightMargin(0.14); gPad->SetTicks(1,1);
        h_clustXY->Draw("COLZ");
        c_cluster->cd(4);
        gPad->SetLeftMargin(0.12); gPad->SetBottomMargin(0.12); gPad->SetTopMargin(0.08); gPad->SetRightMargin(0.05); gPad->SetTicks(1,1);
        set_overlay_y_range({h_clustE_sum});
        h_clustE_sum->Draw();
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
        // Optional peak placeholders retained in summary schema
        // -------------------------
        double s1x_peak = 0.0, s1x_err = 0.0;
        double s1y_peak = 0.0, s1y_err = 0.0;
        double s2x_peak = 0.0, s2x_err = 0.0;
        double s2y_peak = 0.0, s2y_err = 0.0;

        // -------------------------
        // CRITICAL: Write CSV data to per-run temp file BEFORE summary print
        // Bash will reconstruct global CSV from successful per-run files
        // -------------------------
        {
        RunSummaryRow summary_row;
        summary_row.run = run;
        summary_row.accumulated_charge_mC = accumulated_charge_mC;
        summary_row.hel_pos_charge_mC = hel_pos_charge_mC;
        summary_row.hel_neg_charge_mC = hel_neg_charge_mC;
        summary_row.current_mean_uA = run_current_mean;
        summary_row.beam_time_s = beam_time;
        summary_row.total_entries = n_total;
        summary_row.pass_hms = n_pass_hms;
        summary_row.pass_hms_nps = n_selected_for_analysis;
        summary_row.total_coin_entries = bg.n_coin_raw;
        summary_row.estimated_time_accidentals = bg.n_accidentals;
        summary_row.chi2_ndf_comb_bg = res.chi2_ndf;
        summary_row.pi0_mu_MeV = res.mu_MeV;
        summary_row.pi0_sigma_MeV = res.sigma_MeV;
        summary_row.pi0_signal_counts = res.signal_counts;
        summary_row.mmiss_p_mean_GeV = mmiss_p_mean;
        summary_row.mmiss_p_sigma_GeV = mmiss_p_sigma;
        summary_row.s1x_peak = s1x_peak;
        summary_row.s1x_err = s1x_err;
        summary_row.s1y_peak = s1y_peak;
        summary_row.s1y_err = s1y_err;
        summary_row.s2x_peak = s2x_peak;
        summary_row.s2x_err = s2x_err;
        summary_row.s2y_peak = s2y_peak;
        summary_row.s2y_err = s2y_err;
        summary_row.run_status = run_status;

        const std::string kin_safe = sanitize_token(kinematic.Data());
        TString per_run_csv = Form("/tmp/nps_csv_%s_run_%d_pid_%d_%lld.txt",
                       kin_safe.c_str(),
                       run,
                       gSystem->GetPid(),
                       static_cast<Long64_t>(time(nullptr)));
        ofstream fg(per_run_csv.Data(), ios::out);
        if (fg.is_open()) {
            write_summary_csv_row(fg, summary_row);
            fg.flush();
            fg.close();
            cout << "[CSV_WRITTEN] " << per_run_csv << "\n";
        } else {
            logmsg(WARN, Form("Could not write per-run CSV %s", per_run_csv.Data()));
        }

        if (write_global_csv) {
            ofstream fglobal(global_csv.Data(), ios::app);
            if (fglobal.is_open()) {
                write_summary_csv_row(fglobal, summary_row);
                fglobal.flush();
                fglobal.close();
            } else {
                logmsg(WARN, Form("Could not append global CSV %s (errno=%d, %s)",
                                  global_csv.Data(), errno, std::strerror(errno)));
            }
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
        RunSummaryRow txt_row;
        txt_row.run = run;
        txt_row.accumulated_charge_mC = accumulated_charge_mC;
        txt_row.hel_pos_charge_mC = hel_pos_charge_mC;
        txt_row.hel_neg_charge_mC = hel_neg_charge_mC;
        txt_row.total_entries = n_total;
        txt_row.pass_hms = n_pass_hms;
        txt_row.pass_hms_nps = n_selected_for_analysis;
        txt_row.total_coin_entries = bg.n_coin_raw;
        txt_row.estimated_time_accidentals = bg.n_accidentals;
        txt_row.chi2_ndf_comb_bg = res.chi2_ndf;
        txt_row.pi0_mu_MeV = res.mu_MeV;
        txt_row.pi0_sigma_MeV = res.sigma_MeV;
        txt_row.pi0_signal_counts = res.signal_counts;
        txt_row.mmiss_p_mean_GeV = mmiss_p_mean;
        txt_row.mmiss_p_sigma_GeV = mmiss_p_sigma;
        TString txtout = Form("%s/summary_run%d.txt", outSummaryDir.Data(), run);
        write_summary_txt_file(txtout, txt_row, bg.n_accidentals_err);
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

        // Batch cleanup: windowed mass histograms (diag/side/full)
        for (auto *h: h_mgg_diag) if (h) safe_delete(h);
        for (auto *h: h_mgg_hor) if (h) safe_delete(h);
        for (auto *h: h_mgg_ver) if (h) safe_delete(h);

        for (auto& plot : cut_debug_plots) {
            if (plot.hist) safe_delete(plot.hist);
        }
        for (auto& plot : cut_debug_2d_plots) {
            if (plot.hist) safe_delete(plot.hist);
        }

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

    sw_total.Stop();
    logmsg(INFO, Form("ALL RUNS finished. Total runtime: %f s (real)", sw_total.RealTime()));
    
    if (write_global_csv) {
        logmsg(INFO, Form("Wrote global summary CSV to %s", global_csv.Data()));
    }
    
    // Exit cleanly to avoid ROOT cleanup segfaults
    gApplication->Terminate(0);
}
