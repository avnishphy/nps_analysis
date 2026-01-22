// // use Yaopeng's DataBase_production_runs_newBCMOffset.txt file for now.


// #ifndef NPS_LIVETIMES_H
// #define NPS_LIVETIMES_H

// // nps_livetimes.h
// // Compute livetimes from a per-run ROOT file containing trees "TSH" (scalers) and "T" (events).

// #include <TFile.h>
// #include <TTree.h>
// #include <TString.h>
// #include <TSystem.h>

// #include <iostream>
// #include <fstream>
// #include <cmath>
// #include <iomanip>

// namespace nps { namespace live {

// // Structure holding raw counts and computed livetimes
// struct LivetimeResult {
//     int run = -1;

//     // raw counts from scaler trees / TDC counts
//     double t_edtm_accepted = 0.0;
//     double t_hTRIG4_all_accepted = 0.0;
//     double t_hTRIG4_phy_accepted = 0.0;

//     double trig_accp_total = 0.0;       // H.hL1ACCP.scaler (accepted triggers)
//     double scaler_hTRIG4_total = 0.0;   // H.hTRIG4.scaler (trigger scaler)
//     double scaler_edtm_total = 0.0;     // H.EDTM.scaler (EDTM scaler)

//     double scaler_edtm_no_cut = 0.0;    // EDTM total without current cut
//     double trig_accp_no_cut = 0.0;      // accepted triggers w/o current cut

//     // beam-on fractions (current cut applied)
//     double beam_on_percent_edtm = 0.0;
//     double beam_on_percent_trig_accp = 0.0;

//     // computed livetime estimates (unitless: 0..1 typical)
//     double TLT_livetime_edtm = 0.0;     // total livetime from EDTM
//     double CLTA_livetime_tsh = 0.0;     // computer livetime all (TSH method)
//     double CLTA_livetime_tdc = 0.0;     // computer livetime all (TDC counts)
//     double CLTP_livetime_tdc = 0.0;     // computer livetime physics (TDC)

//     // metadata
//     double expected_current_uA = 0.0;
//     double current_tolerance_uA = 1.0;
//     double ps_factor = 1.0;

//     // convenient print
//     void print(std::ostream &out = std::cout) const {
//         out << std::fixed << std::setprecision(6);
//         out << "[LivetimeResult] run: " << run << "\n";
//         out << "  expected_current = " << expected_current_uA << " ± " << current_tolerance_uA << " uA\n";
//         out << "  scaler_edtm_total = " << scaler_edtm_total << "  scaler_edtm_no_cut = " << scaler_edtm_no_cut << "\n";
//         out << "  scaler_hTRIG4_total = " << scaler_hTRIG4_total << "  trig_accp_total = " << trig_accp_total << "\n";
//         out << "  t_edtm_accepted = " << t_edtm_accepted << "  t_hTRIG4_all_accepted = " << t_hTRIG4_all_accepted << "  t_hTRIG4_phy_accepted = " << t_hTRIG4_phy_accepted << "\n";
//         out << "  beam_on_percent_edtm = " << beam_on_percent_edtm << "  beam_on_percent_trig_accp = " << beam_on_percent_trig_accp << "\n";
//         out << "  TLT_livetime_edtm = " << TLT_livetime_edtm << "\n";
//         out << "  CLTA_livetime_tsh = " << CLTA_livetime_tsh << "  CLTA_livetime_tdc = " << CLTA_livetime_tdc << "  CLTP_livetime_tdc = " << CLTP_livetime_tdc << "\n";
//     }

//     // CSV line (header: see append_livetime_csv)
//     TString csv_line() const {
//         TString s = TString::Format("%d,%.6f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
//             run,
//             expected_current_uA,
//             t_edtm_accepted,
//             t_hTRIG4_all_accepted,
//             t_hTRIG4_phy_accepted,
//             trig_accp_total,
//             scaler_hTRIG4_total,
//             scaler_edtm_total,
//             beam_on_percent_edtm,
//             beam_on_percent_trig_accp,
//             TLT_livetime_edtm,
//             CLTA_livetime_tsh,
//             CLTA_livetime_tdc,
//             CLTP_livetime_tdc);
//         return s;
//     }
// };

// // Internal helper: safe branch lookup & SetBranchAddress for double arrays of length 1
// static inline bool set_branch_address_double(TTree* tree, const char* name, Double_t* addr) {
//     if (!tree) return false;
//     TBranch *b = tree->GetBranch(name);
//     if (!b) {
//         // Branch not found: try with dots replaced by underscores (some trees changed)
//         std::string alt(name);
//         for (char &c : alt) if (c == '.') c = '_';
//         b = tree->GetBranch(alt.c_str());
//         if (!b) return false;
//     }
//     tree->SetBranchStatus(name, 1);
//     tree->SetBranchAddress(name, addr);
//     return true;
// }

// // compute_livetimes_from_Ttrees:
// // Given pointer to the scaler tree (TSH) and event T tree (T), compute scalers / TDC counts
// // expected_current_uA: apply current cut when summing scalers (tolerance default 1.0 uA).
// // ps_factor: prescale factor to convert rates as in your python reference (default 1.0).
// static inline LivetimeResult compute_livetimes_from_Ttrees(TTree* tsh_tree, TTree* t_tree,
//                                                            int run_number = -1,
//                                                            double expected_current_uA = 0.0,
//                                                            double tolerance_uA = 1.0,
//                                                            double ps_factor = 1.0)
// {
//     LivetimeResult R;
//     R.run = run_number;
//     R.expected_current_uA = expected_current_uA;
//     R.current_tolerance_uA = tolerance_uA;
//     R.ps_factor = ps_factor;

//     if (!tsh_tree) {
//         std::cerr << "[nps::live] ERROR: tsh_tree is null\n";
//         return R;
//     }
//     if (!t_tree) {
//         std::cerr << "[nps::live] WARNING: t_tree is null (TDC counting will be zero)\n";
//     }

//     // Prepare single-value containers for SetBranchAddress
//     Double_t H_BCM4A_scalerCurrent = 0.0;
//     Double_t H_EDTM_scaler = 0.0;
//     Double_t H_hTRIG4_scaler = 0.0;
//     Double_t H_hL1ACCP_scaler = 0.0;
//     Double_t H_hEL_REAL_scaler = 0.0;

//     // T tree scalers (TDC raw times)
//     Double_t T_hms_hEDTM_tdcTimeRaw = 0.0;
//     Double_t T_hms_hTRIG4_tdcTimeRaw = 0.0;

//     // Set branch addresses if available
//     bool has_BCM4A = set_branch_address_double(tsh_tree, "H.BCM4A.scalerCurrent", &H_BCM4A_scalerCurrent);
//     bool has_EDTM  = set_branch_address_double(tsh_tree, "H.EDTM.scaler", &H_EDTM_scaler);
//     bool has_hTRIG4= set_branch_address_double(tsh_tree, "H.hTRIG4.scaler", &H_hTRIG4_scaler);
//     bool has_hL1ACCP=set_branch_address_double(tsh_tree, "H.hL1ACCP.scaler", &H_hL1ACCP_scaler);
//     bool has_hEL_REAL=set_branch_address_double(tsh_tree, "H.hEL_REAL.scaler", &H_hEL_REAL_scaler);

//     if (!has_EDTM && !has_hTRIG4 && !has_hL1ACCP) {
//         std::cerr << "[nps::live] ERROR: No required scaler branches found in TSH tree.\n";
//         return R;
//     }

//     if (t_tree) {
//         // T tree branch names as in your reference: T.hms.hEDTM_tdcTimeRaw / T.hms.hTRIG4_tdcTimeRaw
//         set_branch_address_double(t_tree, "T.hms.hEDTM_tdcTimeRaw", &T_hms_hEDTM_tdcTimeRaw);
//         set_branch_address_double(t_tree, "T.hms.hTRIG4_tdcTimeRaw", &T_hms_hTRIG4_tdcTimeRaw);
//     }

//     // Iterate over TSH scaler entries and compute differences between consecutive entries
//     Long64_t nentries_tsh = tsh_tree->GetEntries();
//     if (nentries_tsh <= 0) {
//         std::cerr << "[nps::live] WARNING: TSH tree has zero entries\n";
//         return R;
//     }

//     // previous values to compute deltas
//     double prev_EDTM = 0.0;
//     double prev_hTRIG4 = 0.0;
//     double prev_accp = 0.0;
//     double prev_hEL_REAL = 0.0;
//     bool first = true;

//     for (Long64_t i = 0; i < nentries_tsh; ++i) {
//         tsh_tree->GetEntry(i);
//         if (first) {
//             prev_EDTM = H_EDTM_scaler;
//             prev_hTRIG4 = H_hTRIG4_scaler;
//             prev_accp = H_hL1ACCP_scaler;
//             prev_hEL_REAL = H_hEL_REAL_scaler;
//             first = false;
//             continue;
//         }

//         double delta_EDTM = H_EDTM_scaler - prev_EDTM;
//         double delta_hTRIG4 = H_hTRIG4_scaler - prev_hTRIG4;
//         double delta_accp = H_hL1ACCP_scaler - prev_accp;
//         double delta_hEL_REAL = H_hEL_REAL_scaler - prev_hEL_REAL;

//         // accumulate no-cut totals
//         R.scaler_edtm_no_cut += (delta_EDTM > 0 ? delta_EDTM : 0);
//         R.trig_accp_no_cut += (delta_accp > 0 ? delta_accp : 0);

//         // apply current cut (if branch available, otherwise assume always accept)
//         bool current_ok = true;
//         if (has_BCM4A) {
//             current_ok = (std::fabs(H_BCM4A_scalerCurrent - expected_current_uA) <= tolerance_uA);
//         }

//         if (current_ok) {
//             if (delta_EDTM > 0) R.scaler_edtm_total += delta_EDTM;
//             if (delta_hTRIG4 > 0) R.scaler_hTRIG4_total += delta_hTRIG4;
//             if (delta_accp > 0) R.trig_accp_total += delta_accp;
//             if (delta_hEL_REAL > 0) R.scaler_edtm_total += 0.0; // keep for symmetry (hEL_REAL not used directly)
//             // note: we keep separate scalers but only use those needed below
//         }

//         prev_EDTM = H_EDTM_scaler;
//         prev_hTRIG4 = H_hTRIG4_scaler;
//         prev_accp = H_hL1ACCP_scaler;
//         prev_hEL_REAL = H_hEL_REAL_scaler;
//     }

//     // compute beam-on percentages (avoid division by zero)
//     if (R.scaler_edtm_no_cut > 0.0) R.beam_on_percent_edtm = (R.scaler_edtm_total / R.scaler_edtm_no_cut);
//     else R.beam_on_percent_edtm = 0.0;

//     if (R.trig_accp_no_cut > 0.0) R.beam_on_percent_trig_accp = (R.trig_accp_total / R.trig_accp_no_cut);
//     else R.beam_on_percent_trig_accp = 0.0;

//     // T tree counts (TDC) — count how many events have nonzero raw times
//     if (t_tree) {
//         Long64_t nentries_t = t_tree->GetEntries();
//         for (Long64_t i = 0; i < nentries_t; ++i) {
//             t_tree->GetEntry(i);
//             if (T_hms_hEDTM_tdcTimeRaw != 0.0) R.t_edtm_accepted += 1.0;
//             if (T_hms_hTRIG4_tdcTimeRaw != 0.0) R.t_hTRIG4_all_accepted += 1.0;
//             if ( (T_hms_hTRIG4_tdcTimeRaw != 0.0) && (T_hms_hEDTM_tdcTimeRaw == 0.0) ) R.t_hTRIG4_phy_accepted += 1.0;
//         }
//     }

//     // Compute livetime estimates. Use safe guards for zero denominators.
//     // CLTA (TSH): trig_accp_total / scaler_hTRIG4_total
//     if (R.scaler_hTRIG4_total > 0.0) {
//         R.CLTA_livetime_tsh = (R.trig_accp_total / R.scaler_hTRIG4_total) * ps_factor;
//         // TDC-based computer all (accepted T.hTRIG4 / scaler_hTRIG4)
//         R.CLTA_livetime_tdc = (R.t_hTRIG4_all_accepted / R.scaler_hTRIG4_total) * ps_factor * R.beam_on_percent_trig_accp;
//     } else {
//         R.CLTA_livetime_tsh = 0.0;
//         R.CLTA_livetime_tdc = 0.0;
//     }

//     // CLTP (physics livetime from TDC): t_hTRIG4_phy_accepted / (scaler_hTRIG4_total - scaler_edtm_total)
//     double denom_phys = R.scaler_hTRIG4_total - R.scaler_edtm_total;
//     if (denom_phys > 0.0) {
//         R.CLTP_livetime_tdc = (R.t_hTRIG4_phy_accepted / denom_phys) * ps_factor * R.beam_on_percent_trig_accp;
//     } else {
//         R.CLTP_livetime_tdc = 0.0;
//     }

//     // TLT (EDTM): t_edtm_accepted / (scaler_edtm_total / beam_on_percent_edtm)
//     // Be careful: scaler_edtm_total/beam_on_percent_edtm = effective EDTM delivered during accepted beam-on time.
//     if (R.scaler_edtm_total > 0.0 && R.beam_on_percent_edtm > 0.0) {
//         double edtm_delivered_in_beam_on = (R.scaler_edtm_total / R.beam_on_percent_edtm);
//         if (edtm_delivered_in_beam_on > 0.0) {
//             R.TLT_livetime_edtm = (R.t_edtm_accepted / edtm_delivered_in_beam_on) * ps_factor;
//         } else {
//             R.TLT_livetime_edtm = 0.0;
//         }
//     } else {
//         R.TLT_livetime_edtm = 0.0;
//     }

//     // store small metadata again
//     R.expected_current_uA = expected_current_uA;
//     R.current_tolerance_uA = tolerance_uA;
//     R.ps_factor = ps_factor;

//     return R;
// }

// // compute_livetimes_from_tfile:
// // Given an already-open TFile* or a filename, load "TSH" and "T" trees and call the function above.
// // If tfile is null, will attempt to open filename.
// static inline LivetimeResult compute_livetimes_from_tfile(TFile* tfile,
//                                                           int run_number = -1,
//                                                           double expected_current_uA = 0.0,
//                                                           double tolerance_uA = 1.0,
//                                                           double ps_factor = 1.0)
// {
//     LivetimeResult R;
//     if (!tfile) {
//         std::cerr << "[nps::live] ERROR: compute_livetimes_from_tfile requires non-null TFile*\n";
//         return R;
//     }
//     TTree* tsh = nullptr;
//     TTree* t = nullptr;

//     // Try standard names
//     tsh = dynamic_cast<TTree*>(tfile->Get("TSH"));
//     if (!tsh) tsh = dynamic_cast<TTree*>(tfile->Get("tsh")); // try lowercase
//     t = dynamic_cast<TTree*>(tfile->Get("T"));
//     if (!t) t = dynamic_cast<TTree*>(tfile->Get("t"));

//     if (!tsh) {
//         std::cerr << "[nps::live] WARNING: TSH tree not found in file; livetime calculation aborted\n";
//         return R;
//     }

//     R = compute_livetimes_from_Ttrees(tsh, t, run_number, expected_current_uA, tolerance_uA, ps_factor);
//     return R;
// }

// // Convenience wrapper that opens a file path and computes livetimes.
// // Returns a LivetimeResult (run can be set with run_number argument).
// static inline LivetimeResult compute_livetimes_from_filepath(const TString &filepath,
//                                                              int run_number = -1,
//                                                              double expected_current_uA = 0.0,
//                                                              double tolerance_uA = 1.0,
//                                                              double ps_factor = 1.0)
// {
//     LivetimeResult R;
//     if (gSystem->AccessPathName(filepath)) {
//         std::cerr << "[nps::live] ERROR: cannot open file " << filepath.Data() << "\n";
//         return R;
//     }
//     TFile *f = TFile::Open(filepath, "READ");
//     if (!f || f->IsZombie()) {
//         std::cerr << "[nps::live] ERROR: cannot open file " << filepath.Data() << "\n";
//         if (f) { f->Close(); }
//         return R;
//     }
//     R = compute_livetimes_from_tfile(f, run_number, expected_current_uA, tolerance_uA, ps_factor);
//     f->Close();
//     return R;
// }

// // Append header & a single result to CSV (creates file with header if not present).
// static inline bool append_livetime_csv(const TString &csvpath, const LivetimeResult &R) {
//     // header:
//     const char *header = "run,current,"
//                          "t_edtm_accepted,t_hTRIG4_all_accepted,t_hTRIG4_phy_accepted,"
//                          "trig_accp_total,scaler_hTRIG4_total,scaler_edtm_total,"
//                          "beam_on_percent_edtm,beam_on_percent_trig_accp,"
//                          "TLT_livetime_edtm,CLTA_livetime_tsh,CLTA_livetime_tdc,CLTP_livetime_tdc\n";

//     bool file_exists = !gSystem->AccessPathName(csvpath);
//     std::ofstream out(csvpath.Data(), std::ios::app);
//     if (!out.is_open()) {
//         std::cerr << "[nps::live] ERROR: Could not open CSV for append: " << csvpath.Data() << "\n";
//         return false;
//     }
//     if (!file_exists) {
//         out << header;
//     }
//     out << R.csv_line().Data();
//     out.close();
//     return true;
// }

// }} // namespace nps::live

// #endif // NPS_LIVETIMES_H





// ===========================================================================================================
// Reading the livetimes directly from Yaopeng's DataBase_production_runs_newBCMOffset.txt


#pragma once
#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <limits>

namespace nps {
namespace detail_livetime {

// helper to split a string by whitespace
inline std::vector<std::string> split_ws(const std::string &line) {
    std::istringstream ss(line);
    std::vector<std::string> toks;
    std::string tok;
    while (ss >> tok) toks.push_back(tok);
    return toks;
}

// parse a column from the file into a map
inline bool parse_column_map(const std::string &path,
                             const std::string &column_name,
                             std::unordered_map<int,double> &out_map,
                             bool quiet = false)
{
    std::ifstream ifs(path);
    if (!ifs.is_open()) {
        if (!quiet) std::cerr << "nps_livetimes: cannot open " << path << std::endl;
        return false;
    }

    std::string line;
    bool header_found = false;
    std::vector<std::string> header_tokens;
    int idx_run = -1, idx_col = -1;

    while (std::getline(ifs,line)) {
        bool allws = true;
        for (char c:line) if(!isspace((unsigned char)c)) { allws=false; break; }
        if(allws) continue;

        auto toks = split_ws(line);
        if(toks.empty()) continue;

        if(!header_found) {
            if(toks[0]=="RunNo") {
                header_found=true;
                header_tokens=toks;
                for(size_t i=0;i<toks.size();++i){
                    if(header_tokens[i]=="RunNo") idx_run=i;
                    if(header_tokens[i]==column_name) idx_col=i;
                }
                if(idx_run==-1 || idx_col==-1){
                    if(!quiet) std::cerr<<"nps_livetimes: header missing required columns\n";
                    return false;
                }
            }
            continue;
        } else {
            if(toks.size()<=static_cast<size_t>(std::max(idx_run,idx_col))) continue;
            try{
                int run=std::stoi(toks[idx_run]);
                double val=std::stod(toks[idx_col]);
                out_map[run]=val;
            } catch(...){ continue; }
        }
    }

    if(out_map.empty() && !quiet){
        std::cerr<<"nps_livetimes: no data found for column '"<<column_name<<"'\n";
        return false;
    }
    return true;
}

inline const std::string kDBPath_livetime = "config/DataBase_production_runs_newBCMOffset.txt";

} // namespace detail_livetime

inline bool getCPU_LT(int run, double &out_val){
    static std::unordered_map<int,double> s_map;
    static bool s_loaded=false;
    if(!s_loaded){
        s_loaded = detail_livetime::parse_column_map(detail_livetime::kDBPath_livetime,
                                                    "CPU_LT", s_map, false);
    }
    auto it = s_map.find(run);
    if(it==s_map.end()) return false;
    out_val = it->second;
    return true;
}

inline double getCPU_LT_or_nan(int run){
    double v;
    if(getCPU_LT(run,v)) return v;
    return std::numeric_limits<double>::quiet_NaN();
}

} // namespace nps
