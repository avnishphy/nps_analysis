// ============================================================================
// File : nps_analysis_full.C
// Author: Avnish Singh
// Revised: ChatGPT (2025-11-07)
// Purpose: Full diagnostic pipeline for NPS pi0 analysis with robust outputs:
//   - PNG canvases (timing plane, pi0 overlays, mm vs mgg scatter, cluster diagnostics)
//   - per-run ROOT file with histograms, fit canvases, TParameters
//   - per-run CSV summary and appended global CSV
//   - per-run human-readable TXT summary
//
// Dependencies (unchanged): utils.C (logmsg, trim), nps_helper.h (invariant_mass_pi0,
//    fit_pi0_peak, packClusters, nps_spatial_energy_cuts, missing_mass_proton_pi0,
//    missing_mass_dvcs, constants), nps_time_bg.h (estimate_coincidence_background_default,
//    default windows, integral_and_area_TH2).
//
// Usage: in ROOT
//    .x nps_analysis_full.C+()
// or
//    nps_analysis_full("path/to/skim/","output/plots/","config/runlist.txt");
//
// Notes / recommendations:
// - Verify branch names match your tree. The macro checks essential branches but doesn't
//   replace missing/renamed branches automatically.
// - Increase MAX_CLUS only if your events truly require it.
// ============================================================================

#include "utils.C"
#include "nps_helper.h"
#include "nps_time_bg.h"

#include <TFile.h>
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

using namespace std;

// ------------------------------------------------------------------
// Configurable limits & constants
// ------------------------------------------------------------------
constexpr int MAX_CLUS = 20;            // maximum cluster array size handled
constexpr double DEFAULT_TIME_THRESH_NS = 2.0;
constexpr double DEFAULT_TIME_WINDOW_WRT_150 = 10.0; // 150ns+-(DEFAULT_TIME_WINDOW_WRT_150)
constexpr double EBEAM_DEFAULT = 10.538;  // GeV (set to experimental beam energy)
constexpr double MM_PROTON_SIGMA = 0.05;  // GeV, used for DVCS veto
constexpr int NPRINT_PROGRESS = 16384;

// ------------------------------------------------------------------
// Read run list (supports comments beginning '#')
// ------------------------------------------------------------------
vector<int> readRunList(const string &fname) {
    vector<int> runs;
    ifstream in(fname);
    if (!in.is_open()) {
        logmsg(WARN, "Cannot open runlist: " + fname);
        return runs;
    }
    string line;
    while (getline(in, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        try { runs.push_back(stoi(line)); }
        catch (...) { logmsg(WARN, "Skipping invalid runlist line: " + line); }
    }
    return runs;
}

// ------------------------------------------------------------------
// Choose best pi0 pair: closest to pi0 mass AND within time threshold.
// Returns indices into clus arrays (ia, ib) or (-1,-1).
// ------------------------------------------------------------------
pair<int,int> choose_best_pair_closest_pi0(const vector<int> &good_idx,
                                           double *clusE, double *clusX, double *clusY, double *clusT,
                                           double z_nps,
                                           double target = nps::kPi0Mass_GeV,
                                           double time_thresh_ns = DEFAULT_TIME_THRESH_NS)
{
    int best_i = -1, best_j = -1;
    double best_diff = 1e9;
    double best_totalE = -1.0;
    const int n = static_cast<int>(good_idx.size());
    if (n < 2) return { -1, -1 };

    for (int a = 0; a < n; ++a) {
        for (int b = a + 1; b < n; ++b) {
            const int ia = good_idx[a];
            const int ib = good_idx[b];
            if (ia < 0 || ib < 0) continue;
            const double dt = fabs(clusT[ia] - clusT[ib]);
            if (dt > time_thresh_ns) continue;
            const double m = nps::invariant_mass_pi0(clusE[ia], clusE[ib],
                                                    clusX[ia], clusX[ib],
                                                    clusY[ia], clusY[ib],
                                                    z_nps);
            const double d = fabs(m - target);
            if (d < best_diff - 1e-12) {
                best_diff = d;
                best_i = ia; best_j = ib;
                best_totalE = clusE[ia] + clusE[ib];
            } else if (fabs(d - best_diff) < 1e-12) {
                double tot = clusE[ia] + clusE[ib];
                if (tot > best_totalE) {
                    best_i = ia; best_j = ib; best_totalE = tot;
                }
            }
        }
    }
    return {best_i, best_j};
}

// ------------------------------------------------------------------
// Helper: write CSV header for global summary
// ------------------------------------------------------------------
void write_global_csv_header(const TString &path) {
    ofstream f(path.Data(), ios::out);
    f << "run,total_events,pass_hms,n_selected,coin_raw,acc_est,acc_err,"
      << "pi0_mean,pi0_mean_err,pi0_sigma,pi0_sigma_err,pi0_signal,pi0_signal_err,pi0_bg,pi0_bg_err,pi0_chi2,pi0_ndf,"
      << "coin_roi_yield,acc_roi_yield,acc_est_full_scaled,signal_sub_full,signal_sub_full_err,acc_est_hv_avg,signal_sub_hv_avg,signal_sub_hv_avg_err\n";
    f.close();
}

// ------------------------------------------------------------------
// Main macro: single-file, full diagnostics
// ------------------------------------------------------------------
void nps_analysis(const TString &skimDir_in="output/skimmed/",
                       const TString &outPlotDir_in="output/plots/",
                       const TString &runlistFile="config/runlist_x60_4b.txt",
                       const double Ebeam = EBEAM_DEFAULT)
{
    TStopwatch sw_total; sw_total.Start();
    logmsg(INFO, "=========== NPS π0 FULL diagnostic analysis ===========");

    TString skimDir = skimDir_in.EndsWith("/") ? skimDir_in : skimDir_in + "/";
    TString outPlotDir = outPlotDir_in.EndsWith("/") ? outPlotDir_in : outPlotDir_in + "/";
    gSystem->mkdir(outPlotDir, true);
    gSystem->mkdir("output", true);

    vector<int> runs = readRunList(runlistFile.Data());
    if (runs.empty()) {
        logmsg(ERROR, std::string("No runs found in runlist: ") + runlistFile.Data());
        return;
    }

    // prepare global CSV
    TString global_csv = outPlotDir + "/summary_all_runs.csv";
    write_global_csv_header(global_csv);

    // per-run loop
    for (int run : runs) {
        TStopwatch sw_run; sw_run.Start();
        TString infile = Form("%sskim_run%d.root", skimDir.Data(), run);
        if (gSystem->AccessPathName(infile)) {
            logmsg(WARN, Form("Skipping run %d: file not found: %s", run, infile.Data()));
            continue;
        }
        logmsg(INFO, Form("Processing run %d (file: %s)", run, infile.Data()));

        TFile *f = TFile::Open(infile, "READ");
        if (!f || f->IsZombie()) { logmsg(ERROR, Form("Error opening %s", infile.Data())); if (f) f->Close(); continue; }

        TTree *T = dynamic_cast<TTree*>(f->Get("T"));
        if (!T) { logmsg(ERROR, Form("Tree 'T' not found in %s", infile.Data())); f->Close(); continue; }

        // -------------------------
        // Branch variables (safe sizes)
        // -------------------------
        Double_t HgtrX=0, HgtrY=0, HgtrTh=0, HgtrPh=0, hdelta=0, HgtrP=0, hreactz=0, hcernpeSum=0, hcaletotnorm=0;
        Double_t HgtrPx=0, HgtrPy=0, HgtrPz=0;
        Double_t edtmtdc=0;
        static Double_t clusE[MAX_CLUS], clusX[MAX_CLUS], clusY[MAX_CLUS], clusT[MAX_CLUS];
        Double_t nclust_dbl = 0;
        Double_t BCM2_scalerCurrent = 0;

        // set branch status / addresses (enable only needed branches)
        T->SetBranchStatus("*", 0);
        // adjust these names if your tree uses other labels
        if (T->GetBranch("H.gtr.x")) { T->SetBranchStatus("H.gtr.x", 1); T->SetBranchAddress("H.gtr.x", &HgtrX); }
        if (T->GetBranch("H.gtr.y")) { T->SetBranchStatus("H.gtr.y", 1); T->SetBranchAddress("H.gtr.y", &HgtrY); }
        if (T->GetBranch("H.gtr.p")) { T->SetBranchStatus("H.gtr.p", 1); T->SetBranchAddress("H.gtr.p", &HgtrP); }
        if (T->GetBranch("H.gtr.px")) { T->SetBranchStatus("H.gtr.px", 1); T->SetBranchAddress("H.gtr.px", &HgtrPx); }
        if (T->GetBranch("H.gtr.py")) { T->SetBranchStatus("H.gtr.py", 1); T->SetBranchAddress("H.gtr.py", &HgtrPy); }
        if (T->GetBranch("H.gtr.pz")) { T->SetBranchStatus("H.gtr.pz", 1); T->SetBranchAddress("H.gtr.pz", &HgtrPz); }
        if (T->GetBranch("H.gtr.dp")) { T->SetBranchStatus("H.gtr.dp", 1); T->SetBranchAddress("H.gtr.dp", &hdelta); }
        if (T->GetBranch("H.gtr.th")) { T->SetBranchStatus("H.gtr.th", 1); T->SetBranchAddress("H.gtr.th", &HgtrTh); }
        if (T->GetBranch("H.gtr.ph")) { T->SetBranchStatus("H.gtr.ph", 1); T->SetBranchAddress("H.gtr.ph", &HgtrPh); }
        if (T->GetBranch("H.react.z")) { T->SetBranchStatus("H.react.z", 1); T->SetBranchAddress("H.react.z", &hreactz); }
        if (T->GetBranch("H.cer.npeSum")) { T->SetBranchStatus("H.cer.npeSum", 1); T->SetBranchAddress("H.cer.npeSum", &hcernpeSum); }
        if (T->GetBranch("H.cal.etotnorm")) { T->SetBranchStatus("H.cal.etotnorm", 1); T->SetBranchAddress("H.cal.etotnorm", &hcaletotnorm); }
        if (T->GetBranch("T.hms.hEDTM_tdcTimeRaw")) { T->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1); T->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc); }
        if (T->GetBranch("H.BCM2.scalerCurrent")) { T->SetBranchStatus("H.BCM2.scalerCurrent", 1); T->SetBranchAddress("H.BCM2.scalerCurrent", &BCM2_scalerCurrent); }

        if (T->GetBranch("NPS.cal.nclust")) { T->SetBranchStatus("NPS.cal.nclust", 1); T->SetBranchAddress("NPS.cal.nclust", &nclust_dbl); }
        if (T->GetBranch("NPS.cal.clusE")) { T->SetBranchStatus("NPS.cal.clusE", 1); T->SetBranchAddress("NPS.cal.clusE", &clusE); }
        if (T->GetBranch("NPS.cal.clusX")) { T->SetBranchStatus("NPS.cal.clusX", 1); T->SetBranchAddress("NPS.cal.clusX", &clusX); }
        if (T->GetBranch("NPS.cal.clusY")) { T->SetBranchStatus("NPS.cal.clusY", 1); T->SetBranchAddress("NPS.cal.clusY", &clusY); }
        if (T->GetBranch("NPS.cal.clusT")) { T->SetBranchStatus("NPS.cal.clusT", 1); T->SetBranchAddress("NPS.cal.clusT", &clusT); }

        Long64_t nentries = T->GetEntries();
        cout << "Run " << run << " entries: " << nentries << endl;
        
        double run_current = -1.0;

        if (T->GetBranch("H.BCM2.scalerCurrent")) {
            double BCM2_scalerCurrent = 0;
            T->SetBranchStatus("H.BCM2.scalerCurrent", 1);
            T->SetBranchAddress("H.BCM2.scalerCurrent", &BCM2_scalerCurrent);

            TH1D hCurrent("hCurrent", "Beam Current;I (uA);Counts", 200, 0, 100);
            Long64_t nentries = T->GetEntries();

            for (Long64_t i = 0; i < nentries; i++) {
                T->GetEntry(i);
                if (BCM2_scalerCurrent > 2.0)
                    hCurrent.Fill(BCM2_scalerCurrent);
            }

            if (hCurrent.GetEntries() > 0) {
                int maxBin = hCurrent.GetMaximumBin();
                run_current = hCurrent.GetXaxis()->GetBinCenter(maxBin); // ✅ correct
            } else {
                run_current = 0.0;
            }
        }




        // -------------------------
        // Histograms (unique names per run)
        // -------------------------
        auto name = [&](const char* base){ return TString::Format("%s_run%d", base, run); };

        TH1D *h_nclusters = new TH1D(name("h_nclusters"), "NPS clusters per event;N_{clus};Events", 21, -0.5, 20.5);
        TH1D *h_clustE = new TH1D(name("h_clustE"), "Cluster energy;E_{clus} [GeV];Counts", 200, 0.0, 8.0);
        TH1D *h_clustT = new TH1D(name("h_clustT"), "Cluster time; t [ns];Counts", 200, 120, 180);
        TH2D *h_clustE_vs_T = new TH2D(name("h_clustE_vs_T"), "Cluster E vs t; t [ns];E [GeV]", 200, 120, 180, 200, 0, 8);

        TH1D *h_mpi0_all = new TH1D(name("h_mpi0_all"), "Invariant mass (all best-pairs);M_{#gamma#gamma} [GeV];Events", 200, 0.0, 0.4);
        TH1D *h_mpi0_2 = new TH1D(name("h_mpi0_2"), "Invariant mass (2-cluster);M [GeV];Events", 200, 0.0, 0.4);
        TH1D *h_mpi0_3 = new TH1D(name("h_mpi0_3"), "Invariant mass (3-cluster best pair);M [GeV];Events", 200, 0.0, 0.4);
        TH1D *h_mpi0_4 = new TH1D(name("h_mpi0_4"), "Invariant mass (4-cluster best pair);M [GeV];Events", 200, 0.0, 0.4);

        TH1D *h_mmiss_2 = new TH1D(name("h_mmiss_2"), "Missing mass (2-cluster);M_{miss} [GeV];Events", 200, 0.0, 2.0);
        TH1D *h_mmiss_3 = new TH1D(name("h_mmiss_3"), "Missing mass (3-cluster);M_{miss} [GeV];Events", 200, 0.0, 2.0);
        TH1D *h_mmiss_4 = new TH1D(name("h_mmiss_4"), "Missing mass (4-cluster);M_{miss} [GeV];Events", 200, 0.0, 2.0);
        TH1D *h_mmiss_dvcs = new TH1D(name("h_mmiss_dvcs"), "Missing mass (DVCS-like single photon);M_{miss} [GeV];Events", 200, 0.0, 2.0);

        const double t_min = 140.0, t_max = 160.0;
        const int nbins_t = 100;
        TH2D *h_t1_t2 = new TH2D(name("h_t1_t2"), "t1 (y) vs t2 (x);t2 [ns];t1 [ns]", nbins_t, t_min, t_max, nbins_t, t_min, t_max);
        TH1D *h_t1_proj = new TH1D(name("h_t1_proj"), "t1 projection; t1 [ns];Entries", nbins_t, t_min, t_max);
        TH1D *h_t2_proj = new TH1D(name("h_t2_proj"), "t2 projection; t2 [ns];Entries", nbins_t, t_min, t_max);

        // coin / acc histos and per-box mgg histograms
        TH1D *h_m_pi0_coin = new TH1D(name("h_m_pi0_coin"), "Pi0 mass (Coincidence);M [GeV];Counts", 200, 0.0, 0.4);
        TH1D *h_m_pi0_acc = new TH1D(name("h_m_pi0_acc"), "Pi0 mass (Accidentals - outside coin);M [GeV];Counts", 200, 0.0, 0.4);

        // Per-window mgg histograms (diagonal, horizontal, vertical, full boxes)
        vector<pair<double,double>> diag_windows = nps_bg::default_diag_windows();
        vector<pair<double,double>> side_windows = nps_bg::default_side_windows();
        auto coin_win = nps_bg::default_coin_window();
        auto full1_t1 = nps_bg::default_full_acc1_t1();
        auto full1_t2 = nps_bg::default_full_acc1_t2();
        auto full2_t1 = nps_bg::default_full_acc2_t1();
        auto full2_t2 = nps_bg::default_full_acc2_t2();

        vector<TH1D*> h_mgg_diag;
        for (size_t i=0;i<diag_windows.size();++i) {
            h_mgg_diag.push_back(new TH1D(name(TString::Format("h_mgg_diag%d", (int)i)).Data(),
                                         TString::Format("m_{#gamma#gamma} diag %d;M [GeV];Counts",(int)i).Data(),
                                         200,0.0,0.4));
        }
        vector<TH1D*> h_mgg_hor; // horizontal sideboxes (t1=coin, t2 in side)
        vector<TH1D*> h_mgg_ver; // vertical sideboxes (t2=coin, t1 in side)
        for (size_t i=0;i<side_windows.size();++i) {
            h_mgg_hor.push_back(new TH1D(name(TString::Format("h_mgg_hor%d",(int)i)).Data(),
                                         TString::Format("m_{#gamma#gamma} hor %d;M [GeV];Counts",(int)i).Data(),200,0.0,0.4));
            h_mgg_ver.push_back(new TH1D(name(TString::Format("h_mgg_ver%d",(int)i)).Data(),
                                         TString::Format("m_{#gamma#gamma} ver %d;M [GeV];Counts",(int)i).Data(),200,0.0,0.4));
        }
        TH1D *h_mgg_full1 = new TH1D(name("h_mgg_full1"), "m_{#gamma#gamma} fullbox1;M [GeV];Counts", 200,0.0,0.4);
        TH1D *h_mgg_full2 = new TH1D(name("h_mgg_full2"), "m_{#gamma#gamma} fullbox2;M [GeV];Counts", 200,0.0,0.4);

        // mm vs mgg scatter
        TH2D *h_mmiss_vs_mgg = new TH2D(name("h_mmiss_vs_mgg"), "M_{miss} vs M_{#gamma#gamma};M_{#gamma#gamma} [GeV];M_{miss} [GeV]",
                                       200,0.0,0.4,200,0.0,2.0);

        // bookkeeping counters
        Long64_t n_total = 0, n_pass_hms = 0, n_ge2_hms = 0;
        Long64_t n_mult2_hms = 0, n_mult3_hms = 0, n_mult4_hms = 0;
        Long64_t n_mult2_hms_nps = 0, n_mult3_hms_nps = 0, n_mult4_hms_nps = 0;
        Long64_t n_dvcs_flagged = 0, n_selected_for_analysis = 0;

        // -------------------------
        // Event loop
        // -------------------------
        for (Long64_t ev=0; ev<nentries; ++ev) {
            if ((ev & (NPRINT_PROGRESS-1)) == 0) {
                cout << " run " << run << " event " << ev << " / " << nentries << "\r" << flush;
            }
            T->GetEntry(ev);
            ++n_total;

            // HMS electron selection using your helper
            if (!nps::hms_electron_cuts(edtmtdc, hdelta, HgtrTh, HgtrPh, hcernpeSum, hcaletotnorm, hreactz)) continue;
            ++n_pass_hms;

            // safe cast of cluster count
            int nclust = static_cast<int>(lrint(nclust_dbl));
            if (nclust < 0) nclust = 0;
            if (nclust > MAX_CLUS) nclust = MAX_CLUS;
            h_nclusters->Fill(nclust);
            if (nclust < 2) continue;

            ++n_ge2_hms;
            if (nclust == 2) ++n_mult2_hms;
            else if (nclust == 3) ++n_mult3_hms;
            else if (nclust == 4) ++n_mult4_hms;

            // per-cluster diagnostics
            for (int i=0;i<nclust;++i) {
                h_clustE->Fill(clusE[i]);
                h_clustT->Fill(clusT[i]);
                h_clustE_vs_T->Fill(clusT[i], clusE[i]);
            }

            // pack clusters (your helper) - returns number after packing
            const int n_after = nps::packClusters(clusE, clusX, clusY, clusT, nclust);

            // select good clusters by spatial & energy cuts
            vector<int> good_idx; good_idx.reserve(8);
            for (int i=0;i<n_after;++i) {
                if (nps::nps_spatial_energy_cuts(clusE[i], clusX[i], clusY[i], clusT[i], DEFAULT_TIME_WINDOW_WRT_150))
                    good_idx.push_back(i);
            }

            if (good_idx.size() < 2) continue;
            if (good_idx.size() > 4) {
                sort(good_idx.begin(), good_idx.end(), [&](int a,int b){ return clusE[a] > clusE[b]; });
                good_idx.resize(4);
            }

            if (good_idx.size() == 2) ++n_mult2_hms_nps;
            else if (good_idx.size() == 3) ++n_mult3_hms_nps;
            else if (good_idx.size() == 4) ++n_mult4_hms_nps;

            // DVCS veto (single-photon events)
            const double px_e = HgtrPx, py_e = HgtrPy, pz_e = HgtrPz;
            const double p_e_mom = sqrt(max(0.0, px_e*px_e + py_e*py_e + pz_e*pz_e));
            const double Ee = sqrt(max(0.0, p_e_mom*p_e_mom + nps::kElectronMass_GeV*nps::kElectronMass_GeV));
            bool flagged_dvcs = false;
            for (int idx : good_idx) {
                double mm_single = nps::missing_mass_dvcs(Ebeam, Ee, px_e, py_e, pz_e, clusE[idx], clusX[idx], clusY[idx], nps::kDefaultZ_NPS_cm, -17.51);
                if (fabs(mm_single - nps::kProtonMass_GeV) < 3.0 * MM_PROTON_SIGMA) {
                    flagged_dvcs = true;
                    h_mmiss_dvcs->Fill(mm_single);
                    break;
                }
            }
            if (flagged_dvcs) { ++n_dvcs_flagged; continue; }

            // select pair for pi0
            int sel_i=-1, sel_j=-1;
            if (good_idx.size() == 2) { sel_i = good_idx[0]; sel_j = good_idx[1]; }
            else {
                auto pr = choose_best_pair_closest_pi0(good_idx, clusE, clusX, clusY, clusT, nps::kDefaultZ_NPS_cm, nps::kPi0Mass_GeV, DEFAULT_TIME_THRESH_NS);
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

            // missing mass
            const double mm_p = nps::missing_mass_proton_pi0(Ebeam, Ee, px_e, py_e, pz_e,
                                                            clusE[sel_i], clusE[sel_j],
                                                            clusX[sel_i], clusY[sel_i],
                                                            clusX[sel_j], clusY[sel_j],
                                                            nps::kDefaultZ_NPS_cm, -17.51);
            if (good_idx.size()==2) h_mmiss_2->Fill(mm_p);
            else if (good_idx.size()==3) h_mmiss_3->Fill(mm_p);
            else if (good_idx.size()==4) h_mmiss_4->Fill(mm_p);

            h_mmiss_vs_mgg->Fill(mgg, mm_p);
            ++n_selected_for_analysis;

            // determine which timing box this event belongs to (coin/diag/side/full)
            const double coin_time = 0.5*(t1 + t2);
            bool in_coin = (coin_time > coin_win.first && coin_time < coin_win.second);
            if (in_coin) h_m_pi0_coin->Fill(mgg);
            else h_m_pi0_acc->Fill(mgg);

            // check diagonal windows
            bool matched_diag=false;
            for (size_t i=0;i<diag_windows.size();++i) {
                const auto &w = diag_windows[i];
                if (t1> w.first && t1 < w.second && t2 > w.first && t2 < w.second) {
                    h_mgg_diag[i]->Fill(mgg);
                    matched_diag=true;
                    break;
                }
            }

            // horizontal sidebands (t1 in coin, t2 in side windows)
            bool matched_hor=false, matched_ver=false;
            for (size_t i=0;i<side_windows.size();++i) {
                const auto &w = side_windows[i];
                if ( (t1 > coin_win.first && t1 < coin_win.second) && (t2 > w.first && t2 < w.second) ) {
                    h_mgg_hor[i]->Fill(mgg); matched_hor=true; break;
                }
            }
            // vertical sidebands (t2 in coin, t1 in side windows)
            for (size_t i=0;i<side_windows.size();++i) {
                const auto &w = side_windows[i];
                if ( (t2 > coin_win.first && t2 < coin_win.second) && (t1 > w.first && t1 < w.second) ) {
                    h_mgg_ver[i]->Fill(mgg); matched_ver=true; break;
                }
            }

            // full accidental boxes
            if ( (t1 > full1_t1.first && t1 < full1_t1.second && t2 > full1_t2.first && t2 < full1_t2.second) ) {
                h_mgg_full1->Fill(mgg);
            }
            if ( (t1 > full2_t1.first && t1 < full2_t1.second && t2 > full2_t2.first && t2 < full2_t2.second) ) {
                h_mgg_full2->Fill(mgg);
            }

        } // end event loop

        cout << endl;

        // -------------------------
        // Summaries & background estimate
        // -------------------------
        nps_bg::CoincidenceBGResult bg = nps_bg::estimate_coincidence_background_default(h_t1_t2);
        logmsg(INFO, Form("Run %d: coin_raw=%.1f estimated_accidentals=%.3f +- %.3f", run, bg.n_coin_raw, bg.n_accidentals, bg.n_accidentals_err));

        // Fit pi0 peak on h_mpi0_all (keeps your helper which can draw + return FitResult)
        nps::FitResult fitres = nps::fit_pi0_peak(h_mpi0_all, 0.02, 0.30, true, outPlotDir, run);

        // Save timing plane canvas with boxes & annotations
        gStyle->SetOptStat(0);
        TCanvas *c_t12 = new TCanvas(name("c_t12"), "t1 vs t2", 1000,900);
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
        TBox *b_full2 = new TBox(full2_t1.first, full2_t2.first, full2_t1.second, full2_t2.second); b_full2->SetLineColor(kOrange+7); b_full2->SetLineWidth(2); b_full2->SetLineStyle(4); b_full2->SetFillStyle(0); b_full2->Draw("same");

        TLatex txt; txt.SetNDC(); txt.SetTextSize(0.025);
        txt.DrawLatex(0.02, 0.96, Form("Run %d  Coin raw = %.0f  Acc est = %.3f #pm %.3f", run, bg.n_coin_raw, bg.n_accidentals, bg.n_accidentals_err));
        txt.DrawLatex(0.02, 0.92, Form("Coin area (ns^{2}) = %.3f", bg.area_coin));

        double y_text=0.88;
        for (size_t i=0;i<diag_windows.size();++i) {
            auto pr = nps_bg::integral_and_area_TH2(h_t1_t2, diag_windows[i].first, diag_windows[i].second, diag_windows[i].first, diag_windows[i].second);
            double raw = pr.first; double area = pr.second;
            double norm = (area>0) ? raw * (bg.area_coin / area) : 0.0;
            txt.DrawLatex(0.02, y_text, Form("Diag[%g,%g] raw=%.0f norm=%.2f", diag_windows[i].first, diag_windows[i].second, raw, norm)); y_text -= 0.025;
        }

        y_text -= 0.01; txt.DrawLatex(0.02, y_text, "Horizontal sidebands (t1 in coin):"); y_text -= 0.025;
        for (size_t i=0;i<side_windows.size();++i) {
            auto pr = nps_bg::integral_and_area_TH2(h_t1_t2, side_windows[i].first, side_windows[i].second, coin_win.first, coin_win.second);
            double raw = pr.first; double area = pr.second;
            double norm = (area>0) ? raw * (bg.area_coin / area) : 0.0;
            txt.DrawLatex(0.02, y_text, Form("t2 [%g,%g] raw=%.0f norm=%.2f", side_windows[i].first, side_windows[i].second, raw, norm)); y_text -= 0.025;
        }

        y_text -= 0.01; txt.DrawLatex(0.02, y_text, "Vertical sidebands (t2 in coin):"); y_text -= 0.025;
        for (size_t i=0;i<side_windows.size();++i) {
            auto pr = nps_bg::integral_and_area_TH2(h_t1_t2, coin_win.first, coin_win.second, side_windows[i].first, side_windows[i].second);
            double raw = pr.first; double area = pr.second;
            double norm = (area>0) ? raw * (bg.area_coin / area) : 0.0;
            txt.DrawLatex(0.02, y_text, Form("t1 [%g,%g] raw=%.0f norm=%.2f", side_windows[i].first, side_windows[i].second, raw, norm)); y_text -= 0.025;
        }

        // full boxes
        auto pr_full1 = nps_bg::integral_and_area_TH2(h_t1_t2, full1_t2.first, full1_t2.second, full1_t1.first, full1_t1.second);
        auto pr_full2 = nps_bg::integral_and_area_TH2(h_t1_t2, full2_t1.first, full2_t1.second, full2_t2.first, full2_t2.second);
        double raw_f1 = pr_full1.first, area_f1 = pr_full1.second;
        double raw_f2 = pr_full2.first, area_f2 = pr_full2.second;
        double norm_f1 = (area_f1>0) ? raw_f1 * (bg.area_coin / area_f1) : 0.0;
        double norm_f2 = (area_f2>0) ? raw_f2 * (bg.area_coin / area_f2) : 0.0;
        txt.DrawLatex(0.02, y_text, Form("Full box1 raw=%.0f norm=%.2f", raw_f1, norm_f1)); y_text -= 0.03;
        txt.DrawLatex(0.02, y_text, Form("Full box2 raw=%.0f norm=%.2f", raw_f2, norm_f2)); y_text -= 0.03;

        TLegend *leg = new TLegend(0.65, 0.62, 0.92, 0.88); leg->SetBorderSize(0); leg->SetFillColor(0);
        leg->AddEntry(b_coin, "Coincidence box", "l");
        if (!diag_boxes.empty()) leg->AddEntry(diag_boxes.front(), "Diagonal sidebands", "l");
        if (!hor_boxes.empty()) leg->AddEntry(hor_boxes.front(), "Horizontal sidebands", "l");
        if (!ver_boxes.empty()) leg->AddEntry(ver_boxes.front(), "Vertical sidebands", "l");
        leg->AddEntry(b_full1, "Full accidental boxes", "l");
        leg->Draw("same");

        gPad->Update();
        c_t12->SaveAs(Form("%s/t1t2_run%d.png", outPlotDir.Data(), run));

        // -------------------------
        // Pi0 overlay: coin vs acc
        // -------------------------
        TCanvas *c_pi0 = new TCanvas(name("c_pi0"), "Pi0 inv mass coin vs acc", 900,600);
        h_mpi0_all->SetLineColor(kBlack); h_mpi0_all->SetLineWidth(2);
        h_m_pi0_coin->SetLineColor(kBlue); h_m_pi0_coin->SetLineWidth(2);
        h_m_pi0_acc->SetLineColor(kRed); h_m_pi0_acc->SetLineWidth(2); h_m_pi0_acc->SetLineStyle(2);
        h_mpi0_all->Draw("HIST");
        h_m_pi0_coin->Draw("HIST SAME");
        h_m_pi0_acc->Draw("HIST SAME");
        TLegend *l2 = new TLegend(0.55,0.68,0.88,0.88); l2->SetBorderSize(0); l2->SetFillColor(0);
        l2->AddEntry(h_mpi0_all,"All selected #pi^{0} candidates","l");
        l2->AddEntry(h_m_pi0_coin,"Coincidence [149,151] ns","l");
        l2->AddEntry(h_m_pi0_acc,"Outside coincidence (accidentals)","l");
        l2->Draw();
        c_pi0->SaveAs(Form("%s/pi0_mgg_coin_acc_run%d.png", outPlotDir.Data(), run));

        // -------------------------
        // Compute ROI integrals and multiple accidental estimates
        // -------------------------
        const double roi_low = 0.12, roi_high = 0.14;
        auto hist_integral_and_err = [&](TH1D* h)->pair<double,double> {
            int b_lo = h->FindBin(roi_low), b_hi = h->FindBin(roi_high);
            double err = 0.0;
            double val = h->IntegralAndError(b_lo, b_hi, err);
            return {val, err};
        };

        auto coin_roi = hist_integral_and_err(h_m_pi0_coin);
        auto acc_roi = hist_integral_and_err(h_m_pi0_acc);
        double coin_yield = coin_roi.first, coin_yield_err = coin_roi.second;
        double acc_yield = acc_roi.first, acc_yield_err = acc_roi.second;

        // Method A: scale full accidental boxes' mgg histos using area ratios (full1 & full2 then average)
        // Compute full boxes area integrals from TH2 and scale their mgg ROI integrals
        pair<double,double> full1_mgg = hist_integral_and_err(h_mgg_full1);
        pair<double,double> full2_mgg = hist_integral_and_err(h_mgg_full2);

        double acc_est_full1 = 0.0, acc_est_full1_err=0.0;
        if (area_f1 > 0.0) {
            double scale = bg.area_coin / area_f1;
            acc_est_full1 = full1_mgg.first * scale;
            acc_est_full1_err = full1_mgg.second * scale;
        }
        double acc_est_full2 = 0.0, acc_est_full2_err=0.0;
        if (area_f2 > 0.0) {
            double scale = bg.area_coin / area_f2;
            acc_est_full2 = full2_mgg.first * scale;
            acc_est_full2_err = full2_mgg.second * scale;
        }
        // combine
        double acc_est_full = 0.0, acc_est_full_err = 0.0;
        int nfull = 0;
        if (area_f1>0) { acc_est_full += acc_est_full1; acc_est_full_err += acc_est_full1_err*acc_est_full1_err; ++nfull; }
        if (area_f2>0) { acc_est_full += acc_est_full2; acc_est_full_err += acc_est_full2_err*acc_est_full2_err; ++nfull; }
        if (nfull>0) { acc_est_full /= nfull; acc_est_full_err = sqrt(acc_est_full_err)/nfull; }

        // Method B: average of horizontal and vertical sidebands scaled by area ratios
        double acc_sum_hv = 0.0, acc_err2_hv = 0.0; int nhv=0;
        for (size_t i=0;i<side_windows.size();++i) {
            // horizontal
            auto pr_h = nps_bg::integral_and_area_TH2(h_t1_t2, side_windows[i].first, side_windows[i].second, coin_win.first, coin_win.second);
            double area_h = pr_h.second;
            auto hist_h = h_mgg_hor[i];
            auto valh = hist_integral_and_err(hist_h);
            if (area_h > 0) {
                double est = valh.first * (bg.area_coin / area_h);
                double est_err = valh.second * (bg.area_coin / area_h);
                acc_sum_hv += est; acc_err2_hv += est_err*est_err; ++nhv;
            }
            // vertical
            auto pr_v = nps_bg::integral_and_area_TH2(h_t1_t2, coin_win.first, coin_win.second, side_windows[i].first, side_windows[i].second);
            double area_v = pr_v.second;
            auto hist_v = h_mgg_ver[i];
            auto valv = hist_integral_and_err(hist_v);
            if (area_v > 0) {
                double est = valv.first * (bg.area_coin / area_v);
                double est_err = valv.second * (bg.area_coin / area_v);
                acc_sum_hv += est; acc_err2_hv += est_err*est_err; ++nhv;
            }
        }
        double acc_est_hv_avg=0.0, acc_est_hv_avg_err=0.0;
        if (nhv>0) { acc_est_hv_avg = acc_sum_hv / nhv; acc_est_hv_avg_err = sqrt(acc_err2_hv) / nhv; }

        // Method C: diagonal boxes average
        double acc_sum_diag=0.0, acc_err2_diag=0.0; int ndiag=0;
        for (size_t i=0;i<diag_windows.size();++i) {
            auto pr = nps_bg::integral_and_area_TH2(h_t1_t2, diag_windows[i].first, diag_windows[i].second, diag_windows[i].first, diag_windows[i].second);
            double area = pr.second;
            auto val = hist_integral_and_err(h_mgg_diag[i]);
            if (area>0) {
                double est = val.first * (bg.area_coin / area);
                double est_err = val.second * (bg.area_coin / area);
                acc_sum_diag += est; acc_err2_diag += est_err*est_err; ++ndiag;
            }
        }
        double acc_est_diag = 0.0, acc_est_diag_err=0.0;
        if (ndiag>0) { acc_est_diag = acc_sum_diag/ndiag; acc_est_diag_err = sqrt(acc_err2_diag)/ndiag; }

        // Choose a primary accidental estimate (I keep full-box estimate as primary by default if available)
        double acc_primary = acc_est_full, acc_primary_err = acc_est_full_err;
        string acc_method = "full_boxes_avg";
        if (acc_primary <= 0.0 && acc_est_hv_avg>0) { acc_primary = acc_est_hv_avg; acc_primary_err = acc_est_hv_avg_err; acc_method = "hor_ver_avg"; }
        if (acc_primary <= 0.0 && acc_est_diag>0) { acc_primary = acc_est_diag; acc_primary_err = acc_est_diag_err; acc_method = "diag_avg"; }

        // compute subtracted signal and error
        double signal_sub_full = coin_yield - acc_primary;
        double signal_sub_full_err = sqrt(coin_yield_err*coin_yield_err + acc_primary_err*acc_primary_err);

        // also compute subtracted using hv average specifically
        double signal_sub_hv = coin_yield - acc_est_hv_avg;
        double signal_sub_hv_err = sqrt(coin_yield_err*coin_yield_err + acc_est_hv_avg_err*acc_est_hv_avg_err);

        // -------------------------
        // Print results to console
        // -------------------------
        cout << "===== Run " << run << " summary =====\n";
        cout << " Total events: " << n_total << "\n";
        cout << " HMS-passed: " << n_pass_hms << "\n";
        cout << " Selected for pi0 analysis: " << n_selected_for_analysis << "\n";
        cout << " Coin raw (timing plane): " << bg.n_coin_raw << "\n";
        cout << " Estimated accidentals (bg): " << bg.n_accidentals << " +- " << bg.n_accidentals_err << "\n";
        cout << " Accidental primary method: " << acc_method << "\n";
        cout << Form("ROI [%.3f,%.3f] coin_yield = %.3f +- %.3f", roi_low, roi_high, coin_yield, coin_yield_err) << "\n";
        cout << Form("ROI outside coin acc_yield = %.3f +- %.3f", acc_yield, acc_yield_err) << "\n";
        cout << Form("Acc est (full boxes avg) = %.3f +- %.3f", acc_est_full, acc_est_full_err) << "\n";
        cout << Form("Acc est (hor+ver avg) = %.3f +- %.3f", acc_est_hv_avg, acc_est_hv_avg_err) << "\n";
        cout << Form("Acc est (diag avg) = %.3f +- %.3f", acc_est_diag, acc_est_diag_err) << "\n";
        cout << Form("Signal (coin - acc_primary) = %.3f +- %.3f", signal_sub_full, signal_sub_full_err) << "\n";
        cout << Form("Signal (coin - hv_avg) = %.3f +- %.3f", signal_sub_hv, signal_sub_hv_err) << "\n";
        cout << " Fit results (from fit_pi0_peak): mean=" << fitres.mean << " +- " << fitres.mean_err
             << " sigma=" << fitres.sigma << " +- " << fitres.sigma_err << " signal=" << fitres.gauss_integral << " +- " << fitres.gauss_integral_err << "\n";

        // -------------------------
        // Save per-run ROOT file (all histos + TParameters + canvases)
        // -------------------------
        TString outf = Form("%s/diagnostics_run%d.root", outPlotDir.Data(), run);
        TFile *fout = TFile::Open(outf, "RECREATE");
        if (fout && !fout->IsZombie()) {
            // write histograms
            h_nclusters->Write();
            h_clustE->Write(); h_clustT->Write(); h_clustE_vs_T->Write();
            h_mpi0_all->Write(); h_mpi0_2->Write(); h_mpi0_3->Write(); h_mpi0_4->Write();
            h_mmiss_2->Write(); h_mmiss_3->Write(); h_mmiss_4->Write(); h_mmiss_dvcs->Write();
            h_t1_t2->Write("h_t1_t2", TObject::kOverwrite); h_t1_proj->Write(); h_t2_proj->Write();
            h_m_pi0_coin->Write(); h_m_pi0_acc->Write();
            for (auto *h: h_mgg_diag) h->Write();
            for (auto *h: h_mgg_hor) h->Write();
            for (auto *h: h_mgg_ver) h->Write();
            h_mgg_full1->Write(); h_mgg_full2->Write();
            h_mmiss_vs_mgg->Write();

            // write background and fit scalars (TParameter)
            TParameter<double>("coin_raw", bg.n_coin_raw).Write();
            TParameter<double>("accidental_est", bg.n_accidentals).Write();
            TParameter<double>("accidental_err", bg.n_accidentals_err).Write();
            TParameter<double>("coin_area_ns2", bg.area_coin).Write();

            TParameter<double>("pi0_mean_GeV", fitres.mean).Write();
            TParameter<double>("pi0_mean_err_GeV", fitres.mean_err).Write();
            TParameter<double>("pi0_sigma_GeV", fitres.sigma).Write();
            TParameter<double>("pi0_sigma_err_GeV", fitres.sigma_err).Write();
            TParameter<double>("pi0_signal", fitres.gauss_integral).Write();
            TParameter<double>("pi0_signal_err", fitres.gauss_integral_err).Write();
            TParameter<double>("pi0_bg", fitres.bg_integral).Write();
            TParameter<double>("pi0_bg_err", fitres.bg_integral_err).Write();
            TParameter<double>("pi0_chi2", fitres.chi2).Write();
            TParameter<int>("pi0_ndf", fitres.ndf).Write();

            // write canvases
            c_t12->Write(c_t12->GetName(), TObject::kOverwrite);
            c_pi0->Write(c_pi0->GetName(), TObject::kOverwrite);
            // if fit created a canvas, write it too
            if (!fitres.canvas_name.empty()) {
                TCanvas *cfit = dynamic_cast<TCanvas*>(gROOT->FindObject(fitres.canvas_name.c_str()));
                if (cfit) cfit->Write(fitres.canvas_name.c_str(), TObject::kOverwrite);
            }

            fout->Close();
            logmsg(INFO, Form("Wrote per-run ROOT diagnostics to %s", outf.Data()));
        } else {
            logmsg(WARN, Form("Could not create ROOT output %s", outf.Data()));
        }

        // Save PNGs already done; also write TXT summary and CSV per-run
        TString txtout = Form("%s/summary_run%d.txt", outPlotDir.Data(), run);
        ofstream ftxt(txtout.Data());
        ftxt << "Run " << run << " summary\n";
        ftxt << "Total events: " << n_total << "\n";
        ftxt << "HMS passed: " << n_pass_hms << "\n";
        ftxt << "Selected for pi0 analysis: " << n_selected_for_analysis << "\n";
        ftxt << Form("Coin raw (timing plane): %.1f\n", bg.n_coin_raw);
        ftxt << Form("Estimated accidentals (bg): %.3f +- %.3f\n", bg.n_accidentals, bg.n_accidentals_err);
        ftxt << "Accidental estimates:\n";
        ftxt << Form("  full boxes avg = %.3f +- %.3f\n", acc_est_full, acc_est_full_err);
        ftxt << Form("  hor+ver avg = %.3f +- %.3f\n", acc_est_hv_avg, acc_est_hv_avg_err);
        ftxt << Form("  diag avg = %.3f +- %.3f\n", acc_est_diag, acc_est_diag_err);
        ftxt << Form("ROI [%.3f,%.3f] coin_yield = %.3f +- %.3f\n", roi_low, roi_high, coin_yield, coin_yield_err);
        ftxt << Form("ROI outside coin acc_yield = %.3f +- %.3f\n", acc_yield, acc_yield_err);
        ftxt << Form("Signal (coin - acc_primary=%s) = %.3f +- %.3f\n", acc_method.c_str(), signal_sub_full, signal_sub_full_err);
        ftxt << "Fit results:\n";
        ftxt << Form("  mean = %.6f +- %.6f GeV\n", fitres.mean, fitres.mean_err);
        ftxt << Form("  sigma = %.6f +- %.6f GeV\n", fitres.sigma, fitres.sigma_err);
        ftxt << Form("  signal = %.3f +- %.3f\n", fitres.gauss_integral, fitres.gauss_integral_err);
        ftxt.close();

        // per-run CSV
        TString csv_run = Form("%s/summary_run%d.csv", outPlotDir.Data(), run);
        ofstream fc(csv_run.Data());
        fc << "run,run_current,total_events,pass_hms,n_selected,coin_raw,acc_estimate_t1_t2,acc_est,acc_err,pi0_mean,pi0_mean_err,pi0_sigma,pi0_sigma_err,pi0_signal,pi0_signal_err,pi0_bg,pi0_bg_err,pi0_chi2,pi0_ndf,coin_roi_yield,coin_roi_err,acc_roi_yield,acc_roi_err,acc_est_full,acc_est_full_err,acc_est_hv_avg,acc_est_hv_avg_err,acc_est_diag,acc_est_diag_err,signal_sub_full,signal_sub_full_err,signal_sub_hv,signal_sub_hv_err\n";
        fc << run << "," << run_current << "," << n_total << "," << n_pass_hms << "," << n_selected_for_analysis << "," << bg.n_coin_raw << "," << bg.n_accidentals << ","<< acc_primary << "," << acc_primary_err << ",";
        fc << fitres.mean << "," << fitres.mean_err << "," << fitres.sigma << "," << fitres.sigma_err << "," << fitres.gauss_integral << "," << fitres.gauss_integral_err << "," << fitres.bg_integral << "," << fitres.bg_integral_err << "," << fitres.chi2 << "," << fitres.ndf << ",";
        fc << coin_yield << "," << coin_yield_err << "," << acc_yield << "," << acc_yield_err << ",";
        fc << acc_est_full << "," << acc_est_full_err << "," << acc_est_hv_avg << "," << acc_est_hv_avg_err << "," << acc_est_diag << "," << acc_est_diag_err << ",";
        fc << signal_sub_full << "," << signal_sub_full_err << "," << signal_sub_hv << "," << signal_sub_hv_err << "\n";
        fc.close();

        // append to global CSV
        ofstream fglob(global_csv.Data());
        fglob << "run, run_current,total_events,pass_hms,n_selected,coin_raw,acc_estimate_t1_t2,acc_est,acc_err,pi0_mean,pi0_mean_err,pi0_sigma,pi0_sigma_err,pi0_signal,pi0_signal_err,pi0_bg,pi0_bg_err,pi0_chi2,pi0_ndf,coin_roi_yield,coin_roi_err,acc_roi_yield,acc_roi_err,acc_est_full,acc_est_full_err,acc_est_hv_avg,acc_est_hv_avg_err,acc_est_diag,acc_est_diag_err,signal_sub_full,signal_sub_full_err,signal_sub_hv,signal_sub_hv_err\n";
        // fglob.close();

        // ofstream fglob(global_csv.Data(), ios::app);
        fglob << run << "," << run_current << "," << n_total << "," << n_pass_hms << "," << n_selected_for_analysis << "," << bg.n_coin_raw << "," << bg.n_accidentals << "," << acc_primary << "," << acc_primary_err << ",";
        fglob << fitres.mean << "," << fitres.mean_err << "," << fitres.sigma << "," << fitres.sigma_err << "," << fitres.gauss_integral << "," << fitres.gauss_integral_err << "," << fitres.bg_integral << "," << fitres.bg_integral_err << "," << fitres.chi2 << "," << fitres.ndf << ",";
        fglob << coin_yield << "," << coin_yield_err << "," << acc_yield << "," << acc_yield_err << ",";
        fglob << acc_est_full << "," << acc_est_full_err << "," << acc_est_hv_avg << "," << acc_est_hv_avg_err << "," << acc_est_diag << "," << acc_est_diag_err << ",";
        fglob << signal_sub_full << "," << signal_sub_full_err << "," << signal_sub_hv << "," << signal_sub_hv_err << "\n";
        fglob.close();


        // -------------------------
        // Save extra canvases
        // -------------------------
        // mm vs mgg
        TCanvas *c_mmiss_mgg = new TCanvas(name("c_mmiss_mgg"), "Mmiss vs Mgg", 900,700);
        h_mmiss_vs_mgg->Draw("COLZ");
        c_mmiss_mgg->SaveAs(Form("%s/mmiss_vs_mgg_run%d.png", outPlotDir.Data(), run));

        // cluster diagnostics
        TCanvas *c_cluster = new TCanvas(name("c_cluster"), "Cluster diagnostics", 1000,600);
        c_cluster->Divide(2,1);
        c_cluster->cd(1); h_clustE->Draw();
        c_cluster->cd(2); h_clustT->Draw();
        c_cluster->SaveAs(Form("%s/cluster_E_T_run%d.png", outPlotDir.Data(), run));

        // -------------------------
        // Clean up (delete created objects to avoid collisions)
        // -------------------------
        delete h_nclusters; delete h_clustE; delete h_clustT; delete h_clustE_vs_T;
        delete h_mpi0_all; delete h_mpi0_2; delete h_mpi0_3; delete h_mpi0_4;
        delete h_mmiss_2; delete h_mmiss_3; delete h_mmiss_4; delete h_mmiss_dvcs;
        delete h_t1_t2; delete h_t1_proj; delete h_t2_proj;
        delete h_m_pi0_coin; delete h_m_pi0_acc;
        for (auto *h: h_mgg_diag) delete h;
        for (auto *h: h_mgg_hor) delete h;
        for (auto *h: h_mgg_ver) delete h;
        delete h_mgg_full1; delete h_mgg_full2;
        delete h_mmiss_vs_mgg;

        delete c_t12; delete c_pi0; delete c_mmiss_mgg; delete c_cluster;
        delete b_coin; for (auto *b: diag_boxes) delete b; for (auto *b: hor_boxes) delete b; for (auto *b: ver_boxes) delete b;
        delete b_full1; delete b_full2; delete leg; delete l2;

        f->Close();
        sw_run.Stop();
        logmsg(INFO, Form("Run %d finished. Runtime: %f s (real)", run, sw_run.RealTime()));
    } // end runs

    sw_total.Stop();
    logmsg(INFO, Form("ALL RUNS finished. Total runtime: %f s (real)", sw_total.RealTime()));
    logmsg(INFO, Form("Wrote global summary CSV to %s", (outPlotDir + "/summary_all_runs.csv").Data()));
}

