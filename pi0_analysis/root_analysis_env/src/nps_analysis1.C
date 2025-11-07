// ============================================================================
// Author: Avnish Singh
// Revised: ChatGPT (2025-11-06) - handles 2/3/4-cluster events, DVCS flagging
// ============================================================================

#include "utils.C"  // user-provided: logmsg(), trim(), etc.

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TROOT.h>

#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <cstring>

#include "nps_helper.h"
#include "nps_time_bg.h"

using namespace std;

// ============================================================
// Runlist reader
// ============================================================
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

// ============================================================
// Helper: choose best pair from indices that produce invariant mass closest to pi0 mass
// returns pair indices (i,j) relative to good_idx vector; -1,-1 if not found
// pair<int,int> choose_best_pair_closest_pi0(const vector<int> &good_idx,
//                                            double *clusE, double *clusX, double *clusY,
//                                            double z_nps, double target = nps::kPi0Mass_GeV)
// {
//     int best_i=-1, best_j=-1;
//     double best_diff = 1e9;
//     const int n = (int)good_idx.size();
//     for (int a = 0; a < n; ++a) {
//         for (int b = a+1; b < n; ++b) {
//             int ia = good_idx[a], ib = good_idx[b];
//             double m = nps::invariant_mass_pi0(clusE[ia], clusE[ib],
//                                               clusX[ia], clusX[ib],
//                                               clusY[ia], clusY[ib],
//                                               z_nps);
//             double d = std::fabs(m - target);
//             if (d < best_diff) {
//                 best_diff = d;
//                 best_i = ia; best_j = ib;
//             }
//         }
//     }
//     return {best_i, best_j};
// }

// Choose the best photon pair (indices into clus arrays) whose invariant mass is
// closest to the target π0 mass AND whose cluster times differ by <= time_thresh_ns.
// - good_idx: vector of candidate cluster indices (indices into clusE/clusX/clusY/clusT)
// - clusE, clusX, clusY, clusT: arrays of cluster properties
// - z_nps: distance to NPS (cm), used by invariant_mass_pi0
// - target: target invariant mass (GeV), default = nps::kPi0Mass_GeV
// - time_thresh_ns: maximum allowed |t_i - t_j| in ns (default 2.0)
// Returns pair<int,int> = {index_i, index_j} (indices into the clus arrays), or {-1,-1}
pair<int,int> choose_best_pair_closest_pi0(const vector<int> &good_idx,
                                           double *clusE, double *clusX, double *clusY, double *clusT,
                                           double z_nps,
                                           double target = nps::kPi0Mass_GeV,
                                           double time_thresh_ns = 2.0)
{
    int best_i = -1, best_j = -1;
    double best_diff = 1e9;
    double best_totalE = -1.0; // tie-breaker: prefer pair with larger total energy
    const int n = static_cast<int>(good_idx.size());
    if (n < 2) return { -1, -1 };

    for (int a = 0; a < n; ++a) {
        for (int b = a + 1; b < n; ++b) {
            const int ia = good_idx[a];
            const int ib = good_idx[b];

            // time difference requirement: enforce clusters be within time_thresh_ns
            const double dt = std::fabs(clusT[ia] - clusT[ib]);
            if (dt > time_thresh_ns) continue; // skip pairs not within timing window

            // compute invariant mass for this pair
            const double m = nps::invariant_mass_pi0(clusE[ia], clusE[ib],
                                                    clusX[ia], clusX[ib],
                                                    clusY[ia], clusY[ib],
                                                    z_nps);
            const double d = std::fabs(m - target);

            if (d < best_diff) {
                best_diff = d;
                best_i = ia;
                best_j = ib;
                best_totalE = clusE[ia] + clusE[ib];
            } else if (std::fabs(d - best_diff) < 1e-12) {
                // tie-breaker: choose the pair with larger summed energy
                double totalE = clusE[ia] + clusE[ib];
                if (totalE > best_totalE) {
                    best_i = ia;
                    best_j = ib;
                    best_totalE = totalE;
                }
            }
        }
    }

    // If no pair passed the timing requirement, best_i/j remain -1
    return { best_i, best_j };
}


// ============================================================
// Main analysis macro
// ============================================================
void nps_analysis1(const TString &skimDir_in="output/skimmed/",
                  const TString &outPlotDir_in="output/plots/",
                  const TString &runlistFile="config/runlist_x60_4b.txt")
{
    TStopwatch sw_total; sw_total.Start();
    logmsg(INFO, "=========== NPS π0 diagnostic analysis ===========");

    // make sure directories end with '/'
    TString skimDir = skimDir_in.EndsWith("/") ? skimDir_in : skimDir_in + "/";
    TString outPlotDir = outPlotDir_in.EndsWith("/") ? outPlotDir_in : outPlotDir_in + "/";
    gSystem->mkdir(outPlotDir, true);
    gSystem->mkdir("output", true);

    vector<int> runs = readRunList(runlistFile.Data());
    if (runs.empty()) {
        logmsg(ERROR, "No runs found!");
        return;
    }

    // ---------- analysis constants (tune these) ----------
    const double DNPS_CM = nps::kDefaultZ_NPS_cm;
    const double MERGE_SPACE2 = 50.0;          // squared (cm^2)
    const double MERGE_TIME_NS = 2.0;          // ns
    const double EBEAM = 10.538;               // GeV: set to your beam energy
    const double MM_PROTON_SIGMA = 0.05;       // GeV: approximate missing-mass sigma for proton (tune)
    const int NPRINT_PROGRESS = 16384;

    // loop runs
    for (int run : runs) {
        TString infile = Form("%sskim_run%d.root", skimDir.Data(), run);
        if (gSystem->AccessPathName(infile)) {
            logmsg(WARN, Form("Skipping run %d: file not found: %s", run, infile.Data()));
            continue;
        }

        logmsg(INFO, Form("Processing run %d", run));
        TFile *f = TFile::Open(infile, "READ");
        if (!f || f->IsZombie()) { logmsg(ERROR, Form("Error opening file for run %d", run)); if (f) f->Close(); continue; }

        TTree *T = (TTree*)f->Get("T");
        if (!T) { logmsg(ERROR, Form("Tree not found in run %d", run)); f->Close(); continue; }

        // -------------------------
        // Branch variables
        // -------------------------
        Double_t HgtrX=0, HgtrY=0, HgtrTh=0, HgtrPh=0, hdelta=0, HgtrP=0, hreactz=0, hcernpeSum=0, hcaletotnorm=0;
        Double_t HgtrPx=0, HgtrPy=0, HgtrPz=0;
        Double_t edtmtdc=0;
        Double_t clusE[10000], clusX[10000], clusY[10000], clusT[10000];
        Double_t nclust_dbl = 0;

        // --------- set branch status / addresses (only needed branches) ----------
        T->SetBranchStatus("*", 0);

        // HMS branches (use correct names from your tree!)
        T->SetBranchStatus("H.gtr.x", 1); T->SetBranchAddress("H.gtr.x", &HgtrX);
        T->SetBranchStatus("H.gtr.y", 1); T->SetBranchAddress("H.gtr.y", &HgtrY);
        T->SetBranchStatus("H.gtr.p", 1); T->SetBranchAddress("H.gtr.p", &HgtrP);
        T->SetBranchStatus("H.gtr.px", 1); T->SetBranchAddress("H.gtr.px", &HgtrPx);
        T->SetBranchStatus("H.gtr.py", 1); T->SetBranchAddress("H.gtr.py", &HgtrPy);
        T->SetBranchStatus("H.gtr.pz", 1); T->SetBranchAddress("H.gtr.pz", &HgtrPz);
        T->SetBranchStatus("H.gtr.dp", 1); T->SetBranchAddress("H.gtr.dp", &hdelta);
        T->SetBranchStatus("H.gtr.th", 1); T->SetBranchAddress("H.gtr.th", &HgtrTh);
        T->SetBranchStatus("H.gtr.ph", 1); T->SetBranchAddress("H.gtr.ph", &HgtrPh);
        T->SetBranchStatus("H.react.z", 1); T->SetBranchAddress("H.react.z", &hreactz);
        T->SetBranchStatus("H.cer.npeSum", 1); T->SetBranchAddress("H.cer.npeSum", &hcernpeSum); 
        T->SetBranchStatus("H.cal.etotnorm", 1); T->SetBranchAddress("H.cal.etotnorm", &hcaletotnorm); 
        T->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1); T->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);

        // NPS branches
        T->SetBranchStatus("NPS.cal.nclust", 1); T->SetBranchAddress("NPS.cal.nclust", &nclust_dbl);
        T->SetBranchStatus("NPS.cal.clusE", 1); T->SetBranchAddress("NPS.cal.clusE", &clusE);
        T->SetBranchStatus("NPS.cal.clusX", 1); T->SetBranchAddress("NPS.cal.clusX", &clusX);
        T->SetBranchStatus("NPS.cal.clusY", 1); T->SetBranchAddress("NPS.cal.clusY", &clusY);
        T->SetBranchStatus("NPS.cal.clusT", 1); T->SetBranchAddress("NPS.cal.clusT", &clusT);

        Long64_t nentries = T->GetEntries();
        cout << "Run " << run << " entries: " << nentries << endl;

        // -------------- book some diagnostic histograms --------------
        TH1D *h_nclusters = new TH1D("h_nclusters", "NPS clusters per event;N_{clus};Events", 21, -0.5, 20.5);

        // per-multiplicity invariant mass histograms
        TH1D *h_mpi0_2 = new TH1D("h_mpi0_2", "Invariant mass (2-cluster);M_{#gamma#gamma} [GeV];Events", 200, 0.0, 0.4);
        TH1D *h_mpi0_3 = new TH1D("h_mpi0_3", "Invariant mass (3-cluster best-pair);M_{#gamma#gamma} [GeV];Events", 200, 0.0, 0.4);
        TH1D *h_mpi0_4 = new TH1D("h_mpi0_4", "Invariant mass (4-cluster best-pair);M_{#gamma#gamma} [GeV];Events", 200, 0.0, 0.4);

        // per-multiplicity missing-mass histograms
        TH1D *h_mmiss_2 = new TH1D("h_mmiss_2", "Missing mass (2-cluster);M_{miss} [GeV];Events", 200, 0.0, 2.0);
        TH1D *h_mmiss_3 = new TH1D("h_mmiss_3", "Missing mass (3-cluster best-pair);M_{miss} [GeV];Events", 200, 0.0, 2.0);
        TH1D *h_mmiss_4 = new TH1D("h_mmiss_4", "Missing mass (4-cluster best-pair);M_{miss} [GeV];Events", 200, 0.0, 2.0);

        // DVCS histos
        TH1D *h_mmiss_dvcs = new TH1D("h_mmiss_dvcs", "Missing mass (DVCS-like single photon);M_{miss} [GeV];Events", 200, 0.0, 2.0);
        TH1D *h_mpi0_all = new TH1D("h_mpi0_all", "Invariant mass (all selected best-pairs);M_{#gamma#gamma} [GeV];Events", 200, 0.0, 0.4);
        TH1D *h_m_pi0_coin = new TH1D("h_m_pi0_coin", "Invariant Mass #pi^{0} (Coincidence Window);M_{#gamma#gamma} [GeV];Counts", 200, 0.0, 0.4);
        TH1D *h_m_pi0_acc = new TH1D("h_m_pi0_acc", "Invariant Mass #pi^{0} (Outside Coincidence Window);M_{#gamma#gamma} [GeV];Counts", 200, 0.0, 0.4);

        // choose binning & range to cover sidebands and coincidence window; tune as needed
        const double t_min = 140.0; // ns
        const double t_max = 160.0; // ns
        const int nbins_t = 100;     // fine binning for timing

        TH2D *h_t1_t2 = new TH2D("h_t1_t2", "t1 (y) vs t2 (x);t2 [ns];t1 [ns]", nbins_t, t_min, t_max, nbins_t, t_min, t_max);

        // Optionally make 1D projections (filled later), create them now:
        TH1D *h_t1_proj = new TH1D("h_t1_proj", "t1 projection; t1 [ns];Entries", nbins_t, t_min, t_max);
        TH1D *h_t2_proj = new TH1D("h_t2_proj", "t2 projection; t2 [ns];Entries", nbins_t, t_min, t_max);




        // counters & multiplicity bookkeeping
        Long64_t n_total = 0;
        Long64_t n_pass_hms = 0;
        Long64_t n_ge2_hms = 0;
        Long64_t n_mult2_hms = 0, n_mult3_hms = 0, n_mult4_hms = 0;
        Long64_t n_ge2_hms_nps = 0;
        Long64_t n_mult2_hms_nps = 0, n_mult3_hms_nps = 0, n_mult4_hms_nps = 0;
        Long64_t n_dvcs_flagged = 0;
        Long64_t n_selected_for_analysis = 0;

        // -------------- event loop --------------
        for (Long64_t ev = 0; ev < nentries; ++ev) {
            if ((ev & (NPRINT_PROGRESS-1)) == 0) {
                cout << "  run " << run << " event " << ev << " / " << nentries << "\r" << flush;
            }
            T->GetEntry(ev);
            ++n_total;

            // Basic HMS electron selection
            // we used HgtrX/HgtrY placeholders for npeSum/etot in branch setup; adjust if you use proper variables
            if (!nps::hms_electron_cuts(edtmtdc, hdelta, HgtrTh, HgtrPh, hcernpeSum, hcaletotnorm, hreactz)) continue;
            ++n_pass_hms;

            // Guard nclust range and cast
            const int nclust = std::max(0, std::min(10000, static_cast<int>(std::lrint(nclust_dbl))));
            h_nclusters->Fill(nclust);
            if (nclust < 2) continue;
            // if (nclust > 4) continue;

            ++n_ge2_hms; // counts >=2 cluster events
            if (nclust == 2) ++n_mult2_hms;
            else if (nclust == 3) ++n_mult3_hms;
            else if (nclust == 4) ++n_mult4_hms;

            // merge small splits and pack
            // if (nclust > 1) nps::mergeClusters(clusE, clusX, clusY, clusT, nclust, MERGE_SPACE2, MERGE_TIME_NS);
            const int n_after = nps::packClusters(clusE, clusX, clusY, clusT, nclust);

            // find good clusters passing fiducial and energy cuts
            vector<int> good_idx;
            good_idx.reserve(8);
            for (int i = 0; i < n_after; ++i) {
                if (nps::nps_spatial_energy_cuts(clusE[i], clusX[i], clusY[i], clusT[i], 10.0)) good_idx.push_back(i); // time_diff=2.0 (150ns+-2.0ns); last argument
            }

            // Only consider up to 4 clusters for this extended analysis (you can increase if desired)
            if (good_idx.size() < 2) continue;
            if (good_idx.size() > 4) {
                // keep only leading 4 by energy (to limit combinatorics) - sort indices by energy desc
                sort(good_idx.begin(), good_idx.end(), [&](int a, int b){ return clusE[a] > clusE[b]; });
                good_idx.resize(4);
            }

            // Multiplicity bookkeeping (based on good_idx size)
            if (good_idx.size() == 2) ++n_mult2_hms_nps;
            else if (good_idx.size() == 3) ++n_mult3_hms_nps;
            else if (good_idx.size() == 4) ++n_mult4_hms_nps;

            // DVCS rejection step: test each individual cluster as a single photon candidate.
            // if any single-cluster missing mass (with electron) is near the proton mass within 3*sigma,
            // flag the whole event as DVCS-like and exclude it from pi0 analysis.
            // Need outgoing electron momentum components:
            const double px_e = HgtrPx;
            const double py_e = HgtrPy;
            const double pz_e = HgtrPz;
            const double p_e_mom = std::sqrt(std::max(0.0, px_e*px_e + py_e*py_e + pz_e*pz_e));
            const double Ee = std::sqrt(std::max(0.0, p_e_mom*p_e_mom + nps::kElectronMass_GeV*nps::kElectronMass_GeV));

            bool flagged_dvcs = false;
            for (int idx : good_idx) {
                const double mm_single = nps::missing_mass_dvcs(EBEAM, Ee, px_e, py_e, pz_e,
                                                               clusE[idx], clusX[idx], clusY[idx],
                                                               DNPS_CM, -17.51);
                
                if (std::fabs(mm_single - nps::kProtonMass_GeV) < 3.0 * MM_PROTON_SIGMA) {
                    flagged_dvcs = true;
                    h_mmiss_dvcs->Fill(mm_single);
                    break;
                }
            }
            if (flagged_dvcs) {++n_dvcs_flagged; continue; }  // event given to dvcs missing mass histo only if it lies within 3sigma of proton.
            

            // At this point event passed DVCS veto and will be used for pi0 analysis.
            // Choose best pair (for 2-cluster it's trivial; for 3/4 pick best pair closest to pi0 mass).
            int sel_i = -1, sel_j = -1;
            if (good_idx.size() == 2) {
                sel_i = good_idx[0]; sel_j = good_idx[1];
            } else {
                auto pr = choose_best_pair_closest_pi0(good_idx, clusE, clusX, clusY, clusT, DNPS_CM, nps::kPi0Mass_GeV);
                sel_i = pr.first; sel_j = pr.second;
                if (sel_i < 0 || sel_j < 0) continue; // safety
            }

            // inside the event loop, after you decide sel_i, sel_j (the two clusters used for pi0 analysis)
            // fill the 2D timing histogram using cluster times clusT[sel_i], clusT[sel_j].
            // note: user previously used t1 as y-axis and t2 as x-axis. We'll treat clusT[sel_i] as t1 and clusT[sel_j] as t2
            // you can choose ordering; keep consistent with helper assumptions.

            if (sel_i >= 0 && sel_j >= 0) {
                const double t1 = clusT[sel_i]; // y-axis
                const double t2 = clusT[sel_j]; // x-axis
                // fill 2D hist
                h_t1_t2->Fill(t2, t1);
                // project bins accumulation for quick checks
                h_t1_proj->Fill(t1);
                h_t2_proj->Fill(t2);
            }

            // compute invariant mass for the selected pair and fill appropriate histo
            const double m_pi0_sel = nps::invariant_mass_pi0(clusE[sel_i], clusE[sel_j],
                                                             clusX[sel_i], clusX[sel_j],
                                                             clusY[sel_i], clusY[sel_j],
                                                             DNPS_CM);
            h_mpi0_all->Fill(m_pi0_sel);
            if (good_idx.size() == 2) h_mpi0_2->Fill(m_pi0_sel);
            else if (good_idx.size() == 3) h_mpi0_3->Fill(m_pi0_sel);
            else if (good_idx.size() == 4) h_mpi0_4->Fill(m_pi0_sel);

            // compute missing mass for proton for the selected pair
            const double mm_p = nps::missing_mass_proton_pi0(EBEAM, Ee, px_e, py_e, pz_e,
                                                             clusE[sel_i], clusE[sel_j],
                                                             clusX[sel_i], clusY[sel_i],
                                                             clusX[sel_j], clusY[sel_j],
                                                             DNPS_CM, -17.51);
            if (good_idx.size() == 2) h_mmiss_2->Fill(mm_p);
            else if (good_idx.size() == 3) h_mmiss_3->Fill(mm_p);
            else if (good_idx.size() == 4) h_mmiss_4->Fill(mm_p);

            ++n_selected_for_analysis;

            //---------------------------------------------
            // Define coincidence window for NPS cluster timing
            //---------------------------------------------
            const double coin_lo = 149.0;
            const double coin_hi = 151.0;

            // Use the average of the two selected NPS cluster times as the coincidence time
            const double coin_time = 0.5 * (clusT[sel_i] + clusT[sel_j]);

            //---------------------------------------------
            // Fill invariant-mass histograms for timing comparison
            //---------------------------------------------
            if (coin_time > coin_lo && coin_time < coin_hi) {
                // Inside coincidence window
                h_m_pi0_coin->Fill(m_pi0_sel);
            } else {
                // Outside coincidence window (accidentals)
                h_m_pi0_acc->Fill(m_pi0_sel);
            }

        } // end event loop

        cout << endl;

        // -------------------------
        // Print and save summary statistics
        // -------------------------
        const double denom = (n_ge2_hms > 0) ? static_cast<double>(n_ge2_hms) : 1.0;
        const double pct2 = 100.0 * static_cast<double>(n_mult2_hms) / denom;
        const double pct3 = 100.0 * static_cast<double>(n_mult3_hms) / denom;
        const double pct4 = 100.0 * static_cast<double>(n_mult4_hms) / denom;
        const double pct234 = 100.0 * static_cast<double>(n_mult2_hms + n_mult3_hms + n_mult4_hms) / denom;

        cout << "===== Run " << run << " multiplicity summary (among >=2-cluster events) =====\n";
        cout << "  total events processed: " << n_total << "\n";
        cout << "  events passing HMS selection: " << n_pass_hms << "\n";
        cout << "  ==== cluster multiplicity summary passing the HMS cuts ====  " << "\n";
        cout << "    events with >=2 NPS clusters: " << n_ge2_hms << "\n";
        cout << "    2-cluster events: " << n_mult2_hms << " (" << pct2 << " % of >=2)\n";
        cout << "    3-cluster events: " << n_mult3_hms << " (" << pct3 << " % of >=2)\n";
        cout << "    4-cluster events: " << n_mult4_hms << " (" << pct4 << " % of >=2)\n";
        cout << "    (2+3+4) combined: " << (n_mult2_hms+n_mult3_hms+n_mult4_hms) << " (" << pct234 << " % of >=2)\n";
        cout << "  ============================================================  " << "\n";
        cout << "  ==== cluster multiplicity summary passing the HMS cuts and NPS cuts ====  " << "\n";
        // cout << "    events with >=2 NPS clusters: " << n_ge2 << "\n";
        cout << "    2-cluster events: " << n_mult2_hms_nps << "\n";
        cout << "    3-cluster events: " << n_mult3_hms_nps << "\n";
        cout << "    4-cluster events: " << n_mult4_hms_nps << "\n";
        // cout << "    (2+3+4) combined: " << (n_mult2+n_mult3+n_mult4) << " (" << pct234 << " % of >=2)\n";
        cout << "  ============================================================  " << "\n";
        cout << "  DVCS-flagged events: " << n_dvcs_flagged << "\n";
        cout << "  Events selected for π0 analysis after HMS cuts, NPS cuts, and DVCS veto: " << n_selected_for_analysis << "\n";
        cout << "  ============================================================  " << "\n";

        nps_bg::CoincidenceBGResult bg = nps_bg::estimate_coincidence_background_default(h_t1_t2);
        // Print a concise summary to both cout and logmsg
        // std::cout << "========== Coincidence timing BG estimate (run " << run << ") =============\n";
        // std::cout << bg.summary() << std::endl;

        logmsg(INFO, Form("Run %d: coin_raw=%.1f estimated_accidentals=%.3f +- %.3f",
            run, bg.n_coin_raw, bg.n_accidentals, bg.n_accidentals_err));

        // Fit and print result; draw and save PNG in outPlotDir. Provide run number for naming.
        nps::FitResult fitres = nps::fit_pi0_peak(h_mpi0_all, 0.02, 0.30, true, outPlotDir, run);


        // Optionally draw rectangles on the 2D hist to visualize the boxes (one-time canvas)
        // Save a canvas PNG of the 2D plot with coin box overlay:
        // TCanvas *c_t12 = new TCanvas(Form("c_t12_run%d", run), "t1 vs t2", 800, 700);
        // h_t1_t2->Draw("COLZ");
        // gPad->Update();

        // // draw coincidence rectangle
        // double cx1 = 149.0, cx2 = 151.0, cy1 = 149.0, cy2 = 151.0;
        // TLine *r1 = new TLine(cx1, cy1, cx2, cy1); r1->SetLineColor(kRed); r1->SetLineWidth(2); r1->Draw("same");
        // TLine *r2 = new TLine(cx1, cy2, cx2, cy2); r2->SetLineColor(kRed); r2->SetLineWidth(2); r2->Draw("same");
        // TLine *r3 = new TLine(cx1, cy1, cx1, cy2); r3->SetLineColor(kRed); r3->SetLineWidth(2); r3->Draw("same");
        // TLine *r4 = new TLine(cx2, cy1, cx2, cy2); r4->SetLineColor(kRed); r4->SetLineWidth(2); r4->Draw("same");

        // // add text with numbers
        // TLatex txt; txt.SetNDC(); txt.SetTextSize(0.03);
        // txt.DrawLatex(0.15, 0.92, Form("Coin raw: %.1f, Accidentals est: %.3f #pm %.3f",
        //     bg.n_coin_raw, bg.n_accidentals, bg.n_accidentals_err));

        // c_t12->SaveAs(Form("%s/t1t2_run%d.png", outPlotDir.Data(), run));

        // --- Draw all windows and annotate counts ---
        gStyle->SetOptStat(0);
        TCanvas *c_t12 = new TCanvas(Form("c_t12_run%d", run), "t1 vs t2", 1000, 900);
        h_t1_t2->Draw("COLZ");
        gPad->SetRightMargin(0.15);
        gPad->Update();

        // prepare windows from defaults in nps_helper.h
        auto coin_win = nps_bg::default_coin_window();      // {149,151}
        auto diag_wins = nps_bg::default_diag_windows();    // vector of {141,143},...
        auto side_wins = nps_bg::default_side_windows();    // same as diag by default
        auto full1_t1 = nps_bg::default_full_acc1_t1();     // {153,159}
        auto full1_t2 = nps_bg::default_full_acc1_t2();     // {141,147}
        auto full2_t2 = nps_bg::default_full_acc2_t2();     // {153,159}
        auto full2_t1 = nps_bg::default_full_acc2_t1();     // {141,147}

        // coin box (red)
        double cx_lo = coin_win.first, cx_hi = coin_win.second;
        double cy_lo = coin_win.first, cy_hi = coin_win.second;
        TBox *b_coin = new TBox(cx_lo, cy_lo, cx_hi, cy_hi);
        b_coin->SetLineColor(kRed); b_coin->SetLineWidth(2); b_coin->SetFillStyle(0);
        b_coin->Draw("same");

        // diagonal boxes (magenta)
        std::vector<TBox*> diag_boxes;
        for (auto &w : diag_wins) {
            TBox *b = new TBox(w.first, w.first, w.second, w.second);
            b->SetLineColor(kMagenta); b->SetLineWidth(2); b->SetLineStyle(2); b->SetFillStyle(0);
            b->Draw("same");
            diag_boxes.push_back(b);
        }

        // horizontal sidebands: t1 in coin, t2 in side windows (blue)
        std::vector<TBox*> hor_boxes;
        for (auto &w : side_wins) {
            TBox *b = new TBox(w.first, cy_lo, w.second, cy_hi);
            b->SetLineColor(kBlue); b->SetLineWidth(2); b->SetLineStyle(3); b->SetFillStyle(0);
            b->Draw("same");
            hor_boxes.push_back(b);
        }

        // vertical sidebands: t2 in coin, t1 in side windows (green)
        std::vector<TBox*> ver_boxes;
        for (auto &w : side_wins) {
            TBox *b = new TBox(cx_lo, w.first, cx_hi, w.second);
            b->SetLineColor(kGreen+2); b->SetLineWidth(2); b->SetLineStyle(3); b->SetFillStyle(0);
            b->Draw("same");
            ver_boxes.push_back(b);
        }

        // full accidental boxes (orange / dark yellow)
        TBox *b_full1 = new TBox(full1_t2.first, full1_t1.first, full1_t2.second, full1_t1.second);
        b_full1->SetLineColor(kOrange+1); b_full1->SetLineWidth(2); b_full1->SetLineStyle(4); b_full1->SetFillStyle(0);
        b_full1->Draw("same");

        TBox *b_full2 = new TBox(full2_t1.first, full2_t2.first, full2_t1.second, full2_t2.second);
        b_full2->SetLineColor(kOrange+7); b_full2->SetLineWidth(2); b_full2->SetLineStyle(4); b_full2->SetFillStyle(0);
        b_full2->Draw("same");

        // --- compute & annotate counts for each box using helper integral function ---
        // coin counts already available in bg.n_coin_raw and bg.area_coin
        TLatex txt; txt.SetNDC(); txt.SetTextSize(0.025); txt.SetTextFont(42);

        // print global summary at top-left
        txt.DrawLatex(0.02, 0.96, Form("Run %d  Coin raw = %.0f  Acc est = %.3f #pm %.3f",
                                    run, bg.n_coin_raw, bg.n_accidentals, bg.n_accidentals_err));
        txt.DrawLatex(0.02, 0.92, Form("Coin area (ns^{2}) = %.3f", bg.area_coin));

        // annotate diagonal boxes: print each box raw & normalized (area scaling)
        double y_text = 0.88;
        for (auto &w : diag_wins) {
            auto pr = nps_bg::integral_and_area_TH2(h_t1_t2, w.first, w.second, w.first, w.second);
            double raw = pr.first;
            double area = pr.second;
            double norm = (area>0) ? raw * (bg.area_coin / area) : 0.0;
            txt.DrawLatex(0.02, y_text, Form("Diag [%g,%g] raw=%.0f norm=%.2f", w.first, w.second, raw, norm));
            y_text -= 0.025;
        }

        // horizontal sidebands summary (t1 in coin, t2 in side windows)
        y_text -= 0.01;
        txt.DrawLatex(0.02, y_text, "Horizontal sidebands (t1 in coin):"); y_text -= 0.025;
        for (auto &w : side_wins) {
            auto pr = nps_bg::integral_and_area_TH2(h_t1_t2, w.first, w.second, cy_lo, cy_hi); // x: side, y: coin
            double raw = pr.first; double area = pr.second;
            double norm = (area>0) ? raw * (bg.area_coin / area) : 0.0;
            txt.DrawLatex(0.02, y_text, Form("  t2 [%g,%g]  raw=%.0f  norm=%.2f", w.first, w.second, raw, norm));
            y_text -= 0.025;
        }

        // vertical sidebands summary (t2 in coin, t1 in side windows)
        y_text -= 0.01;
        txt.DrawLatex(0.02, y_text, "Vertical sidebands (t2 in coin):"); y_text -= 0.025;
        for (auto &w : side_wins) {
            auto pr = nps_bg::integral_and_area_TH2(h_t1_t2, cx_lo, cx_hi, w.first, w.second); // x: coin, y: side
            double raw = pr.first; double area = pr.second;
            double norm = (area>0) ? raw * (bg.area_coin / area) : 0.0;
            txt.DrawLatex(0.02, y_text, Form("  t1 [%g,%g]  raw=%.0f  norm=%.2f", w.first, w.second, raw, norm));
            y_text -= 0.025;
        }

        // full accidental boxes summary
        auto pr_full1 = nps_bg::integral_and_area_TH2(h_t1_t2, full1_t2.first, full1_t2.second, full1_t1.first, full1_t1.second);
        auto pr_full2 = nps_bg::integral_and_area_TH2(h_t1_t2, full2_t1.first, full2_t1.second, full2_t2.first, full2_t2.second);
        double raw_f1 = pr_full1.first, area_f1 = pr_full1.second;
        double raw_f2 = pr_full2.first, area_f2 = pr_full2.second;
        double norm_f1 = (area_f1>0) ? raw_f1 * (bg.area_coin / area_f1) : 0.0;
        double norm_f2 = (area_f2>0) ? raw_f2 * (bg.area_coin / area_f2) : 0.0;
        txt.DrawLatex(0.02, y_text, Form("Full_box1 x=[%g,%g] y=[%g,%g] raw=%.0f norm=%.2f",
                                        full1_t2.first, full1_t2.second, full1_t1.first, full1_t1.second, raw_f1, norm_f1));
        y_text -= 0.03;
        txt.DrawLatex(0.02, y_text, Form("Full_box2 x=[%g,%g] y=[%g,%g] raw=%.0f norm=%.2f",
                                        full2_t1.first, full2_t1.second, full2_t2.first, full2_t2.second, raw_f2, norm_f2));
        y_text -= 0.03;

        // add legend
        TLegend *leg = new TLegend(0.65, 0.62, 0.92, 0.88);
        leg->SetBorderSize(0); leg->SetFillColor(0); leg->SetTextSize(0.025);
        leg->AddEntry(b_coin, "Coincidence box", "l");
        leg->AddEntry(diag_boxes.front(), "Diagonal sidebands", "l");
        leg->AddEntry(hor_boxes.front(), "Horizontal sidebands", "l");
        leg->AddEntry(ver_boxes.front(), "Vertical sidebands", "l");
        leg->AddEntry(b_full1, "Full accidental boxes", "l");
        leg->Draw("same");

        // redraw TLatex/legend on top
        gPad->Update();

        // Save canvas (PNG and also write canvas into output ROOT file)
        c_t12->SaveAs(Form("%s/t1t2_run%d.png", outPlotDir.Data(), run));
        
        //---------------------------------------------
        // Save and visualize all π0 invariant mass distributions
        //---------------------------------------------
        TCanvas *c_pi0_coin = new TCanvas("c_pi0_coin", "Pi0 Invariant Mass Distributions", 800, 600);

        // Style settings
        h_mpi0_all->SetLineColor(kBlack);
        h_mpi0_all->SetLineWidth(2);
        h_mpi0_all->SetTitle(";Invariant Mass M_{#gamma#gamma} [GeV];Counts");

        h_m_pi0_coin->SetLineColor(kBlue);
        h_m_pi0_coin->SetLineWidth(2);

        h_m_pi0_acc->SetLineColor(kRed);
        h_m_pi0_acc->SetLineStyle(2);
        h_m_pi0_acc->SetLineWidth(2);

        // Draw histograms
        h_mpi0_all->Draw("HIST");              // draw all first for axis scaling
        h_m_pi0_coin->Draw("HIST SAME");
        h_m_pi0_acc->Draw("HIST SAME");

        // Legend
        TLegend *leg_pi0 = new TLegend(0.55, 0.68, 0.88, 0.88);
        leg_pi0->AddEntry(h_mpi0_all,  "All selected #pi^{0} candidates", "l");
        leg_pi0->AddEntry(h_m_pi0_coin, "Coincidence Window [149,151] ns", "l");
        leg_pi0->AddEntry(h_m_pi0_acc,  "Outside Coincidence Window (Accidentals)", "l");
        leg_pi0->Draw();

        // Save plots
        c_pi0_coin->SaveAs(Form("output/plots/pi0_invmass_all_coin_acc_run%d_10ns.png", run));


        // //---------------------------------------------
        // // Compute and print integral of accidental π0 histogram in ROI [0.12, 0.14] GeV
        // //---------------------------------------------
        // const double roi_low = 0.12;
        // const double roi_high = 0.14;

        // // Find the corresponding bin indices
        // int bin_lo = h_m_pi0_acc->FindBin(roi_low);
        // int bin_hi = h_m_pi0_acc->FindBin(roi_high);

        // // Compute integral and error
        // double acc_int_err = 0.0;
        // double acc_int = h_m_pi0_acc->IntegralAndError(bin_lo, bin_hi, acc_int_err);

        // std::cout << Form("[INFO] Accidental π0 integral (%.3f–%.3f GeV) = %.3f ± %.3f",
        //                 roi_low, roi_high, acc_int, acc_int_err) << std::endl;

        //---------------------------------------------
        // Compute and print integrals for all π0 histograms in [0.12, 0.14] GeV
        //---------------------------------------------
        const double roi_low = 0.12;
        const double roi_high = 0.14;

        auto get_integral = [&](TH1 *hist, const char *label, int color, double ypos) {
            int bin_lo = hist->FindBin(roi_low);
            int bin_hi = hist->FindBin(roi_high);
            double err = 0.0;
            double val = hist->IntegralAndError(bin_lo, bin_hi, err);
            std::cout << Form("[INFO] %s π0 integral (%.3f–%.3f GeV) = %.3f ± %.3f",
                            label, roi_low, roi_high, val, err)
                    << std::endl;

            // Draw annotation on plot (optional)
            TLatex latex;
            latex.SetNDC();
            latex.SetTextColor(color);
            latex.SetTextSize(0.03);
            latex.DrawLatex(0.18, ypos,
                            Form("%s: %.0f ± %.1f", label, val, err));
        };

        get_integral(h_mpi0_all, "All", kBlack, 0.88);
        get_integral(h_m_pi0_coin, "Coinc", kBlue + 1, 0.83);
        get_integral(h_m_pi0_acc, "Accid", kRed + 1, 0.78);


        // -------------------------
        // Save histograms to root file
        // -------------------------
        TString outf = Form("%s/diagnostics_run%d.root", outPlotDir.Data(), run);
        TFile *fout = TFile::Open(outf, "RECREATE");
        if (fout && !fout->IsZombie()) {
            h_nclusters->Write();
            h_mpi0_all->Write();
            h_mpi0_2->Write();
            h_mpi0_3->Write();
            h_mpi0_4->Write();
            h_mmiss_2->Write();
            h_mmiss_3->Write();
            h_mmiss_4->Write();
            h_mmiss_dvcs->Write();
            h_m_pi0_coin->Write();
            h_m_pi0_acc->Write();

            h_t1_t2->Write("h_t1_t2", TObject::kOverwrite);
            h_t1_proj->Write("h_t1_proj", TObject::kOverwrite);
            h_t2_proj->Write("h_t2_proj", TObject::kOverwrite);

            // Save bg numbers in a small TTree or as TParameter objects for easy reading later:
            // use TParameter (ROOT 6) to store scalars:
            // Note: include <TParameter.h> at top of macro if you want this.
            // Example: write estimated accidental as a TParameter<double>
            TParameter<double> p_coin_raw("coin_raw", bg.n_coin_raw);
            TParameter<double> p_acc_est("accidental_est", bg.n_accidentals);
            TParameter<double> p_acc_err("accidental_err", bg.n_accidentals_err);
            p_coin_raw.Write();
            p_acc_est.Write();
            p_acc_err.Write();

            
            // Write fit scalars
            TParameter<double> p_mean("pi0_mean_GeV", fitres.mean);
            TParameter<double> p_mean_err("pi0_mean_err_GeV", fitres.mean_err);
            TParameter<double> p_sigma("pi0_sigma_GeV", fitres.sigma);
            TParameter<double> p_sigma_err("pi0_sigma_err_GeV", fitres.sigma_err);
            TParameter<double> p_signal("pi0_signal", fitres.gauss_integral);
            TParameter<double> p_signal_err("pi0_signal_err", fitres.gauss_integral_err);
            TParameter<double> p_bg("pi0_bg", fitres.bg_integral);
            TParameter<double> p_bg_err("pi0_bg_err", fitres.bg_integral_err);
            TParameter<double> p_chi2("pi0_chi2", fitres.chi2);
            TParameter<int> p_ndf("pi0_ndf", fitres.ndf);
            p_mean.Write("", TObject::kOverwrite);
            p_mean_err.Write("", TObject::kOverwrite);
            p_sigma.Write("", TObject::kOverwrite);
            p_sigma_err.Write("", TObject::kOverwrite);
            p_signal.Write("", TObject::kOverwrite);
            p_signal_err.Write("", TObject::kOverwrite);
            p_bg.Write("", TObject::kOverwrite);
            p_bg_err.Write("", TObject::kOverwrite);
            p_chi2.Write("", TObject::kOverwrite);
            p_ndf.Write("", TObject::kOverwrite);

            // If fit drew a canvas we named it c_pi0_fit_run<run>; write it to the file
            if (!fitres.canvas_name.empty()) {
                TCanvas *c = (TCanvas*)gROOT->FindObject(fitres.canvas_name.c_str());
                if (c) c->Write(fitres.canvas_name.c_str(), TObject::kOverwrite);
            }


            fout->Close();
            logmsg(INFO, Form("Wrote diagnostics to %s", outf.Data()));
        } else {
            logmsg(WARN, Form("Could not create output file: %s", outf.Data()));
        }

        // cleanup
        delete h_nclusters;
        delete h_mpi0_all;
        delete h_mpi0_2;
        delete h_mpi0_3;
        delete h_mpi0_4;
        delete h_mmiss_2;
        delete h_mmiss_3;
        delete h_mmiss_4;
        delete h_mmiss_dvcs;
        delete c_t12;
        delete h_t1_t2;
        delete h_t1_proj;
        delete h_t2_proj;

        f->Close();
    } // end runs loop

    sw_total.Stop();
    logmsg(INFO, Form("Total runtime: %f s (real)", sw_total.RealTime()));
}
