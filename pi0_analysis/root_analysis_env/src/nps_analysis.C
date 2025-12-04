// ============================================================================
// File : nps_analysis_full_improved.C
// Author: Avnish Singh (original) + ChatGPT (improvements, 2025-11-12)
// Purpose: Full diagnostic pipeline for NPS pi0 analysis with robust outputs.
// Dependencies: utils.C, nps_helper.h, nps_time_bg.h, nps_comb_bg.h
// ============================================================================

#include "utils.C"
#include "nps_helper.h"
#include "nps_time_bg.h"
#include "nps_comb_bg.h"

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
#include <TF1.h>

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

using namespace std;

// ------------------------------------------------------------------
// Configurable limits & constants
// ------------------------------------------------------------------
constexpr int MAX_CLUS = 20;
constexpr double DEFAULT_TIME_THRESH_NS = 10.0;
constexpr double DEFAULT_TIME_WINDOW_WRT_150 = 10.0;
constexpr double EBEAM_DEFAULT = 10.538;
constexpr int NPRINT_PROGRESS = 16384;

// ------------------------------------------------------------------
// CSV header writer for global summary (overwrite existing file)
// Columns requested by user:
// run, charge_estimate, current_mean_uA, total_entries, pass_hms,
// pass_hms_nps, estimated_accidentals, chi2_ndf_comb_bg,
// pi0_mu_MeV, pi0_sigma_MeV, pi0_signal_counts, mmiss_p_mean_GeV, mmiss_p_sigma_GeV
// ------------------------------------------------------------------
void write_global_csv_header(const TString &path) {
    ofstream f(path.Data(), ios::out);
    f << "run,accumulated_charge(mC),current_mean_uA,total_entries,pass_hms,pass_hms_nps,estimated_accidentals,chi2_ndf_comb_bg,pi0_mu_MeV,pi0_sigma_MeV,pi0_signal_counts,mmiss_p_mean_GeV,mmiss_p_sigma_GeV,run_current_mode_uA\n";
    f.close();
}

// ------------------------------------------------------------------
// Main macro: single-file, full diagnostics
// ------------------------------------------------------------------
void nps_analysis(const TString &skimDir_in="output/skimmed/",
                   const TString &outPlotDir_in="output/plots/x60_4b",
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

    // prepare global CSV (overwrite header)
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

        // open input file & tree
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
        Double_t BCM2_scalerCurrent = 0, BCM2_scalerCharge = 0, H_1MHz_scalerTime = 0;

        // set branch status / addresses (enable only needed branches)
        T->SetBranchStatus("*", 0);
        auto enable = [&](const char* b, void* addr){ if (T->GetBranch(b)) { T->SetBranchStatus(b,1); T->SetBranchAddress(b, addr); } };
        enable("H.gtr.x", &HgtrX);
        enable("H.gtr.y", &HgtrY);
        enable("H.gtr.p", &HgtrP);
        enable("H.gtr.px", &HgtrPx);
        enable("H.gtr.py", &HgtrPy);
        enable("H.gtr.pz", &HgtrPz);
        enable("H.gtr.dp", &hdelta);
        enable("H.gtr.th", &HgtrTh);
        enable("H.gtr.ph", &HgtrPh);
        enable("H.react.z", &hreactz);
        enable("H.cer.npeSum", &hcernpeSum);
        enable("H.cal.etotnorm", &hcaletotnorm);
        enable("T.hms.hEDTM_tdcTimeRaw", &edtmtdc);
        enable("H.BCM2.scalerCurrent", &BCM2_scalerCurrent);
        enable("H.BCM2.scalerCharge", &BCM2_scalerCharge);
        enable("H.1MHz.scalerTime", &H_1MHz_scalerTime);
        enable("NPS.cal.nclust", &nclust_dbl);
        enable("NPS.cal.clusE", &clusE);
        enable("NPS.cal.clusX", &clusX);
        enable("NPS.cal.clusY", &clusY);
        enable("NPS.cal.clusT", &clusT);

        Long64_t nentries = T->GetEntries();
        cout << "Run " << run << " entries: " << nentries << endl;

        // quick-run current estimate (mode & mean) and charge-like diagnostic
        double run_current_mode = 0.0;
        double run_current_mean = 0.0;
        double charge_estimate_uA_counts = 0.0; // diagnostic sum of per-event current readings
        {
            if (T->GetBranch("H.BCM2.scalerCurrent")) {
                TH1D *hCurrent = new TH1D(Form("hCurrent_run%d",run), "Beam Current;I (uA);Counts", 200, 0, 100);
                double sum = 0.0; long cnt = 0;
                for (Long64_t i = 0; i < nentries; ++i) {
                    T->GetEntry(i);
                    if (BCM2_scalerCurrent > 0) {
                        hCurrent->Fill(BCM2_scalerCurrent);
                        sum += BCM2_scalerCurrent;
                        ++cnt;
                        charge_estimate_uA_counts += BCM2_scalerCurrent;
                    }
                }
                if (hCurrent->GetEntries() > 0) {
                    run_current_mode = hCurrent->GetXaxis()->GetBinCenter(hCurrent->GetMaximumBin());
                    run_current_mean = (cnt>0) ? (sum / double(cnt)) : 0.0;
                } else { run_current_mode = run_current_mean = 0.0; }
                delete hCurrent;
            } else {
                run_current_mode = run_current_mean = 0.0;
                charge_estimate_uA_counts = 0.0;
            }
        }

        // double accumulated_charge_mC = nps::get_accumulated_charge(
        //     T,
        //     nentries,
        //     &H_BCM4A_scalerCurrent,
        //     &H_1MHz_scalerTime,
        //     /*min_current=*/ 2.0,  // µA
        //     run,
        //     true
        // );

        // -------------------------
        // Histograms (unique names per run)
        // -------------------------
        auto name = [&](const char* base){ return TString::Format("%s_run%d", base, run); };

        TH1D *h_nclusters = new TH1D(name("h_nclusters"), "NPS clusters per event;N_{clus};Events", 21, -0.5, 20.5);
        TH1D *h_clustE = new TH1D(name("h_clustE"), "Cluster energy;E_{clus} [GeV];Counts", 200, 0.0, 4.0);
        TH1D *h_clustT = new TH1D(name("h_clustT"), "Cluster time; t [ns];Counts", 200, 120, 180);
        TH2D *h_clustE_vs_T = new TH2D(name("h_clustE_vs_T"), "Cluster E vs t; t [ns];E [GeV]", 200, 120, 180, 200, 0, 8);
        TH2D *h_clustXY = new TH2D(name("h_clustXY"), "Cluster X vs Y; X [cm]; Y [cm]", 30, -30, 30, 36, -36, 36);
        TH1D *h_clustE_sum = new TH1D(name("h_clustE_sum"), "Sum cluster E per event;E_{sum} [GeV];Events", 200, 0, 9);
        TH1D *h_opening_angle = new TH1D(name("h_opening_angle"), "Opening angle between two photons;#theta_{open} [rad];Counts", 200, 0, 0.5);
        TH1D *h_photon_Eratio = new TH1D(name("h_photon_Eratio"), "Photon energy ratio E1/E2;E1/E2;Counts", 100, 0, 5);

        TH1D *h_mpi0_all = new TH1D(name("h_mpi0_all"), "Invariant mass (all best-pairs);M_{#gamma#gamma} [GeV];Events", 200, 0.0, 0.4);
        TH1D *h_mpi0_2 = new TH1D(name("h_mpi0_2"), "Invariant mass (2-cluster);M [GeV];Events", 200, 0.0, 0.4);
        TH1D *h_mpi0_3 = new TH1D(name("h_mpi0_3"), "Invariant mass (3-cluster best pair);M [GeV];Events", 200, 0.0, 0.4);
        TH1D *h_mpi0_4 = new TH1D(name("h_mpi0_4"), "Invariant mass (4-cluster best pair);M [GeV];Events", 200, 0.0, 0.4);

        TH1D *h_mmiss_2 = new TH1D(name("h_mmiss_2"), "Missing mass (2-cluster);M_{miss} [GeV];Events", 200, 0.0, 2.0);
        TH1D *h_mmiss_3 = new TH1D(name("h_mmiss_3"), "Missing mass (3-cluster);M_{miss} [GeV];Events", 200, 0.0, 2.0);
        TH1D *h_mmiss_4 = new TH1D(name("h_mmiss_4"), "Missing mass (4-cluster);M_{miss} [GeV];Events", 200, 0.0, 2.0);
        TH1D *h_mmiss_dvcs = new TH1D(name("h_mmiss_dvcs"), "Missing mass (DVCS-like single photon);M_{miss} [GeV];Events", 200, 0.0, 2.0);

        const double t_min = 140.0, t_max = 160.0;
        const int nbins_t = 200;
        TH2D *h_t1_t2 = new TH2D(name("h_t1_t2"), "t1 (y) vs t2 (x);t2 [ns];t1 [ns]", nbins_t, t_min, t_max, nbins_t, t_min, t_max);
        TH1D *h_t1_proj = new TH1D(name("h_t1_proj"), "t1 projection; t1 [ns];Entries", nbins_t, t_min, t_max);
        TH1D *h_t2_proj = new TH1D(name("h_t2_proj"), "t2 projection; t2 [ns];Entries", nbins_t, t_min, t_max);

        // coin / acc histos and per-box mgg histograms
        TH1D *h_m_pi0_coin = new TH1D(name("h_m_pi0_coin"), "Pi0 mass (Coincidence);M [GeV];Counts", 200, 0.0, 0.4);
        TH1D *h_m_pi0_acc  = new TH1D(name("h_m_pi0_acc"),  "Pi0 mass (Accidentals - outside coin);M [GeV];Counts", 200, 0.0, 0.4);
        TH1D *h_coin_bgsub = nullptr;

        // Per-window mgg histograms (diag/side/full)
        vector<pair<double,double>> diag_windows = nps::default_diag_windows();
        vector<pair<double,double>> side_windows = nps::default_side_windows();
        auto coin_win = nps::default_coin_window();
        auto full1_t1 = nps::default_full_acc1_t1();
        auto full1_t2 = nps::default_full_acc1_t2();
        auto full2_t1 = nps::default_full_acc2_t1();
        auto full2_t2 = nps::default_full_acc2_t2();

        vector<TH1D*> h_mgg_diag;
        for (size_t i=0;i<diag_windows.size();++i) {
            h_mgg_diag.push_back(new TH1D(name(TString::Format("h_mgg_diag%d",(int)i)).Data(),
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

        // mm vs t1/t2
        TH2D *h_mmiss_vs_t1 = new TH2D(name("h_mmiss_vs_t1"), "M_{miss} vs t1; t1 [ns]; M_{miss} [GeV]", nbins_t, t_min, t_max, 200, 0.0, 2.0);
        TH2D *h_mmiss_vs_t2 = new TH2D(name("h_mmiss_vs_t2"), "M_{miss} vs t2; t2 [ns]; M_{miss} [GeV]", nbins_t, t_min, t_max, 200, 0.0, 2.0);

        // bookkeeping counters
        Long64_t n_total = 0, n_pass_hms = 0, n_ge2_hms = 0;
        Long64_t n_mult2_hms = 0, n_mult3_hms = 0, n_mult4_hms = 0;
        Long64_t n_mult2_hms_nps = 0, n_mult3_hms_nps = 0, n_mult4_hms_nps = 0;
        Long64_t n_dvcs_flagged = 0, n_selected_for_analysis = 0;

        // -------------------------
        // Event loop
        // -------------------------
        const Long64_t print_every = std::max( (Long64_t)1000, nentries/100 ); // ~100 prints or at least every 1k
        for (Long64_t ev=0; ev<nentries; ++ev) {
            if ((ev % print_every) == 0) {
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
            double Esum = 0;
            for (int i=0;i<nclust;++i) {
                h_clustE->Fill(clusE[i]);
                h_clustT->Fill(clusT[i]);
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

            if (good_idx.size() < 2) continue;
            if (good_idx.size() > 4) {
                sort(good_idx.begin(), good_idx.end(), [&](int a,int b){ return clusE[a] > clusE[b]; });
                good_idx.resize(4);
            }

            if (good_idx.size() == 2) ++n_mult2_hms_nps;
            else if (good_idx.size() == 3) ++n_mult3_hms_nps;
            else if (good_idx.size() == 4) ++n_mult4_hms_nps;

            // // DVCS veto (single-photon events) - placeholder (user-defined)
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
            // compute opening angle from photon unit vectors
            {
                // photon vector from cluster (x,y,z_nps)
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
            const double mm_p = nps::missing_mass_proton_pi0(Ebeam, Ee, px_e, py_e, pz_e,
                                                            clusE[sel_i], clusE[sel_j],
                                                            clusX[sel_i], clusY[sel_i],
                                                            clusX[sel_j], clusY[sel_j],
                                                            nps::kDefaultZ_NPS_cm, -17.51);
            if (good_idx.size()==2) h_mmiss_2->Fill(mm_p);
            else if (good_idx.size()==3) h_mmiss_3->Fill(mm_p);
            else if (good_idx.size()==4) h_mmiss_4->Fill(mm_p);

            h_mmiss_vs_mgg->Fill(mgg, mm_p);
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
        } // end event loop

        cout << endl;

        // -------------------------
        // Summaries & background estimate
        // -------------------------
        nps::CoincidenceBGResult bg = nps::estimate_coincidence_background_default(h_t1_t2);

        // -------------------------
        // Data-driven accidental subtraction (returns bg-subtracted histogram)
        // -------------------------
        h_coin_bgsub = nps::make_and_subtract_accidentals_data_driven(
            h_m_pi0_coin, bg,
            h_mgg_full1, h_mgg_full2,
            h_mgg_hor, h_mgg_ver, h_mgg_diag,
            h_m_pi0_acc,
            h_t1_t2,
            diag_windows, side_windows,
            outPlotDir, run
        );

        // -------------------------
        // Fit combinatorial BG and subtract (user helper)
        // -------------------------
        TH1D* h_final = nullptr;
        nps::BGSubtractionResult res; // default constructed
        // ensure default zeros
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

            // assign to outer-scope pointer (do NOT redeclare)
            h_final = res.h_final;
        }

        // -------------------------
        // Open per-run ROOT file (and ensure we write everything there)
        // -------------------------
        TString outf = Form("%s/diagnostics_run%d.root", outPlotDir.Data(), run);
        TFile *fout = TFile::Open(outf, "RECREATE");
        if (!fout || fout->IsZombie()) {
            logmsg(WARN, Form("Could not create ROOT output %s", outf.Data()));
        } else {
            fout->cd();

            // write histograms (explicit list)
            h_nclusters->Write();
            h_clustE->Write(); h_clustT->Write(); h_clustE_vs_T->Write();
            h_clustXY->Write(); h_clustE_sum->Write(); h_opening_angle->Write(); h_photon_Eratio->Write();
            h_mpi0_all->Write(); h_mpi0_2->Write(); h_mpi0_3->Write(); h_mpi0_4->Write();
            h_mmiss_2->Write(); h_mmiss_3->Write(); h_mmiss_4->Write(); h_mmiss_dvcs->Write();
            h_t1_t2->Write("h_t1_t2", TObject::kOverwrite); h_t1_proj->Write(); h_t2_proj->Write();
            h_m_pi0_coin->Write(); h_m_pi0_acc->Write();
            if (h_coin_bgsub) h_coin_bgsub->Write("h_pi0_coin_bgsub", TObject::kOverwrite);
            if (h_final) h_final->Write("h_pi0_final", TObject::kOverwrite);

            for (auto *h: h_mgg_diag) if (h) h->Write();
            for (auto *h: h_mgg_hor) if (h) h->Write();
            for (auto *h: h_mgg_ver) if (h) h->Write();
            h_mgg_full1->Write(); h_mgg_full2->Write();
            h_mmiss_vs_mgg->Write();
            h_mmiss_vs_t1->Write(); h_mmiss_vs_t2->Write();

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
            logmsg(INFO, Form("Wrote per-run ROOT diagnostics to %s", outf.Data()));
        }

        // -------------------------
        // Draw and save canvases (PNG) for quick look
        // -------------------------
        gStyle->SetOptStat(0);

        // 1) t1 vs t2 canvas with boxes
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
        TBox *b_full2 = new TBox(full2_t2.first, full2_t1.first, full2_t2.second, full2_t1.second); b_full2->SetLineColor(kOrange+7); b_full2->SetLineWidth(2); b_full2->SetLineStyle(4); b_full2->SetFillStyle(0); b_full2->Draw("same");

        TLatex txt; txt.SetNDC(); txt.SetTextSize(0.025);
        txt.DrawLatex(0.02, 0.96, Form("Run %d  Coin raw = %.0f  Acc est = %.3f #pm %.3f", run, bg.n_coin_raw, bg.n_accidentals, bg.n_accidentals_err));
        txt.DrawLatex(0.02, 0.92, Form("Coin area (ns^{2}) = %.3f", bg.area_coin));
        c_t12->SaveAs(Form("%s/t1t2_run%d.png", outPlotDir.Data(), run));

        // 2) pi0 mass overlay coin/acc and bgsub (+ final if exists)
        TCanvas *c_pi0 = new TCanvas(name("c_pi0"), "Pi0 inv mass coin vs acc", 900,600);
        h_mpi0_all->SetLineColor(kBlack); h_mpi0_all->SetLineWidth(1);
        h_m_pi0_coin->SetLineColor(kBlue); h_m_pi0_coin->SetLineWidth(2);
        // h_m_pi0_acc->SetLineColor(kRed); h_m_pi0_acc->SetLineWidth(2); h_m_pi0_acc->SetLineStyle(2);
        h_mpi0_all->Draw("HIST");
        h_m_pi0_coin->Draw("HIST SAME");
        // h_m_pi0_acc->Draw("HIST SAME");
        if (h_coin_bgsub) { h_coin_bgsub->SetLineColor(kGreen+2); h_coin_bgsub->SetLineWidth(2); h_coin_bgsub->Draw("HIST SAME"); }
        if (h_final)      { h_final->SetLineColor(kMagenta+1); h_final->SetLineWidth(2); h_final->Draw("HIST SAME"); }

        TLegend *l2 = new TLegend(0.45,0.60,0.78,0.88); l2->SetBorderSize(0); l2->SetFillColor(0);
        // l2->AddEntry(h_mpi0_all,"All selected #pi^{0} candidates","l");
        // l2->AddEntry(h_m_pi0_coin,"t1 & t2 within coincidence window","l");
        // // l2->AddEntry(h_m_pi0_acc,"t1 & t2 outside coincidence window","l");
        // if (h_coin_bgsub) l2->AddEntry(h_coin_bgsub,"Coincidence - accidental (bg-sub)","l");
        // if (h_final)      l2->AddEntry(h_final,"Final background-subtracted (fit-sub)","l");
        // l2->Draw();
        // Legend shows analysis flow top→bottom
        l2->SetHeader("Analysis flow");
        l2->SetTextSize(0.030);

        l2->AddEntry(h_mpi0_all,   "1) All #pi^{0} candidates",                     "l");
        l2->AddEntry(h_m_pi0_coin, "2) Within coincidence window (t1 & t2)",      "l");
        // l2->AddEntry(h_m_pi0_acc, "Outside coincidence window (accidentals)",   "l");
        if (h_coin_bgsub) l2->AddEntry(h_coin_bgsub, "3) After timing (accidental) subtraction", "l");
        if (h_final)      l2->AddEntry(h_final,      "4) Final after combinatorial (fit) subtraction", "l");

        l2->Draw();


        c_pi0->SaveAs(Form("%s/pi0_mgg_inside_outside_coin_compare_run%d.png", outPlotDir.Data(), run));

        // 3) mm vs mgg
        TCanvas *c_mmiss_mgg = new TCanvas(name("c_mmiss_mgg"), "Mmiss vs Mgg", 900,700);
        h_mmiss_vs_mgg->Draw("COLZ");
        c_mmiss_mgg->SaveAs(Form("%s/mmiss_vs_mgg_run%d.png", outPlotDir.Data(), run));

        // 4) cluster diagnostics
        TCanvas *c_cluster = new TCanvas(name("c_cluster"), "Cluster diagnostics", 1000,800);
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
            // fallback: use combined mmiss_vs_mgg projection if any entries
            if (h_mmiss_vs_mgg->GetEntries() > 0) {
                TH1D *hproj = h_mmiss_vs_mgg->ProjectionY("hproj_mmiss_tmp");
                mmiss_p_mean = hproj->GetMean();
                mmiss_p_sigma = hproj->GetRMS();
                delete hproj;
            } else {
                mmiss_p_mean = mmiss_p_sigma = 0.0;
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
                ftxt << Form("Run current (mode) = %.3f uA   mean = %.3f uA\n", run_current_mode, run_current_mean);
                ftxt.close();
            } else {
                logmsg(WARN, Form("Could not open TXT summary %s for writing", Form("%s/summary_run%d.txt", outPlotDir.Data(), run)));
            }
        }

        // -------------------------
        // Append to global CSV (requested columns)
        // run, charge_estimate_uA_counts, current_mean_uA, total_entries,
        // pass_hms, pass_hms_nps, estimated_accidentals, chi2/ndf for comb bg,
        // final pi0 mu (MeV), final pi0 sigma (MeV), final signal counts,
        // mmiss_p_mean (GeV), mmiss_p_sigma (GeV), run_current_mode_uA
        // -------------------------
        {
            ofstream fg(global_csv.Data(), ios::app);
            if (fg.is_open()) {
                fg << run << ","                             // run
                //    << std::setprecision(6) << accumulated_charge_mC << ","   // accumulated charge
                   << std::setprecision(6) << run_current_mean << ","           // mean current uA
                   << n_total << ","                                           // total entries
                   << n_pass_hms << ","                                        // pass hms
                   << n_selected_for_analysis << ","                          // pass hms + nps
                   << std::setprecision(6) << bg.n_accidentals << ","         // estimated accidentals
                   << std::setprecision(6) << res.chi2_ndf << ","             // chi2/ndf (comb bg)
                   << std::setprecision(6) << res.mu_MeV << ","               // pi0 mu MeV
                   << std::setprecision(6) << res.sigma_MeV << ","            // pi0 sigma MeV
                   << std::setprecision(3) << res.signal_counts << ","        // final signal counts
                   << std::setprecision(6) << mmiss_p_mean << ","             // mmiss p mean (GeV)
                   << std::setprecision(6) << mmiss_p_sigma << ","           // mmiss p sigma (GeV)
                   << std::setprecision(6) << run_current_mode << "\n";      // run current mode
                fg.close();
            } else {
                logmsg(WARN, Form("Could not append to global CSV %s", global_csv.Data()));
            }
        }

        // -------------------------
        // Clean up (delete created objects to avoid collisions)
        // -------------------------
        delete h_nclusters; delete h_clustE; delete h_clustT; delete h_clustE_vs_T;
        delete h_clustXY; delete h_clustE_sum; delete h_opening_angle; delete h_photon_Eratio;
        delete h_mpi0_all; delete h_mpi0_2; delete h_mpi0_3; delete h_mpi0_4;
        delete h_mmiss_2; delete h_mmiss_3; delete h_mmiss_4; delete h_mmiss_dvcs;
        delete h_t1_t2; delete h_t1_proj; delete h_t2_proj;
        delete h_m_pi0_coin; delete h_m_pi0_acc;
        if (h_coin_bgsub) delete h_coin_bgsub;
        if (h_final) delete h_final;
        for (auto *h: h_mgg_diag) delete h;
        for (auto *h: h_mgg_hor) delete h;
        for (auto *h: h_mgg_ver) delete h;
        delete h_mgg_full1; delete h_mgg_full2;
        delete h_mmiss_vs_mgg;
        delete h_mmiss_vs_t1; delete h_mmiss_vs_t2;

        // delete canvases & legend / boxes
        delete c_t12; delete c_pi0; delete c_mmiss_mgg; delete c_cluster;
        delete l2;
        delete b_coin; delete b_full1; delete b_full2;
        for (auto *b: diag_boxes) delete b;
        for (auto *b: hor_boxes) delete b;
        for (auto *b: ver_boxes) delete b;

        f->Close();
        sw_run.Stop();
        logmsg(INFO, Form("Run %d finished. Runtime: %f s (real)", run, sw_run.RealTime()));
    } // end runs

    sw_total.Stop();
    logmsg(INFO, Form("ALL RUNS finished. Total runtime: %f s (real)", sw_total.RealTime()));
    logmsg(INFO, Form("Wrote global summary CSV to %s", (outPlotDir + "/summary_all_runs.csv").Data()));
}
