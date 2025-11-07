// ============================================================================
// Author: Avnish Singh
// Revised: ChatGPT (2025-11-06) - full working ROOT macro using nps_helper.h
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
// Main analysis macro
// ============================================================
void nps_analysis(const TString &skimDir_in="output/skimmed/",
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
    const double DNPS_CM = 407.0;                // NPS z (cm)
    const double MERGE_SPACE2 = 50.0;            // squared (cm^2)
    const double MERGE_TIME_NS = 2.0;            // ns
    const double EBEAM = 10.538;                   // GeV: set to your beam energy
    const double MM_PROTON_SIGMA = 0.05;         // GeV: approximate missing-mass sigma for proton (tune)

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
        Double_t HgtrX=0, HgtrY=0, HgtrTh=0, HgtrPh=0, hdelta=0, HgtrP=0, hreactz=0;
        Double_t HgtrPx=0, HgtrPy=0, HgtrPz=0;
        Double_t hdcx=0, hdcxp=0, hdcy=0, hdcyp=0, hdchit1=0;
        Double_t cointime=0, HhodStatus=0, starTime=0, fptime=0;
        Double_t hbeta=0, hcalepr=0, hcaletot=0, hcernpe=0;
        Double_t hel1=0, helpred=0, helrep=0, helmps=0, helnqrt=0;
        Double_t trig1tdc=0, trig6tdc=0, edtmtdc=0;
        Double_t clusE[10000], clusX[10000], clusY[10000], clusT[10000];
        Double_t vclusE[10000], vclusT[10000];
        Int_t vclusSize[10000];
        Double_t nclust_dbl = 0;
        Double_t block_e[2000], cluster_ID[2000], block_t[2000];
        Int_t vnclus = 0;
        Double_t ctpi1=0, ctpi2=0;
        Double_t hztar=0;

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
        T->SetBranchStatus("H.cer.npeSum", 1); T->SetBranchAddress("H.cer.npeSum", &hcernpe);
        T->SetBranchStatus("H.cal.etotnorm", 1); T->SetBranchAddress("H.cal.etotnorm", &hcaletot);
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
        TH1D *h_nclusters = new TH1D("h_nclusters", "NPS clusters per event;N_{clus};Events", 20, -0.5, 19.5);
        TH1D *h_mpi0 = new TH1D("h_mpi0", "Invariant mass (pi0 candidates);M_{#gamma#gamma} [GeV];Events", 200, 0.0, 0.4);
        TH1D *h_mmiss = new TH1D("h_mmiss", "Missing mass (proton);M_{miss} [GeV];Events", 200, 0.0, 2.0);
        TH1D *h_mmiss_dvcs = new TH1D("h_mmiss_dvcs", "Missing mass (DVCS-like);M_{miss} [GeV];Events", 200, 0.0, 2.0);

        // Counters
        Long64_t n_total = 0, n_pass_hms = 0, n_two_good = 0, n_dvcs_flagged = 0;

        // -------------- event loop --------------
        for (Long64_t ev = 0; ev < nentries; ++ev) {
            if ((ev & 0x3FFF) == 0) { // print progress every ~16384 events
                cout << "  run " << run << " event " << ev << " / " << nentries << "\r" << flush;
            }
            T->GetEntry(ev);
            ++n_total;

            // Basic HMS electron selection
            if (!nps::hms_electron_cuts(edtmtdc, hdelta, HgtrTh, HgtrPh, hcernpe, hcaletot, hreactz)) continue;
            ++n_pass_hms;

            // Guard actual nclust range and cast
            const int nclust = std::max(0, std::min(10000, static_cast<int>(std::lrint(nclust_dbl))));
            h_nclusters->Fill(nclust);

            // merge clusters in-place to remove small splits
            // if (nclust > 1) {
            //     nps::mergeClusters(clusE, clusX, clusY, clusT, nclust, MERGE_SPACE2, MERGE_TIME_NS);
            // }
            const int n_after = nps::packClusters(clusE, clusX, clusY, clusT, nclust);

            // select clusters passing spatial + energy cuts
            vector<int> good_idx;
            good_idx.reserve(4);
            for (int i = 0; i < n_after; ++i) {
                if (nps::nps_spatial_energy_cuts(clusE[i], clusX[i], clusY[i], clusT[i])) good_idx.push_back(i);
            }

            if (good_idx.size() != 2) continue; // only analyze clean 2-cluster candidates
            ++n_two_good;

            // compute invariant mass of two photons
            const double m_pi0 = nps::invariant_mass_pi0(
                clusE[good_idx[0]], clusE[good_idx[1]],
                clusX[good_idx[0]], clusX[good_idx[1]],
                clusY[good_idx[0]], clusY[good_idx[1]],
                DNPS_CM);
            h_mpi0->Fill(m_pi0);

            // compute outgoing electron energy & angle from HMS branches; HgtrP = momentum (GeV)
            const double p_eout = HgtrP;
            const double Ee = std::sqrt(std::max(0.0, p_eout*p_eout + nps::kElectronMass_GeV*nps::kElectronMass_GeV));
            // HgtrTh likely in radians (check tree). If it's radians, convert to degrees:
            const double px_e = HgtrPx;
            const double py_e = HgtrPy;
            const double pz_e = HgtrPz;

            // compute missing mass for proton (π0 final)
            const double mm_proton = nps::missing_mass_proton_pi0(
                EBEAM, Ee, px_e, py_e, pz_e,
                clusE[good_idx[0]], clusE[good_idx[1]],
                clusX[good_idx[0]], clusY[good_idx[0]],
                clusX[good_idx[1]], clusY[good_idx[1]],
                DNPS_CM, -17.51); // theta_nps_deg: tune for your geometry
            

            // If missing mass is far from proton peak (|Δ| > N sigma), check DVCS possibility:
            if (std::fabs(mm_proton - nps::kProtonMass_GeV) > 3.0 * MM_PROTON_SIGMA) {
                // test each photon individually as a DVCS photon -> compute missing mass_dvcs
                const double mm1 = nps::missing_mass_dvcs(EBEAM, Ee, px_e,
                                                         py_e, pz_e,
                                                         clusE[good_idx[0]],
                                                         clusX[good_idx[0]], clusY[good_idx[0]],
                                                         DNPS_CM, -17.51);
                const double mm2 = nps::missing_mass_dvcs(EBEAM, Ee, px_e,
                                                          py_e, pz_e,
                                                         clusE[good_idx[1]],
                                                         clusX[good_idx[1]], clusY[good_idx[1]],
                                                         DNPS_CM, -17.51);

                // If either candidate yields a missing mass near the proton => DVCS-like
                if (std::fabs(mm1 - nps::kProtonMass_GeV) < 3.0 * MM_PROTON_SIGMA ||
                    std::fabs(mm2 - nps::kProtonMass_GeV) < 3.0 * MM_PROTON_SIGMA) {
                    h_mmiss_dvcs->Fill(mm1);
                    h_mmiss_dvcs->Fill(mm2);
                    ++n_dvcs_flagged;
                    continue; // exclude from pi0 analysis (flagged as DVCS-like)
                }
            }

            h_mmiss->Fill(mm_proton);

            // Otherwise: treat as pi0 candidate (we already filled histos above)
        } // end event loop

        cout << endl;
        logmsg(INFO, Form("Run %d summary: total=%lld, pass_hms=%lld, two_good=%lld, dvcs_flagged=%lld",
                         run, n_total, n_pass_hms, n_two_good, n_dvcs_flagged));

        // write histograms to output file
        TString outf = Form("%s/diagnostics_run%d.root", outPlotDir.Data(), run);
        TFile *fout = TFile::Open(outf, "RECREATE");
        if (fout && !fout->IsZombie()) {
            h_nclusters->Write();
            h_mpi0->Write();
            h_mmiss->Write();
            h_mmiss_dvcs->Write();
            fout->Close();
            logmsg(INFO, Form("Wrote diagnostics to %s", outf.Data()));
        } else {
            logmsg(WARN, Form("Could not create output file: %s", outf.Data()));
        }

        // cleanup
        delete h_nclusters;
        delete h_mpi0;
        delete h_mmiss;
        delete h_mmiss_dvcs;

        f->Close();
    } // end runs loop

    sw_total.Stop();
    logmsg(INFO, Form("Total runtime: %f s (real)", sw_total.RealTime()));
}
