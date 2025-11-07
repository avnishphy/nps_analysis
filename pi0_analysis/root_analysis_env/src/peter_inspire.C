#include "utils.C"  // logmsg(), trim(), etc.

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TROOT.h>

#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <cmath>
#include <utility>
#include <algorithm>
#include <cstdio>

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
// Robust Gaussian peak fiTer
// ============================================================
pair<double,double> fit_peak(TH1D* h, const string &label) {
    if (!h || h->GetEntries() < 10) {
        int ib = h ? h->GetMaximumBin() : 1;
        return {h ? h->GetBinCenter(ib) : 0.0, 0.0};
    }

    double mean = h->GetMean();
    double rms  = h->GetRMS();
    int ibmax   = h->GetMaximumBin();
    double xpeak = h->GetBinCenter(ibmax);
    double binwidth = h->GetBinWidth(1);
    double xmin = h->GetXaxis()->GetXmin();
    double xmax = h->GetXaxis()->GetXmax();

    cout << "\n==== Fit Debug Info for " << label << " ====\n";
    cout << "Entries: " << h->GetEntries() << "\n";
    cout << "Mean: " << mean << "  RMS: " << rms << "\n";
    cout << "xpeak (max bin center): " << xpeak << "  BinWidth: " << binwidth << "\n";
    cout << "X Range: [" << xmin << ", " << xmax << "]\n";

    if (!isfinite(mean) || rms <= 0)
        mean = xpeak, rms = std::max(binwidth, 0.05 * std::fabs(mean == 0 ? 1.0 : mean));

    double fit_xmin = std::max(xmin, mean - 1 * rms);
    double fit_xmax = std::min(xmax, mean + 1 * rms);
    if (fit_xmax - fit_xmin < 5 * binwidth)
        fit_xmin = std::max(xmin, xpeak - 5 * binwidth),
        fit_xmax = std::min(xmax, xpeak + 5 * binwidth);

    TF1 *f = new TF1(Form("gaus_%s_%p", label.c_str(), (void*)h), "gaus", fit_xmin, fit_xmax);
    f->SetParameters(h->GetMaximum(), mean, rms);
    f->SetParLimits(2, 0.2 * binwidth, xmax - xmin);

    int fitStatus = h->Fit(f, "QRN");
    double peak = (fitStatus == 0) ? f->GetParameter(1) : mean;
    double peak_err = (fitStatus == 0) ? f->GetParError(1) : 0.0;

    delete f;
    return {peak, peak_err};
}

// ============================================================
// Neutral pi
// ============================================================
void pi0_analysis_peter(const TString &skimDir_in="output/skimmed/",
                        const TString &outPlotDir_in="output/plots/",
                        const TString &runlistFile="config/runlist_x60_4b.txt") {

    TStopwatch sw_total;
    sw_total.Start();
    logmsg(INFO, "=========== NPS π0 diagnostic analysis (improved) ===========");

    TString skimDir = skimDir_in.EndsWith("/") ? skimDir_in : skimDir_in + "/";
    TString outPlotDir = outPlotDir_in.EndsWith("/") ? outPlotDir_in : outPlotDir_in + "/";
    gSystem->mkdir(outPlotDir, true);
    gSystem->mkdir("output", true);

    vector<int> runs = readRunList(runlistFile.Data());
    if (runs.empty()) {
        logmsg(ERROR, "No runs found!");
        return;
    }

    // Output plain-text files used in original script
    FILE *f2 = fopen("output/skim_pi0.txt", "w");      // skim for pi0
    FILE *f3 = fopen("output/skim_dvcs.txt", "w");     // skim for DVCS-like
    FILE *f6 = fopen("output/topdrawer.plt", "w");     // topdrawer-like plotting script
    FILE *f7 = fopen("output/block_by_block.txt", "w"); // per-block outputs
    if (!f2) { logmsg(WARN, "Cannot open output/skim_pi0.txt for writing"); }
    if (!f3) { logmsg(WARN, "Cannot open output/skim_dvcs.txt for writing"); }
    if (!f6) { logmsg(WARN, "Cannot open output/topdrawer.plt for writing"); }
    if (!f7) { logmsg(WARN, "Cannot open output/block_by_block.txt for writing"); }

    // constants
    const int N_BLOCKS = 1080;
    const double dnps_cm = 407.0;   // NPS distance from target in cm (consistent with clusX/cluY in cm)
    const double pi0_mass = 0.135;  // GeV
    const double merge_space_threshold = 50.0; // squared distance threshold (cm^2), tune as needed
    const double merge_time_threshold = 2.0;   // ns

    // loop runs
    for (int run : runs) {
        TString infile = Form("%sskim_run%d.root", skimDir.Data(), run);
        if (gSystem->AccessPathName(infile)) {
            logmsg(WARN, Form("Skipping run %d: file not found.", run));
            continue;
        }

        logmsg(INFO, Form("Processing run %d", run));
        TFile *f = TFile::Open(infile, "READ");
        if (!f || f->IsZombie()) { logmsg(ERROR, Form("Error opening file for run %d", run)); if(f) f->Close(); continue; }

        TTree *T = (TTree*)f->Get("T");
        if (!T) { logmsg(ERROR, Form("Tree not found in run %d", run)); f->Close(); continue; }

        // --- Branch variables ---
        Double_t HgtrX=0, HgtrTh=0, HgtrPh=0, hdelta=0;
        Double_t hdcx=0, hdcxp=0, hdcy=0, hdcyp=0;
        Double_t evtType=0, HgtrP=0;
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
        // book-keeping arrays (initialize later)
        static Int_t nchist[100], mallhist[100], mbesthist[100];
        static Int_t epihist[100], mvbesthist[100], msbesthist[100];
        static Int_t cthist1[120], cthist6[120], cthist16[120];
        static Int_t cltimehist[120], cltimehist6[100], ctimehist[120];
        static Int_t t16diffh[100], td1h[100], td6h[100], td16h[100];
        static Int_t eclhist[100][20], ctxh[100][20], ctyh[100][20], cteh[100][20];
        static Int_t blktimeh[1080][20], blkmassh[1080][20];
        static Int_t eloblkhist[1080][20], ehiblkhist[1080][20];
        static Int_t mbycolhist[100][30];
        static Int_t clth[100][5];
        static Double_t avmb[1080], sumb[1080], minion[1080], minioner[1080];
        static Double_t avm[30];
        static Double_t counters[10], ethcntr[10];

        // Initialize arrays to zero for this run
        memset(nchist, 0, sizeof(nchist));
        memset(mallhist, 0, sizeof(mallhist));
        memset(mbesthist, 0, sizeof(mbesthist));
        memset(epihist, 0, sizeof(epihist));
        memset(mvbesthist, 0, sizeof(mvbesthist));
        memset(msbesthist, 0, sizeof(msbesthist));
        memset(cthist1, 0, sizeof(cthist1));
        memset(cthist6, 0, sizeof(cthist6));
        memset(cthist16, 0, sizeof(cthist16));
        memset(cltimehist, 0, sizeof(cltimehist));
        memset(cltimehist6, 0, sizeof(cltimehist6));
        memset(ctimehist, 0, sizeof(ctimehist));
        memset(t16diffh, 0, sizeof(t16diffh));
        memset(td1h, 0, sizeof(td1h));
        memset(td6h, 0, sizeof(td6h));
        memset(td16h, 0, sizeof(td16h));
        memset(eclhist, 0, sizeof(eclhist));
        memset(ctxh, 0, sizeof(ctxh));
        memset(ctyh, 0, sizeof(ctyh));
        memset(cteh, 0, sizeof(cteh));
        memset(blktimeh, 0, sizeof(blktimeh));
        memset(blkmassh, 0, sizeof(blkmassh));
        memset(eloblkhist, 0, sizeof(eloblkhist));
        memset(ehiblkhist, 0, sizeof(ehiblkhist));
        memset(mbycolhist, 0, sizeof(mbycolhist));
        memset(clth, 0, sizeof(clth));
        memset(avmb, 0, sizeof(avmb));
        memset(sumb, 0, sizeof(sumb));
        memset(minion, 0, sizeof(minion));
        memset(minioner, 0, sizeof(minioner));
        memset(avm, 0, sizeof(avm));
        memset(counters, 0, sizeof(counters));
        memset(ethcntr, 0, sizeof(ethcntr));

        // branch status and address
        T->SetBranchStatus("*", 0);
        T->SetBranchStatus("H.dc.x_fp", 1); T->SetBranchAddress("H.dc.x_fp",  &hdcx);
        T->SetBranchStatus("H.dc.y_fp", 1); T->SetBranchAddress("H.dc.y_fp",  &hdcy);
        T->SetBranchStatus("H.dc.xp_fp", 1); T->SetBranchAddress("H.dc.xp_fp", &hdcxp);
        T->SetBranchStatus("H.dc.yp_fp", 1); T->SetBranchAddress("H.dc.yp_fp", &hdcyp);
        T->SetBranchStatus("H.gtr.y", 1); T->SetBranchAddress("H.gtr.y", &HgtrX); // original used HgtrX for x
        T->SetBranchStatus("H.gtr.x", 1); T->SetBranchAddress("H.gtr.x", &HgtrX);
        T->SetBranchStatus("H.gtr.p", 1); T->SetBranchAddress("H.gtr.p", &HgtrP);
        T->SetBranchStatus("H.gtr.beta", 1); T->SetBranchAddress("H.gtr.beta", &hbeta);
        T->SetBranchStatus("H.gtr.dp", 1); T->SetBranchAddress("H.gtr.dp", &hdelta);
        T->SetBranchStatus("H.gtr.th", 1); T->SetBranchAddress("H.gtr.th", &HgtrTh);
        T->SetBranchStatus("H.gtr.ph", 1); T->SetBranchAddress("H.gtr.ph", &HgtrPh);
        T->SetBranchStatus("H.react.z", 1); T->SetBranchAddress("H.react.z", &hztar);
        T->SetBranchStatus("H.cal.eprtracknorm", 1); T->SetBranchAddress("H.cal.eprtracknorm", &hcalepr);
        T->SetBranchStatus("H.cal.etottracknorm", 1); T->SetBranchAddress("H.cal.etottracknorm", &hcaletot);
        T->SetBranchStatus("H.cal.etotnorm", 1); T->SetBranchAddress("H.cal.etotnorm", &hcaltot);
        T->SetBranchStatus("H.cer.npeSum", 1); T->SetBranchAddress("H.cer.npeSum", &hcernpe);
        T->SetBranchStatus("H.hod.goodscinhit", 1); T->SetBranchAddress("H.hod.goodscinhit", &HhodStatus);

        T->SetBranchStatus("T.helicity.hel", 1); T->SetBranchAddress("T.helicity.hel",&hel1);
        T->SetBranchStatus("T.helicity.helpred", 1); T->SetBranchAddress("T.helicity.helpred",&helpred);
        T->SetBranchStatus("T.helicity.helrep", 1); T->SetBranchAddress("T.helicity.helrep",&helrep);
        T->SetBranchStatus("T.helicity.mps", 1); T->SetBranchAddress("T.helicity.mps",&helmps);
        T->SetBranchStatus("T.helicity.npqrt", 1); T->SetBranchAddress("T.helicity.nqrt",&helnqrt);

        T->SetBranchStatus("H.hod.starttime", 1); T->SetBranchAddress("H.hod.starttime", &starTime);

        T->SetBranchStatus("CTime.epCoinTime1_ROC1", 1); T->SetBranchAddress("CTime.epCoinTime1_ROC1", &ctpi1);
        T->SetBranchStatus("CTime.epCoinTime2_ROC1", 1); T->SetBranchAddress("CTime.epCoinTime2_ROC1", &ctpi2);
        T->SetBranchStatus("H.hod.fpHitsTime", 1); T->SetBranchAddress("H.hod.fpHitsTime", &fptime);
        T->SetBranchStatus("H.dc.ntrack", 1); T->SetBranchAddress("H.dc.ntrack", &hdchit1);

        T->SetBranchStatus("NPS.cal.nclust", 1); T->SetBranchAddress("NPS.cal.nclust", &nclust_dbl);
        T->SetBranchStatus("NPS.cal.clusE", 1); T->SetBranchAddress("NPS.cal.clusE", &clusE);
        T->SetBranchStatus("NPS.cal.clusX", 1); T->SetBranchAddress("NPS.cal.clusX", &clusX);
        T->SetBranchStatus("NPS.cal.clusY", 1); T->SetBranchAddress("NPS.cal.clusY", &clusY);
        T->SetBranchStatus("NPS.cal.clusT", 1); T->SetBranchAddress("NPS.cal.clusT", &clusT);
        T->SetBranchStatus("NPS.cal.fly.block_clusterID", 1); T->SetBranchAddress("NPS.cal.fly.block_clusterID",&cluster_ID);
        T->SetBranchStatus("NPS.cal.fly.e", 1); T->SetBranchAddress("NPS.cal.fly.e",&block_e);
        T->SetBranchStatus("NPS.cal.fly.goodAdcTdcDiffTime", 1); T->SetBranchAddress("NPS.cal.fly.goodAdcTdcDiffTime",&block_t);
        T->SetBranchStatus("NPS.cal.vtpClusE", 1); T->SetBranchAddress("NPS.cal.vtpClusE",&vclusE);
        T->SetBranchStatus("NPS.cal.vtpClusSize", 1); T->SetBranchAddress("NPS.cal.vtpClusSize",&vclusSize);
        T->SetBranchStatus("NPS.cal.vtpClusTime", 1); T->SetBranchAddress("NPS.cal.vtpClusTime",&vclusT);
        T->SetBranchStatus("NPS.cal.vtpClusY", 1); T->SetBranchAddress("NPS.cal.vtpClusY",&vnclus);
        T->SetBranchStatus("T.hms.npsTRIG1_tdcTimeRaw", 1); T->SetBranchAddress("T.hms.npsTRIG1_tdcTimeRaw",&trig1tdc);
        T->SetBranchStatus("T.hms.npsTRIG6_tdcTimeRaw", 1); T->SetBranchAddress("T.hms.npsTRIG6_tdcTimeRaw",&trig6tdc);
        T->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1); T->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw",&edtmtdc);

        Long64_t nentries = T->GetEntries();
        cout << "Number of entries in the run: " << nentries << endl;

        // --- Define histograms (once per run)
        TH2F* hEblock_low  = new TH2F("hEblock_low",  "Block Energy Occupancy (Low E);Cluster Energy [GeV];Block ID;Counts",
                                    25, 0.0, 0.5, N_BLOCKS, 0.5, N_BLOCKS + 0.5);
        TH2F* hEblock_high = new TH2F("hEblock_high", "Block Energy Occupancy (High E);Cluster Energy [GeV];Block ID;Counts",
                                    20, 0.0, 5.0, N_BLOCKS, 0.5, N_BLOCKS + 0.5);
        TH1F* hNBlocksPerCluster = new TH1F("hNBlocksPerCluster", "Number of Blocks per Cluster;N_{blocks};Counts", 21, -0.5, 20.5);
        TH2F* hClusterXY = new TH2F("hClusterXY", "Cluster Position;X [cm];Y [cm];Counts", 120, -60, 60, 120, -60, 60);
        TH2F* hClusterEvsT = new TH2F("hClusterEvsT", "Cluster Energy vs Time;Time [ns];Energy [GeV];Counts", 200, 0, 400, 100, 0, 5);

        TH1F* h_nclust      = new TH1F("h_nclust",      "Number of Clusters per Event;# Clusters;Counts", 50, -0.5, 49.5);
        TH1F* h_ctime       = new TH1F("h_ctime",       "Coincidence Time (ctpi1);ctpi1 - 50 [ns];Counts", 120, 0, 120);
        TH1F* h_clust_time  = new TH1F("h_clust_time",  "Cluster Time (All);clusT - 105 [ns];Counts", 100, 0, 100);
        TH1F* h_vclust_time = new TH1F("h_vclust_time", "Cluster Time (VTP);clusT - 105 [ns];Counts", 100, 0, 100);

        TH1F *h_pi0_mass = new TH1F("h_pi0_mass", "Invariant Mass of #pi^{0} Candidates; M_{#gamma#gamma} [GeV]; Counts", 150, 0, 0.5);
        TH1F *h_pi0_energy = new TH1F("h_pi0_energy", "Sum of Cluster Energies; E_{#gamma1}+E_{#gamma2} [GeV]; Counts", 100, 0, 10);
        TH1F *h_opening_angle = new TH1F("h_opening_angle", "Opening Angle between Photons; #theta [rad]; Counts", 100, 0, 0.5);
        TH2F *h_mass_vs_angle = new TH2F("h_mass_vs_angle", "M_{#gamma#gamma} vs #theta; #theta [rad]; M [GeV]", 100, 0, 0.5, 150, 0, 0.5);
        TH2F *h_mass_vs_energy = new TH2F("h_mass_vs_energy", "M_{#gamma#gamma} vs E_{tot}; E_{#gamma1}+E_{#gamma2} [GeV]; M [GeV]", 100, 0, 10, 150, 0, 0.5);

        // Extra histograms that might be useful
        TH1F *h_nblocks = new TH1F("h_nblocks", "Number of blocks in event;N_{blocks};Counts", 500, 0, 500);
        TH1F *h_cluster_energy = new TH1F("h_cluster_energy", "Cluster energy distribution;E [GeV];Counts", 200, 0, 10);

        // Event loop
        for (Long64_t evnt = 0; evnt < nentries; ++evnt) {
            if (evnt % 10000 == 0) cout << "Progress: " << evnt << " / " << nentries << " processing..." << endl;
            T->GetEntry(evnt);
            int nclust = static_cast<int>(nclust_dbl);
            // HMS sanity cut
            if (!(edtmtdc < 0.1 && hdelta > -15.0 && hdelta < 15.0)) continue;

            // HMS angular + calorimeter quick print for early events
            if (HgtrTh > -1 && HgtrTh < 1 && HgtrPh > -1 && HgtrPh < 1 && hcaltot > 0.001) {
                if (evnt > 300 && evnt < 320)
                    printf("HMS %lld %3.1f %6.2f %7.4f %7.4f %5.2f %5.2f\n",
                        evnt, evtType, hdelta, HgtrTh, HgtrPh, hcaletot, hcernpe);
            }

            // Merge nearby clusters (spatial + temporal proximity)
            if (nclust > 1) {
                for (int n1 = 0; n1 < nclust - 1; ++n1) {
                    if (clusE[n1] <= 0) continue;
                    for (int n2 = n1 + 1; n2 < nclust; ++n2) {
                        if (clusE[n2] <= 0) continue;
                        double xdiff = clusX[n1] - clusX[n2];
                        double ydiff = clusY[n1] - clusY[n2];
                        double tdiff = fabs(clusT[n1] - clusT[n2]);
                        double space_diff_sq = xdiff * xdiff + ydiff * ydiff;
                        if (space_diff_sq < merge_space_threshold && tdiff < merge_time_threshold) {
                            double E1 = clusE[n1], E2 = clusE[n2];
                            double Etot = E1 + E2;
                            // energy-weighted combination
                            clusX[n1] = (clusX[n1] * E1 + clusX[n2] * E2) / Etot;
                            clusY[n1] = (clusY[n1] * E1 + clusY[n2] * E2) / Etot;
                            clusT[n1] = (clusT[n1] * E1 + clusT[n2] * E2) / Etot;
                            clusE[n1] = Etot;
                            clusE[n2] = 0.0; // mark removed
                            // reassign block IDs from n2 -> n1
                            for (int b = 0; b < N_BLOCKS; ++b) {
                                if (static_cast<int>(cluster_ID[b]) == n2) cluster_ID[b] = n1;
                            }
                        }
                    }
                }
            }

            // Find highest-energy cluster (nhi/ehi) for DVCS-like checks
            int nhi = -1;
            double ehi = 0.0;
            for (int i = 0; i < nclust; ++i) {
                if (clusE[i] > ehi) { ehi = clusE[i]; nhi = i; }
            }

            // Optionally normalize block energies (only if you want)
            const double enorm = 1.35 * 0.93; // keep original value but flag use if needed
            for (int b = 0; b < N_BLOCKS; ++b) {
                if (isfinite(block_e[b])) block_e[b] *= enorm;
            }

            // Block-cluster correlations -> fill 2D occupancy histograms
            int totalBlocksThisEvent = 0;
            for (int b = 0; b < N_BLOCKS; ++b) {
                int icl = static_cast<int>(cluster_ID[b]);
                if (icl < 0 || icl >= nclust) continue;
                double e_block = block_e[b];
                double e_cluster = clusE[icl];
                double t_cluster = clusT[icl];
                if (e_cluster <= 0 || !isfinite(e_block)) continue;
                if (e_block > 0.3 * e_cluster && e_block < 0.8 * e_cluster && t_cluster > 0. && t_cluster < 300.) {
                    if (e_cluster >= 0.1 && e_cluster < 0.5) hEblock_low->Fill(e_cluster, b + 1);
                    if (e_cluster >= 0.25 && e_cluster < 5.0) hEblock_high->Fill(e_cluster, b + 1);
                }
                ++totalBlocksThisEvent;
            }
            h_nblocks->Fill(totalBlocksThisEvent);

            // Cluster-level histograms
            int nBlocksInCluster = 0;
            for (int c = 0; c < nclust; ++c) {
                if (clusE[c] <= 0) continue;
                hClusterXY->Fill(clusX[c], clusY[c]);
                hClusterEvsT->Fill(clusT[c], clusE[c]);
                h_cluster_energy->Fill(clusE[c]);

                // count associated blocks
                nBlocksInCluster = 0;
                for (int b = 0; b < N_BLOCKS; ++b) {
                    if (static_cast<int>(cluster_ID[b]) == c && block_e[b] > 0) nBlocksInCluster++;
                }
                hNBlocksPerCluster->Fill(nBlocksInCluster);
            }

            // Fill ncluster histogram
            if (nclust >= 0 && nclust < 100) h_nclust->Fill(nclust);
            else if (nclust >= 100) h_nclust->Fill(99);

            // Fill coincidence time histogram (ctpi1)
            double ctime_bin = ctpi1 - 50.0;
            if (ctime_bin >= 0 && ctime_bin < 120) h_ctime->Fill(ctime_bin);

            // Cluster times (all)
            for (int c = 0; c < nclust; ++c) {
                if (clusE[c] > 0.3) {
                    double t_bin = clusT[c] - 105.0;
                    if (t_bin >= 0 && t_bin < 100) h_clust_time->Fill(t_bin);
                }
            }

            // VTP cluster times
            for (int c = 0; c < vnclus; ++c) {
                if (vclusE[c] > 400) {
                    double t_bin = vclusT[c] - 105.0;
                    if (t_bin >= 0 && t_bin < 100) h_vclust_time->Fill(t_bin);
                }
            }

            // ---------- π0 reconstruction ----------
            double ebest = 0.0, mbest = 0.0, opening_theta_best = 0.0;
            int n1best = -1, n2best = -1;
            double massmin = 1e6;

            if (nclust > 1) {
                for (int n1 = 0; n1 < nclust - 1; ++n1) {
                    if (clusE[n1] < 0.1) continue;
                    for (int n2 = n1 + 1; n2 < nclust; ++n2) {
                        if (clusE[n2] < 0.1) continue;
                        double tdiff = fabs(clusT[n1] - clusT[n2]);
                        if (tdiff > 3.0) continue; // time_cut
                        double xdiff = clusX[n1] - clusX[n2];
                        double ydiff = clusY[n1] - clusY[n2];
                        double r_sep = sqrt(xdiff * xdiff + ydiff * ydiff); // cm
                        if (r_sep < 0.03) continue; // remove trivial duplicates

                        double opening_theta = r_sep / dnps_cm; // radians (small-angle approx)
                        double sin2_half = pow(sin(opening_theta / 2.0), 2.0);
                        double arg = 4.0 * clusE[n1] * clusE[n2] * sin2_half;
                        double mass = arg > 0.0 ? sqrt(arg) : 0.0;

                        // fill quick integer hist arrays (same binning as your original)
                        int icc = static_cast<int>(mass / pi0_mass * 50.0);
                        icc = std::clamp(icc, 0, 99);
                        mallhist[icc]++;

                        double mdiff = fabs(mass - pi0_mass);
                        if (mdiff < massmin) {
                            massmin = mdiff;
                            ebest = clusE[n1] + clusE[n2];
                            mbest = mass;
                            n1best = n1;
                            n2best = n2;
                            opening_theta_best = opening_theta;
                        }
                    } // n2
                } // n1

                // Fill histograms for best pair
                if (n1best >= 0 && n2best >= 0) {
                    h_pi0_mass->Fill(mbest);
                    h_pi0_energy->Fill(ebest);
                    h_opening_angle->Fill(opening_theta_best);
                    h_mass_vs_angle->Fill(opening_theta_best, mbest);
                    h_mass_vs_energy->Fill(ebest, mbest);
                }
            } // nclust>1

            // Additional counters & hist arrays around best candidate
            if (ebest > 0.0) {
                int icc = static_cast<int>(mbest / 0.135 * 50.0);
                icc = std::clamp(icc, 0, 99);
                mbesthist[icc]++;

                bool t_best_narrow = (clusT[n1best] > 148.0 && clusT[n1best] < 152.0) &&
                                     (clusT[n2best] > 148.0 && clusT[n2best] < 152.0);
                bool e_best_cut = (clusE[n1best] > 0.6 && clusE[n2best] > 0.6);

                if (t_best_narrow) mvbesthist[icc]++;
                if (t_best_narrow && e_best_cut) {
                    msbesthist[icc]++;
                    double xdiff = clusX[n1best] - clusX[n2best];
                    double xav = 0.5 * (clusX[n1best] + clusX[n2best]);
                    if (fabs(xdiff) < 5.0) {
                        int icol = static_cast<int>((xav + 30.0) / 2.0);
                        icol = std::clamp(icol, 0, 30);
                        mbycolhist[icc][icol]++;
                    }
                    int ie = static_cast<int>(ebest * 10.0);
                    ie = std::clamp(ie, 0, 99);
                    if (mbest > 0.11 && mbest < 0.16) epihist[ie]++;

                    // broad time + mass window counters
                    if ((clusT[n1best] > 140.0 && clusT[n1best] < 164.0) &&
                        (clusT[n2best] > 140.0 && clusT[n2best] < 164.0) &&
                        mbest > 0.01 && mbest < 0.16) {

                        ethcntr[0]++;
                        if (clusE[n1best] > 0.6 && clusE[n2best] > 0.6) ethcntr[1]++;
                        if (clusE[n1best] > 0.7 && clusE[n2best] > 0.7) ethcntr[2]++;
                        if (clusE[n1best] > 0.8 && clusE[n2best] > 0.8) ethcntr[3]++;
                        if (clusE[n1best] > 0.9 && clusE[n2best] > 0.9) ethcntr[4]++;
                        if (clusE[n1best] > 1.0 && clusE[n2best] > 1.0) ethcntr[5]++;
                    }
                }
            }

            // Electron in HMS cut and trig logic
            if (hdelta > -18.0 && hdelta < 18.0 && hcaletot > 0.6 && hcernpe > 0.5) {
                counters[2]++;
                int icc = static_cast<int>(ctpi1 - 50.0); icc = std::clamp(icc, 0, 119);

                int ic1 = (n1best >= 0 ? static_cast<int>(clusE[n1best] / 0.05) : 0);
                int ic2 = (n2best >= 0 ? static_cast<int>(clusE[n2best] / 0.05) : 0);
                ic1 = std::clamp(ic1, 0, 99);
                ic2 = std::clamp(ic2, 0, 99);
                int ichi = (n1best >= 0 && n2best >= 0 && clusE[n1best] > clusE[n2best]) ? ic1 : ic2;
                int iclo = (n1best >= 0 && n2best >= 0 && clusE[n1best] > clusE[n2best]) ? ic2 : ic1;

                int j1 = (n1best >= 0 ? static_cast<int>(clusT[n1best] - ctpi1) : 0);
                int j2 = (n2best >= 0 ? static_cast<int>(clusT[n2best] - ctpi1) : 0);
                j1 = std::clamp(j1, 0, 99);
                j2 = std::clamp(j2, 0, 99);

                bool trg1 = (trig1tdc > 31300 && trig1tdc < 35000);
                bool trg6 = (trig6tdc > 31300 && trig6tdc < 35000);

                if (mbest > 0.11 && mbest < 0.16 && nclust > 1) {
                    if (trg6 && !trg1) {
                        eclhist[iclo][0]++; eclhist[ichi][1]++;
                        td1h[j1]++; td16h[j2]++;
                    }
                    if (trg1 && !trg6) {
                        eclhist[iclo][2]++; eclhist[ichi][3]++;
                        td6h[j1]++; td6h[j2]++;
                    }
                }

                if (trg1 && !trg6) cthist1[icc]++;
                if (trg1 && trg6) {
                    cthist16[icc]++;
                    int icc_diff = static_cast<int>((trig1tdc - trig6tdc)/10.0 + 50.0);
                    icc_diff = std::clamp(icc_diff, 0, 99);
                    t16diffh[icc_diff]++;
                }
                if (trg6 && !trg1) cthist6[icc]++;

                if (trg6 && !trg1) {
                    for (int n = 0; n < nclust; ++n) {
                        if (clusE[n] > 0.3) {
                            int tbin = std::clamp(static_cast<int>(clusT[n] - 105.0), 0, 99);
                            cltimehist6[tbin]++;
                        }
                    }
                }

                // cluster times for good pi0 mass
                if (mbest > 0.11 && mbest < 0.16 && nclust > 1) {
                    int t1 = std::clamp(static_cast<int>(clusT[n1best] - 100.0), 0, 99);
                    int t2 = std::clamp(static_cast<int>(clusT[n2best] - 100.0), 0, 99);
                    int tav = std::clamp(static_cast<int>((clusT[n1best] + clusT[n2best]) / 2.0 - 100.0), 0, 99);
                    clth[t1][0]++; clth[t2][1]++; clth[tav][2]++;
                }

                // ctpi2 placeholder; if you have a real ctpi2 fill it from branch
                ctpi2 = 0.0;

                // skim writing for possible pi0 events
                if (ctpi2 > -30.0 && ctpi2 < 30.0) {
                    counters[3]++;
                    auto isValid = [](double v){ return std::isfinite(v); };
                    if (nclust > 1 && mbest > 0.05 &&
                        isValid(clusE[n1best]) && isValid(clusE[n2best]) &&
                        isValid(clusX[n1best]) && isValid(clusX[n2best]) &&
                        isValid(clusY[n1best]) && isValid(clusY[n2best]) &&
                        isValid(clusT[n1best]) && isValid(clusT[n2best])) {

                        counters[4]++;
                        if (f2) {
                            fprintf(f2,
                                "%6.2f %1d %1d %1d %6.2f %7.4f %7.4f %6.2f %7.4f %7.4f %6.2f %6.2f %7.4f %7.4f %7.4f %7.2f %7.2f %7.2f %7.4f %7.2f %7.2f %7.2f %7.3f  %7.3f   \n",
                                ctpi2,(int)hel1+1,(int)helrep+1,(int)helmps+1,
                                hdelta,HgtrTh,HgtrPh,hcernpe,hcalepr,hcaletot,
                                hdcx, hdcy, hdcxp, hdcyp,
                                clusE[n1best],clusX[n1best],clusY[n1best],clusT[n1best],
                                clusE[n2best],clusX[n2best],clusY[n2best],clusT[n2best],
                                mbest,hztar);
                        }

                        // write block info for each cluster
                        for (int j = 0; j < 2; ++j) {
                            int nn = (j == 0) ? n1best : n2best;
                            double etot = clusE[n1best] + clusE[n2best];
                            for (int b = 0; b < N_BLOCKS; ++b) {
                                int icl = static_cast<int>(cluster_ID[b]);
                                if (icl == nn && block_e[b] > 0.0) {
                                    if (f2) fprintf(f2,"%4d %3d %8.4f %8.1f \n", b, icl, block_e[b], block_t[b]);
                                    if (block_e[b] > 0.3 * etot) {
                                        int jj = std::clamp(static_cast<int>(block_t[b] - 142.0), 0, 19);
                                        blktimeh[b][jj]++;
                                        int jj2 = std::clamp(static_cast<int>((mbest - 0.115) * 20.0 / 0.040), 0, 19);
                                        blkmassh[b][jj2]++;
                                    }
                                }
                            }
                        }
                        if (f2) fprintf(f2,"  -1     0.000 \n");
                    }
                }

                // DVCS-like skim
                if (nhi > -1 && ehi > 2.0 && mbest < 0.05) {
                    counters[5]++;
                    if (f3) {
                        fprintf(f3,
                            "%6.2f %1d %1d %1d %6.2f %7.4f %7.4f %6.2f %7.4f %7.4f %6.2f %6.2f %7.4f %7.4f %7.4f %7.2f %7.2f %7.2f %7.3f   \n",
                            ctpi2,(int)hel1+1,(int)helrep+1,(int)helmps+1,
                            hdelta,HgtrTh,HgtrPh,hcernpe,hcalepr,hcaletot,
                            hdcx, hdcy, hdcxp, hdcyp,
                            clusE[nhi],clusX[nhi],clusY[nhi],clusT[nhi],
                            mbest);
                        for (int b = 0; b < N_BLOCKS; ++b) {
                            int icl = static_cast<int>(cluster_ID[b]);
                            if (icl == nhi && block_e[b] > 0.0) fprintf(f3,"%4d %3d %8.4f %8.1f \n", b, icl, block_e[b], block_t[b]);
                        }
                    }
                }

            } // end hdelta/hcaletot/hcernpe cut

        } // end event loop

        // -------------------------
        // End-of-run prints (condensed & structured)
        // -------------------------
        // print summary tables like your original but in cleaner format
        printf("\nRun %d summary (selected hist arrays):\n", run);
        printf("bin nchist mallhist mbesthist mvbesthist msbesthist epihist cth1 cth6 cth16 clt clt6\n");
        for (int icc = 0; icc < 100; ++icc) {
            printf("%2d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d %6d\n",
                   icc,
                   nchist[icc], mallhist[icc], mbesthist[icc],
                   mvbesthist[icc], msbesthist[icc],
                   epihist[icc], cthist1[icc], cthist6[icc],
                   cthist16[icc], cltimehist[icc], cltimehist6[icc]);
        }

        // counters & energy counters
        printf("\nCounters summary:\n");
        for (int i = 0; i < 6; ++i) {
            double denom = (counters[0] != 0.0) ? counters[0] : 1.0;
            printf(" %7.0f  ratio= %7.3f\n", counters[i], counters[i] / denom);
        }
        printf("\nEnergy counters (ethcntr): threshold frac/total\n");
        for (int i = 0; i < 6; ++i) {
            double denom = (ethcntr[0] != 0.0) ? ethcntr[0] : 1.0;
            printf(" %.2f  %7.0f  %7.3f\n", 0.1 * (i + 5), ethcntr[i], ethcntr[i] / denom);
        }

        // compute avm (average mass by column)
        for (int col = 0; col < 30; ++col) {
            avm[col] = 0.0;
            double sum = 0.0;
            for (int j = 40; j < 60; ++j) {
                double mm = 0.135 / 50.0 * (j + 0.5);
                avm[col] += mm * mbycolhist[j][col];
                sum += mbycolhist[j][col];
            }
            if (sum > 0.0) avm[col] /= sum;
        }
        printf("\nColumn averages & counts (columns 1..30):\n");
        for (int i = 0; i < 30; ++i) {
            printf("%3d %.3f ", i, avm[i]);
            for (int j = 43; j < 56; ++j) printf(" %2d", mbycolhist[j][i]);
            printf(" %2d\n", mbycolhist[57][i]);
        }

        // by-block outputs into f7 and compute minion/minioner & blkmass averages
        if (f7) {
            for (int i = 0; i < N_BLOCKS; ++i) {
                fprintf(f7, "%4d ", i);
                for (int j = 0; j < 20; ++j) fprintf(f7, " %4d", eloblkhist[i][j]);
                fprintf(f7, "\n");
                double sum1 = 0., sum2 = 0.;
                for (int j = 0; j < 19; ++j) {
                    double ebin = 0.1 + 0.02 * (j + 0.5);
                    sum1 += eloblkhist[i][j];
                    sum2 += ebin * eloblkhist[i][j];
                }
                minion[i] = 0.0; minioner[i] = 0.0;
                if (sum1 > 20.0) {
                    minion[i] = 1000.0 * sum2 / sum1;
                    minioner[i] = 50.0 / sqrt(sum1);
                }
            }
            // write block hi-energy & timing histos
            for (int i = 0; i < N_BLOCKS; ++i) {
                fprintf(f7, "%4d ", i);
                for (int j = 0; j < 20; ++j) fprintf(f7, " %4d", ehiblkhist[i][j]);
                fprintf(f7, "\n");
            }
            for (int i = 0; i < N_BLOCKS; ++i) {
                fprintf(f7, "%4d ", i);
                for (int j = 0; j < 20; ++j) fprintf(f7, " %4d", blktimeh[i][j]);
                fprintf(f7, "\n");
            }
            for (int i = 0; i < N_BLOCKS; ++i) {
                avmb[i] = 0.; sumb[i] = 0.;
                for (int j = 0; j < 20; ++j) {
                    double mm = 0.115 + 0.040 * (j + 0.5) / 20.0;
                    avmb[i] += mm * blkmassh[i][j];
                    sumb[i] += blkmassh[i][j];
                }
                fprintf(f7, "%4d ", i);
                for (int j = 0; j < 20; ++j) fprintf(f7, " %4d", blkmassh[i][j]);
                fprintf(f7, "\n");
            }
        }

        // topdrawer-like output (mimic original script)
        if (f6) {
            fprintf(f6, "set device postscript \n");
            fprintf(f6, " set window x 1 of 3 y 1 of 2\n");
            fprintf(f6, " title bottom 'M0GG1 (GeV)' \n");
            fprintf(f6, " case         ' XGGX'\n");
            fprintf(f6, " title top 'run %d' \n", run);
            for (int icc = 20; icc < 80; ++icc) fprintf(f6, "%.3f %d\n", 0.135/50.0*(icc+0.5), mbesthist[icc]);
            fprintf(f6, "hist\n");
            for (int icc = 20; icc < 80; ++icc) fprintf(f6, "%.3f %d\n", 0.135/50.0*(icc+0.5), mvbesthist[icc]);
            fprintf(f6, "hist\n");
            for (int icc = 20; icc < 80; ++icc) fprintf(f6, "%.3f %d\n", 0.135/50.0*(icc+0.5), msbesthist[icc]);
            fprintf(f6, "hist\n");
            // more topdrawer content could be added here following original
            fflush(f6);
        }

        // Save root histograms into a run-specific file
        TString outHistFile = Form("output/nps_diagnostics_run%d.root", run);
        TFile *fout = new TFile(outHistFile, "RECREATE");
        hEblock_low->Write();
        hEblock_high->Write();
        hNBlocksPerCluster->Write();
        hClusterXY->Write();
        hClusterEvsT->Write();
        h_nclust->Write();
        h_ctime->Write();
        h_clust_time->Write();
        h_vclust_time->Write();
        h_pi0_mass->Write();
        h_pi0_energy->Write();
        h_opening_angle->Write();
        h_mass_vs_angle->Write();
        h_mass_vs_energy->Write();
        h_nblocks->Write();
        h_cluster_energy->Write();
        fout->Close();

        // visualization example (saving canvases)
        TCanvas* cBlock = new TCanvas("cBlock", "NPS Block Occupancy", 1200, 800);
        cBlock->Divide(1, 2);
        cBlock->cd(1); gPad->SetLogz(); hEblock_low->Draw("COLZ");
        cBlock->cd(2); gPad->SetLogz(); hEblock_high->Draw("COLZ");
        cBlock->SaveAs(Form("output/block_cluster_energy_occupancy_run%d.pdf", run));
        delete cBlock;

        TCanvas* cTime = new TCanvas("cTime", "NPS Timing Diagnostics", 1000, 800);
        cTime->Divide(2,2);
        cTime->cd(1); h_nclust->Draw();
        cTime->cd(2); h_ctime->Draw();
        cTime->cd(3); h_clust_time->Draw();
        cTime->cd(4); h_vclust_time->Draw();
        cTime->SaveAs(Form("output/nps_timing_diagnostics_run%d.pdf", run));
        delete cTime;

        // cleanup hist objects
        delete hEblock_low; delete hEblock_high; delete hNBlocksPerCluster;
        delete hClusterXY; delete hClusterEvsT;
        delete h_nclust; delete h_ctime; delete h_clust_time; delete h_vclust_time;
        delete h_pi0_mass; delete h_pi0_energy; delete h_opening_angle;
        delete h_mass_vs_angle; delete h_mass_vs_energy;
        delete h_nblocks; delete h_cluster_energy;

        // close opened files per-run if you prefer. Keep global f2/f3 for multi-run
        f->Close();
        delete f;
    } // end run-loop

    // close text outputs
    if (f2) fclose(f2);
    if (f3) fclose(f3);
    if (f6) fclose(f6);
    if (f7) fclose(f7);

    sw_total.Stop();
    cout << "Total time: " << sw_total.RealTime() << " s (real), " << sw_total.CpuTime() << " s (cpu)\n";
    logmsg(INFO, "Finished NPS diagnostic analysis.");
}
