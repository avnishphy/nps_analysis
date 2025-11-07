#include "utils.C"  // logmsg(), trim(), etc.

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <cmath>
#include <utility>
#include <algorithm>

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
                             const TString &outPlotDir_in="output/plots/efficiency/",
                             const TString &runlistFile="config/runlist_x60_4b.txt") {

    TStopwatch sw_total;
    sw_total.Start();
    logmsg(INFO, "=========== HMS Tracking Efficiency ===========");

    TString skimDir = skimDir_in.EndsWith("/") ? skimDir_in : skimDir_in + "/";
    TString outPlotDir = outPlotDir_in.EndsWith("/") ? outPlotDir_in : outPlotDir_in + "/";
    gSystem->mkdir(outPlotDir, true);

    vector<int> runs = readRunList(runlistFile.Data());
    if (runs.empty()) {
        logmsg(ERROR, "No runs found!");
        return;
    }



    // // ============================================================ Open output ROOT file
    // TString outfile = outPlotDir + "hms_tracking_efficiency_extended.root";
    // TFile *fout = new TFile(outfile, "RECREATE");
    // fout->mkdir("scaler_fits");  // subdirectory for per-run scaler fits

    // ============================================================ Loop over runs
    for (int run : runs) {
        TString infile = Form("%sskim_run%d.root", skimDir.Data(), run);
        if (gSystem->AccessPathName(infile)) {
            logmsg(WARN, Form("Skipping run %d: file not found.", run));
            continue;
        }

        logmsg(INFO, Form("Processing run %d", run));
        TFile *f = TFile::Open(infile, "READ");
        if (!f || f->IsZombie()) { logmsg(ERROR, Form("Error opening file for run %d", run)); continue; }

        TTree *T = (TTree*)f->Get("T");
        if (!T) { logmsg(ERROR, Form("Tree not found in run %d", run)); f->Close(); continue; }

        gROOT->SetBatch(kTRUE);

        // --- Branch variables ---
        Double_t HgtrX, HgtrTh, HgtrPh, hdelta ; 
        Double_t hdcx,hdcxp,hdcy,hdcyp ; 
        Double_t evType, HgtrP ;
        Double_t cointime, HhodStatus, starTime, fptime, paeronpe,paeropt[99],paeront[99],paeroptd[99],paerontd[99];
        Double_t ntime,xexit,yexit ; 
        Double_t hbeta, hcalepr, hcaletot, hcernpe ;
        Double_t sum1,sum2,counters[10],ethcntr[10];
        Double_t  helmpsa,helnega,helposa,hgtry,pgtry,pcaltot,pgdsc,gdtrk,trkeff,trkeffer,prf, hrf, prfraw, hrfraw, hztar, pztar,trkeffg,trkeffger,trkefft,trkeffter,htof1,htof2 ;
        Double_t  hgdsc,htrkeff,htrkeffer,phbig,ppbig,t1big,t2big,sum3,sum4,sum5,sum6,sum7,sum8 ;
        Double_t hel1,helpred,helrep,helmps,helnqrt,hcaltot ;  
        Int_t  helmpsi,helnegi,helposi,icc,chist[100],chistp[100],chistpt[100];
        Int_t chistpg[100],chisT[100],chistg[100], cltimehist6[100] ;
        Int_t t16diffh[100],td1h[100],td6h[100], td16h[100] ;
        Int_t ctxh[100][20], ctyh[100][20],eclhist[100][20];
        Int_t cteh[100][20] ;
        Int_t ctimehist[120], cltimehist[120], vcltimehist[100] ;
        Int_t evtypecnr[10],ptdcmulth[15],icluster;
        Double_t trig1tdc,trig6tdc,edtmtdc ;
        Double_t  hdchit1,hdchit2,hntrk,ngcorr ;
        Double_t  clusE[10000],clusX[10000],clusY[10000],clusT[10000] ;
        Double_t  vclusE[10000],vclusX[10000],vclusY[10000],vclusT[10000] ;
        Int_t vclusSize[10000] ;
        int vnclus ;
        Double_t nclust,block_e[2000],cluster_ID[2000] ;
        Double_t block_x[2000],block_y[2000],block_t[2000] ;
        Double_t ctpi1,ctpi2,hst,pst,hpst,avm[30];
        Double_t avmb[1080],sumb[1080],minion[1080],minioner[1080];
        Int_t nchist[100], mallhist[100], mbesthist[100] ;
        Int_t epihist[100], mvbesthist[100], msbesthist[100] ;
        Int_t cthist1[120],cthist6[120],cthist16[120],cthist0[120] ;  
        Int_t eoverphist[100],cerhist[100],blktimeh[1080][20] ;  
        Int_t blkmassh[1080][20],mbycolhist[100][30] ;  
        Int_t eloblkhist[1080][20], ehiblkhist[1080][20] ;
        Int_t nevnt=0, trg1,trg6,ic1,ic2,iclo,ichi,clth[100][5];
        Int_t nevntall=0;
        Double_t enorm = 1.35 * 0.93 ; 
        doubel eslope = 0;
        Double_t ctimepi, ctimepinew, ctimeK, ctimep,ctimeraw,tmp ;

        T->SetBranchStatus("*", 0);

        T->SetBranchStatus("H.dc.x_fp", 1); T->SetBranchAddress("H.dc.x_fp",  &hdcx); 
        T->SetBranchStatus("H.dc.y_fp", 1); T->SetBranchAddress("H.dc.y_fp",  &hdcy); 
        T->SetBranchStatus("H.dc.xp_fp", 1); T->SetBranchAddress("H.dc.xp_fp", &hdcxp); 
        T->SetBranchStatus("H.dc.yp_fp", 1); T->SetBranchAddress("H.dc.yp_fp", &hdcyp); 
        T->SetBranchStatus("H.gtr.y", 1); T->SetBranchAddress("H.gtr.y", &hgtry); 
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
        T->SetBranchStatus("H.hod.goodscinhit", 1); T->SetBranchAddress("H.hod.goodscinhit", &hgdsc);

        T->SetBranchStatus("T.helicity.hel", 1); T->SetBranchAddress("T.helicity.hel",&hel1) ;
        T->SetBranchStatus("T.helicity.helpred", 1); T->SetBranchAddress("T.helicity.helpred",&helpred) ;
        T->SetBranchStatus("T.helicity.helrep", 1); T->SetBranchAddress("T.helicity.helrep",&helrep) ; 
        T->SetBranchStatus("T.helicity.mps", 1); T->SetBranchAddress("T.helicity.mps",&helmps) ;
        T->SetBranchStatus("T.helicity.npqrt", 1); T->SetBranchAddress("T.helicity.nqrt",&helnqrt) ;
                            
        T->SetBranchStatus("H.hod.starttime", 1); T->SetBranchAddress("H.hod.starttime", &starTime);       

        T->SetBranchStatus("CTime.epCoinTime1_ROC1", 1); T->SetBranchAddress("CTime.epCoinTime1_ROC1", &ctpi1);
        T->SetBranchStatus("CTime.epCoinTime2_ROC1", 1); T->SetBranchAddress("CTime.epCoinTime2_ROC1", &ctpi2);
        T->SetBranchStatus("H.hod.fpHitsTime", 1); T->SetBranchAddress("H.hod.fpHitsTime", &fptime); 
        T->SetBranchStatus("H.dc.ntrack", 1); T->SetBranchAddress("H.dc.ntrack", &hntrk); 
        
        T->SetBranchStatus("NPS.cal.nclust", 1); T->SetBranchAddress("NPS.cal.nclust", &nclust); 
        T->SetBranchStatus("NPS.cal.clusE", 1); T->SetBranchAddress("NPS.cal.clusE", &clusE); 
        T->SetBranchStatus("NPS.cal.clusX", 1); T->SetBranchAddress("NPS.cal.clusX", &clusX); 
        T->SetBranchStatus("NPS.cal.clusY", 1); T->SetBranchAddress("NPS.cal.clusY", &clusY); 
        T->SetBranchStatus("NPS.cal.clusT", 1); T->SetBranchAddress("NPS.cal.clusT", &clusT); 
        T->SetBranchStatus("NPS.cal.fly.block_clusterID", 1); T->SetBranchAddress("NPS.cal.fly.block_clusterID",&cluster_ID);
        T->SetBranchStatus("NPS.cal.fly.e", 1); T->SetBranchAddress("NPS.cal.fly.e",&block_e) ;
        T->SetBranchStatus("NPS.cal.fly.goodAdcTdcDiffTime", 1); T->SetBranchAddress("NPS.cal.fly.goodAdcTdcDiffTime",&block_t) ;
        T->SetBranchStatus("NPS.cal.vtpClusE", 1); T->SetBranchAddress("NPS.cal.vtpClusE",&vclusE);
        T->SetBranchStatus("NPS.cal.vtpClusSize", 1); T->SetBranchAddress("NPS.cal.vtpClusSize",&vclusSize);
        T->SetBranchStatus("NPS.cal.vtpClusTime", 1); T->SetBranchAddress("NPS.cal.vtpClusTime",&vclusT);
        T->SetBranchStatus("NPS.cal.vtpClusY", 1); T->SetBranchAddress("NPS.cal.vtpClusY",&vnclus);
        T->SetBranchStatus("T.hms.npsTRIG1_tdcTimeRaw", 1); T->SetBranchAddress("T.hms.npsTRIG1_tdcTimeRaw",&trig1tdc);
        T->SetBranchStatus("T.hms.npsTRIG6_tdcTimeRaw", 1); T->SetBranchAddress("T.hms.npsTRIG6_tdcTimeRaw",&trig6tdc);
        T->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1); T->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw",&edtmtdc);

        Long64_t nentries = T->GetEntries();
        cout << "Number of entries in the run: " << endl;
        
        // ============================================================================
        // NPS Analysis Event Loop (in-progress diagnostic build)
        // ============================================================================

        const int N_BLOCKS = 1080;

        // --- Define histograms *before* the event loop (outside normally)
        TH2F* hEblock_low  = new TH2F("hEblock_low",  "Block Energy Occupancy (Low E);Cluster Energy [GeV];Block ID;Counts",
                                    25, 0.0, 0.5, N_BLOCKS, 0.5, N_BLOCKS + 0.5);
        TH2F* hEblock_high = new TH2F("hEblock_high", "Block Energy Occupancy (High E);Cluster Energy [GeV];Block ID;Counts",
                                    20, 0.0, 5.0, N_BLOCKS, 0.5, N_BLOCKS + 0.5);
        TH1F* hNBlocksPerCluster = new TH1F("hNBlocksPerCluster", "Number of Blocks per Cluster;N_{blocks};Counts", 20, 0, 20);
        TH2F* hClusterXY = new TH2F("hClusterXY", "Cluster Position;X [cm];Y [cm];Counts", 100, -100, 100, 100, -100, 100);
        TH2F* hClusterEvsT = new TH2F("hClusterEvsT", "Cluster Energy vs Time;Time [ns];Energy [GeV];Counts", 200, 0, 400, 100, 0, 5);

        TH1F* h_nclust      = new TH1F("h_nclust",      "Number of Clusters per Event;# Clusters;Counts", 100, 0, 100);
        TH1F* h_ctime       = new TH1F("h_ctime",       "Coincidence Time (ctpi1);ctpi1 - 50 [ns];Counts", 120, 0, 120);
        TH1F* h_clust_time  = new TH1F("h_clust_time",  "Cluster Time (All);clusT - 105 [ns];Counts", 100, 0, 100);
        TH1F* h_vclust_time = new TH1F("h_vclust_time", "Cluster Time (VTP);clusT - 105 [ns];Counts", 100, 0, 100);

        TH1F *h_pi0_mass = new TH1F("h_pi0_mass", "Invariant Mass of #pi^{0} Candidates; M_{#gamma#gamma} [GeV]; Counts", 100, 0, 0.3);
        TH1F *h_pi0_energy = new TH1F("h_pi0_energy", "Sum of Cluster Energies; E_{#gamma1}+E_{#gamma2} [GeV]; Counts", 100, 0, 10);
        TH1F *h_opening_angle = new TH1F("h_opening_angle", "Opening Angle between Photons; #theta [rad]; Counts", 100, 0, 0.3);
        TH2F *h_mass_vs_angle = new TH2F("h_mass_vs_angle", "M_{#gamma#gamma} vs #theta; #theta [rad]; M [GeV]", 100, 0, 0.3, 100, 0, 0.3);
        TH2F *h_mass_vs_energy = new TH2F("h_mass_vs_energy", "M_{#gamma#gamma} vs E_{tot}; E_{#gamma1}+E_{#gamma2} [GeV]; M [GeV]", 100, 0, 10, 100, 0, 0.3);


        // ============================================================================
        // Event Loop
        // ============================================================================
        for (int evnt = 0; evnt < nentries; evnt++) {
            if (evnt % 10000 == 0)
                cout << "Progress: " << evnt << " / " << nentries << " processing..." << endl;

            T->GetEntry(evnt);
            evtType = T->GetLeaf("fEvtHdr.fEvtType")->GetValue();

            // HMS sanity cut
            if (edtmtdc < 0.1 && hdelta > -15.0 && hdelta < 15.0) {

                // HMS angular + calorimeter cuts
                if (HgtrTh > -1 && HgtrTh < 1 && HgtrPh > -1 && HgtrPh < 1 && hcaltot > 0.001) {
                    if (evnt > 300 && evnt < 320)
                        printf("HMS %d %3.1f %6.2f %7.4f %7.4f %5.2f %5.2f\n",
                            evnt, evtType, hdelta, HgtrTh, HgtrPh, hcaletot, hcernpe);
                }

                // =========================================================================
                // Merge nearby clusters (spatial + temporal proximity)
                // =========================================================================
                if (nclust > 1) {
                    for (int n1 = 0; n1 < nclust - 1; n1++) {
                        if (clusE[n1] <= 0) continue;

                        for (int n2 = n1 + 1; n2 < nclust; n2++) {
                            if (clusE[n2] <= 0) continue;

                            double xdiff = clusX[n1] - clusX[n2];
                            double ydiff = clusY[n1] - clusY[n2];
                            double tdiff = fabs(clusT[n1] - clusT[n2]);
                            double space_diff_sq = xdiff * xdiff + ydiff * ydiff;

                            // Merge threshold (~7 cm spatial, <2 ns timing)
                            if (space_diff_sq < 50.0 && tdiff < 2.0) { // why are we taking the space distance limit for clusters < 7.5 cm?
                                double E1 = clusE[n1], E2 = clusE[n2];
                                double Etot = E1 + E2;

                                // Energy-weighted merge
                                clusX[n1] = (clusX[n1] * E1 + clusX[n2] * E2) / Etot; // some kind of energy weighted positioning...?
                                clusY[n1] = (clusY[n1] * E1 + clusY[n2] * E2) / Etot;
                                clusT[n1] = (clusT[n1] * E1 + clusT[n2] * E2) / Etot;
                                clusE[n1] = Etot;

                                clusE[n2] = 0.0; // mark as merged/invalid

                                // Reassign blocks; removing the other cluster and giving it first cluster's IDs.
                                for (int i = 0; i < N_BLOCKS; i++) {
                                    if (cluster_ID[i] == n2)
                                        cluster_ID[i] = n1;
                                }
                            }
                        }
                    }
                }

                // =========================================================================
                // Cluster energy normalization (optional) GPT says that it should be only if the energies are uncalibrated.
                // If the energy is already calibrated then we shouldn't use this. 
                // Also, how is the "enorm" and "eslope" calculated?
                // =========================================================================
                double ehi = 0.0;
                int nhi = -1;
                // clusE[i] *= enorm*(1 + eslope*(clusX[i] + 30)); apply this only after confirming why this is applied!!
                for (int i = 0; i < nclust; i++) {
                    if (clusE[i] > ehi) {
                        nhi = i;
                        ehi = clusE[i];
                    }
                }

                // Normalize block energies?? (only if enorm defined)
                for (int i = 0; i < N_BLOCKS; i++) {
                    block_e[i] *= enorm;
                }

                // =========================================================================
                // NPS Block-Cluster Correlations
                // =========================================================================
                for (int i = 0; i < N_BLOCKS; i++) {
                    int icluster = (int)cluster_ID[i];
                    if (icluster < 0 || icluster >= nclust) continue;

                    double e_block = block_e[i];
                    double e_cluster = clusE[icluster];
                    double t_cluster = clusT[icluster];

                    if (e_cluster <= 0 || !isfinite(e_block)) continue;

                    if (e_block > 0.3 * e_cluster && e_block < 0.8 * e_cluster &&
                        t_cluster > 0. && t_cluster < 300.) {

                        if (e_cluster >= 0.1 && e_cluster < 0.5)
                            hEblock_low->Fill(e_cluster, i + 1);

                        if (e_cluster >= 0.25 && e_cluster < 5.0)
                            hEblock_high->Fill(e_cluster, i + 1);
                    }
                }

                // =========================================================================
                // Cluster-level Histograms
                // =========================================================================
                int nBlocksInCluster = 0;
                for (int i = 0; i < nclust; i++) {
                    if (clusE[i] <= 0) continue;

                    hClusterXY->Fill(clusX[i], clusY[i]);
                    hClusterEvsT->Fill(clusT[i], clusE[i]);

                    // count associated blocks
                    nBlocksInCluster = 0;
                    for (int j = 0; j < N_BLOCKS; j++) {
                        if ((int)cluster_ID[j] == i && block_e[j] > 0)
                            nBlocksInCluster++;
                    }
                    hNBlocksPerCluster->Fill(nBlocksInCluster);
                }

                // =========================================================================
                // NPS Debug print for early events
                // =========================================================================
                if (evnt > 10 && evnt < 30 && nclust > 1) {
                    printf("------------------------------------------------------------\n");
                    printf("[Event %d Diagnostics]\n", evnt);
                    printf("ct= %8.2f %8.2f  TDCs: trig1=%8.0f trig6=%8.0f  edtm=%8.0f  start=%8.0f  fp=%8.0f\n",
                        ctpi1, ctpi2, trig1tdc, trig6tdc, edtmtdc, starttime, fptime);
                    printf("evtType= %.1f   nclust= %d   vnclus= %d\n",
                        evtType, nclust, vnclus);

                    if (nclust > 0) {
                        for (int i = 0; i < nclust; i++) {
                            if (clusE[i] > 0.3) {
                                printf("  Cluster %2d | E = %6.3f  X = %6.1f  Y = %6.1f  T = %6.1f\n",
                                    i, clusE[i], clusX[i], clusY[i], clusT[i]);

                                // Loop over all blocks and find those belonging to this cluster
                                int nBlocksInCluster = 0;
                                for (int j = 0; j < N_BLOCKS; j++) {
                                    int icluster = (int)cluster_ID[j];
                                    if (icluster == i && block_e[j] > 0) {
                                        if (nBlocksInCluster == 0)
                                            printf("    [Blocks contributing to cluster %d]\n", i);
                                        
                                        printf("      Block %4d | icluster=%3d  E=%.4f  T=%.4f\n",
                                            j, icluster, block_e[j], block_t[j]);
                                        nBlocksInCluster++;
                                    }
                                }

                                if (nBlocksInCluster == 0)
                                    printf("    [No blocks found for this cluster]\n");
                            }
                        }
                    }
                    printf("------------------------------------------------------------\n");
                }

                // --------------------------------------------------
                // Fill ncluster histogram
                // --------------------------------------------------
                if (nclust >= 0 && nclust < 100) h_nclust->Fill(nclust);
                else if (nclust >= 100) h_nclust->Fill(99);

                // --------------------------------------------------
                // Fill coincidence time histogram
                // --------------------------------------------------
                double ctime_bin = ctpi1 - 50;
                if (ctime_bin >= 0 && ctime_bin < 120) h_ctime->Fill(ctime_bin);

                // --------------------------------------------------
                // Cluster time (all events)
                // --------------------------------------------------
                for (int n1 = 0; n1 < nclust; n1++) {
                    if (clusE[n1] > 0.3) {
                        double t_bin = clusT[n1] - 105;
                        if (t_bin >= 0 && t_bin < 100) h_clust_time->Fill(t_bin);
                    }
                }

                // --------------------------------------------------
                // Cluster time (VTP events)
                // --------------------------------------------------
                for (int n1 = 0; n1 < vnclus; n1++) {
                    if (vclusE[n1] > 400) {
                        double t_bin = vclusT[n1] - 105;
                        if (t_bin >= 0 && t_bin < 100) h_vclust_time->Fill(t_bin);
                    }
                }

                // --------------------------------------------------
                // π⁰ Reconstruction
                // chooses the photon candidates giving mass closest to expected pi0 mass
                // --------------------------------------------------
                const double dnps = 300.0;        // NPS distance from target [cm], depends on kinematics
                const double pi0_mass = 0.135;  // π⁰ mass [GeV/c²]
                const double time_cut = 3.0;    // Cluster time difference [ns]
                const double e_cut = 0.1;       // Minimum cluster energy [GeV]

                int n1best = -1, n2best = -1;
                double ebest = 0.0;
                double mbest = 0.0;
                double_opening_angle_best = 0;
                double massmin = 1e6;  // large initial difference

                if (nclust > 1) {

                    for (int n1 = 0; n1 < nclust - 1; n1++) {
                        if (clusE[n1] < e_cut) continue;

                        for (int n2 = n1 + 1; n2 < nclust; n2++) {
                            if (clusE[n2] < e_cut) continue;

                            double tdiff = fabs(clusT[n1] - clusT[n2]);
                            if (tdiff > time_cut) continue;

                            // Spatial separation in detector plane
                            double xdiff = clusX[n1] - clusX[n2];
                            double ydiff = clusY[n1] - clusY[n2];
                            double r_sep = sqrt(xdiff * xdiff + ydiff * ydiff);
                            if (r_sep < 0.03) continue; // avoid same-cluster double counting

                            // Opening angle (small-angle approx)
                            double opening_theta = r_sep / dnps;  // radians
                            double sin2_half = pow(sin(opening_theta / 2.0), 2);

                            // π⁰ invariant mass (two-photon system)
                            double arg = 4.0 * clusE[n1] * clusE[n2] * sin2_half;
                            double mass = arg > 0 ? sqrt(arg) : 0.0;

                            // Histogram
                            int icc = static_cast<int>(mass / pi0_mass * 50.0);
                            icc = std::clamp(icc, 0, 99);
                            mallhist[icc]++;

                            // Choose best (closest to true π⁰ mass)
                            double mdiff = fabs(mass - pi0_mass);
                            if (mdiff < massmin) {
                                massmin = mdiff;
                                ebest = clusE[n1] + clusE[n2];
                                mbest = mass;
                                n1best = n1;
                                n2best = n2;
                                opening_theta_best = opening_theta;
                            }

                #ifdef DEBUG_PI0
                            if (evnt < 10)
                                printf("n1=%d n2=%d | Xdiff=%6.2f Ydiff=%6.2f | θ=%.3f rad | M=%.3f GeV | ΔM=%.3f\n",
                                    n1, n2, xdiff, ydiff, opening_theta, mass, mdiff);
                #endif
                        }
                    }

                    // Fill best-pair histograms
                    if (n1best >= 0 && n2best >= 0) {
                        h_pi0_mass->Fill(mbest);
                        h_pi0_energy->Fill(ebest);
                        h_opening_angle->Fill(opening_theta_best);
                    }
                }

                // Constants for cuts
                const double t_low = 148.0;
                const double t_high = 152.0;
                const double t_wide_low = 140.0;
                const double t_wide_high = 164.0;
                const double e_cut_hist = 0.6;
                const double xdiff_cut = 5.0;
                const double xav_offset = 30.0;
                const double col_width = 2.0;

                // mbesthist	Mass histogram of best π⁰ candidate
                // mvbesthist	Mass histogram for clusters in narrow time window (~150 ns)
                // msbesthist	Mass histogram for clusters in narrow time + energy cut
                // mbycolhist	Mass vs. detector X-position mapping (to check geometry)
                // epihist	Total π⁰ candidate energy in a mass window around 135 MeV
                // ethcntr	Energy counters for different cluster thresholds in broad time + mass window
                
                // Only proceed if a best π⁰ candidate exists
                if (ebest > 0.0) {

                    bool t_best_narrow = (clusT[n1best] > t_low && clusT[n1best] < t_high) &&
                                        (clusT[n2best] > t_low && clusT[n2best] < t_high);

                    bool e_best_cut = (clusE[n1best] > e_cut_hist && clusE[n2best] > e_cut_hist);

                    // Optional debug print
                    if (evnt > 7990 && evnt < 8000 && e_best_cut)
                        printf("%.3f %.3f %.1f %.1f \n", mbest, ebest, clusT[n1best], clusT[n2best]);

                    // Histogram bin index for mass
                    int icc = static_cast<int>(mbest / 0.135 * 50.0);
                    icc = std::clamp(icc, 0, 99);
                    mbesthist[icc]++;

                    // Narrow time window histogram
                    if (t_best_narrow) mvbesthist[icc]++;

                    // Narrow time + energy cut histogram
                    if (t_best_narrow && e_best_cut) {
                        msbesthist[icc]++;

                        // X difference / column mapping
                        double xdiff = clusX[n1best] - clusX[n2best];
                        double xav = (clusX[n1best] + clusX[n2best]) / 2.0;
                        if (fabs(xdiff) < xdiff_cut) {
                            int icol = static_cast<int>((xav + xav_offset) / col_width);
                            icol = std::clamp(icol, 0, 30);
                            mbycolhist[icc][icol]++;
                        }

                        // π⁰ energy histogram
                        int ie = static_cast<int>(ebest * 10.0);
                        ie = std::clamp(ie, 0, 99);
                        if (mbest > 0.11 && mbest < 0.16) epihist[ie]++;

                        // Broad time + mass window counters
                        if ((clusT[n1best] > t_wide_low && clusT[n1best] < t_wide_high) &&
                            (clusT[n2best] > t_wide_low && clusT[n2best] < t_wide_high) &&
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

                                // Constants
                const double hdelta_low = -18.0;
                const double hdelta_high = 18.0;
                const double hcal_cut = 0.6;
                const double hcern_cut = 0.5;
                const double clusE_min = 0.0;
                const int max_blocks = 1080;
                const double t_low = 148.0;
                const double t_high = 152.0;

                // Electron in HMS cut
                if (hdelta > hdelta_low && hdelta < hdelta_high && hcaletot > hcal_cut && hcernpe > hcern_cut) {

                    counters[2]++;

                    // Basic histogram indices
                    int icc = static_cast<int>(ctpi1 - 50);
                    icc = std::clamp(icc, 0, 119);

                    int ic1 = static_cast<int>(clusE[n1best] / 0.05);
                    int ic2 = static_cast<int>(clusE[n2best] / 0.05);
                    ic1 = std::clamp(ic1, 0, 99);
                    ic2 = std::clamp(ic2, 0, 99);

                    int ichi = (clusE[n1best] > clusE[n2best]) ? ic1 : ic2;
                    int iclo = (clusE[n1best] > clusE[n2best]) ? ic2 : ic1;

                    if (ichi < 0 || iclo < 0) 
                        printf("error %d %d %.0f %.0f\n", iclo, ichi, trig1tdc, trig6tdc);

                    // Cluster time difference
                    int j1 = static_cast<int>(clusT[n1best] - ctpi1);
                    int j2 = static_cast<int>(clusT[n2best] - ctpi1);
                    j1 = std::clamp(j1, 0, 99);
                    j2 = std::clamp(j2, 0, 99);

                    // Trigger flags
                    bool trg1 = (trig1tdc > 31300 && trig1tdc < 35000);
                    bool trg6 = (trig6tdc > 31300 && trig6tdc < 35000);

                    // Histogram filling for best π⁰ candidates
                    if (mbest > 0.11 && mbest < 0.16 && nclust > 1) {
                        if (trg6 && !trg1) {
                            eclhist[iclo][0]++;
                            eclhist[ichi][1]++;
                            td1h[j1]++;
                            td16h[j2]++;
                        }
                        if (trg1 && !trg6) {
                            eclhist[iclo][2]++;
                            eclhist[ichi][3]++;
                            td6h[j1]++;
                            td6h[j2]++;
                        }
                    }

                    // Coincidence time histograms
                    if (trg1 && !trg6) cthist1[icc]++;
                    if (trg1 && trg6) {
                        cthist16[icc]++;
                        int icc_diff = static_cast<int>((trig1tdc - trig6tdc)/10.0 + 50);
                        icc_diff = std::clamp(icc_diff, 0, 99);
                        t16diffh[icc_diff]++;
                    }
                    if (trg6 && !trg1) cthist6[icc]++;

                    // Cluster times for trg6 only
                    if (trg6 && !trg1) {
                        for (int n = 0; n < nclust; n++) {
                            if (clusE[n] > 0.3) {
                                int tbin = std::clamp(static_cast<int>(clusT[n] - 105), 0, 99);
                                cltimehist6[tbin]++;
                            }
                        }
                    }

                    // Cluster times for good π⁰ mass
                    if (mbest > 0.11 && mbest < 0.16 && nclust > 1) {
                        int ichi_best = (clusE[n1best] > clusE[n2best]) ? n1best : n2best;
                        int iclo_best = (clusE[n1best] > clusE[n2best]) ? n2best : n1best;

                        int t1 = std::clamp(static_cast<int>(clusT[n1best] - 100), 0, 99);
                        int t2 = std::clamp(static_cast<int>(clusT[n2best] - 100), 0, 99);
                        int tav = std::clamp(static_cast<int>((clusT[n1best] + clusT[n2best]) / 2 - 100), 0, 99);

                        clth[t1][0]++;
                        clth[t2][1]++;
                        clth[tav][2]++;
                    }

                    // Temporary ctpi2 placeholder
                    ctpi2 = 0.0;

                    // Skim file writing for possible π⁰ events
                    if (ctpi2 > -30.0 && ctpi2 < 30.0) {
                        counters[3]++;
                        auto isValidCluster = [](double val){ return std::isfinite(val); };
                        if (nclust > 1 && mbest > 0.05 &&
                            isValidCluster(clusE[n1best]) && isValidCluster(clusE[n2best]) &&
                            isValidCluster(clusX[n1best]) && isValidCluster(clusX[n2best]) &&
                            isValidCluster(clusY[n1best]) && isValidCluster(clusY[n2best]) &&
                            isValidCluster(clusT[n1best]) && isValidCluster(clusT[n2best])) {

                            counters[4]++;
                            fprintf(f2,
                                "%6.2f %1d %1d %1d %6.2f %7.4f %7.4f %6.2f %7.4f %7.4f %6.2f %6.2f %7.4f %7.4f %7.4f %7.2f %7.2f %7.2f %7.4f %7.2f %7.2f %7.2f %7.3f  %7.3f   \n",
                                ctpi2,(int)hel1+1,(int)helrep+1,(int)helmps+1,
                                hdelta,HgtrTh,HgtrPh,hcernpe,hcalepr,hcaletot,
                                hdcx, hdcy, hdcxp, hdcyp,
                                clusE[n1best],clusX[n1best],clusY[n1best],clusT[n1best],
                                clusE[n2best],clusX[n2best],clusY[n2best],clusT[n2best],
                                mbest,hztar);

                            // Loop over blocks for both clusters
                            for (int j = 0; j < 2; j++) {
                                int nn = (j == 0) ? n1best : n2best;
                                double etot = clusE[n1best] + clusE[n2best];

                                for (int i = 0; i < max_blocks; i++) {
                                    int icluster = static_cast<int>(cluster_ID[i]);
                                    if (icluster == nn && block_e[i] > 0.0) {
                                        fprintf(f2,"%4d %3d %8.4f %8.1f \n", i, icluster, block_e[i], block_t[i]);

                                        if (block_e[i] > 0.3 * etot) {
                                            int jj = std::clamp(static_cast<int>(block_t[i] - 142), 0, 19);
                                            blktimeh[i][jj]++;

                                            jj = std::clamp(static_cast<int>((mbest - 0.115) * 20.0 / 0.040), 0, 19);
                                            blkmassh[i][jj]++;
                                        }
                                    }
                                }
                            }
                            fprintf(f2,"  -1     0.000 \n");
                        }
                    }

                    // DVCS skim writing (optional)
                    if (nhi > -1 && ehi > 2.0 && mbest < 0.05) {
                        counters[5]++;
                        fprintf(f3,
                            "%6.2f %1d %1d %1d %6.2f %7.4f %7.4f %6.2f %7.4f %7.4f %6.2f %6.2f %7.4f %7.4f %7.4f %7.2f %7.2f %7.2f %7.3f   \n",
                            ctpi2,(int)hel1+1,(int)helrep+1,(int)helmps+1,
                            hdelta,HgtrTh,HgtrPh,hcernpe,hcalepr,hcaletot,
                            hdcx, hdcy, hdcxp, hdcyp,
                            clusE[nhi],clusX[nhi],clusY[nhi],clusT[nhi],
                            mbest);

                        for (int i = 0; i < max_blocks; i++) {
                            int icluster = static_cast<int>(cluster_ID[i]);
                            if (icluster == nhi && block_e[i] > 0.0)
                                fprintf(f3,"%4d %3d %8.4f %8.1f \n", i, icluster, block_e[i], block_t[i]);
                        }
                    }
                }


                
                






            } // end HMS + NPS cut
        } // end event loop

        for(icc=0 ; icc<100 ; icc++){
    printf("%2d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d \n",
      icc,nchist[icc],mallhist[icc],mbesthist[icc],
      mvbesthist[icc],msbesthist[icc],
      epihist[icc],cthist1[icc],cthist6[icc],
	   cthist16[icc],cltimehist[icc],cltimehist6[icc]);
	   //ctimehist[icc],cltimehist[icc],vcltimehist[icc]) ;
  }
  /*  for(icc=0 ; icc<100 ; icc++){
    printf("%2d %4d %4d %4d %4d %4d %4d %4d %5d  \n",
	   icc,t16diffh[icc],td1h[icc],td16h[icc],td6h[icc],
           eclhist[icc][0],eclhist[icc][1],eclhist[icc][2],
           eclhist[icc][3]) ;
	   }*/
  /*  for(int i=0 ; i<100 ; i++){
    printf("%2d %4d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d\n",
  i,ctxh[i][0],ctxh[i][1],ctxh[i][2], ctxh[i][3], ctxh[i][4], ctxh[i][5],
    ctxh[i][6], ctxh[i][7], ctxh[i][8], ctxh[i][9], ctxh[i][10],ctxh[i][11],
    ctxh[i][12],ctxh[i][13],ctxh[i][14],ctxh[i][15],ctxh[i][16],ctxh[i][17],
    ctxh[i][18],ctxh[i][19]) ;
  }
  for(int i=0 ; i<100 ; i++){
    printf("%2d %4d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d\n",
  i,ctyh[i][0],ctyh[i][1],ctyh[i][2], ctyh[i][3], ctyh[i][4], ctyh[i][5],
    ctyh[i][6], ctyh[i][7], ctyh[i][8], ctyh[i][9], ctyh[i][10],ctyh[i][11],
    ctyh[i][12],ctyh[i][13],ctyh[i][14],ctyh[i][15],ctyh[i][16],ctyh[i][17],
    ctyh[i][18],ctyh[i][19]) ;
  }

  for(int i=0 ; i<100 ; i++){
    printf("%2d %4d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d %2d\n",
  i,cteh[i][0],cteh[i][1],cteh[i][2], cteh[i][3], cteh[i][4], cteh[i][5],
    cteh[i][6], cteh[i][7], cteh[i][8], cteh[i][9], cteh[i][10],cteh[i][11],
    cteh[i][12],cteh[i][13],cteh[i][14],cteh[i][15],cteh[i][16],cteh[i][17],
    cteh[i][18],cteh[i][19]) ;
  }
  for(icc=0 ; icc<99 ; icc++){
    printf(" %2d %5d %5d %5d %5d  \n",icc,clth[icc][0],
	   clth[icc][1],clth[icc][2],clth[icc][3]) ;
  }
  */

  for(icc=0 ; icc<6 ; icc++){
    printf(" %7.0f %7.3f \n",counters[icc],
	   counters[icc]/counters[0]) ;
  }
  for(icc=0 ; icc<6 ; icc++){
    printf(" %.2f %7.0f %7.3f \n",0.1*(icc+5),ethcntr[icc],
	   ethcntr[icc]/ethcntr[0]) ;
  }
  for(int i=0 ; i<30 ; i++){
    avm[i]=0 ;
    double sum = 0. ; 
    for(int j=40 ; j<60 ; j++) {
      double mm = 0.135/50.*(j+0.5) ;
      avm[i] = avm[i] + mm * mbycolhist[j][i] ; 
      sum += mbycolhist[j][i] ; 
    }
    if(sum > 0.) avm[i] = avm[i]/sum ;
    printf("%3d %.3f",i,avm[i]) ;
    for(int j=43 ; j<56 ; j++) printf(" %2d",mbycolhist[j][i]) ;
    printf(" %2d \n",mbycolhist[57][i]) ;
  }


// by block outputs
// low energy clusters by block
  for(int i=0 ; i<1080 ; i++){
    fprintf(f7,"%4d ",i) ;
    for(int j=0 ; j<19 ; j++) fprintf(f7," %2d",eloblkhist[i][j]) ;
    fprintf(f7," %2d \n",eloblkhist[i][19]) ;
    double sum1 = 0.;
    double sum2 = 0.;
    for(int j=0 ; j<19 ; j++) {
      double ebin = 0.1 + 0.02*(j+05) ;
      sum1 += eloblkhist[i][j] ;
      sum2 += ebin * eloblkhist[i][j] ;
    }
    minion[i]=0. ;
    minioner[i]=0. ;
    if(sum1 > 20) {
      minion[i] = 1000. * sum2/sum1 ;
      minioner[i]= 50./sqrt(sum1) ;
    }
  }
  for(int i=0 ; i<1080 ; i++){
    fprintf(f7,"%4d ",i) ;
    for(int j=0 ; j<19 ; j++) fprintf(f7," %2d",ehiblkhist[i][j]) ;
    fprintf(f7," %2d \n",ehiblkhist[i][19]) ;
  }
  for(int i=0 ; i<1080 ; i++){
    fprintf(f7,"%4d ",i) ;
    for(int j=0 ; j<19 ; j++) fprintf(f7," %2d",blktimeh[i][j]) ;
    fprintf(f7," %2d \n",blktimeh[i][19]) ;
    //printf(" %4d %2d \n",i,blktimeh[i][19]) ;
  }
  for(int i=0 ; i<1080 ; i++){
    avmb[i]=0. ;
    sumb[i]=0. ;
    for(int j=0 ; j<20 ; j++) {
      double mm = 0.115 + 0.040 * (j+0.5)/20. ;
      avmb[i] = avmb[i] + mm * blkmassh[i][j] ; 
      sumb[i] += blkmassh[i][j] ; 
    }
    fprintf(f7,"%4d ",i) ;
    for(int j=0 ; j<19 ; j++) 
      fprintf(f7," %2d",blkmassh[i][j]) ;
    fprintf(f7," %2d \n",blkmassh[i][19]) ;
    //printf(" %4d %2d \n",i,blktimeh[i][19]) ;
  }

//topdrawer plots
  fprintf(f6,"set device postscript \n") ;
// pi0 mass
  fprintf(f6," set window x 1 of 3 y 1 of 2\n") ;
  fprintf(f6," title bottom 'M0GG1 (GeV)' \n") ;
  fprintf(f6," case         ' XGGX'\n") ;
  fprintf(f6," title top 'run %d' \n",runNumber) ;
  for(icc=20 ; icc<80 ; icc++){
    fprintf(f6,"%.3f %d\n",0.135/50.*(icc+0.5),mbesthist[icc]) ;
  }
  fprintf(f6,"hist\n") ;
  for(icc=20 ; icc<80 ; icc++){
    fprintf(f6,"%.3f %d\n",0.135/50.*(icc+0.5),mvbesthist[icc]) ;
  }
  fprintf(f6,"hist\n") ;
  for(icc=20 ; icc<80 ; icc++){
    fprintf(f6,"%.3f %d\n",0.135/50.*(icc+0.5),msbesthist[icc]) ;
  }
  fprintf(f6,"hist\n") ;
// pi0 energy
  fprintf(f6," set window x 2 of 3 y 1 of 2\n") ;
  fprintf(f6," title bottom 'E0GG1 (GeV)' \n") ;
  fprintf(f6," case         ' XGGX'\n") ;
  fprintf(f6," title top 'run %d' \n",runNumber) ;
  for(icc=1 ; icc<80 ; icc++){
    fprintf(f6,"%.3f %d\n",0.1*(icc-0.5),epihist[icc]) ;
  }
  fprintf(f6,"hist\n") ;
// Cointime
  fprintf(f6," set window x 1 of 3 y 2 of 2\n") ;
  fprintf(f6," title bottom 'Cluster time (nsec)' \n") ;
  fprintf(f6," set scale y log \n") ;
  fprintf(f6," title top 'black: lo E  blue: hi E' \n") ;
  for(icc=1 ; icc<100 ; icc++){
    if(clth[icc][0]>1)
    fprintf(f6,"%.3f %d\n",1.0*(icc-0.5)+100.,clth[icc][0]) ;
  }
  fprintf(f6,"hist\n") ;
  fprintf(f6," set color blue \n") ;
  for(icc=1 ; icc<119 ; icc++){
    if(clth[icc][1]>1)
    fprintf(f6,"%.3f %d\n",1.0*(icc-0.5)+100.,clth[icc][1]) ;
  }
  fprintf(f6,"hist\n") ;

  fprintf(f6," set color white \n") ;
  fprintf(f6," set window x 2 of 3 y 2 of 2\n") ;
  fprintf(f6," title bottom 'Coin time (nsec)' \n") ;
  fprintf(f6," set scale y log \n") ;
  fprintf(f6," title top 'average of 2 clusters' \n") ;
  for(icc=1 ; icc<100 ; icc++){
    if(clth[icc][2]>1)
    fprintf(f6,"%.3f %d\n",1.0*(icc-0.5)+55.,clth[icc][2]) ;
  }
  fprintf(f6,"hist\n") ;
  /*  fprintf(f6," set color blue \n") ;
  for(icc=1 ; icc<119 ; icc++){
    fprintf(f6,"%.3f %d\n",1.0*(icc-0.5)+55.,clth[icc][3]) ;
    }*/
  fprintf(f6,"hist\n") ;
  fprintf(f6," set color white \n") ;
  fprintf(f6," set window x 3 of 3 y 2 of 2\n") ;
  fprintf(f6," title bottom 'X NPS (cm)' \n") ;
  fprintf(f6," set limits y 0.12 0.15 \n") ;
  fprintf(f6," title top 'pi0 peak with N0=%.2f  slope=%.2f' \n",enorm,eslope*60.) ;
  for(icc=0 ; icc<30 ; icc++){
    fprintf(f6,"%.3f %.4f\n",-30. + 2.*(icc-0.5),avm[icc]) ;
  }
  fprintf(f6,"hist\n") ;
  fprintf(f6," set color white \n") ;
  fprintf(f6," set window x 3 of 3 y 1 of 2\n") ;
  fprintf(f6," title bottom 'block number' \n") ;
  fprintf(f6," set limits y 100. 500. \n") ;
  fprintf(f6," set sym 9O size 0.5 ; set bar size 0. \n") ;
  fprintf(f6," set order x y dy \n") ;
  fprintf(f6," title top 'min ion peak (MeV)' \n") ;
  for(icc=0 ; icc<1080 ; icc++){
    /*    if(sumb[icc]>30.){
    double avval = avmb[icc]/sumb[icc] ;
    double averr = 0.020 / sqrt(sumb[icc]) ;
    fprintf(f6,"%d %.4f %.4f\n",icc,avval,averr) ; 
    } */
    fprintf(f6,"%d %.4f %.4f\n",icc,minion[icc],minioner[icc]) ;
  }
  fprintf(f6,"plot\n") ;

  // plot mass versus column

  fprintf(f6,"new frame\n") ;
  int ix=0 ; int iy = 4 ;
  for(icc=5 ; icc<30 ; icc++){
    ix = ix + 1; 
    if(ix>6) {ix=1 ; iy = iy-1 ; }
    fprintf(f6," set window x %d of 6 y %d of 4\n",ix,iy) ;
    fprintf(f6," title bottom 'M0GG1 (GeV)' \n") ;
    fprintf(f6," case         ' XGGX'\n") ;
    fprintf(f6," title top 'column %d' \n",icc+1) ;
    for(int j=38 ; j<62 ; j++) { 
      fprintf(f6,"%.4f %4d \n",0.135/50.*(j+0.5),mbycolhist[j][icc]) ;
    }
    fprintf(f6,"hist\n") ;
    fprintf(f6," 0.135 0. ; 0.135 1000000. ; join dash\n") ;
  }

  fprintf(f2,"  0.00  0.00 0.00\n") ;

  fclose(f2) ;
  //  fclose(f5) ;
  fclose(f6) ;
  } else {
      printf("--------\n") ;
      printf("run %d not FOUND\n",runNumber) ;
      printf("--------\n") ;

    gSystem->Exit(1) ;
  } 

  gROOT
->SetBatch(kFALSE);

  gSystem->Exit(1) ; 

        // ============================================================================
        // Save Results
        // ============================================================================
        TFile* fout = new TFile("output/nps_diagnostics.root", "RECREATE");
        hEblock_low->Write();
        hEblock_high->Write();
        hNBlocksPerCluster->Write();
        hClusterXY->Write();
        hClusterEvsT->Write();
        fout->Close();

        // Visualization example
        TCanvas* cBlock = new TCanvas("cBlock", "NPS Block Occupancy", 1200, 800);
        cBlock->Divide(1, 2);
        cBlock->cd(1); gPad->SetLogz(); hEblock_low->Draw("COLZ");
        cBlock->cd(2); gPad->SetLogz(); hEblock_high->Draw("COLZ");
        cBlock->SaveAs("output/block_cluster_energy_occupancy.pdf");

        TCanvas* cTime = new TCanvas("cTime", "NPS Timing Diagnostics", 1000, 800);
        cTime->Divide(2,2);

        cTime->cd(1);
        h_nclust->Draw();

        cTime->cd(2);
        h_ctime->Draw();

        cTime->cd(3);
        h_clust_time->Draw();

        cTime->cd(4);
        h_vclust_time->Draw();

        cTime->SaveAs("output/nps_timing_diagnostics.pdf");

        // Optionally save histograms to ROOT file
        TFile* fout = new TFile("output/nps_timing_diagnostics.root", "RECREATE");
        h_nclust->Write();
        h_ctime->Write();
        h_clust_time->Write();
        h_vclust_time->Write();
        fout->Close();


    
    
    
    }


}
