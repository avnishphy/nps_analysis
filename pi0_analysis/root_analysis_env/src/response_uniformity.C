#include "utils.C" // for readRunList(), logmsg(), trim(), etc.

#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TLegend.h>

#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

// ============================================================
// Runlist reader (keeps your original behaviour but with trim)
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
// response_uniformity: mean cluster energy vs X and vs Y
// - per-run: 2D histos, profile, linear fit (central region)
// - global: accumulate profiles and fit
// - output: per-run root + pdf, global root + pdf, CSV summary
// ============================================================
void response_uniformity(const TString &skimDir_in="output/skimmed/",
                         const TString &outDir_in="output/plots/response/",
                         const TString &runlistFile="config/runlist_x60_4b.txt")
{
    TStopwatch sw; sw.Start();
    logmsg(INFO, "======= NPS response uniformity (E vs X and E vs Y) =======");

    TString skimDir = skimDir_in.EndsWith("/") ? skimDir_in : skimDir_in + "/";
    TString outDir  = outDir_in.EndsWith("/") ? outDir_in : outDir_in + "/";
    gSystem->mkdir(outDir, true);

    // Read runs
    vector<int> runs = readRunList(runlistFile.Data());
    if (runs.empty()) { logmsg(ERROR, "No runs found in runlist"); return; }

    // Global accumulators:
    // profiles (mean E vs position) and 2D scatter to show density
    const Int_t nXBins = 30;   // X profile bins (adapt to geometry)
    const Double_t xMin = -30.0, xMax = 30.0;
    const Int_t nYBins = 36;   // Y profile bins
    const Double_t yMin = -36.0, yMax = 36.0;
    const Int_t nEBins = 100;  // energy histogram resolution

    TProfile *pEvsX_global = new TProfile("pEvsX_global", "Mean cluster E vs X (global);X (cm);Mean E (GeV)",
                                          nXBins, xMin, xMax);
    TProfile *pEvsY_global = new TProfile("pEvsY_global", "Mean cluster E vs Y (global);Y (cm);Mean E (GeV)",
                                          nYBins, yMin, yMax);

    TH2D *h2EvsX_global = new TH2D("h2EvsX_global", "E vs X (global);X (cm);E (GeV)",
                                  nXBins, xMin, xMax, nEBins, 0.0, 4.0);
    TH2D *h2EvsY_global = new TH2D("h2EvsY_global", "E vs Y (global);Y (cm);E (GeV)",
                                  nYBins, yMin, yMax, nEBins, 0.0, 4.0);

    // cluster energy spectrum and multiplicity per event (global)
    TH1D *hE_global = new TH1D("hE_global", "Cluster E (global);E (GeV);counts", nEBins, 0.0, 4.0);
    TH1D *hNclust_global = new TH1D("hNclust_global", "Clusters per event (global);nclusters;events", 21, -0.5, 20.5);

    // CSV summary header (two axes)
    TString csvFile = outDir + "response_uniformity_summary.csv";
    ofstream csv(csvFile.Data());
    csv << "Run,"
        << "slopeX_GeV_per_cm,slopeX_err,interceptX_GeV,interceptX_err,chi2_ndf_X,"
        << "slopeY_GeV_per_cm,slopeY_err,interceptY_GeV,interceptY_err,chi2_ndf_Y,"
        << "entries_profileX,entries_profileY,total_clusters\n";

    // constants for branch sizes
    const int MAX_CLUST = 256; // should be >= maximum clusters in your skims

    // Loop runs
    for (int run : runs) {
        TString infile = Form("%sskim_run%d.root", skimDir.Data(), run);
        if (gSystem->AccessPathName(infile)) {
            logmsg(WARN, Form("Skipping run %d: file not found (%s)", run, infile.Data()));
            continue;
        }

        logmsg(INFO, Form("Processing run %d", run));
        TFile *f = TFile::Open(infile, "READ");
        if (!f || f->IsZombie()) { logmsg(ERROR, Form("Error opening file for run %d", run)); continue; }

        TTree *T = (TTree*)f->Get("T");
        if (!T) { logmsg(ERROR, Form("Tree 'T' not found in run %d", run)); f->Close(); continue; }

        // --- prepare branch variables
        // Common issue: branches can be int or double in different skims;
        // we read as double where appropriate and cast safely.
        double nclust_d = 0.0;
        double clusE_arr[MAX_CLUST]; memset(clusE_arr, 0, sizeof(clusE_arr));
        double clusX_arr[MAX_CLUST]; memset(clusX_arr, 0, sizeof(clusX_arr));
        double clusY_arr[MAX_CLUST]; memset(clusY_arr, 0, sizeof(clusY_arr));

        // Disable everything then enable only what we need for speed
        T->SetBranchStatus("*", 0);
        // Try multiple possible branch names robustly if needed (here we assume the names you listed)
        if (T->GetBranch("NPS.cal.nclust")) {
            T->SetBranchStatus("NPS.cal.nclust", 1);
            T->SetBranchAddress("NPS.cal.nclust", &nclust_d);
        } else if (T->GetBranch("NPS.cal.nclust_i")) {
            T->SetBranchStatus("NPS.cal.nclust_i", 1);
            T->SetBranchAddress("NPS.cal.nclust_i", &nclust_d);
        } else {
            // If missing, try reading other metadata or continue with a warning
            logmsg(WARN, Form("Run %d: branch NPS.cal.nclust not found — events may be skipped", run));
            // Still attempt to proceed; many skims store nclust implicitly by vector size — but TTree arrays are expected.
        }

        if (T->GetBranch("NPS.cal.clusE")) {
            T->SetBranchStatus("NPS.cal.clusE", 1);
            T->SetBranchAddress("NPS.cal.clusE", clusE_arr);
        } else { logmsg(ERROR, Form("Run %d: required branch NPS.cal.clusE not found — skipping run", run)); f->Close(); continue; }

        if (T->GetBranch("NPS.cal.clusX")) {
            T->SetBranchStatus("NPS.cal.clusX", 1);
            T->SetBranchAddress("NPS.cal.clusX", clusX_arr);
        } else { logmsg(WARN, Form("Run %d: branch NPS.cal.clusX not found", run)); }

        if (T->GetBranch("NPS.cal.clusY")) {
            T->SetBranchStatus("NPS.cal.clusY", 1);
            T->SetBranchAddress("NPS.cal.clusY", clusY_arr);
        } else { logmsg(WARN, Form("Run %d: branch NPS.cal.clusY not found", run)); }

        // --- per-run histos / profiles
        TH2D *h2EvsX = new TH2D(Form("h2EvsX_run%d", run),
                                Form("Run %d: E vs X;X (cm);E (GeV)", run), nXBins, xMin, xMax, nEBins, 0.0, 4.0);
        TH2D *h2EvsY = new TH2D(Form("h2EvsY_run%d", run),
                                Form("Run %d: E vs Y;Y (cm);E (GeV)", run), nYBins, yMin, yMax, nEBins, 0.0, 4.0);

        TProfile *pEvsX = new TProfile(Form("pEvsX_run%d", run),
                                       Form("Run %d: Mean E vs X;X (cm);Mean E (GeV)", run), nXBins, xMin, xMax);
        TProfile *pEvsY = new TProfile(Form("pEvsY_run%d", run),
                                       Form("Run %d: Mean E vs Y;Y (cm);Mean E (GeV)", run), nYBins, yMin, yMax);

        TH1D *hE_run = new TH1D(Form("hE_run%d", run), Form("Cluster E run %d;E (GeV);counts", run), nEBins, 0.0, 4.0);
        TH1D *hNclust_run = new TH1D(Form("hNclust_run%d", run), Form("Clusters per event run %d;Nclusters;events", run), 21, -0.5, 20.5);

        Long64_t nentries = T->GetEntries();
        if (nentries == 0) {
            logmsg(WARN, Form("No entries for run %d", run));
            // cleanup objects
            delete h2EvsX; delete h2EvsY; delete pEvsX; delete pEvsY; delete hE_run; delete hNclust_run;
            f->Close();
            continue;
        }

        Long64_t totalClustersThisRun = 0;
        // Event loop
        for (Long64_t ev=0; ev<nentries; ++ev) {
            T->GetEntry(ev);

            int nclust = (int) round(nclust_d); // robust cast
            if (nclust <= 0) {
                hNclust_run->Fill(0);
                continue;
            }
            hNclust_run->Fill(nclust);

            // select the highest-energy cluster in the event
            double ehi = -1.0;
            int idx_hi = -1;
            for (int c = 0; c < nclust && c < MAX_CLUST; ++c) {
                double e = clusE_arr[c];
                if (!isfinite(e)) continue;
                if (e > ehi) { ehi = e; idx_hi = c; }
            }
            if (idx_hi < 0) continue;

            double e = clusE_arr[idx_hi];
            double x = clusX_arr[idx_hi];
            double y = clusY_arr[idx_hi]; // <-- fixed: use Y branch here
            if (!isfinite(x) || !isfinite(y) || !isfinite(e)) continue;

            // Fill per-run histos and global accumulators
            h2EvsX->Fill(x, e);
            pEvsX->Fill(x, e);
            h2EvsX_global->Fill(x, e);
            pEvsX_global->Fill(x, e);

            h2EvsY->Fill(y, e);
            pEvsY->Fill(y, e);
            h2EvsY_global->Fill(y, e);
            pEvsY_global->Fill(y, e);

            hE_run->Fill(e);
            hE_global->Fill(e);

            ++totalClustersThisRun;
        } // end event loop

        // --- Fit profile X and Y with linear polynomial over bins with sufficient stats
        auto fitProfileLinear = [&](TProfile *prof, Int_t minEntriesPerBin = 5) {
            int nb = prof->GetNbinsX();
            int firstBin = 1, lastBin = nb;
            // find first/last bin with > minEntriesPerBin
            for (int b=1; b<=nb; ++b) {
                if (prof->GetBinEntries(b) > minEntriesPerBin) { firstBin = b; break; }
            }
            for (int b=nb; b>=1; --b) {
                if (prof->GetBinEntries(b) > minEntriesPerBin) { lastBin = b; break; }
            }
            double xmin = prof->GetXaxis()->GetBinCenter(firstBin) - prof->GetXaxis()->GetBinWidth(1)/2.0;
            double xmax = prof->GetXaxis()->GetBinCenter(lastBin)  + prof->GetXaxis()->GetBinWidth(1)/2.0;
            if (xmin >= xmax) { // fallback: full range
                xmin = prof->GetXaxis()->GetXmin();
                xmax = prof->GetXaxis()->GetXmax();
            }
            TF1 *f = new TF1("pol1_fit", "pol1", xmin, xmax);
            int fitStat = prof->Fit(f, "QR"); // quiet, use range
            double slope = f->GetParameter(1);
            double slope_err = f->GetParError(1);
            double intercept = f->GetParameter(0);
            double intercept_err = f->GetParError(0);
            double chi2ndf = (f->GetNDF()>0) ? f->GetChisquare()/f->GetNDF() : -1.0;
            return tuple<TF1*, double, double, double, double, double>(f, slope, slope_err, intercept, intercept_err, chi2ndf);
        };

        TF1 *fitX = nullptr;
        double slopeX=0, slopeXerr=0, interceptX=0, interceptXerr=0, chi2X=-1;
        tie(fitX, slopeX, slopeXerr, interceptX, interceptXerr, chi2X) = fitProfileLinear(pEvsX);

        TF1 *fitY = nullptr;
        double slopeY=0, slopeYerr=0, interceptY=0, interceptYerr=0, chi2Y=-1;
        tie(fitY, slopeY, slopeYerr, interceptY, interceptYerr, chi2Y) = fitProfileLinear(pEvsY);

        // --- Save per-run ROOT with objects
        TString runOutRoot = outDir + Form("response_run%d.root", run);
        TFile fout(runOutRoot, "RECREATE");
        h2EvsX->Write();
        pEvsX->Write();
        if (fitX) fitX->Write("fitX");
        h2EvsY->Write();
        pEvsY->Write();
        if (fitY) fitY->Write("fitY");
        hE_run->Write();
        hNclust_run->Write();
        fout.Close();

        // --- Save per-run PDFs (2 pages: colz + profile) for X and Y
        {
            TCanvas c1("c1", "E vs X", 900, 700);
            h2EvsX->Draw("COLZ");
            c1.SaveAs(outDir + Form("EvsX_run%d_colz.pdf", run));

            TCanvas c2("c2", "Profile X", 900, 700);
            pEvsX->Draw();
            if (fitX) { fitX->SetLineColor(kRed); fitX->Draw("same"); }
            c2.SaveAs(outDir + Form("EvsX_run%d_profile.pdf", run));
        }
        {
            TCanvas c3("c3", "E vs Y", 900, 700);
            h2EvsY->Draw("COLZ");
            c3.SaveAs(outDir + Form("EvsY_run%d_colz.pdf", run));

            TCanvas c4("c4", "Profile Y", 900, 700);
            pEvsY->Draw();
            if (fitY) { fitY->SetLineColor(kRed); fitY->Draw("same"); }
            c4.SaveAs(outDir + Form("EvsY_run%d_profile.pdf", run));
        }

        // --- Log & CSV (concise)
        logmsg(INFO, Form("Run %d: clusters=%lld; slopeX=%.6e±%.6e GeV/cm chi2/ndf=%.2f; slopeY=%.6e±%.6e GeV/cm chi2/ndf=%.2f",
                          run, (long long)totalClustersThisRun,
                          slopeX, slopeXerr, chi2X,
                          slopeY, slopeYerr, chi2Y));
        csv << run << ","
            << slopeX << "," << slopeXerr << "," << interceptX << "," << interceptXerr << "," << chi2X << ","
            << slopeY << "," << slopeYerr << "," << interceptY << "," << interceptYerr << "," << chi2Y << ","
            << (long long)pEvsX->GetEntries() << "," << (long long)pEvsY->GetEntries() << "," << totalClustersThisRun << "\n";

        // cleanup per-run heap objects
        if (fitX) delete fitX;
        if (fitY) delete fitY;
        delete h2EvsX;
        delete h2EvsY;
        delete pEvsX;
        delete pEvsY;
        delete hE_run;
        delete hNclust_run;
        f->Close();
        delete f;
    } // end runs loop

    // --- Global fit for X and Y profiles
    TF1 *fitGlobalX = new TF1("fitGlobalX", "pol1", pEvsX_global->GetXaxis()->GetXmin(), pEvsX_global->GetXaxis()->GetXmax());
    pEvsX_global->Fit(fitGlobalX, "QR");
    double slopeGX = fitGlobalX->GetParameter(1), slopeGXerr = fitGlobalX->GetParError(1);
    double interceptGX = fitGlobalX->GetParameter(0), interceptGXerr = fitGlobalX->GetParError(0);
    double chi2ndfGX = (fitGlobalX->GetNDF()>0) ? fitGlobalX->GetChisquare()/fitGlobalX->GetNDF() : -1.0;

    TF1 *fitGlobalY = new TF1("fitGlobalY", "pol1", pEvsY_global->GetXaxis()->GetXmin(), pEvsY_global->GetXaxis()->GetXmax());
    pEvsY_global->Fit(fitGlobalY, "QR");
    double slopeGY = fitGlobalY->GetParameter(1), slopeGYerr = fitGlobalY->GetParError(1);
    double interceptGY = fitGlobalY->GetParameter(0), interceptGYerr = fitGlobalY->GetParError(0);
    double chi2ndfGY = (fitGlobalY->GetNDF()>0) ? fitGlobalY->GetChisquare()/fitGlobalY->GetNDF() : -1.0;

    // Save global ROOT and PDF plots
    TString globalRoot = outDir + "response_global.root";
    TFile foutG(globalRoot, "RECREATE");
    h2EvsX_global->Write();
    pEvsX_global->Write();
    fitGlobalX->Write("fitGlobalX");
    h2EvsY_global->Write();
    pEvsY_global->Write();
    fitGlobalY->Write("fitGlobalY");
    hE_global->Write();
    hNclust_global->Write();
    foutG.Close();

    // Global PDFs
    {
        TCanvas cG1("cG1", "E vs X (global)", 1000, 800);
        h2EvsX_global->Draw("COLZ");
        cG1.SaveAs(outDir + "EvsX_global_colz.pdf");

        TCanvas cG2("cG2", "Profile X (global)", 1000, 800);
        pEvsX_global->Draw();
        fitGlobalX->SetLineColor(kRed); fitGlobalX->Draw("same");
        cG2.SaveAs(outDir + "EvsX_global_profile.pdf");
    }
    {
        TCanvas cG3("cG3", "E vs Y (global)", 1000, 800);
        h2EvsY_global->Draw("COLZ");
        cG3.SaveAs(outDir + "EvsY_global_colz.pdf");

        TCanvas cG4("cG4", "Profile Y (global)", 1000, 800);
        pEvsY_global->Draw();
        fitGlobalY->SetLineColor(kRed); fitGlobalY->Draw("same");
        cG4.SaveAs(outDir + "EvsY_global_profile.pdf");
    }

    // Final logging
    logmsg(INFO, Form("Global X fit: slope=%.6e ± %.6e  intercept=%.3f ± %.3f  chi2/ndf=%.2f",
                      slopeGX, slopeGXerr, interceptGX, interceptGXerr, chi2ndfGX));
    logmsg(INFO, Form("Global Y fit: slope=%.6e ± %.6e  intercept=%.3f ± %.3f  chi2/ndf=%.2f",
                      slopeGY, slopeGYerr, interceptGY, interceptGYerr, chi2ndfGY));

    csv.close();

    // cleanup global objects
    delete fitGlobalX;
    delete fitGlobalY;
    delete h2EvsX_global;
    delete pEvsX_global;
    delete h2EvsY_global;
    delete pEvsY_global;
    delete hE_global;
    delete hNclust_global;

    sw.Stop();
    logmsg(INFO, Form("Done. Time: %.2f s. Summary CSV: %s", sw.RealTime(), csvFile.Data()));
}
