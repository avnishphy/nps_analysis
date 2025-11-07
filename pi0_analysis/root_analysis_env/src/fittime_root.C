// fittime_root.C
// Enhanced NPS fittime macro
// - reads block_by_block.txt (from pi0_analysis_peter.C)
// - supports two offset methods (mean and Gaussian-fit)
// - produces layout maps (30x36), comparisons, overlayed timing spectra,
//   records fit failures and writes .out/.param/.summary/.root/.pdf
//
// Author: adapted for user's workflow (Avnish + ChatGPT), 2025-11-05
// ----------------------------------------------------------------------------

#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLegend.h>
#include <TText.h>
#include <TSystem.h>

#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;

void fittime_root(const TString &blockFile = "output/block_by_block.txt",
                  const TString &outDir   = "output/fittime/",
                  Int_t nBlocks = 1080,
                  Int_t nBins   = 20,
                  Int_t nRows   = 30,      // NPS layout: 30 rows
                  Int_t nCols   = 36,      // NPS layout: 36 columns (30x36=1080)
                  Int_t minCountsForOffset = 10, // min counts required to compute offsets
                  Double_t idealTime = 150.0     // target time [ns]
                  )
{
    // --- style
    gStyle->SetOptStat(0);
    cout << "\n=== fittime_root (enhanced) ===\n";
    cout << "Input: " << blockFile << "\nOutput dir: " << outDir << "\n";

    // Ensure output directory exists
    TString outdir = outDir;
    if (!outdir.EndsWith("/")) outdir += "/";
    gSystem->mkdir(outdir, kTRUE);

    // Output filenames
    TString outRootName = outdir + "fittime_results.root";
    TString outPdfName  = outdir + "fittime_offsets.pdf";
    TString outTxtName  = outdir + "fittime.out";
    TString paramName   = outdir + "tdc_offset.param";
    TString summaryName = outdir + "fittime_summary.txt";

    // --- read input file
    ifstream fin(blockFile.Data());
    if (!fin.is_open()) {
        cerr << "ERROR: cannot open " << blockFile << "\n";
        return;
    }
    cout << "Reading block file and extracting timing histograms..." << endl;

    // containers
    vector<vector<int>> blktimeh(nBlocks, vector<int>(nBins, 0));
    vector<int> ntot(nBlocks, 0);

    // skip first 2*nBlocks lines (eloblkhist + ehiblkhist as produced by pi0_analysis_peter.C)
    string line;
    int skipLines = 2 * nBlocks;
    int lineCount = 0;
    while (lineCount < skipLines && getline(fin, line)) ++lineCount;
    if (lineCount < skipLines) {
        cerr << "Warning: file shorter than expected while skipping first sections (" << lineCount << "/" << skipLines << ")\n";
    }

    // read next nBlocks lines (timing histograms)
    int readBlocks = 0;
    while (readBlocks < nBlocks && getline(fin, line)) {
        if (line.find_first_not_of(" \t\r\n") == string::npos) continue;
        istringstream iss(line);
        int idx;
        if (!(iss >> idx)) continue;
        int ib = idx;
        // allow 1-based index mapping
        if (ib < 0 || ib >= nBlocks) {
            if (idx >= 1 && idx <= nBlocks) ib = idx - 1;
            else ib = readBlocks; // fallback
        }
        for (int j = 0; j < nBins; ++j) {
            int v=0;
            if (!(iss >> v)) v = 0;
            blktimeh[ib][j] = v;
            ntot[ib] += v;
        }
        ++readBlocks;
    }
    fin.close();

    cout << "Read " << readBlocks << " timing histograms (expected " << nBlocks << ").\n";

    // --- compute offsets by two methods
    vector<double> mean_time(nBlocks, 0.0);
    vector<double> offset_mean(nBlocks, 0.0);
    vector<double> offset_fit(nBlocks, 0.0);
    vector<int> fitStatus(nBlocks, 0);       // 0 = not attempted/skipped, 1 = success, -1 = failed
    vector<double> fit_mean(nBlocks, 0.0);
    vector<double> fit_sigma(nBlocks, 0.0);
    vector<double> fit_chi2(nBlocks, 0.0);

    // absolute center for bin k (0-based): center = 140.0 + k + 0.5
    auto binCenterAbs = [](int k)->double { return 140.0 + (double)k + 0.5; };

    int npos = 0, nneg = 0;
    for (int ib = 0; ib < nBlocks; ++ib) {
        long long sumN = 0;
        double sumT = 0.0;
        for (int k = 0; k < nBins; ++k) {
            int cnt = blktimeh[ib][k];
            sumN += cnt;
            sumT += cnt * binCenterAbs(k);
        }
        if (sumN > 0) mean_time[ib] = sumT / (double)sumN;
        else mean_time[ib] = 0.0;
        if (sumN >= minCountsForOffset) {
            offset_mean[ib] = idealTime - mean_time[ib]; // offset to shift mean_time -> idealTime
        } else {
            offset_mean[ib] = 0.0;
        }
        if (offset_mean[ib] > 0) ++npos;
        if (offset_mean[ib] < 0) ++nneg;
    }

    // --- Fit method: fit gaussian to each block's timing histogram (if enough counts)
    // We'll reuse a temporary TH1D for fitting to avoid creating 1080 histos.
    TH1D *htemp = new TH1D("htemp","tmp", nBins, 140.0, 140.0 + nBins); // bins [140,140+20)
    // Note: bin i covers [140+i, 141+i), center at 140+i+0.5 as intended.

    for (int ib = 0; ib < nBlocks; ++ib) {
        // skip low-stat blocks
        if (ntot[ib] < minCountsForOffset) {
            fitStatus[ib] = 0; // skipped due to low stats
            continue;
        }
        // fill temporary hist
        htemp->Reset();
        for (int k = 0; k < nBins; ++k) {
            int cnt = blktimeh[ib][k];
            // Fill bin-by-bin: fill content directly to avoid Poisson fluctuations
            htemp->SetBinContent(k+1, cnt); // TH1 bins are 1-based
            htemp->SetBinError(k+1, sqrt((double)cnt));
        }
        // find max bin to choose fit window
        int ibmax = htemp->GetMaximumBin();
        double center_guess = htemp->GetBinCenter(ibmax);
        double amp_guess = htemp->GetBinContent(ibmax);
        double sigma_guess = 1.0; // ns, sensible initial

        // define fit range: center +/- 3 bins (~3 ns)
        double fit_lo = max(140.0, center_guess - 3.0);
        double fit_hi = min(140.0 + nBins, center_guess + 3.0);

        TF1 *fgaus = new TF1("fgaus","gaus", fit_lo, fit_hi);
        fgaus->SetParameters(amp_guess, center_guess, sigma_guess);
        // set reasonable parameter limits
        fgaus->SetParLimits(2, 0.2, 5.0); // sigma between 0.2 and 5 ns

        // perform fit quietly ("Q") and robust ("R") without drawing ("N")
        int fitRes = htemp->Fit(fgaus, "QRN"); // returns 0 on success typically

        double mean_f = fgaus->GetParameter(1);
        double sigma_f = fgaus->GetParameter(2);
        double amp_f = fgaus->GetParameter(0);
        double mean_err = fgaus->GetParError(1);
        double sigma_err = fgaus->GetParError(2);
        double chi2 = fgaus->GetChisquare();
        int ndf = fgaus->GetNDF();

        // decide success criteria: fitRes==0 and mean within [140,160] and sigma reasonable
        bool success = (fitRes == 0) && (mean_f >= 140.0-2.0) && (mean_f <= 160.0+2.0) && (sigma_f > 0.0 && sigma_f < 10.0);

        if (success) {
            fitStatus[ib] = 1;
            fit_mean[ib] = mean_f;
            fit_sigma[ib] = sigma_f;
            fit_chi2[ib] = (ndf>0) ? chi2/ndf : chi2;
            offset_fit[ib] = idealTime - mean_f;
        } else {
            fitStatus[ib] = -1;  // failed fit
            fit_mean[ib] = mean_f;
            fit_sigma[ib] = sigma_f;
            fit_chi2[ib] = (ndf>0) ? chi2/ndf : chi2;
            offset_fit[ib] = 0.0; // do not apply fit offset by default
        }
        delete fgaus;
    }
    delete htemp;

    // collect fit failure indices
    vector<int> fitFailures;
    for (int ib = 0; ib < nBlocks; ++ib) if (fitStatus[ib] == -1) fitFailures.push_back(ib);

    // --- Build global timing histograms (before, mean-corrected, fit-corrected)
    double tmin = 140.0, tmax = 160.0;
    int tBins = 100; // fine binning for overlay plots
    TH1D *hBefore = new TH1D("hBefore", "NPS cluster timing (no offset);Time (ns);Counts", tBins, tmin, tmax);
    TH1D *hMeanAfter = new TH1D("hMeanAfter", "NPS cluster timing (mean offset applied);Time (ns);Counts", tBins, tmin, tmax);
    TH1D *hFitAfter  = new TH1D("hFitAfter",  "NPS cluster timing (fit offset applied);Time (ns);Counts", tBins, tmin, tmax);

    // Also produce histogram contributions for failed fit blocks (before & mean-corrected)
    TH1D *hFailedBefore = new TH1D("hFailedBefore", "Timing from fit-failed blocks (before);Time (ns);Counts", tBins, tmin, tmax);
    TH1D *hFailedMeanAfter = new TH1D("hFailedMeanAfter", "Timing from fit-failed blocks (mean-applied);Time (ns);Counts", tBins, tmin, tmax);

    // Fill histograms by summing block contributions at absolute bin centers
    for (int ib = 0; ib < nBlocks; ++ib) {
        for (int k = 0; k < nBins; ++k) {
            int cnt = blktimeh[ib][k];
            if (cnt <= 0) continue;
            double center_abs = binCenterAbs(k); // 140 + k + 0.5

            // before
            hBefore->Fill(center_abs, (double)cnt);

            // mean method
            double corrected_mean = center_abs + offset_mean[ib];
            if (corrected_mean >= tmin && corrected_mean <= tmax) hMeanAfter->Fill(corrected_mean, (double)cnt);

            // fit method (only apply if fit succeeded)
            if (fitStatus[ib] == 1) {
                double corrected_fit = center_abs + offset_fit[ib];
                if (corrected_fit >= tmin && corrected_fit <= tmax) hFitAfter->Fill(corrected_fit, (double)cnt);
            } else {
                // record failed-block contributions separately
                hFailedBefore->Fill(center_abs, (double)cnt);
                double corrected_mean_failed = center_abs + offset_mean[ib];
                if (corrected_mean_failed >= tmin && corrected_mean_failed <= tmax)
                    hFailedMeanAfter->Fill(corrected_mean_failed, (double)cnt);
            }
        }
    }

    // --- Create block-layout histograms (30x36) for mean, fit and diff
    TH2D *hLayoutMean = new TH2D("hLayoutMean", "NPS layout: offset (mean method);col;row;offset (ns)",
                                nCols, 0.5, nCols+0.5, nRows, 0.5, nRows+0.5);
    TH2D *hLayoutFit  = new TH2D("hLayoutFit",  "NPS layout: offset (fit method);col;row;offset (ns)",
                                nCols, 0.5, nCols+0.5, nRows, 0.5, nRows+0.5);
    TH2D *hLayoutDiff = new TH2D("hLayoutDiff", "NPS layout: fit - mean (ns);col;row;diff (ns)",
                                nCols, 0.5, nCols+0.5, nRows, 0.5, nRows+0.5);

    // Fill layout maps; mapping chosen: block index -> row = ib / nCols, col = ib % nCols
    for (int ib = 0; ib < nBlocks; ++ib) {
        int row = ib / nCols; // 0..nRows-1
        int col = ib % nCols; // 0..nCols-1
        int binx = col + 1;
        int biny = row + 1;
        hLayoutMean->SetBinContent(binx, biny, offset_mean[ib]);
        hLayoutFit->SetBinContent(binx, biny, offset_fit[ib]);
        hLayoutDiff->SetBinContent(binx, biny, (offset_fit[ib] - offset_mean[ib]));
    }

    // --- Offset distributions and comparisons
    TH1D *hOffsetMean = new TH1D("hOffsetMean", "Offset distribution (mean method);offset (ns);blocks", 120, -12, 12);
    TH1D *hOffsetFit  = new TH1D("hOffsetFit",  "Offset distribution (fit method);offset (ns);blocks", 120, -12, 12);
    TH1D *hOffsetDiff = new TH1D("hOffsetDiff", "Offset difference (fit - mean);diff (ns);blocks", 120, -6, 6);

    for (int ib = 0; ib < nBlocks; ++ib) {
        if (ntot[ib] >= minCountsForOffset) {
            hOffsetMean->Fill(offset_mean[ib]);
            if (fitStatus[ib] == 1) hOffsetFit->Fill(offset_fit[ib]);
            if (fitStatus[ib] == 1) hOffsetDiff->Fill(offset_fit[ib] - offset_mean[ib]);
        }
    }

    // Graph mean vs fit for blocks where fit succeeded
    vector<double> vMean;
    vector<double> vFit;
    vector<int>    vIdx;
    for (int ib = 0; ib < nBlocks; ++ib) {
        if (fitStatus[ib] == 1 && ntot[ib] >= minCountsForOffset) {
            vMean.push_back(offset_mean[ib]);
            vFit.push_back(offset_fit[ib]);
            vIdx.push_back(ib);
        }
    }
    TGraph *gMeanVsFit = nullptr;
    if (!vMean.empty()) {
        gMeanVsFit = new TGraph(vMean.size());
        for (size_t i = 0; i < vMean.size(); ++i) gMeanVsFit->SetPoint(i, vMean[i], vFit[i]);
        gMeanVsFit->SetTitle("Mean vs Fit offsets;offset_mean (ns);offset_fit (ns)");
    }

    // --- Prepare canvases and save multi-page PDF
    // Page order:
    //  1) timing heatmap (block vs time)
    //  2) layout mean
    //  3) layout fit
    //  4) layout diff
    //  5) offset histograms (mean, fit)
    //  6) mean vs fit scatter + diff hist
    //  7) Timing overlays (no offset, mean, fit)
    //  8) failed fit summary (hist + small table)
    // Save as multi-page PDF (use "(" and ")" convention)

    // heatmap (block vs time) - similar to previous, but use absolute time axis
    TH2D *hCounts2D = new TH2D("hCounts2D", "Timing histogram per block;Block;Time (ns)",
                               nBlocks, -0.5, nBlocks - 0.5, nBins, 140.0 - 0.5, 140.0 + nBins + 0.5);

    for (int ib = 0; ib < nBlocks; ++ib) {
        for (int k = 0; k < nBins; ++k) {
            int cnt = blktimeh[ib][k];
            if (cnt > 0) {
                double center_abs = binCenterAbs(k);
                hCounts2D->Fill(ib, center_abs, (double)cnt);
            }
        }
    }

    // Canvas 1: overview heatmap + offsets graph + offset hist
    TCanvas *c1 = new TCanvas("c1", "fittime overview", 1400, 900);
    c1->Divide(2,2);

    c1->cd(1);
    hCounts2D->SetTitle("Timing histogram per block;block index;time (ns)");
    hCounts2D->GetZaxis()->SetTitle("counts");
    hCounts2D->Draw("COLZ");
    gPad->SetLogz(1);

    c1->cd(2);
    hOffsetMean->SetLineColor(kBlue); hOffsetMean->Draw();
    hOffsetFit->SetLineColor(kRed); hOffsetFit->SetLineWidth(2); hOffsetFit->Draw("SAME");
    TLegend *leg = new TLegend(0.6,0.7,0.88,0.88);
    leg->AddEntry(hOffsetMean,"Mean offsets","l");
    leg->AddEntry(hOffsetFit,"Fit offsets (successful fits)","l");
    leg->SetBorderSize(0); leg->Draw();

    c1->cd(3);
    if (gMeanVsFit) {
        gMeanVsFit->SetMarkerStyle(20);
        gMeanVsFit->Draw("AP");
        // draw diagonal
        TLine *diag = new TLine(-12,-12,12,12); diag->SetLineStyle(2); diag->SetLineColor(kGray+2); diag->Draw("same");
    } else {
        TH1D *hEmpty = new TH1D("hEmpty","No fit-successful blocks",10,0,1);
        hEmpty->Draw();
    }

    c1->cd(4);
    hOffsetDiff->Draw();

    // Canvas 2: layout maps
    TCanvas *c2 = new TCanvas("c2","layout maps", 1400, 900);
    c2->Divide(3,1);

    c2->cd(1);
    hLayoutMean->SetContour(99);
    hLayoutMean->Draw("COLZ");
    TText t1(0.15, 0.92, Form("Mean offset (target %.1f ns)", idealTime));
    t1.SetNDC(); t1.Draw();

    c2->cd(2);
    hLayoutFit->SetContour(99);
    hLayoutFit->Draw("COLZ");
    TText t2(0.15, 0.92, "Fit offset (successful fits only)");
    t2.SetNDC(); t2.Draw();

    c2->cd(3);
    hLayoutDiff->SetContour(99);
    hLayoutDiff->Draw("COLZ");
    TText t3(0.15, 0.92, "Fit - Mean (ns)");
    t3.SetNDC(); t3.Draw();

    // Canvas 3: timing overlay before / mean / fit
    TCanvas *c3 = new TCanvas("c3","timing overlay", 1200, 600);
    hBefore->SetLineColor(kBlack); hBefore->SetLineWidth(1);
    hMeanAfter->SetLineColor(kBlue); hMeanAfter->SetLineWidth(2);
    hFitAfter->SetLineColor(kRed); hFitAfter->SetLineWidth(2);

    hBefore->Draw();
    hMeanAfter->Draw("SAME");
    hFitAfter->Draw("SAME");

    TLegend *legTim = new TLegend(0.6,0.6,0.92,0.85);
    legTim->AddEntry(hBefore,"Raw (no offset)","l");
    legTim->AddEntry(hMeanAfter,"Mean-offset applied","l");
    legTim->AddEntry(hFitAfter,"Fit-offset applied (fit success only)","l");
    legTim->SetBorderSize(0);
    legTim->Draw();

    // Canvas 4: failed fit diagnostics
    TCanvas *c4 = new TCanvas("c4","fit failures", 1000, 600);
    c4->Divide(2,1);

    c4->cd(1);
    TH1D *hFailCount = new TH1D("hFailCount", "Counts in fit-failed blocks (before);Time (ns);Counts", tBins, tmin, tmax);
    hFailCount->Add(hFailedBefore);
    hFailCount->Draw();

    c4->cd(2);
    TH1D *hFailMeanAfter = new TH1D("hFailMeanAfter", "Fit-failed blocks, mean-shifted;Time (ns);Counts", tBins, tmin, tmax);
    hFailMeanAfter->Add(hFailedMeanAfter);
    hFailMeanAfter->Draw();

    // Save multi-page PDF with ordered pages
    cout << "Writing PDF pages to " << outPdfName << " ..." << endl;

    // Force updates before saving
    c1->Update();
    c2->Update();
    c3->Update();
    c4->Update();

    // Page 1 (start the PDF)
    c1->SaveAs((outPdfName + "(").Data());  

    // Page 2
    c2->SaveAs(outPdfName.Data());

    // Page 3
    c3->SaveAs(outPdfName.Data());

    // Page 4
    c4->SaveAs(outPdfName.Data());

    // Close the PDF stream
    TCanvas *cend = new TCanvas("cend","end",200,200);
    cend->Update();
    cend->SaveAs((outPdfName + ")").Data());
    delete cend;

    cout << "PDF written successfully to " << outPdfName << endl;


    // --- Write ROOT file with objects
    cout << "Writing ROOT file: " << outRootName << "\n";
    TFile froot(outRootName, "RECREATE");
    hCounts2D->Write();
    hLayoutMean->Write();
    hLayoutFit->Write();
    hLayoutDiff->Write();
    hOffsetMean->Write();
    hOffsetFit->Write();
    hOffsetDiff->Write();
    if (gMeanVsFit) gMeanVsFit->Write("gMeanVsFit");
    hBefore->Write();
    hMeanAfter->Write();
    hFitAfter->Write();
    hFailedBefore->Write();
    hFailedMeanAfter->Write();
    hFailCount->Write();
    hFailMeanAfter->Write();
    froot.Close();

    // --- Write text outputs: fittime.out, tdc_offset.param, fittime_summary.txt
    ofstream fout(outTxtName.Data());
    ofstream fparam(paramName.Data());
    ofstream fsummary(summaryName.Data());
    if (!fout.is_open() || !fparam.is_open() || !fsummary.is_open()) {
        cerr << "ERROR: cannot open one of the text outputs for writing.\n";
    } else {
        fout << "# block ntot mean_time offset_mean offset_fit fitStatus(1=ok,-1=bad,0=skipped) fit_mean fit_sigma fit_chi2\n";
        for (int ib = 0; ib < nBlocks; ++ib) {
            fout << setw(4) << ib
                 << setw(8) << ntot[ib]
                 << setw(10) << fixed << setprecision(2) << mean_time[ib]
                 << setw(10) << fixed << setprecision(3) << offset_mean[ib]
                 << setw(10) << fixed << setprecision(3) << offset_fit[ib]
                 << setw(4) << fitStatus[ib]
                 << setw(10) << fixed << setprecision(3) << fit_mean[ib]
                 << setw(10) << fixed << setprecision(3) << fit_sigma[ib]
                 << setw(10) << fixed << setprecision(3) << fit_chi2[ib]
                 << "\n";
        }
        fout.close();

        // tdc_offset.param: write mean offsets by default (you can choose fit offsets if you prefer)
        // format: 108 rows x 10 cols -> 1080 numbers total, f5.2 format approximated
        const int colsPerRow = 10;
        int nRows = (nBlocks + colsPerRow - 1) / colsPerRow;
        for (int r = 0; r < nRows; ++r) {
            int k1 = r * colsPerRow;
            int k2 = std::min(k1 + colsPerRow - 1, nBlocks - 1);
            for (int k = k1; k <= k2; ++k) {
                double val = offset_mean[k]; // writing mean-method offsets by default
                fparam << fixed << setprecision(2) << setw(5) << val;
                if (k != k2) fparam << ",";
            }
            if (k2 - k1 + 1 < colsPerRow) {
                for (int p = 0; p < colsPerRow - (k2 - k1 + 1); ++p) fparam << ", 0.00";
            }
            fparam << "\n";
        }
        fparam.close();

        // summary
        fsummary << "# fittime summary\n";
        fsummary << "# input file: " << blockFile << "\n";
        fsummary << "# blocks: " << nBlocks << "  bins: " << nBins << "\n";
        fsummary << "# minCountsForOffset: " << minCountsForOffset << "  idealTime: " << idealTime << "\n";
        fsummary << "# total blocks read: " << readBlocks << "\n";
        fsummary << "# fit failures: " << fitFailures.size() << "\n";
        if (!fitFailures.empty()) {
            fsummary << "# failed block indices (count = " << fitFailures.size() << "):\n";
            for (size_t i = 0; i < fitFailures.size(); ++i) {
                fsummary << fitFailures[i] << ((i%20==19) ? "\n" : " ");
            }
            fsummary << "\n";
        }
        fsummary.close();
    }

    cout << "Wrote outputs in " << outdir << "\n";
    cout << "Fit failures: " << fitFailures.size() << " blocks (listed in " << summaryName << ")\n";
    cout << "Done.\n";

    // clean up (delete canvases and graphs to free memory)
    delete c1;
    delete c2;
    delete c3;
    delete c4;
    delete hCounts2D;
    delete hLayoutMean;
    delete hLayoutFit;
    delete hLayoutDiff;
    delete hOffsetMean;
    delete hOffsetFit;
    delete hOffsetDiff;
    if (gMeanVsFit) delete gMeanVsFit;
    delete hBefore;
    delete hMeanAfter;
    delete hFitAfter;
    delete hFailedBefore;
    delete hFailedMeanAfter;
    delete hFailCount;
    delete hFailMeanAfter;
}
