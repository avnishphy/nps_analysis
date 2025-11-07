// fittime_root.C
// Enhanced + optimized NPS fittime macro
// - reads block_by_block.txt (from pi0_analysis_peter.C)
// - supports two offset methods (mean and Gaussian-fit)
// - produces 30x36 layout maps (block numbers shown), comparisons,
//   overlayed timing spectra, records fit failures and writes .out/.param/.summary/.root/.pdf
//
// Author: adapted for user's workflow (Avnish + ChatGPT), 2025-11-05
// Improvements:
//  - block numbers printed on layout (lower-left origin, left->right then up)
//  - colorbar shown and margins adjusted for COLZ plots
//  - minimal but informative console output only
//  - all 1D histos drawn as lines ("HIST") so overlays are clear
//  - multi-page PDF saving fixed and robust
//  - writes both mean-based and fit-based param files
// ----------------------------------------------------------------------------

#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLine.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TBox.h>
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

void fittime_root1(const TString &blockFile = "output/block_by_block.txt",
                  const TString &outDir   = "output/fittime/",
                  Int_t nBlocks = 1080,
                  Int_t nBins   = 20,
                  Int_t nRows   = 30,      // NPS layout: 30 rows (vertical)
                  Int_t nCols   = 36,      // NPS layout: 36 columns (horizontal)
                  Int_t minCountsForOffset = 10, // min counts required to compute offsets
                  Double_t idealTime = 150.0     // target time [ns]
                  )
{
    // --- style / global settings
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis); // nicer palette
    gStyle->SetTitleOffset(1.1, "Z");
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
    TString paramMeanName = outdir + "tdc_offset.param";        // mean-based offsets
    TString paramFitName  = outdir + "tdc_offset_fit.param";   // fit-based offsets
    TString summaryName = outdir + "fittime_summary.txt";

    // --- read input file (block_by_block)
    ifstream fin(blockFile.Data());
    if (!fin.is_open()) {
        cerr << "ERROR: cannot open " << blockFile << "\n";
        return;
    }

    cout << "Reading '" << blockFile << "' and extracting timing histograms...\n";

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
        // allow 1-based index mapping if present
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

    // --- helper: absolute center for bin k (0-based): center = 140.0 + k + 0.5
    auto binCenterAbs = [](int k)->double { return 140.0 + double(k) + 0.5; };

    // --- compute mean-based offsets
    vector<double> mean_time(nBlocks, 0.0);
    vector<double> offset_mean(nBlocks, 0.0);
    int npos_mean = 0, nneg_mean = 0;
    for (int ib=0; ib<nBlocks; ++ib) {
        long long sumN = 0;
        double sumT = 0.0;
        for (int k=0; k<nBins; ++k) {
            const int cnt = blktimeh[ib][k];
            sumN += cnt;
            sumT += cnt * binCenterAbs(k);
        }
        mean_time[ib] = (sumN>0) ? (sumT / (double)sumN) : 0.0;
        if (sumN >= minCountsForOffset) offset_mean[ib] = idealTime - mean_time[ib];
        else offset_mean[ib] = 0.0;
        if (offset_mean[ib] > 0) ++npos_mean;
        if (offset_mean[ib] < 0) ++nneg_mean;
    }

    // --- fit-based offsets (Gaussian fit) - one TH1 reused for speed
    vector<double> offset_fit(nBlocks, 0.0);
    vector<int> fitStatus(nBlocks, 0);      // 0 skipped, 1 success, -1 failed
    vector<double> fit_mean(nBlocks, 0.0), fit_sigma(nBlocks, 0.0), fit_chi2(nBlocks, 0.0);
    TH1D *htmp = new TH1D("htmp","tmp bin hist", nBins, 140.0, 140.0 + nBins); // edges [140,160)
    htmp->SetDirectory(nullptr);

    cout << "Fitting blocks (fit attempted only when ntot >= " << minCountsForOffset << ")...\n";
    int fitAttempts = 0, fitSuccess = 0, fitFails = 0;
    for (int ib=0; ib<nBlocks; ++ib) {
        if (ntot[ib] < minCountsForOffset) {
            fitStatus[ib] = 0; // skip
            continue;
        }
        // fill histogram by direct bin content (no Poisson refill)
        htmp->Reset();
        for (int k=0; k<nBins; ++k) {
            int cnt = blktimeh[ib][k];
            htmp->SetBinContent(k+1, cnt);
            htmp->SetBinError(k+1, (cnt>0) ? sqrt((double)cnt) : 1.0);
        }

        // auto choose fit window around maximum
        int ibmax = htmp->GetMaximumBin();
        double center_guess = htmp->GetBinCenter(ibmax);
        double amp_guess = htmp->GetBinContent(ibmax);
        double sigma_guess = 1.0;

        double fit_lo = max(140.0, center_guess - 3.0);
        double fit_hi = min(140.0 + nBins, center_guess + 3.0);

        TF1 *fgaus = new TF1("fgaus","gaus", fit_lo, fit_hi);
        fgaus->SetParameters(amp_guess, center_guess, sigma_guess);
        fgaus->SetParLimits(2, 0.2, 5.0);

        // perform fit quietly
        ++fitAttempts;
        int fitRes = htmp->Fit(fgaus, "QRN"); // quiet, robust, no draw

        double mean_f = fgaus->GetParameter(1);
        double sigma_f = fgaus->GetParameter(2);
        double chi2 = fgaus->GetChisquare();
        int ndf = fgaus->GetNDF();

        bool success = (fitRes == 0) && (mean_f >= 138.0) && (mean_f <= 162.0) && (sigma_f > 0.0 && sigma_f < 10.0);
        if (success) {
            fitStatus[ib] = 1;
            fit_mean[ib] = mean_f;
            fit_sigma[ib] = sigma_f;
            fit_chi2[ib] = (ndf>0) ? (chi2/ndf) : chi2;
            offset_fit[ib] = idealTime - mean_f;
            ++fitSuccess;
        } else {
            fitStatus[ib] = -1; // fit failed
            fit_mean[ib] = mean_f;
            fit_sigma[ib] = sigma_f;
            fit_chi2[ib] = (ndf>0) ? (chi2/ndf) : chi2;
            offset_fit[ib] = 0.0;
            ++fitFails;
        }
        delete fgaus;
        // modest progress printing
        if ((ib%120)==0) cout << "  fitted block " << ib << " / " << nBlocks << "\n";
    }
    delete htmp;

    // collect fit failure indices
    vector<int> fitFailures;
    for (int ib=0; ib<nBlocks; ++ib) if (fitStatus[ib] == -1) fitFailures.push_back(ib);

    // --- Build global timing histograms (before, mean-corrected, fit-corrected)
    double tmin = 140.0, tmax = 160.0;
    int tBins = 20;
    TH1D *hBefore = new TH1D("hBefore", "NPS cluster timing (no offset);Time (ns);Counts", tBins, tmin, tmax);
    TH1D *hMeanAfter = new TH1D("hMeanAfter", "NPS cluster timing (mean offset applied);Time (ns);Counts", tBins, tmin, tmax);
    TH1D *hFitAfter  = new TH1D("hFitAfter",  "NPS cluster timing (fit offset applied);Time (ns);Counts", tBins, tmin, tmax);
    TH1D *hFailedBefore = new TH1D("hFailedBefore", "Timing from fit-failed blocks (before);Time (ns);Counts", tBins, tmin, tmax);
    TH1D *hFailedMeanAfter = new TH1D("hFailedMeanAfter", "Timing from fit-failed blocks (mean-shifted);Time (ns);Counts", tBins, tmin, tmax);

    for (int ib=0; ib<nBlocks; ++ib) {
        for (int k=0; k<nBins; ++k) {
            int cnt = blktimeh[ib][k];
            if (cnt == 0) continue;

            double binWidthInput = (tmax - tmin) / nBins;
            double center_abs = tmin + (k + 0.5) * binWidthInput;

            hBefore->Fill(center_abs, cnt);

            double cm = center_abs + offset_mean[ib];
            if (cm >= tmin && cm <= tmax) hMeanAfter->Fill(cm, cnt);

            if (fitStatus[ib] == 1) {
                double cf = center_abs + offset_fit[ib];
                if (cf >= tmin && cf <= tmax) hFitAfter->Fill(cf, cnt);
            } else {
                hFailedBefore->Fill(center_abs, cnt);
                double cmf = center_abs + offset_mean[ib];
                if (cmf >= tmin && cmf <= tmax) hFailedMeanAfter->Fill(cmf, cnt);
            }
        }
    }

// normalize by bin width for counts/ns
double bw = hBefore->GetBinWidth(1);
hBefore->Scale(1.0 / bw);
hMeanAfter->Scale(1.0 / bw);
hFitAfter->Scale(1.0 / bw);
hFailedBefore->Scale(1.0 / bw);
hFailedMeanAfter->Scale(1.0 / bw);


    // ensure 1D histos draw as lines
    hBefore->SetLineColor(kBlack); hBefore->SetLineWidth(1);
    hMeanAfter->SetLineColor(kBlue); hMeanAfter->SetLineWidth(2);
    hFitAfter->SetLineColor(kRed); hFitAfter->SetLineWidth(2);
    hFailedBefore->SetLineColor(kMagenta); hFailedBefore->SetLineWidth(1);
    hFailedMeanAfter->SetLineColor(kGreen+2); hFailedMeanAfter->SetLineWidth(1);

    // --- layout maps (2D) and fill them
    TH2D *hLayoutMean = new TH2D("hLayoutMean", "NPS layout: offset (mean method);col;row;offset (ns)",
                                nCols, 0.5, nCols+0.5, nRows, 0.5, nRows+0.5);
    TH2D *hLayoutFit  = new TH2D("hLayoutFit",  "NPS layout: offset (fit method);col;row;offset (ns)",
                                nCols, 0.5, nCols+0.5, nRows, 0.5, nRows+0.5);
    TH2D *hLayoutDiff = new TH2D("hLayoutDiff", "NPS layout: fit - mean (ns);col;row;diff (ns)",
                                nCols, 0.5, nCols+0.5, nRows, 0.5, nRows+0.5);

    // map block index -> col,row such that index 0 is lower-left, increase to right then up
    for (int ib=0; ib<nBlocks; ++ib) {
        int row = ib / nCols;   // 0 = bottom row
        int col = ib % nCols;   // 0 = left-most column
        int binx = col + 1;     // TH2 bins start at 1
        int biny = row + 1;
        hLayoutMean->SetBinContent(binx, biny, offset_mean[ib]);
        hLayoutFit->SetBinContent(binx, biny, offset_fit[ib]);
        hLayoutDiff->SetBinContent(binx, biny, (offset_fit[ib] - offset_mean[ib]));
    }

    // --- offset histograms
    TH1D *hOffsetMean = new TH1D("hOffsetMean", "Offset distribution (mean);offset (ns);blocks", 120, -12, 12);
    TH1D *hOffsetFit  = new TH1D("hOffsetFit",  "Offset distribution (fit);offset (ns);blocks", 120, -12, 12);
    TH1D *hOffsetDiff = new TH1D("hOffsetDiff", "Offset difference (fit - mean);diff (ns);blocks", 120, -6, 6);

    for (int ib=0; ib<nBlocks; ++ib) {
        if (ntot[ib] >= minCountsForOffset) {
            hOffsetMean->Fill(offset_mean[ib]);
            if (fitStatus[ib]==1) hOffsetFit->Fill(offset_fit[ib]);
            if (fitStatus[ib]==1) hOffsetDiff->Fill(offset_fit[ib] - offset_mean[ib]);
        }
    }
    hOffsetMean->SetLineColor(kBlue); hOffsetMean->SetLineWidth(2);
    hOffsetFit->SetLineColor(kRed);  hOffsetFit->SetLineWidth(2);
    hOffsetDiff->SetLineColor(kViolet); hOffsetDiff->SetLineWidth(2);

    // Mean vs Fit scatter
    vector<double> vMean, vFit;
    for (int ib=0; ib<nBlocks; ++ib) {
        if (fitStatus[ib]==1 && ntot[ib]>=minCountsForOffset) {
            vMean.push_back(offset_mean[ib]);
            vFit.push_back(offset_fit[ib]);
        }
    }
    TGraph *gMeanVsFit = nullptr;
    if (!vMean.empty()) {
        gMeanVsFit = new TGraph((int)vMean.size());
        for (size_t i=0;i<vMean.size();++i) gMeanVsFit->SetPoint(i, vMean[i], vFit[i]);
        gMeanVsFit->SetMarkerStyle(20);
        gMeanVsFit->SetMarkerSize(0.8);
    }

    // --- Build block-vs-time heatmap (2D): use absolute time axis [140..160]
    TH2D *hCounts2D = new TH2D("hCounts2D", "Timing histogram per block;Block;Time (ns)",
                               nBlocks, -0.5, nBlocks - 0.5, nBins, 140.0 - 0.5, 140.0 + nBins + 0.5);
    for (int ib=0; ib<nBlocks; ++ib) {
        for (int k=0; k<nBins; ++k) {
            int cnt = blktimeh[ib][k];
            if (cnt>0) hCounts2D->Fill(ib, binCenterAbs(k), (double)cnt);
        }
    }

    // --- Prepare canvases (we will save them to a multi-page PDF)
    // Page order:
    //  1) block vs time heatmap
    //  2) layout mean (with block numbers and failed-fit marks)
    //  3) layout fit (with block numbers and failed-fit marks)
    //  4) layout diff
    //  5) offset histograms overlay (mean & fit) + diff
    //  6) mean vs fit scatter
    //  7) timing overlays (before / mean / fit)
    //  8) failed-fit diagnostics (hist + list)
    //
    // We keep canvases alive until after the final SaveAs("pdf)")

    // Canvas 1 - overview
    TCanvas *c1 = new TCanvas("c1","overview",1400,900);
    c1->Divide(2,2);

    // heatmap panel (top-left)
    c1->cd(1);
    gPad->SetRightMargin(0.14);
    hCounts2D->GetZaxis()->SetTitle("counts");
    hCounts2D->Draw("COLZ");
    gPad->SetLogz(1);
    // nice ticks
    hCounts2D->GetXaxis()->SetNdivisions(510);
    hCounts2D->GetYaxis()->SetNdivisions(510);

    // offset distributions (top-right)
    c1->cd(2);
    hOffsetMean->Draw("HIST");
    hOffsetFit->Draw("HIST SAME");
    TLegend *leg = new TLegend(0.6,0.6,0.88,0.88);
    leg->AddEntry(hOffsetMean,"Mean offsets","l");
    leg->AddEntry(hOffsetFit,"Fit offsets (successful fits)","l");
    leg->SetBorderSize(0); leg->Draw();

    // mean vs fit (bottom-left)
    c1->cd(3);
    if (gMeanVsFit) {
        gMeanVsFit->Draw("AP");
        TLine *d = new TLine(-12,-12,12,12); d->SetLineStyle(2); d->SetLineColor(kGray+2); d->Draw("same");
    } else {
        TH1D *hEmpty = new TH1D("hEmpty","No fit-successful blocks",10,0,1);
        hEmpty->Draw();
    }

    // diff histogram (bottom-right)
    c1->cd(4);
    hOffsetDiff->Draw("HIST");

    // Canvas 2 - layout mean with block numbers + failed marks
    TCanvas *c2 = new TCanvas("c2","layout_mean",1400,900);
    c2->Divide(1,1);
    c2->cd();
    gPad->SetRightMargin(0.14);
    hLayoutMean->Draw("COLZ");
    // draw block numbers on top of layout
    TLatex latex;
    latex.SetTextAlign(22); // center
    latex.SetTextFont(42);
    latex.SetTextSize(0.017); // small
    for (int ib=0; ib<nBlocks; ++ib) {
        int row = ib / nCols;
        int col = ib % nCols;
        double x = double(col) + 1.0; // bin center x
        double y = double(row) + 1.0; // bin center y
        double val = hLayoutMean->GetBinContent(int(x), int(y));
        // choose contrasting color depending on bin value
        if (fabs(val) > 6.0) latex.SetTextColor(kWhite); else latex.SetTextColor(kBlack);
        latex.DrawLatex(x, y, Form("%d", ib));
    }
    // mark failed-fit blocks with a red 'X' and a small box
    for (int f : fitFailures) {
        int row = f / nCols;
        int col = f % nCols;
        double x = double(col) + 1.0;
        double y = double(row) + 1.0;
        latex.SetTextColor(kRed);
        latex.SetTextSize(0.02);
        latex.DrawLatex(x, y, "X");
    }
    // label
    TLatex tmean; tmean.SetNDC(); tmean.SetTextSize(0.03); tmean.DrawLatex(0.02,0.96, Form("Mean offsets (target %.1f ns)", idealTime));

    // Canvas 3 - layout fit
    TCanvas *c3 = new TCanvas("c3","layout_fit",1400,900);
    c3->cd();
    gPad->SetRightMargin(0.14);
    hLayoutFit->Draw("COLZ");
    latex.SetTextSize(0.017); latex.SetTextColor(kBlack);
    for (int ib=0; ib<nBlocks; ++ib) {
        int row = ib / nCols;
        int col = ib % nCols;
        double x = double(col) + 1.0;
        double y = double(row) + 1.0;
        latex.SetTextColor(kBlack);
        latex.DrawLatex(x, y, Form("%d", ib));
    }
    // mark failed fit blocks (these will be empty in fit map)
    for (int f : fitFailures) {
        int row = f / nCols;
        int col = f % nCols;
        double x = double(col) + 1.0;
        double y = double(row) + 1.0;
        latex.SetTextColor(kRed);
        latex.SetTextSize(0.02);
        latex.DrawLatex(x, y, "X");
    }
    TLatex tfit; tfit.SetNDC(); tfit.SetTextSize(0.03); tfit.DrawLatex(0.02,0.96, "Fit offsets (successful fits only)");

    // Canvas 4 - layout diff (fit - mean)
    TCanvas *c4 = new TCanvas("c4","layout_diff",1400,900);
    c4->cd();
    gPad->SetRightMargin(0.14);
    hLayoutDiff->Draw("COLZ");
    TLatex td; td.SetNDC(); td.SetTextSize(0.03); td.DrawLatex(0.02,0.96, "Fit - Mean (ns)");

    // Canvas 5 - timing overlays
    TCanvas *c5 = new TCanvas("c5","timing_overlays",1200,700);
    c5->cd();
    hBefore->Draw("HIST");
    hMeanAfter->Draw("HIST SAME");
    hFitAfter->Draw("HIST SAME");
    TLegend *legT = new TLegend(0.65,0.6,0.92,0.85);
    legT->AddEntry(hBefore,"Raw (no offset)","l");
    legT->AddEntry(hMeanAfter,"Mean-offset applied","l");
    legT->AddEntry(hFitAfter,"Fit-offset applied (fit success only)","l");
    legT->SetBorderSize(0); legT->Draw();

    // Canvas 6 - failed-fit diagnostics
    TCanvas *c6 = new TCanvas("c6","failed_fit",1000,700);
    c6->Divide(2,1);
    c6->cd(1);
    hFailedBefore->Draw("HIST");
    c6->cd(2);
    hFailedMeanAfter->Draw("HIST");

    // Canvas 7 - offset histograms zoomed
    TCanvas *c7 = new TCanvas("c7","offset_hists",1000,800);
    c7->cd();
    hOffsetMean->Draw("HIST");
    hOffsetFit->Draw("HIST SAME");
    TLegend *leg2 = new TLegend(0.6,0.7,0.88,0.88);
    leg2->AddEntry(hOffsetMean,"Mean","l");
    leg2->AddEntry(hOffsetFit,"Fit","l");
    leg2->SetBorderSize(0); leg2->Draw();

    // ---------- Save multi-page PDF properly ----------
    cout << "Writing multi-page PDF: " << outPdfName << " ...\n";
    // ensure drawn
    c1->Update();
    c2->Update();
    c3->Update();
    c4->Update();
    c5->Update();
    c6->Update();
    c7->Update();

    // Start PDF
    c1->SaveAs((outPdfName + "(").Data());
    // other pages (same filename adds pages)
    c2->SaveAs(outPdfName.Data());
    c3->SaveAs(outPdfName.Data());
    c4->SaveAs(outPdfName.Data());
    c5->SaveAs(outPdfName.Data());
    c6->SaveAs(outPdfName.Data());
    c7->SaveAs(outPdfName.Data());
    // close multi-page PDF
    TCanvas cend("cend","end",200,200);
    cend.SaveAs((outPdfName + ")").Data());
    cout << "PDF saved to " << outPdfName << " (multi-page)\n";

    // --- Write ROOT file
    cout << "Writing ROOT file: " << outRootName << "\n";
    TFile fout(outRootName, "RECREATE");
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
    hOffsetDiff->Write();
    fout.Close();

    // --- Write text outputs: fittime.out, tdc_offset.param (mean), tdc_offset_fit.param (fit), fittime_summary.txt
    ofstream foutTxt(outTxtName.Data());
    ofstream foutMean(paramMeanName.Data());
    ofstream foutFit(paramFitName.Data());
    ofstream foutSummary(summaryName.Data());

    if (!foutTxt.is_open() || !foutMean.is_open() || !foutFit.is_open() || !foutSummary.is_open()) {
        cerr << "ERROR: could not open one or more text outputs for writing.\n";
    } else {
        foutTxt << "# block ntot mean_time offset_mean offset_fit fitStatus fit_mean fit_sigma fit_chi2\n";
        for (int ib=0; ib<nBlocks; ++ib) {
            foutTxt << setw(4) << ib
                    << setw(8) << ntot[ib]
                    << setw(10) << fixed << setprecision(2) << mean_time[ib]
                    << setw(10) << fixed << setprecision(3) << offset_mean[ib]
                    << setw(10) << fixed << setprecision(3) << offset_fit[ib]
                    << setw(5) << fitStatus[ib]
                    << setw(10) << fixed << setprecision(3) << fit_mean[ib]
                    << setw(10) << fixed << setprecision(3) << fit_sigma[ib]
                    << setw(10) << fixed << setprecision(3) << fit_chi2[ib]
                    << "\n";
        }
        foutTxt.close();

        // write param files: 108 rows x 10 cols (~f5.2, comma separated)
        const int colsPerRow = 10;
        int nRows = (nBlocks + colsPerRow - 1) / colsPerRow;
        for (int r=0;r<nRows;++r) {
            int k1 = r*colsPerRow;
            int k2 = min(k1+colsPerRow-1, nBlocks-1);
            for (int k=k1; k<=k2; ++k) {
                foutMean << fixed << setprecision(2) << setw(5) << offset_mean[k];
                if (k!=k2) foutMean << ",";
            }
            if (k2-k1+1 < colsPerRow) {
                for (int p=0; p<colsPerRow-(k2-k1+1); ++p) foutMean << ", 0.00";
            }
            foutMean << "\n";
        }
        foutMean.close();

        // fit-based param file - write offsets where fit succeeded, otherwise 0.00
        for (int r=0;r<nRows;++r) {
            int k1 = r*colsPerRow;
            int k2 = min(k1+colsPerRow-1, nBlocks-1);
            for (int k=k1;k<=k2;++k) {
                double v = (fitStatus[k]==1) ? offset_fit[k] : 0.0;
                foutFit << fixed << setprecision(2) << setw(5) << v;
                if (k!=k2) foutFit << ",";
            }
            if (k2-k1+1 < colsPerRow) {
                for (int p=0;p<colsPerRow-(k2-k1+1);++p) foutFit << ", 0.00";
            }
            foutFit << "\n";
        }
        foutFit.close();

        // summary
        foutSummary << "# fittime summary\n";
        foutSummary << "# input file: " << blockFile << "\n";
        foutSummary << "# blocks: " << nBlocks << "  bins: " << nBins << "\n";
        foutSummary << "# minCountsForOffset: " << minCountsForOffset << "  idealTime: " << idealTime << "\n";
        foutSummary << "# total blocks read: " << readBlocks << "\n";
        foutSummary << "# mean offsets: +ve=" << npos_mean << "  -ve=" << nneg_mean << "\n";
        foutSummary << "# fit attempts: " << fitAttempts << "  success: " << fitSuccess << "  fails: " << fitFails << "\n";
        foutSummary << "# fit failures list (" << fitFailures.size() << "):\n";
        for (size_t i=0;i<fitFailures.size();++i) {
            foutSummary << fitFailures[i] << ((i%20==19) ? "\n" : " ");
        }
        foutSummary << "\n";
        foutSummary.close();
    }

    // final console summary
    cout << "\n=== Summary ===\n";
    cout << "Blocks read:           " << readBlocks << "\n";
    cout << "Blocks with mean offset (ntot >= " << minCountsForOffset << "): " << npos_mean + nneg_mean << "\n";
    cout << "  mean offsets: +ve=" << npos_mean << "  -ve=" << nneg_mean << "\n";
    cout << "Fit attempts:          " << fitAttempts << "  success: " << fitSuccess << "  fails: " << fitFails << "\n";
    cout << "Fit failures recorded: " << fitFailures.size() << " (listed in " << summaryName << ")\n";
    cout << "Outputs written to: " << outdir << "\n";
    cout << "  " << outRootName << "\n";
    cout << "  " << outPdfName << "\n";
    cout << "  " << outTxtName << "\n";
    cout << "  " << paramMeanName << "  (mean-based)\n";
    cout << "  " << paramFitName << "   (fit-based: only successful fits written, rest 0.00)\n";
    cout << "  " << summaryName << "\n";
    cout << "=== Done ===\n";

    // delete histos/graphs (ROOT file has stored them into fout root already)
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
    delete hOffsetDiff;

        // Cleanup: do not aggressively delete histos/canvases if you want to inspect them in interactive session.
    // But free the canvases we created to avoid leaks when used in batch runs:
    delete c1;
    delete c2;
    delete c3;
    delete c4;
    delete c5;
    delete c6;
    delete c7;

}
