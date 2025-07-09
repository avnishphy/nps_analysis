// Compile with:
// g++ stabilityCheck.cpp `root-config --cflags --libs` -o stabilityCheck

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TF1.h>
#include <TLegend.h>
#include <TStyle.h>

#include <iostream>
#include <vector>
#include <string>
#include <map>
#include <fstream>
#include <filesystem>
#include <regex>
#include <cmath>

namespace fs = std::filesystem;

bool pass_cuts(double edt, double dp, double et, double npe, double th, double ph, double beta) {
    return (edt < 0.1 && std::abs(dp) <= 8.5 && et > 0.9 && npe > 1.0 &&
            std::abs(th) <= 0.09 && std::abs(ph) <= 0.09 && std::abs(beta - 1) < 0.5);
}

void plot_graph(const std::vector<int>& runs, const std::vector<double>& values, const std::vector<double>& errors,
                const std::string& title, const std::string& yLabel, const std::string& filename, int markerStyle, int color) {
    std::vector<double> run_d(runs.begin(), runs.end());
    TCanvas* c = new TCanvas(filename.c_str(), title.c_str(), 800, 600);
    TGraphErrors* graph = new TGraphErrors(runs.size(), run_d.data(), values.data(), nullptr, errors.data());
    graph->SetTitle((title + ";Run Number;" + yLabel).c_str());
    graph->SetMarkerStyle(markerStyle);
    graph->SetMarkerColor(color);
    graph->SetLineColor(color);
    graph->Draw("APL");

    c->SetGrid();
    c->SaveAs((filename + ".pdf").c_str());
}

int main() {
    std::string dir_path = "/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/updated";
    std::regex file_pattern(R"(nps_hms_coin_(\d+)_(\d+)_1_-1\.root)");
    std::map<int, TChain*> run_chains;

    for (const auto& entry : fs::directory_iterator(dir_path)) {
        std::string filename = entry.path().filename();
        std::smatch match;
        if (std::regex_match(filename, match, file_pattern)) {
            int run = std::stoi(match[1]);
            std::string full_path = entry.path().string();
            if (run_chains.find(run) == run_chains.end()) {
                run_chains[run] = new TChain("T");
            }
            run_chains[run]->Add(full_path.c_str());
        }
    }

    // Results to write
    std::ofstream csv_out("stability_peaks.csv");
    csv_out << "Run,Entries,Etotnorm_Peak,Etotnorm_Err,Etotnorm_Sigma,Etotnorm_SigmaErr,Etotnorm_Chi2NDF,"
            << "Etottracknorm_Peak,Etottracknorm_Err,Etottracknorm_Sigma,Etottracknorm_SigmaErr,Etottracknorm_Chi2NDF,"
            << "Etracknorm_Peak,Etracknorm_Err,Etracknorm_Sigma,Etracknorm_SigmaErr,Etracknorm_Chi2NDF\n";

    std::vector<int> runs;
    std::vector<double> etotnorm_peaks, etotnorm_errors;
    std::vector<double> etottracknorm_peaks, etottracknorm_errors;
    std::vector<double> etracknorm_peaks, etracknorm_errors;

    for (const auto& [run, chain] : run_chains) {
        chain->SetBranchStatus("*", 0);

        double edt, dp, beta, npe, th, ph;
        double etot, etotnorm, etrack, etracknorm, etottracknorm;

        chain->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
        chain->SetBranchStatus("H.gtr.dp", 1);
        chain->SetBranchStatus("H.gtr.beta", 1);
        chain->SetBranchStatus("H.cer.npeSum", 1);
        chain->SetBranchStatus("H.gtr.th", 1);
        chain->SetBranchStatus("H.gtr.ph", 1);
        chain->SetBranchStatus("H.cal.etot", 1);
        chain->SetBranchStatus("H.cal.etotnorm", 1);
        chain->SetBranchStatus("H.cal.etottracknorm", 1);
        chain->SetBranchStatus("H.cal.etrack", 1);
        chain->SetBranchStatus("H.cal.etracknorm", 1);

        chain->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edt);
        chain->SetBranchAddress("H.gtr.dp", &dp);
        chain->SetBranchAddress("H.gtr.beta", &beta);
        chain->SetBranchAddress("H.cer.npeSum", &npe);
        chain->SetBranchAddress("H.gtr.th", &th);
        chain->SetBranchAddress("H.gtr.ph", &ph);
        chain->SetBranchAddress("H.cal.etot", &etot);
        chain->SetBranchAddress("H.cal.etotnorm", &etotnorm);
        chain->SetBranchAddress("H.cal.etottracknorm", &etottracknorm);
        chain->SetBranchAddress("H.cal.etrack", &etrack);
        chain->SetBranchAddress("H.cal.etracknorm", &etracknorm);

        TH1D h_etotnorm(Form("h_etotnorm_%d", run), "", 100, 0.5, 1.5);
        TH1D h_etottracknorm(Form("h_etottracknorm_%d", run), "", 100, 0.5, 1.5);
        TH1D h_etracknorm(Form("h_etracknorm_%d", run), "", 100, 0.5, 1.5);

        Long64_t nEntries = chain->GetEntries();
        int passed = 0;

        for (Long64_t i = 0; i < nEntries; ++i) {
            chain->GetEntry(i);
            if (pass_cuts(edt, dp, etotnorm, npe, th, ph, beta)) {
                h_etotnorm.Fill(etotnorm);
                h_etottracknorm.Fill(etottracknorm);
                h_etracknorm.Fill(etracknorm);
                passed++;
            }
        }

        if (passed < 50) {
            std::cerr << "Run " << run << " has insufficient entries after cuts.\n";
            continue;
        }

        auto fit_gaussian = [](TH1D& h) {
            TF1 fit("gaus", "gaus", 0.7, 1.3);
            h.Fit(&fit, "RQ0");
            double mu = fit.GetParameter(1);
            double mu_err = fit.GetParError(1);
            double sigma = fit.GetParameter(2);
            double sigma_err = fit.GetParError(2);
            double chi2ndf = fit.GetChisquare() / fit.GetNDF();
            return std::tuple{mu, mu_err, sigma, sigma_err, chi2ndf};
        };

        auto [etn_peak, etn_err, etn_sig, etn_sigerr, etn_chi2] = fit_gaussian(h_etotnorm);
        auto [ett_peak, ett_err, ett_sig, ett_sigerr, ett_chi2] = fit_gaussian(h_etottracknorm);
        auto [etk_peak, etk_err, etk_sig, etk_sigerr, etk_chi2] = fit_gaussian(h_etracknorm);

        // Save values
        runs.push_back(run);
        etotnorm_peaks.push_back(etn_peak);
        etotnorm_errors.push_back(etn_err);
        etottracknorm_peaks.push_back(ett_peak);
        etottracknorm_errors.push_back(ett_err);
        etracknorm_peaks.push_back(etk_peak);
        etracknorm_errors.push_back(etk_err);

        csv_out << run << "," << passed << ","
                << etn_peak << "," << etn_err << "," << etn_sig << "," << etn_sigerr << "," << etn_chi2 << ","
                << ett_peak << "," << ett_err << "," << ett_sig << "," << ett_sigerr << "," << ett_chi2 << ","
                << etk_peak << "," << etk_err << "," << etk_sig << "," << etk_sigerr << "," << etk_chi2 << "\n";

        std::cout << "Run " << run << ": "
                  << "Etotnorm peak = " << etn_peak << " ± " << etn_err << ", "
                  << "Etottracknorm peak = " << ett_peak << ", "
                  << "Etracknorm peak = " << etk_peak << std::endl;
    }

    csv_out.close();

    // === Plotting ===
    plot_graph(runs, etotnorm_peaks, etotnorm_errors,
               "Etotnorm Peak vs Run", "Peak of H.cal.etotnorm", "etotnorm_peak_vs_run", 20, kBlack);
    plot_graph(runs, etottracknorm_peaks, etottracknorm_errors,
               "Etottracknorm Peak vs Run", "Peak of H.cal.etottracknorm", "etottracknorm_peak_vs_run", 21, kBlue + 2);
    plot_graph(runs, etracknorm_peaks, etracknorm_errors,
               "Etracknorm Peak vs Run", "Peak of H.cal.etracknorm", "etracknorm_peak_vs_run", 22, kRed + 1);

    return 0;
}
