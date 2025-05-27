#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <string>
#include <numeric>
#include <regex>
#include <filesystem>

using namespace std;
namespace fs = std::filesystem;

TF1* FitPeak(TH1D* hist) {
    TF1* gaus = new TF1("gaus", "gaus", 0.7, 1.5);
    hist->Fit(gaus, "Q", "", 0.7, 1.5);
    return gaus;
}

map<int, double> ReadCurrentCSV(const string& csv_path) {
    map<int, double> run_current;
    ifstream file(csv_path);
    string line;
    getline(file, line); // skip header

    while (getline(file, line)) {
        stringstream ss(line);
        string run_str, current_str;
        getline(ss, run_str, ',');
        getline(ss, current_str, ',');

        int run = stoi(run_str);
        double current = stod(current_str);
        run_current[run] = current;
    }
    return run_current;
}

void cal_peak_vs_current() {
    string root_dir = "/lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/luminosity_all";
    string current_csv = "/w/hallc-scshelf2102/nps/singhav/nps_analysis/luminosity_analysis/livetime_results.csv";
    string output_csv = "/w/hallc-scshelf2102/nps/singhav/nps_analysis/luminosity_analysis/cal_peak_vs_current_cpp.csv";

    map<int, vector<string>> run_files;
    for (const auto& entry : fs::directory_iterator(root_dir)) {
        string fname = entry.path().filename();
        smatch match;
        regex rgx("nps_hms_coin_(\\d+)_\\d+_1_-1\\.root");
        if (regex_search(fname, match, rgx)) {
            int run = stoi(match[1]);
            run_files[run].push_back(entry.path());
        }
    }

    map<int, double> run_current = ReadCurrentCSV(current_csv);

    vector<double> v_run, v_current, v_etotnorm, v_etottracknorm;

    ofstream out(output_csv);
    out << "run,current,etotnorm_peak,etottracknorm_peak\n";

    for (auto& [run, files] : run_files) {
        cout << "Processing run " << run << " with " << files.size() << " segment(s)..." << endl;
        vector<double> norm_peaks, track_peaks;

        for (const auto& fpath : files) {
            TFile* f = TFile::Open(fpath.c_str());
            if (!f || f->IsZombie()) {
                cerr << "Failed to open " << fpath << endl;
                continue;
            }

            TTree* T = (TTree*)f->Get("T");
            double_t etotnorm = 0, etottracknorm = 0;
            T->SetBranchStatus("*", 0);
            T->SetBranchStatus("H.cal.etotnorm", 1);
            T->SetBranchStatus("H.cal.etottracknorm", 1);
            T->SetBranchAddress("H.cal.etotnorm", &etotnorm);
            T->SetBranchAddress("H.cal.etottracknorm", &etottracknorm);

            TH1D* h1 = new TH1D("h1", "H.cal.etotnorm", 100, 0.7, 1.5);
            TH1D* h2 = new TH1D("h2", "H.cal.etottracknorm", 100, 0.7, 1.5);

            for (Long64_t i = 0; i < T->GetEntries(); ++i) {
                T->GetEntry(i);
                if (etotnorm > 0.7) h1->Fill(etotnorm);
                if (etottracknorm > 0.7) h2->Fill(etottracknorm);
            }

            TF1* fit1 = FitPeak(h1);
            TF1* fit2 = FitPeak(h2);

            if (fit1 && fit2) {
                norm_peaks.push_back(fit1->GetParameter(1));
                track_peaks.push_back(fit2->GetParameter(1));
            }

            delete h1; delete h2;
            f->Close();
        }

        if (!norm_peaks.empty() && !track_peaks.empty()) {
            double avg_norm = accumulate(norm_peaks.begin(), norm_peaks.end(), 0.0) / norm_peaks.size();
            double avg_track = accumulate(track_peaks.begin(), track_peaks.end(), 0.0) / track_peaks.size();

            double current = run_current.count(run) ? run_current[run] : -1;
            out << run << "," << current << "," << avg_norm << "," << avg_track << "\n";

            cout << "Run " << run << ": current=" << current << ", avg_norm=" << avg_norm << ", avg_track=" << avg_track << endl;

            if (current > 0) {
                v_run.push_back(run);
                v_current.push_back(current);
                v_etotnorm.push_back(avg_norm);
                v_etottracknorm.push_back(avg_track);
            }
        }
    }

    out.close();
    cout << "Results saved to " << output_csv << endl;

    // ======= Plotting ======= //
    gStyle->SetOptStat(0);

    TCanvas* c1 = new TCanvas("c1", "etotnorm peak vs current", 800, 600);
    TGraph* g1 = new TGraph(v_current.size(), &v_current[0], &v_etotnorm[0]);
    g1->SetTitle("H.cal.etotnorm Peak vs Beam Current;Beam Current (#muA);etotnorm Peak Position");
    g1->SetMarkerStyle(20);
    g1->SetMarkerColor(kBlue);
    g1->Draw("AP");
    c1->SaveAs("etotnorm_peak_vs_current.png");

    TCanvas* c2 = new TCanvas("c2", "etottracknorm peak vs current", 800, 600);
    TGraph* g2 = new TGraph(v_current.size(), &v_current[0], &v_etottracknorm[0]);
    g2->SetTitle("H.cal.etottracknorm Peak vs Beam Current;Beam Current (#muA);etottracknorm Peak Position");
    g2->SetMarkerStyle(21);
    g2->SetMarkerColor(kRed);
    g2->Draw("AP");
    c2->SaveAs("etottracknorm_peak_vs_current.png");

    cout << "Plots saved as PNGs." << endl;
}
