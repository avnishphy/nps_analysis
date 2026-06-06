#include <filesystem>
#include <set>
#include <iomanip>

//
// Compilation:
//   g++ -O2 -std=c++17 -o check_edtm_tdc check_edtm_tdc.cxx `root-config --cflags --libs`
//
// # Usage:
// #   ./check_edtm_tdc "<run_list>" <root_dir> <output_pdf>
// #
// # <run_list> can be:
// #   - Individual runs: 4387 4390 4400
// #   - Ranges: 4387-4390
// #   - Comma or space separated: 4387,4390 4400-4402
// #   - Any combination: "4387-4390,4400 4500-4502"
// #
// # Example:
// #   ./check_edtm_tdc "4387-4390,4400 4500-4502" /path/to/rootfiles out.pdf


#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>      // For gROOT
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TROOT.h>
#include <sstream>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <map>
#include <filesystem>
#include <iomanip>

namespace fs = std::filesystem;

// Scan directory ONCE and build run→filelist map
std::map<int, std::vector<std::string>> build_run_file_map(const std::string& dir, const std::vector<int>& runs) {
    std::map<int, std::vector<std::string>> run_files;
    std::set<int> run_set(runs.begin(), runs.end());
    for (const auto& entry : fs::directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;
        std::string fname = entry.path().filename().string();
        if (fname.find(".root") == std::string::npos) continue;
        for (int run : run_set) {
            if (fname.find(std::to_string(run)) != std::string::npos) {
                run_files[run].push_back(entry.path().string());
            }
        }
    }
    // Sort file lists for each run
    for (auto& kv : run_files) std::sort(kv.second.begin(), kv.second.end());
    return run_files;
}
    // (Removed stray/erroneous line)

double find_main_peak(const TH1D* h) {
    int maxbin = h->GetMaximumBin();
    return h->GetBinCenter(maxbin);
}
// (Removed unused/erroneous get_files_for_run declaration)

// Simple progress bar (single-threaded, safe for bash parallelization)
class ProgressBar {
    int total;
    int width;
public:
    ProgressBar(int total_, int width_ = 50) : total(total_), width(width_) {}
    void print(int current) {
        int pos = (total > 0) ? (current * width) / total : 0;
        std::cout << "\r[";
        for (int i = 0; i < width; ++i) std::cout << (i < pos ? "=" : " ");
        std::cout << "] " << current << "/" << total << std::flush;
        if (current == total) std::cout << std::endl;
    }
};

// Process a single run: build TChain, fill hist, return peak
struct RunResult {
    int run;
    double peak;
    TH1D* hist;
};

RunResult process_run(int run, const std::vector<std::string>& files) {
    TChain chain("T");
    for (const auto& f : files) chain.Add(f.c_str());
    auto* hist = new TH1D(Form("h%d", run), Form("Run %d;T.hms.hEDTM_tdcTimeRaw;Counts", run), 2000, 0, 20000);
    hist->SetDirectory(nullptr);
    TTreeReader reader(&chain);
    TTreeReaderValue<double> tdc(reader, "T.hms.hEDTM_tdcTimeRaw");
    while (reader.Next()) {
        if (*tdc > 0) hist->Fill(*tdc);
    }
    double peak = find_main_peak(hist);
    return {run, peak, hist};
}

void plot_results(const std::vector<RunResult>& results, const std::string& out_pdf) {
    double global_max_peak = 0;
    for (const auto& r : results) if (r.peak > global_max_peak) global_max_peak = r.peak;
    double xmax = global_max_peak + 10000;
    gStyle->SetOptStat(0);
    gROOT->SetBatch(kTRUE);
    std::unique_ptr<TCanvas> c(new TCanvas("c", "EDTM TDC Spectra", 1200, 800));
    c->SetLogy();
    std::unique_ptr<TLegend> leg(new TLegend(0.7, 0.7, 0.95, 0.95));
    int color = 1;
    for (size_t i = 0; i < results.size(); ++i) {
        auto* h = results[i].hist;
        h->GetXaxis()->SetRangeUser(0, xmax);
        h->SetLineColor(color++);
        h->SetLineWidth(2);
        if (i == 0) h->Draw("hist");
        else h->Draw("hist same");
        leg->AddEntry(h, Form("Run %d (peak=%.0f)", results[i].run, results[i].peak), "l");
    }
    leg->Draw();
    c->Modified();
    c->Update();
    c->Print(out_pdf.c_str(), "pdf");
    for (const auto& r : results) delete r.hist;
}
// (Removed stray/erroneous for loop)
// Parse run list string (supports ranges, commas, spaces)
std::vector<int> parse_run_list(const std::string& run_str) {
    std::vector<int> runs;
    std::string s = run_str;
    std::replace(s.begin(), s.end(), ',', ' ');
    std::istringstream iss(s);
    std::string token;
    while (iss >> token) {
        size_t dash = token.find('-');
        if (dash != std::string::npos) {
            int start = std::stoi(token.substr(0, dash));
            int end = std::stoi(token.substr(dash + 1));
            for (int r = start; r <= end; ++r) runs.push_back(r);
        } else {
            runs.push_back(std::stoi(token));
        }
    }
    std::sort(runs.begin(), runs.end());
    runs.erase(std::unique(runs.begin(), runs.end()), runs.end());
    return runs;
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <run_list> <root_dir> <output_pdf>\n";
        std::cerr << "Example: " << argv[0] << " \"4387-4390,4400 4500-4502\" /path/to/rootfiles out.pdf\n";
        return 1;
    }
    std::string run_str = argv[1];
    std::string root_dir = argv[2];
    std::string out_pdf = argv[3];
    std::vector<int> runs = parse_run_list(run_str);
    auto run_file_map = build_run_file_map(root_dir, runs);
    std::vector<RunResult> results;
    results.reserve(runs.size());
    ProgressBar pbar(runs.size());
    int count = 0;
    for (int run : runs) {
        auto it = run_file_map.find(run);
        if (it == run_file_map.end() || it->second.empty()) {
            std::cerr << "No files for run " << run << std::endl;
            ++count;
            pbar.print(count);
            continue;
        }
        results.push_back(process_run(run, it->second));
        ++count;
        pbar.print(count);
    }
    if (!results.empty()) plot_results(results, out_pdf);
    return 0;
}