#include <TCanvas.h>
#include <TGraph.h>
#include <TApplication.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <sys/stat.h>

bool file_exists(const std::string& filename) {
    struct stat buffer;
    return (stat(filename.c_str(), &buffer) == 0);
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: ./YieldScaler <start_run> <end_run>\n";
        return 1;
    }

    int run_start = std::atoi(argv[1]);
    int run_end = std::atoi(argv[2]);

    std::string base_dir = "/w/hallc-scshelf2102/nps/nps-ana/REPORT_OUTPUT_pass1/COIN/SKIM/";
    std::ofstream debug("yield_scaler_debug.txt");

    std::vector<double> currents, yields;
    std::vector<int> runs;

    for (int run = run_start; run <= run_end; ++run) {
        std::string filename = base_dir + "skim_NPS_HMS_report_" + std::to_string(run) + "_-1.report";
        if (!file_exists(filename)) {
            debug << "WARNING: File not found: " << filename << "\n";
            continue;
        }

        std::ifstream infile(filename);
        std::string line;
        int line_number = 0;
        double bcm_charge = -1, bcm_current = -1, edtm_trig = -1;
        int ps_index = -1;
        double htrig_val = -1;
        std::map<int, int> ps_factors;
        std::map<int, double> htrigs;

        while (std::getline(infile, line)) {
            ++line_number;

            if (line_number == 40 && line.find("BCM4A Beam Cut Current") != std::string::npos)
                sscanf(line.c_str(), "BCM4A Beam Cut Current: %lf uA", &bcm_current);

            if (line_number == 47 && line.find("BCM4A Beam Cut Charge") != std::string::npos)
                sscanf(line.c_str(), "BCM4A Beam Cut Charge: %lf uC", &bcm_charge);

            if (line_number >= 52 && line_number <= 61) {
                for (int i = 1; i <= 6; ++i) {
                    std::string key = "Ps" + std::to_string(i) + "_factor";
                    if (line.find(key) != std::string::npos) {
                        int val;
                        sscanf(line.c_str(), (key + " = %d").c_str(), &val);
                        ps_factors[i] = val;
                    }
                }
            }

            if (line_number >= 109 && line_number <= 114) {
                for (int i = 1; i <= 6; ++i) {
                    std::string key = "hTRIG" + std::to_string(i);
                    if (line.find(key) != std::string::npos) {
                        double val;
                        sscanf(line.c_str(), (key + " : %lf").c_str(), &val);
                        htrigs[i] = val;
                    }
                }
            }

            if (line_number == 120 && line.find("EDTM Triggers") != std::string::npos)
                sscanf(line.c_str(), "EDTM Triggers           : %lf", &edtm_trig);
        }

        infile.close();

        for (auto& kv : ps_factors) {
            if (kv.second > 0) {
                ps_index = kv.first;
                break;
            }
        }

        if (ps_index == -1 || htrigs.find(ps_index) == htrigs.end()) {
            debug << "Run " << run << ": Missing PS or hTRIG info.\n";
            continue;
        }

        htrig_val = htrigs[ps_index];
        double number_scaler = htrig_val - edtm_trig;
        double yield = (bcm_charge > 0) ? number_scaler / bcm_charge : -1;

        debug << "Run: " << run << "\n";
        debug << "Beam Current: " << bcm_current << " uA\n";
        debug << "Charge: " << bcm_charge << " uC\n";
        debug << "hTRIG" << ps_index << ": " << htrig_val << "\n";
        debug << "EDTM Triggers: " << edtm_trig << "\n";
        debug << "Number Scaler: " << number_scaler << "\n";
        debug << "Yield Scaler: " << yield << "\n\n";

        runs.push_back(run);
        currents.push_back(bcm_current);
        yields.push_back(yield);
    }

    debug.close();

    TApplication app("app", &argc, argv);
    TCanvas* c = new TCanvas("c", "Yield vs. Beam Current", 800, 600);
    TGraph* graph = new TGraph(currents.size());

    for (size_t i = 0; i < currents.size(); ++i)
        graph->SetPoint(i, currents[i], yields[i]);

    graph->SetTitle("Yield Scaler vs. Beam Current;Beam Current (uA);Yield Scaler (counts/uC)");
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(kBlue);
    graph->Draw("AP");

    c->SaveAs("yield_scaler_vs_current.png");

    return 0;
}
