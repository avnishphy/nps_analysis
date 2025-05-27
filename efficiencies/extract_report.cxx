#include <TSystemDirectory.h>
#include <TSystemFile.h>
#include <TList.h>
#include <TRegexp.h>
#include <TString.h>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>
#include <sstream>
#include <iomanip>

double extract_value(const std::string &line) {
    std::istringstream iss(line);
    std::string word;
    while (iss >> word) {
        try {
            return std::stod(word);
        } catch (...) {
            continue;
        }
    }
    return -999;
}

void extract_efficiencies_summary() {
    TString dir_path = "/w/hallc-scshelf2102/nps/nps-ana/REPORT_OUTPUT_pass1/COIN/SKIM";
    TSystemDirectory dir("report_dir", dir_path);
    TList *files = dir.GetListOfFiles();
    TRegexp pattern("skim_NPS_HMS_report_[0-9]+_-1\\.report");

    std::map<int, std::map<TString, std::vector<double>>> run_data;

    std::map<TString, TString> targets = {
        {"BCM4A Charge", "charge"},
        {"PS factor", "ps_factor"},
        {"HMS_TRACK_EFF", "tracking_eff"},
        {"HGCER Efficiency", "hgc_eff"},
        {"HODO Efficiency", "hod_eff"},
        {"HMS EDTM Rate", "rate_edtm"},
        {"HMS PHYSICS Rate", "rate_phys"},
        {"HMS3/4 Rate", "rate_hms34"}
    };

    if (!files) {
        std::cerr << "No files found in the directory." << std::endl;
        return;
    }

    TIter next(files);
    TSystemFile *file;
    while ((file = (TSystemFile *)next())) {
        TString fname = file->GetName();
        if (!file->IsDirectory() && fname.Contains(pattern)) {
            int run = fname.ReplaceAll("skim_NPS_HMS_report_", "")
                          .ReplaceAll("_-1.report", "")
                          .Atoi();
            TString full_path = dir_path + "/" + fname;

            std::ifstream fin(full_path.Data());
            if (!fin.is_open()) {
                std::cerr << "Cannot open: " << full_path << std::endl;
                continue;
            }

            std::vector<std::string> lines;
            std::string line;
            while (getline(fin, line)) {
                lines.push_back(line);
            }

            // Trigger Times (lines 132–137)
            for (int i = 132; i <= 137 && i < lines.size(); ++i) {
                double val = extract_value(lines[i]);
                if (val >= 0)
                    run_data[run]["trig_time"].push_back(val);
            }

            // HMS Live Time (line 86)
            if (lines.size() > 86) {
                double val = extract_value(lines[86]);
                if (val >= 0)
                    run_data[run]["lt_hms"].push_back(val);
            }

            // EDTM Live Time (line 125)
            if (lines.size() > 125) {
                double val = extract_value(lines[125]);
                if (val >= 0)
                    run_data[run]["lt_edtm"].push_back(val);
            }

            // Keyword scan
            for (auto &l : lines) {
                for (auto &[key, label] : targets) {
                    if (l.find(key.Data()) != std::string::npos) {
                        double val = extract_value(l);
                        if (val >= 0)
                            run_data[run][label].push_back(val);
                        break;
                    }
                }
            }
        }
    }

    // Write to CSV
    std::ofstream fout("/w/hallc-scshelf2102/nps/singhav/nps_analysis/efficiencies/report_summary.csv");
    fout << "Run,Avg_PS_Factor,Avg_Trig_Time,Charge,Tracking_Eff,HGC_Eff,Hod_Eff,LT_EDTM,LT_HMS,Rate_EDTM,Rate_PHYS,Rate_HMS34\n";

    for (const auto &[run, entries] : run_data) {
        fout << run;

        auto avg = [&entries](const TString &key) {
            const auto &vec = entries.at(key);
            double sum = 0;
            for (auto v : vec) sum += v;
            return (vec.empty() ? -1.0 : sum / vec.size());
        };
        

        fout << "," << avg("ps_factor");
        fout << "," << avg("trig_time");
        fout << "," << avg("charge");
        fout << "," << avg("tracking_eff");
        fout << "," << avg("hgc_eff");
        fout << "," << avg("hod_eff");
        fout << "," << avg("lt_edtm");
        fout << "," << avg("lt_hms");
        fout << "," << avg("rate_edtm");
        fout << "," << avg("rate_phys");
        fout << "," << avg("rate_hms34") << "\n";
    }

    fout.close();
    std::cout << "✅ Summary written to report_summary.csv\n";
}
