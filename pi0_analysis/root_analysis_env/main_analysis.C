// main_analysis.C
// Master script that integrates all modules and executes the pipeline.
// It loads config files, builds run lists, and calls modules in order.
// Usage:
// root> .L main_analysis.C+
// root> main_analysis();

#include "utils.h"

// Forward declarations for macros (they must be loaded/compiled first)
bool skim_data(const std::vector<int>&, const std::string&, const std::string&, const std::string&, const std::string&);
bool pid_analysis(const std::vector<int>&, const std::string&, const std::string&, const std::string&, const std::string&, const std::string&);
bool efficiency_calc(const std::vector<int>&, const std::string&, const std::string&, const std::string&, const std::string&, const std::string&);
bool invariant_mass(const std::vector<int>&, const std::string&, const std::string&, const std::string&, const std::string&, const std::string&);
bool background_subtract(const std::vector<int>&, const std::string&, const std::string&, const std::string&);

#include "TStopwatch.h"
#include <vector>
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

void main_analysis() {
    // Flags to enable/disable stages
    bool RUN_SKIM = true;
    bool RUN_PID = true;
    bool RUN_EFF = true;
    bool RUN_MASS = true;
    bool RUN_BKG = true;

    // Paths and config
    std::string cfg_dir = "./config";
    std::string branch_file = cfg_dir + "/branches_to_read.txt";
    std::string cuts_file = cfg_dir + "/cuts.conf";
    std::string runlist_file = cfg_dir + "/runlist.txt";

    // Read runlist
    vector<int> runs;
    ifstream rin(runlist_file.c_str());
    if(!rin.is_open()) {
        logmsg(ERROR, "Cannot open runlist: " + runlist_file);
        return;
    }
    string line;
    while(getline(rin,line)) {
        line = trim(line);
        if(line.empty() || line[0]=='#') continue;
        runs.push_back(atoi(line.c_str()));
    }
    rin.close();

    if(runs.empty()) {
        logmsg(ERROR, "No runs found in runlist; exiting.");
        return;
    }

    // Read simple config keys
    auto kv = readKeyValueConfig(cuts_file);
    std::string raw_pattern = (kv.find("DEFAULT_INPUT_DIR")!=kv.end()) ? (kv["DEFAULT_INPUT_DIR"] + "/nps_%d_*.root") : std::string("/data/nps_%d_*.root");
    std::string skim_out_dir = (kv.find("DEFAULT_OUTPUT_DIR")!=kv.end()) ? kv["DEFAULT_OUTPUT_DIR"] : std::string("output/skimmed");
    std::string plots_dir = "output/plots";
    std::string logs_dir = "output/logs";

    // Print summary
    logmsg(INFO, "Starting main_analysis with runs:");
    for(auto r : runs) std::cout << r << " ";
    std::cout << std::endl;

    TStopwatch totalTimer; totalTimer.Start();

    if(RUN_SKIM) {
        TStopwatch t; t.Start();
        logmsg(INFO, "==== RUNNING SKIM ====");
        bool ok = skim_data(runs, raw_pattern, skim_out_dir, branch_file, "T");
        t.Stop();
        logmsg(INFO, TString::Format("SKIM stage done in %.2f s", t.RealTime()).Data());
        if(!ok) logmsg(WARN, "Skim stage returned false");
    }

    if(RUN_PID) {
        TStopwatch t; t.Start();
        logmsg(INFO, "==== RUNNING PID ====");
        bool ok = pid_analysis(runs, skim_out_dir, plots_dir, branch_file, cuts_file, "skim");
        t.Stop();
        logmsg(INFO, TString::Format("PID stage done in %.2f s", t.RealTime()).Data());
        if(!ok) logmsg(WARN, "PID stage returned false");
    }

    if(RUN_EFF) {
        TStopwatch t; t.Start();
        logmsg(INFO, "==== RUNNING EFFICIENCY ====");
        bool ok = efficiency_calc(runs, skim_out_dir, plots_dir, branch_file, cuts_file, "skim");
        t.Stop();
        logmsg(INFO, TString::Format("EFF stage done in %.2f s", t.RealTime()).Data());
        if(!ok) logmsg(WARN, "EFF stage returned false");
    }

    if(RUN_MASS) {
        TStopwatch t; t.Start();
        logmsg(INFO, "==== RUNNING MASS RECONSTRUCTION ====");
        bool ok = invariant_mass(runs, skim_out_dir, plots_dir, branch_file, cuts_file, "skim");
        t.Stop();
        logmsg(INFO, TString::Format("MASS stage done in %.2f s", t.RealTime()).Data());
        if(!ok) logmsg(WARN, "MASS stage returned false");
    }

    if(RUN_BKG) {
        TStopwatch t; t.Start();
        logmsg(INFO, "==== RUNNING BACKGROUND SUBTRACTION ====");
        bool ok = background_subtract(runs, plots_dir, cuts_file, "mass_run");
        t.Stop();
        logmsg(INFO, TString::Format("BKG stage done in %.2f s", t.RealTime()).Data());
        if(!ok) logmsg(WARN, "BKG stage returned false");
    }

    totalTimer.Stop();
    logmsg(INFO, TString::Format("main_analysis completed in %.2f s", totalTimer.RealTime()).Data());
}
