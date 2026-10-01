// this code computes the boiling correction for the NPS pi0 analysis. It is based on the code in boil_corr.cxx, but has been modified to work with the NPS pi0 analysis framework. The code reads in the data and simulation trees, applies the necessary cuts, and computes the boiling correction as a function of beam current. The output is a ROOT file containing the boiling correction histogram and fit parameters.
// it selects the good scalers thorugh the use of Christine's good event selection, uses those events to compute the EDTM livetimes, and finally the yields vs current. The carbon target is used to compute the boiling correction, which is then applied to the LH2 and LD2 targets. The output is a ROOT file containing the boiling correction histogram and fit parameters.

#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TStyle.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <vector>       

namespace fs = std::filesystem;

struct ConfigRow{
    std::string prescale_token;
    double ps_val;
    std::string target;
    double edtm_live_time;
};

static std::map<int, ConfigRow> build_config_lookup(const std::string& cfg_path){
    std::ifstream in(cfg_path);
    if(!in) throw std::runtime_error("Failed to open config file: " + cfg_path);
    std::string line;
    //read header
    if(!std::getline(in, line)) throw std::runtime_error("Config csv is empty: " + cfg_path);
    std::istringstream iss_header(line);
    std::vector<std::string> header;
    std::string col;
    while(std::getline(iss_header, col, ',')) header.push_back(col);
    // find relevant columns
    int ps_col = find_prescale_column(header);
    int run_col = -1, target_col = -1, edtm_live_col = -1;
    for(size_t i = 0; i < header.size(); ++i){
        std:: string h = header[i];
        std::transform(h.begin(), h.end(), h.begin(), ::tolower);
        if (run_col == -1 && h.find("run_number") != std::string::npos) run_col = i;
        if (target_col == -1 && h.find("target") != std::string::npos) target_col = i;
        if (edtm_live_col == -1 && h.find("NewGen_EDTM_livetime") != std::string::npos) edtm_live_col = i; 
}

if(run_col == -1) throw std::runtime_error("Could not find run number column in config csv.");
// parse rows
std::map<int, ConfigRow> lookup;
while (std::getline(in, line)){
    std::istringstream iss(line);
    
}