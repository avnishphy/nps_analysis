#pragma once
#include <fstream>
#include <stdexcept>
#include <string>
#include <cstdio>

// Written before processing and atomically replaced only after all outputs close.
// A stale ROOT file cannot override a failed or interrupted attempt.
inline void write_nps_run_status(const std::string& path, int run,
        bool success, const std::string& stage, bool fit_valid=false,
        int minimizer_status=-1, int covariance_status=-1,
        std::string reason="", int attempts=0, bool at_boundary=false,
        const std::string& classification="interior_valid") {
    for (auto& c : reason) if (c==',' || c=='\n' || c=='\r') c=';';
    const std::string tmp=path+".tmp";
    std::ofstream out(tmp);
    out << "run,success,stage,fit_valid,minimizer_status,covariance_status,reason,attempts,at_boundary,classification\n"
        << run << ',' << success << ',' << stage << ',' << fit_valid << ','
        << minimizer_status << ',' << covariance_status << ',' << reason << ','
        << attempts << ',' << at_boundary << ',' << classification << '\n';
    out.close();
    if (!out || std::rename(tmp.c_str(),path.c_str()) != 0)
        throw std::runtime_error("Cannot publish run status: "+path);
}
