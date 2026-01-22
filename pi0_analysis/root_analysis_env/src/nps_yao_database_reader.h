#pragma once
/* nps__yao_database_reader.h
   Header-only utility to read config/DataBase_production_runs_newBCMOffset.txt
   and provide CPU_LT, Ave_BeamCurr, and Charge_tot for a given run number.

   Usage example:
     #include "nps__yao_database_reader.h"

     double charge = nps::getChargeTot_or_nan(1572);  // returns uC or NaN
     double cpu_lt = nps::getCPU_LT_or_nan(1572);     // returns CPU_LT or NaN
     double ave_curr = nps::getAve_BeamCurr_or_nan(1572); // returns µA or NaN
*/

#include <string>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <limits>

namespace nps {

namespace detail {

inline std::vector<std::string> split_ws(const std::string &line) {
    std::istringstream ss(line);
    std::vector<std::string> toks;
    std::string tok;
    while (ss >> tok) toks.push_back(tok);
    return toks;
}

inline bool parse_column_map(const std::string &path,
                             const std::string &column_name,
                             std::unordered_map<int,double> &out_map)
{
    std::ifstream ifs(path);
    if (!ifs.is_open()) {
        std::cerr << "nps_database_reader: cannot open " << path << std::endl;
        return false;
    }

    std::string line;
    bool header_found = false;
    int idx_run = -1, idx_col = -1;

    while (std::getline(ifs, line)) {
        bool allws = true;
        for (char c : line) if (!isspace((unsigned char)c)) { allws = false; break; }
        if (allws) continue;

        std::vector<std::string> toks = split_ws(line);
        if (toks.empty()) continue;

        if (!header_found) {
            if (toks[0] == "RunNo") {
                header_found = true;
                for (size_t i = 0; i < toks.size(); ++i) {
                    if (toks[i] == "RunNo") idx_run = static_cast<int>(i);
                    if (toks[i] == column_name) idx_col = static_cast<int>(i);
                }
                if (idx_run == -1 || idx_col == -1) {
                    std::cerr << "nps_database_reader: header found but column '"
                              << column_name << "' not present.\n";
                    return false;
                }
            }
            continue;
        } else {
            if (toks.size() <= static_cast<size_t>(std::max(idx_run, idx_col))) continue;
            try {
                int run = std::stoi(toks[idx_run]);
                double val = std::stod(toks[idx_col]);
                out_map[run] = val;
            } catch (...) {
                continue;
            }
        }
    }

    return true;
}

inline const std::string kDBPath = "config/DataBase_production_runs_newBCMOffset.txt";

} // namespace detail

// --- Public API: only _or_nan functions ---

inline double getChargeTot_or_nan(int run) {
    static std::unordered_map<int,double> s_map;
    static bool loaded = detail::parse_column_map(detail::kDBPath, "Charge_tot", s_map);
    auto it = s_map.find(run);
    return it != s_map.end() ? it->second : std::numeric_limits<double>::quiet_NaN();
}

inline double getCPU_LT_or_nan(int run) {
    static std::unordered_map<int,double> s_map;
    static bool loaded = detail::parse_column_map(detail::kDBPath, "CPU_LT", s_map);
    auto it = s_map.find(run);
    return it != s_map.end() ? it->second : std::numeric_limits<double>::quiet_NaN();
}

inline double getAve_BeamCurr_or_nan(int run) {
    static std::unordered_map<int,double> s_map;
    static bool loaded = detail::parse_column_map(detail::kDBPath, "Ave_BeamCurr", s_map);
    auto it = s_map.find(run);
    return it != s_map.end() ? it->second : std::numeric_limits<double>::quiet_NaN();
}

inline double getBeam_Time_or_nan(int run) {
    static std::unordered_map<int,double> s_map;
    static bool loaded = detail::parse_column_map(detail::kDBPath, "Beam_Time", s_map);
    auto it = s_map.find(run);
    return it != s_map.end() ? it->second : std::numeric_limits<double>::quiet_NaN();
}

} // namespace nps
