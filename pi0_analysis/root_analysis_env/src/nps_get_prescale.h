// nps_get_prescale.h
#pragma once

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iostream>
#include <regex>
#include <algorithm>
#include <cctype>

#include "TString.h"

namespace nps {
namespace detail {

// trim helpers
static inline std::string trim_copy(const std::string &s) {
    size_t a = 0;
    while (a < s.size() && std::isspace((unsigned char)s[a])) ++a;
    size_t b = s.size();
    while (b > a && std::isspace((unsigned char)s[b-1])) --b;
    return s.substr(a, b-a);
}

// robust CSV line splitter (handles quoted fields)
static inline std::vector<std::string> split_csv_line(const std::string &line) {
    std::vector<std::string> out;
    std::string cur;
    bool in_quotes = false;
    for (size_t i = 0; i < line.size(); ++i) {
        char c = line[i];
        if (c == '"') {
            if (in_quotes && i+1 < line.size() && line[i+1] == '"') {
                cur.push_back('"'); // escaped quote
                ++i;
            } else {
                in_quotes = !in_quotes;
            }
        } else if (c == ',' && !in_quotes) {
            out.push_back(trim_copy(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    out.push_back(trim_copy(cur));
    return out;
}

// Parse prescale token using regex "ps<l>=<r>" (e.g. "ps4=2").
// Returns true if matched and sets l_out and r_out.
static inline bool parse_ps_token_lr(const std::string &tok, int &l_out, int &r_out) {
    try {
        // Regex: optional whitespace, 'ps' (case-insensitive), digits (group1), optional spaces, '=', optional spaces, digits (group2)
        static const std::regex re(R"(\s*ps\s*(\d+)\s*=\s*(\d+)\s*)", std::regex::icase);
        std::smatch m;
        if (std::regex_search(tok, m, re) && m.size() >= 3) {
            l_out = std::stoi(m[1].str());
            r_out = std::stoi(m[2].str());
            return true;
        }
    } catch (const std::regex_error &e) {
        std::cerr << "[get_prescale] regex_error in parse_ps_token_lr: " << e.what() << "  token='" << tok << "'\n";
    } catch (...) {
        // fall through
    }
    return false;
}

// fallback: parse first integer in token
static inline bool parse_first_int(const std::string &tok, int &out) {
    std::string digits;
    for (char c : tok) if (std::isdigit((unsigned char)c)) digits.push_back(c);
    if (digits.empty()) return false;
    try { out = std::stoi(digits); return true; } catch(...) { return false; }
}

static inline unsigned long long compute_prescale_from_r(int r) {
    if (r <= 0) return 1ULL;
    if (r > 63) {
        std::cerr << "[get_prescale] WARNING: r (" << r << ") > 63; capping to 63 to avoid overflow\n";
        r = 63;
    }
    // compute 2^(r-1) + 1
    return (1ULL << (r - 1)) + 1ULL;
}

// parse CSV input stream and compute prescale for run_number
static inline unsigned long long parse_csv_stream_for_prescale(std::istream &is, int run_number, const std::string &source_label) {
    std::string header;
    if (!std::getline(is, header)) {
        std::cerr << "[get_prescale] ERROR: empty CSV (" << source_label << ")\n";
        return 1ULL;
    }
    auto headers = split_csv_line(header);
    int idx_run = -1, idx_prescale = -1;
    for (size_t i = 0; i < headers.size(); ++i) {
        std::string h = headers[i];
        std::string hl = h;
        std::transform(hl.begin(), hl.end(), hl.begin(), ::tolower);
        if (hl == "run_number" || hl == "run") idx_run = (int)i;
        if (hl == "prescale") idx_prescale = (int)i;
    }
    if (idx_run < 0 || idx_prescale < 0) {
        std::cerr << "[get_prescale] ERROR: CSV (" << source_label << ") missing run_number or prescale column\n";
        return 1ULL;
    }

    std::string line;
    while (std::getline(is, line)) {
        if (line.empty()) continue;
        auto toks = split_csv_line(line);
        if ((int)toks.size() <= std::max(idx_run, idx_prescale)) continue;
        std::string runstr = trim_copy(toks[idx_run]);
        if (runstr.empty()) continue;
        int run = 0;
        try { run = std::stoi(runstr); } catch(...) { continue; }
        if (run != run_number) continue;

        // found matching run
        std::string ps_tok = trim_copy(toks[idx_prescale]);
        std::cout << "[get_prescale] run=" << run << " prescale_token='" << ps_tok << "'\n";

        int l=0, r=0;
        if (parse_ps_token_lr(ps_tok, l, r)) {
            std::cout << "[get_prescale] parsed token: ps" << l << "=" << r << "  (using right-hand value R=" << r << " to compute prescale)\n";
            unsigned long long val = compute_prescale_from_r(r);
            std::cout << "[get_prescale] computed prescale = " << val << "\n";
            return val;
        }

        // fallback: try parse first integer as r (right-hand)
        if (parse_first_int(ps_tok, r)) {
            std::cout << "[get_prescale] fallback parsed integer r=" << r << " from token\n";
            unsigned long long val = compute_prescale_from_r(r);
            std::cout << "[get_prescale] computed prescale (fallback) = " << val << "\n";
            return val;
        }

        std::cout << "[get_prescale] WARNING: could not parse prescale token -> returning 1\n";
        return 1ULL;
    }

    std::cout << "[get_prescale] run " << run_number << " not found in CSV (" << source_label << ") -> returning 1\n";
    return 1ULL;
}

} // namespace detail

// -----------------------------------------------------------------------------
// Public API
// -----------------------------------------------------------------------------

// Primary API: accepts runlist file + run number and optional explicit csvPath.
// If csvPath is non-empty and readable -> parse it directly.
// Otherwise fall back to scanning runlist for a referenced CSV and then to main.csv fallback.
static inline unsigned long long get_prescale_from_runlist(const std::string &runlistFilePath,
                                                           int run_number,
                                                           const std::string &csvPath = std::string())
{
    // if explicit csvPath provided, try it first
    if (!csvPath.empty()) {
        std::ifstream csvifs(csvPath);
        if (csvifs.is_open()) {
            std::cout << "[get_prescale] using provided CSV: " << csvPath << "\n";
            return detail::parse_csv_stream_for_prescale(csvifs, run_number, csvPath);
        } else {
            std::cerr << "[get_prescale] WARNING: provided CSV '" << csvPath << "' cannot be opened. Falling back.\n";
        }
    }

    // try to open runlist file
    std::ifstream ifs(runlistFilePath);
    if (!ifs.is_open()) {
        std::cerr << "[get_prescale] ERROR: cannot open runlistFile '" << runlistFilePath << "'\n";
        return 1ULL;
    }

    // peek first non-empty line
    std::string first_line;
    std::streampos start_pos = ifs.tellg();
    while (std::getline(ifs, first_line)) {
        if (!first_line.empty()) break;
    }
    ifs.clear();
    ifs.seekg(0, std::ios::beg);

    std::string lower = first_line;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

    auto contains = [&](const std::string &hay, const std::string &needle){ return hay.find(needle) != std::string::npos; };

    // if file itself is CSV
    if (!first_line.empty() && contains(first_line, ",") && contains(lower, "prescale")) {
        std::cout << "[get_prescale] treating '" << runlistFilePath << "' as CSV\n";
        return detail::parse_csv_stream_for_prescale(ifs, run_number, runlistFilePath);
    }

    // otherwise, scan the runlist for a .csv path
    std::string line;
    std::string found_csv;
    while (std::getline(ifs, line)) {
        if (line.empty()) continue;
        if (line.find(".csv") != std::string::npos) {
            size_t pos = line.find(".csv");
            size_t left = pos;
            while (left > 0 && line[left-1] != '"' && line[left-1] != '\'' && line[left-1] != ' ' && line[left-1] != '\t') --left;
            size_t right = pos + 4;
            while (right < line.size() && line[right] != '"' && line[right] != '\'' && line[right] != ' ' && line[right] != '\t' && line[right] != ',') ++right;
            std::string candidate = detail::trim_copy(line.substr(left, right-left));
            if (!candidate.empty() && (candidate.front() == '"' || candidate.front() == '\'')) candidate = candidate.substr(1);
            if (!candidate.empty() && (candidate.back() == '"' || candidate.back() == '\'')) candidate.pop_back();
            if (!candidate.empty()) { found_csv = candidate; break; }
        }
    }

    if (!found_csv.empty()) {
        std::cout << "[get_prescale] found CSV reference in runlist: '" << found_csv << "'\n";
        std::ifstream csvifs(found_csv);
        if (csvifs.is_open()) return detail::parse_csv_stream_for_prescale(csvifs, run_number, found_csv);
        std::cerr << "[get_prescale] WARNING: referenced CSV '" << found_csv << "' cannot be opened\n";
    }

    // fallback: main.csv in same directory
    size_t slash = runlistFilePath.find_last_of("/\\");
    std::string dir = (slash == std::string::npos) ? std::string(".") : runlistFilePath.substr(0, slash);
    std::string main_csv = dir + "/main.csv";
    std::ifstream mainifs(main_csv);
    if (mainifs.is_open()) {
        std::cout << "[get_prescale] using fallback CSV: " << main_csv << "\n";
        return detail::parse_csv_stream_for_prescale(mainifs, run_number, main_csv);
    }

    std::cout << "[get_prescale] no CSV found; returning default prescale=1\n";
    return 1ULL;
}

// TString overload
static inline unsigned long long get_prescale_from_runlist(const TString &runlistFilePath,
                                                           int run_number,
                                                           const TString &csvPath = TString())
{
    return get_prescale_from_runlist(std::string(runlistFilePath.Data()), run_number,
                                     std::string(csvPath.Data()));
}

} // namespace nps
