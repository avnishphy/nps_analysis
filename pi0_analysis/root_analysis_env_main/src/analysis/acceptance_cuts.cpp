// acceptance_cuts.cpp
// Loads and applies acceptance cuts from INI-style config sections.

#include "acceptance_cuts.h"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {

using SectionMap = std::map<std::string, std::map<std::string, std::string>>;

std::string trim_copy(const std::string& s) {
    size_t first = 0;
    while (first < s.size() && std::isspace(static_cast<unsigned char>(s[first]))) {
        ++first;
    }

    size_t last = s.size();
    while (last > first && std::isspace(static_cast<unsigned char>(s[last - 1]))) {
        --last;
    }

    return s.substr(first, last - first);
}

std::string to_lower_copy(const std::string& s) {
    std::string out = s;
    std::transform(out.begin(), out.end(), out.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return out;
}

std::string strip_inline_comment(const std::string& value) {
    bool in_single_quote = false;
    bool in_double_quote = false;

    for (size_t i = 0; i < value.size(); ++i) {
        const char c = value[i];
        if (c == '\'' && !in_double_quote) {
            in_single_quote = !in_single_quote;
            continue;
        }
        if (c == '"' && !in_single_quote) {
            in_double_quote = !in_double_quote;
            continue;
        }

        if (!in_single_quote && !in_double_quote && (c == '#' || c == ';')) {
            if (i == 0 || std::isspace(static_cast<unsigned char>(value[i - 1]))) {
                return trim_copy(value.substr(0, i));
            }
        }
    }

    return trim_copy(value);
}

std::vector<std::string> split_csv_tokens(const std::string& value) {
    std::vector<std::string> out;
    std::stringstream ss(value);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        tok = trim_copy(tok);
        if (!tok.empty()) {
            out.push_back(tok);
        }
    }
    return out;
}

std::string bool_to_text(bool v) {
    return v ? "true" : "false";
}

bool parse_bool_strict(const std::string& raw,
                       const std::string& section,
                       const std::string& key,
                       const std::string& cfg_file) {
    const std::string v = to_lower_copy(trim_copy(raw));
    if (v == "1" || v == "true" || v == "yes" || v == "on") return true;
    if (v == "0" || v == "false" || v == "no" || v == "off") return false;

    throw std::runtime_error(
        "Invalid boolean value for key '" + key + "' in section [" + section +
        "] of " + cfg_file + ": '" + raw + "'");
}

double parse_double_strict(const std::string& raw,
                           const std::string& section,
                           const std::string& key,
                           const std::string& cfg_file) {
    try {
        size_t pos = 0;
        const double val = std::stod(raw, &pos);
        if (pos != raw.size()) {
            throw std::runtime_error("trailing characters");
        }
        return val;
    } catch (...) {
        throw std::runtime_error(
            "Invalid numeric value for key '" + key + "' in section [" + section +
            "] of " + cfg_file + ": '" + raw + "'");
    }
}

const std::string& require_key(const std::map<std::string, std::string>& merged,
                               const std::string& section,
                               const std::string& key,
                               const std::string& cfg_file) {
    auto it = merged.find(to_lower_copy(key));
    if (it == merged.end()) {
        throw std::runtime_error(
            "Missing required key '" + key + "' in effective section [" + section +
            "] of " + cfg_file);
    }
    return it->second;
}

void require_range(double lo,
                   double hi,
                   const std::string& label,
                   const std::string& cfg_file,
                   const std::string& section) {
    if (!(lo < hi)) {
        throw std::runtime_error(
            "Invalid range for " + label + " in section [" + section + "] of " +
            cfg_file + ": min must be < max");
    }
}

void require_positive(double value,
                      const std::string& label,
                      const std::string& cfg_file,
                      const std::string& section) {
    if (!(value > 0.0)) {
        throw std::runtime_error(
            "Invalid value for " + label + " in section [" + section + "] of " +
            cfg_file + ": must be > 0");
    }
}

bool is_auto_token(const std::string& raw) {
    return to_lower_copy(trim_copy(raw)) == "auto";
}

SectionMap parse_config_file(const std::string& config_file) {
    std::ifstream in(config_file.c_str());
    if (!in.is_open()) {
        throw std::runtime_error("Cannot open acceptance config file: " + config_file);
    }

    SectionMap sections;
    std::string current_section;
    std::string line;
    int line_no = 0;

    while (std::getline(in, line)) {
        ++line_no;
        const std::string stripped = trim_copy(line);
        if (stripped.empty()) continue;
        if (stripped[0] == '#' || stripped[0] == ';') continue;

        if (stripped.front() == '[' && stripped.back() == ']') {
            current_section = trim_copy(stripped.substr(1, stripped.size() - 2));
            if (current_section.empty()) {
                throw std::runtime_error("Empty section name at " + config_file +
                                         ":" + std::to_string(line_no));
            }
            sections[current_section];
            continue;
        }

        if (current_section.empty()) {
            throw std::runtime_error(
                "Key/value line before any section at " + config_file + ":" +
                std::to_string(line_no));
        }

        const size_t eq_pos = stripped.find('=');
        if (eq_pos == std::string::npos) {
            throw std::runtime_error("Expected key=value at " + config_file + ":" +
                                     std::to_string(line_no));
        }

        std::string key = trim_copy(stripped.substr(0, eq_pos));
        std::string value = strip_inline_comment(stripped.substr(eq_pos + 1));
        if (key.empty()) {
            throw std::runtime_error("Empty key at " + config_file + ":" +
                                     std::to_string(line_no));
        }
        if (value.empty()) {
            throw std::runtime_error("Empty value for key '" + key + "' at " +
                                     config_file + ":" + std::to_string(line_no));
        }

        sections[current_section][to_lower_copy(key)] = value;
    }

    return sections;
}

void merge_section_recursive(const SectionMap& sections,
                             const std::string& section,
                             std::map<std::string, std::string>& merged,
                             std::vector<std::string>& stack,
                             const std::string& config_file) {
    const auto it_section = sections.find(section);
    if (it_section == sections.end()) {
        throw std::runtime_error("Section [" + section + "] was requested but not found in " +
                                 config_file);
    }

    if (std::find(stack.begin(), stack.end(), section) != stack.end()) {
        std::ostringstream oss;
        oss << "Inheritance cycle detected in " << config_file << ": ";
        for (const auto& s : stack) oss << s << " -> ";
        oss << section;
        throw std::runtime_error(oss.str());
    }

    stack.push_back(section);

    const auto it_inherit = it_section->second.find("inherit");
    if (it_inherit != it_section->second.end()) {
        const std::vector<std::string> parents = split_csv_tokens(it_inherit->second);
        for (const auto& parent : parents) {
            merge_section_recursive(sections, parent, merged, stack, config_file);
        }
    }

    for (const auto& kv : it_section->second) {
        if (kv.first == "inherit") continue;
        merged[kv.first] = kv.second;
    }

    stack.pop_back();
}

std::string list_available_kin_sections(const SectionMap& sections) {
    std::vector<std::string> names;
    names.reserve(sections.size());
    for (const auto& kv : sections) {
        if (to_lower_copy(kv.first) == "global") continue;
        names.push_back(kv.first);
    }
    std::sort(names.begin(), names.end());

    std::ostringstream oss;
    for (size_t i = 0; i < names.size(); ++i) {
        if (i) oss << ", ";
        oss << names[i];
    }
    return oss.str();
}

} // namespace

AcceptanceCuts::AcceptanceCuts(const std::string& config_file,
                               const std::string& kin_setting)
    : source_file_(config_file), kin_setting_(kin_setting) {
    if (trim_copy(kin_setting).empty()) {
        throw std::runtime_error("AcceptanceCuts requires a non-empty kinematic setting");
    }

    const SectionMap sections = parse_config_file(config_file);
    if (sections.empty()) {
        throw std::runtime_error("Acceptance config has no sections: " + config_file);
    }

    if (sections.find(kin_setting) == sections.end()) {
        throw std::runtime_error(
            "Acceptance cuts section [" + kin_setting + "] not found in " +
            config_file + ". Available kin sections: " + list_available_kin_sections(sections));
    }

    std::map<std::string, std::string> merged;
    std::vector<std::string> stack;

    const auto it_global = sections.find("global");
    if (it_global != sections.end()) {
        merge_section_recursive(sections, "global", merged, stack, config_file);
    }
    merge_section_recursive(sections, kin_setting, merged, stack, config_file);

    hms_data_.edtm_tdc_max = parse_double_strict(require_key(merged, kin_setting, "hms_data.edtm_tdc_max", config_file),
                                                 kin_setting, "hms_data.edtm_tdc_max", config_file);
    hms_data_.react_z_min = parse_double_strict(require_key(merged, kin_setting, "hms_data.react_z_min", config_file),
                                                kin_setting, "hms_data.react_z_min", config_file);
    hms_data_.react_z_max = parse_double_strict(require_key(merged, kin_setting, "hms_data.react_z_max", config_file),
                                                kin_setting, "hms_data.react_z_max", config_file);
    hms_data_.delta_min = parse_double_strict(require_key(merged, kin_setting, "hms_data.delta_min", config_file),
                                              kin_setting, "hms_data.delta_min", config_file);
    hms_data_.delta_max = parse_double_strict(require_key(merged, kin_setting, "hms_data.delta_max", config_file),
                                              kin_setting, "hms_data.delta_max", config_file);
    hms_data_.gtr_th_min = parse_double_strict(require_key(merged, kin_setting, "hms_data.gtr_th_min", config_file),
                                               kin_setting, "hms_data.gtr_th_min", config_file);
    hms_data_.gtr_th_max = parse_double_strict(require_key(merged, kin_setting, "hms_data.gtr_th_max", config_file),
                                               kin_setting, "hms_data.gtr_th_max", config_file);
    hms_data_.gtr_ph_min = parse_double_strict(require_key(merged, kin_setting, "hms_data.gtr_ph_min", config_file),
                                               kin_setting, "hms_data.gtr_ph_min", config_file);
    hms_data_.gtr_ph_max = parse_double_strict(require_key(merged, kin_setting, "hms_data.gtr_ph_max", config_file),
                                               kin_setting, "hms_data.gtr_ph_max", config_file);
    hms_data_.cer_npe_sum_min = parse_double_strict(require_key(merged, kin_setting, "hms_data.cer_npe_sum_min", config_file),
                                                    kin_setting, "hms_data.cer_npe_sum_min", config_file);
    hms_data_.cal_etotnorm_min = parse_double_strict(require_key(merged, kin_setting, "hms_data.cal_etotnorm_min", config_file),
                                                     kin_setting, "hms_data.cal_etotnorm_min", config_file);
    hms_data_.cal_etotnorm_max = parse_double_strict(require_key(merged, kin_setting, "hms_data.cal_etotnorm_max", config_file),
                                                     kin_setting, "hms_data.cal_etotnorm_max", config_file);

    hms_simulation_.react_z_min = parse_double_strict(require_key(merged, kin_setting, "hms_sim.react_z_min", config_file),
                                                      kin_setting, "hms_sim.react_z_min", config_file);
    hms_simulation_.react_z_max = parse_double_strict(require_key(merged, kin_setting, "hms_sim.react_z_max", config_file),
                                                      kin_setting, "hms_sim.react_z_max", config_file);
    hms_simulation_.delta_min = parse_double_strict(require_key(merged, kin_setting, "hms_sim.delta_min", config_file),
                                                    kin_setting, "hms_sim.delta_min", config_file);
    hms_simulation_.delta_max = parse_double_strict(require_key(merged, kin_setting, "hms_sim.delta_max", config_file),
                                                    kin_setting, "hms_sim.delta_max", config_file);
    hms_simulation_.gtr_th_min = parse_double_strict(require_key(merged, kin_setting, "hms_sim.gtr_th_min", config_file),
                                                     kin_setting, "hms_sim.gtr_th_min", config_file);
    hms_simulation_.gtr_th_max = parse_double_strict(require_key(merged, kin_setting, "hms_sim.gtr_th_max", config_file),
                                                     kin_setting, "hms_sim.gtr_th_max", config_file);
    hms_simulation_.gtr_ph_min = parse_double_strict(require_key(merged, kin_setting, "hms_sim.gtr_ph_min", config_file),
                                                     kin_setting, "hms_sim.gtr_ph_min", config_file);
    hms_simulation_.gtr_ph_max = parse_double_strict(require_key(merged, kin_setting, "hms_sim.gtr_ph_max", config_file),
                                                     kin_setting, "hms_sim.gtr_ph_max", config_file);

    nps_cluster_.energy_min = parse_double_strict(require_key(merged, kin_setting, "nps_cluster.energy_min", config_file),
                                                  kin_setting, "nps_cluster.energy_min", config_file);
    nps_cluster_.x_min = parse_double_strict(require_key(merged, kin_setting, "nps_cluster.x_min", config_file),
                                             kin_setting, "nps_cluster.x_min", config_file);
    nps_cluster_.x_max = parse_double_strict(require_key(merged, kin_setting, "nps_cluster.x_max", config_file),
                                             kin_setting, "nps_cluster.x_max", config_file);
    nps_cluster_.y_min = parse_double_strict(require_key(merged, kin_setting, "nps_cluster.y_min", config_file),
                                             kin_setting, "nps_cluster.y_min", config_file);
    nps_cluster_.y_max = parse_double_strict(require_key(merged, kin_setting, "nps_cluster.y_max", config_file),
                                             kin_setting, "nps_cluster.y_max", config_file);
    nps_cluster_.time_center_ns = parse_double_strict(require_key(merged, kin_setting, "nps_cluster.time_center_ns", config_file),
                                                      kin_setting, "nps_cluster.time_center_ns", config_file);

    {
        const std::string raw = require_key(merged, kin_setting, "nps_cluster.time_halfwidth_ns", config_file);
        if (is_auto_token(raw)) {
            nps_cluster_.has_time_halfwidth_override = false;
        } else {
            nps_cluster_.has_time_halfwidth_override = true;
            nps_cluster_.time_halfwidth_ns_override =
                parse_double_strict(raw, kin_setting, "nps_cluster.time_halfwidth_ns", config_file);
        }
    }

    {
        const std::string raw = require_key(merged, kin_setting, "pairing.time_diff_max_ns", config_file);
        if (is_auto_token(raw)) {
            pairing_.has_time_diff_max_override = false;
        } else {
            pairing_.has_time_diff_max_override = true;
            pairing_.time_diff_max_ns_override =
                parse_double_strict(raw, kin_setting, "pairing.time_diff_max_ns", config_file);
        }
    }

    {
        const std::string raw = require_key(merged, kin_setting, "timing.use_shifted_sidebands", config_file);
        if (is_auto_token(raw)) {
            timing_.has_shifted_sidebands_override = false;
        } else {
            timing_.has_shifted_sidebands_override = true;
            timing_.use_shifted_sidebands_override =
                parse_bool_strict(raw, kin_setting, "timing.use_shifted_sidebands", config_file);
        }
    }

    weighted_.mmiss_corr_min = parse_double_strict(require_key(merged, kin_setting, "weighted.mmiss_corr_min", config_file),
                                                   kin_setting, "weighted.mmiss_corr_min", config_file);
    weighted_.mmiss_corr_max = parse_double_strict(require_key(merged, kin_setting, "weighted.mmiss_corr_max", config_file),
                                                   kin_setting, "weighted.mmiss_corr_max", config_file);

    require_range(hms_data_.react_z_min, hms_data_.react_z_max,
                  "hms_data.react_z", config_file, kin_setting);
    require_range(hms_data_.delta_min, hms_data_.delta_max,
                  "hms_data.delta", config_file, kin_setting);
    require_range(hms_data_.gtr_th_min, hms_data_.gtr_th_max,
                  "hms_data.gtr_th", config_file, kin_setting);
    require_range(hms_data_.gtr_ph_min, hms_data_.gtr_ph_max,
                  "hms_data.gtr_ph", config_file, kin_setting);
    require_range(hms_data_.cal_etotnorm_min, hms_data_.cal_etotnorm_max,
                  "hms_data.cal_etotnorm", config_file, kin_setting);

    require_range(hms_simulation_.react_z_min, hms_simulation_.react_z_max,
                  "hms_sim.react_z", config_file, kin_setting);
    require_range(hms_simulation_.delta_min, hms_simulation_.delta_max,
                  "hms_sim.delta", config_file, kin_setting);
    require_range(hms_simulation_.gtr_th_min, hms_simulation_.gtr_th_max,
                  "hms_sim.gtr_th", config_file, kin_setting);
    require_range(hms_simulation_.gtr_ph_min, hms_simulation_.gtr_ph_max,
                  "hms_sim.gtr_ph", config_file, kin_setting);

    require_range(nps_cluster_.x_min, nps_cluster_.x_max,
                  "nps_cluster.x", config_file, kin_setting);
    require_range(nps_cluster_.y_min, nps_cluster_.y_max,
                  "nps_cluster.y", config_file, kin_setting);
    require_positive(nps_cluster_.energy_min,
                     "nps_cluster.energy_min", config_file, kin_setting);
    if (nps_cluster_.has_time_halfwidth_override) {
        require_positive(nps_cluster_.time_halfwidth_ns_override,
                         "nps_cluster.time_halfwidth_ns", config_file, kin_setting);
    }
    if (pairing_.has_time_diff_max_override) {
        require_positive(pairing_.time_diff_max_ns_override,
                         "pairing.time_diff_max_ns", config_file, kin_setting);
    }

    require_range(weighted_.mmiss_corr_min, weighted_.mmiss_corr_max,
                  "weighted.mmiss_corr", config_file, kin_setting);
}

double AcceptanceCuts::resolved_cluster_time_halfwidth_ns(double fallback_halfwidth_ns) const noexcept {
    return nps_cluster_.has_time_halfwidth_override
               ? nps_cluster_.time_halfwidth_ns_override
               : fallback_halfwidth_ns;
}

double AcceptanceCuts::resolved_pair_time_diff_max_ns(double fallback_pair_time_diff_ns) const noexcept {
    return pairing_.has_time_diff_max_override
               ? pairing_.time_diff_max_ns_override
               : fallback_pair_time_diff_ns;
}

bool AcceptanceCuts::resolved_use_shifted_sidebands(bool fallback_use_shifted) const noexcept {
    return timing_.has_shifted_sidebands_override
               ? timing_.use_shifted_sidebands_override
               : fallback_use_shifted;
}

bool AcceptanceCuts::pass_hms_data(double edtm_tdc,
                                   double h_delta,
                                   double h_gtr_th,
                                   double h_gtr_ph,
                                   double h_cer_npe_sum,
                                   double h_cal_etotnorm,
                                   double h_react_z) const noexcept {
    if (edtm_tdc > hms_data_.edtm_tdc_max) return false;
    if (h_react_z < hms_data_.react_z_min || h_react_z > hms_data_.react_z_max) return false;
    if (h_delta < hms_data_.delta_min || h_delta > hms_data_.delta_max) return false;
    if (h_gtr_th < hms_data_.gtr_th_min || h_gtr_th > hms_data_.gtr_th_max) return false;
    if (h_gtr_ph < hms_data_.gtr_ph_min || h_gtr_ph > hms_data_.gtr_ph_max) return false;
    if (h_cer_npe_sum < hms_data_.cer_npe_sum_min) return false;
    if (h_cal_etotnorm < hms_data_.cal_etotnorm_min || h_cal_etotnorm > hms_data_.cal_etotnorm_max) return false;
    return true;
}

bool AcceptanceCuts::pass_hms_simulation(double h_delta,
                                         double h_gtr_th,
                                         double h_gtr_ph,
                                         double h_react_z) const noexcept {
    if (h_react_z < hms_simulation_.react_z_min || h_react_z > hms_simulation_.react_z_max) return false;
    if (h_delta < hms_simulation_.delta_min || h_delta > hms_simulation_.delta_max) return false;
    if (h_gtr_th < hms_simulation_.gtr_th_min || h_gtr_th > hms_simulation_.gtr_th_max) return false;
    if (h_gtr_ph < hms_simulation_.gtr_ph_min || h_gtr_ph > hms_simulation_.gtr_ph_max) return false;
    return true;
}

bool AcceptanceCuts::pass_nps_cluster(double clusE,
                                      double clusX,
                                      double clusY,
                                      double clusT,
                                      double resolved_time_halfwidth_ns) const noexcept {
    if (clusE <= nps_cluster_.energy_min) return false;
    if (!(clusX > nps_cluster_.x_min && clusX < nps_cluster_.x_max)) return false;
    if (!(clusY > nps_cluster_.y_min && clusY < nps_cluster_.y_max)) return false;

    const double t_lo = nps_cluster_.time_center_ns - resolved_time_halfwidth_ns;
    const double t_hi = nps_cluster_.time_center_ns + resolved_time_halfwidth_ns;
    if (!(clusT > t_lo && clusT < t_hi)) return false;

    return true;
}

bool AcceptanceCuts::pass_weighted_exclusive(double mmiss_corr) const noexcept {
    return (mmiss_corr >= weighted_.mmiss_corr_min && mmiss_corr <= weighted_.mmiss_corr_max);
}

std::string AcceptanceCuts::summary() const {
    std::ostringstream oss;
    oss << "Acceptance cuts loaded from " << source_file_ << "\n";
    oss << "  kin_setting: " << kin_setting_ << "\n";

    oss << "  hms_data.edtm_tdc_max=" << hms_data_.edtm_tdc_max << "\n";
    oss << "  hms_data.react_z=[" << hms_data_.react_z_min << ", " << hms_data_.react_z_max << "]\n";
    oss << "  hms_data.delta=[" << hms_data_.delta_min << ", " << hms_data_.delta_max << "]\n";
    oss << "  hms_data.gtr_th=[" << hms_data_.gtr_th_min << ", " << hms_data_.gtr_th_max << "]\n";
    oss << "  hms_data.gtr_ph=[" << hms_data_.gtr_ph_min << ", " << hms_data_.gtr_ph_max << "]\n";
    oss << "  hms_data.cer_npe_sum_min=" << hms_data_.cer_npe_sum_min << "\n";
    oss << "  hms_data.cal_etotnorm=[" << hms_data_.cal_etotnorm_min << ", " << hms_data_.cal_etotnorm_max << "]\n";

    oss << "  hms_sim.react_z=[" << hms_simulation_.react_z_min << ", " << hms_simulation_.react_z_max << "]\n";
    oss << "  hms_sim.delta=[" << hms_simulation_.delta_min << ", " << hms_simulation_.delta_max << "]\n";
    oss << "  hms_sim.gtr_th=[" << hms_simulation_.gtr_th_min << ", " << hms_simulation_.gtr_th_max << "]\n";
    oss << "  hms_sim.gtr_ph=[" << hms_simulation_.gtr_ph_min << ", " << hms_simulation_.gtr_ph_max << "]\n";

    oss << "  nps_cluster.energy_min=" << nps_cluster_.energy_min << "\n";
    oss << "  nps_cluster.x=[" << nps_cluster_.x_min << ", " << nps_cluster_.x_max << "]\n";
    oss << "  nps_cluster.y=[" << nps_cluster_.y_min << ", " << nps_cluster_.y_max << "]\n";
    oss << "  nps_cluster.time_center_ns=" << nps_cluster_.time_center_ns << "\n";
    oss << "  nps_cluster.time_halfwidth_ns="
        << (nps_cluster_.has_time_halfwidth_override ? std::to_string(nps_cluster_.time_halfwidth_ns_override) : "auto")
        << "\n";
    oss << "  pairing.time_diff_max_ns="
        << (pairing_.has_time_diff_max_override ? std::to_string(pairing_.time_diff_max_ns_override) : "auto")
        << "\n";
    oss << "  timing.use_shifted_sidebands="
        << (timing_.has_shifted_sidebands_override ? bool_to_text(timing_.use_shifted_sidebands_override) : "auto")
        << "\n";

    oss << "  weighted.mmiss_corr=[" << weighted_.mmiss_corr_min << ", " << weighted_.mmiss_corr_max << "]";
    return oss.str();
}
