#pragma once

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "eff_math_utils.h"
#include "efficiency_types.h"

namespace effstuff {

struct ConfigCsvData {
  std::vector<RunConfigRow> rows;
  int malformed_row_count = 0;
  std::vector<std::string> malformed_messages;
};

inline std::vector<std::string> parse_csv_line(const std::string& line) {
  std::vector<std::string> fields;
  std::string current;
  bool in_quotes = false;

  for (size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];

    if (c == '"') {
      if (in_quotes && (i + 1) < line.size() && line[i + 1] == '"') {
        current.push_back('"');
        ++i;
      } else {
        in_quotes = !in_quotes;
      }
      continue;
    }

    if (c == ',' && !in_quotes) {
      fields.push_back(current);
      current.clear();
      continue;
    }

    current.push_back(c);
  }

  fields.push_back(current);
  return fields;
}

inline int find_header_index(const std::vector<std::string>& headers,
                             const std::string& needle) {
  for (size_t i = 0; i < headers.size(); ++i) {
    if (iequals(trim_copy(headers[i]), needle)) {
      return static_cast<int>(i);
    }
  }
  return -1;
}

inline bool parse_int_strict(const std::string& text, int& value_out) {
  const std::string s = trim_copy(text);
  if (s.empty()) {
    return false;
  }

  char* end_ptr = nullptr;
  const long parsed = std::strtol(s.c_str(), &end_ptr, 10);
  if (end_ptr == s.c_str() || *end_ptr != '\0') {
    return false;
  }

  value_out = static_cast<int>(parsed);
  return true;
}

inline bool load_config_csv(const std::string& csv_path,
                            ConfigCsvData& out,
                            std::string& error_out) {
  std::ifstream in(csv_path);
  if (!in) {
    error_out = "Could not open config CSV: " + csv_path;
    return false;
  }

  std::string header_line;
  if (!std::getline(in, header_line)) {
    error_out = "Config CSV is empty: " + csv_path;
    return false;
  }

  std::vector<std::string> headers = parse_csv_line(header_line);
  const int run_idx = find_header_index(headers, "run_number");
  const int kin_idx = find_header_index(headers, "Kin_old");
  const int type_idx = find_header_index(headers, "Type");
  const int prescale_idx = find_header_index(headers, "prescale");

  if (run_idx < 0 || kin_idx < 0 || type_idx < 0 || prescale_idx < 0) {
    error_out = "Config CSV is missing one or more required headers: run_number, Kin_old, Type, prescale";
    return false;
  }

  const int max_required_index = std::max(std::max(run_idx, kin_idx), std::max(type_idx, prescale_idx));

  std::string line;
  int line_number = 1;
  while (std::getline(in, line)) {
    ++line_number;

    if (trim_copy(line).empty()) {
      continue;
    }

    const std::vector<std::string> fields = parse_csv_line(line);
    if (static_cast<int>(fields.size()) <= max_required_index) {
      ++out.malformed_row_count;
      out.malformed_messages.push_back("line " + std::to_string(line_number) + ": missing required columns");
      continue;
    }

    int run_number = 0;
    if (!parse_int_strict(fields[run_idx], run_number)) {
      ++out.malformed_row_count;
      out.malformed_messages.push_back("line " + std::to_string(line_number) + ": invalid run_number '" + trim_copy(fields[run_idx]) + "'");
      continue;
    }

    RunConfigRow row;
    row.run_number = run_number;
    row.kin_old = trim_copy(fields[kin_idx]);
    row.run_type = trim_copy(fields[type_idx]);
    row.prescale_token = trim_copy(fields[prescale_idx]);

    if (row.kin_old.empty()) {
      ++out.malformed_row_count;
      out.malformed_messages.push_back("line " + std::to_string(line_number) + ": empty Kin_old");
      continue;
    }

    out.rows.push_back(row);
  }

  if (out.rows.empty()) {
    error_out = "No valid rows were parsed from config CSV: " + csv_path;
    return false;
  }

  return true;
}

inline std::vector<std::string> collect_unique_kinematics(const std::vector<RunConfigRow>& rows) {
  std::set<std::string> unique;
  for (const auto& row : rows) {
    unique.insert(row.kin_old);
  }
  return std::vector<std::string>(unique.begin(), unique.end());
}

inline std::vector<std::string> split_comma_list(const std::string& text) {
  std::vector<std::string> out;
  std::stringstream ss(text);
  std::string token;
  while (std::getline(ss, token, ',')) {
    const std::string cleaned = trim_copy(token);
    if (!cleaned.empty()) {
      out.push_back(cleaned);
    }
  }
  return out;
}

inline std::vector<std::string> default_allowed_types() {
  return {"production", "Production"};
}

inline bool type_is_allowed(const std::string& run_type,
                            const std::vector<std::string>& allowed_types) {
  for (const auto& allowed : allowed_types) {
    if (run_type == allowed) {
      return true;
    }
  }
  return false;
}

inline std::string sanitize_for_filename(const std::string& text) {
  std::string out;
  out.reserve(text.size());
  for (unsigned char c : text) {
    if (std::isalnum(c) != 0) {
      out.push_back(static_cast<char>(c));
    } else {
      out.push_back('_');
    }
  }
  return out;
}

inline std::string csv_quote(const std::string& text) {
  std::string escaped;
  escaped.reserve(text.size() + 2);
  escaped.push_back('"');
  for (char c : text) {
    if (c == '"') {
      escaped.push_back('"');
      escaped.push_back('"');
    } else {
      escaped.push_back(c);
    }
  }
  escaped.push_back('"');
  return escaped;
}

}  // namespace effstuff
