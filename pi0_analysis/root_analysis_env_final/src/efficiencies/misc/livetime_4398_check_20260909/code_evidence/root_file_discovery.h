#pragma once

#include <algorithm>
#include <filesystem>
#include <string>
#include <vector>

#include "efficiency_types.h"

namespace effstuff {

namespace fs = std::filesystem;

inline bool starts_with(const std::string& text, const std::string& prefix) {
  return text.rfind(prefix, 0) == 0;
}

inline bool ends_with(const std::string& text, const std::string& suffix) {
  return text.size() >= suffix.size() &&
         text.compare(text.size() - suffix.size(), suffix.size(), suffix) == 0;
}

inline int parse_segment_from_filename(const std::string& filepath, int run_number) {
  const std::string name = fs::path(filepath).filename().string();
  const std::string prefix = "nps_hms_coin_" + std::to_string(run_number) + "_";
  const std::string suffix = "_1_-1.root";

  if (!starts_with(name, prefix) || !ends_with(name, suffix)) {
    return -1;
  }

  const size_t start = prefix.size();
  const size_t end = name.size() - suffix.size();
  if (end <= start) {
    return -1;
  }

  try {
    return std::stoi(name.substr(start, end - start));
  } catch (...) {
    return -1;
  }
}

inline std::vector<std::string> find_run_files_in_directory(const std::string& directory,
                                                            int run_number) {
  std::vector<std::pair<int, std::string>> segment_and_file;

  if (!fs::exists(directory)) {
    return {};
  }

  const std::string prefix = "nps_hms_coin_" + std::to_string(run_number) + "_";
  const std::string suffix = "_1_-1.root";

  for (const auto& entry : fs::directory_iterator(directory)) {
    if (!entry.is_regular_file()) {
      continue;
    }

    const std::string name = entry.path().filename().string();
    if (!starts_with(name, prefix) || !ends_with(name, suffix)) {
      continue;
    }

    const int seg = parse_segment_from_filename(entry.path().string(), run_number);
    segment_and_file.push_back({seg, entry.path().string()});
  }

  std::sort(segment_and_file.begin(), segment_and_file.end(),
            [](const auto& a, const auto& b) {
              if (a.first == b.first) {
                return a.second < b.second;
              }
              return a.first < b.first;
            });

  std::vector<std::string> files;
  files.reserve(segment_and_file.size());
  for (const auto& item : segment_and_file) {
    files.push_back(item.second);
  }

  return files;
}

inline LocatedRunFiles locate_run_files_prefer_updated(int run_number,
                                                       const std::string& updated_dir,
                                                       const std::string& production_dir) {
  LocatedRunFiles out;

  const auto updated_files = find_run_files_in_directory(updated_dir, run_number);
  if (!updated_files.empty()) {
    out.files = updated_files;
    out.source = "updated";
    out.segment_count = static_cast<int>(updated_files.size());
    return out;
  }

  const auto production_files = find_run_files_in_directory(production_dir, run_number);
  if (!production_files.empty()) {
    out.files = production_files;
    out.source = "production";
    out.segment_count = static_cast<int>(production_files.size());
    return out;
  }

  out.source = "none";
  out.segment_count = 0;
  return out;
}

}  // namespace effstuff
