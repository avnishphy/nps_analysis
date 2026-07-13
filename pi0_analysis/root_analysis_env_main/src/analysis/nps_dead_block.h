#ifndef NPS_DEAD_BLOCK_H
#define NPS_DEAD_BLOCK_H

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace nps_dead_block {

struct Mask {
    static constexpr int kNCols = 30;
    static constexpr int kNRows = 36;
    static constexpr int kNBlocks = kNCols * kNRows;
    static constexpr double kXMin = -30.0;
    static constexpr double kXMax = 30.0;
    static constexpr double kYMin = -36.0;
    static constexpr double kYMax = 36.0;

    std::set<int> dead_blocks;
    std::set<int> excluded_blocks;

    bool empty() const noexcept { return excluded_blocks.empty(); }

    static bool block_to_col_row(int block, int& col, int& row) noexcept {
        if (block < 0 || block >= kNBlocks) return false;
        col = block % kNCols;
        row = block / kNCols;
        return true;
    }

    static int col_row_to_block(int col, int row) noexcept {
        if (col < 0 || col >= kNCols || row < 0 || row >= kNRows) return -1;
        return row * kNCols + col;
    }

    static int xy_to_block(double x, double y) noexcept {
        if (!std::isfinite(x) || !std::isfinite(y)) return -1;
        if (x < kXMin || x >= kXMax || y < kYMin || y >= kYMax) return -1;

        const double cell_w = (kXMax - kXMin) / static_cast<double>(kNCols);
        const double cell_h = (kYMax - kYMin) / static_cast<double>(kNRows);
        int col = static_cast<int>(std::floor((x - kXMin) / cell_w));
        int row = static_cast<int>(std::floor((y - kYMin) / cell_h));

        if (col < 0) col = 0;
        if (row < 0) row = 0;
        if (col >= kNCols) col = kNCols - 1;
        if (row >= kNRows) row = kNRows - 1;

        return col_row_to_block(col, row);
    }

    bool rejects_xy(double x, double y) const noexcept {
        const int block = xy_to_block(x, y);
        return block >= 0 && excluded_blocks.count(block) != 0;
    }
};

struct RunRange {
    int start_run = 0;
    int end_run = 0;
    std::set<int> dead_blocks;
};

inline std::string trim(std::string s) {
    auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [&](unsigned char c) { return !is_space(c); }));
    s.erase(std::find_if(s.rbegin(), s.rend(), [&](unsigned char c) { return !is_space(c); }).base(), s.end());
    return s;
}

inline std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return s;
}

inline std::vector<std::string> split_csv_quoted(const std::string& line) {
    std::vector<std::string> fields;
    std::string cur;
    bool in_quotes = false;

    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (c == '"') {
            if (in_quotes && (i + 1 < line.size()) && line[i + 1] == '"') {
                cur.push_back('"');
                ++i;
            } else {
                in_quotes = !in_quotes;
            }
        } else if (c == ',' && !in_quotes) {
            fields.push_back(trim(cur));
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    fields.push_back(trim(cur));
    return fields;
}

inline bool parse_int_cell(const std::vector<std::string>& row, int idx, int& out) {
    if (idx < 0 || idx >= static_cast<int>(row.size())) return false;
    const std::string cell = trim(row[idx]);
    if (cell.empty()) return false;

    try {
        size_t pos = 0;
        const int value = std::stoi(cell, &pos, 10);
        if (pos != cell.size()) return false;
        out = value;
        return true;
    } catch (...) {
        return false;
    }
}

inline void add_dead_block_and_neighbors(Mask& mask, int block) {
    int col = -1;
    int row = -1;
    if (!Mask::block_to_col_row(block, col, row)) return;

    mask.dead_blocks.insert(block);
    for (int dr = -1; dr <= 1; ++dr) {
        for (int dc = -1; dc <= 1; ++dc) {
            const int neighbor = Mask::col_row_to_block(col + dc, row + dr);
            if (neighbor >= 0) mask.excluded_blocks.insert(neighbor);
        }
    }
}

inline bool load_run_ranges(const std::string& csv_path,
                            std::vector<RunRange>& ranges,
                            std::string& error) {
    ranges.clear();

    std::ifstream in(csv_path);
    if (!in.is_open()) {
        error = "cannot open dead-block CSV: " + csv_path;
        return false;
    }

    std::string line;
    if (!std::getline(in, line)) {
        error = "dead-block CSV is empty: " + csv_path;
        return false;
    }

    if (line.size() >= 3 &&
        static_cast<unsigned char>(line[0]) == 0xEF &&
        static_cast<unsigned char>(line[1]) == 0xBB &&
        static_cast<unsigned char>(line[2]) == 0xBF) {
        line.erase(0, 3);
    }

    const std::vector<std::string> header = split_csv_quoted(line);
    int idx_start = -1;
    int idx_end = -1;
    int idx_blocks = -1;

    for (size_t i = 0; i < header.size(); ++i) {
        const std::string h = lower(trim(header[i]));
        if (h == "start run") idx_start = static_cast<int>(i);
        else if (h == "end run") idx_end = static_cast<int>(i);
        else if (h == "blocks (dead/missing)") idx_blocks = static_cast<int>(i);
    }

    if (idx_start < 0 || idx_end < 0 || idx_blocks < 0) {
        error = "dead-block CSV is missing Start run, End run, or Blocks (Dead/Missing) header";
        return false;
    }

    while (std::getline(in, line)) {
        if (trim(line).empty()) continue;
        const std::vector<std::string> row = split_csv_quoted(line);

        int start_run = 0;
        int end_run = 0;
        if (!parse_int_cell(row, idx_start, start_run) || !parse_int_cell(row, idx_end, end_run)) {
            continue;
        }
        if (end_run < start_run) std::swap(start_run, end_run);

        RunRange range;
        range.start_run = start_run;
        range.end_run = end_run;

        for (int i = idx_blocks; i < static_cast<int>(row.size()); ++i) {
            int block = -1;
            if (parse_int_cell(row, i, block) && block >= 0 && block < Mask::kNBlocks) {
                range.dead_blocks.insert(block);
            }
        }

        if (!range.dead_blocks.empty()) ranges.push_back(range);
    }

    return true;
}

inline Mask mask_for_run(const std::vector<RunRange>& ranges, int run) {
    Mask mask;
    for (const auto& range : ranges) {
        if (run < range.start_run || run > range.end_run) continue;
        for (int block : range.dead_blocks) {
            add_dead_block_and_neighbors(mask, block);
        }
    }
    return mask;
}

inline std::string join_blocks(const std::set<int>& blocks, int max_items = 18) {
    std::ostringstream out;
    int count = 0;
    for (int block : blocks) {
        if (count > 0) out << ", ";
        if (count >= max_items) {
            out << "...";
            break;
        }
        out << block;
        ++count;
    }
    return count > 0 ? out.str() : std::string("none");
}

}  // namespace nps_dead_block

#endif  // NPS_DEAD_BLOCK_H
