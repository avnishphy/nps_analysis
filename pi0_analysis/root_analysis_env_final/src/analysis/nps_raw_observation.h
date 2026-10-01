#ifndef NPS_RAW_OBSERVATION_H
#define NPS_RAW_OBSERVATION_H

#include <cstdint>
#include <utility>
#include <vector>

namespace nps_raw_observation {

enum RegionBit : std::uint32_t {
    kPrompt = 1u << 0,
    kDiagonal = 1u << 1,
    kHorizontal = 1u << 2,
    kVertical = 1u << 3,
    kFull1 = 1u << 4,
    kFull2 = 1u << 5
};

enum Category : int {
    kOutside = 0,
    kPromptCategory = 1,
    kDiagonalCategory = 2,
    kHorizontalCategory = 3,
    kVerticalCategory = 4,
    kFull1Category = 5,
    kFull2Category = 6,
    kAmbiguous = 7
};

struct Classification {
    std::uint32_t region_mask = 0;
    int category = kOutside;
};

inline bool in_open_window(double value, const std::pair<double, double>& window)
{
    return value > window.first && value < window.second;
}

inline Classification classify(
    double t1,
    double t2,
    const std::pair<double, double>& prompt_window,
    const std::vector<std::pair<double, double>>& diagonal_windows,
    const std::vector<std::pair<double, double>>& side_windows,
    const std::pair<double, double>& full1_t1,
    const std::pair<double, double>& full1_t2,
    const std::pair<double, double>& full2_t1,
    const std::pair<double, double>& full2_t2)
{
    Classification result;

    if (in_open_window(t1, prompt_window) && in_open_window(t2, prompt_window)) {
        result.region_mask |= kPrompt;
    }
    for (const auto& window : diagonal_windows) {
        if (in_open_window(t1, window) && in_open_window(t2, window)) {
            result.region_mask |= kDiagonal;
            break;
        }
    }
    for (const auto& window : side_windows) {
        if (in_open_window(t1, prompt_window) && in_open_window(t2, window)) {
            result.region_mask |= kHorizontal;
            break;
        }
    }
    for (const auto& window : side_windows) {
        if (in_open_window(t2, prompt_window) && in_open_window(t1, window)) {
            result.region_mask |= kVertical;
            break;
        }
    }
    if (in_open_window(t1, full1_t1) && in_open_window(t2, full1_t2)) {
        result.region_mask |= kFull1;
    }
    if (in_open_window(t1, full2_t1) && in_open_window(t2, full2_t2)) {
        result.region_mask |= kFull2;
    }

    if (result.region_mask == 0) {
        result.category = kOutside;
    } else if ((result.region_mask & (result.region_mask - 1u)) != 0) {
        result.category = kAmbiguous;
    } else if (result.region_mask == kPrompt) {
        result.category = kPromptCategory;
    } else if (result.region_mask == kDiagonal) {
        result.category = kDiagonalCategory;
    } else if (result.region_mask == kHorizontal) {
        result.category = kHorizontalCategory;
    } else if (result.region_mask == kVertical) {
        result.category = kVerticalCategory;
    } else if (result.region_mask == kFull1) {
        result.category = kFull1Category;
    } else if (result.region_mask == kFull2) {
        result.category = kFull2Category;
    }

    return result;
}

inline const char* category_name(int category)
{
    switch (category) {
        case kOutside: return "outside";
        case kPromptCategory: return "prompt";
        case kDiagonalCategory: return "diagonal";
        case kHorizontalCategory: return "horizontal";
        case kVerticalCategory: return "vertical";
        case kFull1Category: return "full1";
        case kFull2Category: return "full2";
        case kAmbiguous: return "ambiguous";
        default: return "invalid";
    }
}

}  // namespace nps_raw_observation

#endif
