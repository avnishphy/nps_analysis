#include "../src/analysis/nps_raw_observation.h"

#include <cassert>
#include <utility>
#include <vector>

int main()
{
    using namespace nps_raw_observation;
    const std::pair<double, double> prompt{149.0, 151.0};
    const std::vector<std::pair<double, double>> side{
        {141.0, 143.0}, {143.0, 145.0}, {145.0, 147.0},
        {153.0, 155.0}, {155.0, 157.0}, {157.0, 159.0}};
    const std::pair<double, double> full1_t1{153.0, 159.0};
    const std::pair<double, double> full1_t2{141.0, 147.0};
    const std::pair<double, double> full2_t1{141.0, 147.0};
    const std::pair<double, double> full2_t2{153.0, 159.0};

    auto category = [&](double t1, double t2) {
        return classify(t1, t2, prompt, side, side,
                        full1_t1, full1_t2, full2_t1, full2_t2);
    };

    assert(category(150.0, 150.0).category == kPromptCategory);
    assert(category(142.0, 142.0).category == kDiagonalCategory);
    assert(category(150.0, 142.0).category == kHorizontalCategory);
    assert(category(142.0, 150.0).category == kVerticalCategory);
    assert(category(156.0, 144.0).category == kFull1Category);
    assert(category(144.0, 156.0).category == kFull2Category);
    assert(category(148.0, 152.0).category == kOutside);

    // Existing histogram fills use strict open boundaries.
    assert(category(149.0, 150.0).category == kOutside);
    assert(category(150.0, 141.0).category == kOutside);

    for (double t1 : {142.0, 144.0, 146.0, 150.0, 154.0, 156.0, 158.0}) {
        for (double t2 : {142.0, 144.0, 146.0, 150.0, 154.0, 156.0, 158.0}) {
            const auto result = category(t1, t2);
            assert(result.category >= kOutside && result.category <= kAmbiguous);
            assert(result.category != kAmbiguous);
            assert((result.region_mask & ~std::uint32_t{63}) == 0);
        }
    }

    return 0;
}
