#include <cmath>
#include <iostream>
#include <vector>

#include "../src/analysis/nps_2d_mass_cut.h"

int main() {
    std::vector<nps2d::Point> points;
    Long64_t id = 0;

    // Anchored anti-correlated signal component.
    for (int ix = -14; ix <= 14; ++ix) {
        const double x = 0.13498 + 0.00035 * ix;
        for (int iy = -14; iy <= 14; ++iy) {
            const double residual_y = 0.004 * iy;
            const double y = 0.938 - 15.0 * (x - 0.13498) + residual_y;
            const double radius2 = std::pow(ix / 7.0, 2) + std::pow(iy / 7.0, 2);
            const double weight = 4.0 * std::exp(-0.5 * radius2);
            points.push_back({id++, x, y, weight});
        }
    }

    // Stronger remote ridge: old global-maximum logic can lock onto this.
    for (int ix = -20; ix <= 20; ++ix) {
        const double x = 0.122 + 0.00045 * ix;
        for (int iy = -18; iy <= 18; ++iy) {
            const double y = 1.30 - 12.0 * (x - 0.122) + 0.006 * iy;
            const double weight = 5.0 * std::exp(-0.5 * (std::pow(ix / 10.0, 2) +
                                                         std::pow(iy / 10.0, 2)));
            points.push_back({id++, x, y, weight});
        }
    }

    nps2d::Config cfg;
    cfg.write_debug = false;
    const nps2d::Result result = nps2d::evaluate_mass_cuts(points, cfg);
    if (!result.params.valid || !result.params.ellipse_valid) {
        std::cerr << "mass-cut fit invalid\n";
        return 1;
    }
    if (std::abs(result.params.mean_mpi0 - cfg.seed_mpi0) > cfg.max_model_mpi0_offset ||
        std::abs(result.params.mean_mmiss - cfg.seed_mmiss) > cfg.max_model_mmiss_offset) {
        std::cerr << "ellipse model escaped physical anchor\n";
        return 2;
    }
    if (result.params.fit_subset_total_fraction > cfg.auto_max_core_total_fraction) {
        std::cerr << "fit subset exceeded leak limit\n";
        return 3;
    }
    if (result.params.fit_subset_total_fraction <= 0.10 || result.params.ellipse_growth_steps <= 0) {
        std::cerr << "fit subset did not grow from density seed\n";
        return 5;
    }
    if (result.params.core_total_fraction <= 0.0 ||
        result.params.core_total_fraction >= result.params.ellipse_total_fraction) {
        std::cerr << "reported core is not nested inside ellipse\n";
        return 6;
    }
    if (result.params.peak_mmiss < cfg.seed_mmiss - cfg.seed_mmiss_half_width ||
        result.params.peak_mmiss > cfg.seed_mmiss + cfg.seed_mmiss_half_width) {
        std::cerr << "peak escaped seed window\n";
        return 4;
    }

    std::cout << "PASS peak=(" << result.params.peak_mpi0 << ","
              << result.params.peak_mmiss << ") core="
              << result.params.core_total_fraction << " fit_subset="
              << result.params.fit_subset_total_fraction << " mean=("
              << result.params.mean_mpi0 << "," << result.params.mean_mmiss
              << ") mcd_valid=" << result.params.mcd_valid << "\n";
    return 0;
}
