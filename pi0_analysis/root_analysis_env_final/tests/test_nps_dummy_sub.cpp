// Standalone regression test; also a small read-only ROOT validation driver.
#include "../src/analysis/nps_dummy_sub.h"
#include <TFile.h>
#include <TTree.h>
#include <algorithm>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>

namespace d = nps::dummy;

void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
void near(double actual, double expected) {
    check(std::abs(actual - expected) <= 1e-11 * std::max(1.0, std::abs(expected)),
          "numeric mismatch");
}
template<class F> void rejects(F action) {
    bool caught = false;
    try { action(); } catch (const std::invalid_argument&) { caught = true; }
    check(caught, "invalid input accepted");
}

void unit_tests() {
    check(d::window(-1) == d::Window::Upstream, "upstream sign");
    check(d::window(1) == d::Window::Downstream, "downstream sign");
    check(d::window(0) == d::Window::Unassigned, "zero split");
    near(d::material_weight(-1), 1.0 / 8.467);
    near(d::material_weight(1), 1.0 / 4.256);
    near(d::dummy_subtraction_weight(1, 0, 2), 0);
    near(d::effective_charge_mC(10000, 2, .8, .9, .95), 3.42);
    // Estimated livetime can slightly exceed one in the real CSV.
    near(d::effective_charge_mC(1000, 1, 1.002, 1, 1), 1.002);

    d::Config cfg;
    cfg.upstream_ratio = 2; cfg.downstream_ratio = 4;
    d::WeightedSum lh2, up, down;
    lh2.add(2); lh2.add(3); up.add(4); down.add(2);
    const auto result = d::subtract(lh2, up, down, 10, 2, cfg);
    near(result.subtracted.sumw, -.75); // Negative bins must survive.
    near(result.subtracted.sumw2, 1.1925);
    near(result.background.sumw, 1.25);
    near(result.subtracted.stat_error(), std::sqrt(1.1925));
    d::WeightedSum signed_events;
    signed_events.add(d::production_event_weight(2, 10));
    signed_events.add(d::production_event_weight(3, 10));
    signed_events.add(d::dummy_subtraction_weight(4, -1, 2, cfg));
    signed_events.add(d::dummy_subtraction_weight(2, 1, 2, cfg));
    near(signed_events.sumw, result.subtracted.sumw);
    near(signed_events.sumw2, result.subtracted.sumw2);
    near(d::dummy_subtraction_weight(-2, 1, 2, cfg), .25);
    // Splitting one exposure into two runs must not double its yield.
    near(d::production_event_weight(4, 2),
         d::production_event_weight(2, 1 + 1) * 2);
    cfg.ytar_split_cm = 1;
    check(d::window(.5, cfg) == d::Window::Upstream, "custom split");
    const double nan = std::numeric_limits<double>::quiet_NaN();
    rejects([&] { d::window(nan); });
    rejects([&] { d::production_event_weight(nan, 1); });
    rejects([&] { d::subtract(lh2, up, down, 1, 0); });
    rejects([&] { d::effective_charge_mC(1000, 0, 1, 1, 1); });
    cfg.upstream_ratio = -1;
    rejects([&] { d::material_weight(-1, cfg); });
    up.sumw2 = -1;
    rejects([&] { d::subtract(lh2, up, down, 1, 1); });
    std::cout << "PASS unit tests\n";
}

int main(int argc, char** argv) {
    try {
        unit_tests();
        if (argc == 1) return 0;
        check(argc == 6, "usage: test BASE Q_LH2 Q_DUMMY_4260 Q_DUMMY_4261 LH2_RUN");
        const double q_lh2 = std::stod(argv[2]);
        const double q_dummy = std::stod(argv[3]) + std::stod(argv[4]);
        // argv[5] is the configurable LH2 run, enabling a different small fixture.
        const int runs[] = {std::stoi(argv[5]), 4260, 4261};
        d::WeightedSum lh2, up, down, signed_events;
        for (int i = 0; i < 3; ++i) {
            const std::string path = std::string(argv[1]) + "/output/KinC_x60_4b/root/diagnostics_run"
                                     + std::to_string(runs[i]) + ".root";
            std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
            check(file && !file->IsZombie(), "cannot open fixture");
            auto* tree = dynamic_cast<TTree*>(file->Get("physics"));
            check(tree != nullptr, "missing physics tree");
            double ytar = 0, weight = 0;
            tree->SetBranchStatus("*", 0);
            tree->SetBranchStatus("ytar", 1);
            tree->SetBranchStatus("pi0_weight", 1);
            check(tree->SetBranchAddress("ytar", &ytar) >= 0, "missing ytar");
            check(tree->SetBranchAddress("pi0_weight", &weight) >= 0, "missing pi0_weight");
            for (Long64_t j = 0; j < tree->GetEntries(); ++j) {
                check(tree->GetEntry(j) > 0, "entry read failed");
                if (i == 0) {
                    lh2.add(weight);
                    signed_events.add(d::production_event_weight(weight, q_lh2));
                } else {
                    const auto side = d::window(ytar);
                    if (side == d::Window::Upstream) up.add(weight);
                    if (side == d::Window::Downstream) down.add(weight);
                    signed_events.add(d::dummy_subtraction_weight(weight, ytar, q_dummy));
                }
            }
        }
        const auto result = d::subtract(lh2, up, down, q_lh2, q_dummy);
        near(signed_events.sumw, result.subtracted.sumw);
        near(signed_events.sumw2, result.subtracted.sumw2);
        std::cout << std::setprecision(17) << "RESULT " << result.lh2.sumw << ' '
                  << result.background.sumw << ' ' << result.subtracted.sumw << ' '
                  << result.subtracted.sumw2 << '\n';
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
