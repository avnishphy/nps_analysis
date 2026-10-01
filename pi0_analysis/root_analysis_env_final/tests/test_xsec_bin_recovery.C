// Deterministic response fixtures: exercise production retry/export code without
// event files. Expose analysis state only in this test translation unit.
#include "../src/xsec_extract/xsec_root.h"
#include "../src/xsec_extract/xsec_response.h"
#define private public
#include "../src/xsec_extract/xsec_analysis.h"
#undef private
#define main xsec_production_main
#include "../src/xsec_extract/excl_xsec_pi0_analysis_no_simc_model.C"
#undef main

namespace {
void require(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
AnalysisConfig config() {
    AnalysisConfig c;
    c.n_q2 = 2; c.n_xb = 2; c.n_tprime = 2; c.n_phi = 12;
    c.tgt_contam = 1.; c.write_png = c.write_pdf = false;
    return c;
}
void fixture(ExclPi0XSecAnalysis& a, bool feed_excluded = false) {
    const int ns = a.cfg.n_q2 * a.cfg.n_xb * a.cfg.n_tprime;
    a.slices.resize(ns); a.truth_moments.resize(ns + 6);
    a.phi_edges.resize(a.cfg.n_phi + 1);
    a.q2_edges = {4., 5., 6.}; a.xb_edges_by_q2 = {{.3, .4, .5}, {.3, .4, .5}};
    a.tprime_edges = {-1., -.5, 0.};
    for (int p = 0; p <= a.cfg.n_phi; ++p) a.phi_edges[p] = 2 * TMath::Pi() * p / a.cfg.n_phi;
    a.migration_response.assign(ns * a.cfg.n_phi, std::vector<nps_xsec::ResponseCell>(ns + 6));
    for (int b = 0; b < ns; ++b) {
        auto& m = a.truth_moments[b];
        m.weight = 1.; m.q2 = 5.; m.xb = .4; m.t = -.5; m.tprime = -.3; m.epsilon = .7;
        a.slices[b].phi.resize(a.cfg.n_phi);
        for (int p = 0; p < a.cfg.n_phi; ++p) {
            const double phi = .5 * (a.phi_edges[p] + a.phi_edges[p + 1]);
            auto& row = a.migration_response[b * a.cfg.n_phi + p];
            row[b].basis = {1., std::cos(phi), std::cos(2 * phi)};
            // Groups 2 and 3 mix, so successful recovery must retain their
            // off-diagonal covariance and one simultaneous fit.
            if (b % 4 >= 2) {
                const int other = b % 4 == 2 ? b + 1 : b - 1;
                row[other].basis = {.2, .2 * std::cos(phi), .2 * std::cos(2 * phi)};
            }
            // Group 0 becomes a free nuisance after its reco rows are removed.
            // Distinct angular response makes those nuisance columns identifiable.
            if (feed_excluded && b % 4 == 1)
                row[b - 1].basis = {.2 * (1 + .3 * std::cos(3 * phi)), .03 * std::sin(phi), .03 * std::sin(2 * phi)};
            auto& obs = a.slices[b].phi[p];
            for (int t = 0; t < ns; ++t)
                obs.data += row[t].basis[0] * (10. + t) + row[t].basis[1] * .4 + row[t].basis[2] * -.2;
            obs.data_sumw2 = 1.;
        }
    }
    a.compute_ratios_and_xsec();
}
void unsupported_row(ExclPi0XSecAnalysis& a, int group) {
    // Failure in the second t' bin must remove both t' slices in this group.
    auto& row = a.migration_response[(4 + group) * a.cfg.n_phi + 10];
    for (auto& cell : row) cell = nps_xsec::ResponseCell{};
}
void verify(ExclPi0XSecAnalysis& a, int retained) {
    require(a.successful_fit_groups == retained, "wrong number of retained Q2/xB groups");
    for (size_t b = 0; b < a.slices.size(); ++b) {
        const auto& s = a.slices[b];
        if (a.retained_fit_groups[b % 4]) {
            require(s.fit_xsec.ok, "retained bin has no fit");
            require(std::abs(s.fit_xsec.sigmaU - (10. + b)) < 1e-8, "wrong recovered U");
            require(std::abs(s.fit_xsec.sigmaTL - .4) < 1e-8, "wrong recovered LT");
            require(s.fit_xsec.ndf == a.migration_fit.ndf, "fit was not global");
        } else {
            require(!s.fit_xsec.ok && !std::isfinite(s.fit_xsec.sigmaU), "failed bin published as a cross section");
            require(!std::isfinite(s.phi[0].sim) && !std::isfinite(s.phi[0].xsec), "failed prediction was not cleared");
            require(!s.fit_failure_reason.empty(), "missing failure reason");
        }
    }
}
void export_fixture(ExclPi0XSecAnalysis& a, const fs::path& out) {
    fs::create_directories(out); a.cfg.out_dir = out.string();
    a.cfg.out_csv = (out / "summary.csv").string();
    a.cfg.out_slice_csv = (out / "slices.csv").string();
    a.fout = TFile::Open((out / "results.root").c_str(), "RECREATE");
    a.write_migration_results(); a.write_csv(); a.write_slice_csv();
    if (out.filename() == "partial" || out.filename() == "all_failed") {
        a.cfg.write_pdf = true;
        a.cfg.out_all_plots_pdf = (out / "plots.pdf").string();
        a.make_slice_plots();
        a.make_sigma_vs_tprime_plots();
    }
}
}

int main(int argc, char** argv) {
    {
        AnalysisConfig c;
        c.phi_bin_edges[1] = .13;
        c.tprime_bin_edges[1] = -.62;
        c.q2_bin_edges = {3., 3.7, 5.}; c.n_q2 = 2;
        c.xb_bin_edges_by_q2 = {{.25, .31, .455}, {.25, .39, .455}};
        c.n_xb = 2;
        validate_xsec_binning(c);
        ExclPi0XSecAnalysis a(c);
        a.build_binning();
        require(a.phi_edges == c.phi_bin_edges &&
                a.tprime_edges == c.tprime_bin_edges &&
                a.q2_edges == c.q2_bin_edges &&
                a.xb_edges_by_q2 == c.xb_bin_edges_by_q2,
                "response extractor did not use configured edges");
        c.phi_bin_edges[1] = c.phi_bin_edges[0];
        bool rejected = false;
        try { validate_xsec_binning(c); }
        catch (const std::runtime_error&) { rejected = true; }
        require(rejected, "unordered configured edges were accepted");
    }
    gROOT->SetBatch(kTRUE);
    const fs::path out = argc > 1 ? argv[1] : "/tmp/nps_xsec_recovery_test";
    {
        ExclPi0XSecAnalysis a(config()); fixture(a); a.fit_slices(); verify(a, 4);
        require(!a.fit_fallback && a.fit_attempts.size() == 1, "healthy global fit changed path");
        export_fixture(a, out / "healthy");
    }
    {
        ExclPi0XSecAnalysis a(config()); fixture(a); unsupported_row(a, 0);
        a.fit_slices(); verify(a, 3);
        require(a.fit_attempts.size() == 2, "unsupported group did not recover in one retry");
        auto col = [&](int block) {
            return 3 * (std::find(a.active_truth_blocks.begin(), a.active_truth_blocks.end(), block) - a.active_truth_blocks.begin());
        };
        const size_t np = a.migration_fit.parameters.size();
        require(std::abs(a.migration_fit.covariance[col(2) * np + col(3)]) > 1e-5, "cross-bin covariance lost");
        export_fixture(a, out / "partial");
    }
    {
        ExclPi0XSecAnalysis a(config()); fixture(a, true); unsupported_row(a, 0);
        a.fit_slices(); verify(a, 3);
        require(std::find(a.active_truth_blocks.begin(), a.active_truth_blocks.end(), 0) != a.active_truth_blocks.end(),
                "excluded-bin truth feed-in was dropped");
        export_fixture(a, out / "nuisance");
    }
    {
        ExclPi0XSecAnalysis a(config()); fixture(a);
        // Nonzero but rank-deficient response in group 1: two identical columns.
        for (int b : {1, 5}) for (int p = 0; p < a.cfg.n_phi; ++p)
            a.migration_response[b * a.cfg.n_phi + p][b].basis[2] = a.migration_response[b * a.cfg.n_phi + p][b].basis[1];
        unsupported_row(a, 0); a.fit_slices(); verify(a, 2);
        require(a.retained_fit_groups[2] && a.retained_fit_groups[3], "wrong rank-failure subset");
        export_fixture(a, out / "rank");
    }
    {
        ExclPi0XSecAnalysis a(config()); fixture(a);
        a.truth_moments[0].weight = 0;
        for (auto& row : a.migration_response) row[0] = nps_xsec::ResponseCell{};
        a.compute_ratios_and_xsec(); // unit target factor; empty truth must not abort here
        a.fit_slices(); verify(a, 3);
    }
    {
        auto c = config(); c.tgt_contam = .5;
        ExclPi0XSecAnalysis a(c); fixture(a); unsupported_row(a, 0);
        const double before = a.slices[1].phi[0].data;
        a.fit_slices();
        require(a.slices[1].phi[0].data == before, "target correction applied again on retry");
        require(std::abs(a.slices[1].fit_xsec.sigmaU - 22.) < 1e-8, "target correction not preserved");
    }
    {
        ExclPi0XSecAnalysis a(config()); fixture(a);
        for (int g = 0; g < 4; ++g) unsupported_row(a, g);
        a.fit_slices(); verify(a, 0); export_fixture(a, out / "all_failed");
    }
    {
        auto c = config(); c.mc_max_iterations = 0; // force exhaustion in every attempted subset
        ExclPi0XSecAnalysis a(c); fixture(a); a.fit_slices(); verify(a, 0);
        require(a.fit_attempts.size() == 15, "not all nonempty subsets were attempted");
    }
    std::cout << "PASS global recovery, covariance, nuisance feed-in, rank failure, empty truth, all-failed exports, convergence exhaustion\n";
}
