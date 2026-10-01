// Synthetic selected-event coverage test; no production files or PARTONS.
// Compile with ROOT flags after sourcing the Hall C/NPS environment.
#include "../src/xsec_extract/xsec_root.h"
#include "../src/xsec_extract/xsec_response.h"
#define private public
#include "../src/xsec_extract/xsec_analysis.h"
#undef private
#define main xsec_coverage_production_main
#include "../src/xsec_extract/excl_xsec_pi0_analysis_no_simc_model.C"
#undef main

namespace {
void require_coverage(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

std::vector<double> histogram_snapshot(const TH2D& h) {
    std::vector<double> result{h.GetEntries(), h.GetXaxis()->GetXmin(), h.GetXaxis()->GetXmax(),
        h.GetYaxis()->GetXmin(), h.GetYaxis()->GetXmax(), h.GetMinimumStored(), h.GetMaximumStored()};
    for (int bin = 0; bin < h.GetNcells(); ++bin) {
        result.push_back(h.GetBinContent(bin));
        result.push_back(h.GetBinError(bin));
    }
    return result;
}

void check_counts(const TH2D& h, long long expected) {
    require_coverage(h.GetEntries() == expected, "coverage entry count differs from response");
    require_coverage(h.Integral() == expected, "regular coverage bins lost selected events");
    require_coverage(h.Integral(0, h.GetNbinsX() + 1, 0, h.GetNbinsY() + 1) == expected,
                     "coverage has nonzero underflow/overflow");
}

void fixture(const fs::path& directory, bool graphics, bool diagnostics) {
    AnalysisConfig c;
    c.n_tprime = c.n_q2 = c.n_xb = 1; c.n_phi = 4;
    c.q2_min = 5.; c.q2_max = 6.; c.xb_min = .5; c.xb_max = .6;
    c.tprime_min = -1.; c.tprime_max = 0.;
    c.out_dir = directory.string(); c.out_all_plots_pdf = "combined.pdf";
    c.write_pdf = c.write_png = graphics; c.diagnostics = diagnostics;
    fs::create_directories(directory);
    ExclPi0XSecAnalysis a(c);
    a.q2_edges = {c.q2_min, c.q2_max};
    a.xb_edges = {c.xb_min, c.xb_max}; a.xb_edges_by_q2 = {a.xb_edges};
    a.tprime_edges = {c.tprime_min, c.tprime_max};
    a.phi_edges = {0., .5 * TMath::Pi(), TMath::Pi(), 1.5 * TMath::Pi(), 2. * TMath::Pi()};
    a.has_vertex_kinematics = true;
    a.init_storage();
    const auto fill = [&](double generated_q2, double generated_xb, int exclusive = 1,
                          float mass = .94f, float reco_q2 = 5.5f) {
        const double W = std::sqrt(c.mp * c.mp + generated_q2 * (1. / generated_xb - 1.));
        const double positive_minus_t = -nps_xsec::forward_t(generated_q2, W, c.mp, c.mpi0) + .2;
        a.fill_sim_event(reco_q2, -.9f, -.6f, .55f, .4f, 2.e-8f, exclusive,
                         mass, 1.e-8f, generated_q2, W, positive_minus_t, .3f, 0, false, 2.3f);
    };
    fill(5.5, .55); // Published truth cell.
    // Both generated axes extend beyond the initial half-span margins. Their
    // reconstructed coordinates remain inside selection: these are guard events.
    fill(2., .25);
    fill(9., .85);
    fill(5.5, .55, 0);       // Wrong generated channel.
    fill(5.5, .55, 1, 2.f);  // Reconstructed missing mass rejected.
    fill(5.5, .55, 1, .94f, 7.f); // Reconstructed Q2 rejected.

    long long response_events = 0;
    for (const auto& row : a.migration_response)
        for (const auto& cell : row) response_events += cell.events;
    require_coverage(response_events == 3, "fixture response selection changed");
    require_coverage(a.truth_moments[0].events == 1 && a.truth_moments[3].events == 1 &&
                     a.truth_moments[4].events == 1, "outside-range truth did not reach Q2 guards");
    const auto original_response = a.migration_response;
    a.migration_fit.parameters = {12., -3., 1.};
    const auto original_parameters = a.migration_fit.parameters;
    a.fout = TFile::Open((directory / "coverage.root").c_str(), "RECREATE");
    require_coverage(a.fout && !a.fout->IsZombie(), "cannot create coverage fixture output");
    gStyle->SetPalette(55); gStyle->SetNumberContours(37);
    std::vector<int> palette;
    for (int i = 0; i < gStyle->GetNumberOfColors(); ++i) palette.push_back(gStyle->GetColorPalette(i));

    if (!diagnostics) {
        require_coverage(!a.h_migration_vertex_q2_xb && !a.h_migration_reco_q2_xb,
                         "disabled diagnostics allocated coverage histograms");
        a.make_migration_coverage_plots();
        require_coverage(!a.fout->GetDirectory("migration"), "disabled coverage wrote ROOT objects");
        return;
    }
    check_counts(*a.h_migration_vertex_q2_xb, response_events);
    check_counts(*a.h_migration_reco_q2_xb, response_events);
    auto* generated = a.h_migration_vertex_q2_xb.get();
    require_coverage(generated->GetXaxis()->GetXmin() <= .25 && generated->GetXaxis()->GetXmax() > .85 &&
                     generated->GetYaxis()->GetXmin() <= 2. && generated->GetYaxis()->GetXmax() > 9.,
                     "generated coverage axes did not extend for feed-in");
    const auto before_generated = histogram_snapshot(*generated);
    const auto before_reco = histogram_snapshot(*a.h_migration_reco_q2_xb);
    a.make_migration_coverage_plots();
    require_coverage(histogram_snapshot(*generated) == before_generated &&
                     histogram_snapshot(*a.h_migration_reco_q2_xb) == before_reco,
                     "display density conversion mutated stored count histograms");
    require_coverage(a.migration_fit.parameters == original_parameters, "coverage changed fit parameters");
    for (size_t r = 0; r < original_response.size(); ++r)
        for (size_t b = 0; b < original_response[r].size(); ++b) {
            const auto& before = original_response[r][b]; const auto& after = a.migration_response[r][b];
            require_coverage(before.events == after.events && before.basis == after.basis &&
                             before.covariance == after.covariance, "coverage changed response cells");
        }
    require_coverage(gStyle->GetNumberContours() == 37 &&
                     gStyle->GetNumberOfColors() == static_cast<int>(palette.size()), "coverage changed global style");
    for (size_t i = 0; i < palette.size(); ++i)
        require_coverage(gStyle->GetColorPalette(i) == palette[i], "coverage changed global palette colors");
    for (const char* name : {"vertex_q2_xb_selected", "reco_q2_xb_selected"}) {
        auto* persisted = dynamic_cast<TH2D*>(a.fout->Get((std::string("migration/") + name).c_str()));
        require_coverage(persisted != nullptr, "missing persisted coverage histogram");
        check_counts(*persisted, response_events);
    }
    require_coverage(a.fout->Get("migration/coverage_semantics") != nullptr, "missing coverage semantics");
    a.close_combined_pdf();
    for (const char* extension : {".pdf", ".png"})
        require_coverage(fs::exists(directory / "migration" / (std::string("vertex_reco_q2_xb_coverage") + extension)) == graphics,
                         "coverage graphics switch/output contract broken");
}
} // namespace

int main(int argc, char** argv) {
    gROOT->SetBatch(kTRUE);
    TH1::SetDefaultSumw2(kTRUE);
    const fs::path out = argc > 1 ? argv[1] : "/tmp/nps_migration_coverage_check";
    fixture(out / "graphics", true, true);
    fixture(out / "root_only", false, true);
    fixture(out / "disabled", true, false);
    std::cout << "PASS coverage selected-event accounting, guard feed-in, axis extension, count persistence, "
                 "display isolation, palette restoration, and graphics switches\n";
}
