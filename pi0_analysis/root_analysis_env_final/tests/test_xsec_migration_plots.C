// Synthetic normalization and output-contract test; no production inputs.
// Compile with root-config --cflags --libs in the Hall C/NPS environment.
// Usage: test_xsec_migration_plots [scratch_output_directory]
#include "../src/xsec_extract/xsec_root.h"
#include "../src/xsec_extract/xsec_response.h"
#define private public
#include "../src/xsec_extract/xsec_analysis.h"
#undef private
#include "../src/xsec_extract/xsec_plot_global.h"
#include "../src/xsec_extract/xsec_plot_migration.h"
#include "../src/xsec_extract/xsec_output.h"

namespace {
void require_migration(bool ok, const char* explanation) {
    if (!ok) throw std::runtime_error(explanation);
}
bool close_migration(double a, double b) { return std::abs(a - b) < 1e-14; }

void fixture(const fs::path& directory, bool graphics, bool diagnostics) {
    AnalysisConfig c;
    c.n_tprime = 3; c.n_q2 = c.n_xb = 1; c.n_phi = 2;
    c.out_dir = directory.string();
    c.write_pdf = c.write_png = graphics; c.diagnostics = diagnostics;
    fs::create_directories(directory);
    ExclPi0XSecAnalysis a(c);
    a.slices.resize(3); a.truth_moments.resize(9);
    a.migration_response.assign(6, std::vector<nps_xsec::ResponseCell>(9));
    const auto fill = [&](int r, int b, int n) {
        auto& cell = a.migration_response[r][b];
        cell.events = n;
        cell.basis = {double(n), double(n) * (r % 2 ? -.8 : .7), double(n) * (r > 1 ? .3 : -.4)};
    };
    // Two supported published truth blocks, one empty published block, one
    // populated guard, five empty guards, and one unsupported reconstructed
    // slice. Unequal phi counts exercise aggregation BEFORE normalization.
    fill(0, 0, 15); fill(1, 0, 5); fill(2, 0, 4); fill(3, 0, 1);
    fill(0, 1, 1); fill(1, 1, 4); fill(2, 1, 5); fill(3, 1, 15);
    fill(0, 3, 2); fill(3, 3, 3);
    const auto original_response = a.migration_response;
    a.migration_fit.parameters = {99., -72., 31.};
    const auto original_fit = a.migration_fit.parameters;
    a.fout = TFile::Open((directory / "fixture.root").c_str(), "RECREATE");
    gStyle->SetPalette(55); gStyle->SetNumberContours(37);
    std::vector<int> original_palette;
    for (int i = 0; i < gStyle->GetNumberOfColors(); ++i)
        original_palette.push_back(gStyle->GetColorPalette(i));
    a.make_migration_plots();
    require_migration(gStyle->GetNumberContours() == 37, "global contour count changed");
    require_migration(gStyle->GetNumberOfColors() == static_cast<int>(original_palette.size()), "palette length changed");
    for (size_t i = 0; i < original_palette.size(); ++i)
        require_migration(gStyle->GetColorPalette(i) == original_palette[i], "palette color changed");
    require_migration(a.migration_fit.parameters == original_fit, "diagnostics changed fitted coefficients");
    for (size_t r = 0; r < original_response.size(); ++r)
        for (size_t b = 0; b < original_response[r].size(); ++b) {
            const auto& before = original_response[r][b];
            const auto& after = a.migration_response[r][b];
            require_migration(after.events == before.events && after.basis == before.basis &&
                              after.covariance == before.covariance, "diagnostics changed response cells");
        }
    if (!diagnostics) {
        require_migration(a.fout->GetDirectory("migration") == nullptr, "disabled diagnostics wrote migration objects");
        require_migration(!fs::exists(directory / "migration"), "disabled diagnostics wrote plots");
        return;
    }
    const auto get = [&](const char* name) {
        auto* h = dynamic_cast<TH2D*>(a.fout->Get((std::string("migration/") + name).c_str()));
        require_migration(h != nullptr, "missing migration ROOT histogram");
        return h;
    };
    auto* n = get("selected_counts");
    auto* pr = get("p_reco_given_truth_selected");
    auto* pt = get("p_truth_given_reco_selected");
    require_migration(n->GetNbinsX() == 9 && n->GetNbinsY() == 3, "full truth-block indices not preserved");
    require_migration(n->Integral() == 55., "selected event count not conserved");
    require_migration(n->GetBinContent(1, 1) == 20. && n->GetBinContent(1, 2) == 5., "wrong phi count aggregation");
    require_migration(close_migration(pr->GetBinContent(1, 1), .8), "wrong P(reco|truth) direction");
    require_migration(close_migration(pt->GetBinContent(1, 1), 20. / 27.), "guard excluded from P(truth|reco) denominator");
    for (int b = 1; b <= 9; ++b) {
        double total = 0.;
        for (int s = 1; s <= 3; ++s) total += pr->GetBinContent(b, s);
        require_migration(close_migration(total, b == 1 || b == 2 || b == 4 ? 1. : 0.), "unsupported truth column or column sum incorrect");
    }
    for (int s = 1; s <= 3; ++s) {
        double total = 0.;
        for (int b = 1; b <= 9; ++b) total += pt->GetBinContent(b, s);
        require_migration(close_migration(total, s < 3 ? 1. : 0.), "unsupported reco row or row sum incorrect");
    }
    const char* components[] = {"response_U", "response_LT", "response_TT"};
    for (int component = 0; component < 3; ++component) {
        auto* h = get(components[component]);
        require_migration(h->GetNbinsX() == 9 && h->GetNbinsY() == 6, "response dimensions changed");
        for (int r = 0; r < 6; ++r) for (int b = 0; b < 9; ++b)
            require_migration(h->GetBinContent(b + 1, r + 1) == original_response[r][b].basis[component],
                              "raw signed response was normalized or transposed");
    }
    require_migration(a.fout->Get("migration/definitions") != nullptr, "normalization metadata missing");
    for (const char* stem : {"response_components", "selected_migration_support", "selected_migration_fractions"})
        for (const char* extension : {".pdf", ".png"})
            require_migration(fs::exists(directory / "migration" / (std::string(stem) + extension)) == graphics,
                              "graphics switch/output contract broken");
}
}

int main(int argc, char** argv) {
    gROOT->SetBatch(true);
    const fs::path out = argc > 1 ? argv[1] : "/tmp/nps_xsec_migration_plots_test";
    fixture(out / "no_graphics", false, true);
    fixture(out / "graphics", true, true);
    fixture(out / "disabled", false, false);
    std::cout << "PASS migration counts/fractions, raw signed response, guards/zero support, unchanged fit, ROOT/graphics switches, palette restoration\n";
}
