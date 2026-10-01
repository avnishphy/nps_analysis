#include <cassert>
#include <string>
#include "TROOT.h"
#include "../src/analysis/nps_plot_diagnostics.h"

int main(int argc, char** argv) {
    gROOT->SetBatch(true);
    npsplot::Diagnostics diagnostics;
    TH1D before("h_test_run1", "test", 20, 0., 2.);
    before.SetDirectory(nullptr);
    diagnostics.register_after(&before);
    for (double x : {-0.1, 0.5, 1.0, 2.5}) {
        before.Fill(x);
        diagnostics.stage_pass = x >= 0.4 && x <= 1.1;
        diagnostics.fill(&before, x, 1.0);
    }
    assert(before.GetEntries() == 4);
    assert(before.GetBinContent(0) == 1 && before.GetBinContent(21) == 1);
    assert(diagnostics.after.at(&before)->GetEntries() == 2);
    const auto bins = diagnostics.display_bins(&before, 0.2, 1.8);
    assert(before.GetXaxis()->GetBinLowEdge(bins.first) <= 0.2);
    assert(before.GetXaxis()->GetBinUpEdge(bins.second) >= 1.8);
    // Sparse cluster occupancy must not zoom away the detector edges.
    TH1D cluster_x("h_cut_cluster_x_run1", "x", 34, -34., 34.);
    TH1D cluster_y("h_cut_cluster_y_run1", "y", 40, -40., 40.);
    for (TH1D* h : {&cluster_x, &cluster_y}) {
        h->Fill(1.);
        const auto full = diagnostics.display_bins(h, -10., 10.);
        assert(full.first == 1 && full.second == h->GetNbinsX());
        assert(h->GetXaxis()->GetBinWidth(1) == 2.);
    }
    TH2D cluster_xy("h_debug_cluster_xy_run1", "NPS cluster position;x [cm];y [cm]",
                    34, -34., 34., 40, -40., 40.);
    cluster_xy.SetDirectory(nullptr);
    diagnostics.register_after(&cluster_xy);
    for (double x : {-35., -33., -33., -33., 0., 33., 35.}) {
        cluster_xy.Fill(x, 1.);
        diagnostics.stage_pass = x == 0.;
        diagnostics.fill(&cluster_xy, x, 1., 1.);
    }
    const std::string output_dir = argc > 1 ? argv[1] : "/tmp";
    TCanvas maps("maps", "Separate before/after regression", 1600, 1000);
    maps.Divide(2, 2);
    for (int side = 0; side < 2; ++side) {
        maps.cd(side + 1);
        gPad->SetLeftMargin(.12); gPad->SetRightMargin(.15);
        gPad->SetTopMargin(.16); gPad->SetBottomMargin(.12);
        auto* view = diagnostics.draw_2d(&cluster_xy, side == 1);
        assert(view != &cluster_xy && view != diagnostics.after2d.at(&cluster_xy));
        assert(view->GetEntries() == (side ? 1. : 7.));
        assert(view->GetXaxis()->GetBinLowEdge(view->GetXaxis()->GetFirst()) == -34.);
        assert(view->GetXaxis()->GetBinUpEdge(view->GetXaxis()->GetLast()) == 34.);
        assert(view->GetYaxis()->GetBinLowEdge(view->GetYaxis()->GetFirst()) == -40.);
        assert(view->GetYaxis()->GetBinUpEdge(view->GetYaxis()->GetLast()) == 40.);
        assert(view->GetXaxis()->GetBinWidth(1) == 2. && view->GetYaxis()->GetBinWidth(1) == 2.);
        assert(view->GetMinimum() == 0. && view->GetMaximum() == (side ? 1. : 3.));
        int histograms = 0;
        for (TObject* object : *gPad->GetListOfPrimitives())
            if (object->InheritsFrom(TH2::Class())) ++histograms;
        assert(histograms == 1);
        TLatex label; label.SetNDC(); label.SetTextSize(.04);
        label.DrawLatex(.15,.96,"TEST: NPS cluster x vs y");
        label.SetTextSize(.03);
        label.DrawLatex(.15,.90,side ? "After NPS cluster cuts" : "Before NPS cluster cuts");
        label.DrawLatex(.65,.90,Form("Entries %.0f",view->GetEntries()));
    }
    assert(cluster_xy.GetEntries() == 7. && cluster_xy.Integral(0,35,0,41) == 7.);
    assert(cluster_xy.GetBinContent(0,cluster_xy.GetYaxis()->FindBin(1.)) == 1.);
    assert(cluster_xy.GetBinContent(35,cluster_xy.GetYaxis()->FindBin(1.)) == 1.);
    maps.cd(3);
    gPad->SetTopMargin(.16); gPad->SetLeftMargin(.12); gPad->SetRightMargin(.05);
    before.SetStats(0); before.SetTitle(""); before.SetMaximum(2.);
    before.Draw("HIST"); diagnostics.draw_after(&before,"HMS cuts");
    TLatex caption; caption.SetNDC(); caption.SetTextSize(.04);
    caption.DrawLatex(.15,.96,"TEST: 1D legend layout");
    caption.SetTextSize(.03); caption.DrawLatex(.15,.90,"Cut: 0.4 < x < 1.1");
    maps.SaveAs((output_dir + "/separate_2d_test.png").c_str());
    // The after population may be empty; its view must still match before.
    TH2D sparse("h_debug_test_run1", "test;x;y", 100, -50., 50., 100, -50., 50.);
    sparse.SetDirectory(nullptr); diagnostics.register_after(&sparse); sparse.Fill(4.,6.,3.);
    maps.cd(4);
    auto* pre = diagnostics.draw_2d(&sparse,false);
    const int xfirst = pre->GetXaxis()->GetFirst(), ylast = pre->GetYaxis()->GetLast();
    gPad->Clear();
    auto* post = diagnostics.draw_2d(&sparse,true);
    assert(post->GetXaxis()->GetFirst() == xfirst && post->GetYaxis()->GetLast() == ylast);
    assert(post->GetMaximum() == 3. && post->GetEntries() == 0.);
    assert(sparse.GetXaxis()->GetFirst() == 1 && sparse.GetYaxis()->GetLast() == 100);
    nps2d::Params params;
    params.valid = params.ellipse_valid = true;
    params.fit_subset_bins = 4; params.fit_subset_total_fraction = .000303725;
    params.cov_det = 1e-12;
    nps2d::Config cfg;
    assert(npsplot::ellipse_failure(params,cfg).find("too few") != std::string::npos);
    // Same raw coordinates, mixed selectors and unequal weights. No histogram
    // clipping/reweighting may manufacture before/after distributions.
    std::vector<nps2d::Point> points = {{0,.135,.93,1.},{1,.135,.93,2.},{2,.135,.93,3.}};
    const auto original = points;
    nps2d::Result result;
    result.pass_ellipse = {0,1,0}; result.pass_mcd = {1,1,0};
    cfg.output_dir = argc > 1 ? argv[1] : "/tmp";
    cfg.tag = "mass_cut_diagnostic_test";
    npsplot::draw_mass_comparison(points, {1,0,1}, result, cfg, "ellipse fit unavailable");
    for (std::size_t i=0; i<points.size(); ++i) {
        assert(points[i].mpi0 == original[i].mpi0 && points[i].mmiss == original[i].mmiss);
        assert(points[i].weight == original[i].weight);
    }
    TFile file((cfg.output_dir+"/"+cfg.tag+".root").c_str());
    auto* all = file.Get<TH2D>("h_mass_compare_all");
    auto* decor = file.Get<TH2D>("h_mass_compare_decorrelation");
    assert(all && decor && all->Integral() == 6. && decor->Integral() == 4.);
    assert(all->GetBinError(all->FindBin(.135,.93)) == std::sqrt(14.));
    auto* mass_canvas = file.Get<TCanvas>("mass_cut_comparison");
    assert(mass_canvas);
    for (int pad = 1; pad <= 2; ++pad) {
        int histograms = 0;
        for (TObject* object : *mass_canvas->GetPad(pad)->GetListOfPrimitives())
            if (object->InheritsFrom(TH2::Class())) ++histograms;
        assert(histograms == 1);
    }
    std::cout << "PASS: separate 2D populations, independent NPS colors, shared spatial axes, original counts/flow/weights and selector fallback\n";
}
