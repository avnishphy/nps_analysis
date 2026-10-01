#pragma once

// Defurne Fig. 3.9 (printed p.51) motivates displaying vertex kinematics for
// events accepted at reconstruction. That DIS figure uses vertex beam energy
// as its color. The present input has no such field: our adaptation explicitly
// uses selected MC counts and the configured Q2/xB boundaries instead.
#include "xsec_plot_global.h"
#include <TExec.h>

inline void ExclPi0XSecAnalysis::make_migration_coverage_plots() {
    if (!cfg.diagnostics || !h_migration_vertex_q2_xb || !h_migration_reco_q2_xb) return;
    fout->cd();
    auto* directory = fout->GetDirectory("migration");
    if (!directory) directory = fout->mkdir("migration");
    directory->cd();
    h_migration_vertex_q2_xb->Write();
    h_migration_reco_q2_xb->Write();
    TObjString("Adaptation of Defurne Fig.3.9: selected-MC counts, not vertex beam energy; "
               "dashed lines are reconstructed analysis bin boundaries, not the HRS acceptance contour; "
               "both stored histograms contain exactly the events used to build the response, including guards; "
               "display-only clones are divided by their xB-Q2 bin area for comparable color densities.")
        .Write("coverage_semantics");
    fout->cd();
    if (!cfg.write_pdf && !cfg.write_png) return;
    fs::create_directories(fs::path(cfg.out_dir) / "migration");

    std::vector<int> old_palette(gStyle->GetNumberOfColors());
    for (size_t i = 0; i < old_palette.size(); ++i) old_palette[i] = gStyle->GetColorPalette(i);
    const int old_contours = gStyle->GetNumberContours();
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(100);
    TCanvas canvas("c_migration_coverage", "Generated and reconstructed selected-MC coverage", 1500, 800);
    canvas.Divide(2, 1);
    // ROOT can extend truth axes farther than the reco axes, merging bins.
    // Compare density, not counts per potentially different bin area. Stored
    // histograms remain raw counts for an exact event-accounting check.
    auto density = [](TH2D* source, const char* name) {
        auto h = std::unique_ptr<TH2D>(static_cast<TH2D*>(source->Clone(name)));
        h->SetDirectory(nullptr);
        h->Scale(1.0, "width");
        h->SetTitle(""); // TLatex below provides the single panel heading.
        h->GetZaxis()->SetTitle("Selected MC entries / (#Deltax_{B} #DeltaQ^{2})");
        return h;
    };
    auto vertex_density = density(h_migration_vertex_q2_xb.get(), "vertex_coverage_density");
    auto reco_density = density(h_migration_reco_q2_xb.get(), "reco_coverage_density");
    const double maximum = std::max({1.0, vertex_density->GetMaximum(), reco_density->GetMaximum()});
    const double minimum = std::max(1e-12, maximum * 1e-4);
    std::array<TH2D*, 2> maps{vertex_density.get(), reco_density.get()};
    std::vector<std::unique_ptr<TLine>> boundaries;
    std::vector<std::unique_ptr<TExec>> palettes;
    for (int pad = 0; pad < 2; ++pad) {
        canvas.cd(pad + 1);
        style_current_pad(0.13, 0.17, 0.12, 0.18);
        gPad->SetLogz();
        auto* h = maps[pad];
        h->SetStats(false);
        h->SetMinimum(minimum);
        h->SetMaximum(maximum);
        h->Draw("AXIS");
        // ROOT paints pads again during PDF export. Retain a per-pad palette
        // command so another figure's signed palette cannot leak into this one.
        auto palette = std::make_unique<TExec>(("migration_coverage_palette_" + std::to_string(pad)).c_str(),
                                               "gStyle->SetPalette(kBird);");
        palette->Draw(); palettes.push_back(std::move(palette));
        h->Draw("COLZ SAME");
        style_hist_axes(h, 0.042, 0.032);
        h->GetZaxis()->SetTitleOffset(1.3);
        h->GetZaxis()->SetTitleSize(0.035);
        h->GetZaxis()->SetLabelSize(0.028);
        auto line = [&](double x1, double y1, double x2, double y2) {
            auto edge = std::make_unique<TLine>(x1, y1, x2, y2);
            edge->SetLineColor(kBlack); edge->SetLineStyle(2); edge->SetLineWidth(2);
            edge->Draw(); boundaries.push_back(std::move(edge));
        };
        // Include the outer analysis rectangle and the frozen internal bin
        // edges. xB edges may differ by Q2 interval, so draw each row locally.
        for (double q : q2_edges) line(cfg.xb_min, q, cfg.xb_max, q);
        for (int iq = 0; iq < cfg.n_q2; ++iq)
            for (double x : xb_edges_by_q2[iq]) line(x, q2_edges[iq], x, q2_edges[iq + 1]);
        TLatex note; note.SetNDC(); note.SetTextFont(42); note.SetTextSize(0.031);
        note.DrawLatex(0.13, 0.96, pad == 0 ? "Generated Q^{2}, x_{B}: accepted MC" : "Reconstructed Q^{2}, x_{B}: same MC");
        note.SetTextSize(0.025);
        note.DrawLatex(0.13, 0.91, "Dashed: reconstructed analysis bin boundaries");
        note.DrawLatex(0.13, 0.87, "Color: MC density (log); not vertex beam energy");
    }
    canvas.Update();
    write_canvas_pdf_png(&canvas, (fs::path(cfg.out_dir) / "migration" / "vertex_reco_q2_xb_coverage").string());
    if (!old_palette.empty()) gStyle->SetPalette(static_cast<int>(old_palette.size()), old_palette.data());
    gStyle->SetNumberContours(old_contours);
}
