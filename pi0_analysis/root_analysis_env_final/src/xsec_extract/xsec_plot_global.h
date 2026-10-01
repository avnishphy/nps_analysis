#pragma once

// Global QA, missing-mass shape overlays and PDF lifecycle; display scaling never enters the fit.
#include "xsec_plot_style.h"
#include "xsec_analysis.h"
#include <spawn.h>
#include <sys/wait.h>
#include <unistd.h>

extern char** environ;

// Keep legacy per-canvas filenames and append the same canvas to the optional
// multipage PDF. Plot generation never changes the fitted response matrix.
inline void ExclPi0XSecAnalysis::write_canvas_pdf_png(TCanvas* c, const std::string& base) {
    if (!c || (!cfg.write_pdf && !cfg.write_png)) return;
    if (cfg.write_pdf) {
        const std::string page = base + ".pdf";
        c->SaveAs(page.c_str());
        append_to_combined_pdf(page);
    }
    if (cfg.write_png) c->SaveAs((base + ".png").c_str());
}

// The numeric edges are shared by the truth and reconstructed grids. The
// caller's title/axis must identify which grid its observable belongs to.
inline std::string ExclPi0XSecAnalysis::format_slice_bin_label(int it, int iq, int ix) const {
    if (it < 0 || it >= cfg.n_tprime || iq < 0 || iq >= cfg.n_q2 || ix < 0 || ix >= cfg.n_xb) {
        return "Invalid bin index";
    }

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3)
        << "it=" << it << " iq=" << iq << " ix=" << ix
        << " | Q^{2}[" << q2_edges[iq] << "," << q2_edges[iq + 1] << "]"
        << " x_{B}[" << xb_edges_by_q2[iq][ix] << "," << xb_edges_by_q2[iq][ix + 1] << "]"
        << " t'[" << tprime_edges[it] << "," << tprime_edges[it + 1] << "]";
    return oss.str();
}

// Annotate the current pad with grid indices and explicit numeric boundaries.
inline void ExclPi0XSecAnalysis::draw_slice_bin_label(int it, int iq, int ix, double y_ndc) const {
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextFont(42);
    lat.SetTextSize(0.026);
    lat.DrawLatex(0.11, y_ndc, format_slice_bin_label(it, iq, ix).c_str());
}

// Resolve a relative multipage filename beneath this run's output directory.
inline void ExclPi0XSecAnalysis::init_combined_pdf() {
    if (!cfg.write_pdf) return;
    fs::path p(cfg.out_all_plots_pdf);
    if (p.is_relative()) p = fs::path(cfg.out_dir) / p;
    combined_pdf_path = p.string();
    combined_pdf_pages.clear();
}

// Record the already-rendered page. Interleaving two live ROOT PDF streams
// caused TPDF object-state warnings and made page rendering order fragile.
inline void ExclPi0XSecAnalysis::append_to_combined_pdf(const std::string& page_path) {
    if (!cfg.write_pdf) return;
    if (combined_pdf_path.empty()) init_combined_pdf();
    if (combined_pdf_path.empty()) return;
    combined_pdf_pages.push_back(page_path);
}

// Merge complete single-page PDFs without reopening ROOT's PDF driver. argv
// avoids shell interpretation of output paths; a temporary target protects an
// existing combined PDF if merging fails.
inline void ExclPi0XSecAnalysis::close_combined_pdf() {
    if (!cfg.write_pdf || combined_pdf_pages.empty() || combined_pdf_path.empty()) return;
    const fs::path destination(combined_pdf_path);
    if (destination.has_parent_path()) fs::create_directories(destination.parent_path());
    const fs::path temporary = destination.string() + ".merge." + std::to_string(getpid()) + ".pdf";
    if (combined_pdf_pages.size() == 1) {
        fs::copy_file(combined_pdf_pages.front(), temporary, fs::copy_options::overwrite_existing);
    } else {
        std::vector<char*> arguments;
        arguments.push_back(const_cast<char*>("pdfunite"));
        for (const auto& page : combined_pdf_pages)
            arguments.push_back(const_cast<char*>(page.c_str()));
        const std::string temp_name = temporary.string();
        arguments.push_back(const_cast<char*>(temp_name.c_str()));
        arguments.push_back(nullptr);
        pid_t pid = -1;
        const int launch = posix_spawnp(&pid, "pdfunite", nullptr, nullptr, arguments.data(), environ);
        if (launch != 0) die("Cannot launch pdfunite for combined plot PDF: error " + std::to_string(launch));
        int status = 0;
        if (waitpid(pid, &status, 0) != pid || !WIFEXITED(status) || WEXITSTATUS(status) != 0)
            die("pdfunite failed while building combined plot PDF: " + combined_pdf_path);
    }
    fs::rename(temporary, destination);
    combined_pdf_pages.clear();
    log("Wrote combined plot PDF: " + combined_pdf_path);
}

// Reconstructed occupancy and shape QA use the original model-weighted MC.
// Area matching here affects detached histogram clones only: these overlays
// are not absolute forward-fit predictions and do not normalize the extraction.
inline void ExclPi0XSecAnalysis::make_global_plots() {
    if (!cfg.diagnostics || (!cfg.write_pdf && !cfg.write_png)) return;
    fs::create_directories(fs::path(cfg.out_dir) / "global");

    // Global overlays should use the same target-contamination normalization
    // as the slice-level extraction to avoid apparent data/sim mismatches.
    const double contam_scale = (std::isfinite(cfg.tgt_contam) && cfg.tgt_contam > 0.0)
                                    ? (1.0 / cfg.tgt_contam)
                                    : 1.0;
    auto clone_scaled_data_hist = [&](const TH1D* src, const char* name) {
        TH1D* h = dynamic_cast<TH1D*>(src->Clone(name));
        if (h) {
            h->SetDirectory(nullptr);
            h->Scale(contam_scale);
        }
        return std::unique_ptr<TH1D>(h);
    };

    auto h_q2_data_corr = clone_scaled_data_hist(h_q2_data.get(), "h_q2_data_corr");
    auto h_xb_data_corr = clone_scaled_data_hist(h_xb_data.get(), "h_xb_data_corr");
    auto h_tprime_data_corr = clone_scaled_data_hist(h_tprime_data.get(), "h_tprime_data_corr");
    auto h_phi_data_corr = clone_scaled_data_hist(h_phi_data.get(), "h_phi_data_corr");
    h_q2_data_corr->GetXaxis()->SetTitle("Q^{2}_{reco} [GeV^{2}]");
    h_xb_data_corr->GetXaxis()->SetTitle("x_{B,reco}");
    h_tprime_data_corr->GetXaxis()->SetTitle("t'_{reco} [GeV^{2}]");
    h_phi_data_corr->GetXaxis()->SetTitle("#phi_{reco} [rad]");
    for (TH2D* h : {h_q2_xb_data.get(), h_q2_xb_sim.get()}) {
        h->GetXaxis()->SetTitle("Q^{2}_{reco} [GeV^{2}]");
        h->GetYaxis()->SetTitle("x_{B,reco}");
    }
    for (TH2D* h : {h_tprime_phi_data.get(), h_tprime_phi_sim.get()}) {
        h->GetXaxis()->SetTitle("t'_{reco} [GeV^{2}]");
        h->GetYaxis()->SetTitle("#phi_{reco} [rad]");
    }

    auto clone_area_scaled_sim_hist = [](const TH1D* src, const char* name, const TH1D* data_ref) {
        TH1D* h = dynamic_cast<TH1D*>(src->Clone(name));
        if (h) {
            h->SetDirectory(nullptr);
            const double sim_int = h->Integral();
            const double data_int = data_ref ? data_ref->Integral() : 0.0;
            if (sim_int > 0.0 && data_int > 0.0) h->Scale(data_int / sim_int);
        }
        return std::unique_ptr<TH1D>(h);
    };

    auto h_q2_sim_shape = clone_area_scaled_sim_hist(h_q2_sim.get(), "h_q2_sim_shape", h_q2_data_corr.get());
    auto h_xb_sim_shape = clone_area_scaled_sim_hist(h_xb_sim.get(), "h_xb_sim_shape", h_xb_data_corr.get());
    auto h_tprime_sim_shape = clone_area_scaled_sim_hist(h_tprime_sim.get(), "h_tprime_sim_shape", h_tprime_data_corr.get());
    auto h_phi_sim_shape = clone_area_scaled_sim_hist(h_phi_sim.get(), "h_phi_sim_shape", h_phi_data_corr.get());

    auto draw_shape_overlay = [](TH1D* hdata, TH1D* hsim, const char* panel_title) {
        style_current_pad(0.13,0.04,0.12,0.22);
        hdata->SetTitle(panel_title);
        hdata->SetLineWidth(2);
        hdata->SetLineColor(kBlack);
        hsim->SetLineWidth(2);
        hsim->SetLineColor(kRed);
        hdata->GetYaxis()->SetTitle("Weighted yield / bin (shape only)");
        const double ymax = std::max(hdata->GetMaximum(), hsim->GetMaximum());
        hdata->SetMinimum(0.0);
        hdata->SetMaximum((ymax > 0.0) ? 1.18 * ymax : 1.0);
        hdata->Draw("hist");
        hsim->Draw("hist same");
        hdata->Draw("hist same");
        auto leg = make_compact_legend(0.48, 0.81, 0.93, 0.92, 0.032);
        leg->AddEntry(hdata, "Data / target factor", "l");
        leg->AddEntry(hsim, "SIMC rescaled to data area", "l");
        leg->Draw();
    };

    TCanvas c1("c1", "global", 1200, 900);
    c1.Divide(2,2);

    c1.cd(1);
    draw_shape_overlay(h_q2_data_corr.get(), h_q2_sim_shape.get(), "Reconstructed Q^{2}: shape comparison");

    c1.cd(2);
    draw_shape_overlay(h_xb_data_corr.get(), h_xb_sim_shape.get(), "Reconstructed x_{B}: shape comparison");

    c1.cd(3);
    draw_shape_overlay(h_tprime_data_corr.get(), h_tprime_sim_shape.get(), "Reconstructed t': shape comparison");

    c1.cd(4);
    draw_shape_overlay(h_phi_data_corr.get(), h_phi_sim_shape.get(), "Reconstructed #phi: shape comparison");

    c1.Update();
    write_canvas_pdf_png(&c1, (fs::path(cfg.out_dir) / "global" / "global_1d_distributions").string());

    // --- DEBUG: Q2:xB bin grid overlay ---
    TCanvas c2("c2", "Q2:xB binning", 800, 700);
    h_q2_xb_data->Draw("COLZ");
    // Draw Q2 and xB bin edges
    for (size_t i = 1; i < q2_edges.size() - 1; ++i) {
        TLine* l = new TLine(q2_edges[i], cfg.xb_min, q2_edges[i], cfg.xb_max);
        l->SetLineColor(kBlue+2); l->SetLineStyle(2); l->SetLineWidth(3); l->Draw();
    }
    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        const double qlo = q2_edges[iq];
        const double qhi = q2_edges[iq + 1];
        const auto& xrow = xb_edges_by_q2[iq];
        for (size_t i = 1; i < xrow.size() - 1; ++i) {
            TLine* l = new TLine(qlo, xrow[i], qhi, xrow[i]);
            l->SetLineColor(kRed+2); l->SetLineStyle(2); l->SetLineWidth(3); l->Draw();
        }
    }
    c2.Update();
    write_canvas_pdf_png(&c2, (fs::path(cfg.out_dir) / "global" / "q2_xb_binning_debug").string());

    // --- DEBUG: configured t' bin edges ---
    TCanvas c3("c3", "t' binning", 800, 700);
    h_tprime_data_corr->GetYaxis()->SetTitle("Tgt-corrected weighted counts");
    h_tprime_data_corr->SetMinimum(0.0);
    h_tprime_data_corr->SetMaximum(1.15 * std::max(1e-12, h_tprime_data_corr->GetMaximum()));
    h_tprime_data_corr->Draw("hist");
    for (size_t i = 1; i < tprime_edges.size() - 1; ++i) {
        TLine* l = new TLine(tprime_edges[i], 0, tprime_edges[i], h_tprime_data_corr->GetMaximum());
        l->SetLineColor(kGreen+2); l->SetLineStyle(2); l->SetLineWidth(3); l->Draw();
    }
    c3.Update();
    write_canvas_pdf_png(&c3, (fs::path(cfg.out_dir) / "global" / "tprime_binning_debug").string());

    // --- DEBUG: phi bin edges ---
    TCanvas c4("c4", "phi binning", 800, 700);
    h_phi_data_corr->GetYaxis()->SetTitle("Tgt-corrected weighted counts");
    h_phi_data_corr->SetMinimum(0.0);
    h_phi_data_corr->SetMaximum(1.15 * std::max(1e-12, h_phi_data_corr->GetMaximum()));
    h_phi_data_corr->Draw("hist");
    for (size_t i = 1; i < phi_edges.size() - 1; ++i) {
        TLine* l = new TLine(phi_edges[i], 0, phi_edges[i], h_phi_data_corr->GetMaximum());
        l->SetLineColor(kMagenta+2); l->SetLineStyle(2); l->SetLineWidth(3); l->Draw();
    }
    c4.Update();
    write_canvas_pdf_png(&c4, (fs::path(cfg.out_dir) / "global" / "phi_binning_debug").string());

    // --- Existing 2D plots ---
    TCanvas c5("c5", "2d", 1200, 500);
    c5.Divide(2,1);
    c5.cd(1); style_current_pad(0.12,0.16,0.12,0.10);
    h_q2_xb_data->SetTitle("Data weighted yield before target divisor");h_q2_xb_data->Draw("COLZ");
    c5.cd(2); style_current_pad(0.12,0.16,0.12,0.10);
    h_tprime_phi_data->SetTitle("Data weighted yield before target divisor");h_tprime_phi_data->Draw("COLZ");
    c5.Update();
    write_canvas_pdf_png(&c5, (fs::path(cfg.out_dir) / "global" / "data_occupancy_2d").string());

    TCanvas c6("c6", "2d_sim", 1200, 500);
    c6.Divide(2,1);
    c6.cd(1); style_current_pad(0.12,0.16,0.12,0.10);
    h_q2_xb_sim->SetTitle("SIMC response weight (full_weight/sigcm)");h_q2_xb_sim->Draw("COLZ");
    c6.cd(2); style_current_pad(0.12,0.16,0.12,0.10);
    h_tprime_phi_sim->SetTitle("SIMC response weight (full_weight/sigcm)");h_tprime_phi_sim->Draw("COLZ");
    c6.Update();
    write_canvas_pdf_png(&c6, (fs::path(cfg.out_dir) / "global" / "simc_occupancy_2d").string());
}

// Missing-mass overlays diagnose reconstructed candidates before the mass cut.
// Their separate display-only area scale is never reused by the response fit.
inline void ExclPi0XSecAnalysis::make_mmiss_comparison_plots() {
    // Plots span 0-2.5 GeV even though only the configured mass window enters
    // the cross-section. Comparison remains visual; no background-channel
    // model or subtraction is applied.
    if (!cfg.diagnostics || (!cfg.write_pdf && !cfg.write_png)) return;
    const fs::path plot_dir = fs::path(cfg.out_dir) / "mmiss_comparison";
    fs::create_directories(plot_dir);
    for (int it = 0; it < cfg.n_tprime; ++it)
        for (int iq = 0; iq < cfg.n_q2; ++iq)
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                auto& m = mmiss_slices[static_cast<size_t>(slice_index(it, iq, ix))];
                const std::string suffix = "t" + std::to_string(it) + "_q" +
                    std::to_string(iq) + "_x" + std::to_string(ix);
                auto clone = [&](TH1D* src, const std::string& name, double factor) {
                    auto h = std::unique_ptr<TH1D>(dynamic_cast<TH1D*>(src->Clone(name.c_str())));
                    h->SetDirectory(nullptr);
                    h->Scale(factor);
                    return h;
                };
                auto he = clone(m.exclusive.get(), "plot_mmiss_exclusive_" + suffix,
                                m.exclusive_normalization_scale);
                TCanvas c(("c_mmiss_" + suffix).c_str(), "Missing-mass comparison", 1100, 700);
                c.SetLeftMargin(0.12); c.SetRightMargin(0.03);
                m.data->SetMarkerStyle(20); m.data->SetMarkerSize(0.65);
                m.data->SetLineColor(kBlack);
                he->SetLineColor(kBlue+1); he->SetLineWidth(2); he->SetLineStyle(2);
                const double ymax = 1.35 * std::max({m.data->GetMaximum(), he->GetMaximum(), 1e-6});
                m.data->SetMinimum(0.0); m.data->SetMaximum(ymax);
                m.data->SetTitle(("Reconstructed bin, full missing mass: " + format_slice_bin_label(it, iq, ix)).c_str());
                m.data->Draw("E1");
                he->Draw("HIST SAME");
                m.data->Draw("E1 SAME");
                TLine cut_lo(cfg.mmiss_lower_gev, 0.0, cfg.mmiss_lower_gev, ymax);
                TLine cut_hi(cfg.mmiss_upper_gev, 0.0, cfg.mmiss_upper_gev, ymax);
                TLegend leg(0.58, 0.69, 0.96, 0.89);
                leg.SetBorderSize(0); leg.SetFillStyle(0);
                leg.AddEntry(m.data.get(), "Data: all candidates", "lep");
                leg.AddEntry(he.get(), "Exclusive SIMC: area-scaled", "l");
                if (cfg.mmiss_select == "window") {
                    for (TLine* cut : {&cut_lo, &cut_hi}) {
                        cut->SetLineColor(kMagenta+2); cut->SetLineStyle(3);
                        cut->SetLineWidth(2); cut->Draw();
                    }
                    std::ostringstream cut_label;
                    cut_label << std::fixed << std::setprecision(2) << cfg.mmiss_lower_gev
                              << "-" << cfg.mmiss_upper_gev << " GeV fit window";
                    leg.AddEntry(&cut_hi, cut_label.str().c_str(), "l");
                }
                leg.Draw();
                // Keep the explanatory text within the left half of the pad;
                // long notes would cover the right-hand legend and hide data.
                TLatex note; note.SetNDC(); note.SetTextSize(0.025);
                note.DrawLatex(0.13, 0.88, Form("Display scale: data / SIMC integral in %.2f-%.2f GeV",
                                                 cfg.mmiss_lower_gev,cfg.mmiss_upper_gev));
                note.DrawLatex(0.13, 0.83, "Visual comparison only; no background subtraction");
                if (cfg.mmiss_select != "window")
                    note.DrawLatex(0.13, 0.78, ("Selected by " + cfg.mmiss_select +
                        "; see mmiss_selection_2d").c_str());
                c.Update();
                write_canvas_pdf_png(&c, (plot_dir / ("mmiss_" + suffix)).string());
            }
}

// Four views make the event decision visible for both measured and generated-exclusive
// reconstructed samples. The data axes are mmiss_all versus mpi0_all.
inline void ExclPi0XSecAnalysis::make_mass_selection_plot() {
    if (!h_mass_data_all || !cfg.diagnostics || (!cfg.write_pdf && !cfg.write_png)) return;
    TCanvas canvas("c_mass_selection", "2D mass selection", 1300, 1050);
    canvas.Divide(2, 2);
    TH2D* views[] = {h_mass_data_all.get(), h_mass_data_selected.get(),
                     h_mass_sim_all.get(), h_mass_sim_selected.get()};
    const char* titles[] = {"Data: all candidates", "Data: selected",
                            "Generated-exclusive SIMC: all reconstructed",
                            "Generated-exclusive SIMC: selected"};
    std::unique_ptr<TGraph> boundary;
    if (cfg.mmiss_select == "mcd" || cfg.mmiss_select == "ellipse") {
        const auto& g = mass_geometry;
        const double trace = g.cov_xx + g.cov_yy;
        const double delta = std::sqrt(0.25 * (g.cov_xx - g.cov_yy) *
            (g.cov_xx - g.cov_yy) + g.cov_xy * g.cov_xy);
        const double angle = 0.5 * std::atan2(2.0 * g.cov_xy, g.cov_xx - g.cov_yy);
        const double major = std::sqrt(g.d2_cut * (0.5 * trace + delta));
        const double minor = std::sqrt(g.d2_cut * (0.5 * trace - delta));
        boundary = std::make_unique<TGraph>(241);
        for (int i = 0; i <= 240; ++i) {
            const double a = 2.0 * TMath::Pi() * i / 240.0;
            const double u = major * std::cos(a), v = minor * std::sin(a);
            boundary->SetPoint(i, g.mean_x + u * std::cos(angle) - v * std::sin(angle),
                                g.mean_y + u * std::sin(angle) + v * std::cos(angle));
        }
        boundary->SetLineColor(kRed + 1);
        boundary->SetLineWidth(3);
    }
    for (int i = 0; i < 4; ++i) {
        canvas.cd(i + 1);
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.15);
        gPad->SetBottomMargin(0.12);
        views[i]->SetTitle(Form("%s (%lld plotted events);%s [GeV];%s [GeV]",
                                titles[i], static_cast<long long>(views[i]->Integral()),
                                i < 2 ? "mpi0_all" : "mpi0",
                                i < 2 ? "mmiss_all" : "mmiss"));
        views[i]->Draw("COLZ");
        if (boundary) boundary->Draw("L SAME");
        TLatex label;
        label.SetNDC(true);
        label.SetTextSize(0.035);
        label.DrawLatex(0.15, 0.87, ("Cut: " + cfg.mmiss_select).c_str());
    }
    canvas.Update();
    write_canvas_pdf_png(&canvas,
        (fs::path(cfg.out_dir) / "mmiss_selection_2d").string());
}

// Keep corrected data Mx on its own page. The next page compares the same
// reconstructed observable in data (mmiss_all) and SIMC (mmiss). These are
// unweighted diagnostic counts; no display scaling enters the fit.
inline void ExclPi0XSecAnalysis::make_mmiss_selection_1d_plot() {
    if (!h_mmiss_corr_data_all || !cfg.diagnostics ||
        (!cfg.write_pdf && !cfg.write_png)) return;

    struct MassPanel {
        const TH1D* all;
        const TH1D* selected;
        std::string title;
        std::string axis;
        std::string all_label;
        std::string selected_label;
        bool zoom;
    };
    auto render_page = [&](const char* name, const char* caption,
                           const std::vector<MassPanel>& panels) {
        TCanvas canvas(name, caption, 1400, panels.size() == 2 ? 600 : 1100);
        canvas.Divide(2, static_cast<int>(panels.size() / 2));
        std::vector<std::unique_ptr<TH1D>> all(panels.size()), selected(panels.size());
        std::vector<std::unique_ptr<TLegend>> legends(panels.size());
        for (size_t i = 0; i < panels.size(); ++i) {
            const auto& panel = panels[i];
            canvas.cd(static_cast<int>(i + 1));
            gPad->SetLeftMargin(0.12);
            gPad->SetRightMargin(0.04);
            gPad->SetBottomMargin(0.13);
            all[i].reset(static_cast<TH1D*>(panel.all->Clone(
                (std::string(name) + "_all_" + std::to_string(i)).c_str())));
            selected[i].reset(static_cast<TH1D*>(panel.selected->Clone(
                (std::string(name) + "_selected_" + std::to_string(i)).c_str())));
            all[i]->SetDirectory(nullptr);
            selected[i]->SetDirectory(nullptr);
            const double low = panel.zoom ? 0.65 : 0.0;
            const double high = panel.zoom ? 1.15 : 2.5;
            for (TH1D* h : {all[i].get(), selected[i].get()})
                h->GetXaxis()->SetRangeUser(low, high);
            const std::string title = panel.title + (panel.zoom ? ": peak zoom;" : ": full spectrum;") +
                panel.axis + ";Events / 20 MeV";
            all[i]->SetTitle(title.c_str());
            all[i]->SetLineColor(kBlue + 1);
            all[i]->SetLineWidth(2);
            selected[i]->SetLineColor(kOrange + 7);
            selected[i]->SetLineWidth(2);
            selected[i]->SetFillColorAlpha(kOrange + 1, 0.45);
            gPad->SetLogy(!panel.zoom);
            all[i]->SetMinimum(panel.zoom ? 0.0 : 0.8);
            // Scale from the visible bins so an outlying full-range bin
            // cannot hide the selected peak in a zoom panel.
            const int first_bin = all[i]->FindFixBin(low);
            const int last_bin = all[i]->FindFixBin(std::nextafter(high, low));
            double visible_max = 1.0;
            for (int bin = first_bin; bin <= last_bin; ++bin)
                visible_max = std::max({visible_max,
                    all[i]->GetBinContent(bin), selected[i]->GetBinContent(bin)});
            all[i]->SetMaximum(1.35 * visible_max);
            all[i]->Draw("HIST");
            selected[i]->Draw("HIST SAME");
            all[i]->Draw("HIST SAME");
            legends[i] = std::make_unique<TLegend>(
                panel.zoom ? 0.13 : 0.55, 0.75, panel.zoom ? 0.54 : 0.94, 0.89);
            legends[i]->SetBorderSize(0);
            legends[i]->SetFillStyle(0);
            legends[i]->AddEntry(all[i].get(), panel.all_label.c_str(), "l");
            legends[i]->AddEntry(selected[i].get(), panel.selected_label.c_str(), "lf");
            legends[i]->Draw();
        }
        canvas.Update();
        write_canvas_pdf_png(&canvas, (fs::path(cfg.out_dir) / name).string());
    };

    const std::string data_selected = "Data: combined " + cfg.mmiss_select + " selected";
    const std::string sim_selected = "SIMC: data " + cfg.mmiss_select + " selected";
    render_page("mmiss_selection_1d", "Data corrected Mx selection", {
        {h_mmiss_corr_data_all.get(), h_mmiss_corr_data_selected.get(),
         "Data corrected Mx", "mmiss_all_corr [GeV]", "All data candidates", data_selected, false},
        {h_mmiss_corr_data_all.get(), h_mmiss_corr_data_selected.get(),
         "Data corrected Mx", "mmiss_all_corr [GeV]", "All data candidates", data_selected, true},
    });
    render_page("mmiss_reconstructed_comparison", "Reconstructed data and SIMC Mx", {
        {h_mmiss_reco_data_all.get(), h_mmiss_reco_data_selected.get(),
         "Data reconstructed Mx", "mmiss_all [GeV]", "All data candidates", data_selected, false},
        {h_mmiss_reco_data_all.get(), h_mmiss_reco_data_selected.get(),
         "Data reconstructed Mx", "mmiss_all [GeV]", "All data candidates", data_selected, true},
        {h_mmiss_sim_all.get(), h_mmiss_sim_selected.get(),
         "SIMC reconstructed Mx", "mmiss [GeV]", "Generated-exclusive SIMC", sim_selected, false},
        {h_mmiss_sim_all.get(), h_mmiss_sim_selected.get(),
         "SIMC reconstructed Mx", "mmiss [GeV]", "Generated-exclusive SIMC", sim_selected, true},
    });
}

// Display the generated-bin reference epsilon used to reconstruct the angular
// cross-section curves. The fit itself integrates event-level vertex epsilon
// in every matrix element; these plotted averages do not replace that integral.
inline void ExclPi0XSecAnalysis::make_epsilon_plots() {
    if (!cfg.diagnostics || (!cfg.write_pdf && !cfg.write_png)) return;
    fs::create_directories(fs::path(cfg.out_dir) / "global");

    TCanvas c_eps_t("c_eps_t", "epsilon_vs_tprime", 1200, 900);
    TH1D hframe("hframe_eps_t", "Generated-bin reference #epsilon;#LTt'_{gen}#GT [GeV^{2}];#epsilon", 100, cfg.tprime_min, cfg.tprime_max);
    hframe.SetMinimum(0.0);
    hframe.SetMaximum(1.05);
    hframe.Draw("AXIS");

    std::vector<std::unique_ptr<TGraphErrors>> graphs;
    TLegend leg(0.52, 0.58, 0.88, 0.88);
    leg.SetBorderSize(0);
    leg.SetFillStyle(0);

    const int colors[] = {kBlue + 1, kRed + 1, kGreen + 2, kMagenta + 2, kOrange + 7, kCyan + 2};

    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        for (int ix = 0; ix < cfg.n_xb; ++ix) {
            std::vector<double> xt, yeps, ey;
            xt.reserve(cfg.n_tprime);
            yeps.reserve(cfg.n_tprime);
            ey.reserve(cfg.n_tprime);

            for (int it = 0; it < cfg.n_tprime; ++it) {
                const SliceResult& s = slice(it, iq, ix);
                if (!s.fit_xsec.ok || !std::isfinite(s.epsilon)) continue;
                xt.push_back(s.mean_tprime_vertex_sim);
                yeps.push_back(s.epsilon);
                ey.push_back(0.0);
            }

            if (xt.empty()) continue;

            auto g = std::make_unique<TGraphErrors>(static_cast<int>(xt.size()), xt.data(), yeps.data(), nullptr, ey.data());
            int series = iq * cfg.n_xb + ix;
            g->SetLineColor(colors[series % 6]);
            g->SetMarkerColor(colors[series % 6]);
            g->SetMarkerStyle(20 + (series % 10));
            g->SetLineWidth(2);
            g->Draw("PL SAME");

            std::ostringstream lbl;
            lbl << "Q^{2}[" << std::fixed << std::setprecision(2) << q2_edges[iq] << "," << q2_edges[iq + 1]
                << "], x_{B}[" << xb_edges_by_q2[iq][ix] << "," << xb_edges_by_q2[iq][ix + 1] << "]";
            leg.AddEntry(g.get(), lbl.str().c_str(), "lp");
            graphs.emplace_back(std::move(g));
        }
    }

    leg.Draw();
    TLatex lat;
    lat.SetNDC(true);
    lat.SetTextSize(0.032);
    lat.DrawLatex(0.12, 0.93, Form("E_{beam} = %.3f GeV", cfg.ebeam));
    c_eps_t.Update();
    write_canvas_pdf_png(&c_eps_t, (fs::path(cfg.out_dir) / "global" / "epsilon_vs_tprime_by_q2_xb").string());

    const int nslices = cfg.n_tprime * cfg.n_q2 * cfg.n_xb;
    TH1D h_eps_slice("h_eps_slice", "Reference #epsilon by generated bin;generated (t',Q^{2},x_{B}) bin index;#epsilon", nslices, 0.5, nslices + 0.5);
    int bin = 1;
    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                const SliceResult& s = slice(it, iq, ix);
                h_eps_slice.SetBinContent(bin, s.epsilon);
                std::ostringstream bl;
                bl << std::fixed << std::setprecision(2)
                   << "t[" << tprime_edges[it] << "," << tprime_edges[it + 1] << "] "
                   << "Q2[" << q2_edges[iq] << "," << q2_edges[iq + 1] << "] "
                   << "xB[" << xb_edges_by_q2[iq][ix] << "," << xb_edges_by_q2[iq][ix + 1] << "]";
                h_eps_slice.GetXaxis()->SetBinLabel(bin, bl.str().c_str());
                ++bin;
            }
        }
    }

    TCanvas c_eps_idx("c_eps_idx", "epsilon_by_slice", 1300, 700);
    h_eps_slice.SetMinimum(0.0);
    h_eps_slice.SetMaximum(1.05);
    h_eps_slice.SetStats(0);
    h_eps_slice.SetMarkerStyle(20);
    h_eps_slice.SetLineWidth(2);
    h_eps_slice.GetXaxis()->SetLabelSize(0.018);
    h_eps_slice.GetXaxis()->LabelsOption("v");
    h_eps_slice.Draw("P HIST");
    c_eps_idx.SetBottomMargin(0.22);
    c_eps_idx.Update();
    write_canvas_pdf_png(&c_eps_idx, (fs::path(cfg.out_dir) / "global" / "epsilon_by_slice_index").string());
}
