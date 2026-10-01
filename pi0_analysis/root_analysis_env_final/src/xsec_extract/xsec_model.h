#pragma once

// Optional GK predictions and comparison plots. Predictions do not determine data fit normalization.
#include "xsec_plot_global.h"

// Evaluate an external theory reference after fitting. Each reference uses
// response-weighted generated-origin means for the matching truth bin. GK is
// not integrated over that bin or folded through the detector; its point value
// must not be interpreted as a fitted bin average or a normalization constraint.
// In particular, epsilon(mean Q2,mean xB) need not equal the response-weighted
// mean epsilon used for the fitted angular display. Their difference is not
// a normalization bug: the reference model and extracted bin constant are
// different estimands until a finite-bin/model-variation study connects them.
inline void ExclPi0XSecAnalysis::compute_partons_projection() {
    if (!cfg.partons_projection) return;
    if (successful_fit_groups == 0) {
        warn("Skipping PARTONS: no successful global fit; writing failure diagnostics");
        return;
    }
#ifdef NPS_ENABLE_PARTONS
    // PARTONS requires physical generated t, not t'=t-tmin. The moments
    // aggregate all accepted reconstructed destinations of each truth bin.
    nps_partons_pi0::Model model;
    model.initialize(cfg.partons_executable_path, cfg.partons_warmups, cfg.partons_calls);
    int projected = 0;
    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                SliceResult& s = slice(it, iq, ix);
                if (!s.fit_xsec.ok || s.truth_response_sum <= 0.0) continue;
                try {
                    const auto p = model.predict(s.mean_q2_vertex_sim, s.mean_xb_vertex_sim,
                                                 s.mean_t_sim, cfg.ebeam);
                    if (!p.valid) {
                        warn("PARTONS returned an unavailable pi0 prediction (including a zero from the GK tmin guard) for slice " +
                             std::to_string(it) + "/" + std::to_string(iq) + "/" +
                             std::to_string(ix));
                        continue;
                    }
                    s.partons_ok = true;
                    s.partons_epsilon = p.epsilon;
                    s.partons_electron_flux_xbq2 = p.electron_flux_xbq2;
                    s.partons_sigmaU = p.sigmaU;
                    s.partons_sigmaLT = p.sigmaLT;
                    s.partons_sigmaTT = p.sigmaTT;
                    ++projected;
                } catch (const std::exception& e) {
                    warn("PARTONS pi0 projection failed for slice " +
                         std::to_string(it) + "/" + std::to_string(iq) + "/" +
                         std::to_string(ix) + ": " + e.what());
                }
            }
        }
    }
    if (projected == 0) die("PARTONS produced no valid pi0 slice predictions.");
    log("PARTONS GK06/GPDGK19 projected " + std::to_string(projected) +
        " generated bins; theory evaluated at response-weighted generated means, not bin averages.");
#else
    die("--partons requires compilation with NPS_ENABLE_PARTONS and the native PARTONS libraries.");
#endif
}

// Display bin-constant fitted coefficients at their generated reference means.
// Horizontal markers carry no kinematic measurement uncertainty; bin intervals
// are recorded in the CSV. Vertical errors are marginal fit uncertainties and
// do not display the nonzero cross-bin covariance or common target systematic.
inline void ExclPi0XSecAnalysis::make_partons_projection_plots() {
    if (!cfg.partons_projection || (!cfg.write_pdf && !cfg.write_png)) return;
    fs::create_directories(fs::path(cfg.out_dir) / "partons_projection");

    // Raw fitted and PARTONS coefficients are microbarn/MeV^2 in the SIMC
    // CM convention. A factor 10^9 gives readable nb/GeV^2 on the plots;
    // CSV and ROOT values remain raw so a reviewer can audit the conversion.
    constexpr double to_nb_per_gev2 = 1.0e9;
    const std::array<const char*, 3> names{"sigmaU", "sigmaLT", "sigmaTT"};
    const std::array<const char*, 3> titles{"#sigma_{U} = #sigma_{T} + #varepsilon#sigma_{L}",
                                            "#sigma_{LT}", "#sigma_{TT}"};
    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        for (int ix = 0; ix < cfg.n_xb; ++ix) {
            for (int term = 0; term < 3; ++term) {
                std::vector<double> data_x, data_ex_low, data_ex_high, data_y, data_ey;
                std::vector<double> model_x, model_y;
                double xmin = std::numeric_limits<double>::infinity();
                double xmax = -std::numeric_limits<double>::infinity();
                double ymin = std::numeric_limits<double>::infinity();
                double ymax = -std::numeric_limits<double>::infinity();
                for (int it = 0; it < cfg.n_tprime; ++it) {
                    const SliceResult& s = slice(it, iq, ix);
                    if (!s.fit_xsec.ok) continue;
                    const double center = -s.mean_tprime_vertex_sim;
                    const double fitted = (term == 0 ? s.fit_xsec.sigmaU :
                                           term == 1 ? s.fit_xsec.sigmaTL : s.fit_xsec.sigmaTT) * to_nb_per_gev2;
                    const double fitted_err = to_nb_per_gev2*std::max(0.,
                        plot_coefficient_error(slice_index(it,iq,ix),term));
                    data_x.push_back(center);
                    data_ex_low.push_back(std::max(0.,center+tprime_edges[it+1]));
                    data_ex_high.push_back(std::max(0.,-tprime_edges[it]-center));
                    data_y.push_back(fitted); data_ey.push_back(fitted_err);
                    xmin = std::min(xmin, -tprime_edges[it + 1]);
                    xmax = std::max(xmax, -tprime_edges[it]);
                    ymin = std::min(ymin, fitted - fitted_err);
                    ymax = std::max(ymax, fitted + fitted_err);
                    if (s.partons_ok) {
                        const double xm = -s.mean_tprime_vertex_sim;
                        const double ym = (term == 0 ? s.partons_sigmaU :
                                           term == 1 ? s.partons_sigmaLT : s.partons_sigmaTT) * to_nb_per_gev2;
                        model_x.push_back(xm); model_y.push_back(ym);
                        xmin = std::min(xmin, xm); xmax = std::max(xmax, xm);
                        ymin = std::min(ymin, ym); ymax = std::max(ymax, ym);
                    }
                }
                if (data_x.empty()) continue;
                const double xspan = std::max(1e-4, xmax - xmin);
                const double yspan = std::max(1.0, ymax - ymin);
                std::ostringstream tag;
                tag << names[term] << "_q" << iq << "_x" << ix;
                TCanvas canvas(("c_partons_" + tag.str()).c_str(), "PARTONS pi0 projection", 1100, 820);
                canvas.cd();
                style_current_pad(0.16, 0.04, 0.13, 0.18);
                TH1D frame(("h_partons_frame_" + tag.str()).c_str(),
                           ";-t'_{gen} [GeV^{2}];Structure function [nb/GeV^{2}]", 100,
                           xmin - 0.08 * xspan, xmax + 0.08 * xspan);
                frame.SetMinimum(ymin - 0.18 * yspan);
                frame.SetMaximum(ymax + 0.24 * yspan);
                frame.Draw();
                style_hist_axes(&frame, 0.045, 0.04);
                frame.GetXaxis()->CenterTitle();
                frame.GetXaxis()->SetNdivisions(505);
                if (ymin<0. && ymax>0.) {
                    auto zero=new TLine(frame.GetXaxis()->GetXmin(),0.,frame.GetXaxis()->GetXmax(),0.);
                    zero->SetLineStyle(2);zero->SetLineColor(kGray+2);zero->Draw();
                }

                TGraphAsymmErrors fitted_graph(static_cast<int>(data_x.size()), data_x.data(),
                    data_y.data(), data_ex_low.data(), data_ex_high.data(),data_ey.data(),data_ey.data());
                fitted_graph.SetMarkerStyle(20); fitted_graph.SetMarkerSize(1.0);
                fitted_graph.SetMarkerColor(kBlack); fitted_graph.SetLineColor(kBlack);
                fitted_graph.Draw("PZ SAME");
                TGraph model_graph(static_cast<int>(model_x.size()), model_x.data(), model_y.data());
                model_graph.SetMarkerStyle(25); model_graph.SetMarkerSize(1.2);
                model_graph.SetMarkerColor(kMagenta + 2); model_graph.SetLineColor(kMagenta + 2);
                if (!model_x.empty()) model_graph.Draw("P SAME");
                // Keep the legend below the kinematics and projection caveats.
                // The fitted points are displayed in nb/GeV^2 after conversion.
                auto legend = make_compact_legend(0.48, 0.69, 0.84, 0.79, 0.034);
                legend->AddEntry(&fitted_graph, positivity_boundary_active ?
                    (positive_toys_successful ? "Positive fit (refit-toy SD)" : "Positive fit (toy SD unavailable)") :
                    "Generated-bin coefficient", "p");
                if (!model_x.empty()) legend->AddEntry(&model_graph, "GK06/GPDGK19 at means", "p");
                legend->Draw();

                TLatex text;
                text.SetNDC(true); text.SetTextFont(42); text.SetTextAlign(22);
                text.SetTextSize(0.037);
                text.DrawLatex(0.50, 0.965, titles[term]);
                text.SetTextSize(0.029);
                text.DrawLatex(0.50, 0.925,
                    Form("Generated Q^{2} #in [%.3f, %.3f] GeV^{2}, x_{B} #in [%.3f, %.3f]",
                         q2_edges[iq], q2_edges[iq + 1],
                         xb_edges_by_q2[iq][ix], xb_edges_by_q2[iq][ix + 1]));
                text.SetTextSize(0.025);
                text.DrawLatex(0.50, 0.886, "Vertex means; horizontal bin spans; GK not detector-folded");
                text.DrawLatex(0.50, 0.850,Form("Objective: %s; positivity: %s; GK error omitted%s",
                    cfg.fit_objective.c_str(),positivity_boundary_active ? "boundary" :
                    (cfg.positive_xsec ? "on" : "off"),term==1 ? "; check LT #phi sign" : ""));
                canvas.Update();
                write_canvas_pdf_png(&canvas,
                    (fs::path(cfg.out_dir) / "partons_projection" /
                     ("partons_" + tag.str() + "_vs_minus_tprime")).string());
            }
        }
    }
}
