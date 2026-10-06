#pragma once

// Reconstructed-yield residuals and generated-bin cross-section plots, with explicitly distinct axes.
#include "xsec_plot_global.h"
#include "xsec_helicity.h"

// Pads 1, 2 and 4 show reconstructed observations. Pad 3 instead displays the
// generated-bin angular function and separately marked residual-corrected
// experimental points. Identical grid indices do not imply identical events.
inline void ExclPi0XSecAnalysis::make_slice_plots() {
    if (!cfg.write_pdf && !cfg.write_png) return;
    fs::create_directories(fs::path(cfg.out_dir) / "slices");
    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                const SliceResult& s = slice(it, iq, ix);
                std::ostringstream tag;
                tag << "t" << it << "_q" << iq << "_x" << ix;
                TCanvas c(("c_"+tag.str()).c_str(), tag.str().c_str(), 1600, 1100);
                c.Divide(2, 2, 0.01, 0.02);

                if (!s.fit_xsec.ok) {
                    for (int pad = 1; pad <= 4; ++pad) {
                        c.cd(pad);
                        draw_pad_message("Excluded from global fit", "See fit_status.csv for the failure reason.");
                        if (pad == 1) draw_slice_bin_label(it, iq, ix);
                    }
                    c.Update();
                    write_canvas_pdf_png(&c, (fs::path(cfg.out_dir) / "slices" / ("slice_" + tag.str())).string());
                    continue;
                }

                std::vector<double> phi_c(cfg.n_phi), ydata(cfg.n_phi), yerr(cfg.n_phi), ysim(cfg.n_phi), ysimerr(cfg.n_phi);
                for (int ip = 0; ip < cfg.n_phi; ++ip) {
                    phi_c[ip] = 0.5 * (phi_edges[ip] + phi_edges[ip+1]);
                    const PhiBin& pb = s.phi[ip];
                    ydata[ip] = pb.data; yerr[ip] = std::sqrt(std::max(0.0, pb.data_sumw2));
                    ysim[ip] = pb.sim;
                    const int row=slice_index(it,iq,ix)*cfg.n_phi+ip;
                    const double toy_var=positivity_boundary_active ?
                        plot_parameter_variance(response_design[row]) : pb.sim_sumw2;
                    ysimerr[ip] = std::isfinite(toy_var) ? std::sqrt(std::max(0.0,toy_var)) : 0.0;
                }

                c.cd(1);
                style_current_pad(0.13, 0.04, 0.12, 0.10);
                auto gdata = new TGraphErrors(cfg.n_phi, phi_c.data(), ydata.data(), nullptr, yerr.data());
                auto gsim  = new TGraphErrors(cfg.n_phi, phi_c.data(), ysim.data(), nullptr, ysimerr.data());
                gdata->SetTitle("Reconstructed-bin forward-fit yields;#phi_{reco} [rad];Weighted counts");
                gdata->SetMarkerStyle(20); gdata->SetLineColor(kBlack); gdata->SetMarkerColor(kBlack);
                gsim->SetMarkerStyle(24); gsim->SetLineColor(kRed); gsim->SetMarkerColor(kRed);
                gdata->Draw("AP");
                style_graph_axes(gdata, 0.046, 0.04);
                gsim->Draw("P SAME");
                auto l1 = make_compact_legend(0.56, 0.76, 0.90, 0.89, 0.030);
                l1->AddEntry(gdata,"Data (tgt corrected)","p");
                l1->AddEntry(gsim, positivity_boundary_active ?
                            (positive_toys_successful ? "Positive fit (refit-toy SD)" : "Positive fit (toy SD unavailable)") :
                            "Global forward prediction", "p"); l1->Draw();
                draw_slice_bin_label(it, iq, ix);

                c.cd(2);
                style_current_pad(0.13, 0.04, 0.12, 0.10);
                std::vector<double> ratio_phi, ratio_y, ratio_err;
                double ratio_min=1., ratio_max=1.;
                for (int ip=0;ip<cfg.n_phi;++ip) {
                    const auto& pb=s.phi[ip];
                    if (!std::isfinite(pb.ratio) || !std::isfinite(pb.ratio_err) || pb.ratio_err<0.) continue;
                    ratio_phi.push_back(phi_c[ip]);ratio_y.push_back(pb.ratio);ratio_err.push_back(pb.ratio_err);
                    ratio_min=std::min(ratio_min,pb.ratio-pb.ratio_err);
                    ratio_max=std::max(ratio_max,pb.ratio+pb.ratio_err);
                }
                if (ratio_phi.empty()) {
                    draw_pad_message("Residual ratio unavailable", "No row has a positive forward prediction.");
                } else {
                    auto grat = new TGraphErrors(static_cast<int>(ratio_phi.size()),
                                                 ratio_phi.data(), ratio_y.data(), nullptr, ratio_err.data());
                    grat->SetTitle("Reconstructed-yield residual;#phi_{reco} [rad];Data / forward prediction");
                    grat->SetMarkerStyle(20);grat->SetMarkerColor(kBlack);grat->SetLineColor(kBlack);
                    const double span=std::max(0.2,ratio_max-ratio_min);
                    grat->SetMinimum(ratio_min-0.18*span);
                    grat->SetMaximum(ratio_max+0.28*span);
                    grat->Draw("AP");
                    grat->GetXaxis()->SetLimits(cfg.phi_min,cfg.phi_max);
                    style_graph_axes(grat, 0.046, 0.04);
                    auto unity=new TLine(cfg.phi_min,1.,cfg.phi_max,1.);
                    unity->SetLineColor(kGray+2);unity->SetLineStyle(2);unity->Draw();
                    grat->Draw("P SAME");
                    TLatex fit_note;
                    fit_note.SetNDC(true);fit_note.SetTextFont(42);fit_note.SetTextSize(0.028);
                    fit_note.DrawLatex(0.15,0.84,Form("Global %s = %.1f/%.0f",
                                                      "Objective / nominal DOF",
                                                      s.fit_xsec.chi2,s.fit_xsec.ndf));
                    fit_note.SetTextSize(0.024);
                    fit_note.DrawLatex(0.15,0.79,Form("Objective: %s; positivity: %s",
                        cfg.fit_objective.c_str(),positivity_boundary_active ? "boundary" :
                        (cfg.positive_xsec ? "on" : "off")));
                    fit_note.DrawLatex(0.15,0.74,Form("%zu/%d ratios; bars use data sumw^{2} only",
                        ratio_phi.size(),cfg.n_phi));
                }

                c.cd(3);
                style_current_pad(0.13, 0.04, 0.12, 0.10);
                if (s.fit_xsec.ok && s.fit_xsec.p.size() >= 3) {
                    constexpr double display_scale=1.0e9; // microbarn/MeV^2 to nb/GeV^2
                    const bool fit_tlp = (has_helicity && has_sim_helicity && s.fit_xsec.p.size() >= 4);
                    const double eps = clamp(s.epsilon, 0.0, 1.0);
                    // Evaluate the virtual-photon CM angular cross section at
                    // this generated bin's reference epsilon. No electron
                    // flux multiplies this cross-section display; that flux is
                    // already in the response used to predict measured yields.
                    const double inv2pi = 1.0 / (2.0 * TMath::Pi());
                    const double k_lt = std::sqrt(std::max(0.0, 2.0 * eps * (1.0 + eps)));
                    const double k_tlp = std::sqrt(std::max(0.0, 2.0 * eps * (1.0 - eps)));

                    std::vector<double> xsec(cfg.n_phi), xsecerr(cfg.n_phi);
                    std::vector<double> point_phi, point_sigma, point_error;
                    for (int ip = 0; ip < cfg.n_phi; ++ip) {
                        const PhiBin& pb = s.phi[ip];
                        xsec[ip] = pb.xsec*display_scale;
                        std::vector<double> gradient(migration_fit.parameters.size(),0.);
                        const auto found=std::find(active_truth_blocks.begin(),active_truth_blocks.end(),
                            slice_index(it,iq,ix));
                        if(found!=active_truth_blocks.end()) {
                            const auto basis=sigma_model_basis(phi_c[ip],s.epsilon,0.,0);
                            for(int term=0;term<3;++term)
                                gradient[3*(found-active_truth_blocks.begin())+term]=basis[term];
                        }
                        const double plot_var=plot_parameter_variance(gradient);
                        xsecerr[ip]=std::isfinite(plot_var) ? display_scale*std::sqrt(plot_var) : 0.;
                        const int row=slice_index(it,iq,ix)*cfg.n_phi+ip;
                        if (row<static_cast<int>(experimental_points.size()) &&
                            (experimental_points[row].status=="ok" ||
                             experimental_points[row].status=="central_only_boundary" ||
                             experimental_points[row].status=="central_only_joint")) {
                            point_phi.push_back(phi_c[ip]);
                            point_sigma.push_back(experimental_points[row].sigma_exp*display_scale);
                            const double point_sd=plot_point_error(row);
                            point_error.push_back(std::isfinite(point_sd) ? point_sd*display_scale : 0.);
                        }
                    }

                    const int ncurve = 361;
                    std::vector<double> ph_curve(ncurve);
                    std::vector<double> y_tot(ncurve), y_u(ncurve), y_lt(ncurve), y_tt(ncurve), y_tlp(ncurve, 0.0);
                    for (int i = 0; i < ncurve; ++i) {
                        const double ph = cfg.phi_min + (cfg.phi_max - cfg.phi_min) * (static_cast<double>(i) / static_cast<double>(ncurve - 1));
                        ph_curve[i] = ph;
                        y_u[i] = display_scale*inv2pi * s.fit_xsec.sigmaU;
                        y_lt[i] = display_scale*inv2pi * k_lt * s.fit_xsec.sigmaTL * std::cos(ph);
                        y_tt[i] = display_scale*inv2pi * eps * s.fit_xsec.sigmaTT * std::cos(2.0 * ph);
                        if (fit_tlp) y_tlp[i] = display_scale*inv2pi * k_tlp * s.fit_xsec.sigmaTLp * std::sin(ph);
                        y_tot[i] = y_u[i] + y_lt[i] + y_tt[i] + (fit_tlp ? y_tlp[i] : 0.0);
                    }
                    double ymin = std::numeric_limits<double>::infinity();
                    double ymax = -std::numeric_limits<double>::infinity();
                    for (int ip = 0; ip < cfg.n_phi; ++ip) {
                        if (!std::isfinite(xsec[ip])) continue;
                        ymin = std::min(ymin, xsec[ip] - std::max(0.0, xsecerr[ip]));
                        ymax = std::max(ymax, xsec[ip] + std::max(0.0, xsecerr[ip]));
                    }
                    for (int i = 0; i < ncurve; ++i) {
                        ymin = std::min({ymin, y_tot[i], y_u[i], y_lt[i], y_tt[i]});
                        ymax = std::max({ymax, y_tot[i], y_u[i], y_lt[i], y_tt[i]});
                        if (fit_tlp) {
                            ymin = std::min(ymin, y_tlp[i]);
                            ymax = std::max(ymax, y_tlp[i]);
                        }
                    }
                    for (size_t i=0;i<point_phi.size();++i) {
                        ymin=std::min(ymin,point_sigma[i]-point_error[i]);
                        ymax=std::max(ymax,point_sigma[i]+point_error[i]);
                    }
                    if (!std::isfinite(ymin) || !std::isfinite(ymax)) {
                        ymin = -1.0;
                        ymax = 1.0;
                    }
                    const double span = std::max(1e-12, ymax - ymin);

                    auto gx = new TGraphErrors(cfg.n_phi, phi_c.data(), xsec.data(), nullptr, xsecerr.data());
                    gx->SetTitle("Virtual-photon angular cross section;#phi_{gen} [rad];d^{2}#sigma/(dt d#phi) [nb/(GeV^{2} rad)]");
                    gx->SetMarkerStyle(20);
                    gx->SetMarkerColor(kBlack);
                    gx->SetLineColor(kBlack);
                    gx->SetMinimum(ymin - 0.14 * span);
                    gx->SetMaximum(ymax + 0.32 * span);
                    gx->Draw("AP");
                    style_graph_axes(gx, 0.046, 0.04);

                    auto g_tot = new TGraphErrors(ncurve, ph_curve.data(), y_tot.data(), nullptr, nullptr);
                    auto g_u = new TGraphErrors(ncurve, ph_curve.data(), y_u.data(), nullptr, nullptr);
                    auto g_lt = new TGraphErrors(ncurve, ph_curve.data(), y_lt.data(), nullptr, nullptr);
                    auto g_tt = new TGraphErrors(ncurve, ph_curve.data(), y_tt.data(), nullptr, nullptr);
                    TGraphErrors* g_tlp = nullptr;
                    if (fit_tlp) g_tlp = new TGraphErrors(ncurve, ph_curve.data(), y_tlp.data(), nullptr, nullptr);

                    g_tot->SetLineColor(kBlue + 1); g_tot->SetLineWidth(3);
                    g_u->SetLineColor(kBlack); g_u->SetLineStyle(2); g_u->SetLineWidth(2);
                    g_lt->SetLineColor(kRed + 1); g_lt->SetLineStyle(7); g_lt->SetLineWidth(2);
                    g_tt->SetLineColor(kGreen + 2); g_tt->SetLineStyle(9); g_tt->SetLineWidth(2);
                    if (g_tlp) {
                        g_tlp->SetLineColor(kMagenta + 2);
                        g_tlp->SetLineStyle(5);
                        g_tlp->SetLineWidth(2);
                    }

                    g_tot->Draw("L SAME");
                    g_u->Draw("L SAME");
                    g_lt->Draw("L SAME");
                    g_tt->Draw("L SAME");
                    if (g_tlp) g_tlp->Draw("L SAME");

                    TGraphErrors* g_points=nullptr;
                    if (!point_phi.empty()) {
                        g_points=new TGraphErrors(static_cast<int>(point_phi.size()),
                            point_phi.data(),point_sigma.data(),nullptr,point_error.data());
                        g_points->SetMarkerStyle(25);g_points->SetMarkerColor(kMagenta+2);
                        g_points->SetLineColor(kMagenta+2);g_points->Draw("P SAME");
                    }

                    auto l3 = make_compact_legend(0.63, 0.58, 0.91, 0.88, 0.030);
                    l3->AddEntry(gx, positivity_boundary_active ?
                                 (positive_toys_successful ? "Fit at #phi centers (refit-toy SD)" : "Fit at #phi centers (toy SD unavailable)") :
                                 "Fit at #phi centers #pm marginal error", "p");
                    l3->AddEntry(g_tot, "Fit at mean vertex #varepsilon", "l");
                    if (g_points) l3->AddEntry(g_points,!cfg.joint_plot_input.empty() ?
                        "Residual-corrected (central only)" : positivity_boundary_active ?
                        (positive_toys_successful ? "Residual-corrected (refit-toy SD)" : "Residual-corrected (toy SD unavailable)") :
                        "Residual-corrected points","p");
                    l3->AddEntry(g_u, "#sigma_{U} term", "l");
                    l3->AddEntry(g_lt, "#sigma_{LT} cos#phi term", "l");
                    l3->AddEntry(g_tt, "#sigma_{TT} cos2#phi term", "l");
                    if (g_tlp) l3->AddEntry(g_tlp, "#sigma_{TL'} sin#phi term", "l");
                    l3->Draw();

                    auto comp = new TPaveText(0.13, 0.70, 0.58, 0.86, "NDC");
                    comp->SetBorderSize(0);
                    comp->SetFillColorAlpha(kWhite, 0.85);
                    comp->SetTextAlign(12);
                    comp->SetTextFont(42);
                    comp->SetTextSize(0.032);
                    if (positivity_boundary_active) {
                        comp->AddText(Form("#sigma_{U}=%.3g, #sigma_{LT}=%.3g, #sigma_{TT}=%.3g",
                                           display_scale*s.fit_xsec.sigmaU, display_scale*s.fit_xsec.sigmaTL,
                                           display_scale*s.fit_xsec.sigmaTT));
                        comp->AddText(Form("Positivity boundary: %d conditional refit toys",positive_toys_successful));
                    } else if (fit_tlp) {
                        comp->AddText(Form("#sigma_{U}=%.3g#pm%.2g, #sigma_{LT}=%.3g#pm%.2g",
                                           display_scale*s.fit_xsec.sigmaU, display_scale*std::max(0.0, s.fit_xsec.sigmaU_err),
                                           display_scale*s.fit_xsec.sigmaTL, display_scale*std::max(0.0, s.fit_xsec.sigmaTL_err)));
                        comp->AddText(Form("#sigma_{TT}=%.3g#pm%.2g, #sigma_{TL'}=%.3g#pm%.2g",
                                           display_scale*s.fit_xsec.sigmaTT, display_scale*std::max(0.0, s.fit_xsec.sigmaTT_err),
                                           display_scale*s.fit_xsec.sigmaTLp, display_scale*std::max(0.0, s.fit_xsec.sigmaTLp_err)));
                    } else {
                        comp->AddText(Form("#sigma_{U}=%.3g#pm%.2g, #sigma_{LT}=%.3g#pm%.2g",
                                           display_scale*s.fit_xsec.sigmaU, display_scale*std::max(0.0, s.fit_xsec.sigmaU_err),
                                           display_scale*s.fit_xsec.sigmaTL, display_scale*std::max(0.0, s.fit_xsec.sigmaTL_err)));
                        comp->AddText(Form("#sigma_{TT}=%.3g#pm%.2g",
                                           display_scale*s.fit_xsec.sigmaTT, display_scale*std::max(0.0, s.fit_xsec.sigmaTT_err)));
                    }
                    comp->AddText(Form("nb/GeV^{2}; #LT#varepsilon_{vertex}#GT=%.3f",s.epsilon));
                    comp->Draw();
                } else {
                    draw_pad_message("Fit unavailable", "No valid global migration fit for this generated bin.");
                }

                c.cd(4);
                style_current_pad(0.13, 0.04, 0.12, 0.10);
                if (has_helicity) {
                    std::vector<double> asym_phi, asym, asymerr;
                    for (int ip = 0; ip < cfg.n_phi; ++ip) {
                        const PhiBin& pb = s.phi[ip];
                        const auto diagnostic = nps_xsec::helicity_yield_asymmetry(
                            pb.data_plus, pb.data_minus, pb.data_plus_sumw2, pb.data_minus_sumw2);
                        if (!diagnostic.valid) continue;
                        asym_phi.push_back(phi_c[ip]);
                        asym.push_back(diagnostic.value);
                        asymerr.push_back(diagnostic.error);
                    }
                    if (asym.empty()) {
                        draw_pad_message("Helicity diagnostic unavailable",
                            "No #phi bin has a valid two-helicity yield estimate.");
                    } else {
                        auto ga = new TGraphErrors(static_cast<int>(asym.size()), asym_phi.data(), asym.data(), nullptr, asymerr.data());
                        ga->SetTitle("Reconstructed helicity yield diagnostic;#phi_{reco} [rad];(Y_{+}-Y_{-})/(Y_{+}+Y_{-})");
                        ga->SetMarkerStyle(20);
                        ga->SetMarkerColor(kBlack);
                        ga->SetLineColor(kBlack);
                        ga->Draw("AP");
                        style_graph_axes(ga, 0.046, 0.04);

                        // The global fit has only helicity-even coefficients.
                        // Raw yield balance has no helicity-charge or beam-
                        // polarization correction and cannot determine LT'.
                        TLatex note;
                        note.SetNDC(true); note.SetTextFont(42); note.SetTextSize(0.030);
                        note.DrawLatex(0.13, 0.83, "Yield balance only; LT' is not fitted");
                        note.DrawLatex(0.13, 0.77, "No helicity-charge or polarization correction");
                    }
                } else {
                    draw_pad_message("No helicity branch", "TL' is not extracted for this input.");
                }

                c.Update();
                write_canvas_pdf_png(&c, (fs::path(cfg.out_dir) / "slices" / ("slice_" + tag.str())).string());
            }
        }
    }
}

// Plot bin-constant U/LT/TT parameters at response-weighted generated means.
// A marker's x location is a reference coordinate, not a bin-centering
// correction. Marginal vertical errors omit cross-bin covariance and the
// separately exported common target-contamination systematic.
inline void ExclPi0XSecAnalysis::make_sigma_vs_tprime_plots() {
    if (!cfg.write_pdf && !cfg.write_png) return;
    fs::create_directories(fs::path(cfg.out_dir) / "sigma_vs_tprime");

    constexpr double display_scale=1.0e9; // microbarn/MeV^2 to nb/GeV^2
    const bool have_tlp = false;

    auto update_range = [](double y, double ey, double& ymin, double& ymax) {
        if (!std::isfinite(y)) return;
        const double err = (std::isfinite(ey) && ey > 0.0) ? ey : 0.0;
        ymin = std::min(ymin, y - err);
        ymax = std::max(ymax, y + err);
    };

    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        for (int ix = 0; ix < cfg.n_xb; ++ix) {
            std::vector<double> x, ex_low, ex_high;
            std::vector<double> sigma_u, sigma_u_err;
            std::vector<double> sigma_tl, sigma_tl_err;
            std::vector<double> sigma_tt, sigma_tt_err;
            std::vector<double> sigma_tlp, sigma_tlp_err;

            x.reserve(cfg.n_tprime);
            ex_low.reserve(cfg.n_tprime); ex_high.reserve(cfg.n_tprime);
            sigma_u.reserve(cfg.n_tprime);
            sigma_u_err.reserve(cfg.n_tprime);
            sigma_tl.reserve(cfg.n_tprime);
            sigma_tl_err.reserve(cfg.n_tprime);
            sigma_tt.reserve(cfg.n_tprime);
            sigma_tt_err.reserve(cfg.n_tprime);
            sigma_tlp.reserve(cfg.n_tprime);
            sigma_tlp_err.reserve(cfg.n_tprime);

            for (int it = 0; it < cfg.n_tprime; ++it) {
                const SliceResult& s = slice(it, iq, ix);
                if (!s.fit_xsec.ok) continue;

                const double reference=-s.mean_tprime_vertex_sim;
                x.push_back(reference);
                // Horizontal whiskers show the generated bin span, not an
                // uncertainty on the accepted-response-weighted reference.
                ex_low.push_back(std::max(0.,reference+tprime_edges[it+1]));
                ex_high.push_back(std::max(0.,-tprime_edges[it]-reference));

                sigma_u.push_back(display_scale*s.fit_xsec.sigmaU);
                sigma_u_err.push_back(display_scale*std::max(0.0,plot_coefficient_error(slice_index(it,iq,ix),0)));
                sigma_tl.push_back(display_scale*s.fit_xsec.sigmaTL);
                sigma_tl_err.push_back(display_scale*std::max(0.0,plot_coefficient_error(slice_index(it,iq,ix),1)));
                sigma_tt.push_back(display_scale*s.fit_xsec.sigmaTT);
                sigma_tt_err.push_back(display_scale*std::max(0.0,plot_coefficient_error(slice_index(it,iq,ix),2)));
                sigma_tlp.push_back(display_scale*s.fit_xsec.sigmaTLp);
                sigma_tlp_err.push_back(display_scale*std::max(0.0, s.fit_xsec.sigmaTLp_err));
            }

            if (x.empty()) continue;

            std::ostringstream tag;
            tag << "q" << iq << "_x" << ix;
            TCanvas c(("c_sigma_vs_t_" + tag.str()).c_str(), "sigma_vs_tprime", 1100, 820);
            c.cd();
            style_current_pad(0.16, 0.04, 0.13, 0.18);
            double xmin = std::numeric_limits<double>::infinity();
            double xmax = -std::numeric_limits<double>::infinity();
            double ymin = std::numeric_limits<double>::infinity();
            double ymax = -std::numeric_limits<double>::infinity();
            for (size_t i = 0; i < x.size(); ++i) {
                xmin = std::min(xmin, x[i] - ex_low[i]);
                xmax = std::max(xmax, x[i] + ex_high[i]);
                update_range(sigma_u[i], sigma_u_err[i], ymin, ymax);
                update_range(sigma_tl[i], sigma_tl_err[i], ymin, ymax);
                update_range(sigma_tt[i], sigma_tt_err[i], ymin, ymax);
                if (have_tlp) update_range(sigma_tlp[i], sigma_tlp_err[i], ymin, ymax);
            }

            if (!(std::isfinite(xmin) && std::isfinite(xmax) && std::isfinite(ymin) && std::isfinite(ymax))) {
                continue;
            }

            const double xspan = std::max(1e-6, xmax - xmin);
            const double yspan = std::max(1.0, ymax - ymin);
            TH1D h_coeff_frame(("h_coeff_frame_" + tag.str()).c_str(),
                               ";-t'_{gen} [GeV^{2}];Structure function [nb/GeV^{2}]",
                               100,
                               xmin - 0.08 * xspan,
                               xmax + 0.08 * xspan);
            h_coeff_frame.SetMinimum(ymin - 0.18 * yspan);
            h_coeff_frame.SetMaximum(ymax + 0.24 * yspan);
            h_coeff_frame.Draw();
            style_hist_axes(&h_coeff_frame, 0.045, 0.04);
            h_coeff_frame.GetXaxis()->CenterTitle();
            h_coeff_frame.GetXaxis()->SetNdivisions(505);
            if (ymin<0. && ymax>0.) {
                auto zero=new TLine(h_coeff_frame.GetXaxis()->GetXmin(),0.,
                                    h_coeff_frame.GetXaxis()->GetXmax(),0.);
                zero->SetLineStyle(2);zero->SetLineColor(kGray+2);zero->Draw();
            }

            auto g_u = new TGraphAsymmErrors(static_cast<int>(x.size()), x.data(), sigma_u.data(), ex_low.data(), ex_high.data(), sigma_u_err.data(), sigma_u_err.data());
            std::vector<double> no_x_error(x.size(),0.);
            auto g_tl = new TGraphAsymmErrors(static_cast<int>(x.size()), x.data(), sigma_tl.data(), no_x_error.data(), no_x_error.data(), sigma_tl_err.data(), sigma_tl_err.data());
            auto g_tt = new TGraphAsymmErrors(static_cast<int>(x.size()), x.data(), sigma_tt.data(), no_x_error.data(), no_x_error.data(), sigma_tt_err.data(), sigma_tt_err.data());
            g_u->SetLineColor(kBlack); g_u->SetMarkerColor(kBlack); g_u->SetMarkerStyle(20); g_u->SetLineWidth(2);
            g_tl->SetLineColor(kRed + 1); g_tl->SetMarkerColor(kRed + 1); g_tl->SetMarkerStyle(21); g_tl->SetLineWidth(2);
            g_tt->SetLineColor(kBlue + 1); g_tt->SetMarkerColor(kBlue + 1); g_tt->SetMarkerStyle(22); g_tt->SetLineWidth(2);

            g_u->Draw("PZ SAME");
            g_tl->Draw("PZ SAME");
            g_tt->Draw("PZ SAME");

            TGraphAsymmErrors* g_tlp = nullptr;
            if (have_tlp) {
                g_tlp = new TGraphAsymmErrors(static_cast<int>(x.size()), x.data(), sigma_tlp.data(), no_x_error.data(), no_x_error.data(), sigma_tlp_err.data(), sigma_tlp_err.data());
                g_tlp->SetLineColor(kGreen + 2);
                g_tlp->SetMarkerColor(kGreen + 2);
                g_tlp->SetMarkerStyle(23);
                g_tlp->SetLineWidth(2);
                g_tlp->Draw("PZ SAME");
            }

            auto leg_coeff = make_compact_legend(0.72, 0.62, 0.90, 0.78, 0.037);
            leg_coeff->AddEntry(g_u, "#sigma_{U}", "p");
            leg_coeff->AddEntry(g_tl, "#sigma_{LT}", "p");
            leg_coeff->AddEntry(g_tt, "#sigma_{TT}", "p");
            if (g_tlp) leg_coeff->AddEntry(g_tlp, "#sigma_{TL'}", "p");
            leg_coeff->Draw();

            TLatex lat1;
            lat1.SetNDC(true);
            lat1.SetTextFont(42);
            lat1.SetTextAlign(22);
            lat1.SetTextSize(0.040);
            lat1.DrawLatex(0.50, 0.965, "Generated-bin coefficients from global migration fit");
            lat1.SetTextSize(0.030);
            lat1.DrawLatex(0.50, 0.925,
                           Form("Generated Q^{2} #in [%.3f, %.3f] GeV^{2}, x_{B} #in [%.3f, %.3f]",
                                q2_edges[iq], q2_edges[iq + 1], xb_edges_by_q2[iq][ix], xb_edges_by_q2[iq][ix + 1]));
            lat1.SetTextSize(0.026);
            lat1.DrawLatex(0.50, 0.885, positivity_boundary_active ?
                (positive_toys_successful ?
                    Form("Vertex means; U bars show bin spans; refit-toy SD (N=%d)",positive_toys_successful) :
                    "Vertex means; U bars show bin spans; toy SD unavailable") :
                "Vertex means; U bars show bin spans; marginal errors only");
            lat1.SetTextSize(0.023);
            lat1.DrawLatex(0.50,0.852,Form("Objective: %s; positivity: %s; target error separate",
                cfg.fit_objective.c_str(),positivity_boundary_active ? "boundary" :
                (cfg.positive_xsec ? "on" : "off")));

            c.Update();
            write_canvas_pdf_png(&c, (fs::path(cfg.out_dir) / "sigma_vs_tprime" / ("sigma_terms_vs_minus_tprime_" + tag.str())).string());
        }
    }
}
