#pragma once

// Serialize results, covariance and reproducibility metadata; close files owned by the analysis.
#include "xsec_analysis.h"
#include "xsec_migration_output.h"

inline void ExclPi0XSecAnalysis::write_results() {
    fout->cd();
    write_migration_results();
    // Save bin edges
    TVectorD vphi(phi_edges.size()), vt(tprime_edges.size()), vq(q2_edges.size()), vx(xb_edges.size());
    TMatrixD vxb2d(cfg.n_q2, cfg.n_xb + 1);
    for (size_t i = 0; i < phi_edges.size(); ++i) vphi[i] = phi_edges[i];
    for (size_t i = 0; i < tprime_edges.size(); ++i) vt[i] = tprime_edges[i];
    for (size_t i = 0; i < q2_edges.size(); ++i) vq[i] = q2_edges[i];
    for (size_t i = 0; i < xb_edges.size(); ++i) vx[i] = xb_edges[i];
    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        for (int ix = 0; ix <= cfg.n_xb; ++ix) {
            vxb2d(iq, ix) = xb_edges_by_q2[iq][ix];
        }
    }
    vphi.Write("phi_edges");
    vt.Write("tprime_edges");
    vq.Write("q2_edges");
    vx.Write("xb_edges");
    vxb2d.Write("xb_edges_by_q2");

    h_q2_data->Write();
    h_q2_sim->Write();
    h_xb_data->Write();
    h_xb_sim->Write();
    h_tprime_data->Write();
    h_tprime_sim->Write();
    h_phi_data->Write();
    h_phi_sim->Write();
    h_q2_xb_data->Write();
    h_q2_xb_sim->Write();
    h_tprime_phi_data->Write();
    h_tprime_phi_sim->Write();
    if (h_mass_data_all) {
        h_mass_data_all->Write(); h_mass_data_selected->Write();
        h_mass_sim_all->Write(); h_mass_sim_selected->Write();
        h_mmiss_corr_data_all->Write(); h_mmiss_corr_data_selected->Write();
        h_mmiss_reco_data_all->Write(); h_mmiss_reco_data_selected->Write();
        h_mmiss_sim_all->Write(); h_mmiss_sim_selected->Write();
    }
    for (const auto& m : mmiss_slices) {
        m.data->Write(); m.exclusive->Write();
    }

    TParameter<Long64_t>("n_data_total", cutflow.n_data_total).Write();
    TParameter<Long64_t>("n_data_pass", cutflow.n_data_pass).Write();
    TParameter<Long64_t>("n_data_inrange", cutflow.n_data_inrange).Write();
    TParameter<Long64_t>("n_sim_total", cutflow.n_sim_total).Write();
    TParameter<Long64_t>("n_sim_pass", cutflow.n_sim_pass).Write();
    TParameter<Long64_t>("n_sim_inrange", cutflow.n_sim_inrange).Write();
    TParameter<Long64_t>("n_sim_zero_sigcm", cutflow.n_sim_zero_sigcm).Write();
    TParameter<Long64_t>("n_sim_vertex_fallback", cutflow.n_sim_vertex_fallback).Write();

    if (cfg.partons_projection) {
        // Save every slice's model point beside its extracted coefficients.
        // ROOT and CSV retain raw SIMC sigcm units; only plot axes convert
        // to nb/GeV^2. A zero partons_ok marks an uncomputed/invalid point.
        TTree model_tree("partons_projection", "GK06/GPDGK19 pi0 points at response-weighted generated-bin means");
        int it = 0, iq = 0, ix = 0, model_ok = 0;
        double q2 = 0, xb = 0, t = 0, tprime = 0, eps = 0, flux = 0;
        double fit_u = 0, fit_lt = 0, fit_tt = 0, model_u = 0, model_lt = 0, model_tt = 0;
        model_tree.Branch("it", &it); model_tree.Branch("iq", &iq); model_tree.Branch("ix", &ix);
        model_tree.Branch("partons_ok", &model_ok);
        model_tree.Branch("Q2", &q2); model_tree.Branch("xB", &xb);
        model_tree.Branch("t", &t); model_tree.Branch("tprime", &tprime);
        model_tree.Branch("epsilon", &eps); model_tree.Branch("electron_flux_xbq2", &flux);
        model_tree.Branch("fit_sigmaU", &fit_u); model_tree.Branch("fit_sigmaLT", &fit_lt);
        model_tree.Branch("fit_sigmaTT", &fit_tt);
        model_tree.Branch("partons_sigmaU", &model_u); model_tree.Branch("partons_sigmaLT", &model_lt);
        model_tree.Branch("partons_sigmaTT", &model_tt);
        for (it = 0; it < cfg.n_tprime; ++it) {
            for (iq = 0; iq < cfg.n_q2; ++iq) {
                for (ix = 0; ix < cfg.n_xb; ++ix) {
                    const SliceResult& s = slice(it, iq, ix);
                    model_ok = s.partons_ok ? 1 : 0;
                    q2 = s.mean_q2_vertex_sim; xb = s.mean_xb_vertex_sim;
                    t = s.mean_t_sim; tprime = s.mean_tprime_vertex_sim;
                    eps = s.partons_epsilon; flux = s.partons_electron_flux_xbq2;
                    fit_u = s.fit_xsec.sigmaU; fit_lt = s.fit_xsec.sigmaTL;
                    fit_tt = s.fit_xsec.sigmaTT;
                    model_u = s.partons_sigmaU; model_lt = s.partons_sigmaLT;
                    model_tt = s.partons_sigmaTT;
                    model_tree.Fill();
                }
            }
        }
        model_tree.Write();
    }

    // write slice fit summaries as TObjString for portability
    std::ostringstream meta;
    meta << std::setprecision(std::numeric_limits<double>::max_digits10);
    meta << "tprime_edges:";
    for (double x : tprime_edges) meta << " " << x;
    meta << "\nq2_edges:";
    for (double x : q2_edges) meta << " " << x;
    meta << "\nxb_edges:";
    for (double x : xb_edges) meta << " " << x;
    meta << "\nxb_edges_by_q2:";
    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        meta << "\n  iq=" << iq << ":";
        for (double x : xb_edges_by_q2[iq]) meta << " " << x;
    }
    meta << "\nphi_edges:";
    for (double x : phi_edges) meta << " " << x;
    meta << "\nconfigured_kinematic=" << cfg.configured_kinematic;
    meta << "\ninput_data_file=" << cfg.data_file;
    meta << "\ninput_simc_file=" << cfg.simc_file;
    meta << "\nmodel_xsec_branch=" << (has_model_xsec ? model_xsec_branch : "NONE");
    meta << "\nresponse_mode=global_generated_to_reconstructed_forward_fit";
    meta << "\nfit_subset_recovery=" << (fit_fallback ? "yes" : "no");
    meta << "\nfit_retained_q2_xb_groups=" << successful_fit_groups;
    meta << "\nfit_subset_policy=exclude_whole_reco_Q2_xB_groups_retain_all_truth_feed_in_as_free_nuisance";
    meta << "\nfit_subset_selection=largest_supported_converged_subset_first_tie_by_increasing_excluded_group_index";
    meta << "\nfit_subset_caveat=data_dependent_exclusion_requires_closure_covariance_is_conditional_on_selected_subset";
    meta << "\nresponse_equation=y_reco_r=sum_truth_b_component_a_A_rba_sigma_ba";
    meta << "\nresponse_columns=active_truth_block_times_U_LT_TT_integrated_event_basis_no_probability_normalization";
    meta << "\ntruth_region_guards=six_disjoint_exterior_faces_tprime_below_above_then_Q2_below_above_then_xB_below_above";
    meta << "\ntruth_guard_treatment=only_populated_guards_have_free_U_LT_TT_no_fixed_generator_background_no_edge_clamping";
    meta << "\ntruth_guard_limitation=piecewise_constant_exterior_shapes_require_guard_variation_and_closure";
    meta << "\nfit_objective=" << cfg.fit_objective;
    meta << "\nscaled_poisson_upstream_pi0_weight=fixed_background_correction";
    meta << "\nscaled_poisson_mc_stat=fixed_response";
    meta << "\nscaled_poisson_reference_min_neff=20";
    meta << "\nscaled_poisson_empty_scale_choice=" << cfg.scaled_empty_scale;
    meta << "\nscaled_poisson_minuit_status=" << scaled_minuit_status;
    meta << "\nscaled_poisson_covariance_status=" << scaled_covariance_status;
    meta << "\nscaled_poisson_edm=" << scaled_edm;
    meta << "\nscaled_poisson_calls=" << scaled_calls;
    meta << "\nfit_method=" << (cfg.fit_objective=="scaled-poisson" ? "Minuit2_scaled_Poisson_fixed_response_continuous_phi_positivity" : (cfg.positive_xsec ? "full_rank_SVD_then_continuous_angular_positivity_constrained_WLS" :
        "full_rank_whitened_column_normalized_SVD_no_regularization_no_positivity_clipping"));
    meta << "\nfit_positive_xsec=" << (cfg.positive_xsec ? "yes" : "no");
    meta << "\nfit_positivity_scope=all_active_truth_blocks_all_phi_epsilon_up_to_observed_block_maximum";
    meta << "\nfit_positivity_boundary_active=" << (positivity_boundary_active ? "yes" : "no");
    meta << "\nfit_positivity_zero_allowed=yes_no_artificial_positive_floor";
    meta << "\nfit_variance_mode=" << (cfg.fit_objective=="scaled-poisson" ? "ignored_fixed_response" : cfg.fit_variance_mode);
    meta << "\nfit_variance=" << (cfg.fit_objective=="scaled-poisson" ?
        "not_used_fixed_response_scaled_Poisson" : (cfg.fit_variance_mode == "data" ?
        "data_sumw2_after_target_divisor_squared_Eq_5_23" :
        "data_sumw2_plus_iterated_Poissonized_MC_event_outer_products_extension"));
    meta << "\nfit_covariance=" << (positivity_boundary_active ? "unavailable_boundary_constrained_estimator_NaN" :
        cfg.fit_objective=="scaled-poisson" ? "Minuit2_Hessian_physical_transformation_conditional_fixed_response" :
        "conditional_known_variance_inverse_information_no_chi2_ndf_rescaling");
    meta << "\nfit_curvature_inverse=" << (cfg.fit_objective=="scaled-poisson" ?
        "Minuit2_Hessian_physical_transform_diagnostic_not_boundary_interval" :
        "unconstrained_inverse_information_at_final_weights_diagnostic_not_boundary_covariance");
    meta << "\npositive_plot_errors=" << (cfg.fit_objective=="scaled-poisson" ?
        "unavailable_at_boundary_no_Gaussian_refit_toys" :
        "conditional_constrained_refit_toy_sampling_SD_not_coverage_interval");
    meta << "\npositive_refit_toys_requested=256";
    meta << "\npositive_refit_toys_successful=" << positive_toys_successful;
    meta << "\npositive_refit_toy_limitations=frozen_response_and_binning_final_input_row_variances_MC_effect_as_Gaussian_row_noise_target_and_model_systematics_separate";
    meta << "\nfitted_prediction_variance=" << (positivity_boundary_active ? "unavailable_boundary_constrained_estimator_NaN" :
        cfg.fit_objective=="scaled-poisson" ? "X_Cminuit_Xt_diag_fixed_response_parameter_only" :
        (cfg.fit_variance_mode=="data" ? "XCXt_diag_fixed_response" :
         "XCXt_diag_plus_Vmc_minus_2_Hdiag_Vmc_conditional_first_order_same_fit"));
    meta << "\nfit_ndf_interpretation=" << (cfg.positive_xsec ?
        "rows_minus_parameters_nominal_only_no_boundary_chi_square_calibration" : "rows_minus_parameters");
    meta << "\nlegacy_ratio_error=data_error_over_fitted_prediction_residual_display_only_not_ratio_uncertainty";
    meta << "\nfit_mc_covariance_limitation=fixed_generated_count_multinomial_correlations_not_included";
    meta << "\nfit_target_covariance=separate_rank_one_sigma_i_sigma_j_times_target_relative_error_squared";
    meta << "\nfit_rank_tolerance=" << cfg.rank_tolerance;
    meta << "\nfit_mc_max_iterations=" << cfg.mc_max_iterations;
    meta << "\nfit_mc_tolerance=" << cfg.mc_fit_tolerance;
    meta << "\nfit_mc_iterations=" << mc_iterations;
    meta << "\nfit_mc_converged=" << (cfg.fit_objective=="scaled-poisson" || cfg.fit_variance_mode=="data" ? "not_applicable" :
        (mc_converged ? "yes" : "no"));
    meta << "\nfit_excluded_supported_zero_variance_rows=" << omitted_zero_variance_rows;
    meta << "\nfit_empty_row_policy=" << (cfg.fit_objective=="scaled-poisson" ?
        "include_supported_zero_yield_reference_scale_from_block_then_slice_then_global" :
        "exclude_zero_observed_variance_no_pseudocounts_requires_sparse_bin_closure");
    meta << "\nfit_nonphysical_truth_bins=" << nonphysical_truth_bins;
    meta << "\nfit_chi2_and_ndf_scope=" << (cfg.fit_objective=="scaled-poisson" ?
        "deviance_and_nominal_rows_minus_parameters_global_not_reduced_chi_square" :
        "global_shared_by_all_slice_rows_not_independent_per_slice_fits");
    meta << "\nlegacy_phi_csv_semantics=data_sim_ratio_are_reconstructed_residuals_xsec_is_generated_bin_fitted_curve_at_phi_center";
    meta << "\nexperimental_point_equation=sigma_fit_reference_times_(data_minus_all_other_truth_phi_cell_contributions)_over_corresponding_truth_phi_cell_contribution";
    meta << "\nexperimental_point_correspondence=same_tprime_Q2_xB_phi_grid_indices_at_reconstruction_and_vertex";
    meta << "\nexperimental_point_reference=accepted_vertex_phi_cell_response_weighted_Q2_xB_tprime_epsilon_and_phi_bin_center_function_evaluation_not_phi_average";
    meta << "\nexperimental_point_covariance=" << (cfg.fit_objective=="scaled-poisson" ?
        "unavailable_central_postfit_points_only" : (positivity_boundary_active ?
        "unavailable_boundary_constrained_estimator" :
        (cfg.fit_variance_mode == "data" ? "analytic_same_data_fixed_response_fixed_variance" :
         "analytic_same_data_plus_first_order_Poissonized_MC_response_and_reference_epsilon_fixed_final_variances")));
    meta << "\nexperimental_point_covariance_limitations=" << (cfg.fit_objective=="scaled-poisson" ?
        "no_analytic_point_covariance_central_postfit_diagnostics_only" :
        "conditional_weighted_data_variance_fixed_binning_and_final_GLS_variances_no_upstream_signal_weight_or_fixed_Ngen_or_detector_systematics");
    meta << "\nexperimental_point_target_covariance=separate_rank_one_common_divisor";
    meta << "\neta_correction=omitted_eta_equals_one";
    meta << "\noutput_cross_section=virtual_photon_d2sigma_dt_dphi_sigcm_units_not_four_fold_electron";
    meta << "\nU_approximation=bin_constant_T_plus_epsilon_L_with_event_varying_epsilon_no_Rosenbluth_separation";
    meta << "\nlegacy_mean_t_sim_and_epsilon=generated_origin_response_weighted_means";
    meta << "\ntruth_reference_means=accepted_reconstruction_response_weighted_generated_coordinates_not_bin_centering_correction";
    meta << "\ndiamond_xb_q2_vertices=" << xsec_diamond_vertices_text(cfg);
    meta << "\nexclusive_selection=" << cfg.mmiss_select;
    meta << "\nsim_reconstructed_mass_selection="
         << ((cfg.mmiss_select == "mcd" || cfg.mmiss_select == "ellipse") ?
             "combined_data_geometry" : "configured_missing_mass_window");
    if (cfg.mmiss_select == "mcd" || cfg.mmiss_select == "ellipse")
        meta << "\ncombined_mass_cut_file=" << cfg.mmiss_cut_file;
    meta << "\nmmiss_lower_gev=" << cfg.mmiss_lower_gev;
    meta << "\nmmiss_upper_gev=" << cfg.mmiss_upper_gev;
    meta << "\ntarget_contam_factor=" << cfg.tgt_contam;
    meta << "\ntarget_contam_factor_err=" << cfg.tgt_contam_err;
    meta << "\ntarget_contam_usage=data_yield_divided_by_factor";
    meta << "\nresponse_absolute_scale=physical_SIMC_normalization_only_no_data_area_matching";
    meta << "\npartons_projection=" << (cfg.partons_projection ? "yes" : "no");
    if (cfg.partons_projection) {
        meta << "\npartons_model=DVMPProcessGK06_DVMPCFFGK06_GPDGK19_LO";
        meta << "\npartons_mc_warmups=" << cfg.partons_warmups;
        meta << "\npartons_mc_calls=" << cfg.partons_calls;
        meta << "\npartons_kinematics=response_weighted_generated_bin_mean_Q2_xB_signed_physical_t";
        meta << "\npartons_electron_flux_xbq2=FULL_Hand_Gamma_equals_2pi_times_native_GK06_response_prefactor";
        meta << "\npartons_conversion=electron_nb_divided_by_FULL_Hand_Gamma_times_1e_minus9_then_Fourier_coefficients_times_2pi";
        meta << "\npartons_phi_projection=three_points_0_pi_over_2_pi_to_sigmaU_sigmaLT_sigmaTT";
        meta << "\npartons_warning=point_predictions_no_bin_or_detector_folding_MC_integration_error_not_propagated_LT_phi_sign_requires_check_zero_GK_tmin_guard_marked_unavailable";
    }
    meta << "\nvertex_kinematics_available=" << (has_vertex_kinematics ? "yes" : "no");
    meta << "\nvertex_kinematics_policy=mandatory_no_reconstructed_fallback";
    meta << "\nvertex_source=" << cfg.vertex_simc_file;
    meta << "\nvertex_epsilon_source=original_exclusive_h10_Q2i_Wi_hsxptari_hsyptari_SIMC_electron_angle";
    meta << "\nvertex_electron_hms_theta_deg=" << cfg.hms_theta_deg;
    meta << "\nvertex_raw_sigcm_matched_events=" << n_vertex_matched;
    meta << "\nvertex_lookup_note=raw_h10_entry_is_smeared_exclusive_event_id_sigcm_checked_for_each_selected_event";
    meta << "\nvertex_t_convention=SIMC_ti_positive_minus_t_converted_to_signed_physical_t_for_PARTONS";
    meta << "\nvertex_response_note=generated_Q2_xB_tprime_select_truth_block_vertex_epsilon_phi_define_basis_reconstructed_coordinates_select_row";
    meta << "\nbackground_note=full_data_mmiss_compared_to_area_normalized_smeared_generated_exclusive_SIMC_no_background_simulation_or_subtraction";
    meta << "\nmmiss_comparison_note=bin_residuals_and_pulls_are_diagnostics_user_selects_background_onset_boundary";
    meta << "\nfit_helicity_mode=unpolarized_total_only";
    meta << "\nfit_basis=(1/(2pi))*{sigmaU,sqrt(2eps(1+eps))*sigmaLT*cos(phi),eps*sigmaTT*cos(2phi)}";
    meta << "\ngamma_flux_usage=calculated_and_written_only_not_multiplied_in_fit";
    meta << "\ngamma_flux_note=full_weight_over_sigcm_retains_siglab_over_sigcm_equals_davejac_times_gtpr_times_fac";
    meta << "\nsim_weight_mode=full_weight_over_sigcm";
    meta << "\nhelicity_data=" << (has_helicity ? "yes" : "no");
    meta << "\nhelicity_sim=" << (has_sim_helicity ? "yes" : "no");
    TObjString(meta.str().c_str()).Write("analysis_metadata");
}

inline void ExclPi0XSecAnalysis::write_csv() {
    std::ofstream out(cfg.out_csv);
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    out << "it,iq,ix,ip,phi_lo,phi_hi,phi_center,q2_mean_data,xb_mean_data,tprime_mean_data,q2_mean_sim,xb_mean_sim,tprime_mean_sim,"
        << "epsilon,gamma_flux,data,data_err,sim,sim_err,ratio,ratio_err,xsec,xsec_err,xsec_sys_tgt,"
        << "q2_mean_truth,xb_mean_truth,tprime_mean_truth,t_mean_truth,truth_response_weight,fit_scope,fit_xsec_ok,fit_failure_reason,fit_objective,mc_stat_treatment\n";
    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                const SliceResult& s = slice(it, iq, ix);
                for (int ip = 0; ip < cfg.n_phi; ++ip) {
                    const PhiBin& pb = s.phi[ip];
                    double phi_lo = phi_edges[ip];
                    double phi_hi = phi_edges[ip+1];
                    double phi_c = 0.5 * (phi_lo + phi_hi);
                    out << it << "," << iq << "," << ix << "," << ip << ","
                        << phi_lo << "," << phi_hi << "," << phi_c << ","
                        << pb.mean_q2_data << "," << pb.mean_xb_data << "," << pb.mean_tprime_data << ","
                        << pb.mean_q2_sim << "," << pb.mean_xb_sim << "," << pb.mean_tprime_sim << ","
                        << s.epsilon << "," << s.gamma_flux << ","
                        << pb.data << "," << std::sqrt(std::max(0.0, pb.data_sumw2)) << ","
                        << pb.sim << "," << (s.fit_xsec.ok && !positivity_boundary_active ? std::sqrt(std::max(0.0, pb.sim_sumw2)) : std::numeric_limits<double>::quiet_NaN()) << ","
                        << pb.ratio << "," << pb.ratio_err << ","
                        << pb.xsec << "," << pb.xsec_err << "," << pb.xsec_sys_tgt << ","
                        << s.mean_q2_vertex_sim << "," << s.mean_xb_vertex_sim << "," << s.mean_tprime_vertex_sim << ","
                        << s.mean_t_sim << "," << s.truth_response_sum << ',' << s.fit_scope << ','
                        << s.fit_xsec.ok << ',' << xsec_csv_quote(s.fit_failure_reason) << ',' << cfg.fit_objective << ','
                        << (cfg.fit_objective=="scaled-poisson" ? "fixed_response" : cfg.fit_variance_mode) << '\n';
                }
            }
        }
    }
}

inline void ExclPi0XSecAnalysis::write_slice_csv() {
    std::ofstream out(cfg.out_slice_csv);
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    out << "it,iq,ix,tprime_lo,tprime_hi,tprime_center,q2_lo,q2_hi,xb_lo,xb_hi,"
        << "mean_q2_data,mean_xb_data,mean_tprime_data,mean_q2_sim,mean_xb_sim,mean_tprime_sim,"
        << "epsilon,gamma_flux,has_model_xsec,sumw_data,sumw_sim,"
        << "fit_ratio_ok,fit_ratio_chi2,fit_ratio_ndf,fit_ratio_A,fit_ratio_Aerr,fit_ratio_B,fit_ratio_Berr,fit_ratio_C,fit_ratio_Cerr,"
        << "fit_xsec_ok,fit_xsec_chi2,fit_xsec_ndf,fit_xsec_sigmaU,fit_xsec_sigmaUerr,fit_xsec_sigmaTL,fit_xsec_sigmaTLerr,fit_xsec_sigmaTT,fit_xsec_sigmaTTerr,"
        << "fit_asym_ok,fit_asym_chi2,fit_asym_ndf,fit_asym_sigmaTLp,fit_asym_sigmaTLperr,"
        // Preserve historical column names, but fitted coefficients now belong
        // to generated bins. Reconstructed means remain residual diagnostics;
        // vertex means and mean_t_sim define the generated-origin model point.
        << "mean_t_sim,mean_q2_vertex_sim,mean_xb_vertex_sim,partons_ok,partons_epsilon,partons_electron_flux_xbq2,"
        << "partons_sigmaU,partons_sigmaLT,partons_sigmaTT,mean_tprime_vertex_sim,truth_response_weight,"
        << "fit_scope,fit_global_rank,fit_global_condition,fit_mc_iterations,fit_mc_converged,fit_failure_reason,fit_objective,fit_objective_value,mc_stat_treatment\n";
    for (int it = 0; it < cfg.n_tprime; ++it) {
        for (int iq = 0; iq < cfg.n_q2; ++iq) {
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                const SliceResult& s = slice(it, iq, ix);
                double tlo = tprime_edges[it], thi = tprime_edges[it+1];
                double qlo = q2_edges[iq], qhi = q2_edges[iq+1];
                double xlo = xb_edges_by_q2[iq][ix], xhi = xb_edges_by_q2[iq][ix+1];
                double tcenter = 0.5 * (tlo + thi);
                out << it << "," << iq << "," << ix << ","
                    << tlo << "," << thi << "," << tcenter << ","
                    << qlo << "," << qhi << ","
                    << xlo << "," << xhi << ","
                    << s.mean_q2_data << "," << s.mean_xb_data << "," << s.mean_tprime_data << ","
                    << s.mean_q2_sim << "," << s.mean_xb_sim << "," << s.mean_tprime_sim << ","
                    << s.epsilon << "," << s.gamma_flux << "," << (s.has_model_xsec ? 1 : 0) << ","
                    << s.sumw_data << "," << s.sumw_sim << ","
                    << (s.fit_ratio.ok ? 1 : 0) << "," << s.fit_ratio.chi2 << "," << s.fit_ratio.ndf << ","
                    << (s.fit_ratio.p.size() > 0 ? s.fit_ratio.p[0] : 0.0) << "," << (s.fit_ratio.perr.size() > 0 ? s.fit_ratio.perr[0] : 0.0) << ","
                    << (s.fit_ratio.p.size() > 1 ? s.fit_ratio.p[1] : 0.0) << "," << (s.fit_ratio.perr.size() > 1 ? s.fit_ratio.perr[1] : 0.0) << ","
                    << (s.fit_ratio.p.size() > 2 ? s.fit_ratio.p[2] : 0.0) << "," << (s.fit_ratio.perr.size() > 2 ? s.fit_ratio.perr[2] : 0.0) << ","
                    << (s.fit_xsec.ok ? 1 : 0) << "," << (cfg.fit_objective=="scaled-poisson" ? std::numeric_limits<double>::quiet_NaN() : s.fit_xsec.chi2) << "," << s.fit_xsec.ndf << ","
                    << s.fit_xsec.sigmaU << "," << s.fit_xsec.sigmaU_err << ","
                    << s.fit_xsec.sigmaTL << "," << s.fit_xsec.sigmaTL_err << ","
                    << s.fit_xsec.sigmaTT << "," << s.fit_xsec.sigmaTT_err << ","
                    << (s.fit_asym.ok ? 1 : 0) << "," << s.fit_asym.chi2 << "," << s.fit_asym.ndf << ","
                    << s.fit_asym.sigmaTLp << "," << s.fit_asym.sigmaTLp_err << ","
                    << s.mean_t_sim << "," << s.mean_q2_vertex_sim << ","
                    << s.mean_xb_vertex_sim << "," << (s.partons_ok ? 1 : 0) << ","
                    << s.partons_epsilon << "," << s.partons_electron_flux_xbq2 << ","
                    << s.partons_sigmaU << "," << s.partons_sigmaLT << ","
                    << s.partons_sigmaTT << "," << s.mean_tprime_vertex_sim << "," << s.truth_response_sum
                    << ',' << s.fit_scope << ',' << s.fit_rank << ',' << s.fit_condition << ','
                    << s.fit_mc_iterations << ',' << (cfg.fit_variance_mode=="data"?-1:(s.fit_xsec.ok && mc_converged?1:0))
                    << ',' << xsec_csv_quote(s.fit_failure_reason) << ',' << cfg.fit_objective << ',' << s.fit_xsec.chi2
                    << ',' << (cfg.fit_objective=="scaled-poisson" ? "fixed_response" : cfg.fit_variance_mode) << '\n';
            }
        }
    }
}

inline void ExclPi0XSecAnalysis::write_mmiss_comparison_csv() {
    const fs::path path = fs::path(cfg.out_dir) / "mmiss_exclusive_comparison.csv";
    std::ofstream out(path);
    out << std::setprecision(10);
    out << "it,iq,ix,mmiss_lo,mmiss_hi,data_yield,data_error,"
           "sim_exclusive_yield,sim_exclusive_error,exclusive_normalization_scale,"
           "scaled_sim_exclusive_yield,scaled_sim_exclusive_error,data_minus_exclusive,pull,in_xsec_region\n";
    for (int it = 0; it < cfg.n_tprime; ++it)
        for (int iq = 0; iq < cfg.n_q2; ++iq)
            for (int ix = 0; ix < cfg.n_xb; ++ix) {
                const auto& m = mmiss_slices[static_cast<size_t>(slice_index(it, iq, ix))];
                for (int bin = 1; bin <= m.data->GetNbinsX(); ++bin) {
                    const double data = m.data->GetBinContent(bin);
                    const double data_err = m.data->GetBinError(bin);
                    const double sim = m.exclusive->GetBinContent(bin);
                    const double sim_err = m.exclusive->GetBinError(bin);
                    const double scaled = m.exclusive_normalization_scale * sim;
                    const double scaled_err = m.exclusive_normalization_scale * sim_err;
                    const double residual = data - scaled;
                    const double variance = data_err * data_err + scaled_err * scaled_err;
                    const double pull = variance > 0.0 ? residual / std::sqrt(variance) : 0.0;
                    out << it << ',' << iq << ',' << ix << ','
                        << m.data->GetXaxis()->GetBinLowEdge(bin) << ','
                        << m.data->GetXaxis()->GetBinUpEdge(bin) << ','
                        << data << ',' << data_err << ',' << sim << ',' << sim_err << ','
                        << m.exclusive_normalization_scale << ',' << scaled << ',' << scaled_err << ','
                        << residual << ',' << pull << ','
                        << (cfg.mmiss_select == "window" ?
                            (m.data->GetBinCenter(bin) >= cfg.mmiss_lower_gev &&
                             m.data->GetBinCenter(bin) <= cfg.mmiss_upper_gev ? 1 : 0) :
                            -1) << '\n'; // 2D selection cannot be inferred from a 1D bin
                }
            }
}

inline void ExclPi0XSecAnalysis::cleanup() {

    close_combined_pdf();

    if (fout) { fout->Write(); fout->Close(); fout = nullptr; }
    if (f_sim) { f_sim->Close(); f_sim = nullptr; }
    if (f_data) { f_data->Close(); f_data = nullptr; }
    if (f_vertex) { f_vertex->Close(); f_vertex = nullptr; }
}
