#pragma once

// Fit-independent export after event accumulation and the target correction.
// No active-fit subset, fitted coefficient, or positivity result is required.
#include "xsec_analysis.h"

inline void ExclPi0XSecAnalysis::write_joint_inputs() {
    const auto open = [&](const char* name) {
        std::ofstream out(fs::path(cfg.out_dir) / name);
        if (!out) die(std::string("Cannot create joint input: ") + name);
        out << std::setprecision(std::numeric_limits<double>::max_digits10);
        return out;
    };
    const auto finish = [&](std::ofstream& out) {
        out.flush();
        if (!out) die("Failed writing joint input");
    };
    auto bins = open("joint_input_slices.csv");
    bins << "it,iq,ix,tprime_lo,tprime_hi,q2_lo,q2_hi,xb_lo,xb_hi\n";
    auto reco = open("migration_reco_rows.csv");
    reco << "reco_row,it,iq,ix,ip,phi_lo,phi_hi,data,data_variance\n";
    for (int it=0; it<cfg.n_tprime; ++it)
        for (int iq=0; iq<cfg.n_q2; ++iq)
            for (int ix=0; ix<cfg.n_xb; ++ix) {
                const int b = slice_index(it, iq, ix);
                bins << it << ',' << iq << ',' << ix << ',' << tprime_edges[it] << ','
                     << tprime_edges[it+1] << ',' << q2_edges[iq] << ',' << q2_edges[iq+1]
                     << ',' << xb_edges_by_q2[iq][ix] << ',' << xb_edges_by_q2[iq][ix+1] << '\n';
                for (int ip=0; ip<cfg.n_phi; ++ip) {
                    const auto& p = slices[b].phi[ip];
                    reco << b*cfg.n_phi+ip << ',' << it << ',' << iq << ',' << ix << ','
                         << ip << ',' << phi_edges[ip] << ',' << phi_edges[ip+1] << ','
                         << p.data << ',' << p.data_sumw2 << '\n';
                }
            }
    auto truth = open("migration_truth_blocks.csv");
    truth << "truth_block,events,epsilon_max\n";
    for (size_t b=0; b<truth_moments.size(); ++b)
        truth << b << ',' << truth_moments[b].events << ',' << truth_moments[b].epsilon_max << '\n';
    auto response = open("migration_response_cells.csv");
    response << "reco_row,truth_block,events,basis_U,basis_LT,basis_TT";
    for (const char* a : {"U", "LT", "TT"})
        for (const char* b : {"U", "LT", "TT"}) response << ",cov_" << a << '_' << b;
    response << '\n';
    for (size_t r=0; r<migration_response.size(); ++r)
        for (size_t b=0; b<migration_response[r].size(); ++b) {
            const auto& cell = migration_response[r][b];
            response << r << ',' << b << ',' << cell.events;
            for (double v : cell.basis) response << ',' << v;
            for (double v : cell.covariance) response << ',' << v;
            response << '\n';
        }
    for (auto* stream : {&bins, &reco, &truth, &response}) finish(*stream);

    // Written last: this marker distinguishes a completed preparation from
    // a fitted extraction and from interrupted preparation artifacts.
    auto meta = open("joint_input_metadata.txt");
    meta << "input_stage=prepared_joint_inputs_v1\nfit_objective=not_run\n";
    meta << "configured_kinematic=" << cfg.configured_kinematic
         << "\ninput_data_file=" << cfg.data_file << "\ninput_simc_file=" << cfg.simc_file
         << "\nvertex_source=" << cfg.vertex_simc_file
         << "\nvertex_epsilon_source=original_exclusive_h10_Q2i_Wi_hsxptari_hsyptari_SIMC_electron_angle"
         << "\nvertex_electron_hms_theta_deg=" << cfg.hms_theta_deg
         << "\nexclusive_selection=" << cfg.mmiss_select
         << "\nsim_reconstructed_mass_selection=" << (cfg.mmiss_select == "window" ?
             "configured_missing_mass_window" : "combined_data_geometry")
         << "\ncombined_mass_cut_file=" << cfg.mmiss_cut_file
         << "\nmmiss_lower_gev=" << cfg.mmiss_lower_gev << "\nmmiss_upper_gev=" << cfg.mmiss_upper_gev
         << "\ntarget_contam_factor=" << cfg.tgt_contam << "\ntarget_contam_factor_err=" << cfg.tgt_contam_err
         << "\ntarget_contam_usage=data_yield_divided_by_factor"
         << "\nresponse_absolute_scale=physical_SIMC_normalization_only_no_data_area_matching"
         << "\nsim_weight_mode=full_weight_over_sigcm\npartons_projection=no"
         << "\ndiamond_xb_q2_vertices=" << xsec_diamond_vertices_text(cfg);
    meta << "\ntprime_edges:";
    for (double v : tprime_edges) meta << ' ' << v;
    meta << "\nq2_edges:";
    for (double v : q2_edges) meta << ' ' << v;
    meta << "\nphi_edges:";
    for (double v : phi_edges) meta << ' ' << v;
    meta << "\nxb_edges_by_q2:";
    for (size_t iq=0; iq<xb_edges_by_q2.size(); ++iq) {
        meta << "\n  iq=" << iq << ':';
        for (double v : xb_edges_by_q2[iq]) meta << ' ' << v;
    }
    meta << '\n';
    finish(meta);
}
