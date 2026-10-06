#pragma once

// Event cache for independent reconstructed/truth/publication binning.
// Existing selection, run normalization and raw-SIMC matching are reused.
#include "xsec_analysis.h"

inline std::string forward_json_string(const std::string& value) {
    std::ostringstream out;
    out << '"';
    for (unsigned char c : value) {
        if (c == '"' || c == '\\') out << '\\' << c;
        else if (c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
        else out << c;
    }
    out << '"';
    return out.str();
}

inline void ExclPi0XSecAnalysis::begin_forward_inputs() {
    for (const char* name : {"data_events.csv", "mc_events.csv", "forward_cache_manifest.json",
                             "data_events.csv.partial", "mc_events.csv.partial",
                             "forward_cache_manifest.json.partial"})
        if (fs::exists(fs::path(cfg.out_dir)/name))
            die("Forward event export refuses existing/partial cache: " + (fs::path(cfg.out_dir)/name).string());
    forward_data_stream.open((fs::path(cfg.out_dir)/"data_events.csv.partial").string());
    forward_mc_stream.open((fs::path(cfg.out_dir)/"mc_events.csv.partial").string());
    if (!forward_data_stream || !forward_mc_stream) die("Cannot open forward event cache");
    forward_data_stream << std::setprecision(17)
        << "event_id,run_number,q2,xb,tprime,phi,weight\n";
    forward_mc_stream << std::setprecision(17)
        << "event_id,reco_q2,reco_xb,reco_tprime,reco_phi,truth_q2,truth_xb,truth_tprime,truth_phi,epsilon,base_weight,nominal_weight\n";
}

inline void ExclPi0XSecAnalysis::finish_forward_inputs() {
    forward_data_stream.flush(); forward_mc_stream.flush();
    if (!forward_data_stream || !forward_mc_stream) die("Failed writing forward event cache");
    forward_data_stream.close(); forward_mc_stream.close();
    if (!forward_data_count || !forward_mc_count) die("Empty selected data or MC event cache");
    std::ofstream out((fs::path(cfg.out_dir)/"forward_cache_manifest.json.partial").string());
    out << std::setprecision(17) << "{\n  \"schema_version\": 2,\n  \"complete\": true,\n"
        << "  \"data_file\": " << forward_json_string(cfg.data_file) << ",\n"
        << "  \"simc_file\": " << forward_json_string(cfg.simc_file) << ",\n"
        << "  \"vertex_simc_file\": " << forward_json_string(cfg.vertex_simc_file) << ",\n"
        << "  \"data_root_uuid\": " << forward_json_string(f_data->GetUUID().AsString()) << ",\n"
        << "  \"simc_root_uuid\": " << forward_json_string(f_sim->GetUUID().AsString()) << ",\n"
        << "  \"vertex_root_uuid\": " << forward_json_string(f_vertex->GetUUID().AsString()) << ",\n"
        << "  \"data_events\": " << forward_data_count << ",\n  \"mc_events\": " << forward_mc_count << ",\n"
        << "  \"target_divisor\": " << cfg.tgt_contam << ",\n  \"target_divisor_error\": " << cfg.tgt_contam_err << ",\n"
        << "  \"data_weight_units\": \"yield_per_mC_after_target_division\",\n"
        << "  \"mc_base_definition\": \"full_weight/sigcm; original generated normalization retained\",\n"
        << "  \"mc_nominal_definition\": \"full_weight; fixed Q2/xB feed-in row contribution\",\n"
        << "  \"rectangular_kinematic_cuts_applied\": false,\n"
        << "  \"mass_selection\": " << forward_json_string(cfg.mmiss_select) << ",\n"
        << "  \"mass_cut_file\": " << forward_json_string(cfg.mmiss_cut_file) << ",\n"
        << "  \"mmiss_lower_gev\": " << cfg.mmiss_lower_gev << ",\n  \"mmiss_upper_gev\": " << cfg.mmiss_upper_gev << ",\n"
        << "  \"hms_theta_deg\": " << cfg.hms_theta_deg << ",\n  \"mp_gev\": " << cfg.mp << ",\n  \"mpi0_gev\": " << cfg.mpi0 << ",\n"
        << "  \"diamond_xb_q2_vertices\": [";
    for (size_t i=0;i<cfg.diamond_xb_q2_vertices.size();++i) {
        if (i) out << ',';
        out << '[' << cfg.diamond_xb_q2_vertices[i][0] << ',' << cfg.diamond_xb_q2_vertices[i][1] << ']';
    }
    out << "],\n  \"upstream_weight_uncertainty_included\": false,\n"
        << "  \"mc_sampling\": \"accepted records with original event_id; failed trials encoded in full_weight normalization\"\n}\n";
    out.flush();
    if (!out) die("Failed writing forward event manifest");
    out.close();
    fs::rename(fs::path(cfg.out_dir)/"data_events.csv.partial", fs::path(cfg.out_dir)/"data_events.csv");
    fs::rename(fs::path(cfg.out_dir)/"mc_events.csv.partial", fs::path(cfg.out_dir)/"mc_events.csv");
    // The completion manifest is committed last; partial exports are not usable.
    fs::rename(fs::path(cfg.out_dir)/"forward_cache_manifest.json.partial", fs::path(cfg.out_dir)/"forward_cache_manifest.json");
}
