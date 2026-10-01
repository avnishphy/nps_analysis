#pragma once

// Read-only input validation and generated-event linkage. Invalid physics inputs must fail explicitly.
#include "xsec_analysis.h"

inline void ExclPi0XSecAnalysis::load_input() {
    // Read-only files: normalization comes from their production contracts.
    f_sim = TFile::Open(cfg.simc_file.c_str(), "READ");
    if (!f_sim || f_sim->IsZombie()) die("Cannot open SIMC input file.");

    f_data = TFile::Open(cfg.data_file.c_str(), "READ");
    if (!f_data || f_data->IsZombie()) die("Cannot open data input file.");

    t_sim = dynamic_cast<TTree*>(f_sim->Get(cfg.simc_tree.c_str()));
    t_data = dynamic_cast<TTree*>(f_data->Get(cfg.data_tree.c_str()));
    if (!t_sim) die("Cannot find SIMC tree.");
    if (!t_data) die("Cannot find data tree.");
    for(const char* name:{"Q2","t","tmin","xB","phi","pi0_weight","scale","charge_uC","run_number","mmiss_all","W"})
        if(!t_data->GetBranch(name)) die("Data missing required branch: "+std::string(name));
    for(const char* name:{"Q2","t","tmin","xB","phi","full_weight","is_exclusive","mmiss","W","sigcm"})
        if(!t_sim->GetBranch(name)) die("SIMC missing required branch: "+std::string(name));
    if (cfg.mmiss_select != "window") {
        if (!t_data->GetBranch("mpi0_all") || !t_data->GetBranch("mmiss_all_corr") ||
            !t_sim->GetBranch("mpi0"))
            die("Combined mass selection diagnostics require data mpi0_all/mmiss_all_corr and SIMC mpi0 branches.");
        const char* flag = cfg.mmiss_select == "mcd" ?
            "is_exclusive_mcd_combined" : "is_exclusive_ellipse_combined";
        if (!t_data->GetBranch(flag)) die("Data missing selected mass-cut flag: " + std::string(flag));
    }
    if (cfg.vertex_simc_file.empty())
        die("Vertex epsilon requires the matching original exclusive h10; provide --vertex_simc_file.");
    f_vertex = TFile::Open(cfg.vertex_simc_file.c_str(), "READ");
    if (!f_vertex || f_vertex->IsZombie()) die("Cannot open --vertex_simc_file: " + cfg.vertex_simc_file);
    t_vertex = dynamic_cast<TTree*>(f_vertex->Get("h10"));
    if (!t_vertex) die("--vertex_simc_file must contain an original SIMC h10 tree.");
    for (const char* branch : {"Q2i", "Wi", "ti", "phipqi", "sigcm", "hsxptari", "hsyptari"})
        if (!t_vertex->GetBranch(branch)) die("--vertex_simc_file missing h10 branch: " + std::string(branch));
    if (!t_sim->GetBranch("event_id")) die("--vertex_simc_file requires event_id in the smeared tree.");
    vertex_from_raw_simc = true;
}

inline void ExclPi0XSecAnalysis::load_mass_cut() {
    if (cfg.mmiss_select == "mcd" || cfg.mmiss_select == "ellipse") {
        mass_geometry = read_xsec_mass_geometry(cfg.mmiss_cut_file, cfg.mmiss_select);
        log("Loaded " + cfg.mmiss_select + " geometry from " + cfg.mmiss_cut_file);
    }
}

inline void ExclPi0XSecAnalysis::detect_optional_branches() {
    has_helicity = (t_data->GetBranch("helicity") != nullptr);
    has_sim_helicity = (t_sim->GetBranch("helicity") != nullptr);

    for (const auto& b : cfg.model_xsec_candidates) {
        if (t_sim->GetBranch(b.c_str())) {
            has_model_xsec = true;
            model_xsec_branch = b;
            break;
        }
    }

    if (!has_model_xsec) {
        die("No SIMC CM model cross-section branch 'sigcm' found. This extraction requires de-modeling with full_weight/sigcm.");
    }
    // The SIMC channel flag is generated truth. Reconstructed mass selection
    // is the window or the measured combined-data 2D geometry.
    if (!t_data->GetBranch("mmiss_all") || !t_sim->GetBranch("mmiss") ||
        !t_sim->GetBranch("is_exclusive"))
        die("Extraction requires data mmiss_all and SIMC mmiss/is_exclusive branches.");
    if (!t_data->GetBranch("mmiss_all") || !t_sim->GetBranch("mmiss"))
        die("Missing-mass diagnostics require data mmiss_all and SIMC mmiss branches.");
    has_vertex_kinematics = vertex_from_raw_simc;
    log("Raw exclusive h10 supplies generated coordinates and vertex epsilon; selected event_id and sigcm must match.");
    if (!has_helicity) warn("No helicity branch found in data tree. TL' will not be optimized.");
    if (has_helicity && !has_sim_helicity) {
        warn("Data helicity exists but SIMC helicity branch is missing. TL' optimization is disabled.");
    }
}
