#pragma once

// Build and freeze shared generated/reconstructed bin edges before forming the response.
#include "xsec_analysis.h"

inline void ExclPi0XSecAnalysis::build_binning() {
    // The two methods read identical fixed Q2, xB, phi and tprime edges.
    // The ratio method uses its separately configured physical-t edges.
    phi_edges = cfg.phi_bin_edges;
    tprime_edges = cfg.tprime_bin_edges;
    q2_edges = cfg.q2_bin_edges;
    xb_edges_by_q2 = cfg.xb_bin_edges_by_q2;
    // Retain the legacy global xB metadata field. Per-Q2 rows are authoritative.
    xb_edges = xb_edges_by_q2.front();

    auto log_edges = [&](const char* name, const std::vector<double>& edges) {
        std::ostringstream out;
        out << name << " edges:";
        for (double edge : edges) out << " " << edge;
        log(out.str());
    };
    log_edges("phi", phi_edges);
    log_edges("tprime", tprime_edges);
    log_edges("Q2", q2_edges);
    for (int iq = 0; iq < cfg.n_q2; ++iq) {
        std::string name = "xB Q2 bin " + std::to_string(iq);
        log_edges(name.c_str(), xb_edges_by_q2[iq]);
    }
}
