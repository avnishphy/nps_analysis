#pragma once

// Exact combined-data 2D mass geometry used for the smeared SIMC selection.
// The data decision itself comes from the stored combined event flag.
#include "xsec_config.h"

struct XsecMassGeometry {
    double mean_x = 0, mean_y = 0;
    double cov_xx = 0, cov_xy = 0, cov_yy = 0;
    double d2_cut = 0;
    double x_min = 0, x_max = 0, y_min = 0, y_max = 0;

    bool contains(double x, double y) const {
        if (!std::isfinite(x) || !std::isfinite(y) ||
            x < x_min || x >= x_max || y < y_min || y >= y_max) return false;
        const double det = cov_xx * cov_yy - cov_xy * cov_xy;
        const double dx = x - mean_x, dy = y - mean_y;
        return (cov_yy * dx * dx - 2.0 * cov_xy * dx * dy +
                cov_xx * dy * dy) / det <= d2_cut;
    }
};

inline XsecMassGeometry read_xsec_mass_geometry(const std::string& path,
                                                 const std::string& mode) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot read combined mass-cut metadata: " + path);
    std::map<std::string, std::string> values;
    std::string line;
    while (std::getline(in, line)) {
        const auto pos = line.find('=');
        if (pos != std::string::npos)
            values[line.substr(0, pos)] = line.substr(pos + 1);
    }
    auto number = [&](const std::string& key) {
        const auto it = values.find(key);
        if (it == values.end())
            throw std::runtime_error("Combined mass-cut metadata lacks " + key +
                " in " + path + "; regenerate it with combine_analysis_branches.py.");
        size_t used = 0;
        const double value = std::stod(it->second, &used);
        if (used != it->second.size() || !std::isfinite(value))
            throw std::runtime_error("Invalid combined mass-cut metadata " + key + " in " + path);
        return value;
    };
    if (values["tag"] != "combined_2d_mass_cut")
        throw std::runtime_error("Unexpected combined mass-cut metadata tag in " + path);
    const std::string prefix = mode == "mcd" ? "mcd_" : "";
    if (number(prefix == "mcd_" ? "mcd_valid" : "ellipse_valid") != 1.0)
        throw std::runtime_error("Requested " + mode + " mass cut is invalid in " + path);
    XsecMassGeometry g;
    g.mean_x = number(prefix + "mean_mpi0");
    g.mean_y = number(prefix + "mean_mmiss");
    g.cov_xx = number(prefix + "cov_mpi0_mpi0");
    g.cov_xy = number(prefix + "cov_mpi0_mmiss");
    g.cov_yy = number(prefix + "cov_mmiss_mmiss");
    g.d2_cut = number(prefix == "mcd_" ? "mcd_d2_cut" : "ellipse_d2_cut");
    g.x_min = number("mpi0_min");
    g.x_max = number("mpi0_max");
    g.y_min = number("mmiss_min");
    g.y_max = number("mmiss_max");
    const double det = g.cov_xx * g.cov_yy - g.cov_xy * g.cov_xy;
    if (!(g.cov_xx > 0 && g.cov_yy > 0 && det > 0 && g.d2_cut > 0 &&
          g.x_min < g.x_max && g.y_min < g.y_max))
        throw std::runtime_error("Invalid combined " + mode + " mass-cut geometry in " + path);
    return g;
}
