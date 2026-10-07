#pragma once

// Exact no-SIMC-model estimator adapter for local calibration toys.  The
// production response solver and finite-MC covariance contraction are reused
// directly; this header owns only the same outer feasible-GLS iteration used
// by xsec_fit.h.  Campaign code must prove central numerical parity before
// accepting any toys.
#include "xsec_positive_solver.h"
#include "xsec_response.h"
#include <numeric>

namespace nps_xsec {

struct DirectCalibrationFit {
    LinearSolution fit;
    std::vector<double> final_variance;
    std::vector<double> minima;
    std::vector<double> boundary_tolerances;
    std::vector<double> feasibility_tolerances;
    int mc_iterations = 0;
    int positivity_iterations = 0;
    bool mc_converged = false;
    bool boundary_active = false;
};

inline DirectCalibrationFit fit_direct_calibration_problem(
    const std::vector<std::vector<double>>& design,
    const std::vector<double>& subtracted_yield,
    const std::vector<double>& data_variance,
    const std::vector<std::vector<ResponseCell>>& response_cells,
    const std::vector<double>& fixed_feedin_mc_variance,
    const std::vector<double>& epsilon_max,
    const std::string& variance_mode,
    double rank_tolerance,
    int max_iterations,
    double fit_tolerance) {
    if (design.size() != response_cells.size() || design.size() != fixed_feedin_mc_variance.size())
        throw std::runtime_error("Direct calibration fit: response row-count mismatch");
    if (design.empty() || epsilon_max.empty() || design.front().size() != 3 * epsilon_max.size())
        throw std::runtime_error("Direct calibration fit: parameter/block dimension mismatch");
    if (variance_mode != "data" && variance_mode != "finite-mc")
        throw std::runtime_error("Direct calibration fit: variance mode must be data or finite-mc");
    if (max_iterations < 1 || !(fit_tolerance > 0. && fit_tolerance < 1.))
        throw std::runtime_error("Direct calibration fit: invalid finite-MC controls");

    std::vector<int> blocks(epsilon_max.size());
    std::iota(blocks.begin(), blocks.end(), 0);
    DirectCalibrationFit output;
    auto solve = [&](const std::vector<double>& variance) {
        auto result = solve_positive_response(design, subtracted_yield, variance,
                                              epsilon_max, rank_tolerance);
        output.minima = result.minima;
        output.boundary_tolerances = result.boundary_tolerances;
        output.feasibility_tolerances = result.feasibility_tolerances;
        output.positivity_iterations = result.iterations;
        output.boundary_active = result.boundary_active;
        return result.fit;
    };

    output.final_variance = data_variance;
    output.fit = solve(output.final_variance);
    output.mc_converged = variance_mode == "data";
    for (; variance_mode == "finite-mc" && output.mc_iterations < max_iterations;) {
        ++output.mc_iterations;
        std::vector<double> next_variance = data_variance;
        for (size_t row = 0; row < design.size(); ++row)
            next_variance[row] += mc_prediction_variance(
                response_cells[row], blocks, output.fit.parameters) +
                fixed_feedin_mc_variance[row];
        auto next = solve(next_variance);
        double parameter_change = 0., variance_change = 0.;
        const size_t np = next.parameters.size();
        for (size_t column = 0; column < np; ++column) {
            const double scale = std::max(std::abs(next.parameters[column]),
                std::sqrt(next.covariance[column * np + column]));
            parameter_change = std::max(parameter_change,
                std::abs(next.parameters[column] - output.fit.parameters[column]) / scale);
        }
        for (size_t row = 0; row < next_variance.size(); ++row)
            variance_change = std::max(variance_change,
                std::abs(next_variance[row] - output.final_variance[row]) / next_variance[row]);
        output.fit = std::move(next);
        output.final_variance = std::move(next_variance);
        if (std::max(parameter_change, variance_change) < fit_tolerance) {
            output.mc_converged = true;
            break;
        }
    }
    if (!output.mc_converged)
        throw std::runtime_error("Direct calibration fit: finite-MC iteration did not converge");
    return output;
}

} // namespace nps_xsec
