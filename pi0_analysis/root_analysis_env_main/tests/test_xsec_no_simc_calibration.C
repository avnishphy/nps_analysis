#include "../src/xsec_extract/xsec_no_simc_calibration.h"

#include <iostream>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
bool close(double first, double second, double tolerance = 2e-10) {
    return std::abs(first - second) <= tolerance * std::max(1., std::abs(second));
}
std::vector<double> predict(const std::vector<std::vector<double>>& design,
                            const std::vector<double>& parameters) {
    std::vector<double> result(design.size());
    for (size_t row = 0; row < design.size(); ++row)
        for (size_t column = 0; column < parameters.size(); ++column)
            result[row] += design[row][column] * parameters[column];
    return result;
}
}

int main() {
    std::vector<std::vector<double>> design(12, std::vector<double>(6));
    for (size_t row = 0; row < design.size(); ++row) design[row][row % 6] = 1.;
    const std::vector<double> truth{3., -.4, .2, 7., .8, -.3};
    const auto yield = predict(design, truth);
    const std::vector<double> variance(design.size(), 1.);
    const std::vector<double> fixed(design.size(), 0.);
    std::vector<std::vector<nps_xsec::ResponseCell>> cells(
        design.size(), std::vector<nps_xsec::ResponseCell>(2));
    const auto data = nps_xsec::fit_direct_calibration_problem(
        design, yield, variance, cells, fixed, {.5, .7}, "data", 1e-10, 100, 1e-6);
    require(data.fit.rank == 6 && data.mc_converged && data.mc_iterations == 0,
            "data-variance direct fit diagnostics failed");
    for (size_t index = 0; index < truth.size(); ++index)
        require(close(data.fit.parameters[index], truth[index]), "central parity failed");

    for (size_t row = 0; row < design.size(); ++row)
        for (size_t block = 0; block < 2; ++block)
            for (int component = 0; component < 3; ++component)
                cells[row][block].covariance[3 * component + component] = 1e-4;
    const auto finite = nps_xsec::fit_direct_calibration_problem(
        design, yield, variance, cells, fixed, {.5, .7}, "finite-mc", 1e-10, 100, 1e-8);
    require(finite.mc_converged && finite.mc_iterations > 0,
            "finite-MC iteration did not run to convergence");
    for (size_t index = 0; index < truth.size(); ++index)
        require(close(finite.fit.parameters[index], truth[index]), "finite-MC parity failed");

    std::vector<std::vector<double>> boundary_design(6, std::vector<double>(3));
    for (size_t row = 0; row < boundary_design.size(); ++row)
        boundary_design[row][row % 3] = 1.;
    std::vector<std::vector<nps_xsec::ResponseCell>> boundary_cells(
        boundary_design.size(), std::vector<nps_xsec::ResponseCell>(1));
    const auto boundary = nps_xsec::fit_direct_calibration_problem(
        boundary_design, predict(boundary_design, {-1., 0., 0.}),
        std::vector<double>(6, 1.), boundary_cells, std::vector<double>(6, 0.),
        {.5}, "data", 1e-10, 20, 1e-6);
    require(boundary.boundary_active && close(boundary.fit.parameters[0], 0.),
            "continuous positivity boundary was not enforced");

    bool singular_rejected = false;
    for (auto& row : boundary_design) row[1] = row[0];
    try {
        nps_xsec::fit_direct_calibration_problem(
            boundary_design, std::vector<double>(6, 0.), std::vector<double>(6, 1.),
            boundary_cells, std::vector<double>(6, 0.), {.5}, "data", 1e-10, 20, 1e-6);
    } catch (const std::runtime_error&) { singular_rejected = true; }
    require(singular_rejected, "singular calibration toy was not rejected");
    std::cout << "PASS no-SIMC calibration adapter parity, iteration, boundary, and singular rejection\n";
}
