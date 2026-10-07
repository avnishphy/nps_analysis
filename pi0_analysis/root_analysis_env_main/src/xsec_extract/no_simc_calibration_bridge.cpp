#include "xsec_no_simc_calibration.h"

#include <algorithm>
#include <cstring>

namespace {
void copy_error(const std::string& message, char* output, int capacity) {
    if (!output || capacity <= 0) return;
    const size_t count = std::min(message.size(), static_cast<size_t>(capacity - 1));
    std::memcpy(output, message.data(), count);
    output[count] = '\0';
}
}

extern "C" int nps_no_simc_calibration_fit(
    int nr, int np, const double* design, const double* subtracted_yield,
    const double* data_variance, const double* response_covariance,
    const double* fixed_feedin_mc_variance, const double* epsilon_max,
    int finite_mc, double rank_tolerance, int max_iterations, double fit_tolerance,
    double* parameters, double* covariance, double* final_variance,
    double* minima, double* boundary_tolerances, double* feasibility_tolerances,
    int* rank, double* condition, double* chi2, int* ndf,
    int* mc_iterations, int* positivity_iterations, int* boundary_active,
    char* error, int error_capacity) {
    try {
        if (nr <= 0 || np <= 0 || np % 3 != 0)
            throw std::runtime_error("invalid bridge dimensions");
        const int nb = np / 3;
        std::vector<std::vector<double>> matrix(nr, std::vector<double>(np));
        std::vector<double> y(subtracted_yield, subtracted_yield + nr);
        std::vector<double> variance(data_variance, data_variance + nr);
        std::vector<double> fixed(fixed_feedin_mc_variance,
                                  fixed_feedin_mc_variance + nr);
        std::vector<double> eps(epsilon_max, epsilon_max + nb);
        std::vector<std::vector<nps_xsec::ResponseCell>> cells(
            nr, std::vector<nps_xsec::ResponseCell>(nb));
        for (int row = 0; row < nr; ++row) {
            std::copy(design + static_cast<size_t>(row) * np,
                      design + static_cast<size_t>(row + 1) * np,
                      matrix[row].begin());
            for (int block = 0; block < nb; ++block)
                std::copy(response_covariance +
                              (static_cast<size_t>(row) * nb + block) * 9,
                          response_covariance +
                              (static_cast<size_t>(row) * nb + block + 1) * 9,
                          cells[row][block].covariance.begin());
        }
        const auto result = nps_xsec::fit_direct_calibration_problem(
            matrix, y, variance, cells, fixed, eps,
            finite_mc ? "finite-mc" : "data", rank_tolerance,
            max_iterations, fit_tolerance);
        std::copy(result.fit.parameters.begin(), result.fit.parameters.end(), parameters);
        std::copy(result.fit.covariance.begin(), result.fit.covariance.end(), covariance);
        std::copy(result.final_variance.begin(), result.final_variance.end(), final_variance);
        std::copy(result.minima.begin(), result.minima.end(), minima);
        std::copy(result.boundary_tolerances.begin(), result.boundary_tolerances.end(),
                  boundary_tolerances);
        std::copy(result.feasibility_tolerances.begin(), result.feasibility_tolerances.end(),
                  feasibility_tolerances);
        *rank = static_cast<int>(result.fit.rank);
        *condition = result.fit.condition;
        *chi2 = result.fit.chi2;
        *ndf = result.fit.ndf;
        *mc_iterations = result.mc_iterations;
        *positivity_iterations = result.positivity_iterations;
        *boundary_active = result.boundary_active ? 1 : 0;
        copy_error("", error, error_capacity);
        return 0;
    } catch (const std::exception& exception) {
        copy_error(exception.what(), error, error_capacity);
        return 1;
    }
}
