// Shared LT/TT with independent U across kinematic settings. The Python driver
// validates source presets and writes this program's compact, versioned input.
// Each setting contributes separate yield-per-mC rows and its own one-mC SIMC
// response. Stacking rows shares hadronic coefficients without summing or
// exposure-rescaling independent measurements.
#include "xsec_response.h"
#include "xsec_joint_positive_solver.h"

#include <TFile.h>
#include <TObjString.h>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

struct JointCell {
    std::array<double,3> basis{};
    std::array<double,9> covariance{};
};

static double mc_variance(const std::vector<JointCell>& cells,
                          const std::vector<int>& blocks,
                          const std::vector<std::array<int,3>>& columns,
                          const std::vector<double>& parameters) {
    long double total=0;
    for (size_t b=0; b<blocks.size(); ++b) {
        if (columns[b][0] < 0) continue;
        const auto& cell=cells[blocks[b]];
        long double sum=0, correction=0, absolute=0;
        for (int i=0; i<3; ++i) for (int j=0; j<3; ++j) {
            const long double term=static_cast<long double>(parameters[columns[b][i]])*
                cell.covariance[3*i+j]*parameters[columns[b][j]];
            absolute+=std::abs(term);
            const long double increment=term-correction;
            const long double next=sum+increment;
            correction=(next-sum)-increment;
            sum=next;
        }
        const long double roundoff=128.L*std::numeric_limits<double>::epsilon()*absolute;
        if (!std::isfinite(sum) || !std::isfinite(absolute) || sum < -roundoff)
            throw std::runtime_error("Joint xsec fit: invalid finite-MC covariance quadratic form");
        total+=std::max(0.L,sum);
    }
    const double result=static_cast<double>(total);
    if (!std::isfinite(result))
        throw std::runtime_error("Joint xsec fit: finite-MC variance overflow");
    return result;
}

struct JointRow {
    int setting = -1, local_row = -1;
    double data = 0, data_variance = 0;
    std::vector<JointCell> cells;
};

static void require(bool condition, const std::string& message) {
    if (!condition) throw std::runtime_error("Joint xsec fit: " + message);
}

static std::ofstream output(const fs::path& directory, const char* name) {
    std::ofstream stream(directory / name);
    require(static_cast<bool>(stream), std::string("cannot write ") + name);
    stream << std::setprecision(std::numeric_limits<double>::max_digits10);
    return stream;
}

int main(int argc, char** argv) {
    try {
        require(argc == 8, "usage: joint_solver problem.txt out_dir data|finite-mc positive_0_or_1 rank_tolerance max_iterations convergence_tolerance");
        const fs::path problem = argv[1], out_dir = argv[2];
        const std::string variance_mode = argv[3];
        require(variance_mode == "data" || variance_mode == "finite-mc", "invalid variance mode");
        const int positive_arg = std::stoi(argv[4]);
        require(positive_arg == 0 || positive_arg == 1, "positive flag must be 0 or 1");
        const bool positive = positive_arg == 1;
        const double rank_tolerance = std::stod(argv[5]);
        const int max_iterations = std::stoi(argv[6]);
        const double convergence_tolerance = std::stod(argv[7]);
        require(max_iterations > 0 && std::isfinite(convergence_tolerance) &&
                    convergence_tolerance > 0 && convergence_tolerance < 1,
                "invalid finite-MC iteration settings");

        std::ifstream input(problem);
        require(static_cast<bool>(input), "cannot open validated problem");
        std::string version;
        int nsettings = 0, nt = 0, nq = 0, nx = 0, nphi = 0, nr = 0, nb = 0;
        input >> version >> nsettings >> nt >> nq >> nx >> nphi >> nr >> nb;
        const int published = nt * nq * nx;
        require(version == "joint_xsec_v3" && nsettings >= 2 && nt > 0 &&
                    nq > 0 && nx > 0 && nphi >= 3 &&
                    nr == nsettings * published * nphi && nb == published + 6,
                "invalid problem dimensions");
        std::vector<std::vector<double>> epsilon_max(nsettings, std::vector<double>(nb));
        for (auto& setting : epsilon_max)
            for (double& value : setting) input >> value;
        require(static_cast<bool>(input), "cannot read epsilon maxima");

        std::vector<JointRow> rows(static_cast<size_t>(nr));
        for (int r = 0; r < nr; ++r) {
            auto& row = rows[r];
            input >> row.setting >> row.local_row >> row.data >> row.data_variance;
            require(static_cast<bool>(input) && row.setting == r / (published * nphi) &&
                        row.local_row == r % (published * nphi) &&
                        std::isfinite(row.data) && std::isfinite(row.data_variance) &&
                        row.data_variance >= 0,
                    "invalid or misordered reconstructed row " + std::to_string(r));
            row.cells.resize(static_cast<size_t>(nb));
            for (auto& cell : row.cells) {
                for (double& value : cell.basis) input >> value;
                for (double& value : cell.covariance) input >> value;
                require(static_cast<bool>(input), "truncated response at row " + std::to_string(r));
                for (double value : cell.basis)
                    require(std::isfinite(value), "nonfinite response basis");
                for (double value : cell.covariance)
                    require(std::isfinite(value), "nonfinite MC covariance");
            }
        }
        std::string trailing;
        require(!(input >> trailing), "unexpected trailing problem data");

        // A zero observed variance supplies no Gaussian fit weight. Keep the
        // row in diagnostics, but never invent a pseudocount.
        std::vector<int> fit_rows;
        int omitted_zero_variance = 0;
        for (int r = 0; r < nr; ++r) {
            const auto& row = rows[r];
            bool support = false;
            for (const auto& cell : row.cells)
                for (double value : cell.basis) support = support || value != 0.0;
            if (!support && (row.data != 0.0 || row.data_variance > 0.0))
                throw std::runtime_error("Joint xsec fit: data outside SIMC support at setting " +
                    std::to_string(row.setting) + " row " + std::to_string(row.local_row));
            if (support && row.data_variance == 0.0) ++omitted_zero_variance;
            if (support && row.data_variance > 0.0) fit_rows.push_back(r);
        }

        // Published blocks always retain their coefficients. A guard
        // enters only when it feeds a fitted row in at least one setting.
        std::vector<int> active_blocks;
        for (int b = 0; b < nb; ++b) {
            bool supported = b < published;
            if (!supported)
                for (int r : fit_rows)
                    for (double value : rows[r].cells[b].basis)
                        supported = supported || value != 0.0;
            if (supported) active_blocks.push_back(b);
        }
        // Each block has U_s for each setting and common LT/TT. A guard U_s
        // with no fitted-row support is omitted, never shared with another U.
        std::vector<std::vector<std::array<int,3>>> columns(nsettings,
            std::vector<std::array<int,3>>(active_blocks.size(), {{-1,-1,-1}}));
        struct Parameter { int block, setting, term; };
        std::vector<Parameter> labels;
        std::vector<std::array<int,3>> positive_columns;
        std::vector<std::pair<int,int>> positive_blocks;
        std::vector<double> active_epsilon;
        for (size_t j = 0; j < active_blocks.size(); ++j) {
            const int b = active_blocks[j];
            for (int s = 0; s < nsettings; ++s) {
                bool supported = b < published;
                for (int r : fit_rows)
                    if (rows[r].setting == s)
                        for (double value : rows[r].cells[b].basis)
                            supported = supported || value != 0.0;
                if (supported) {
                    columns[s][j][0] = static_cast<int>(labels.size());
                    labels.push_back({b, s, 0});
                }
            }
            for (int term = 1; term < 3; ++term) {
                for (int s = 0; s < nsettings; ++s)
                    columns[s][j][term] = static_cast<int>(labels.size());
                labels.push_back({b, -1, term});
            }
            for (int s = 0; s < nsettings; ++s) {
                if (columns[s][j][0] < 0) continue;
                const double e = epsilon_max[s][b];
                require(std::isfinite(e) && e > 0.0 && e < 1.0,
                        "active setting/truth block has no valid event epsilon maximum: " +
                        std::to_string(s) + "/" + std::to_string(b));
                active_epsilon.push_back(e);
                positive_columns.push_back(columns[s][j]);
                positive_blocks.emplace_back(s, b);
            }
        }
        const size_t npar = labels.size();
        require(fit_rows.size() > npar, "too few independent rows for per-setting U and shared LT/TT");
        const auto design_row = [&](const JointRow& row) {
            std::vector<double> design(npar);
            for (size_t j = 0; j < active_blocks.size(); ++j) {
                const auto& index = columns[row.setting][j];
                if (index[0] < 0) continue;
                for (int term = 0; term < 3; ++term)
                    design[index[term]] = row.cells[active_blocks[j]].basis[term];
            }
            return design;
        };
        std::vector<std::vector<double>> design;
        std::vector<double> data, data_variance;
        for (int r : fit_rows) {
            design.push_back(design_row(rows[r]));
            data.push_back(rows[r].data);
            data_variance.push_back(rows[r].data_variance);
        }

        nps_xsec::PositiveSolution positive_state;
        auto solve = [&](const std::vector<double>& variance) {
            if (!positive)
                return nps_xsec::solve_weighted_response(design, data, variance, rank_tolerance);
            positive_state = nps_xsec::solve_joint_positive_response(
                design, data, variance, active_epsilon, positive_columns, rank_tolerance);
            return positive_state.fit;
        };
        std::vector<double> fit_variance = data_variance;
        auto fit = solve(fit_variance);
        bool converged = variance_mode == "data";
        int iterations = 0;
        for (; variance_mode == "finite-mc" && iterations < max_iterations;) {
            ++iterations;
            std::vector<double> next_variance = data_variance;
            for (size_t i = 0; i < fit_rows.size(); ++i)
                next_variance[i] += mc_variance(
                    rows[fit_rows[i]].cells, active_blocks,
                    columns[rows[fit_rows[i]].setting], fit.parameters);
            auto next = solve(next_variance);
            double parameter_change = 0.0, variance_change = 0.0;
            for (size_t j = 0; j < npar; ++j) {
                const double scale = std::max(
                    std::abs(next.parameters[j]), std::sqrt(next.covariance[j * npar + j]));
                parameter_change = std::max(parameter_change,
                    std::abs(next.parameters[j] - fit.parameters[j]) / scale);
            }
            for (size_t i = 0; i < fit_rows.size(); ++i)
                variance_change = std::max(variance_change,
                    std::abs(next_variance[i] - fit_variance[i]) / next_variance[i]);
            fit = std::move(next);
            fit_variance = std::move(next_variance);
            if (std::max(parameter_change, variance_change) < convergence_tolerance) {
                converged = true;
                break;
            }
        }
        require(converged, "finite-MC variance iteration did not converge");
        require(fit.rank == npar, "shared response is not full rank");

        // A constrained-boundary curvature is only a diagnostic. It is not a
        // Gaussian covariance or a confidence interval for the estimate.
        const bool boundary = positive && positive_state.boundary_active;
        const double nan = std::numeric_limits<double>::quiet_NaN();
        auto parameters = output(out_dir, "joint_parameters.csv");
        parameters << "parameter_index,truth_block,region,it,iq,ix,component,value,error_stat_plus_mc,epsilon_max,setting_index\n";
        for (size_t col = 0; col < npar; ++col) {
            const auto& label = labels[col];
            const int b = label.block;
            const bool pub = b < published;
            const int it = pub ? b / (nq * nx) : -1;
            const int iq = pub ? (b / nx) % nq : -1;
            const int ix = pub ? b % nx : -1;
            const std::string region = pub ? "published" : nps_xsec::guard_name(b - published);
            const double err = boundary ? nan : std::sqrt(fit.covariance[col * npar + col]);
            double e = 0;
            for (int s = 0; s < nsettings; ++s)
                if (label.setting < 0 || label.setting == s) e = std::max(e, epsilon_max[s][b]);
            parameters << col << ',' << b << ',' << region << ',' << it << ','
                       << iq << ',' << ix << ',' << (label.term == 0 ? "U" : label.term == 1 ? "LT" : "TT")
                       << ',' << fit.parameters[col] << ',' << err << ',' << e
                       << ',' << label.setting << '\n';
        }

        auto covariance = output(out_dir, "joint_covariance.csv");
        covariance << "parameter_i,parameter_j,stat_plus_mc_covariance,unconstrained_curvature_inverse\n";
        for (size_t i = 0; i < npar; ++i)
            for (size_t j = 0; j < npar; ++j)
                covariance << i << ',' << j << ','
                           << (boundary ? nan : fit.covariance[i * npar + j]) << ','
                           << fit.covariance[i * npar + j] << '\n';

        auto used_rows = output(out_dir, "joint_rows.csv");
        used_rows << "setting_index,reco_row,fit_index,data,data_variance,variance_used,mc_variance_final,prediction,residual,pull\n";
        std::vector<double> chi2_by_setting(static_cast<size_t>(nsettings), 0.0);
        int fit_index = 0;
        for (int r = 0; r < nr; ++r) {
            const auto& row = rows[r];
            const auto x = design_row(row);
            const double prediction = std::inner_product(x.begin(), x.end(), fit.parameters.begin(), 0.0);
            const double mc_final =
                mc_variance(row.cells, active_blocks, columns[row.setting], fit.parameters);
            const bool included = fit_index < static_cast<int>(fit_rows.size()) && fit_rows[fit_index] == r;
            const double variance = included ? fit_variance[fit_index] : nan;
            const double residual = row.data - prediction;
            const double pull = included ? residual / std::sqrt(variance) : nan;
            if (included) chi2_by_setting[row.setting] += pull * pull;
            used_rows << row.setting << ',' << row.local_row << ','
                      << (included ? fit_index : -1) << ',' << row.data << ','
                      << row.data_variance << ',' << variance << ',' << mc_final << ','
                      << prediction << ',' << residual << ',' << pull << '\n';
            if (included) ++fit_index;
        }
        auto setting_chi2 = output(out_dir, "joint_setting_chi2.csv");
        setting_chi2 << "setting_index,chi2_contribution\n";
        for (int s = 0; s < nsettings; ++s)
            setting_chi2 << s << ',' << chi2_by_setting[s] << '\n';

        auto positivity = output(out_dir, "joint_positivity.csv");
        positivity << "truth_block,epsilon_max,minimum_response_bracket,boundary_tolerance,feasibility_tolerance,setting_index\n";
        for (size_t j = 0; j < positive_columns.size(); ++j) {
            const auto& index = positive_columns[j];
            const double minimum = static_cast<double>(nps_xsec::positive_detail::minimum(
                fit.parameters[index[0]], fit.parameters[index[1]], fit.parameters[index[2]],
                active_epsilon[j]).first);
            positivity << positive_blocks[j].second << ',' << active_epsilon[j] << ',' << minimum << ','
                       << (positive ? positive_state.boundary_tolerances[j] : nan) << ','
                       << (positive ? positive_state.feasibility_tolerances[j] : nan)
                       << ',' << positive_blocks[j].first << '\n';
        }
        auto summary = output(out_dir, "joint_summary.txt");
        summary << "fit_components=U_per_setting,LT_shared,TT_shared\n"
                << "positivity_epsilon_range=zero_to_observed_maximum_per_setting_and_truth_block\n"
                << "settings=" << nsettings << '\n'
                << "rows_used=" << fit_rows.size() << '\n'
                << "rows_available=" << nr << '\n'
                << "omitted_supported_zero_variance_rows=" << omitted_zero_variance << '\n'
                << "active_truth_blocks=" << active_blocks.size() << '\n'
                << "parameters=" << npar << '\n'
                << "rank=" << fit.rank << '\n'
                << "condition=" << fit.condition << '\n'
                << "chi2=" << fit.chi2 << '\n'
                << "ndf=" << fit.ndf << '\n'
                << "variance_mode=" << variance_mode << '\n'
                << "mc_iterations=" << iterations << '\n'
                << "positive_xsec=" << (positive ? 1 : 0) << '\n'
                << "positivity_boundary_active=" << (boundary ? 1 : 0) << '\n'
                << "covariance_status=" << (boundary ? "unavailable_boundary_constrained_estimate" :
                                            "conditional_inverse_information") << '\n'
                << "input_normalization=each_setting_yield_per_mC_and_one_mC_SIMC_response_stacked_as_separate_rows\n"
                << "target_factor_uncertainty=not_propagated\n";
        // Mirror the numerical result in ROOT for downstream analysis. The
        // manifest and CSVs retain block labels, setting paths, and row maps.
        TFile root((out_dir / "joint_xsec_output.root").c_str(), "RECREATE");
        require(!root.IsZombie(), "cannot create joint ROOT output");
        TVectorD root_parameters(static_cast<int>(npar));
        TVectorD root_singular(static_cast<int>(fit.singular_values.size()));
        TMatrixD root_covariance(static_cast<int>(npar), static_cast<int>(npar));
        TMatrixD root_curvature(static_cast<int>(npar), static_cast<int>(npar));
        for (size_t i=0; i<npar; ++i) {
            root_parameters[static_cast<int>(i)] = fit.parameters[i];
            for (size_t j=0; j<npar; ++j) {
                const double value=fit.covariance[i*npar+j];
                root_covariance(static_cast<int>(i),static_cast<int>(j))=boundary?nan:value;
                root_curvature(static_cast<int>(i),static_cast<int>(j))=value;
            }
        }
        for (size_t i=0; i<fit.singular_values.size(); ++i)
            root_singular[static_cast<int>(i)]=fit.singular_values[i];
        root_parameters.Write("joint_parameters");
        root_covariance.Write("joint_covariance_stat_plus_mc");
        root_curvature.Write("joint_unconstrained_curvature_inverse");
        root_singular.Write("joint_scaled_singular_values");
        std::ostringstream root_meta;
        root_meta << "method=shared_LT_TT_independent_U_forward_response_fit\n"
                  << "settings=" << nsettings << "\nrows_used=" << fit_rows.size()
                  << "\nrank=" << fit.rank << "\nchi2=" << fit.chi2
                  << "\nndf=" << fit.ndf << "\nvariance_mode=" << variance_mode
                  << "\npositive_xsec=" << positive
                  << "\npositivity_boundary_active=" << boundary
                  << "\npositivity_epsilon_range=zero_to_observed_maximum_per_setting_and_truth_block\n"
                  << "target_factor_uncertainty=not_propagated";
        TObjString(root_meta.str().c_str()).Write("analysis_metadata");
        root.Write();
        require(!root.TestBit(TFile::kWriteError), "failed to write joint ROOT output");
        root.Close();
        require(static_cast<bool>(parameters) && static_cast<bool>(covariance) &&
                    static_cast<bool>(used_rows) && static_cast<bool>(setting_chi2) &&
                    static_cast<bool>(positivity) && static_cast<bool>(summary),
                "failed while writing joint outputs");
        std::cout << "Joint fit complete: settings=" << nsettings
                  << " rows=" << fit_rows.size() << " parameters=" << npar
                  << " rank=" << fit.rank << " chi2/ndf=" << fit.chi2 << '/' << fit.ndf
                  << " boundary=" << boundary << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FATAL] " << error.what() << '\n';
        return 1;
    }
}
