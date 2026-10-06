// Reuse the reference steering and all its event/response/output machinery.
// Only the entry point selects the model mapping inside the shared fit layer.
// Keep the reference source byte-for-byte unchanged; its main is unused here.
#define main nps_unused_reference_main
#include "excl_xsec_pi0_analysis_no_simc_model.C"
#undef main

namespace {

bool nonempty_environment(const char* name) {
    const char* value = std::getenv(name);
    return value && value[0];
}

bool shared_option_takes_value(const std::string& option) {
    static const std::set<std::string> options = {
        "--kin", "--target", "--output-base", "--root-dir", "--data-file",
        "--sim-file", "--vertex_simc_file", "--vertex-simc-file", "--out-dir",
        "--out-root", "--out-csv", "--out-slice-csv", "--all-plots-pdf",
        "--mmiss-lower", "--mmiss-upper", "--mmiss_select", "--mmiss-select",
        "--mmiss-cut-file", "--target-contam", "--target-contam-err", "--ebeam",
        "--partons-warmups", "--partons-calls", "--svd-rank-tolerance",
        "--mc-max-iterations", "--mc-fit-tolerance", "--fit-variance",
        "--fit-objective", "--scaled-empty-scale", "--joint-plot-input"
    };
    return options.count(option) != 0;
}

std::string model_output_name(std::string name) {
    const auto position = name.find("no_simc_model");
    if (position != std::string::npos)
        name.replace(position, 13, "simc_model");
    return name;
}

} // namespace

int main(int argc, char** argv) {
    try {
        nps_xsec::ProxyOptions options;
        std::vector<char*> args{argv[0]};
        bool out_dir_explicit = nonempty_environment("NPS_XSEC_OUT_DIR");
        std::array<bool, 4> output_explicit = {
            nonempty_environment("NPS_XSEC_OUT_ROOT"),
            nonempty_environment("NPS_XSEC_OUT_CSV"),
            nonempty_environment("NPS_XSEC_OUT_SLICE_CSV"),
            nonempty_environment("NPS_XSEC_ALL_PLOTS_PDF")
        };
        size_t positional_count = 0;

        for (int i = 1; i < argc; ++i) {
            const std::string key = argv[i];
            auto next = [&]() {
                if (++i >= argc) throw std::runtime_error("Missing value for " + key);
                return std::string(argv[i]);
            };

            if (key == "--model") {
                throw std::runtime_error(
                    "Physics-model selection removed; use --mode simc_model in the pipeline");
            } else if (key == "--fit-strategy") {
                options.fit_strategy = next();
                if (options.fit_strategy != "joint_minuit" &&
                    options.fit_strategy != "staged_feasible")
                    throw std::runtime_error(
                        "--fit-strategy requires joint_minuit or staged_feasible");
            } else if (key == "--model-initial") {
                std::stringstream stream(next());
                std::string token;
                for (size_t j = 0; j < nps_xsec::kModelParameters; ++j) {
                    if (!std::getline(stream, token, ','))
                        throw std::runtime_error(
                            "--model-initial needs N_U,DeltaB_U,N_LT,N_TT");
                    options.initial[j] = std::stod(token);
                }
                if (std::getline(stream, token, ','))
                    throw std::runtime_error("Too many initial parameters");
            } else if (key == "--model-fix-u-slope") {
                options.fix_u_slope = true;
            } else if (key == "--model-starts") {
                options.starts = std::stoi(next());
            } else if (key == "--model-max-evaluations") {
                options.max_calls = std::stoi(next());
            } else if (key == "--model-max-iterations") {
                options.max_iterations = std::stoi(next());
            } else if (key == "--model-tolerance") {
                options.tolerance = std::stod(next());
            } else if (key == "--fixed-default-model" || key == "--model-free" ||
                       key == "--normalize_mmiss" || key == "--normalize-mmiss" ||
                       key == "--simc-yield-scale") {
                throw std::runtime_error(
                    "Obsolete ratio-method option: " + key +
                    "; the proxy uses the reference absolute response");
            } else {
                if (key == "--out-dir") out_dir_explicit = true;
                if (key == "--out-root") output_explicit[0] = true;
                if (key == "--out-csv") output_explicit[1] = true;
                if (key == "--out-slice-csv") output_explicit[2] = true;
                if (key == "--all-plots-pdf") output_explicit[3] = true;

                args.push_back(argv[i]);
                if (!key.empty() && key[0] == '-' && shared_option_takes_value(key)) {
                    if (i + 1 < argc) args.push_back(argv[++i]);
                } else if (key.empty() || key[0] != '-') {
                    if (positional_count < 2) output_explicit[positional_count] = true;
                    ++positional_count;
                }
                if (key == "--help" || key == "-h") {
                    std::cout
                        << "SigParam2021-inspired pi0: --model-initial "
                           "N_U,DeltaB_U,N_LT,N_TT --model-starts N (default 6)\n"
                        << "Pivoted U slope in [-20,20] GeV^-2; N_U>0; LT/TT signed. "
                           "--model-fix-u-slope fixes slope=0 for validation.\n"
                        << "--fit-strategy joint_minuit (default) or staged_feasible "
                           "(Gaussian positive-xsec diagnostic central values; no confidence errors).\n";
                }
            }
        }

        if (options.starts < 1 || options.max_calls < 1 || options.max_iterations < 1 ||
            !(options.tolerance > 0.0))
            throw std::runtime_error("Invalid model numerical controls");
        for (double parameter : options.initial)
            if (!std::isfinite(parameter))
                throw std::runtime_error("Nonfinite initial model parameter");

        gROOT->SetBatch(kTRUE);
        TH1::SetDefaultSumw2(kTRUE);
        gStyle->SetOptStat(0);
        auto cfg = parse_xsec_config(static_cast<int>(args.size()), args.data());
        const auto& definition = nps_xsec::xsec_model();
        if (cfg.model_identifier != definition.id)
            throw std::runtime_error(
                "Configuration must identify the current model: " + definition.id);
        for (size_t j = 0; j < nps_xsec::kModelParameters; ++j)
            if (options.initial[j] < definition.lower[j] ||
                options.initial[j] > definition.upper[j])
                throw std::runtime_error(
                    "Initial parameter outside model bounds: " + definition.names[j]);

        // The shared parser intentionally owns all reference defaults. Map
        // only defaults (never explicit paths) into the model namespace.
        if (!out_dir_explicit) {
            fs::path directory(cfg.out_dir);
            const std::string leaf = directory.filename().string();
            if (leaf == "xsec")
                cfg.out_dir = (directory.parent_path() / "xsec_simc_model").string();
            else if (leaf == "output_pi0_xsec_no_simc_model")
                cfg.out_dir = (directory.parent_path() /
                               "output_pi0_xsec_simc_model").string();
        }
        std::array<std::string*, 4> output_paths = {
            &cfg.out_root, &cfg.out_csv, &cfg.out_slice_csv, &cfg.out_all_plots_pdf
        };
        for (size_t index = 0; index < output_paths.size(); ++index) {
            if (output_explicit[index]) continue;
            const fs::path original(*output_paths[index]);
            *output_paths[index] =
                (fs::path(cfg.out_dir) / model_output_name(original.filename().string())).string();
        }
        if (cfg.verbose) {
            std::cout << "SIMC-model resolved output configuration:\n"
                      << "  out dir:   " << cfg.out_dir << '\n'
                      << "  out root:  " << cfg.out_root << '\n'
                      << "  out csv:   " << cfg.out_csv << '\n'
                      << "  out slice: " << cfg.out_slice_csv << '\n'
                      << "  plots pdf: " << cfg.out_all_plots_pdf << std::endl;
        }

        if (cfg.prepare_forward_inputs || !cfg.joint_plot_input.empty())
            throw std::runtime_error(
                "Model entry point requires a detector-level model fit; "
                "use reference entry for forward-cache or joint-plot operations");

        ExclPi0XSecAnalysis analysis(cfg);
        analysis.enable_proxy(options);
        analysis.Run();
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "[FATAL] " << error.what() << '\n';
        return 1;
    }
}
