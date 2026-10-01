#pragma once

// CLI and environment parsing. Reject obsolete or ambiguous physics settings instead of silently changing them.
#include "xsec_config.h"

static std::string getenv_or_default(const char* key, const std::string& fallback = "") {
    const char* v = std::getenv(key);
    if (!v || !v[0]) return fallback;
    return std::string(v);
}

static std::string trim_copy_cli(const std::string& s) {
    size_t first = 0;
    while (first < s.size() && std::isspace(static_cast<unsigned char>(s[first]))) ++first;
    size_t last = s.size();
    while (last > first && std::isspace(static_cast<unsigned char>(s[last - 1]))) --last;
    return s.substr(first, last - first);
}

static std::string sanitize_token_cli(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (char c : s) {
        if (std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-') out.push_back(c);
        else out.push_back('_');
    }
    return out;
}

static std::string require_next_arg(int& i, int argc, char** argv, const std::string& opt_name) {
    if (i + 1 >= argc) throw std::runtime_error("Missing value for option " + opt_name);
    ++i;
    return std::string(argv[i]);
}

static void print_usage(const char* prog) {
    std::cout << "Usage: " << prog << R"( [options]

Primary options:
  --kin <Kin_old>                 Derive standard paths only; physics defaults
                                  remain those in xsec_config.h
  --target <name>                 Combined target token (default: LH2)
  --output-base <path>            Output base directory (default: $NPS_OUTPUT_BASE or output)
  --root-dir <path>               Root directory for combined/sim files (default from kin)
  --data-file <path>              Combined data ROOT file
  --sim-file <path>               Simulation ROOT file
  --vertex_simc_file <file|dir>   Required original exclusive SIMC h10 ROOT file or worksim directory
  --out-dir <path>                Output directory for xsec products
  --out-root <path>               Output ROOT file
  --out-csv <path>                Output CSV summary
  --out-slice-csv <path>          Output per-slice CSV summary
  --all-plots-pdf <path>          Combined plots PDF path

Binning/range options:
  --mmiss-lower <GeV>            Shared data/exclusive-SIMC lower bound (default 0.80)
  --mmiss-upper <GeV>            Shared data/exclusive-SIMC upper bound (default 1.10)
  --mmiss_select <mode>          mcd or ellipse (default: window)
  --mmiss-cut-file <path>        Combined-data geometry text (default: beside data ROOT)
  --target-contam <factor>       Data yield divisor (default 0.584; use 1 to omit)
  --target-contam-err <factor>   Absolute divisor uncertainty (default 0.014; use 0 with 1)
  --partons                      Project native GK06/GPDGK19 pi0 response points
  --partons-warmups <int>        DVMP CFF MC warm-up calls (default 10000)
  --partons-calls <int>          DVMP CFF MC calls (default 100000)
  Bin edges are fixed in xsec_config.h (phi, tprime, Q2, and xB per Q2).

Behavior options:
  --prepare-joint-inputs          Export joint-fit inputs only; no individual fit or plots
  --prepare-forward-inputs        Export selected events before rectangular cuts; no fit or plots
  --joint-plot-input <file>        Render imported joint coefficients; no individual fit or toys
  --ebeam <float>                 Beam energy
  --svd-rank-tolerance <float>    Relative singular-value rank cutoff (default 1e-10)
  --mc-max-iterations <int>       MC-variance fit iteration limit (default 30)
  --mc-fit-tolerance <float>      Relative MC-variance fit convergence (default 1e-6)
  --fit-variance <data|finite-mc> Eq. 5.23 data-only reference or finite-MC extension (default)
  --fit-objective <gaussian|scaled-poisson>  Fit statistic (default gaussian)
  --scaled-empty-scale <auto|slice|global>  Empty-row reference sensitivity
  --positive-xsec                 Require full angular cross section >= 0; LT/TT remain signed
  --no-positive-xsec              Disable positivity (default; overrides environment)
  --quiet                         Disable verbose logging
  --no-diagnostics                Omit QA, missing-mass, epsilon and migration plots; keep core fit plots
  --no-png                        Disable all PNG plots
  --no-pdf                        Disable individual and combined PDF plots
  --help                          Show this message

Compatibility:
  positional arg1                 out_root (legacy)
  positional arg2                 out_csv (legacy)

Positivity environment:
  NPS_XSEC_POSITIVE_XSEC=0|1       Default off; last explicit on/off flag wins
)";
}


inline AnalysisConfig parse_xsec_config(int argc, char** argv) {
        AnalysisConfig cfg;
        cfg.partons_executable_path = argv[0]; // PARTONS locates properties beside this binary.

        std::string cli_kin = getenv_or_default("NPS_XSEC_KIN", "");
        std::string cli_target = getenv_or_default("NPS_XSEC_TARGET", "LH2");
        std::string cli_output_base = getenv_or_default("NPS_OUTPUT_BASE", "output");
        std::string cli_root_dir = getenv_or_default("NPS_XSEC_ROOT_DIR", "");
        std::string positive_mode = getenv_or_default("NPS_XSEC_POSITIVE_XSEC", "0");

        bool mmiss_select_set = false;
        bool explicitly_disabled_positivity = false;
        bool data_file_set = false;
        bool sim_file_set = false;
        bool out_dir_set = false;
        bool out_root_set = false;
        bool out_csv_set = false;
        bool out_slice_csv_set = false;
        bool out_all_plots_pdf_set = false;

        auto apply_env_path = [&](const char* key, std::string& dst, bool& set_flag) {
            std::string v = getenv_or_default(key, "");
            if (!v.empty()) {
                dst = v;
                set_flag = true;
            }
        };

        apply_env_path("NPS_XSEC_DATA_FILE", cfg.data_file, data_file_set);
        apply_env_path("NPS_XSEC_SIM_FILE", cfg.simc_file, sim_file_set);
        cfg.vertex_simc_file = getenv_or_default("NPS_XSEC_VERTEX_SIMC_FILE", "");
        apply_env_path("NPS_XSEC_OUT_DIR", cfg.out_dir, out_dir_set);
        apply_env_path("NPS_XSEC_OUT_ROOT", cfg.out_root, out_root_set);
        apply_env_path("NPS_XSEC_OUT_CSV", cfg.out_csv, out_csv_set);
        apply_env_path("NPS_XSEC_OUT_SLICE_CSV", cfg.out_slice_csv, out_slice_csv_set);
        apply_env_path("NPS_XSEC_ALL_PLOTS_PDF", cfg.out_all_plots_pdf, out_all_plots_pdf_set);

        int positional_idx = 0;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--help" || arg == "-h") { print_usage(argv[0]); std::exit(0); }
            if (arg == "--kin") { cli_kin = require_next_arg(i, argc, argv, arg); continue; }
            if (arg == "--target") { cli_target = require_next_arg(i, argc, argv, arg); continue; }
            if (arg == "--output-base") { cli_output_base = require_next_arg(i, argc, argv, arg); continue; }
            if (arg == "--root-dir") { cli_root_dir = require_next_arg(i, argc, argv, arg); continue; }
            if (arg == "--data-file") { cfg.data_file = require_next_arg(i, argc, argv, arg); data_file_set = true; continue; }
            if (arg == "--sim-file") { cfg.simc_file = require_next_arg(i, argc, argv, arg); sim_file_set = true; continue; }
            if (arg == "--vertex_simc_file" || arg == "--vertex-simc-file") {
                cfg.vertex_simc_file = require_next_arg(i, argc, argv, arg); continue;
            }
            if (arg == "--out-dir") { cfg.out_dir = require_next_arg(i, argc, argv, arg); out_dir_set = true; continue; }
            if (arg == "--out-root") { cfg.out_root = require_next_arg(i, argc, argv, arg); out_root_set = true; continue; }
            if (arg == "--out-csv") { cfg.out_csv = require_next_arg(i, argc, argv, arg); out_csv_set = true; continue; }
            if (arg == "--out-slice-csv") { cfg.out_slice_csv = require_next_arg(i, argc, argv, arg); out_slice_csv_set = true; continue; }
            if (arg == "--all-plots-pdf") { cfg.out_all_plots_pdf = require_next_arg(i, argc, argv, arg); out_all_plots_pdf_set = true; continue; }
            if (arg == "--mmiss-lower") { cfg.mmiss_lower_gev = std::stod(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--mmiss-upper") { cfg.mmiss_upper_gev = std::stod(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--mmiss_select" || arg == "--mmiss-select") {
                cfg.mmiss_select = require_next_arg(i, argc, argv, arg);
                mmiss_select_set = true;
                continue;
            }
            if (arg == "--mmiss-cut-file") { cfg.mmiss_cut_file = require_next_arg(i, argc, argv, arg); continue; }
            if (arg == "--target-contam") { cfg.tgt_contam = std::stod(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--target-contam-err") { cfg.tgt_contam_err = std::stod(require_next_arg(i, argc, argv, arg)); continue; }
            // A data/SIMC area match would use the measured yield to set its
            // own response normalization and override the target correction.
            // Reject historical commands explicitly rather than reinterpret them.
            if (arg == "--normalize_mmiss" || arg == "--normalize-mmiss")
                throw std::runtime_error(arg + " was removed: use the physical --target-contam correction; data/SIMC area matching is not an absolute cross-section normalization.");
            if (arg == "--partons") { cfg.partons_projection = true; continue; }
            if (arg == "--partons-warmups") { cfg.partons_warmups = std::stoi(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--partons-calls") { cfg.partons_calls = std::stoi(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--ebeam") { cfg.ebeam = std::stod(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--svd-rank-tolerance") { cfg.rank_tolerance = std::stod(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--mc-max-iterations") { cfg.mc_max_iterations = std::stoi(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--mc-fit-tolerance") { cfg.mc_fit_tolerance = std::stod(require_next_arg(i, argc, argv, arg)); continue; }
            if (arg == "--fit-variance") { cfg.fit_variance_mode = require_next_arg(i, argc, argv, arg); continue; }
            if (arg == "--fit-objective") { cfg.fit_objective = require_next_arg(i, argc, argv, arg); continue; }
            if (arg == "--scaled-empty-scale") { cfg.scaled_empty_scale = require_next_arg(i, argc, argv, arg); continue; }
            if (arg == "--positive-xsec") { positive_mode = "1"; continue; }
            if (arg == "--no-positive-xsec") { positive_mode = "0"; explicitly_disabled_positivity = true; continue; }
            if (arg == "--quiet") { cfg.verbose = false; continue; }
            if (arg == "--prepare-joint-inputs") { cfg.prepare_joint_inputs = true; continue; }
            if (arg == "--prepare-forward-inputs") { cfg.prepare_forward_inputs = true; continue; }
            if (arg == "--joint-plot-input") { cfg.joint_plot_input = require_next_arg(i, argc, argv, arg); continue; }
            if (arg == "--no-diagnostics") { cfg.diagnostics = false; continue; }
            if (arg == "--no-png") { cfg.write_png = false; continue; }
            if (arg == "--no-pdf") { cfg.write_pdf = false; continue; }

            if (!arg.empty() && arg[0] == '-') throw std::runtime_error("Unknown option: " + arg);
            if (positional_idx == 0) { cfg.out_root = arg; out_root_set = true; }
            else if (positional_idx == 1) { cfg.out_csv = arg; out_csv_set = true; }
            else throw std::runtime_error("Unexpected positional argument: " + arg);
            ++positional_idx;
        }

        // Validate the effective setting after explicit flags override the
        // environment, matching the shell wrapper's precedence.
        if (positive_mode != "0" && positive_mode != "1")
            throw std::runtime_error("NPS_XSEC_POSITIVE_XSEC must be 0 or 1 (or override with --positive-xsec/--no-positive-xsec).");
        cfg.positive_xsec = positive_mode == "1";
        if (cfg.prepare_forward_inputs) {
            if (cfg.prepare_joint_inputs || !cfg.joint_plot_input.empty() || cfg.partons_projection)
                throw std::runtime_error("Forward event export cannot be combined with joint stages or PARTONS.");
            cfg.write_png = cfg.write_pdf = cfg.diagnostics = false;
        }
        if (cfg.prepare_joint_inputs) {
            if (!cfg.joint_plot_input.empty())
                throw std::runtime_error("Joint preparation and joint plotting are separate stages.");
            if (cfg.partons_projection)
                throw std::runtime_error("Joint input preparation does not compute PARTONS projections.");
            cfg.write_png = cfg.write_pdf = cfg.diagnostics = false;
        }
        if (cfg.scaled_empty_scale!="auto" && cfg.scaled_empty_scale!="slice" && cfg.scaled_empty_scale!="global")
            throw std::runtime_error("--scaled-empty-scale must be auto, slice, or global.");
        if (cfg.fit_objective != "gaussian" && cfg.fit_objective != "scaled-poisson")
            throw std::runtime_error("--fit-objective must be gaussian or scaled-poisson.");
        // The scaled likelihood always uses a physical angular cross section;
        // an explicit negative-physics request is incompatible with it.
        if (cfg.fit_objective == "scaled-poisson") {
            if (explicitly_disabled_positivity || (positive_mode == "0" && getenv_or_default("NPS_XSEC_POSITIVE_XSEC","").size()))
                throw std::runtime_error("scaled-poisson requires angular positivity.");
            cfg.positive_xsec = true;
        }

        cli_kin = trim_copy_cli(cli_kin);
        cli_target = trim_copy_cli(cli_target);
        cli_output_base = trim_copy_cli(cli_output_base);
        cli_root_dir = trim_copy_cli(cli_root_dir);
        if (cli_target.empty()) cli_target = "LH2";
        if (cli_output_base.empty()) cli_output_base = "output";

        if (!cli_kin.empty() && cli_kin != cfg.configured_kinematic) {
            std::cerr << "[WARN] --kin " << cli_kin << " changes paths only; physics defaults remain "
                      << cfg.configured_kinematic
                      << " from xsec_config.h unless explicitly overridden.\n";
        }

        if (!cli_kin.empty()) {
            const std::string kin_safe = sanitize_token_cli(cli_kin);
            if (kin_safe.empty()) throw std::runtime_error("Invalid --kin value after sanitization: '" + cli_kin + "'");
            if (cli_root_dir.empty()) cli_root_dir = (fs::path(cli_output_base) / kin_safe / "root").string();
            if (!out_dir_set) { cfg.out_dir = (fs::path(cli_output_base) / kin_safe / "xsec").string(); out_dir_set = true; }
        }

        if (!cli_root_dir.empty()) {
            if (!data_file_set) {
                cfg.data_file = (fs::path(cli_root_dir) / ("combined_branches_" + sanitize_token_cli(cli_target) + ".root")).string();
                data_file_set = true;
            }
            if (!sim_file_set) {
                const fs::path smeared = fs::path(cli_root_dir) / "simc_pi0_analysis_output_smeared.root";
                cfg.simc_file = smeared.string();
                sim_file_set = true;
            }
            if (!out_dir_set) { cfg.out_dir = (fs::path(cli_root_dir).parent_path() / "xsec").string(); out_dir_set = true; }
        }

        if (!out_dir_set || trim_copy_cli(cfg.out_dir).empty()) cfg.out_dir = "output_pi0_xsec_no_simc_model";
        if (!out_root_set || trim_copy_cli(cfg.out_root).empty()) cfg.out_root = (fs::path(cfg.out_dir) / "excl_xsec_pi0_analysis_no_simc_model_output.root").string();
        if (!out_csv_set || trim_copy_cli(cfg.out_csv).empty()) cfg.out_csv = (fs::path(cfg.out_dir) / "excl_xsec_pi0_analysis_no_simc_model_summary.csv").string();
        if (!out_slice_csv_set || trim_copy_cli(cfg.out_slice_csv).empty()) cfg.out_slice_csv = (fs::path(cfg.out_dir) / "excl_xsec_pi0_analysis_no_simc_model_slice_summary.csv").string();
        if (!out_all_plots_pdf_set || trim_copy_cli(cfg.out_all_plots_pdf).empty()) cfg.out_all_plots_pdf = (fs::path(cfg.out_dir) / "all_generated_plots_no_simc_model.pdf").string();

        cfg.vertex_simc_file = trim_copy_cli(cfg.vertex_simc_file);
        if (!cfg.vertex_simc_file.empty()) {
            fs::path vertex_path(cfg.vertex_simc_file);
            if (fs::is_directory(vertex_path)) {
                // Directory mode targets the exclusive generated channel.
                if (cli_kin.empty())
                    throw std::runtime_error("--vertex_simc_file directory mode requires --kin.");
                std::string kin_token = sanitize_token_cli(cli_kin);
                if (kin_token.rfind("KinC_", 0) == 0) kin_token.erase(0, 5);
                vertex_path /= "nps_excl_pi0_" + kin_token + ".root";
            }
            if (!fs::is_regular_file(vertex_path))
                throw std::runtime_error("Original exclusive SIMC ROOT file not found: " + vertex_path.string());
            cfg.vertex_simc_file = vertex_path.string();
        }

        if (trim_copy_cli(cfg.data_file).empty()) throw std::runtime_error("Data file path is empty. Provide --data-file or --kin/--root-dir.");
        if (trim_copy_cli(cfg.simc_file).empty()) throw std::runtime_error("Simulation file path is empty. Provide --sim-file or --kin/--root-dir.");
        if (!fs::exists(cfg.data_file)) throw std::runtime_error("Data ROOT file not found: " + cfg.data_file);
        if (!fs::exists(cfg.simc_file)) throw std::runtime_error("Simulation ROOT file not found: " + cfg.simc_file);
        if ((mmiss_select_set && cfg.mmiss_select != "mcd" && cfg.mmiss_select != "ellipse") ||
            (!mmiss_select_set && cfg.mmiss_select != "window"))
            throw std::runtime_error("--mmiss_select must be mcd or ellipse.");
        if (cfg.mmiss_cut_file.empty()) {
            const fs::path data(cfg.data_file);
            cfg.mmiss_cut_file = (data.parent_path() /
                (data.stem().string() + "_combined_2d_mass_cut_debug.txt")).string();
        }
        if (!(std::isfinite(cfg.mmiss_lower_gev) && std::isfinite(cfg.mmiss_upper_gev) &&
              cfg.mmiss_lower_gev >= 0.0 && cfg.mmiss_lower_gev < cfg.mmiss_upper_gev))
            throw std::runtime_error("Missing-mass bounds must be finite with 0 <= lower < upper.");
        if (!(std::isfinite(cfg.tgt_contam) && cfg.tgt_contam > 0.0))
            throw std::runtime_error("--target-contam must be positive and finite.");
        if (!(std::isfinite(cfg.tgt_contam_err) && cfg.tgt_contam_err >= 0.0))
            throw std::runtime_error("--target-contam-err must be nonnegative and finite.");
        // These tolerances detect unresolved response directions and failure to
        // converge; they must not be interpreted as physics regularization.
        if (!(std::isfinite(cfg.rank_tolerance) && cfg.rank_tolerance > 0.0 && cfg.rank_tolerance < 1.0))
            throw std::runtime_error("--svd-rank-tolerance must be finite and strictly between 0 and 1.");
        if (cfg.mc_max_iterations <= 0)
            throw std::runtime_error("--mc-max-iterations must be positive.");
        if (cfg.fit_variance_mode != "data" && cfg.fit_variance_mode != "finite-mc")
            throw std::runtime_error("--fit-variance must be data or finite-mc.");
        if (!(std::isfinite(cfg.mc_fit_tolerance) && cfg.mc_fit_tolerance > 0.0 && cfg.mc_fit_tolerance < 1.0))
            throw std::runtime_error("--mc-fit-tolerance must be finite and strictly between 0 and 1.");
        if (cfg.partons_warmups <= 0 || cfg.partons_calls <= 0)
            throw std::runtime_error("PARTONS MC warmups and calls must be positive.");
#ifndef NPS_ENABLE_PARTONS
        if (cfg.partons_projection)
            throw std::runtime_error("--partons requires a native PARTONS-linked executable.");
#endif
        validate_xsec_binning(cfg);

        if (cfg.prepare_joint_inputs) {
            std::cout << "Joint input preparation (no individual fit):\n"
                      << "  data: " << cfg.data_file << '\n'
                      << "  SIMC: " << cfg.simc_file << '\n'
                      << "  vertex SIMC: " << cfg.vertex_simc_file << '\n'
                      << "  selection: " << cfg.mmiss_select << '\n'
                      << "  target divisor: " << cfg.tgt_contam << '\n'
                      << "  input tables: " << cfg.out_dir << std::endl;
            return cfg;
        }

        std::cout << "XSec input configuration:" << std::endl
                  << "  data file: " << cfg.data_file << std::endl
                  << "  sim file:  " << cfg.simc_file << std::endl
                  << "  vertex SIMC: " << (cfg.vertex_simc_file.empty() ? "off" : cfg.vertex_simc_file) << std::endl
                  << "  selection: " << cfg.mmiss_select << " (generated-exclusive SIMC only)" << std::endl
                  << "  SIMC mass cut: " << (cfg.mmiss_select == "window" ?
                      "configured missing-mass window" : "combined-data geometry") << std::endl
                  << "  target factor: " << cfg.tgt_contam << " +/- " << cfg.tgt_contam_err << " (data divided by factor)" << std::endl
                  << "  SVD rank tolerance: " << cfg.rank_tolerance << std::endl
                  << "  MC variance iterations: " << cfg.mc_max_iterations
                  << ", convergence tolerance: " << cfg.mc_fit_tolerance << std::endl
                  << "  fit variance: " << (cfg.fit_objective=="scaled-poisson" ? "ignored (fixed SIMC response)" : cfg.fit_variance_mode) << std::endl
                  << "  fit objective: " << cfg.fit_objective << std::endl
                  << "  scaled empty scale: " << cfg.scaled_empty_scale << std::endl
                  << "  angular positivity: " << (cfg.positive_xsec ? "on (zero allowed; LT/TT signed)" : "off") << std::endl
                  << "  PARTONS projection: " << (cfg.partons_projection ? "on" : "off")
                  << " (warmups " << cfg.partons_warmups << ", calls " << cfg.partons_calls << ")" << std::endl
                  << "  out dir:   " << cfg.out_dir << std::endl
                  << "  out root:  " << cfg.out_root << std::endl
                  << "  out csv:   " << cfg.out_csv << std::endl
                  << "  out slice: " << cfg.out_slice_csv << std::endl
                  << "  plots pdf: " << cfg.out_all_plots_pdf << std::endl;

    return cfg;
}
