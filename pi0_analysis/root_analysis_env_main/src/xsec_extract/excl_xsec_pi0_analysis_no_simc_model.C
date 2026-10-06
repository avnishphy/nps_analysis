// Cross-section steering file: stage ordering is intentionally visible here.
// Implementations are grouped by responsibility in xsec_*.h; the standalone
// translation unit keeps the existing ROOT/C++ executable and wrapper interface.
// The wrapper generates xsec_config.h from a JSON preset in xsec_config/.
#include "xsec_config.h"
#include "xsec_input.h"
#include "xsec_forward_inputs.h"
#include "xsec_binning.h"
#include "xsec_accumulation.h"
#include "xsec_fit.h"
#include "xsec_scaled_poisson.h"
#include "xsec_experimental_points.h"
#include "xsec_positive_toys.h"
#include "xsec_plot_global.h"
#include "xsec_plot_slices.h"
#include "xsec_plot_migration.h"
#include "xsec_plot_migration_coverage.h"
#include "xsec_model.h"
#include "xsec_output.h"
#include "xsec_joint_inputs.h"
#include "xsec_joint_plot.h"
#include "xsec_cli.h"

void ExclPi0XSecAnalysis::Run() {
    // Keep mode-specific stages from mutating an output directory or reading
    // full event trees before their incompatibility is diagnosed. The model
    // entry point also checks these at its CLI boundary; retain the invariant
    // here for the standalone no-model entry and programmatic callers.
    if (!event_model() && cfg.prepare_joint_inputs)
        die("Joint M0 preparation requires the SIMC-model entry point.");
    if (event_model() && (cfg.prepare_forward_inputs || !cfg.joint_plot_input.empty()))
        die("SIMC-model extraction does not support forward-cache or joint-plot stages.");
    validate_xsec_binning(cfg);
    // 1. Resolve the measured sample and its matched generated MC. Binning is
    // frozen before response construction and reused for both coordinate sets.
    if (cfg.prepare_joint_inputs && fs::exists(cfg.out_dir) && !fs::is_empty(cfg.out_dir))
        die("Joint input preparation requires a new or empty output directory: " + cfg.out_dir);
    fs::create_directories(cfg.out_dir);
    if (!cfg.prepare_joint_inputs && !cfg.prepare_forward_inputs) init_combined_pdf();
    load_input();
    load_mass_cut();
    detect_optional_branches();
    if (cfg.prepare_forward_inputs) {
        begin_forward_inputs();
        fill_from_trees();
        finish_forward_inputs();
        log("Forward event cache complete; no binned fit or plotting was run.");
        return;
    }
    build_binning();
    init_storage();
    // 2. Accumulate physical event weights, keeping reconstruction and origin
    // separate. Area matching below is restricted to missing-mass shape QA.
    fill_from_trees();
    compute_mmiss_area_scales();
    // 3. Apply the independent target correction to data, then fit all
    // published blocks jointly. Low-tprime feed-in is fitted; other exterior
    // faces retain the shared fixed nominal contribution.
    compute_ratios_and_xsec();
    if (cfg.prepare_joint_inputs) {
        write_joint_inputs();
        log("Joint inputs prepared; no individual fit was run.");
        return;
    }
    if (cfg.joint_plot_input.empty()) fit_slices();
    else load_joint_plot_input();
    compute_experimental_points();
    if (cfg.joint_plot_input.empty()) compute_positive_toy_errors();
    // 4. GK is a comparison evaluated after extraction. It neither supplies
    // fitted coefficients nor changes the detector response/normalization.
    compute_partons_projection();

    fout = TFile::Open(cfg.out_root.c_str(), "RECREATE");
    if (!fout || fout->IsZombie())
        die("Cannot open output ROOT file.");

    // 5. Plot and persist both the fitted quantities and the complete response
    // problem. Retain covariance, guard mapping and failure/diagnostic flags.
    make_global_plots();
    make_mmiss_comparison_plots();
    make_mass_selection_plot();
    make_mmiss_selection_1d_plot();
    make_epsilon_plots();
    make_slice_plots();
    make_sigma_vs_tprime_plots();
    make_partons_projection_plots();
    // Migration figures inspect the response itself, including exterior
    // feed-in. They do not rescale A or alter the fitted coefficients.
    make_migration_plots();
    make_migration_coverage_plots();
    close_combined_pdf();
    if (!cfg.joint_plot_input.empty()) {
        // The authoritative fit/covariance remain in the parent joint output.
        // This file contains only plotting diagnostics and imported-fit slices.
        fout->cd();
        TNamed provenance("joint_plot_source", cfg.joint_plot_input.c_str());
        provenance.Write();
        write_experimental_points();
        write_slice_csv();
        write_mmiss_comparison_csv();
        fout->Write();
        log("Joint-result plots complete; no individual fit or refit toys were run.");
        return;
    }
    write_results();
    write_experimental_points();
    write_csv();
    write_slice_csv();
    write_scaled_poisson_diagnostics();
    write_mmiss_comparison_csv();
    // cleanup() is now handled by the destructor; do not call here to avoid double-free.

    if (successful_fit_groups == 0)
        die("No Q2/xB subset produced a valid global fit; failure diagnostics were written to " + cfg.out_dir);
    log(fit_fallback ? "Analysis complete with excluded Q2/xB bins; see fit_status.csv." : "Analysis complete.");
}

int main(int argc, char **argv) {
    try {
        // Batch mode prevents graphics initialization from requiring a display.
        gROOT->SetBatch(kTRUE);
        TH1::SetDefaultSumw2(kTRUE);
        gStyle->SetOptStat(0);
        const AnalysisConfig cfg = parse_xsec_config(argc, argv);
        ExclPi0XSecAnalysis analysis(cfg);
        analysis.Run();
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "[FATAL] " << error.what() << std::endl;
        return 1;
    }
}
