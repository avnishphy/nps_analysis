#pragma once

// Analysis state and stage interfaces. Histograms belong to one analysis instance, never globals.
#include "xsec_physics.h"
#include "xsec_response.h"
#include "xsec_mass_cut.h"
#include "xsec_proxy_solver.h"

class ExclPi0XSecAnalysis {
    // Test-only access to archived responses; no production input bypass.
    friend struct ProxyValidationAccess;

public:
    ~ExclPi0XSecAnalysis() {
        // Reset histograms before closing files
        h_q2_data.reset();
        h_q2_sim.reset();
        h_xb_data.reset();
        h_xb_sim.reset();
        h_tprime_data.reset();
        h_tprime_sim.reset();
        h_phi_data.reset();
        h_phi_sim.reset();
        h_q2_xb_data.reset();
        h_q2_xb_sim.reset();
        h_tprime_phi_data.reset();
        h_tprime_phi_sim.reset();
        // Close and null files
        cleanup();
        t_sim = nullptr;
        t_data = nullptr;
    }
    explicit ExclPi0XSecAnalysis(const AnalysisConfig& c) : cfg(c) {}
    void Run();
    void enable_proxy(const nps_xsec::ProxyOptions& options) {
        model_fit_mode=true; proxy_options=options;
    }

private:
    bool model_fit_mode=false;
    bool synthetic_validation=false;
    nps_xsec::ProxyOptions proxy_options;
    nps_xsec::ProxyResult proxy_result;
    std::vector<nps_xsec::ModelEvent> model_events;
    nps_xsec::ProxyProblem model_all_rows;
    std::vector<nps_xsec::ProxyProblem::RowEvaluation> model_rows;
    bool event_model() const {return model_fit_mode;}
    double model_row_variance(size_t r) const {
        if(positivity_boundary_active)return std::numeric_limits<double>::quiet_NaN();
        const auto& j=model_rows.at(r).jacobian;double v=0.;
        for(size_t a=0;a<j.size();++a)for(size_t b=0;b<j.size();++b)v+=j[a]*proxy_result.covariance[a*j.size()+b]*j[b];
        return v;
    }
    nps_xsec::ModelEvaluation model_curve(double tp) const;
    void write_sigparam_diagnostics();
    void fit_proxy_subset(const std::vector<bool>& groups);
    void write_proxy_diagnostics(const nps_xsec::ProxyProblem& problem);
    AnalysisConfig cfg;
    CutFlow cutflow;
    XsecMassGeometry mass_geometry;

    TFile* f_sim = nullptr;
    TFile* f_data = nullptr;
    TFile* f_vertex = nullptr;
    TTree* t_sim = nullptr;
    TTree* t_data = nullptr;
    TTree* t_vertex = nullptr;
    TFile* fout = nullptr;

    bool has_helicity = false;
    bool has_sim_helicity = false;
    bool has_model_xsec = false;
    bool has_vertex_kinematics = false;
    bool vertex_from_raw_simc = false;
    long long n_vertex_matched = 0;
    std::string model_xsec_branch;

    std::string combined_pdf_path;
    std::ofstream forward_data_stream, forward_mc_stream;
    unsigned long long forward_event_id = 0;
    long long forward_data_count = 0, forward_mc_count = 0;
    void begin_forward_inputs();
    void finish_forward_inputs();
    std::vector<std::string> combined_pdf_pages;

    std::vector<double> phi_edges, tprime_edges, q2_edges, xb_edges;
    std::vector<std::vector<double>> xb_edges_by_q2;
    std::vector<SliceResult> slices;
    // Response row = reconstructed slice*n_phi+phi. Truth blocks consist of
    // published bins followed by six exterior regions; only tprime_below is
    // a fitted exterior block. Q2/xB feed-in remains in fixed row offsets.
    std::vector<std::vector<nps_xsec::ResponseCell>> migration_response;
    // Sparse cell key = truth_block*n_phi + vertex_phi_bin. All cells share
    // their parent block's three coefficients; no new fit parameters enter.
    std::vector<std::map<int, nps_xsec::ResponseCell>> truth_phi_response;
    std::vector<nps_xsec::TruthMoments> truth_phi_moments;
    std::vector<ExperimentalPoint> experimental_points;
    std::vector<double> experimental_point_covariance;
    std::vector<nps_xsec::TruthMoments> truth_moments;
    std::vector<int> active_truth_blocks, fixed_truth_blocks, fit_rows;
    std::vector<double> fixed_feedin_prediction, fixed_feedin_mc_variance;
    std::vector<std::vector<double>> response_design;
    nps_xsec::LinearSolution migration_fit;
    // Boundary-constrained estimates do not have the unconstrained Gaussian
    // covariance. Retain that inverse curvature only as a named diagnostic.
    std::vector<double> fit_curvature_inverse;
    // Conditional constrained-refit toy spread, kept distinct from the
    // unavailable Gaussian covariance when a positivity boundary is active.
    std::vector<double> positive_toy_parameter_covariance;
    std::vector<double> positive_toy_point_covariance;
    int positive_toys_successful = 0;
    std::vector<double> positivity_boundary_tolerances, positivity_feasibility_tolerances;
    bool positivity_boundary_active = false;
    int positivity_iterations = 0;
    std::vector<double> fit_variance;
    std::vector<XsecScaledRow> scaled_rows;
    std::map<int, XsecWeightMoments> run_weight_moments;
    int scaled_minuit_status = -1, scaled_covariance_status = -1;
    double scaled_edm = std::numeric_limits<double>::quiet_NaN();
    unsigned int scaled_calls = 0;
    void fit_scaled_poisson_subset(const std::vector<bool>& groups);
    void finalize_fit_subset(const std::vector<bool>& groups);
    void write_scaled_poisson_diagnostics();
    int mc_iterations=0, omitted_zero_variance_rows=0, nonphysical_truth_bins=0;
    bool mc_converged=false;
    struct FitAttempt {
        std::vector<bool> groups; // iq*n_xb+ix; all t' bins retained together.
        bool ok = false;
        std::string reason;
    };
    std::vector<FitAttempt> fit_attempts;
    std::vector<bool> retained_fit_groups;
    bool fit_fallback = false;
    int successful_fit_groups = 0;
    // Full-range missing-mass spectra are diagnostics only: all data
    // candidates versus area-normalized smeared generated-exclusive SIMC.
    // No background simulation or subtraction enters.
    std::vector<MissingMassSliceDiagnostic> mmiss_slices;
    std::unique_ptr<TH2D> h_mass_data_all, h_mass_data_selected;
    std::unique_ptr<TH2D> h_mass_sim_all, h_mass_sim_selected;
    std::unique_ptr<TH1D> h_mmiss_corr_data_all, h_mmiss_corr_data_selected;
    std::unique_ptr<TH1D> h_mmiss_reco_data_all, h_mmiss_reco_data_selected;
    std::unique_ptr<TH1D> h_mmiss_sim_all, h_mmiss_sim_selected;


    // Global QA histograms
    std::unique_ptr<TH1D> h_q2_data, h_q2_sim, h_xb_data, h_xb_sim, h_tprime_data, h_tprime_sim, h_phi_data, h_phi_sim;
    std::unique_ptr<TH2D> h_q2_xb_data, h_q2_xb_sim, h_tprime_phi_data, h_tprime_phi_sim;
    // The same accepted MC events viewed at the vertex and reconstruction.
    // Counts diagnose migration support, not cross sections or efficiencies.
    std::unique_ptr<TH2D> h_migration_vertex_q2_xb, h_migration_reco_q2_xb;

    int slice_index(int it, int iq, int ix) const {
        return (it * cfg.n_q2 + iq) * cfg.n_xb + ix;
    }

    SliceResult& slice(int it, int iq, int ix) {
        return slices[slice_index(it, iq, ix)];
    }
    const SliceResult& slice(int it, int iq, int ix) const {
        return slices[slice_index(it, iq, ix)];
    }

    void log(const std::string& s) const { if (cfg.verbose) std::cout << "[INFO] " << s << "\n"; }
    void warn(const std::string& s) const { std::cerr << "[WARN] " << s << "\n"; }
    [[noreturn]] void die(const std::string& s) const { throw std::runtime_error(s); }

    void load_input();
    void load_mass_cut();
    void detect_optional_branches();
    void build_binning();
    void init_storage();
    void fill_mmiss_data_diagnostic(double q2, double t, double tmin, double xb, double phi,
                                    double mmiss, double pi0_weight, float scale,
                                    double charge_uC, double total_charge_uC);
    void fill_mmiss_sim_diagnostic(float q2, float t, float tmin, float xb, float phi,
                                   float mmiss, float full_weight, int exclusive);
    void compute_mmiss_area_scales();
    void make_mmiss_comparison_plots();
    void make_mass_selection_plot();
    void make_mmiss_selection_1d_plot();
    void write_mmiss_comparison_csv();
    void fill_from_trees();
    void compute_ratios_and_xsec();
    void fit_slices();
    void compute_experimental_points();
    void compute_positive_toy_errors();
    double plot_parameter_variance(const std::vector<double>& gradient) const;
    double plot_coefficient_error(int truth_block, int term) const;
    double plot_point_error(int row) const;
    void write_experimental_points();
    void fit_global_subset(const std::vector<bool>& groups);
    void invalidate_fit_results();
    void compute_partons_projection();
    void make_global_plots();
    void make_migration_plots();
    void make_migration_coverage_plots();
    void make_epsilon_plots();
    void make_slice_plots();
    void make_sigma_vs_tprime_plots();
    void make_partons_projection_plots();
    void write_results();
    void write_joint_inputs();
    void load_joint_plot_input();
    void write_migration_results();
    void write_csv();
    void write_slice_csv();
    void cleanup();
    void init_combined_pdf();
    void append_to_combined_pdf(const std::string& page_path);
    void close_combined_pdf();
    std::string format_slice_bin_label(int it, int iq, int ix) const;
    void draw_slice_bin_label(int it, int iq, int ix, double y_ndc = 0.92) const;

    void fill_data_event(double q2, double t, double tmin, double xb, double phi, double pi0_weight, float scale, double charge_uC, double total_charge_uC, double mmiss_all, double mpi0_all, int exclusive_flag, int helicity, bool use_helicity, double W, int run_number);
    void fill_sim_event(float q2,
                        float t,
                        float tmin,
                        float xb,
                        float phi,
                        float full_weight,
                        int is_exclusive,
                        float mmiss,
                        float mpi0,
                        float model_xsec,
                        float vertex_q2,
                        float vertex_W,
                        float vertex_t,
                        float vertex_phi,
                        float vertex_hsxptari,
                        float vertex_hsyptari,
                        int helicity,
                        bool use_helicity,
                        float W);

    bool slice_passes_kin(const float q2, const float xb, const double tprime) const {
        return (q2 >= cfg.q2_min && q2 <= cfg.q2_max &&
                xb >= cfg.xb_min && xb <= cfg.xb_max &&
                tprime >= cfg.tprime_min && tprime <= cfg.tprime_max &&
                xsec_inside_diamond(cfg, xb, q2));
    }

    bool selected_data(double mmiss, double mpi0, int exclusive_flag) const {
        if (cfg.mmiss_select != "window")
            return exclusive_flag != 0; // exact decision stored by the combined-data producer
        return std::isfinite(mmiss) && mmiss >= cfg.mmiss_lower_gev &&
               mmiss <= cfg.mmiss_upper_gev;
    }
    bool selected_sim(int generated_exclusive, double mmiss, double mpi0) const {
        // Only exclusive generated MC has a usable pi0 sigcm denominator.
        if (!generated_exclusive) return false;
        if (cfg.mmiss_select != "window")
            return mass_geometry.contains(mpi0, mmiss);
        return std::isfinite(mmiss) && mmiss >= cfg.mmiss_lower_gev &&
               mmiss <= cfg.mmiss_upper_gev;
    }

    double calc_tprime(float t, float tmin) const { return static_cast<double>(t) - static_cast<double>(tmin); }

    double calc_w(float q2, float xb) const {
        return std::sqrt(std::max(0.0, q2_xb_to_w2(q2, xb, cfg.mp)));
    }

    void accumulate_global_histograms(float q2, float xb, double tprime, double phi, double weight_data, double weight_sim);
    void write_canvas_pdf_png(TCanvas* c, const std::string& base);
};
