// acceptance_cuts.h
// Loads and applies kinematic-dependent acceptance cuts from an INI-style config.
// No external parser dependency is required.

#ifndef ACCEPTANCE_CUTS_H
#define ACCEPTANCE_CUTS_H

#include <string>
#include <limits>

struct HMSDataCuts {
    double edtm_tdc_max = 0.1;
    double react_z_min = -4.0;
    double react_z_max = 4.0;
    double delta_min = -15.0;
    double delta_max = 15.0;
    double gtr_th_min = -0.1;
    double gtr_th_max = 0.1;
    double gtr_ph_min = -0.04;
    double gtr_ph_max = 0.04;
    double cer_npe_sum_min = 1.5;
    double cal_etotnorm_min = 0.7;
    double cal_etotnorm_max = 1.2;
};

struct HMSSimulationCuts {
    double react_z_min = -4.0;
    double react_z_max = 4.0;
    double delta_min = -15.0;
    double delta_max = 15.0;
    double gtr_th_min = -0.1;
    double gtr_th_max = 0.1;
    double gtr_ph_min = -0.04;
    double gtr_ph_max = 0.04;
};

struct NPSClusterCuts {
    double energy_min = 0.8;
    double x_min = -30.0;
    double x_max = 30.0;
    double y_min = -36.0;
    double y_max = 36.0;
    double time_center_ns = 150.0;

    // If false, caller fallback is used.
    bool has_time_halfwidth_override = false;
    double time_halfwidth_ns_override = std::numeric_limits<double>::quiet_NaN();
};

struct PairingCuts {
    // If false, caller fallback is used.
    bool has_time_diff_max_override = false;
    double time_diff_max_ns_override = std::numeric_limits<double>::quiet_NaN();
};

struct TimingWindowCuts {
    // If false, caller fallback is used.
    bool has_shifted_sidebands_override = false;
    bool use_shifted_sidebands_override = false;
};

struct WeightedCuts {
    double mmiss_corr_min = 0.8;
    double mmiss_corr_max = 1.1;
};

class AcceptanceCuts {
public:
    // Throws std::runtime_error if file, section, or values are missing/invalid.
    AcceptanceCuts(const std::string& config_file, const std::string& kin_setting);

    const std::string& source_file() const noexcept { return source_file_; }
    const std::string& kin_setting() const noexcept { return kin_setting_; }

    const HMSDataCuts& hms_data() const noexcept { return hms_data_; }
    const HMSSimulationCuts& hms_simulation() const noexcept { return hms_simulation_; }
    const NPSClusterCuts& nps_cluster() const noexcept { return nps_cluster_; }
    const PairingCuts& pairing() const noexcept { return pairing_; }
    const TimingWindowCuts& timing() const noexcept { return timing_; }
    const WeightedCuts& weighted() const noexcept { return weighted_; }

    double resolved_cluster_time_halfwidth_ns(double fallback_halfwidth_ns) const noexcept;
    double resolved_pair_time_diff_max_ns(double fallback_pair_time_diff_ns) const noexcept;
    bool resolved_use_shifted_sidebands(bool fallback_use_shifted) const noexcept;

    bool pass_hms_data(double edtm_tdc,
                       double h_delta,
                       double h_gtr_th,
                       double h_gtr_ph,
                       double h_cer_npe_sum,
                       double h_cal_etotnorm,
                       double h_react_z) const noexcept;

    bool pass_hms_simulation(double h_delta,
                             double h_gtr_th,
                             double h_gtr_ph,
                             double h_react_z) const noexcept;

    bool pass_nps_cluster(double clusE,
                          double clusX,
                          double clusY,
                          double clusT,
                          double resolved_time_halfwidth_ns) const noexcept;

    bool pass_weighted_exclusive(double mmiss_corr) const noexcept;

    std::string summary() const;

private:
    std::string source_file_;
    std::string kin_setting_;
    HMSDataCuts hms_data_;
    HMSSimulationCuts hms_simulation_;
    NPSClusterCuts nps_cluster_;
    PairingCuts pairing_;
    TimingWindowCuts timing_;
    WeightedCuts weighted_;
};

#endif // ACCEPTANCE_CUTS_H
