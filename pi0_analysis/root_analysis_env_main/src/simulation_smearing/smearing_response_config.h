#ifndef NPS_SMEARING_RESPONSE_CONFIG_H
#define NPS_SMEARING_RESPONSE_CONFIG_H

// Fixed detector-response settings shared by the smearing fitter and producer.
// Fitted coefficients and each program's RNG policy remain local/runtime inputs.
namespace SmearingResponseConfig {

constexpr double NONPOSITIVE_CLAMP = 1e-6;

// Producer-facing interpolated maps are sampled onto this fixed grid.  The
// fitter's producer-response diagnostic must use the same binning and lookup
// convention as simc_pi0_analysis.C (FindBin followed by in-range clamping).
constexpr int INTERPOLATED_MAP_NBINS_X = 100;
constexpr int INTERPOLATED_MAP_NBINS_Y = 100;
constexpr double SECTION_BOUNDARY_TOLERANCE_CM = 1e-9;
constexpr double PRODUCER_FALLBACK_SIGMA = 0.05;

// Photon response mean: E_mean = a + b*E_safe + c*ln(E_safe/1 GeV).
// In scalar-mu mode the fitter fixes a=c=0, leaving E_mean=b*E_safe.
constexpr bool ENABLE_ENERGY_DEPENDENT_MU = true;
constexpr const char* ENERGY_MEAN_MODEL = "a_plus_bE_plus_clnE";
constexpr double MU_ENERGY_MIN_GEV = 0.2;

// Per-photon response PDF; sigma is a standard deviation for Gaussian and
// a FWHM coefficient for Landau.
constexpr int SMEAR_SHAPE_GAUSSIAN = 0;
constexpr int SMEAR_SHAPE_LANDAU = 1;
constexpr int ENERGY_SMEAR_SHAPE = SMEAR_SHAPE_GAUSSIAN;

constexpr double LANDAU_FWHM_TO_SCALE = 4.0;
constexpr double LANDAU_SCALE_DENOMINATOR_MIN = 1e-9;
constexpr bool ENABLE_TRUNCATED_LANDAU = true;
constexpr int LANDAU_MAX_REDRAWS = 64;
constexpr double LANDAU_MAX_FWHM_ABOVE_MPV = 8.0;

// sigma_E = sigma*sqrt(E) when true. Otherwise use the three-term model with
// res_A/res_B/res_C loaded from the fit artifact or these fixed defaults.
constexpr bool USE_SIMPLE_STOCHASTIC_MODEL = true;
constexpr double RESOLUTION_A_DEFAULT = 0.97;
constexpr double RESOLUTION_B_DEFAULT = 1.1;
constexpr double RESOLUTION_C_DEFAULT = 1.14;

// Position response: additive Gaussian x/y shifts, optionally scaled as
// sigma_pos(E)=sigma_pos_ref*sqrt(E0/E_safe).
constexpr bool ENABLE_POSITION_SMEARING = true;
constexpr bool ENABLE_ENERGY_DEPENDENT_SIGMA_POS = false;
constexpr double SIGMA_POS_ENERGY_E0_GEV = 2.0;
constexpr double SIGMA_POS_ENERGY_MIN_GEV = 0.2;

// Optional common multiplicative correction applied to both photon energies.
constexpr bool ENABLE_MGG_LINEAR_ENERGY_CORRECTION = false;
constexpr double MGG_LINEAR_SLOPE = -40.0;
constexpr double MGG_LINEAR_PIVOT_GEV = 0.135;
constexpr bool MGG_LINEAR_USE_INVERSE = false;
constexpr double MGG_LINEAR_INVERSE_DENOMINATOR_MIN = 1e-9;
constexpr double MGG_LINEAR_FACTOR_MIN = 0.85;
constexpr double MGG_LINEAR_FACTOR_MAX = 1.15;

// Optional p_e scaling used only by missing-mass reconstruction.
constexpr bool ENABLE_ELECTRON_MOMENTUM_SCALING = false;

}  // namespace SmearingResponseConfig

#endif  // NPS_SMEARING_RESPONSE_CONFIG_H
