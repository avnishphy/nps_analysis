#ifndef NPS_DUMMY_SUB_H
#define NPS_DUMMY_SUB_H

// Dummy-window subtraction adapted from elastic_analysis/yield_plots.cxx.
// Header-only C++11, usable directly or through ROOT.gInterpreter.Declare.
//
// Inputs are the already selected physics-tree ytar (= H.gtr.y, cm) and
// pi0_weight. Do not reapply the old elastic W/PID cuts to these pi0 events.
// Use the SAME physics selection for LH2 and dummy, including any mass cut.
//
// For each target separately, pool the effective exposure of accepted runs:
//   Qeff = sum_r [HEL_charge_after_cut_uC / 1000 * NewGen_EDTM_livetime
//                 * HMS_tracking_eff * HMS_hodo_3of4_eff / ps_value].
// ps_value is the decoded prescale FACTOR, not the ps4/ps6 register setting.
// Add each included run once, even if it contributes zero selected events.
// Then, for each bin (or an integrated selection):
//   Y = sum_LH2(w)/Qeff_LH2
//       - sum_dummy_up(w)/(8.467*Qeff_dummy)
//       - sum_dummy_down(w)/(4.256*Qeff_dummy).
// This is the legacy effective-charge prescription; all charges here are mC.
// The same full dummy exposure normalizes BOTH windows (no factor of two).
//
// Example:
//   nps::dummy::WeightedSum lh2, up, down;
//   lh2.add(pi0_weight);  // selected LH2 event
//   // Route dummy pi0_weight into up/down using window(ytar).
//   auto result = nps::dummy::subtract(lh2, up, down, q_lh2, q_dummy);
// Alternatively fill a histogram with production_event_weight(...) and
// negative dummy_subtraction_weight(...); call Sumw2 BEFORE filling ROOT TH1.
//
// Do NOT multiply these weights by the combiner's per-run `scale`: that would
// normalize twice. That script currently sums per-run yields, not pooled
// yields. Adapting its normalization/branch schema is a separate integration.
// Keep signed weights and negative subtracted bins. Errors below are sum(w^2)
// statistics with fixed pi0 weights, exposures and material ratios; correlated
// fit, charge, efficiency, livetime and thickness uncertainties are NOT included.

#include <cmath>
#include <stdexcept>
#include <string>

namespace nps {
namespace dummy {

namespace detail {
inline void finite(double value, const char* name) {
    if (!std::isfinite(value))
        throw std::invalid_argument(std::string("nps_dummy_sub: nonfinite ") + name);
}
inline void positive(double value, const char* name) {
    finite(value, name);
    if (value <= 0.0)
        throw std::invalid_argument(std::string("nps_dummy_sub: nonpositive ") + name);
}
} // namespace detail

struct Config {
    // Dummy / LH2 aluminum thickness ratios, NOT multiplicative corrections.
    // Legacy downstream uses the tip (4.256); wall alternative was 3.711.
    // These legacy values are configurable, not a new geometry calibration.
    double upstream_ratio = 8.467;
    double downstream_ratio = 4.256;
    double ytar_split_cm = 0.0;

    void validate() const {
        detail::positive(upstream_ratio, "upstream_ratio");
        detail::positive(downstream_ratio, "downstream_ratio");
        detail::finite(ytar_split_cm, "ytar_split_cm");
    }
};

enum class Window { Unassigned = 0, Upstream = -1, Downstream = 1 };

inline Window window(double ytar_cm, const Config& config = Config()) {
    config.validate();
    detail::finite(ytar_cm, "ytar_cm");
    // Follow the old event loop, whose histogram titles have reversed signs.
    if (ytar_cm < config.ytar_split_cm) return Window::Upstream;
    if (ytar_cm > config.ytar_split_cm) return Window::Downstream;
    return Window::Unassigned; // Exact split excluded, as in the old loop.
}

inline double material_weight(double ytar_cm, const Config& config = Config()) {
    const Window side = window(ytar_cm, config);
    if (side == Window::Upstream) return 1.0 / config.upstream_ratio;
    if (side == Window::Downstream) return 1.0 / config.downstream_ratio;
    return 0.0;
}

// Sum this result across accepted runs of ONE target. No PID factor is added:
// the current combine stage uses tracking * hodo efficiency only.
inline double effective_charge_mC(double charge_uC, double ps_value,
                                  double livetime, double tracking_eff,
                                  double hodo_3of4_eff) {
    detail::positive(charge_uC, "charge_uC");
    detail::positive(ps_value, "ps_value");
    detail::positive(livetime, "livetime");
    detail::positive(tracking_eff, "tracking_eff");
    detail::positive(hodo_3of4_eff, "hodo_3of4_eff");
    const double result = (charge_uC / 1000.0) * livetime * tracking_eff
                          * hodo_3of4_eff / ps_value;
    detail::positive(result, "effective_charge_mC");
    return result;
}

inline double production_event_weight(double pi0_weight, double lh2_charge_mC) {
    detail::finite(pi0_weight, "pi0_weight");
    detail::positive(lh2_charge_mC, "lh2_charge_mC");
    const double result = pi0_weight / lh2_charge_mC;
    detail::finite(result, "production_event_weight");
    return result;
}

// Signed contribution to LH2-minus-dummy, not the positive background yield.
inline double dummy_subtraction_weight(double pi0_weight, double ytar_cm,
                                        double dummy_charge_mC,
                                        const Config& config = Config()) {
    const double result = -production_event_weight(pi0_weight, dummy_charge_mC)
                          * material_weight(ytar_cm, config);
    detail::finite(result, "dummy_subtraction_weight");
    return result;
}

// Raw weighted count and its statistical variance. Can represent one bin or
// a whole selection. When importing histograms use sumw2 = bin_error^2;
// never substitute sqrt(sumw) for a weighted statistical error.
struct WeightedSum {
    double sumw = 0.0;
    double sumw2 = 0.0;

    void validate() const {
        detail::finite(sumw, "sumw");
        detail::finite(sumw2, "sumw2");
        if (sumw2 < 0.0)
            throw std::invalid_argument("nps_dummy_sub: negative sumw2");
    }

    void add(double weight) {
        detail::finite(weight, "weight");
        WeightedSum next = *this;
        next.sumw += weight;
        next.sumw2 += weight * weight;
        next.validate();
        *this = next;
    }

    double stat_error() const {
        validate();
        return std::sqrt(sumw2);
    }
};

inline WeightedSum scaled(const WeightedSum& input, double factor) {
    input.validate();
    detail::finite(factor, "scale factor");
    WeightedSum output;
    output.sumw = input.sumw * factor;
    output.sumw2 = input.sumw2 * factor * factor;
    output.validate();
    return output;
}

struct Result {
    WeightedSum lh2;
    WeightedSum upstream;
    WeightedSum downstream;
    WeightedSum background;
    WeightedSum subtracted;
};

// Independent LH2 and disjoint upstream/downstream dummy event samples.
// All three inputs must be raw pi0-weight sums, not charge-normalized yields.
inline Result subtract(const WeightedSum& lh2, const WeightedSum& upstream,
                       const WeightedSum& downstream, double lh2_charge_mC,
                       double dummy_charge_mC, const Config& config = Config()) {
    config.validate();
    detail::positive(lh2_charge_mC, "lh2_charge_mC");
    detail::positive(dummy_charge_mC, "dummy_charge_mC");
    Result result;
    result.lh2 = scaled(lh2, 1.0 / lh2_charge_mC);
    result.upstream = scaled(upstream, (1.0 / dummy_charge_mC) / config.upstream_ratio);
    result.downstream = scaled(downstream, (1.0 / dummy_charge_mC) / config.downstream_ratio);
    result.background.sumw = result.upstream.sumw + result.downstream.sumw;
    result.background.sumw2 = result.upstream.sumw2 + result.downstream.sumw2;
    result.subtracted.sumw = result.lh2.sumw - result.background.sumw;
    result.subtracted.sumw2 = result.lh2.sumw2 + result.background.sumw2;
    result.background.validate();
    result.subtracted.validate();
    return result;
}

} // namespace dummy
} // namespace nps

#endif // NPS_DUMMY_SUB_H
