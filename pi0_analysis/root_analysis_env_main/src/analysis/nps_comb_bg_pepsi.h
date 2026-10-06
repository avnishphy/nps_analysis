// nps_comb_bg_pepsi.h
//
// PEPSI-motivated combinatorial-background subtraction for the two-photon
// invariant-mass spectrum.
//
// Physics motivation
// ------------------
// The PEPSI study supplied with this analysis (Peter Boseted elog: https://hallcweb.jlab.org/elogs/NPS-RG1a-Analysis/34) shows that wrong photon pairs
// (photons originating from two different pi0 mesons in the same event) form
// a broad continuum which is approximately flat/slowly varying at low mass
// and then turns down rapidly near the pi0 mass.  A free high-order polynomial
// does not encode that behavior: when it is fitted in two disconnected
// sidebands it can oscillate under the excluded peak and can become negative.
//
// No numerical PEPSI histogram was supplied, so this file does NOT claim to be
// a direct Monte-Carlo template fit.  It implements a positive analytic proxy
// for the PEPSI shape, the three-parameter Fermi (logistic) turn-off
//
//                 A
//   B(m) = ------------------- ,
//          1 + exp((m-mt)/w)
//
// where A is the low-mass continuum level, mt is the turn-off position, and
// w controls its width.  The formula is positive for A >= 0 and has only three
// parameters, making it substantially more stable than a quartic polynomial
// for per-run fits.  Once a machine-readable PEPSI template is available, it
// should supersede this analytic proxy (or be used for a closure test).
//
// Interface compatibility
// -----------------------
// The public BGSubtractionResult and FitCombinatorialBGAndSubtract(...) API are
// intentionally identical to nps_comb_bg.h.  The old `poly_order` argument is
// retained as `legacy_poly_order` solely so nps_analysis_main.C needs only an
// include change.  It does not alter the Fermi model; it remains in output file
// names to preserve the existing plot/output discovery contract.
//
// Statistical notes
// -----------------
// * The input has already had accidental coincidences subtracted, so its bins
//   can be negative.  The sideband fit therefore uses the supplied Gaussian
//   bin uncertainties rather than a Poisson likelihood.
// * The fitted-parameter covariance is propagated point-by-point into the
//   subtracted histogram.  Those prediction errors are correlated between
//   mass bins; TH1D stores only their diagonal part.  An integrated-yield
//   uncertainty must use the full covariance, not IntegralAndError alone.
// * The pull panel contains sideband points only.  The pi0 signal window is not
//   a background goodness-of-fit region.

#ifndef NPS_COMB_BG_PEPSI_H
#define NPS_COMB_BG_PEPSI_H

#include <TH1D.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TBox.h>
#include <TPaveText.h>
#include <TMatrixDSym.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TDecompChol.h>
#include <TFile.h>
#include <TNamed.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TStyle.h>
#include <RVersion.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

namespace nps {

// Kept byte-for-byte compatible at the field/interface level with the result
// type in nps_comb_bg.h, so no downstream analysis code needs to change.
struct BGSubtractionResult {
    TH1D* h_final = nullptr;  // caller owns this detached histogram
    double chi2_ndf = -1.0;
    double mu_MeV = 0.0;
    double sigma_MeV = 0.0;
    double signal_counts = 0.0;
    bool success = false;
    int minimizer_status = -1;
    int covariance_status = -1;
    int attempts = 0;
    double chi2 = 0.0;
    int ndf = 0;
    double edm = 0.0;
    unsigned int calls = 0;
    bool at_boundary = false;
    bool zero_background = false;
    bool active_shape_limits = false;
    double amplitude = 0.0;
    double background_integral = 0.0;
    double zero_chi2 = 0.0;
    double boundary_score_upper = 0.0;
    std::string failure_reason = "invalid_input";
};

// For chi2(A,theta)=chi2(0)-2*A*q(theta)+A*A*d(theta), A>=0,
// zero is the global optimum iff q(theta)<=0 for every allowed shape.
// Bound q on rectangles using exact corner extrema of the logistic function.
// Refine ambiguous rectangles deterministically; never infer zero from a
// small fitted amplitude or a failed covariance. The tolerance bounds only
// floating-point summation error, in units of q, not amplitude.
inline bool CertifyPepsiZero(const TGraphErrors& graph, double tlo, double thi,
        double wlo, double whi, double& upper, double& chi0) {
    struct Cell { double tl,th,wl,wh; int depth; };
    std::vector<Cell> cells{{tlo,thi,wlo,whi,0}};
    double scale=0; chi0=0; upper=-std::numeric_limits<double>::infinity();
    for (int i=0;i<graph.GetN();++i) {
        double x,y; graph.GetPoint(i,x,y); const double e=graph.GetErrorY(i);
        scale+=std::abs(y/(e*e)); chi0+=(y/e)*(y/e);
    }
    const double tol=32*std::numeric_limits<double>::epsilon()*graph.GetN()*scale;
    const double objective_precision=64*std::numeric_limits<double>::epsilon()*graph.GetN()*std::max(1.,chi0);
    int visited=0;
    while (!cells.empty()) {
        const auto cell=cells.back(); cells.pop_back();
        double bound=0, score=0, dlow=0, dmid=0;
        for(int i=0;i<graph.GetN();++i) {
            double x,y; graph.GetPoint(i,x,y); const double e=graph.GetErrorY(i);
            const double c=y/(e*e);
            double low=1,high=0;
            for(double t:{cell.tl,cell.th}) for(double w:{cell.wl,cell.wh}) {
                const double f=1/(1+std::exp((x-t)/w));
                low=std::min(low,f); high=std::max(high,f);
            }
            bound+=c*(c>=0?high:low);
            const double mid=1/(1+std::exp((x-(cell.tl+cell.th)/2)/std::sqrt(cell.wl*cell.wh)));
            score+=c*mid; dlow+=low*low/(e*e); dmid+=mid*mid/(e*e);
        }
        if (bound<=tol || (dlow>0 && bound*bound/dlow<=objective_precision)) {
            upper=std::max(upper,bound); continue;
        }
        if ((score>tol && score*score/dmid>objective_precision) || ++visited>65536 || cell.depth>=40) {
            upper=bound; return false;
        }
        if ((cell.th-cell.tl)/(thi-tlo) >= std::log(cell.wh/cell.wl)/std::log(whi/wlo)) {
            const double m=(cell.tl+cell.th)/2;
            cells.push_back({cell.tl,m,cell.wl,cell.wh,cell.depth+1});
            cells.push_back({m,cell.th,cell.wl,cell.wh,cell.depth+1});
        } else {
            const double m=std::sqrt(cell.wl*cell.wh);
            cells.push_back({cell.tl,cell.th,cell.wl,m,cell.depth+1});
            cells.push_back({cell.tl,cell.th,m,cell.wh,cell.depth+1});
        }
    }
    return true;
}

inline bool ValidPepsiFit(const TFitResultPtr& fit) {
    if (!fit.Get() || !fit->IsValid() || int(fit) != 0 ||
        fit->CovMatrixStatus() != 3 || fit->Ndf() <= 0 ||
        !std::isfinite(fit->Chi2()) || !std::isfinite(fit->Edm())) return false;
    const auto cov = fit->GetCovarianceMatrix();
    for (unsigned int i=0; i<fit->NPar(); ++i) {
        if (!std::isfinite(fit->Parameter(i)) || !std::isfinite(fit->ParError(i))) return false;
        for (unsigned int j=0; j<fit->NPar(); ++j)
            if (!std::isfinite(cov(i,j))) return false;
    }
    std::vector<unsigned int> free;
    for(unsigned int i=0;i<fit->NPar();++i) if(!fit->IsParameterFixed(i)) free.push_back(i);
    if(free.empty()) return false;
    TMatrixDSym free_cov(free.size());
    for(unsigned int i=0;i<free.size();++i) for(unsigned int j=0;j<free.size();++j)
        free_cov(i,j)=cov(free[i],free[j]);
    TDecompChol chol(free_cov);
    return chol.Decompose();
}

inline bool InPepsiSideband(double x,
                            double left_lo, double left_hi,
                            double right_lo, double right_hi)
{
    return (x >= left_lo && x <= left_hi) ||
           (x >= right_lo && x <= right_hi);
}

// Convert the selected histogram bins into a graph.  TGraphErrors is used
// because the accidental-subtracted spectrum is not a Poisson-distributed
// non-negative histogram.  Existing Sumw2 uncertainties are retained.
inline TGraphErrors* MakePepsiSidebandGraph(const TH1D* h,
                                             double left_lo, double left_hi,
                                             double right_lo, double right_hi)
{
    if (!h) return nullptr;

    std::vector<double> xs, ys, exs, eys;
    xs.reserve(h->GetNbinsX());
    ys.reserve(h->GetNbinsX());
    exs.reserve(h->GetNbinsX());
    eys.reserve(h->GetNbinsX());

    for (int bin = 1; bin <= h->GetNbinsX(); ++bin) {
        const double x = h->GetXaxis()->GetBinCenter(bin);
        if (!InPepsiSideband(x, left_lo, left_hi, right_lo, right_hi)) continue;

        const double y = h->GetBinContent(bin);
        double error = h->GetBinError(bin);
        // This fallback is used only when the source histogram has no usable
        // uncertainty.  For a properly Sumw2-enabled accidental subtraction,
        // the first branch should normally supply the error.
        if (!(error > 0.0) || !std::isfinite(error))
            error = (y > 0.0) ? std::sqrt(y) : 1.0;

        xs.push_back(x);
        ys.push_back(y);
        exs.push_back(0.0);
        eys.push_back(error);
    }

    if (xs.empty()) return nullptr;
    auto* graph = new TGraphErrors(static_cast<int>(xs.size()));
    graph->SetName("g_pepsi_sideband_points");
    for (std::size_t i = 0; i < xs.size(); ++i) {
        graph->SetPoint(static_cast<int>(i), xs[i], ys[i]);
        graph->SetPointError(static_cast<int>(i), exs[i], eys[i]);
    }
    return graph;
}

// Propagate the full 3x3 fit covariance to B(x).  For q=1/(1+exp(z)),
// dB/dA=q, dB/dmt=A*q*(1-q)/w, and
// dB/dw=A*q*(1-q)*(x-mt)/w^2.
inline double PepsiFermiErrorAtX(const TF1& function,
                                 const TMatrixDSym& covariance,
                                 double x)
{
    if (covariance.GetNrows() < 3) return 0.0;

    const double amplitude = function.GetParameter(0);
    const double turn_mass = function.GetParameter(1);
    const double width = function.GetParameter(2);
    if (!(width > 0.0)) return 0.0;

    const double z = (x - turn_mass) / width;
    double q = 0.0;
    if (z > 50.0) q = std::exp(-z);
    else if (z < -50.0) q = 1.0;
    else q = 1.0 / (1.0 + std::exp(z));

    const double common = amplitude * q * (1.0 - q);
    const double gradient[3] = {
        q,
        common / width,
        common * (x - turn_mass) / (width * width)
    };

    double variance = 0.0;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            variance += gradient[i] * covariance(i, j) * gradient[j];

    // Tiny negative values can appear from floating-point cancellation.
    return (variance > 0.0 && std::isfinite(variance))
        ? std::sqrt(variance) : 0.0;
}

inline BGSubtractionResult FitCombinatorialBGAndSubtract(
    TH1D* h_coin_bgsub,
    const char* outDir = "",
    int run = -1,
    int legacy_poly_order = 2,
    double left_lo = 0.01, double left_hi = 0.10,
    double right_lo = 0.15, double right_hi = 0.40,
    bool draw = true)
{
    BGSubtractionResult result;

    if (!h_coin_bgsub) {
        std::cerr << "[nps::PEPSICombBG] ERROR: null input histogram\n";
        return result;
    }
    if (!(left_lo < left_hi && left_hi < right_lo && right_lo < right_hi)) {
        std::cerr << "[nps::PEPSICombBG] ERROR: sideband windows must satisfy "
                  << "left_lo < left_hi < right_lo < right_hi\n";
        return result;
    }

#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 0, 0)
    if (h_coin_bgsub->GetSumw2N() == 0) h_coin_bgsub->Sumw2();
#else
    h_coin_bgsub->Sumw2();
#endif

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    std::unique_ptr<TGraphErrors> sidebands(
        MakePepsiSidebandGraph(h_coin_bgsub,
                               left_lo, left_hi, right_lo, right_hi));
    if (!sidebands) {
        std::cerr << "[nps::PEPSICombBG] ERROR: no bins in the sideband windows\n";
        return result;
    }

    sidebands->SetMarkerStyle(21);
    sidebands->SetMarkerSize(0.85);
    sidebands->SetMarkerColor(kAzure + 7);

    const double x_min = h_coin_bgsub->GetXaxis()->GetXmin();
    const double x_max = h_coin_bgsub->GetXaxis()->GetXmax();

    // Use the positive mean of the left sideband as the plateau seed.  The
    // fitted amplitude remains free; this seed only helps Minuit converge.
    double left_sum = 0.0;
    int left_count = 0;
    double largest_sideband_value = 0.0;
    for (int i = 0; i < sidebands->GetN(); ++i) {
        double x = 0.0, y = 0.0;
        sidebands->GetPoint(i, x, y);
        largest_sideband_value = std::max(largest_sideband_value, y);
        if (x >= left_lo && x <= left_hi && y > 0.0) {
            left_sum += y;
            ++left_count;
        }
    }
    const double amplitude_seed = (left_count > 0)
        ? left_sum / static_cast<double>(left_count)
        : std::max(1.0, largest_sideband_value);

    const TString function_name = TString::Format("f_bg_pepsi_run%d", run);
    // A formula-based TF1 is used instead of a C++ callback so ROOT can safely
    // serialize and reopen `fitted_bg` in the diagnostic ROOT file.
    std::unique_ptr<TF1> f_bg(new TF1(function_name,
        "[0]/(1.0+exp((x-[1])/[2]))", x_min, x_max));
    f_bg->SetParNames("A", "m_turn", "width");
    f_bg->SetParameters(std::max(1.0, amplitude_seed),
                        std::max(left_hi, right_lo - 0.005), 0.020);

    // Bounds encode only broad physical/numerical requirements.  In
    // particular, mt is allowed inside the excluded peak interval; forcing it
    // into either sideband would bias the background interpolation.
    const double amplitude_upper = std::max(10.0, 20.0 * std::max(1.0, largest_sideband_value));
    const double turn_lower = std::max(x_min, left_hi);
    const double turn_upper = std::min(x_max, std::max(right_lo + 0.07, 0.20));
    const double bin_width = h_coin_bgsub->GetXaxis()->GetBinWidth(1);
    const double width_lower = std::max(0.001, 0.25 * bin_width);
    const double width_upper = std::max(width_lower * 2.0,
                                        std::min(0.10, right_hi - left_lo));
    f_bg->SetParLimits(0, 0.0, amplitude_upper);
    f_bg->SetParLimits(1, turn_lower, turn_upper);
    f_bg->SetParLimits(2, width_lower, width_upper);

    // R: respect TF1 range, Q: quiet, S: return covariance, N: do not attach
    // the temporary fit function to the graph.
    result.zero_background = sidebands->GetN()>3 && CertifyPepsiZero(*sidebands,
        turn_lower,turn_upper,width_lower,width_upper,
        result.boundary_score_upper,result.zero_chi2);
    // Histogram-only replicas need no fit of undefined shape parameters once
    // the global zero prediction has been certified. Nominal diagnostic runs
    // retain the ordinary-fit status for the historical comparison.
    const bool skip_undefined_fit=result.zero_background && !draw && (!outDir || !*outDir);
    TFitResultPtr fit_result(-1);
    if (!skip_undefined_fit) fit_result=sidebands->Fit(f_bg.get(), "RQSN");
    result.attempts = skip_undefined_fit?0:1;
    auto log_attempt = [&]() {
        std::cout << "[COMB_FIT] run=" << run << " attempt=" << result.attempts
                  << " status=" << int(fit_result)
                  << " valid=" << (fit_result.Get() && fit_result->IsValid())
                  << " covariance=" << (fit_result.Get() ? fit_result->CovMatrixStatus() : -1)
                  << " edm=" << (fit_result.Get() ? fit_result->Edm() : -1) << '\n';
        if (outDir && std::string(outDir).size()) {
            gSystem->mkdir(outDir,true);
            const std::string path=std::string(outDir)+"/combinatorial_fit_run"+std::to_string(run)+".csv";
            std::ofstream out(path, result.attempts==1 ? std::ios::out : std::ios::app);
            if (result.attempts==1) out << "run,attempt,accepted,minimizer_status,covariance_status,chi2,ndf,edm,ncalls,A,A_error,m_turn,m_turn_error,width,width_error\n";
            out << std::setprecision(17) << run << ',' << result.attempts << ','
                << ValidPepsiFit(fit_result) << ',' << int(fit_result) << ','
                << (fit_result.Get()?fit_result->CovMatrixStatus():-1) << ','
                << (fit_result.Get()?fit_result->Chi2():-1) << ','
                << (fit_result.Get()?fit_result->Ndf():-1) << ','
                << (fit_result.Get()?fit_result->Edm():-1) << ','
                << (fit_result.Get()?fit_result->NCalls():0);
            for (int i=0;i<3;++i) out << ',' << f_bg->GetParameter(i) << ',' << f_bg->GetParError(i);
            out << '\n';
        }
    };
    if (!skip_undefined_fit) log_attempt();
    // One deterministic restart from the current point; same objective, data,
    // bounds and options. Never resample/retry until a replica passes.
    if (!skip_undefined_fit && !ValidPepsiFit(fit_result)) {
        fit_result = sidebands->Fit(f_bg.get(), "RQSN");
        ++result.attempts;
        log_attempt();
    }
    // A local Minuit minimum very near zero can miss a positive-amplitude
    // branch. Profile A analytically on a deterministic shape grid, then
    // restart once from its best point if it improves the current objective.
    // This is an optimizer seed/check, never a significance or amplitude cut.
    bool active_shape_check_failed=false;
    const double objective_roundoff=64*std::numeric_limits<double>::epsilon()*sidebands->GetN()*std::max(1.,result.zero_chi2);
    const bool healthy_positive_fit=ValidPepsiFit(fit_result) && result.zero_chi2-fit_result->Chi2()>objective_roundoff;
    if (!result.zero_background && !healthy_positive_fit && sidebands->GetN()>3) {
        double best=result.zero_chi2, best_a=0, best_t=turn_lower, best_w=width_lower;
        for(int it=0;it<=32;++it) for(int iw=0;iw<=32;++iw) {
            const double t=turn_lower+(turn_upper-turn_lower)*it/32.;
            const double w=width_lower*std::pow(width_upper/width_lower,iw/32.);
            double q=0,d=0;
            for(int i=0;i<sidebands->GetN();++i) {
                double x,y; sidebands->GetPoint(i,x,y); const double e=sidebands->GetErrorY(i);
                const double f=1/(1+std::exp((x-t)/w)); q+=y*f/(e*e); d+=f*f/(e*e);
            }
            const double a=d>0?std::min(amplitude_upper,std::max(0.,q/d)):0.;
            const double objective=result.zero_chi2-2*a*q+a*a*d;
            if(objective<best) { best=objective; best_a=a; best_t=t; best_w=w; }
        }
        const double precision=64*std::numeric_limits<double>::epsilon()*sidebands->GetN()*std::max(1.,result.zero_chi2);
        if (!fit_result.Get() || fit_result->Chi2()>best+precision) {
            f_bg->SetParameters(best_a,best_t,best_w);
            fit_result=sidebands->Fit(f_bg.get(),"RQSN"); ++result.attempts; log_attempt();
            if (!ValidPepsiFit(fit_result)) {
                fit_result=sidebands->Fit(f_bg.get(),"RQSN"); ++result.attempts; log_attempt();
            }
        }
        // Test active shape limits with inward KKT gradients. Covariance is
        // required only on the remaining free parameters. Recheck the fixed
        // gradients after refitting, because the free optimum may move.
        if (!ValidPepsiFit(fit_result) && best_a>0 && best_a<amplitude_upper &&
            (best_t==turn_lower || best_t==turn_upper || best_w==width_lower || best_w==width_upper)) {
            double gt=0,gw=0,st=0,sw=0;
            for(int i=0;i<sidebands->GetN();++i) {
                double x,y; sidebands->GetPoint(i,x,y); const double e=sidebands->GetErrorY(i);
                const double f=1/(1+std::exp((x-best_t)/best_w));
                const double common=-2*best_a*(y-best_a*f)*f*(1-f)/(e*e);
                const double dt=common/best_w, dw=dt*(x-best_t)/best_w;
                gt+=dt; gw+=dw; st+=std::abs(dt); sw+=std::abs(dw);
            }
            const double rounding=64*std::numeric_limits<double>::epsilon()*sidebands->GetN();
            const bool kkt_t=(best_t==turn_lower || best_t==turn_upper) && (best_t==turn_lower?gt:-gt)>=-rounding*st;
            const bool kkt_w=(best_w==width_lower || best_w==width_upper) && (best_w==width_lower?gw:-gw)>=-rounding*sw;
            if(kkt_t || kkt_w) {
                const double previous=fit_result.Get()?fit_result->Chi2():std::numeric_limits<double>::infinity();
                f_bg->SetParameters(best_a,best_t,best_w);
                if(kkt_t) f_bg->FixParameter(1,best_t);
                if(kkt_w) f_bg->FixParameter(2,best_w);
                fit_result=sidebands->Fit(f_bg.get(),"RQSN"); ++result.attempts; log_attempt();
                gt=gw=st=sw=0;
                const double a=f_bg->GetParameter(0),t=f_bg->GetParameter(1),w=f_bg->GetParameter(2);
                for(int i=0;i<sidebands->GetN();++i) {
                    double x,y; sidebands->GetPoint(i,x,y); const double e=sidebands->GetErrorY(i);
                    const double f=1/(1+std::exp((x-t)/w));
                    const double dt=-2*a*(y-a*f)*f*(1-f)/(e*e*w),dw=dt*(x-t)/w;
                    gt+=dt;gw+=dw;st+=std::abs(dt);sw+=std::abs(dw);
                }
                result.active_shape_limits=ValidPepsiFit(fit_result) && fit_result->Chi2()<=previous+precision &&
                    (!kkt_t || (best_t==turn_lower?gt:-gt)>=-rounding*st) &&
                    (!kkt_w || (best_w==width_lower?gw:-gw)>=-rounding*sw);
                if(!result.active_shape_limits) {
                    // Fixed limits that fail the constrained-optimum check
                    // cannot become an accepted covariance rescue.
                    active_shape_check_failed=true;
                }
            }
        }
    }
    const bool fit_valid = result.zero_background || (!active_shape_check_failed && ValidPepsiFit(fit_result) &&
        f_bg->GetParameter(0)>0 && fit_result->Chi2()<result.zero_chi2);
    result.success = fit_valid;
    result.minimizer_status = int(fit_result);
    result.covariance_status = fit_result.Get() ? fit_result->CovMatrixStatus() : -1;
    result.failure_reason = active_shape_check_failed?"active_shape_limit_check_failed":fit_valid ? "" :
        (!fit_result.Get() ? "missing_fit_result" :
         fit_result->Ndf() <= 0 ? "insufficient_degrees_of_freedom" :
         int(fit_result) == 1 ? "covariance_forced_positive_definite" :
         int(fit_result) != 0 ? "minimizer_failure" :
         fit_result->CovMatrixStatus() != 3 ? "inaccurate_covariance" :
         "invalid_minimum_or_nonpositive_covariance");
    if (fit_result.Get()) {
        result.chi2 = fit_result->Chi2(); result.ndf = fit_result->Ndf();
        result.edm = fit_result->Edm(); result.calls = fit_result->NCalls();
    }
    if (result.zero_background) {
        // No Hessian exists for the absent shape. Preserve Minuit diagnostics
        // above, but publish exactly zero prediction and the nested objective.
        f_bg->SetParameters(0.0,(turn_lower+turn_upper)/2,std::sqrt(width_lower*width_upper));
        for(int i=0;i<3;++i) f_bg->SetParError(i,0.0);
        result.chi2=result.zero_chi2; result.ndf=sidebands->GetN();
        result.at_boundary=true;
    }
    result.amplitude=f_bg->GetParameter(0);
    for (int i=0; i<3; ++i) {
        double lo=0, hi=0; f_bg->GetParLimits(i,lo,hi);
        const double p=f_bg->GetParameter(i);
        result.at_boundary |= std::min(p-lo,hi-p) < 1e-6*(hi-lo);
    }
    if (!fit_valid) {
        std::cerr << "[nps::PEPSICombBG] FAILED run=" << run
                  << " status=" << result.minimizer_status
                  << " covariance=" << result.covariance_status << '\n';
        if (outDir && std::string(outDir).size()) {
            TFile diagnostics((std::string(outDir)+"/failed_combinatorial_fit_run"+
                std::to_string(run)+".root").c_str(),"RECREATE");
            h_coin_bgsub->Write("fit_input");
            sidebands->Write("sideband_points");
            f_bg->Write("rejected_background_parameters");
            if (fit_result.Get()) fit_result->Write("rejected_fit_result");
        }
        return result; // No subtraction histogram or signal weights from a bad fit.
    }

    TMatrixDSym covariance(3);
    covariance.Zero();
    if (fit_valid && !result.zero_background) {
        covariance = fit_result->GetCovarianceMatrix();
    }

    const double chi2 = result.chi2;
    const int ndf = result.ndf;
    result.chi2_ndf = (ndf > 0) ? chi2 / static_cast<double>(ndf) : -1.0;

    std::cout << "[nps::PEPSICombBG] run=" << run
              << " model=fermi_logistic"
              << " legacy_order_tag=" << legacy_poly_order
              << " fit_valid=" << (fit_valid ? "yes" : "no")
              << " chi2/ndf=" << result.chi2_ndf << "\n"
              << "  A=" << f_bg->GetParameter(0)
              << " +/- " << f_bg->GetParError(0) << " events/bin\n"
              << "  m_turn=" << 1000.0 * f_bg->GetParameter(1)
              << " +/- " << 1000.0 * f_bg->GetParError(1) << " MeV\n"
              << "  width=" << 1000.0 * f_bg->GetParameter(2)
              << " +/- " << 1000.0 * f_bg->GetParError(2) << " MeV\n";

    // The returned histogram is detached from every TFile.  The caller owns
    // it, matching the old helper's ownership contract.
    const TString final_name = TString::Format(
        "%s_bgsub_pepsi_run%d", h_coin_bgsub->GetName(), run);
    auto* h_final = static_cast<TH1D*>(h_coin_bgsub->Clone(final_name));
    h_final->SetDirectory(nullptr);
    h_final->SetTitle(TString::Format(
        "%s (PEPSI-motivated combinatorial BG subtracted)",
        h_coin_bgsub->GetTitle()));
#if ROOT_VERSION_CODE >= ROOT_VERSION(6, 0, 0)
    if (h_final->GetSumw2N() == 0) h_final->Sumw2();
#else
    h_final->Sumw2();
#endif

    for (int bin = 1; bin <= h_coin_bgsub->GetNbinsX(); ++bin) {
        const double x = h_coin_bgsub->GetXaxis()->GetBinCenter(bin);
        const double data = h_coin_bgsub->GetBinContent(bin);
        double data_error = h_coin_bgsub->GetBinError(bin);
        if (!(data_error > 0.0) || !std::isfinite(data_error))
            data_error = (data > 0.0) ? std::sqrt(data) : 1.0;

        const double background = f_bg->Eval(x);
        result.background_integral += background;
        const double background_error = PepsiFermiErrorAtX(*f_bg, covariance, x);
        h_final->SetBinContent(bin, data - background);
        h_final->SetBinError(bin,
            std::hypot(data_error, background_error));
    }

    // Preserve the old helper's downstream signal-summary procedure so a
    // comparison changes only the combinatorial-background model.  This
    // Gaussian is diagnostic; h_final itself is the subtraction product used
    // by the main analysis.
    if (!draw && (!outDir || !*outDir)) {
        result.h_final=h_final;
        result.signal_counts=h_final->Integral();
        return result;
    }
    const double exclusion_lo = left_hi;
    const double exclusion_hi = right_lo;
    const int maximum_bin = h_final->GetMaximumBin();
    const double mu_seed = h_final->GetBinCenter(maximum_bin);
    const double signal_amplitude_seed = h_final->GetBinContent(maximum_bin);
    double sigma_seed = h_final->GetRMS();
    if (!(sigma_seed > 0.0) || !std::isfinite(sigma_seed))
        sigma_seed = (exclusion_hi - exclusion_lo) / 6.0;
    sigma_seed = std::min(sigma_seed, 0.020);

    const double signal_lo = std::max(x_min, mu_seed - 3.0 * sigma_seed);
    const double signal_hi = std::min(x_max, mu_seed + 3.0 * sigma_seed);
    TF1 signal_fit(TString::Format("f_sig_pepsi_run%d", run),
                   "gaus", signal_lo, signal_hi);
    signal_fit.SetParameters(signal_amplitude_seed, mu_seed, sigma_seed);
    TFitResultPtr signal_result = h_final->Fit(&signal_fit, "RQSN");
    (void)signal_result;

    const double signal_mu = signal_fit.GetParameter(1);
    const double signal_sigma = std::abs(signal_fit.GetParameter(2));
    result.mu_MeV = signal_mu * 1000.0;
    result.sigma_MeV = signal_sigma * 1000.0;

    const int signal_bin_lo = h_final->FindBin(signal_lo);
    const int signal_bin_hi = h_final->FindBin(signal_hi);
    result.signal_counts = h_final->Integral(signal_bin_lo, signal_bin_hi);
    result.h_final = h_final;
    if (outDir && std::string(outDir).size()) {
        std::ofstream out(std::string(outDir)+"/background_classification_run"+std::to_string(run)+".csv");
        out << "run,classification,amplitude,minimizer_status,covariance_status,chi2,ndf,background_integral,zero_chi2,boundary_score_upper,reason\n"
            << std::setprecision(17) << run << ',' << (result.zero_background?"zero_background":"interior_valid")
            << ',' << result.amplitude << ',' << result.minimizer_status << ',' << result.covariance_status
            << ',' << result.chi2 << ',' << result.ndf << ',' << result.background_integral
            << ',' << result.zero_chi2 << ',' << result.boundary_score_upper << ','
            << (result.zero_background?"global_nonpositive_amplitude_score":"valid_minimum_and_covariance") << '\n';
    }

    if (draw) {
        if (outDir && std::string(outDir).size() > 0) gSystem->mkdir(outDir, true);
        const TString run_dir = (run >= 0)
            ? TString::Format("%s/run_%d", outDir, run)
            : TString::Format("%s/run_all", outDir);
        gSystem->mkdir(run_dir, true);

        const TString canvas_name = (run >= 0)
            ? TString::Format("c_combbg_pepsi_run%d", run)
            : TString("c_combbg_pepsi");
        // Declare the canvas before all primitives and let RAII destroy it
        // last.  ROOT pads retain non-owning pointers to drawn stack objects;
        // deleting the canvas early can therefore leave dangling references.
        std::unique_ptr<TCanvas> canvas(new TCanvas(canvas_name,
            "PEPSI-motivated combinatorial BG fit and subtraction", 1400, 1000));
        canvas->Divide(2, 2);

        auto style_pad = []() {
            gPad->SetLeftMargin(0.12);
            gPad->SetBottomMargin(0.12);
            gPad->SetTopMargin(0.08);
            gPad->SetRightMargin(0.05);
            gPad->SetTicks(1, 1);
        };

        // Pad 1: input spectrum, the selected sidebands, and the extrapolated
        // positive background.  The shaded gap is deliberately not fitted.
        canvas->cd(1);
        style_pad();
        h_coin_bgsub->SetLineColor(kBlack);
        h_coin_bgsub->SetLineWidth(1);
        h_coin_bgsub->SetMarkerStyle(20);
        h_coin_bgsub->SetMarkerSize(0.7);
        h_coin_bgsub->SetMinimum(0.0);
        h_coin_bgsub->SetMaximum(1.20 * std::max(1.0, h_coin_bgsub->GetMaximum()));
        h_coin_bgsub->Draw("HIST");

        TBox left_box(left_lo, gPad->GetUymin(), left_hi, gPad->GetUymax());
        left_box.SetFillColorAlpha(kCyan - 9, 0.12);
        left_box.SetLineColor(kCyan - 9);
        left_box.Draw("SAME");
        TBox right_box(right_lo, gPad->GetUymin(), right_hi, gPad->GetUymax());
        right_box.SetFillColorAlpha(kCyan - 9, 0.12);
        right_box.SetLineColor(kCyan - 9);
        right_box.Draw("SAME");

        f_bg->SetLineColor(kOrange + 7);
        f_bg->SetLineWidth(3);
        f_bg->Draw("SAME");
        sidebands->Draw("P SAME");

        TLine excluded_left(left_hi, gPad->GetUymin(), left_hi, gPad->GetUymax());
        excluded_left.SetLineStyle(2);
        excluded_left.SetLineColor(kGray + 2);
        excluded_left.Draw("SAME");
        TLine excluded_right(right_lo, gPad->GetUymin(), right_lo, gPad->GetUymax());
        excluded_right.SetLineStyle(2);
        excluded_right.SetLineColor(kGray + 2);
        excluded_right.Draw("SAME");

        TLegend legend1(0.51, 0.64, 0.93, 0.88);
        legend1.SetBorderSize(0);
        legend1.SetFillStyle(0);
        legend1.AddEntry(h_coin_bgsub, "Data after accidental subtraction", "l");
        legend1.AddEntry(sidebands.get(), "Sideband points used in fit", "p");
        legend1.AddEntry(f_bg.get(), "PEPSI-motivated Fermi background", "l");
        legend1.Draw();

        // Pad 2: direct before/after comparison.  A detached clone avoids
        // modifying the input histogram's style or ownership downstream.
        canvas->cd(2);
        style_pad();
        std::unique_ptr<TH1D> h_before(static_cast<TH1D*>(h_coin_bgsub->Clone(
            TString::Format("%s_before_pepsi_draw", h_coin_bgsub->GetName()))));
        h_before->SetDirectory(nullptr);
        h_before->SetLineColor(kGray + 2);
        h_before->SetLineWidth(2);
        h_before->SetMarkerStyle(0);
        h_before->SetMinimum(0.0);
        h_before->SetMaximum(1.20 * std::max(
            std::max(1.0, h_before->GetMaximum()), h_final->GetMaximum()));
        h_before->Draw("HIST");
        h_final->SetLineColor(kBlue + 1);
        h_final->SetLineWidth(2);
        h_final->SetMarkerStyle(0);
        h_final->Draw("HIST SAME");
        TLegend legend2(0.50, 0.70, 0.93, 0.88);
        legend2.SetBorderSize(0);
        legend2.SetFillStyle(0);
        legend2.AddEntry(h_before.get(), "Before combinatorial subtraction", "l");
        legend2.AddEntry(h_final, "After PEPSI-motivated subtraction", "l");
        legend2.Draw();

        // Pad 3: diagnostic Gaussian on the subtracted pi0 peak.  This uses
        // the same summary logic as the original helper for fair comparisons.
        canvas->cd(3);
        style_pad();
        h_final->SetMinimum(0.0);
        h_final->SetMaximum(1.20 * std::max(
            std::max(1.0, h_final->GetMaximum()),
            signal_fit.GetMaximum(signal_lo, signal_hi)));
        h_final->Draw("HIST E");
        signal_fit.SetLineColor(kRed + 1);
        signal_fit.SetLineWidth(2);
        signal_fit.Draw("SAME");
        TPaveText signal_text(0.50, 0.64, 0.91, 0.90, "NDC");
        signal_text.SetFillColor(kWhite);
        signal_text.SetBorderSize(1);
        signal_text.SetTextAlign(12);
        signal_text.SetTextFont(42);
        signal_text.SetTextSize(0.038);
        signal_text.AddText(TString::Format("Run %d", run));
        signal_text.AddText(TString::Format("#mu = %.1f #pm %.1f MeV",
            result.mu_MeV, 1000.0 * signal_fit.GetParError(1)));
        signal_text.AddText(TString::Format("#sigma = %.1f #pm %.1f MeV",
            result.sigma_MeV, 1000.0 * signal_fit.GetParError(2)));
        signal_text.AddText(TString::Format("Signal-window counts = %.1f",
            result.signal_counts));
        signal_text.Draw("SAME");

        // Pad 4: true sideband pulls.  Unlike TH1::GetRMS(), the calculation
        // below is the RMS of the pull values themselves.  Fit uncertainty is
        // not added to each denominator because the chi2 minimized above is
        // defined using the measurement errors; the fitted points are already
        // correlated through the fitted parameters.
        canvas->cd(4);
        style_pad();
        gPad->SetGridy(true);
        std::unique_ptr<TGraphErrors> pull_graph(new TGraphErrors(sidebands->GetN()));
        pull_graph->SetName("g_pepsi_sideband_pulls");
        double pull_sum = 0.0;
        double pull_square_sum = 0.0;
        int pull_count = 0;
        for (int i = 0; i < sidebands->GetN(); ++i) {
            double x = 0.0, y = 0.0;
            sidebands->GetPoint(i, x, y);
            const double error = sidebands->GetErrorY(i);
            const double pull = (error > 0.0) ? (y - f_bg->Eval(x)) / error : 0.0;
            pull_graph->SetPoint(i, x, pull);
            pull_graph->SetPointError(i, 0.0, 0.0);
            pull_sum += pull;
            pull_square_sum += pull * pull;
            ++pull_count;
        }
        const double pull_mean = (pull_count > 0) ? pull_sum / pull_count : 0.0;
        const double pull_variance = (pull_count > 0)
            ? std::max(0.0, pull_square_sum / pull_count - pull_mean * pull_mean)
            : 0.0;
        const double pull_rms = std::sqrt(pull_variance);

        auto* pull_frame = gPad->DrawFrame(x_min, -5.5, x_max, 5.5,
            ";M_{#gamma#gamma} [GeV];(data - background)/#sigma_{data}");
        pull_frame->SetTitle("Sideband pulls only");
        pull_graph->SetMarkerStyle(20);
        pull_graph->SetMarkerSize(0.75);
        pull_graph->Draw("P SAME");
        TLine zero_line(x_min, 0.0, x_max, 0.0);
        zero_line.SetLineColor(kRed + 1);
        zero_line.Draw("SAME");
        TPaveText pull_stats(0.54, 0.70, 0.92, 0.84, "NDC");
        pull_stats.SetFillColor(kWhite);
        pull_stats.SetBorderSize(1);
        pull_stats.SetTextAlign(12);
        pull_stats.SetTextFont(42);
        pull_stats.SetTextSize(0.030);
        pull_stats.AddText(TString::Format(
            "Pull mean = %.3f, RMS = %.3f", pull_mean, pull_rms));
        pull_stats.AddText(TString::Format(
            "#chi^{2}/ndf = %.3f", result.chi2_ndf));
        pull_stats.Draw("SAME");

        canvas->Update();
        // Preserve legacy filename patterns so existing PDF collection scripts
        // continue to find this diagnostic without modification.
        const TString png_path = TString::Format(
            "%s/combbg_run%d_order%d_enhanced.png",
            run_dir.Data(), run, legacy_poly_order);
        canvas->SaveAs(png_path);
        std::cout << "[nps::PEPSICombBG] wrote PNG: " << png_path << "\n";

        const TString root_path = TString::Format(
            "%s/combbg_run%d_order%d_results.root",
            run_dir.Data(), run, legacy_poly_order);
        TFile output_file(root_path, "RECREATE");
        if (!output_file.IsOpen() || output_file.IsZombie()) {
            std::cerr << "[nps::PEPSICombBG] ERROR: cannot create "
                      << root_path << "\n";
        } else {
            output_file.mkdir(TString::Format("run_%d", run));
            output_file.cd(TString::Format("run_%d", run));

            std::unique_ptr<TH1D> input_copy(static_cast<TH1D*>(
                h_coin_bgsub->Clone("h_coin_bgsub_input")));
            std::unique_ptr<TH1D> final_copy(static_cast<TH1D*>(
                h_final->Clone("h_bgsub_final")));
            // Clone() attaches histograms to the current directory by default.
            // Detach them so the unique_ptrs, rather than TFile::Close(), own
            // their lifetime; otherwise both would try to delete each clone.
            input_copy->SetDirectory(nullptr);
            final_copy->SetDirectory(nullptr);
            std::unique_ptr<TGraphErrors> sideband_copy(static_cast<TGraphErrors*>(
                sidebands->Clone("g_sideband_points")));
            std::unique_ptr<TGraphErrors> pull_copy(static_cast<TGraphErrors*>(
                pull_graph->Clone("g_sideband_pulls")));
            std::unique_ptr<TF1> background_copy(static_cast<TF1*>(
                f_bg->Clone("fitted_bg")));
            TF1 signal_copy(signal_fit);
            signal_copy.SetName("fitted_gaus_on_bgsub");
            TNamed model_description("background_model",
                "PEPSI-motivated Fermi/logistic: A/(1+exp((m-m_turn)/width))");

            input_copy->Write();
            final_copy->Write();
            sideband_copy->Write();
            pull_copy->Write();
            background_copy->Write();
            signal_copy.Write();
            covariance.Write("background_fit_covariance");
            model_description.Write();
            output_file.Close();
            std::cout << "[nps::PEPSICombBG] wrote diagnostics: "
                      << root_path << "\n";
        }

        // h_before, pull_graph, drawn stack primitives, and finally canvas are
        // destroyed automatically in reverse declaration order.
    }

    return result;
}

}  // namespace nps

#endif  // NPS_COMB_BG_PEPSI_H
