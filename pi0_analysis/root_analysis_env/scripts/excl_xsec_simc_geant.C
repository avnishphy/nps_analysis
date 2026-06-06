

// Forward declaration for find_bin (needed for find_bin_multi)
int find_bin(double x, const std::vector<double>& edges);
// excl_xsec.C
//
// How to run:
//   From a shell:
//     cd /work/hallc/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/scripts
//     root -l -b -q 'excl_xsec.C'
//
//   With explicit input/output arguments:
//     root -l -b -q 'excl_xsec.C("/w/hallc-scshelf2102/nps/singhav/geant4_simc/HallC_NPS/DVCS_evt_gen/DVCS/build/nps_excl_x60_4b_500k.root","/work/hallc/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/production_wfpi0/combined_branches_LH2_wfpi0.root","excl_xsec_output",1.0,1.0,true,0.5)'
//
//   Last argument is optional: min cluster energy cut in GeV (default 0.5).
//
//   Outputs:
//     - ROOT file: <out_dir>/excl_xsec_products.root
//     - PNG plots: <out_dir>/*.png
//
// Purpose:
//   Build the kinematic objects needed for pi0 cross-section extraction and
//   migration corrections using data + simulation files.
//
// Implemented outputs:
//   1) t' = t - t_min histograms for data and simulation.
//   2) 2D histograms of Q2 vs t' for data and simulation.
//   3) Five t' bins chosen from quantiles (approximately equal occupancy).
//   4) Twelve phi bins in each t' bin.
//   5) 4D phase-space histograms in (t', Q2, xB, phi).
//   6) Migration matrices per t' bin (phi_reco vs phi_gen), including
//      physics-weighted matrices for F0, FLT, FTT, FLTp terms.
//
// Physics model used for migration weighting:
//   K_r^n = (psf / 2pi) * (1 / Ngen) * Sum_i [F_N(i) * Gamma_i]
//   with
//     F0    = 1
//     FLT   = sqrt(2*eps*(1+eps)) * cos(phi)
//     FTT   = eps * cos(2phi)
//     FLTp  = h * sqrt(2*eps*(1-eps)) * sin(phi)
//
//   eps = (1 - y - Q2/(4E^2)) / (1 - y + y^2/2 + Q2/(4E^2)), y = nu/E.
//   Gamma = (alpha/(8pi)) * Q2/(M^2 K^2) * (1-xB)/xB^3 * 1/(1-eps)
//   K (Hand convention) = (W^2 - M^2)/(2M)
//
// Notes:
//   - Data tree: physics (combined_branches_LH2_wfpi0.root)
//   - Sim tree : nerd   (nps_excl_x60_4b_500k.root)
//   - Generated SIMC branches: Q2i, Wi, ti, phipqi
//   - Reconstructed quantities are rebuilt from GEANT-reconstructed observables
//     (electron momentum + photon hit geometry), not copied from t/phipq.
//   - Reco kinematics are computed via nps::compute_physics_vars_from_detector
//     from ../src/nps_physics_var.h for consistency with the main analysis code.
//
// Coordinate convention note:
//   This macro follows the same convention used in simc_pi0_analysis.C:
//   - azimuth-like branches (phi, phipq, phipqi) are in radians from atan2,
//     naturally in [-pi, pi], then wrapped to [0, 2pi) for binning.
//   - when available, xB is built from (Q2, nu) from the same physics frame,
//     consistent with nps::compute_physics_vars_from_detector usage.

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THnSparse.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TVector3.h"

#include "../src/nps_physics_var.h"

namespace {

// Configurable flag: use simulation event weights (s_weight) or set all weights to 1.0
constexpr bool USE_SIM_WEIGHT = true; // Set false to ignore s_weight and use unweighted simulation

// Minimum number of events per multidimensional bin
constexpr int MIN_EVENTS_PER_BIN = 50;

// Helper: multidimensional quantile binning with minimum count per bin
struct MultiDimBinning {
	std::vector<double> t_edges;
	std::vector<double> phi_edges;
	std::vector<double> q2_edges;
	std::vector<double> xb_edges;
	int n_t, n_phi, n_q2, n_xb;
};

// Returns bin edges for a variable using quantiles, with a minimum count per bin
std::vector<double> quantile_edges_min_count(const std::vector<double>& input, int n_bins, int min_count) {
	std::vector<double> values;
	values.reserve(input.size());
	for (double v : input) {
		if (std::isfinite(v)) values.push_back(v);
	}
	if (values.empty() || n_bins <= 0) return {};
	std::sort(values.begin(), values.end());
	// Adjust n_bins if not enough events for min_count
	int max_bins = std::max(1, static_cast<int>(values.size() / min_count));
	n_bins = std::min(n_bins, max_bins);
	std::vector<double> edges;
	edges.reserve(static_cast<size_t>(n_bins) + 1);
	edges.push_back(values.front());
	for (int i = 1; i < n_bins; ++i) {
		size_t idx = static_cast<size_t>(i * values.size() / n_bins);
		if (idx >= values.size()) idx = values.size() - 1;
		edges.push_back(values[idx]);
	}
	edges.push_back(values.back());
	// Enforce strict monotonicity
	for (size_t i = 1; i < edges.size(); ++i) {
		if (!(edges[i] > edges[i - 1])) edges[i] = edges[i - 1] + 1e-6;
	}
	return edges;
}

// Helper: find multidimensional bin index
int find_bin_multi(double t, double phi, double q2, double xb,
				   const std::vector<double>& t_edges,
				   const std::vector<double>& phi_edges,
				   const std::vector<double>& q2_edges,
				   const std::vector<double>& xb_edges,
				   int& tbin, int& phibin, int& q2bin, int& xbbin) {
	tbin = find_bin(t, t_edges);
	phibin = find_bin(phi, phi_edges);
	q2bin = find_bin(q2, q2_edges);
	xbbin = find_bin(xb, xb_edges);
	if (tbin < 0 || phibin < 0 || q2bin < 0 || xbbin < 0) return -1;
	// Flattened bin index
	int idx = ((tbin * (phi_edges.size() - 1) + phibin)
			  * (q2_edges.size() - 1) + q2bin)
			  * (xb_edges.size() - 1) + xbbin;
	return idx;
}

constexpr double kPi = TMath::Pi();
constexpr double kTwoPi = 2.0 * TMath::Pi();
constexpr double kAlphaEM = 1.0 / 137.035999084;
constexpr double kMProton = 0.9382720813;   // GeV
constexpr double kMPi0 = 0.1349768;         // GeV
constexpr double kBeamE = 10.538;           // GeV
constexpr int kNTPrimeBins = 5;
constexpr int kNPhiBins = 12;

struct DataEvent {
	double t;
	double tmin;
	double tprime;
	double q2;
	double xb;
	double phi;     // radians in [0,2pi)
	double weight;  // event weight from data (pi0_weight * scale * is_exclusive)
};

struct SimEvent {
	// Reconstructed quantities
	double t_rec;
	double tmin_rec;
	double tprime_rec;
	double q2_rec;
	double xb_rec;
	double phi_rec;

	// Generated quantities (truth)
	double t_gen;
	double tmin_gen;
	double tprime_gen;
	double q2_gen;
	double xb_gen;
	double phi_gen;

	double epsilon_gen;
	double gamma_gen;
	double sim_weight;
};

struct RecoKinematics {
	double q2;
	double w;
	double nu;
	double xb;
	double t;
	double tmin;
	double phi;
	bool valid;
};

double normalize_phi_rad(double phi_input) {
	// In this analysis chain phi variables are atan2 outputs (radians), as in
	// simc_pi0_analysis.C. We only wrap to [0, 2pi) and do not reinterpret
	// the values as degrees.
	double phi = phi_input;
	while (phi < 0.0) phi += kTwoPi;
	while (phi >= kTwoPi) phi -= kTwoPi;
	return phi;
}

double compute_xb_from_q2_nu(double q2, double nu) {
	const double denom = 2.0 * kMProton * nu;
	if (!(denom > 0.0)) return std::numeric_limits<double>::quiet_NaN();
	return q2 / denom;
}

double compute_xb_from_q2_w(double q2, double w) {
	const double denom = w * w + q2 - kMProton * kMProton;
	if (denom <= 0.0) return std::numeric_limits<double>::quiet_NaN();
	return q2 / denom;
}

double compute_tmin(double q2, double w) {
	if (!(q2 > 0.0) || !(w > (kMProton + kMPi0))) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	const double e_gamma = (w * w - kMProton * kMProton - q2) / (2.0 * w);
	const double p_gamma2 = e_gamma * e_gamma + q2;
	if (p_gamma2 <= 0.0) return std::numeric_limits<double>::quiet_NaN();
	const double p_gamma = std::sqrt(p_gamma2);

	const double e_pi = (w * w + kMPi0 * kMPi0 - kMProton * kMProton) / (2.0 * w);
	const double p_pi2 = e_pi * e_pi - kMPi0 * kMPi0;
	if (p_pi2 <= 0.0) return std::numeric_limits<double>::quiet_NaN();
	const double p_pi = std::sqrt(p_pi2);

	// Forward meson angle in gamma*-p CM gives t closest to 0 (t_min in |t| sense).
	const double tmin = kMPi0 * kMPi0 - q2 - 2.0 * e_gamma * e_pi + 2.0 * p_gamma * p_pi;
	return tmin;
}

double compute_epsilon(double q2, double xb) {
	if (!(q2 > 0.0) || !(xb > 0.0) || !(xb < 1.0)) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	const double nu = q2 / (2.0 * kMProton * xb);
	const double y = nu / kBeamE;
	const double q2_over_4e2 = q2 / (4.0 * kBeamE * kBeamE);
	const double num = 1.0 - y - q2_over_4e2;
	const double den = 1.0 - y + 0.5 * y * y + q2_over_4e2;
	if (den <= 0.0) return std::numeric_limits<double>::quiet_NaN();
	const double eps = num / den;
	return eps;
}

double compute_gamma_flux(double q2, double xb, double eps, double w) {
	if (!(q2 > 0.0) || !(xb > 0.0) || !(xb < 1.0)) {
		return std::numeric_limits<double>::quiet_NaN();
	}
	if (!(eps < 1.0)) return std::numeric_limits<double>::quiet_NaN();

	// const double k_hand = (w * w - kMProton * kMProton) / (2.0 * kMProton); // trying to use the beam energy directly
	// if (!(k_hand > 0.0)) return std::numeric_limits<double>::quiet_NaN();

	const double pref = kAlphaEM / (8.0 * kPi);
	const double gamma = pref * (q2 / (kMProton * kMProton * kBeamE * kBeamE)) *
											 ((1.0 - xb) / (xb * xb * xb)) * (1.0 / (1.0 - eps));
	return gamma;
}

RecoKinematics build_reco_kinematics_from_geant(double sc_e_px_mev,
																	 double sc_e_py_mev,
																	 double sc_e_pz_mev,
																	 double g1_e,
																	 double g1_x,
																	 double g1_y,
																	 double g2_e,
																	 double g2_x,
																	 double g2_y) {
	RecoKinematics out{};
	out.valid = false;

	if (!(g1_e > 0.0) || !(g2_e > 0.0)) return out;

	// Match simc_pi0_analysis.C coordinate transform (SIMC -> Hall C):
	// px = sc_e_Py/1000, py = -sc_e_Px/1000, pz = sc_e_Pz/1000.
	const double px_e = sc_e_py_mev / 1000.0;
	const double py_e = -sc_e_px_mev / 1000.0;
	const double pz_e = sc_e_pz_mev / 1000.0;
	const double me = 0.00051099895;
	const double ee = std::sqrt(px_e * px_e + py_e * py_e + pz_e * pz_e + me * me);

	const nps::PhysicsVars phys = nps::compute_physics_vars_from_detector(
			kBeamE,
			ee, px_e, py_e, pz_e,
			g1_e, g1_x, g1_y,
			g2_e, g2_x, g2_y,
			nps::kDefaultZ_NPS_cm,
			-17.51,
			false);

	out.q2 = phys.Q2;
	out.w = phys.W;
	out.nu = phys.nu;
	out.xb = phys.xB;
	out.t = phys.t;
	out.tmin = phys.tmin;
	out.phi = normalize_phi_rad(phys.phi);
	out.valid = std::isfinite(out.q2) && std::isfinite(out.w) && std::isfinite(out.nu) &&
						std::isfinite(out.xb) && std::isfinite(out.t) && std::isfinite(out.tmin) &&
						std::isfinite(out.phi);
	return out;
}

double fn_unpolarized(double, double) {
	return 1.0;
}

double fn_lt(double phi, double eps) {
	if (eps < 0.0) return 0.0;
	return std::sqrt(std::max(0.0, 2.0 * eps * (1.0 + eps))) * std::cos(phi);
}

double fn_tt(double phi, double eps) {
	return eps * std::cos(2.0 * phi);
}

double fn_ltp(double phi, double eps, double helicity) {
	if (eps < 0.0 || eps > 1.0) return 0.0;
	return helicity * std::sqrt(std::max(0.0, 2.0 * eps * (1.0 - eps))) * std::sin(phi);
}

// ...existing code...

std::vector<double> quantile_edges(const std::vector<double>& input, int n_bins) {
	std::vector<double> values;
	values.reserve(input.size());
	for (double v : input) {
		if (std::isfinite(v)) values.push_back(v);
	}

	if (values.empty() || n_bins <= 0) return {};
	std::sort(values.begin(), values.end());

	std::vector<double> edges;
	edges.reserve(static_cast<size_t>(n_bins) + 1);
	edges.push_back(values.front());

	for (int i = 1; i < n_bins; ++i) {
		const double q = static_cast<double>(i) / static_cast<double>(n_bins);
		const size_t idx = static_cast<size_t>(q * (values.size() - 1));
		edges.push_back(values[idx]);
	}
	edges.push_back(values.back());

	// Enforce strict monotonicity if duplicated quantile values appear.
	for (size_t i = 1; i < edges.size(); ++i) {
		if (!(edges[i] > edges[i - 1])) edges[i] = edges[i - 1] + 1e-6;
	}
	return edges;
}

void style_2d(TH2D* h, const char* ztitle) {
	h->GetXaxis()->CenterTitle(true);
	h->GetYaxis()->CenterTitle(true);
	h->GetZaxis()->SetTitle(ztitle);
	h->GetZaxis()->CenterTitle(true);
	h->SetStats(false);
}


}  // namespace

// Move find_bin outside the anonymous namespace so ROOT can resolve it
int find_bin(double x, const std::vector<double>& edges) {
	if (edges.size() < 2) return -1;
	if (x < edges.front() || x > edges.back()) return -1;
	for (size_t i = 0; i + 1 < edges.size(); ++i) {
		const bool in_bin = (x >= edges[i]) && ((x < edges[i + 1]) || (i + 2 == edges.size() && x <= edges[i + 1]));
		if (in_bin) return static_cast<int>(i);
	}
	return -1;
}


void excl_xsec_simc_geant(
	const char* sim_file_gen = "/w/hallc-scshelf2102/nps/singhav/geant4_simc/HallC_NPS/DVCS_evt_gen/DVCS/build/nps_excl_x60_4b_500k.root",
	const char* sim_file_rec = "output/plots/x60_4b/production_wfpi0/simc_pi0_analysis_output_smeared.root",
	const char* data_file = "/work/hallc/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env/output/plots/x60_4b/production_wfpi0/combined_branches_LH2_wfpi0.root",
	const char* out_dir = "./output/plots/x60_4b/production_wfpi0/excl_xsec_output",
	double phase_space_factor = 1.0,
	double beam_helicity = 1.0,
	bool apply_data_exclusive_cut = true,
	double min_cluster_energy_gev = 0.5)
{

	gStyle->SetOptStat(0);
	gStyle->SetPalette(kBird);
	gSystem->mkdir(out_dir, true);

	// 1D histograms for tprime (t')
	TH1D* h_tprime_data_1d = new TH1D("h_tprime_data_1d", "Data t';t' [GeV^{2}];Weighted counts", 120, -4.0, 1.5);
	TH1D* h_tprime_sim_1d = new TH1D("h_tprime_sim_1d", "Sim t';t' [GeV^{2}];Weighted counts", 120, -4.0, 1.5);

	TFile* f_data = TFile::Open(data_file, "READ");
	TFile* f_sim_gen = TFile::Open(sim_file_gen, "READ");
	TFile* f_sim_rec = TFile::Open(sim_file_rec, "READ");

	// --- Use per-event charge_uC branch for normalization ---
	if (!f_data || f_data->IsZombie()) {
		std::cerr << "ERROR: cannot open data file: " << data_file << std::endl;
		return;
	}
	if (!f_sim_gen || f_sim_gen->IsZombie()) {
		std::cerr << "ERROR: cannot open simulation file (gen): " << sim_file_gen << std::endl;
		f_data->Close();
		return;
	}
	if (!f_sim_rec || f_sim_rec->IsZombie()) {
		std::cerr << "ERROR: cannot open simulation file (rec): " << sim_file_rec << std::endl;
		f_data->Close();
		f_sim_gen->Close();
		return;
	}

	TTree* t_data = dynamic_cast<TTree*>(f_data->Get("physics"));
	TTree* t_sim_gen = dynamic_cast<TTree*>(f_sim_gen->Get("nerd"));
	TTree* t_sim_rec = dynamic_cast<TTree*>(f_sim_rec->Get("simulation"));
	if (!t_data) {
		std::cerr << "ERROR: data tree 'physics' not found." << std::endl;
		f_data->Close();
		f_sim_gen->Close();
		f_sim_rec->Close();
		return;
	}
	if (!t_sim_gen) {
		std::cerr << "ERROR: sim tree 'nerd' not found in gen file." << std::endl;
		f_data->Close();
		f_sim_gen->Close();
		f_sim_rec->Close();
		return;
	}
	if (!t_sim_rec) {
		std::cerr << "ERROR: sim tree 'simulation' not found in rec file." << std::endl;
		f_data->Close();
		f_sim_gen->Close();
		f_sim_rec->Close();
		return;
	}

	// ------------------------------
	// Data branch hookup
	// ------------------------------
	double d_t = 0.0, d_tmin = 0.0, d_q2 = 0.0, d_xb = 0.0, d_phi = 0.0, d_pi0_weight = 1.0;
	double d_cluster_e1 = 0.0, d_cluster_e2 = 0.0;
	float d_charge_uC = 1.0f;
	float d_scale = 1.0f;
	int d_is_exclusive = 1;

	t_data->SetBranchAddress("t", &d_t);
	t_data->SetBranchAddress("tmin", &d_tmin);
	t_data->SetBranchAddress("Q2", &d_q2);
	t_data->SetBranchAddress("xB", &d_xb);
	t_data->SetBranchAddress("phi", &d_phi);
	t_data->SetBranchAddress("scale", &d_scale);
	const bool has_data_charge_uC = (t_data->GetBranch("charge_uC") != nullptr);
	if (has_data_charge_uC) {
		t_data->SetBranchAddress("charge_uC", &d_charge_uC);
	} else {
		std::cerr << "ERROR: data branch 'charge_uC' is required for normalization." << std::endl;
		f_data->Close();
		f_sim_gen->Close();
		f_sim_rec->Close();
		return;
	}
	const bool has_data_pi0_weight = (t_data->GetBranch("pi0_weight") != nullptr);
	if (has_data_pi0_weight) {
		t_data->SetBranchAddress("pi0_weight", &d_pi0_weight);
	} else {
		std::cerr << "ERROR: data branch 'pi0_weight' is required for requested event weighting." << std::endl;
		f_data->Close();
		f_sim_gen->Close();
		f_sim_rec->Close();
		return;
	}
	const bool has_data_cluster_e1 = (t_data->GetBranch("cluster_e_1") != nullptr);
	const bool has_data_cluster_e2 = (t_data->GetBranch("cluster_e_2") != nullptr);
	if (has_data_cluster_e1) t_data->SetBranchAddress("cluster_e_1", &d_cluster_e1);
	if (has_data_cluster_e2) t_data->SetBranchAddress("cluster_e_2", &d_cluster_e2);
	const bool has_data_is_exclusive = (t_data->GetBranch("is_exclusive") != nullptr);
	if (has_data_is_exclusive) {
		t_data->SetBranchAddress("is_exclusive", &d_is_exclusive);
	} else {
		std::cerr << "ERROR: data branch 'is_exclusive' is required for exclusive-only weighting." << std::endl;
		f_data->Close();
		f_sim_gen->Close();
		f_sim_rec->Close();
		return;
	}

	// ------------------------------
	// Simulation branch hookup (reco: from smeared file, branch 'simulation')
	// ------------------------------
	Float_t srec_t = 0.0, srec_tmin = 0.0, srec_q2 = 0.0, srec_xb = 0.0, srec_phi = 0.0, srec_full_weight = 1.0f;
	Int_t srec_is_exclusive = 1;
		t_sim_rec->SetBranchAddress("t", &srec_t);
		t_sim_rec->SetBranchAddress("tmin", &srec_tmin);
		t_sim_rec->SetBranchAddress("Q2", &srec_q2);
		t_sim_rec->SetBranchAddress("xB", &srec_xb);
		t_sim_rec->SetBranchAddress("phi", &srec_phi);
		t_sim_rec->SetBranchAddress("is_exclusive", &srec_is_exclusive);
		t_sim_rec->SetBranchAddress("full_weight", &srec_full_weight);
		ULong64_t srec_event_id = 0;
		t_sim_rec->SetBranchAddress("event_id", &srec_event_id);

	// Simulation branch hookup (gen: from nerd tree in gen file)
	float s_q2i = 0.0f, s_wi = 0.0f, s_ti = 0.0f, s_phii = 0.0f;
		t_sim_gen->SetBranchAddress("Q2i", &s_q2i);
		t_sim_gen->SetBranchAddress("Wi", &s_wi);
		t_sim_gen->SetBranchAddress("ti", &s_ti);
		t_sim_gen->SetBranchAddress("phipqi", &s_phii);
		Int_t sgen_evtNb_int = 0;
		t_sim_gen->SetBranchAddress("evtNb", &sgen_evtNb_int);

	// ------------------------------
	// Read and cache selected events
	// ------------------------------
	std::vector<DataEvent> data_events;
	std::vector<SimEvent> sim_events;
	data_events.reserve(static_cast<size_t>(t_data->GetEntries()));
	sim_events.reserve(static_cast<size_t>(std::min(t_sim_rec->GetEntries(), t_sim_gen->GetEntries())));

	// Fill data_events with all events for histogramming, but mark t'<0 for analysis/quantiles
	for (Long64_t i = 0; i < t_data->GetEntries(); ++i) {
	       t_data->GetEntry(i);
	       if (d_is_exclusive != 1) continue;
	       if (has_data_cluster_e1 && has_data_cluster_e2) {
		       if (d_cluster_e1 < min_cluster_energy_gev || d_cluster_e2 < min_cluster_energy_gev) continue;
	       }

	       const double tprime = d_t - d_tmin;
	       const double phi = normalize_phi_rad(d_phi);
	       const double wt_scale = std::isfinite(d_scale) ? static_cast<double>(d_scale) : 1.0;
	       const double wt_pi0 = std::isfinite(d_pi0_weight) ? d_pi0_weight : 1.0;
	       double charge_mC = std::isfinite(d_charge_uC) && d_charge_uC > 0.0 ? static_cast<double>(d_charge_uC) / 1000.0 : 1.0;
	       const double wt = wt_pi0 * wt_scale * static_cast<double>(d_is_exclusive) / charge_mC;
	       if (!(wt > 0.0) || !std::isfinite(wt)) continue;

	       if (!std::isfinite(tprime) || !std::isfinite(d_q2) || !std::isfinite(d_xb) || !std::isfinite(phi)) continue;
	       if (d_q2 <= 0.0 || d_xb <= 0.0 || d_xb >= 1.0) continue;

	       data_events.push_back({d_t, d_tmin, tprime, d_q2, d_xb, phi, wt});
	}

	// Fill sim_events with all events for histogramming, but mark t'<0 for analysis/quantiles
	// Use reco variables from smeared file, gen variables from nerd tree (matched by entry index)
	Long64_t n_sim = std::min(t_sim_rec->GetEntries(), t_sim_gen->GetEntries());
		// Build a map from evtNb to gen event index for fast lookup
		std::unordered_map<ULong64_t, Long64_t> gen_evtid_to_idx;
		for (Long64_t i = 0; i < t_sim_gen->GetEntries(); ++i) {
			t_sim_gen->GetEntry(i);
			gen_evtid_to_idx[static_cast<ULong64_t>(sgen_evtNb_int)] = i;
		}

		Long64_t n_sim_rec = t_sim_rec->GetEntries();
		for (Long64_t i = 0; i < n_sim_rec; ++i) {
			t_sim_rec->GetEntry(i);
			if (srec_is_exclusive != 1) continue;
			// Find matching gen event by event_id
			auto it = gen_evtid_to_idx.find(srec_event_id);
			if (it == gen_evtid_to_idx.end()) continue; // No match
			Long64_t igen = it->second;
			t_sim_gen->GetEntry(igen);

			// Reco variables (from smeared file)
			const double tr = srec_t;
			const double tmin_r = srec_tmin;
			const double tprime_r = tr - tmin_r;
			const double q2r = srec_q2;
			const double xbr = srec_xb;
			const double phir = normalize_phi_rad(srec_phi);
			const double wt = USE_SIM_WEIGHT ? static_cast<double>(srec_full_weight) : 1.0;

			// Gen variables (from nerd tree)
			const double q2g = static_cast<double>(s_q2i);
			const double wg = static_cast<double>(s_wi);
			const double tg = -static_cast<double>(s_ti);
			const double phig = normalize_phi_rad(static_cast<double>(s_phii));
			const double tmin_g = compute_tmin(q2g, wg);
			const double tprime_g = tg - tmin_g;
			const double xbg = compute_xb_from_q2_w(q2g, wg);
			double eps_g = compute_epsilon(q2g, xbg);
			double gam_g = 1.0;

			if (!std::isfinite(tprime_r) || !std::isfinite(tprime_g)) continue;
			if (!std::isfinite(q2r) || !std::isfinite(q2g) || q2r <= 0.0 || q2g <= 0.0) continue;
			if (!std::isfinite(xbr) || !std::isfinite(xbg) || xbr <= 0.0 || xbr >= 1.0 || xbg <= 0.0 || xbg >= 1.0) continue;

			sim_events.push_back({
					tr, tmin_r, tprime_r, q2r, xbr, phir,
					tg, tmin_g, tprime_g, q2g, xbg, phig,
					eps_g, gam_g, wt
			});
		}

	if (data_events.empty() || sim_events.empty()) {
		std::cerr << "ERROR: no valid events after selections. data=" << data_events.size()
							<< " sim=" << sim_events.size() << std::endl;
		f_data->Close();
		f_sim_gen->Close();
		f_sim_rec->Close();
		return;
	}

	std::cout << "Selected events: data=" << data_events.size()
						<< " sim=" << sim_events.size() << std::endl;



	       // ------------------------------

				   // 5 tprime bins based on quantiles, 12 equally spaced phi bins
				   // ------------------------------
				   std::vector<double> tprime_vals;
				   double tprime_min = -1.0, tprime_max = 0.0;
				   for (const auto& e : data_events) {
					   if (e.tprime >= tprime_min && e.tprime <= tprime_max) {
						   tprime_vals.push_back(e.tprime);
					   }
				   }
				   int n_tprime = 5;
				   auto tprime_edges = quantile_edges(tprime_vals, n_tprime);
		       int n_phi = 12;
		       std::vector<double> phi_edges;
		       phi_edges.reserve(n_phi + 1);
		       for (int i = 0; i <= n_phi; ++i) {
			       phi_edges.push_back(i * kTwoPi / n_phi);
		       }


		       std::cout << "tprime quantile binning (" << n_tprime << " bins): ";
		       for (auto v : tprime_edges) std::cout << v << " ";
		       std::cout << "\nphi equally spaced bins (" << n_phi << " bins): ";
		       for (auto v : phi_edges) std::cout << v << " ";
		       std::cout << "\n";

	       // TODO: Cross-section calculation: (yield_data/yield_sim) * sig_bin_center
	       // (sig_bin_center to be defined by user)

	       // ...existing code...

	// ...existing code...
	// ------------------------------
	// Diagnostics: projections of 4D phase-space hists (after h4_* are declared)
	// ...existing code...

	// ------------------------------
	// Book histograms for binning visualization
	TH1D* h_t_data_1d = new TH1D("h_t_data_1d", "Data t;t [GeV^{2}];Weighted counts", 120, -4.0, 0.5);
	TH1D* h_q2_data_1d = new TH1D("h_q2_data_1d", "Data Q^{2};Q^{2} [GeV^{2}];Weighted counts", 120, 4.0, 8.0);
	TH1D* h_xb_data_1d = new TH1D("h_xb_data_1d", "Data x_{B};x_{B};Weighted counts", 120, 0.3, 0.9);
	TH1D* h_phi_data_1d = new TH1D("h_phi_data_1d", "Data #phi;#phi [rad];Weighted counts", 120, 0.0, kTwoPi);

	TH1D* h_t_sim_1d = new TH1D("h_t_sim_1d", "Sim t;t [GeV^{2}];Weighted counts", 120, -4.0, 0.5);
	TH1D* h_q2_sim_1d = new TH1D("h_q2_sim_1d", "Sim Q^{2};Q^{2} [GeV^{2}];Weighted counts", 120, 4.0, 8.0);
	TH1D* h_xb_sim_1d = new TH1D("h_xb_sim_1d", "Sim x_{B};x_{B};Weighted counts", 120, 0.3, 0.9);
	TH1D* h_phi_sim_1d = new TH1D("h_phi_sim_1d", "Sim #phi;#phi [rad];Weighted counts", 120, 0.0, kTwoPi);


	std::vector<TH1D*> h_phi_data;
	std::vector<TH1D*> h_phi_sim;
	h_phi_data.reserve(kNTPrimeBins);
	h_phi_sim.reserve(kNTPrimeBins);

		       for (int ib = 0; ib < kNTPrimeBins; ++ib) {
			       const std::string n1 = "h_phi_data_tbin" + std::to_string(ib);
			       const std::string n2 = "h_phi_sim_tbin" + std::to_string(ib);
			       const std::string title = "#phi in t' bin " + std::to_string(ib + 1) +
									       ";#phi [rad];Weighted counts";
			       h_phi_data.push_back(new TH1D(n1.c_str(), title.c_str(), kNPhiBins, 0.0, kTwoPi));
			       h_phi_sim.push_back(new TH1D(n2.c_str(), title.c_str(), kNPhiBins, 0.0, kTwoPi));
		       }

	// 4D phase-space histograms in (t', Q2, xB, phi)
	// No multidimensional histograms needed for 1D tprime/phi binning
	// (If needed, can add 2D/1D hists for tprime and phi only)


	// ------------------------------

	// Fill data hists for binning visualization
	for (const auto& e : data_events) {
		   h_t_data_1d->Fill(e.t, e.weight);
		   h_q2_data_1d->Fill(e.q2, e.weight);
		   h_xb_data_1d->Fill(e.xb, e.weight);
		   h_phi_data_1d->Fill(e.phi, e.weight);
		   h_tprime_data_1d->Fill(e.tprime, e.weight);
		   // Only fill analysis histograms for t'<0 events
		   if (e.tprime < 0.0) {
			   // Fill phi-in-tprime-bin histogram
			   int tprime_bin = find_bin(e.tprime, tprime_edges);
			   if (tprime_bin >= 0 && tprime_bin < (int)h_phi_data.size()) {
				   h_phi_data[tprime_bin]->Fill(e.phi, e.weight);
			   }
		   }
	}

	// ------------------------------

	// Fill simulation hists for binning visualization
	for (const auto& e : sim_events) {
		   h_t_sim_1d->Fill(e.t_rec, e.sim_weight);
		   h_q2_sim_1d->Fill(e.q2_rec, e.sim_weight);
		   h_xb_sim_1d->Fill(e.xb_rec, e.sim_weight);
		   h_phi_sim_1d->Fill(e.phi_rec, e.sim_weight);
		   h_tprime_sim_1d->Fill(e.tprime_rec, e.sim_weight);
		   // Only fill analysis histograms for t'<0 events
		   if (e.tprime_rec < 0.0 && e.tprime_gen < 0.0) {
			   // Fill phi-in-tprime-bin histogram for sim (use reco tprime)
			   int tprime_bin = find_bin(e.tprime_rec, tprime_edges);
			   if (tprime_bin >= 0 && tprime_bin < (int)h_phi_sim.size()) {
				   h_phi_sim[tprime_bin]->Fill(e.phi_rec, e.sim_weight);
			   }
		   }
	}


		       // ------------------------------
		       // Normalize simulation 1D hists to data integral if not using sim weights
		       // ------------------------------

	       // ------------------------------
		// Save all outputs
		// ------------------------------
		const std::string out_root = std::string(out_dir) + "/excl_xsec_products.root";
		TFile* fout = TFile::Open(out_root.c_str(), "RECREATE");

		       h_t_data_1d->Write();
		       h_q2_data_1d->Write();
		       h_xb_data_1d->Write();
		       h_phi_data_1d->Write();
		       h_t_sim_1d->Write();
		       h_q2_sim_1d->Write();
		       h_xb_sim_1d->Write();
		       h_phi_sim_1d->Write();
		       h_tprime_data_1d->Write();
		       h_tprime_sim_1d->Write();
			// No multidimensional histograms to write

		       // Plot and save 1D binning histograms
		       TCanvas c_t("c_t_binning", "t binning", 900, 700);
		       h_t_data_1d->SetLineColor(kBlue+1); h_t_data_1d->SetLineWidth(2);
			   TPad* pad1_t = new TPad("pad1_t", "pad1_t", 0, 0.4, 1, 1.0);
			   pad1_t->SetBottomMargin(0.04);
			pad1_t->Draw();
			pad1_t->cd();
			h_t_sim_1d->SetLineColor(kRed+1); h_t_sim_1d->SetLineWidth(2);
			double t_max = std::max(h_t_data_1d->GetMaximum(), h_t_sim_1d->GetMaximum());
			h_t_data_1d->SetMaximum(t_max * 1.15);
			   h_t_data_1d->Draw("hist"); h_t_sim_1d->Draw("hist same");
			   double t_data_total = h_t_data_1d->Integral(0, h_t_data_1d->GetNbinsX()+1);
			   double t_sim_total = h_t_sim_1d->Integral(0, h_t_sim_1d->GetNbinsX()+1);
			   double t_ratio_total = (t_sim_total > 0) ? t_data_total / t_sim_total : 0.0;
			   TLatex latex_t;
			   latex_t.SetNDC();
			   latex_t.SetTextSize(0.06);
			   latex_t.DrawLatex(0.15, 0.92, Form("Total ratio: %.3f", t_ratio_total));
			// Draw t bin edges
			       // No t_edges lines to draw
			c_t.cd();
			   TPad* pad2_t = new TPad("pad2_t", "pad2_t", 0, 0.0, 1, 0.4);
			   pad2_t->SetTopMargin(0.04);
			   pad2_t->SetBottomMargin(0.25);
			pad2_t->Draw();
			pad2_t->cd();
			TH1D* h_t_ratio = (TH1D*)h_t_data_1d->Clone("h_t_ratio");
			h_t_ratio->SetTitle("Data/Sim Ratio;t [GeV^{2}];Data/Sim");
			h_t_ratio->Divide(h_t_sim_1d);
			h_t_ratio->SetLineColor(kBlack); h_t_ratio->SetLineWidth(2);
			h_t_ratio->SetMinimum(0);
			h_t_ratio->SetMaximum(2);
			h_t_ratio->Draw("hist");
		       TLegend leg_t(0.65,0.75,0.88,0.89); leg_t.AddEntry(h_t_data_1d,"Data","l"); leg_t.AddEntry(h_t_sim_1d,"Sim","l"); leg_t.Draw();
		       c_t.SaveAs((std::string(out_dir)+"/binning_t.png").c_str());
		       c_t.Write();

			       TCanvas c_q2("c_q2_binning", "Q2 binning", 900, 900);
				   TPad* pad1_q2 = new TPad("pad1_q2", "pad1_q2", 0, 0.4, 1, 1.0);
				   pad1_q2->SetBottomMargin(0.04);
			       pad1_q2->Draw();
			       pad1_q2->cd();
			       h_q2_data_1d->SetLineColor(kBlue+1); h_q2_data_1d->SetLineWidth(2);
			       h_q2_sim_1d->SetLineColor(kRed+1); h_q2_sim_1d->SetLineWidth(2);
			       double q2_max = std::max(h_q2_data_1d->GetMaximum(), h_q2_sim_1d->GetMaximum());
			       h_q2_data_1d->SetMaximum(q2_max * 1.15);
				   h_q2_data_1d->Draw("hist"); h_q2_sim_1d->Draw("hist same");
				   double q2_data_total = h_q2_data_1d->Integral(0, h_q2_data_1d->GetNbinsX()+1);
				   double q2_sim_total = h_q2_sim_1d->Integral(0, h_q2_sim_1d->GetNbinsX()+1);
				   double q2_ratio_total = (q2_sim_total > 0) ? q2_data_total / q2_sim_total : 0.0;
				   TLatex latex_q2;
				   latex_q2.SetNDC();
				   latex_q2.SetTextSize(0.06);
				   latex_q2.DrawLatex(0.15, 0.92, Form("Total ratio: %.3f", q2_ratio_total));
				       // No q2_edges lines to draw
			       c_q2.cd();
				   TPad* pad2_q2 = new TPad("pad2_q2", "pad2_q2", 0, 0.0, 1, 0.4);
				   pad2_q2->SetTopMargin(0.04);
				   pad2_q2->SetBottomMargin(0.25);
			       pad2_q2->Draw();
			       pad2_q2->cd();
			       TH1D* h_q2_ratio = (TH1D*)h_q2_data_1d->Clone("h_q2_ratio");
			       h_q2_ratio->SetTitle("Data/Sim Ratio;Q^{2} [GeV^{2}];Data/Sim");
			       h_q2_ratio->Divide(h_q2_sim_1d);
			       h_q2_ratio->SetLineColor(kBlack); h_q2_ratio->SetLineWidth(2);
			       h_q2_ratio->SetMinimum(0);
			       h_q2_ratio->SetMaximum(2);
			       h_q2_ratio->Draw("hist");
			       TLegend leg_q2(0.65,0.75,0.88,0.89); leg_q2.AddEntry(h_q2_data_1d,"Data","l"); leg_q2.AddEntry(h_q2_sim_1d,"Sim","l"); leg_q2.Draw();
			       c_q2.SaveAs((std::string(out_dir)+"/binning_q2.png").c_str());
			       c_q2.Write();

			       TCanvas c_xb("c_xb_binning", "xb binning", 900, 900);
				   TPad* pad1_xb = new TPad("pad1_xb", "pad1_xb", 0, 0.4, 1, 1.0);
				   pad1_xb->SetBottomMargin(0.04);
			       pad1_xb->Draw();
			       pad1_xb->cd();
			       h_xb_data_1d->SetLineColor(kBlue+1); h_xb_data_1d->SetLineWidth(2);
			       h_xb_sim_1d->SetLineColor(kRed+1); h_xb_sim_1d->SetLineWidth(2);
			       double xb_max = std::max(h_xb_data_1d->GetMaximum(), h_xb_sim_1d->GetMaximum());
			       h_xb_data_1d->SetMaximum(xb_max * 1.15);
				   h_xb_data_1d->Draw("hist"); h_xb_sim_1d->Draw("hist same");
				   double xb_data_total = h_xb_data_1d->Integral(0, h_xb_data_1d->GetNbinsX()+1);
				   double xb_sim_total = h_xb_sim_1d->Integral(0, h_xb_sim_1d->GetNbinsX()+1);
				   double xb_ratio_total = (xb_sim_total > 0) ? xb_data_total / xb_sim_total : 0.0;
				   TLatex latex_xb;
				   latex_xb.SetNDC();
				   latex_xb.SetTextSize(0.06);
				   latex_xb.DrawLatex(0.15, 0.92, Form("Total ratio: %.3f", xb_ratio_total));
				       // No xb_edges lines to draw
			       c_xb.cd();
				   TPad* pad2_xb = new TPad("pad2_xb", "pad2_xb", 0, 0.0, 1, 0.4);
				   pad2_xb->SetTopMargin(0.04);
				   pad2_xb->SetBottomMargin(0.25);
			       pad2_xb->Draw();
			       pad2_xb->cd();
			       TH1D* h_xb_ratio = (TH1D*)h_xb_data_1d->Clone("h_xb_ratio");
			       h_xb_ratio->SetTitle("Data/Sim Ratio;x_{B};Data/Sim");
			       h_xb_ratio->Divide(h_xb_sim_1d);
			       h_xb_ratio->SetLineColor(kBlack); h_xb_ratio->SetLineWidth(2);
			       h_xb_ratio->SetMinimum(0);
			       h_xb_ratio->SetMaximum(2);
			       h_xb_ratio->Draw("hist");
			       TLegend leg_xb(0.65,0.75,0.88,0.89); leg_xb.AddEntry(h_xb_data_1d,"Data","l"); leg_xb.AddEntry(h_xb_sim_1d,"Sim","l"); leg_xb.Draw();
			       c_xb.SaveAs((std::string(out_dir)+"/binning_xb.png").c_str());
			       c_xb.Write();

			       TCanvas c_phi("c_phi_binning", "phi binning", 900, 900);
				   TPad* pad1_phi = new TPad("pad1_phi", "pad1_phi", 0, 0.4, 1, 1.0);
				   pad1_phi->SetBottomMargin(0.04);
			       pad1_phi->Draw();
			       pad1_phi->cd();
			       h_phi_data_1d->SetLineColor(kBlue+1); h_phi_data_1d->SetLineWidth(2);
			       h_phi_sim_1d->SetLineColor(kRed+1); h_phi_sim_1d->SetLineWidth(2);
			       double phi_max = std::max(h_phi_data_1d->GetMaximum(), h_phi_sim_1d->GetMaximum());
			       h_phi_data_1d->SetMaximum(phi_max * 1.15);
				   h_phi_data_1d->Draw("hist"); h_phi_sim_1d->Draw("hist same");
				   double phi_data_total = h_phi_data_1d->Integral(0, h_phi_data_1d->GetNbinsX()+1);
				   double phi_sim_total = h_phi_sim_1d->Integral(0, h_phi_sim_1d->GetNbinsX()+1);
				   double phi_ratio_total = (phi_sim_total > 0) ? phi_data_total / phi_sim_total : 0.0;
				   TLatex latex_phi;
				   latex_phi.SetNDC();
				   latex_phi.SetTextSize(0.06);
				   latex_phi.DrawLatex(0.15, 0.92, Form("Total ratio: %.3f", phi_ratio_total));
			       for (size_t i = 0; i < phi_edges.size(); ++i) {
				       TLine* l = new TLine(phi_edges[i], 0, phi_edges[i], phi_max * 1.15);
				       l->SetLineColor(kGreen+2); l->SetLineStyle(2); l->Draw();
			       }
			       c_phi.cd();
				   TPad* pad2_phi = new TPad("pad2_phi", "pad2_phi", 0, 0.0, 1, 0.4);
				   pad2_phi->SetTopMargin(0.04);
				   pad2_phi->SetBottomMargin(0.25);
			       pad2_phi->Draw();
			       pad2_phi->cd();
			       TH1D* h_phi_ratio = (TH1D*)h_phi_data_1d->Clone("h_phi_ratio");
			       h_phi_ratio->SetTitle("Data/Sim Ratio;#phi [rad];Data/Sim");
			       h_phi_ratio->Divide(h_phi_sim_1d);
			       h_phi_ratio->SetLineColor(kBlack); h_phi_ratio->SetLineWidth(2);
			       h_phi_ratio->SetMinimum(0);
			       h_phi_ratio->SetMaximum(2);
			       h_phi_ratio->Draw("hist");
			       TLegend leg_phi(0.65,0.75,0.88,0.89); leg_phi.AddEntry(h_phi_data_1d,"Data","l"); leg_phi.AddEntry(h_phi_sim_1d,"Sim","l"); leg_phi.Draw();
			       c_phi.SaveAs((std::string(out_dir)+"/binning_phi.png").c_str());
			       c_phi.Write();

			       // Plot and save tprime 1D binning histogram
			       TCanvas c_tprime("c_tprime_binning", "t' binning", 900, 900);
				   TPad* pad1_tprime = new TPad("pad1_tprime", "pad1_tprime", 0, 0.4, 1, 1.0);
				   pad1_tprime->SetBottomMargin(0.04);
			       pad1_tprime->Draw();
			       pad1_tprime->cd();
			       h_tprime_data_1d->SetLineColor(kBlue+1); h_tprime_data_1d->SetLineWidth(2);
			       h_tprime_sim_1d->SetLineColor(kRed+1); h_tprime_sim_1d->SetLineWidth(2);
							   double tprime_plot_max = std::max(h_tprime_data_1d->GetMaximum(), h_tprime_sim_1d->GetMaximum());
							   h_tprime_data_1d->SetMaximum(tprime_plot_max * 1.15);
							   h_tprime_data_1d->Draw("hist"); h_tprime_sim_1d->Draw("hist same");
				   double tprime_data_total = h_tprime_data_1d->Integral(0, h_tprime_data_1d->GetNbinsX()+1);
				   double tprime_sim_total = h_tprime_sim_1d->Integral(0, h_tprime_sim_1d->GetNbinsX()+1);
				   double tprime_ratio_total = (tprime_sim_total > 0) ? tprime_data_total / tprime_sim_total : 0.0;
				   TLatex latex_tprime;
				   latex_tprime.SetNDC();
				   latex_tprime.SetTextSize(0.06);
				   latex_tprime.DrawLatex(0.15, 0.92, Form("Total ratio: %.3f", tprime_ratio_total));
							   // Draw tprime bin edges
							   for (size_t i = 0; i < tprime_edges.size(); ++i) {
								   TLine* l = new TLine(tprime_edges[i], 0, tprime_edges[i], tprime_plot_max * 1.15);
								   l->SetLineColor(kGreen+2); l->SetLineStyle(2); l->Draw();
							   }
			       c_tprime.cd();
				   TPad* pad2_tprime = new TPad("pad2_tprime", "pad2_tprime", 0, 0.0, 1, 0.4);
				   pad2_tprime->SetTopMargin(0.04);
				   pad2_tprime->SetBottomMargin(0.25);
			       pad2_tprime->Draw();
			       pad2_tprime->cd();
			       TH1D* h_tprime_ratio = (TH1D*)h_tprime_data_1d->Clone("h_tprime_ratio");
			       h_tprime_ratio->SetTitle("Data/Sim Ratio;t' [GeV^{2}];Data/Sim");
			       h_tprime_ratio->Divide(h_tprime_sim_1d);
			       h_tprime_ratio->SetLineColor(kBlack); h_tprime_ratio->SetLineWidth(2);
			       h_tprime_ratio->SetMinimum(0);
			       h_tprime_ratio->SetMaximum(2);
			       h_tprime_ratio->Draw("hist");
			       TLegend leg_tprime(0.65,0.75,0.88,0.89); leg_tprime.AddEntry(h_tprime_data_1d,"Data","l"); leg_tprime.AddEntry(h_tprime_sim_1d,"Sim","l"); leg_tprime.Draw();
			       c_tprime.SaveAs((std::string(out_dir)+"/binning_tprime.png").c_str());
			       c_tprime.Write();


				// Plot phi distributions in each tprime bin (overlay data/sim, show total ratio)
				for (int ib = 0; ib < kNTPrimeBins; ++ib) {
					TCanvas c_phi_bin(Form("c_phi_tbin%d", ib+1), Form("#phi in t' bin %d", ib+1), 900, 900);
					TPad* pad1_phi_bin = new TPad(Form("pad1_phi_bin_%d", ib+1), "", 0, 0.4, 1, 1.0);
					pad1_phi_bin->SetBottomMargin(0.04);
					pad1_phi_bin->Draw();
					pad1_phi_bin->cd();
					h_phi_data[ib]->SetLineColor(kBlue+1); h_phi_data[ib]->SetLineWidth(2);
					h_phi_sim[ib]->SetLineColor(kRed+1); h_phi_sim[ib]->SetLineWidth(2);
					double phi_bin_max = std::max(h_phi_data[ib]->GetMaximum(), h_phi_sim[ib]->GetMaximum());
					h_phi_data[ib]->SetMaximum(phi_bin_max * 1.15);
					h_phi_data[ib]->Draw("hist"); h_phi_sim[ib]->Draw("hist same");
					// Draw phi bin edges
					for (size_t i = 0; i < phi_edges.size(); ++i) {
						TLine* l = new TLine(phi_edges[i], 0, phi_edges[i], phi_bin_max * 1.15);
						l->SetLineColor(kGreen+2); l->SetLineStyle(2); l->Draw();
					}
					double phi_data_total = h_phi_data[ib]->Integral(0, h_phi_data[ib]->GetNbinsX()+1);
					double phi_sim_total = h_phi_sim[ib]->Integral(0, h_phi_sim[ib]->GetNbinsX()+1);
					double phi_ratio_total = (phi_sim_total > 0) ? phi_data_total / phi_sim_total : 0.0;
					TLatex latex_phi_bin;
					latex_phi_bin.SetNDC();
					latex_phi_bin.SetTextSize(0.06);
					latex_phi_bin.DrawLatex(0.15, 0.92, Form("Total ratio: %.3f", phi_ratio_total));
					c_phi_bin.cd();
					TPad* pad2_phi_bin = new TPad(Form("pad2_phi_bin_%d", ib+1), "", 0, 0.0, 1, 0.4);
					pad2_phi_bin->SetTopMargin(0.04);
					pad2_phi_bin->SetBottomMargin(0.25);
					pad2_phi_bin->Draw();
					pad2_phi_bin->cd();
					TH1D* h_phi_ratio_bin = (TH1D*)h_phi_data[ib]->Clone(Form("h_phi_ratio_bin_%d", ib+1));
					h_phi_ratio_bin->SetTitle("Data/Sim Ratio;#phi [rad];Data/Sim");
					h_phi_ratio_bin->Divide(h_phi_sim[ib]);
					h_phi_ratio_bin->SetLineColor(kBlack); h_phi_ratio_bin->SetLineWidth(2);
					h_phi_ratio_bin->SetMinimum(0);
					h_phi_ratio_bin->SetMaximum(2);
					h_phi_ratio_bin->Draw("hist");
					c_phi_bin.SaveAs((std::string(out_dir)+Form("/binning_phi_tbin%d.png", ib+1)).c_str());
					c_phi_bin.Write();
				}
			}

