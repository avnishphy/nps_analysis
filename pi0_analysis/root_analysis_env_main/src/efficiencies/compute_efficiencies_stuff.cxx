/*
============================================================
NPS Hall C Efficiency & Livetime Calculation Driver
============================================================

This program computes the following quantities for each run/kinematic:

1. HMS PID Efficiency
	- Numerator: Events passing PID cuts (Cherenkov, Calorimeter, Tracking, Hodoscope)
	- Denominator: Events passing denominator cuts (looser PID, tracking, hodo)
	- Cuts: Event selection (evt_type==1, evnum in accepted_gevnum_ranges), PID cuts

2. HMS Cal/Cer Tagging Efficiencies
	- Numerator/Denominator: Events passing respective tagging cuts
	- Cuts: Event selection, tagging cuts

3. HMS Tracking Efficiency
	- Numerator: Events passing tracking cuts
	- Denominator: Events passing should-track cuts
	- Cuts: Event selection, tracking cuts

4. HMS Hodoscope 3/4 Efficiency
	- Numerator: Events with good hodo hits
	- Denominator: Events passing hodo denominator cuts
	- Cuts: Event selection, hodo cuts

5. NewGen EDTM Livetime
	- Numerator: Current-window TTree events with EDTM TDC > 1 ns and within 500 ns of the peak
	- Denominator: Scaler EDTM count from TSH tree, summed over intervals passing the current cut
	- Diagnostic: Good-event/evcount-gated EDTM counts are written separately, but do not define the livetime

6. Beam Time
	- Computed from TSH tree scaler intervals, summed where current and accepted evcount cuts pass

Event-level good selection uses g.evnum; scaler-level good selection uses the
corresponding TSHelH evcount ranges. Segment results are accumulated for each run/kinematic.

*/
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "config_csv_helper.h"
#include "eff_math_utils.h"
#include "efficiency_types.h"
#include "good_event_selection_helper.h"
#include "hms_cal_eff_tag_cer.h"
#include "hms_cer_eff_tag_cal.h"
#include "hms_hodo_3of4_eff.h"
#include "hms_pid_eff.h"
#include "hms_tracking_eff.h"
#include "newgen_edtm_livetime.h"
#include "prescale_beamtime_helper.h"
#include "root_file_discovery.h"

/*
Build:
	g++ -O3 -std=c++17 -Wall -Wextra -pedantic `root-config --cflags` compute_efficiencies_stuff.cxx -o compute_efficiencies_stuff `root-config --libs`
*/

namespace fs = std::filesystem;

std::ostream* g_log_stream = nullptr;

void log_info(const std::string& msg) {
	if (g_log_stream) {
		(*g_log_stream) << msg << std::endl;
	}
	std::cout << msg << std::endl;
}

void log_debug(const std::string& msg) {
	if (g_log_stream) {
		(*g_log_stream) << msg << std::endl;
	}
}

void log_warn(const std::string& msg) {
	if (g_log_stream) {
		(*g_log_stream) << msg << std::endl;
	}
	std::cerr << msg << std::endl;
}

namespace {

struct ProgramConfig {
	std::string config_csv_path = "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/config/nps_dvcs_all_kins_main.csv";
	std::string updated_root_dir = "/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/updated";
	std::string production_root_dir = "/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/production";
	std::string output_dir = "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/output/efficiency_stuff";

	std::vector<std::string> requested_kinematics;
	bool all_kinematics = false;
	bool non_interactive = false;

	std::vector<std::string> allowed_types = effstuff::default_allowed_types();
	std::optional<int> single_run_filter;
    std::optional<std::string> log_path;
};

void print_usage(const char* prog) {
	std::cerr
			<< "Usage: " << prog << " [options]\n"
			<< "Options:\n"
			<< "  --config <path>             Run configuration CSV path\n"
			<< "  --updated-dir <path>        Updated replay ROOT directory\n"
			<< "  --production-dir <path>     Production replay ROOT directory\n"
			<< "  --output-dir <path>         Output directory for CSV files\n"
			<< "  --kin <Kin_old>             Process one kinematic setting (repeatable)\n"
			<< "  --all-kins                  Process all kinematic settings\n"
			<< "  --types <a,b,c>             Allowed Type values (default: production,Production)\n"
			<< "  --run <run_number>          Restrict processing to a single run\n"
			<< "  --no-interactive            Disable terminal prompts\n"
			<< "  --log <file>                Write debug log to <file> (optional)\n"
			<< "  --help                      Show this message\n";
}

bool parse_args(int argc, char** argv, ProgramConfig& cfg) {
	for (int i = 1; i < argc; ++i) {
		const std::string arg = argv[i];

		auto require_value = [&](const std::string& key) -> std::string {
			if (i + 1 >= argc) {
				throw std::runtime_error("Missing value after " + key);
			}
			return std::string(argv[++i]);
		};

		try {
			if (arg == "--config") {
				cfg.config_csv_path = require_value(arg);
			} else if (arg == "--updated-dir") {
				cfg.updated_root_dir = require_value(arg);
			} else if (arg == "--production-dir") {
				cfg.production_root_dir = require_value(arg);
			} else if (arg == "--output-dir") {
				cfg.output_dir = require_value(arg);
			} else if (arg == "--kin") {
				cfg.requested_kinematics.push_back(require_value(arg));
			} else if (arg == "--all-kins") {
				cfg.all_kinematics = true;
			} else if (arg == "--types") {
				cfg.allowed_types = effstuff::split_comma_list(require_value(arg));
				if (cfg.allowed_types.empty()) {
					throw std::runtime_error("--types must not be empty");
				}
			} else if (arg == "--run") {
				cfg.single_run_filter = std::stoi(require_value(arg));
			} else if (arg == "--log") {
				cfg.log_path = require_value(arg);
			} else if (arg == "--no-interactive") {
				cfg.non_interactive = true;
			} else if (arg == "--help" || arg == "-h") {
				print_usage(argv[0]);
				return false;
			} else {
				throw std::runtime_error("Unknown argument: " + arg);
			}
		} catch (const std::exception& ex) {
			std::cerr << "Argument parsing error: " << ex.what() << "\n";
			print_usage(argv[0]);
			return false;
		}
	}

	if (cfg.all_kinematics && !cfg.requested_kinematics.empty()) {
		std::cerr << "Invalid options: use either --all-kins or --kin, not both.\n";
		return false;
	}

	return true;
}

std::vector<std::string> parse_kinematic_selection_input(const std::string& input, const std::vector<std::string>& available) {

	std::vector<std::string> selected;
	const std::string clean = effstuff::trim_copy(input);
	if (clean.empty()) {
		return selected;
	}

	if (effstuff::iequals(clean, "all")) {
		return available;
	}

	std::vector<std::string> tokens;
	std::stringstream ss(clean);
	std::string token;
	while (std::getline(ss, token, ',')) {
		const std::string t = effstuff::trim_copy(token);
		if (!t.empty()) {
			tokens.push_back(t);
		}
	}

	std::set<std::string> dedup;
	for (const auto& t : tokens) {
		bool matched = false;

		try {
			const int idx = std::stoi(t);
			if (idx >= 1 && idx <= static_cast<int>(available.size())) {
				dedup.insert(available[static_cast<size_t>(idx - 1)]);
				matched = true;
			}
		} catch (...) {
			// Not numeric, continue with name matching.
		}

		if (!matched) {
			for (const auto& kin : available) {
				if (kin == t) {
					dedup.insert(kin);
					matched = true;
					break;
				}
			}
		}
	}

	selected.assign(dedup.begin(), dedup.end());
	return selected;
}

std::vector<std::string> choose_kinematics(const ProgramConfig& cfg, const std::vector<std::string>& available_kinematics) {

	if (available_kinematics.empty()) {
		return {};
	}

	if (cfg.all_kinematics) {
		return available_kinematics;
	}

	if (!cfg.requested_kinematics.empty()) {
		std::set<std::string> dedup;
		for (const auto& requested : cfg.requested_kinematics) {
			if (std::find(available_kinematics.begin(), available_kinematics.end(), requested) == available_kinematics.end()) {
				throw std::runtime_error("Invalid kinematic selection: " + requested);
			}
			dedup.insert(requested);
		}
		return std::vector<std::string>(dedup.begin(), dedup.end());
	}

	if (cfg.non_interactive) {
		throw std::runtime_error("No kinematic selection provided. Use --kin or --all-kins when --no-interactive is enabled.");
	}

	std::cout << "Available kinematic settings (Kin_old):\n";
	for (size_t i = 0; i < available_kinematics.size(); ++i) {
		std::cout << "  " << (i + 1) << ") " << available_kinematics[i] << "\n";
	}
	std::cout << "Select one or more kinematics by number/name (comma-separated), or type 'all': ";

	std::string input;
	std::getline(std::cin, input);
	const auto selected = parse_kinematic_selection_input(input, available_kinematics);
	if (selected.empty()) {
		throw std::runtime_error("Invalid interactive kinematic selection.");
	}

	return selected;
}

bool event_is_selected(unsigned int evt_type, double gevnum, const effstuff::GoodSelectionSummary& selection) {

	const bool is_physics_evt = (evt_type == 1);
	const bool in_range = effstuff::value_in_ranges(
			static_cast<long long>(std::llround(gevnum)), selection.accepted_gevnum_ranges);
	return is_physics_evt && in_range;
}

std::string tree_branch_leaf_type(TTree* tree, const char* branch_name) {
	if (tree == nullptr || branch_name == nullptr) {
		return "";
	}

	TBranch* branch = tree->GetBranch(branch_name);
	if (branch == nullptr) {
		return "";
	}

	TLeaf* leaf = branch->GetLeaf(branch_name);
	if (leaf == nullptr && branch->GetListOfLeaves() != nullptr && branch->GetListOfLeaves()->GetEntries() > 0) {
		leaf = static_cast<TLeaf*>(branch->GetListOfLeaves()->At(0));
	}

	if (leaf == nullptr) {
		return "";
	}

	return leaf->GetTypeName();
}

bool validate_required_branch_type(
		TTree* tree,
		const char* branch_name,
		const char* expected_type,
		const std::string& filepath,
		int& missing_branch_segments) {
	const std::string actual_type = tree_branch_leaf_type(tree, branch_name);
	if (actual_type.empty()) {
		++missing_branch_segments;
		log_warn("[warn segment] Missing branch or leaf type for '" + std::string(branch_name) +
		         "' in file: " + filepath);
		return false;
	}

	if (actual_type != expected_type) {
		++missing_branch_segments;
		log_warn("[warn segment] Branch type mismatch for '" + std::string(branch_name) +
		         "' in file: " + filepath +
		         " (expected " + std::string(expected_type) +
		         ", got " + actual_type + ")");
		return false;
	}

	return true;
}

double compute_edtm_peak_from_hist(const std::map<int, long long>& edtm_hist) {
	if (edtm_hist.empty()) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	const double bin_width = 10.0;
	int max_bin = 0;
	long long max_count = -1;
	for (const auto& kv : edtm_hist) {
		if (kv.second > max_count) {
			max_count = kv.second;
			max_bin = kv.first;
		}
	}

	if (max_count <= 0) {
		return std::numeric_limits<double>::quiet_NaN();
	}

	return (static_cast<double>(max_bin) + 0.5) * bin_width;
}

// Combined pass: collects current-window EDTM pulse times and all efficiency
// accumulators in one read of the T tree. The EDTM peak is determined after all
// segments are scanned; the cached pulse times are then counted around that peak.
void accumulate_segment_histogram_and_metrics(
		const std::string& filepath,
		const effstuff::GoodSelectionSummary& selection,
		const effstuff::HMSPidCuts& pid_cuts,
		std::map<int, long long>& out_hist,
		std::vector<double>& out_current_window_edtm_tdc,
		std::vector<double>& out_selected_current_window_edtm_tdc,
		effstuff::HMSPidEffAccumulator& pid_acc,
		effstuff::HMSCalEffTagCerAccumulator& cal_tag_acc,
		effstuff::HMSCerEffTagCalAccumulator& cer_tag_acc,
		effstuff::HMSTrackingEffAccumulator& tracking_acc,
		effstuff::HMSHodo3of4EffAccumulator& hodo_acc,
		int& missing_branch_segments) {
	TFile file(filepath.c_str(), "READ");
	if (file.IsZombie()) {
		++missing_branch_segments;
		return;
	}

	auto* t = dynamic_cast<TTree*>(file.Get("T"));
	if (t == nullptr) {
		++missing_branch_segments;
		return;
	}

	const bool has_core =
			effstuff::tree_has_branch(t, "T.hms.hEDTM_tdcTimeRaw") &&
			effstuff::tree_has_branch(t, "g.evnum") &&
			effstuff::tree_has_branch(t, "fEvtHdr.fEvtType");
	if (!has_core) {
		++missing_branch_segments;
		return;
	}

	if (!validate_required_branch_type(t, "T.hms.hEDTM_tdcTimeRaw", "Double_t", filepath, missing_branch_segments) ||
			!validate_required_branch_type(t, "g.evnum", "Double_t", filepath, missing_branch_segments) ||
			!validate_required_branch_type(t, "fEvtHdr.fEvtType", "UInt_t", filepath, missing_branch_segments)) {
		return;
	}

	const bool have_pid_branches =
			effstuff::tree_has_branch(t, "H.cer.npeSum") &&
			effstuff::tree_has_branch(t, "H.cal.etotnorm") &&
			effstuff::tree_has_branch(t, "H.dc.ntrack") &&
			effstuff::tree_has_branch(t, "H.gtr.dp") &&
			effstuff::tree_has_branch(t, "H.gtr.th") &&
			effstuff::tree_has_branch(t, "H.gtr.ph");

	const bool have_tracking_branches =
			have_pid_branches &&
			effstuff::tree_has_branch(t, "H.hod.goodscinhit") &&
			effstuff::tree_has_branch(t, "H.hod.betanotrack");

	const bool have_hodo_branches =
			effstuff::tree_has_branch(t, "H.dc.ntrack") &&
			effstuff::tree_has_branch(t, "H.hod.1x.nhits") &&
			effstuff::tree_has_branch(t, "H.hod.1y.nhits") &&
			effstuff::tree_has_branch(t, "H.hod.2x.nhits") &&
			effstuff::tree_has_branch(t, "H.hod.2y.nhits");

	const bool have_event_current =
			effstuff::tree_has_branch(t, "H.BCM4A.scalerCurrent");
	if (!have_event_current) {
		++missing_branch_segments;
	}

	t->SetBranchStatus("*", 0);
	t->SetBranchStatus("T.hms.hEDTM_tdcTimeRaw", 1);
	t->SetBranchStatus("g.evnum", 1);
	t->SetBranchStatus("fEvtHdr.fEvtType", 1);
	if (have_event_current) {
		t->SetBranchStatus("H.BCM4A.scalerCurrent", 1);
	}

	if (have_pid_branches) {
		t->SetBranchStatus("H.cer.npeSum", 1);
		t->SetBranchStatus("H.cal.etotnorm", 1);
		t->SetBranchStatus("H.dc.ntrack", 1);
		t->SetBranchStatus("H.gtr.dp", 1);
		t->SetBranchStatus("H.gtr.th", 1);
		t->SetBranchStatus("H.gtr.ph", 1);
	} else if (have_hodo_branches) {
		t->SetBranchStatus("H.dc.ntrack", 1);
	}

	if (have_tracking_branches) {
		// Tracking efficiency intentionally uses the no-track hodoscope beta tag.
		// Keep H.hod.betanotrack separate from H.dc.ntrack to avoid biasing the
		// tracking denominator with the quantity being measured.
		t->SetBranchStatus("H.hod.goodscinhit", 1);
		t->SetBranchStatus("H.hod.betanotrack", 1);
	}

	if (have_hodo_branches) {
		t->SetBranchStatus("H.hod.1x.nhits", 1);
		t->SetBranchStatus("H.hod.1y.nhits", 1);
		t->SetBranchStatus("H.hod.2x.nhits", 1);
		t->SetBranchStatus("H.hod.2y.nhits", 1);
	}

	double edtm_tdc = 0.0;
	double gevnum = 0.0;
	unsigned int evt_type = 0;  // fEvtHdr.fEvtType is UInt_t in ROOT replay files
	double cer_npe = std::numeric_limits<double>::quiet_NaN();
	double cal_etotnorm = std::numeric_limits<double>::quiet_NaN();
	double ntrack = std::numeric_limits<double>::quiet_NaN();
	double hms_dp = std::numeric_limits<double>::quiet_NaN();
	double hms_th = std::numeric_limits<double>::quiet_NaN();
	double hms_ph = std::numeric_limits<double>::quiet_NaN();
	double good_scin = std::numeric_limits<double>::quiet_NaN();
	double hod_notrack = std::numeric_limits<double>::quiet_NaN();
	double hodo_1x = std::numeric_limits<double>::quiet_NaN();
	double hodo_1y = std::numeric_limits<double>::quiet_NaN();
	double hodo_2x = std::numeric_limits<double>::quiet_NaN();
	double hodo_2y = std::numeric_limits<double>::quiet_NaN();
	double event_current = std::numeric_limits<double>::quiet_NaN();

	t->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtm_tdc);
	t->SetBranchAddress("g.evnum", &gevnum);
	t->SetBranchAddress("fEvtHdr.fEvtType", &evt_type);
	if (have_event_current) {
		t->SetBranchAddress("H.BCM4A.scalerCurrent", &event_current);
	}

	if (have_pid_branches) {
		t->SetBranchAddress("H.cer.npeSum", &cer_npe);
		t->SetBranchAddress("H.cal.etotnorm", &cal_etotnorm);
		t->SetBranchAddress("H.dc.ntrack", &ntrack);
		t->SetBranchAddress("H.gtr.dp", &hms_dp);
		t->SetBranchAddress("H.gtr.th", &hms_th);
		t->SetBranchAddress("H.gtr.ph", &hms_ph);
	} else if (have_hodo_branches) {
		t->SetBranchAddress("H.dc.ntrack", &ntrack);
	}

	if (have_tracking_branches) {
		t->SetBranchAddress("H.hod.goodscinhit", &good_scin);
		t->SetBranchAddress("H.hod.betanotrack", &hod_notrack);
	}

	if (have_hodo_branches) {
		t->SetBranchAddress("H.hod.1x.nhits", &hodo_1x);
		t->SetBranchAddress("H.hod.1y.nhits", &hodo_1y);
		t->SetBranchAddress("H.hod.2x.nhits", &hodo_2x);
		t->SetBranchAddress("H.hod.2y.nhits", &hodo_2y);
	}

	// Read-ahead cache: reduces round-trips on Lustre/network filesystems.
	t->SetCacheSize(64 * 1024 * 1024);
	t->AddBranchToCache("*", kTRUE);
	t->StopCacheLearningPhase();

	const Long64_t n = t->GetEntries();
	for (Long64_t i = 0; i < n; ++i) {
		t->GetEntry(i);

		// Evaluate selection once and reuse for all efficiency accumulators.
		const bool selected = event_is_selected(evt_type, gevnum, selection);
		const bool current_ok = have_event_current
				? effstuff::current_in_selection_window(event_current, selection)
				: false;

		if (current_ok && std::isfinite(edtm_tdc) && edtm_tdc > 1.0) {
			++out_hist[static_cast<int>(edtm_tdc / 10.0)];
			out_current_window_edtm_tdc.push_back(edtm_tdc);
			if (selected) {
				out_selected_current_window_edtm_tdc.push_back(edtm_tdc);
			}
		}

		effstuff::HMSEventContext evt;
		evt.selected = selected;
		evt.edtm_tdc = edtm_tdc;
		evt.cer_npe_sum = cer_npe;
		evt.cal_etotnorm = cal_etotnorm;
		evt.hdc_ntrack = ntrack;
		evt.hms_dp = hms_dp;
		evt.hms_th = hms_th;
		evt.hms_ph = hms_ph;
		evt.hod_goodscinhit = good_scin;
		evt.hod_notrack = hod_notrack;
		evt.hod_1x_nhits = hodo_1x;
		evt.hod_1y_nhits = hodo_1y;
		evt.hod_2x_nhits = hodo_2x;
		evt.hod_2y_nhits = hodo_2y;

		if (have_pid_branches) {
			effstuff::accumulate_HMS_pid_eff(pid_acc, evt, pid_cuts);
			effstuff::accumulate_HMS_cal_eff_tag_cer(cal_tag_acc, evt, pid_cuts);
			effstuff::accumulate_HMS_cer_eff_tag_cal(cer_tag_acc, evt, pid_cuts);
		}

		if (have_tracking_branches) {
			effstuff::accumulate_HMS_tracking_eff(tracking_acc, evt, pid_cuts);
		}

		if (have_hodo_branches) {
			effstuff::accumulate_HMS_hodo_3of4_eff(hodo_acc, evt, pid_cuts);
		}
	}
}

long long count_edtm_peak_events(const std::vector<double>& current_window_edtm_tdc,
                                 double edtm_peak) {
	long long count = 0;
	for (const double edtm_tdc : current_window_edtm_tdc) {
		if (std::fabs(edtm_tdc - edtm_peak) <= 500.0) {
			++count;
		}
	}
	return count;
}

effstuff::RunProcessingRow process_run(
		const effstuff::RunConfigRow& run_cfg,
		const ProgramConfig& cfg,
		const SelectionSettings& selection_settings,
		std::vector<effstuff::SelectionReportRow>* selection_report_rows) {
	effstuff::RunProcessingRow out;
	out.run_number = run_cfg.run_number;
	out.kinematic_setting = run_cfg.kin_old;
	out.run_type = run_cfg.run_type;

	const effstuff::PrescaleInfo prescale = effstuff::build_prescale_info(run_cfg.prescale_token);
	out.prescale_token = prescale.prescale_token;
	out.ps_factor = prescale.ps_factor;
	out.which_TRIG = prescale.which_TRIG;

	const effstuff::LocatedRunFiles located = effstuff::locate_run_files_prefer_updated(
			run_cfg.run_number, cfg.updated_root_dir, cfg.production_root_dir);
	out.file_source_used = located.source;
	out.segment_count_found = located.segment_count;

	if (located.files.empty()) {
		out.run_processing_status = "missing_root_files";
		return out;
	}

	int n_segments_used = 0;
	std::map<int, long long> edtm_hist;
	std::vector<double> current_window_edtm_tdc;
	std::vector<double> selected_current_window_edtm_tdc;

	effstuff::BeamTimeAccumulator beam_acc;
	auto accumulate_charge = [](double& acc, double value) {
		if (!std::isfinite(value)) {
			return;
		}
		if (!std::isfinite(acc)) {
			acc = 0.0;
		}
		acc += value;
	};

	// Declare accumulators here so the combined-pass function can fill them
	// incrementally as each segment is processed (avoids a second file-open loop).
	effstuff::HMSPidCuts pid_cuts;
	effstuff::HMSPidEffAccumulator pid_acc;
	effstuff::HMSCalEffTagCerAccumulator cal_tag_acc;
	effstuff::HMSCerEffTagCalAccumulator cer_tag_acc;
	effstuff::HMSTrackingEffAccumulator tracking_acc;
	effstuff::HMSHodo3of4EffAccumulator hodo_acc;

	for (const auto& file_path : located.files) {
		const int seg_id = effstuff::parse_segment_from_filename(file_path, run_cfg.run_number);
		effstuff::GoodSelectionSummary selection = effstuff::build_good_selection_summary(file_path, selection_settings);

		if (selection_report_rows != nullptr) {
			effstuff::SelectionReportRow report;
			report.run_number = run_cfg.run_number;
			report.kinematic_setting = run_cfg.kin_old;
			report.run_type = run_cfg.run_type;
			report.segment_number = seg_id;
			report.segment_file = file_path;
			report.file_source_used = located.source;
			report.selection_ok = selection.ok;
			report.selection_message = selection.message;
			report.helicity_mode = selection.helicity_mode;
			report.quartet_snap_applied = selection.quartet_snap_applied;
			report.current_window_enabled = selection.current_window_enabled;
			report.mean_current_uA = selection.mean_current_uA;
			report.i0_used_uA = selection.i0_used_uA;
			report.current_min_uA = selection.current_min_uA;
			report.current_max_uA = selection.current_max_uA;
			report.n_scaler_reads_pre = selection.n_scaler_reads_pre;
			report.n_scaler_reads_post = selection.n_scaler_reads_post;
			report.hel0_charge_before_cut_uC = selection.hel0_charge_before_cut_uC;
			report.hel_neg_charge_before_cut_uC = selection.hel_neg_charge_before_cut_uC;
			report.hel_pos_charge_before_cut_uC = selection.hel_pos_charge_before_cut_uC;
			report.hel_charge_before_cut_uC = selection.hel_charge_before_cut_uC;
			report.hel0_charge_after_cut_uC = selection.hel0_charge_after_cut_uC;
			report.hel_neg_charge_after_cut_uC = selection.hel_neg_charge_after_cut_uC;
			report.hel_pos_charge_after_cut_uC = selection.hel_pos_charge_after_cut_uC;
			report.hel_charge_after_cut_uC = selection.hel_charge_after_cut_uC;
			report.evcount_cut = selection.evcount_cut;
			report.evnumber_cut = selection.evnumber_cut;
			report.gevnum_cut = selection.gevnum_cut;
			report.accepted_evcount_ranges = effstuff::ranges_to_string(selection.accepted_evcount_ranges);
			report.accepted_evnumber_ranges = effstuff::ranges_to_string(selection.accepted_evnumber_ranges);
			report.accepted_gevnum_ranges = effstuff::ranges_to_string(selection.accepted_gevnum_ranges);
			selection_report_rows->push_back(std::move(report));
		}

		if (!selection.ok) {
			++out.selection_failed_segments;
			continue;
		}

		accumulate_charge(out.hel_charge_before_cut_uC, selection.hel_charge_before_cut_uC);
		accumulate_charge(out.hel_neg_charge_before_cut_uC, selection.hel_neg_charge_before_cut_uC);
		accumulate_charge(out.hel_pos_charge_before_cut_uC, selection.hel_pos_charge_before_cut_uC);
		accumulate_charge(out.hel_charge_after_cut_uC, selection.hel_charge_after_cut_uC);
		accumulate_charge(out.hel_neg_charge_after_cut_uC, selection.hel_neg_charge_after_cut_uC);
		accumulate_charge(out.hel_pos_charge_after_cut_uC, selection.hel_pos_charge_after_cut_uC);

		effstuff::accumulate_beam_time_and_scaler(file_path, selection, beam_acc);
		accumulate_segment_histogram_and_metrics(
				file_path, selection, pid_cuts, edtm_hist, current_window_edtm_tdc,
				selected_current_window_edtm_tdc,
				pid_acc, cal_tag_acc, cer_tag_acc, tracking_acc, hodo_acc,
				out.missing_branch_segments);

		++n_segments_used;
	}
	out.n_segments = n_segments_used;

	out.negative_dt_intervals = beam_acc.negative_dt_intervals;
	out.non_monotonic_scaler_steps = beam_acc.non_monotonic_scaler_steps;
	out.missing_branch_segments += beam_acc.missing_branch_segments;
	out.beam_time = beam_acc.beam_time;
	out.NewGen_EDTM_den_evcount_gated_value = beam_acc.scaler_edtm_total_evcount_gated;
	out.NewGen_EDTM_den_evcount_gated = beam_acc.used_evcount_range_gate;

	if (out.n_segments == 0) {
		out.run_processing_status = "selection_failed";
		return out;
	}

	const double edtm_peak = compute_edtm_peak_from_hist(edtm_hist);
	out.NewGen_EDTM_peak = edtm_peak;
	if (!std::isfinite(edtm_peak)) {
		log_warn("  [warn run " + std::to_string(run_cfg.run_number) +
		         "] EDTM peak could not be determined from current-window EDTM histogram;"
		         " NewGen_EDTM numerator will remain 0 for this run");
	}

	const long long newgen_edtm_num = std::isfinite(edtm_peak)
			? count_edtm_peak_events(current_window_edtm_tdc, edtm_peak)
			: 0;
	const long long newgen_edtm_num_good_event_gated = std::isfinite(edtm_peak)
			? count_edtm_peak_events(selected_current_window_edtm_tdc, edtm_peak)
			: 0;
	out.NewGen_EDTM_num = newgen_edtm_num;
	out.NewGen_EDTM_den = beam_acc.scaler_edtm_total;
	out.NewGen_EDTM_num_good_event_gated = newgen_edtm_num_good_event_gated;

	out.HMS_pid_eff = effstuff::HMS_pid_eff(pid_acc);
	out.HMS_pid_eff_err = effstuff::binomial_efficiency_error(
			pid_acc.numerator_events,
			pid_acc.denominator_events);
	out.HMS_cal_eff_tag_cer = effstuff::HMS_cal_eff_tag_cer(cal_tag_acc);
	out.HMS_cal_eff_tag_cer_err = effstuff::binomial_efficiency_error(
			cal_tag_acc.numerator_events,
			cal_tag_acc.denominator_events);
	out.HMS_cer_eff_tag_cal = effstuff::HMS_cer_eff_tag_cal(cer_tag_acc);
	out.HMS_cer_eff_tag_cal_err = effstuff::binomial_efficiency_error(
			cer_tag_acc.numerator_events,
			cer_tag_acc.denominator_events);
	out.HMS_tracking_eff = effstuff::HMS_tracking_eff(tracking_acc);
	out.HMS_tracking_eff_err = effstuff::HMS_tracking_eff_err(tracking_acc);
	out.HMS_hodo_3of4_eff = effstuff::HMS_hodo_3of4_eff(hodo_acc);
	out.HMS_hodo_3of4_eff_err = effstuff::binomial_efficiency_error(
			hodo_acc.numerator_events,
			hodo_acc.denominator_events);

	// Warn if any efficiency denominator is zero (would produce NaN)
	if (pid_acc.denominator_events == 0) {
		log_warn("  [warn run " + std::to_string(run_cfg.run_number) +
		         "] HMS PID efficiency denominator is 0 (no events passed PID denominator cuts)");
	}
	if (tracking_acc.should_events == 0) {
		log_warn("  [warn run " + std::to_string(run_cfg.run_number) +
		         "] HMS tracking efficiency denominator is 0 (no events passed tracking should-cuts)");
	}
	if (hodo_acc.denominator_events == 0) {
		log_warn("  [warn run " + std::to_string(run_cfg.run_number) +
		         "] HMS hodo 3/4 efficiency denominator is 0 (no events passed hodo denominator cuts)");
	}

	effstuff::NewGenEdtmAccumulator livetime_acc;
	livetime_acc.numerator_events = newgen_edtm_num;
	livetime_acc.denominator_scaler_edtm = beam_acc.scaler_edtm_total;
	out.NewGen_EDTM_livetime = effstuff::NewGen_EDTM_livetime(livetime_acc, out.ps_factor);
	out.NewGen_EDTM_livetime_err = effstuff::NewGen_EDTM_livetime_err(livetime_acc, out.ps_factor);

	if (beam_acc.scaler_edtm_total <= 0.0) {
		log_warn("  [warn run " + std::to_string(run_cfg.run_number) +
		         "] EDTM scaler denominator is 0 (livetime will be 0 or NaN)");
	}

	{
		std::ostringstream dbg;
		dbg << "  [debug run " << run_cfg.run_number << "]"
		    << " source=" << out.file_source_used
		    << " seg_found=" << out.segment_count_found
		    << " seg_used=" << out.n_segments
		    << " edtm_peak=" << edtm_peak
		    << " sel_failed=" << out.selection_failed_segments
		    << " missing_branch_seg=" << out.missing_branch_segments
		    << " pid=" << pid_acc.numerator_events << "/" << pid_acc.denominator_events
		    << " cal_tag=" << cal_tag_acc.numerator_events << "/" << cal_tag_acc.denominator_events
		    << " cer_tag=" << cer_tag_acc.numerator_events << "/" << cer_tag_acc.denominator_events
		    << " tracking=" << tracking_acc.did_events << "/" << tracking_acc.should_events
		    << " hodo3of4=" << hodo_acc.numerator_events << "/" << hodo_acc.denominator_events
		    << " livetime_num=" << newgen_edtm_num
		    << " livetime_den=" << beam_acc.scaler_edtm_total
		    << " livetime_num_good_event_gated=" << newgen_edtm_num_good_event_gated
		    << " livetime_den_evcount_gated_value=" << beam_acc.scaler_edtm_total_evcount_gated
		    << " livetime_den_evcount_gated=" << (beam_acc.used_evcount_range_gate ? 1 : 0);
		log_debug(dbg.str());
	}

	if (out.missing_branch_segments > 0 || out.selection_failed_segments > 0) {
		out.run_processing_status = "processed_partial";
	} else {
		out.run_processing_status = "processed";
	}

	return out;
}

void write_kinematic_csv(const std::string& output_dir,
												 const std::string& kin,
												 const std::vector<effstuff::RunProcessingRow>& rows) {
	fs::create_directories(output_dir);
	const std::string filename = output_dir + "/efficiency_" + effstuff::sanitize_for_filename(kin) + ".csv";

	std::ofstream out(filename);
	if (!out) {
		throw std::runtime_error("Could not open output file: " + filename);
	}

	out << "run_number,kinematic_setting,run_type,"
			<< "HEL_charge_before_cut_uC,HEL_neg_charge_before_cut_uC,HEL_pos_charge_before_cut_uC,"
			<< "HEL_charge_after_cut_uC,HEL_neg_charge_after_cut_uC,HEL_pos_charge_after_cut_uC,"
			<< "HMS_pid_eff,HMS_cal_eff_tag_cer,HMS_cer_eff_tag_cal,HMS_tracking_eff,HMS_hodo_3of4_eff,"
			<< "NewGen_EDTM_livetime,NewGen_EDTM_num,NewGen_EDTM_den,NewGen_EDTM_peak,"
			<< "NewGen_EDTM_num_good_event_gated,NewGen_EDTM_den_evcount_gated_value,"
			<< "NewGen_EDTM_den_evcount_gated,prescale_token,ps_factor,which_TRIG,beam_time,"
			<< "file_source_used,segment_count_found,n_segments,run_processing_status,"
			<< "missing_branch_segments,selection_failed_segments,negative_dt_intervals,non_monotonic_scaler_steps,"
			<< "HMS_pid_eff_err,HMS_cal_eff_tag_cer_err,HMS_cer_eff_tag_cal_err,"
			<< "HMS_tracking_eff_err,HMS_hodo_3of4_eff_err,NewGen_EDTM_livetime_err\n";

	out << std::fixed << std::setprecision(8);
	for (const auto& row : rows) {
		out << row.run_number << ','
				<< effstuff::csv_quote(row.kinematic_setting) << ','
				<< effstuff::csv_quote(row.run_type) << ','
				<< row.hel_charge_before_cut_uC << ','
				<< row.hel_neg_charge_before_cut_uC << ','
				<< row.hel_pos_charge_before_cut_uC << ','
				<< row.hel_charge_after_cut_uC << ','
				<< row.hel_neg_charge_after_cut_uC << ','
				<< row.hel_pos_charge_after_cut_uC << ','
				<< row.HMS_pid_eff << ','
				<< row.HMS_cal_eff_tag_cer << ','
				<< row.HMS_cer_eff_tag_cal << ','
				<< row.HMS_tracking_eff << ','
				<< row.HMS_hodo_3of4_eff << ','
				<< row.NewGen_EDTM_livetime << ','
				<< row.NewGen_EDTM_num << ','
				<< row.NewGen_EDTM_den << ','
				<< row.NewGen_EDTM_peak << ','
				<< row.NewGen_EDTM_num_good_event_gated << ','
				<< row.NewGen_EDTM_den_evcount_gated_value << ','
				<< (row.NewGen_EDTM_den_evcount_gated ? 1 : 0) << ','
				<< effstuff::csv_quote(row.prescale_token) << ','
				<< row.ps_factor << ','
				<< effstuff::csv_quote(row.which_TRIG) << ','
				<< row.beam_time << ','
				<< effstuff::csv_quote(row.file_source_used) << ','
				<< row.segment_count_found << ','
				<< row.n_segments << ','
				<< effstuff::csv_quote(row.run_processing_status) << ','
				<< row.missing_branch_segments << ','
				<< row.selection_failed_segments << ','
				<< row.negative_dt_intervals << ','
				<< row.non_monotonic_scaler_steps << ','
				<< row.HMS_pid_eff_err << ','
				<< row.HMS_cal_eff_tag_cer_err << ','
				<< row.HMS_cer_eff_tag_cal_err << ','
				<< row.HMS_tracking_eff_err << ','
				<< row.HMS_hodo_3of4_eff_err << ','
				<< row.NewGen_EDTM_livetime_err << '\n';
	}

	std::cout << "[output] Wrote " << filename << " with " << rows.size() << " run rows\n";
}

void write_selection_report_csv(
		const std::string& output_dir,
		const std::string& kin,
		const std::vector<effstuff::SelectionReportRow>& rows) {
	fs::create_directories(output_dir);
	const std::string filename = output_dir + "/selection_report_" + effstuff::sanitize_for_filename(kin) + ".csv";

	std::ofstream out(filename);
	if (!out) {
		throw std::runtime_error("Could not open output file: " + filename);
	}

	out << "run_number,kinematic_setting,run_type,segment_number,segment_file,file_source_used,"
			<< "selection_ok,selection_message,helicity_mode,quartet_snap_applied,current_window_enabled,"
			<< "mean_current_uA,i0_used_uA,current_min_uA,current_max_uA,n_scaler_reads_pre,n_scaler_reads_post,"
			<< "hel0_charge_before_cut_uC,HEL_neg_charge_before_cut_uC,HEL_pos_charge_before_cut_uC,HEL_charge_before_cut_uC,"
			<< "hel0_charge_after_cut_uC,HEL_neg_charge_after_cut_uC,HEL_pos_charge_after_cut_uC,HEL_charge_after_cut_uC,"
			<< "evcount_cut,evnumber_cut,gevnum_cut,"
			<< "accepted_evcount_ranges,accepted_evnumber_ranges,accepted_gevnum_ranges\n";

	out << std::fixed << std::setprecision(8);
	for (const auto& row : rows) {
		out << row.run_number << ','
				<< effstuff::csv_quote(row.kinematic_setting) << ','
				<< effstuff::csv_quote(row.run_type) << ','
				<< row.segment_number << ','
				<< effstuff::csv_quote(row.segment_file) << ','
				<< effstuff::csv_quote(row.file_source_used) << ','
				<< (row.selection_ok ? 1 : 0) << ','
				<< effstuff::csv_quote(row.selection_message) << ','
				<< effstuff::csv_quote(row.helicity_mode) << ','
				<< (row.quartet_snap_applied ? 1 : 0) << ','
				<< (row.current_window_enabled ? 1 : 0) << ','
				<< row.mean_current_uA << ','
				<< row.i0_used_uA << ','
				<< row.current_min_uA << ','
				<< row.current_max_uA << ','
				<< row.n_scaler_reads_pre << ','
				<< row.n_scaler_reads_post << ','
				<< row.hel0_charge_before_cut_uC << ','
				<< row.hel_neg_charge_before_cut_uC << ','
				<< row.hel_pos_charge_before_cut_uC << ','
				<< row.hel_charge_before_cut_uC << ','
				<< row.hel0_charge_after_cut_uC << ','
				<< row.hel_neg_charge_after_cut_uC << ','
				<< row.hel_pos_charge_after_cut_uC << ','
				<< row.hel_charge_after_cut_uC << ','
				<< effstuff::csv_quote(row.evcount_cut) << ','
				<< effstuff::csv_quote(row.evnumber_cut) << ','
				<< effstuff::csv_quote(row.gevnum_cut) << ','
				<< effstuff::csv_quote(row.accepted_evcount_ranges) << ','
				<< effstuff::csv_quote(row.accepted_evnumber_ranges) << ','
				<< effstuff::csv_quote(row.accepted_gevnum_ranges) << '\n';
	}

	std::cout << "[output] Wrote " << filename << " with " << rows.size() << " segment rows\n";
}

void write_global_summary_csv(
		const std::string& output_dir,
		const std::vector<effstuff::KinematicProcessingSummary>& summaries) {
	fs::create_directories(output_dir);
	const std::string path = output_dir + "/efficiency_runs_processed.csv";

	std::ofstream out(path);
	if (!out) {
		throw std::runtime_error("Could not open summary output file: " + path);
	}

	out << "kinematic_setting,total_runs_listed,production_runs,runs_selected_by_type,"
			<< "runs_processed,runs_not_found,malformed_rows_skipped\n";
	for (const auto& summary : summaries) {
		out << effstuff::csv_quote(summary.kinematic_setting) << ','
				<< summary.total_runs_listed << ','
				<< summary.production_runs << ','
				<< summary.runs_selected_by_type << ','
				<< summary.runs_processed << ','
				<< summary.runs_not_found << ','
				<< summary.malformed_rows_skipped << '\n';
	}

	std::cout << "[output] Wrote " << path << "\n";
}

}  // namespace

int main(int argc, char** argv) {
	ProgramConfig cfg;
	if (!parse_args(argc, argv, cfg)) {
		return 1;
	}

	std::ofstream log_file;
	if (cfg.log_path.has_value()) {
		fs::path log_path = fs::path(*cfg.log_path);
		if (log_path.is_relative()) {
			log_path = fs::path(cfg.output_dir) / log_path;
		}

		std::error_code ec;
		if (log_path.has_parent_path()) {
			fs::create_directories(log_path.parent_path(), ec);
		}

		log_file.open(log_path.string(), std::ios::out | std::ios::trunc);
		if (!log_file) {
			std::cerr << "[warn] Could not open log file: " << log_path.string() << "\n";
		} else {
			g_log_stream = &log_file;
			log_debug("=== compute_efficiencies_stuff debug log ===");
			log_debug("config=" + cfg.config_csv_path);
			log_debug("updated_dir=" + cfg.updated_root_dir);
			log_debug("production_dir=" + cfg.production_root_dir);
			log_debug("output_dir=" + cfg.output_dir);
		}
	}

	try {
		effstuff::ConfigCsvData csv_data;
		std::string csv_error;
		if (!effstuff::load_config_csv(cfg.config_csv_path, csv_data, csv_error)) {
			std::cerr << "[error] " << csv_error << "\n";
			return 2;
		}

		const std::vector<std::string> all_kinematics = effstuff::collect_unique_kinematics(csv_data.rows);
		const std::vector<std::string> selected_kinematics = choose_kinematics(cfg, all_kinematics);
		if (selected_kinematics.empty()) {
			std::cerr << "[error] No kinematic settings selected.\n";
			return 3;
		}

		std::cout << "[info] Selected " << selected_kinematics.size() << " kinematic setting(s).\n";
		std::cout << "[info] Allowed Type filters:";
		for (const auto& t : cfg.allowed_types) {
			std::cout << " " << t;
		}
		std::cout << "\n";
		if (g_log_stream) {
			std::ostringstream allowed;
			allowed << "allowed_types=";
			for (size_t i = 0; i < cfg.allowed_types.size(); ++i) {
				if (i > 0) {
					allowed << ",";
				}
				allowed << cfg.allowed_types[i];
			}
			log_debug(allowed.str());
		}

		std::map<std::string, std::vector<effstuff::RunConfigRow>> rows_by_kin;
		for (const auto& row : csv_data.rows) {
			rows_by_kin[row.kin_old].push_back(row);
		}

		const SelectionSettings selection_settings = effstuff::make_default_selection_settings();

		std::vector<effstuff::KinematicProcessingSummary> summaries;
		summaries.reserve(selected_kinematics.size());

		for (const auto& kin : selected_kinematics) {
			log_info("[kinematic] Processing " + kin);
			effstuff::KinematicProcessingSummary summary;
			summary.kinematic_setting = kin;
			summary.malformed_rows_skipped = csv_data.malformed_row_count;

			const auto kin_it = rows_by_kin.find(kin);
			if (kin_it == rows_by_kin.end()) {
				summaries.push_back(summary);
				continue;
			}

			const auto& kin_rows = kin_it->second;

			{
				std::set<int> total_runs;
				std::set<int> production_runs;
				for (const auto& row : kin_rows) {
					total_runs.insert(row.run_number);
					if (row.run_type == "production" || row.run_type == "Production") {
						production_runs.insert(row.run_number);
					}
				}
				summary.total_runs_listed = static_cast<int>(total_runs.size());
				summary.production_runs = static_cast<int>(production_runs.size());
			}

			std::map<int, effstuff::RunConfigRow> selected_run_map;
			for (const auto& row : kin_rows) {
				if (!effstuff::type_is_allowed(row.run_type, cfg.allowed_types)) {
					continue;
				}
				if (cfg.single_run_filter.has_value() && row.run_number != cfg.single_run_filter.value()) {
					continue;
				}

				auto it = selected_run_map.find(row.run_number);
				if (it == selected_run_map.end()) {
					selected_run_map[row.run_number] = row;
				} else if (it->second.prescale_token.empty() && !row.prescale_token.empty()) {
					it->second.prescale_token = row.prescale_token;
				}
			}

			summary.runs_selected_by_type = static_cast<int>(selected_run_map.size());

			std::vector<effstuff::RunProcessingRow> run_rows;
			run_rows.reserve(selected_run_map.size());
			std::vector<effstuff::SelectionReportRow> selection_rows;

			for (const auto& kv : selected_run_map) {
				const auto& run_cfg = kv.second;
				log_info("  [run " + std::to_string(run_cfg.run_number) + "] start");

				auto row = process_run(run_cfg, cfg, selection_settings, &selection_rows);
				if (row.run_processing_status == "processed" || row.run_processing_status == "processed_partial") {
					++summary.runs_processed;
				}
				if (row.run_processing_status == "missing_root_files") {
					++summary.runs_not_found;
				}

				log_info("  [run " + std::to_string(run_cfg.run_number) + "] status=" +
				         row.run_processing_status + ", source=" + row.file_source_used +
				         ", segments=" + std::to_string(row.segment_count_found));
				run_rows.push_back(std::move(row));
			}

			write_kinematic_csv(cfg.output_dir, kin, run_rows);
			write_selection_report_csv(cfg.output_dir, kin, selection_rows);
			summaries.push_back(summary);
		}

		write_global_summary_csv(cfg.output_dir, summaries);
	} catch (const std::exception& ex) {
		if (g_log_stream) {
			log_debug(std::string("[fatal] ") + ex.what());
		}
		std::cerr << "[fatal] " << ex.what() << "\n";
		return 4;
	}

	g_log_stream = nullptr;

	return 0;
}
