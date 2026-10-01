// Compare six run-level livetime definitions on one good-event selection.
// Build: csh -c 'source /group/nps/singhav/setup.csh; g++ -std=c++17 -O2 -Wall -Wextra `root-config --cflags` compare_livetimes.cxx -o compare_livetimes `root-config --libs`'
// Run: ./compare_livetimes RUN [segment0.root ...]
// Omit files to discover updated replay segments, falling back to production.
// Parallel launcher selects kinematics, targets and runs; workers write private CSVs.
// Prescale metadata selects both branches: ps4 -> hTRIG4, ps6 -> hTRIG6.
// Input files must be ordered, non-overlapping segments from ONE replay version.
// CSV appends one row per run to output/efficiency_stuff/compare_livetimes_<Kin_old>.csv.

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "config_csv_helper.h"
#include "good_event_selection_helper.h"
#include "prescale_beamtime_helper.h"
#include "root_file_discovery.h"

namespace {

struct Interval {
  long long lo = 0, hi = 0;  // TSH.evNumber: [lo, hi) latch convention.
  bool selected = false;     // Entire interval has recorded-event coverage and good current.
};

struct Pulse {
  double raw = 0.0, corrected = 0.0;
};
static_assert(sizeof(Pulse) == 16, "Pulse layout must remain memory-compact");

using Histogram = std::map<long long, long long>;

struct Counts {
  double d_all = 0, d_current = 0, d_selected = 0;  // H.EDTM.scaler increments.
  double s_all = 0, s_selected = 0;                 // H.hTRIGn.scaler increments.
  double a_selected = 0, a_all = 0;      // H.hL1ACCP.scaler increments.
  long long edtm_tdc = 0, trig_tdc = 0, physics_tdc = 0;
  long long newgen_num = 0, matched_num = 0;
  std::vector<Pulse> pulses;             // Good events with current and a raw EDTM pulse.
  std::vector<bool> pulse_matched;        // One bit/pulse: belongs to a selected TSH interval.
  Histogram raw_bins;                    // Compact 10-channel histogram; no raw-value copy.
  long long selected_intervals = 0, selected_events = 0;
};

double ratio(double n, double d, double p) {
  return d > 0 ? p * n / d : std::numeric_limits<double>::quiet_NaN();
}

void fill_bin(Histogram& bins, double value, double width) {
  if (std::isfinite(value))
    ++bins[static_cast<long long>(std::floor(value / width))];
}

// Histogram maximum; map order gives a reproducible lower-bin tie break.
double peak(const Histogram& bins, double width) {
  long long best_count = 0;
  double best = std::numeric_limits<double>::quiet_NaN();
  for (const auto& [bin, count] : bins)
    if (count > best_count) { best_count = count; best = (bin + 0.5) * width; }
  return best;
}

// Selected event-number ranges use half-open bounds, as in the good-event helper.
bool whole_interval_in_good_range(long long lo, long long hi,
                                  const effstuff::GoodSelectionSummary& sel) {
  for (const auto& range : sel.accepted_evnumber_ranges)
    if (lo >= range.lo && hi <= range.hi) return true;
  return false;
}

// Read exactly the scaler quantities present in slide 3, plus evNumber for pairing.
std::vector<Interval> read_scalers(TTree* t, const std::string& trigger,
                                   const effstuff::GoodSelectionSummary& sel,
                                   Counts& c, long long first_event, long long last_event) {
  for (const auto& name : {"evNumber", "H.1MHz.scalerTime", "H.BCM4A.scalerCurrent",
                           "H.EDTM.scaler", "H.hL1ACCP.scaler"})
    if (!effstuff::tree_has_branch(t, name)) throw std::runtime_error(std::string("Missing TSH branch: ") + name);
  if (!effstuff::tree_has_branch(t, trigger.c_str())) throw std::runtime_error("Missing TSH trigger: " + trigger);

  double evnum = 0, clock = 0, current = 0, d = 0, s = 0, a = 0;
  t->SetBranchStatus("*", 0);
  for (const auto& name : {"evNumber", "H.1MHz.scalerTime", "H.BCM4A.scalerCurrent",
                           "H.EDTM.scaler", "H.hL1ACCP.scaler"}) t->SetBranchStatus(name, 1);
  t->SetBranchStatus(trigger.c_str(), 1);
  t->SetBranchAddress("evNumber", &evnum);
  t->SetBranchAddress("H.1MHz.scalerTime", &clock);
  t->SetBranchAddress("H.BCM4A.scalerCurrent", &current);
  t->SetBranchAddress("H.EDTM.scaler", &d);
  t->SetBranchAddress("H.hL1ACCP.scaler", &a);
  t->SetBranchAddress(trigger.c_str(), &s);

  std::vector<Interval> intervals;
  if (t->GetEntries() < 2) return intervals;
  t->GetEntry(0);
  double old_evnum = evnum, old_clock = clock, old_d = d, old_s = s, old_a = a;
  for (Long64_t i = 1; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    // A reset or backwards event boundary makes paired counts ambiguous.
    if (clock < old_clock || evnum < old_evnum || d < old_d || s < old_s || a < old_a)
      throw std::runtime_error("Non-monotonic TSH snapshot; inspect this segment");
    const double dd = d - old_d, ds = s - old_s, da = a - old_a;
    c.d_all += dd;
    c.s_all += ds;
    c.a_all += da;
    const bool good_current = effstuff::current_in_selection_window(current, sel);
    if (good_current) c.d_current += dd;  // Current NewGen denominator convention.

    const long long lo = std::llround(old_evnum), hi = std::llround(evnum);
    // Do not use initial/final scaler exposure without corresponding recorded events.
    const bool covered = lo >= first_event && hi <= last_event + 1 && hi > lo;
    const bool good = covered && good_current && whole_interval_in_good_range(lo, hi, sel);
    intervals.push_back({lo, hi, good});
    if (good) {
      c.d_selected += dd; c.s_selected += ds; c.a_selected += da;
      ++c.selected_intervals;
    }
    old_evnum = evnum; old_clock = clock; old_d = d; old_s = s; old_a = a;
  }
  return intervals;
}

void read_events(TTree* t, int trig_number, const effstuff::GoodSelectionSummary& sel,
                 const std::vector<Interval>& intervals, Counts& c) {
  const std::string trig = "T.hms.hTRIG" + std::to_string(trig_number) + "_tdcTimeRaw";
  for (const auto& name : {"g.evnum", "fEvtHdr.fEvtType", "H.BCM4A.scalerCurrent",
                           "T.hms.hEDTM_tdcTimeRaw", "T.hms.hEDTM_tdcTime"})
    if (!effstuff::tree_has_branch(t, name)) throw std::runtime_error(std::string("Missing T branch: ") + name);
  if (!effstuff::tree_has_branch(t, trig.c_str())) throw std::runtime_error("Missing T trigger: " + trig);

  double evnum = 0, current = 0, edtm_raw = 0, edtm_time = 0, trig_raw = 0;
  unsigned int event_type = 0;  // fEvtHdr.fEvtType is UInt_t in production replay.
  t->SetBranchStatus("*", 0);
  for (const auto& name : {"g.evnum", "fEvtHdr.fEvtType", "H.BCM4A.scalerCurrent",
                           "T.hms.hEDTM_tdcTimeRaw", "T.hms.hEDTM_tdcTime"}) t->SetBranchStatus(name, 1);
  t->SetBranchStatus(trig.c_str(), 1);
  t->SetBranchAddress("g.evnum", &evnum);
  t->SetBranchAddress("fEvtHdr.fEvtType", &event_type);
  t->SetBranchAddress("H.BCM4A.scalerCurrent", &current);
  t->SetBranchAddress("T.hms.hEDTM_tdcTimeRaw", &edtm_raw);
  t->SetBranchAddress("T.hms.hEDTM_tdcTime", &edtm_time);
  t->SetBranchAddress(trig.c_str(), &trig_raw);

  size_t j = 0;  // T events and TSH.evNumber boundaries are monotonic within segment.
  long long previous = std::numeric_limits<long long>::min();
  for (Long64_t i = 0; i < t->GetEntries(); ++i) {
    t->GetEntry(i);
    const long long number = std::llround(evnum);
    if (number <= previous) throw std::runtime_error("T event numbers are not strictly increasing");
    previous = number;
    while (j < intervals.size() && number >= intervals[j].hi) ++j;
    const bool matched = j < intervals.size() && number >= intervals[j].lo && intervals[j].selected;
    // Canonical good-event gate: physics event type and selected g.evnum range.
    const bool good_event = event_type == 1 &&
        effstuff::event_value_in_ranges(number, sel.accepted_gevnum_ranges);
    if (!good_event) continue;
    const bool current_good = effstuff::current_in_selection_window(current, sel);
    if (current_good && std::isfinite(edtm_raw) && edtm_raw > 1) {
      c.pulses.push_back({edtm_raw, edtm_time});
      c.pulse_matched.push_back(matched);
      fill_bin(c.raw_bins, edtm_raw, 10.0);
    }
    if (!matched) continue;
    ++c.selected_events;
    const bool edtm = std::isfinite(edtm_raw) && edtm_raw != 0;
    const bool trigger = std::isfinite(trig_raw) && trig_raw != 0;
    if (edtm) ++c.edtm_tdc;
    if (trigger) ++c.trig_tdc;
    if (trigger && !edtm) ++c.physics_tdc;
  }
}

}  // namespace

int main(int argc, char** argv) {
  try {
    if (argc < 2) {
      std::cerr << "Usage: " << argv[0]
                << " RUN [segment.root ...] [--config CSV] [--output-dir DIR]"
                   " [--updated-dir DIR] [--production-dir DIR]"
                   " [--file-source updated|production|explicit]\n";
      return 2;
    }
    const int run = std::stoi(argv[1]);
    // Defaults mirror compute_efficiencies_stuff.cxx; options support isolated workers.
    std::string config_path =
        "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_final/config/nps_dvcs_all_kins_main.csv";
    std::string output_dir =
        "/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_final/output/efficiency_stuff";
    std::string updated_dir = "/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/updated";
    std::string production_dir = "/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/production";
    std::string explicit_source;
    std::vector<std::string> files;
    for (int i = 2; i < argc; ++i) {
      const std::string arg = argv[i];
      if (arg == "--config" || arg == "--output-dir" || arg == "--updated-dir" ||
          arg == "--production-dir" || arg == "--file-source") {
        if (++i >= argc) throw std::runtime_error("Missing value after " + arg);
        const std::string value = argv[i];
        if (arg == "--config") config_path = value;
        else if (arg == "--output-dir") output_dir = value;
        else if (arg == "--updated-dir") updated_dir = value;
        else if (arg == "--production-dir") production_dir = value;
        else explicit_source = value;
      } else if (arg.rfind("--", 0) == 0) {
        throw std::runtime_error("Unknown option: " + arg);
      } else {
        files.push_back(arg);
      }
    }

    // Kin_old, target and prescale come from one metadata row per run.
    effstuff::ConfigCsvData config;
    std::string config_error;
    if (!effstuff::load_config_csv(config_path, config, config_error))
      throw std::runtime_error(config_error);
    std::string kin, target, run_type, prescale_token;
    bool found_run = false;
    for (const auto& row : config.rows) {
      if (row.run_number != run) continue;
      if (found_run && kin != row.kin_old)
        throw std::runtime_error("Run has conflicting Kin_old labels in config");
      if (found_run && target != row.target)
        throw std::runtime_error("Run has conflicting target labels in config");
      if (found_run && prescale_token != row.prescale_token)
        throw std::runtime_error("Run has conflicting prescale tokens in config");
      kin = row.kin_old;
      target = row.target;
      run_type = row.run_type;
      prescale_token = row.prescale_token;
      found_run = true;
    }
    if (!found_run) throw std::runtime_error("Run not found in kinematic config");
    if (target.empty()) throw std::runtime_error("Run has no target in config");
    if (prescale_token.empty()) throw std::runtime_error("Run has no prescale token in config");
    const auto ps = effstuff::build_prescale_info(prescale_token);
    if (!ps.valid || ps.multiple_enabled)
      throw std::runtime_error("Run must have exactly one enabled prescale trigger: " + prescale_token);
    // ps.trig_number selects T.hms.hTRIGn_tdcTimeRaw; ps.which_TRIG selects H.hTRIGn.scaler.
    if (kin.find_first_not_of("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-") != std::string::npos)
      throw std::runtime_error("Kin_old contains unsafe filename characters: " + kin);
    if (!explicit_source.empty() && explicit_source != "updated" &&
        explicit_source != "production" && explicit_source != "explicit")
      throw std::runtime_error("--file-source must be updated, production, or explicit");
    std::string file_source = explicit_source.empty() ? "explicit" : explicit_source;
    if (files.empty()) {
      if (!explicit_source.empty())
        throw std::runtime_error("--file-source requires explicit segment files");
      const auto located = effstuff::locate_run_files_prefer_updated(run, updated_dir, production_dir);
      files = located.files;
      file_source = located.source;
    }
    if (files.empty()) throw std::runtime_error("No replay segments found for run " + std::to_string(run));
    const std::string output_path = std::string(output_dir) + "/compare_livetimes_" + kin + ".csv";

    Counts c;
    const auto settings = effstuff::make_default_selection_settings();
    long long previous_last = std::numeric_limits<long long>::min();
    for (const auto& path : files) {
      const auto sel = effstuff::build_good_selection_summary(path, settings);
      if (!sel.ok || sel.accepted_evnumber_ranges.empty() || sel.accepted_gevnum_ranges.empty())
        throw std::runtime_error("Good-event selection failed for " + path + ": " + sel.message);
      TFile file(path.c_str(), "READ");
      if (file.IsZombie()) throw std::runtime_error("Cannot open " + path);
      auto* te = dynamic_cast<TTree*>(file.Get("T"));
      auto* ts = dynamic_cast<TTree*>(file.Get("TSH"));
      if (!te || !ts || te->GetEntries() == 0) throw std::runtime_error("Missing/empty T or TSH in " + path);
      if (!effstuff::tree_has_branch(te, "g.evnum")) throw std::runtime_error("Missing T.g.evnum in " + path);
      double number = 0;
      te->SetBranchAddress("g.evnum", &number);
      te->GetEntry(0); const long long first = std::llround(number);
      te->GetEntry(te->GetEntries() - 1); const long long last = std::llround(number);
      te->ResetBranchAddresses();
      if (first <= previous_last || last < first)
        throw std::runtime_error("Input segments overlap or are out of event-number order: " + path);
      previous_last = last;
      const auto intervals = read_scalers(ts, ps.which_TRIG, sel, c, first, last);
      read_events(te, ps.trig_number, sel, intervals, c);
    }
    if (c.selected_intervals == 0) throw std::runtime_error("No good, recorded TSH intervals");

    // Current NewGen: 10-channel raw histogram maximum and +/-500 raw channels.
    const double raw_peak = peak(c.raw_bins, 10.0);
    Histogram corrected_bins;
    for (const auto& pulse : c.pulses) {
      if (std::abs(pulse.raw - raw_peak) > 500) continue;
      ++c.newgen_num;
      fill_bin(corrected_bins, pulse.corrected, 0.1);
    }
    // Matched EDTM: dense 0.1-ns corrected peak, median within 1 ns, +/-2 ns.
    const double rough = peak(corrected_bins, 0.1);
    std::vector<double> near;
    for (const auto& pulse : c.pulses)
      if (std::abs(pulse.raw - raw_peak) <= 500 &&
          std::isfinite(pulse.corrected) && std::abs(pulse.corrected - rough) <= 1)
        near.push_back(pulse.corrected);
    std::sort(near.begin(), near.end());
    const double corrected_peak = near.empty() ? std::numeric_limits<double>::quiet_NaN()
        : (near[near.size()/2] + near[(near.size()-1)/2]) / 2;
    for (size_t i = 0; i < c.pulses.size(); ++i)
      if (c.pulse_matched[i] && std::isfinite(c.pulses[i].corrected) &&
          std::abs(c.pulses[i].corrected - corrected_peak) <= 2) ++c.matched_num;

    // Private worker directories prevent concurrent appends to the same CSV.
    std::filesystem::create_directories(output_dir);
    const bool write_header = !std::filesystem::exists(output_path) ||
        std::filesystem::file_size(output_path) == 0;
    const std::string header =
        "run,kin,target,type,file_source,segments,trigger,prescale_factor,selected_intervals,selected_events,"
        "scaler_edtm_all,scaler_edtm_current,scaler_edtm_selected,"
        "scaler_htrig_all,scaler_htrig_selected,"
        "scaler_l1accp_selected,scaler_l1accp_all,beam_on_fraction_l1accp,"
        "tdc_edtm_selected,tdc_htrig_selected,tdc_physics_selected,"
        "newgen_raw_peak_channels,newgen_num_good_current,"
        "matched_corrected_peak_ns,matched_edtm_num,"
        "total_edtm_lt,clta_tsh_lt,clta_tdc_lt,cltp_tdc_lt,newgen_edtm_lt,matched_edtm_lt";
    if (!write_header) {
      std::ifstream existing(output_path);
      std::string old_header;
      std::getline(existing, old_header);
      if (old_header != header) throw std::runtime_error("Existing CSV has different columns: " + output_path);
      std::string line;
      while (std::getline(existing, line))
        if (line.substr(0, line.find(',')) == std::to_string(run))
          throw std::runtime_error("Run already exists in CSV: " + output_path);
    }
    std::ofstream out(output_path, std::ios::app);
    if (!out) throw std::runtime_error("Cannot write CSV: " + output_path);
    out << std::setprecision(12);
    if (write_header) out << header << '\n';
    // Direct interval selection replaces slide 3's proxy beam-on multiplier.
    // It remains in CSV so every scaler value and correction is auditable.
    out << run << ',' << effstuff::csv_quote(kin) << ',' << effstuff::csv_quote(target)
        << ',' << effstuff::csv_quote(run_type) << ',' << file_source << ',' << files.size()
        << ',' << ps.trig_number << ',' << ps.ps_factor << ','
        << c.selected_intervals << ',' << c.selected_events << ',' << c.d_all << ','
        << c.d_current << ',' << c.d_selected << ',' << c.s_all << ',' << c.s_selected << ','
        << c.a_selected << ',' << c.a_all << ','
        << (c.a_all > 0 ? c.a_selected/c.a_all : std::numeric_limits<double>::quiet_NaN()) << ','
        << c.edtm_tdc << ',' << c.trig_tdc << ',' << c.physics_tdc << ','
        << raw_peak << ',' << c.newgen_num << ',' << corrected_peak << ',' << c.matched_num << ','
        << ratio(c.edtm_tdc, c.d_selected, ps.ps_factor) << ','
        << ratio(c.a_selected, c.s_selected, ps.ps_factor) << ','
        << ratio(c.trig_tdc, c.s_selected, ps.ps_factor) << ','
        << ratio(c.physics_tdc, c.s_selected-c.d_selected, ps.ps_factor) << ','
        << ratio(c.newgen_num, c.d_current, ps.ps_factor) << ','
        << ratio(c.matched_num, c.d_selected, ps.ps_factor) << '\n';
    std::cout << "Wrote " << output_path << "\n";
  } catch (const std::exception& e) {
    std::cerr << "Error: " << e.what() << '\n';
    return 1;
  }
}
