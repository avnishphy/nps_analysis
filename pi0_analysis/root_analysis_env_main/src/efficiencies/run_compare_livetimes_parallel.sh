#!/usr/bin/env bash
# Run one ROOT process per run, then merge private worker CSVs by Kin_old.
# First load ROOT: csh -c 'source /group/nps/singhav/setup.csh; ./run_compare_livetimes_parallel.sh --kin KinC_x60_4b --target LH2 --jobs 4'
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exe="${here}/compare_livetimes"
config="/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/config/nps_dvcs_all_kins_main.csv"
updated="/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/updated"
production="/lustre24/expphy/cache/hallc/c-nps/analysis/pass2/replays/production"
output="/w/hallc-scshelf2102/nps/singhav/nps_analysis/pi0_analysis/root_analysis_env_main/output/efficiency_stuff"
plots_dir=""              # Defaults to <output>/plots after --output-dir is parsed.
jobs="$(nproc)"
retries=2                  # Retry only workers killed by SIGKILL/OOM (status 137).
types="production,Production"
kins=() targets=() runs=()
all_kins=0 dry_run=0 no_plots=0

usage() {
  cat <<EOF
Usage: $(basename "$0") [--kin Kin_old ... | --all-kins] [--target name ...] [options]
  --kin NAME           Repeat to select several Kin_old settings (default: all)
  --all-kins           Explicitly select all settings
  --target NAME        Repeat to select several targets (default: all)
  --types a,b          Type filter (default: production,Production)
  --run N              Repeat to restrict runs
  --jobs N             Parallel run workers (default: $(nproc))
  --retries N          Retry status-137 workers N times (default: 2)
  --config FILE        Run metadata CSV
  --updated-dir DIR    Updated replay directory
  --production-dir DIR Production replay directory
  --output-dir DIR     Final CSV directory
  --plots-dir DIR      Combined PDF directory (default: <output-dir>/plots)
  --exe FILE           C++ executable path
  --dry-run            Print selected jobs; read no ROOT files
  --no-plots           Skip overlay plots
EOF
}

# Match existing efficiency launcher options without its interactive prompt.
while (($#)); do
  case "$1" in
    --kin) kins+=("$2"); shift 2 ;;
    --all-kins) all_kins=1; shift ;;
    --target) targets+=("$2"); shift 2 ;;
    --types) types="$2"; shift 2 ;;
    --run) runs+=("$2"); shift 2 ;;
    --jobs) jobs="$2"; shift 2 ;;
    --retries) retries="$2"; shift 2 ;;
    --config) config="$2"; shift 2 ;;
    --updated-dir) updated="$2"; shift 2 ;;
    --production-dir) production="$2"; shift 2 ;;
    --output-dir) output="$2"; shift 2 ;;
    --plots-dir) plots_dir="$2"; shift 2 ;;
    --exe) exe="$2"; shift 2 ;;
    --dry-run) dry_run=1; shift ;;
    --no-plots) no_plots=1; shift ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown option: $1" >&2; usage >&2; exit 2 ;;
  esac
done
if ((all_kins && ${#kins[@]})); then echo "Use --kin or --all-kins" >&2; exit 2; fi
if [[ ! "$jobs" =~ ^[1-9][0-9]*$ ]]; then echo "--jobs must be positive" >&2; exit 2; fi
if [[ ! "$retries" =~ ^[0-9]+$ ]]; then echo "--retries must be a non-negative integer" >&2; exit 2; fi
[[ -f "$config" ]] || { echo "Missing config: $config" >&2; exit 2; }

tmp="$(mktemp -d)"
trap 'rm -rf -- "$tmp"' EXIT  # Only the directory created by this invocation.
kin_csv="$(IFS=,; echo "${kins[*]}")"
target_csv="$(IFS=,; echo "${targets[*]}")"
run_csv="$(IFS=,; echo "${runs[*]}")"

# Python's CSV parser handles quoted fields in the master run metadata.
python3 - "$config" "$kin_csv" "$target_csv" "$types" "$run_csv" > "$tmp/jobs.tsv" <<'PY'
import csv, sys
path, kin_arg, target_arg, type_arg, run_arg = sys.argv[1:]
selected_kins = set(filter(None, kin_arg.split(',')))
selected_targets = {x.casefold() for x in filter(None, target_arg.split(','))}
selected_types = {x.casefold() for x in filter(None, type_arg.split(','))}
selected_runs = set(filter(None, run_arg.split(',')))
jobs = {}
with open(path, newline='') as f:
    reader = csv.DictReader(f)
    required = {'run_number', 'Kin_old', 'target', 'Type', 'prescale'}
    if not required.issubset({x.strip() for x in reader.fieldnames or []}):
        raise SystemExit('Config needs run_number, Kin_old, target, Type, prescale')
    for raw in reader:
        row = {k.strip(): (v or '').strip() for k, v in raw.items() if k}
        run, kin, target, typ = (row[k] for k in ('run_number', 'Kin_old', 'target', 'Type'))
        if not run.isdigit() or not kin or not target: continue
        if selected_kins and kin not in selected_kins: continue
        if selected_targets and target.casefold() not in selected_targets: continue
        if typ.casefold() not in selected_types: continue
        if selected_runs and run not in selected_runs: continue
        identity = (kin, target, row['prescale'])
        if run in jobs and jobs[run] != identity:
            raise SystemExit(f'Conflicting metadata for run {run}')
        jobs[run] = identity
for run in sorted(jobs, key=int):
    kin, target, _ = jobs[run]
    if any('\t' in x or '\n' in x for x in (kin, target)):
        raise SystemExit(f'Unsupported tab/newline in metadata for run {run}')
    print(f'{run}\t{kin}\t{target}')
PY

total="$(wc -l < "$tmp/jobs.tsv")"
if ((total == 0)); then echo "No runs match filters" >&2; exit 1; fi
if ((dry_run)); then cat "$tmp/jobs.tsv"; echo "Selected $total runs"; exit 0; fi
if ((jobs > total)); then jobs="$total"; fi

# Workers need the Hall C ROOT environment loaded before build or execution.
command -v root-config >/dev/null || {
  echo "ROOT unavailable; first source /group/nps/singhav/setup.csh in csh" >&2; exit 2;
}
source_file="${here}/compare_livetimes.cxx"
if [[ ! -x "$exe" || "$source_file" -nt "$exe" ]] ||
   find "$here" -maxdepth 1 -name '*.h' -newer "$exe" -print -quit | grep -q .; then
  echo "[build] $exe"
  g++ -std=c++17 -O2 -Wall -Wextra "$source_file" -o "$exe" $(root-config --cflags --libs)
fi
mkdir -p "$output/compare_livetimes_logs"
export exe config updated production output tmp
export MALLOC_ARENA_MAX="${MALLOC_ARENA_MAX:-2}"  # Limit allocator overhead per ROOT process.
echo "[parallel] $total runs; $jobs workers"

# Each subprocess gets a private CSV directory. A failed run cannot corrupt a final CSV.
cut -f1 "$tmp/jobs.tsv" | xargs -P "$jobs" -I {} bash -c '
  run="$1"
  dir="$tmp/worker_$run"
  mkdir -p "$dir"
  log="$output/compare_livetimes_logs/run${run}.log"
  status=0
  "$exe" "$run" --config "$config" --output-dir "$dir" \
    --updated-dir "$updated" --production-dir "$production" > "$log" 2>&1 || status=$?
  printf "%s\n" "$status" > "$tmp/status_$run"
  exit "$status"
' _ {} || true

# A memory-pressure SIGKILL returns 137. Retry those runs after the main wave,
# when completed workers have released their memory. Other failures are not retried.
for ((attempt = 1; attempt <= retries; ++attempt)); do
  retry_file="$tmp/retry_$attempt.txt"
  : > "$retry_file"
  while IFS=$'\t' read -r run _; do
    status_file="$tmp/status_$run"
    [[ -f "$status_file" && "$(<"$status_file")" == 137 ]] && printf '%s\n' "$run" >> "$retry_file"
  done < "$tmp/jobs.tsv"
  retry_total="$(wc -l < "$retry_file")"
  ((retry_total)) || break
  echo "[retry $attempt/$retries] $retry_total status-137 runs; up to $jobs workers"
  xargs -P "$jobs" -I {} bash -c '
    run="$1"
    dir="$tmp/worker_$run"
    rm -f "$dir"/compare_livetimes_*.csv
    log="$output/compare_livetimes_logs/run${run}.log"
    status=0
    "$exe" "$run" --config "$config" --output-dir "$dir" \
      --updated-dir "$updated" --production-dir "$production" > "$log" 2>&1 || status=$?
    printf "%s\n" "$status" > "$tmp/status_$run"
    exit "$status"
  ' _ {} < "$retry_file" || true
done

# Merge once after every worker exits; keep unselected and previously successful runs.
merge_status=0
python3 - "$tmp" "$output" <<'PY' || merge_status=$?
import csv, os, pathlib, sys, tempfile
tmp, output = map(pathlib.Path, sys.argv[1:])
groups, failed = {}, []
with (tmp/'jobs.tsv').open() as f:
    for run, kin, target in csv.reader(f, delimiter='\t'):
        status_file = tmp/f'status_{run}'
        status = status_file.read_text().strip() if status_file.exists() else 'not_started'
        parts = list((tmp/f'worker_{run}').glob('compare_livetimes_*.csv'))
        if status != '0' or len(parts) != 1:
            failed.append((run, kin, target, status)); continue
        lines = parts[0].read_text().splitlines()
        if len(lines) != 2:
            failed.append((run, kin, target, 'bad_csv')); continue
        header, row = lines
        if kin in groups and groups[kin][0] != header:
            raise SystemExit(f'CSV header mismatch for {kin}')
        groups.setdefault(kin, (header, []))[1].append((int(run), row))
for kin, (header, rows) in groups.items():
    if any(x not in 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_-' for x in kin):
        raise SystemExit(f'Unsafe Kin_old filename: {kin}')
    path = output/f'compare_livetimes_{kin}.csv'
    by_run = {}
    if path.exists():
        old = path.read_text().splitlines()
        if old and old[0] != header:
            raise SystemExit(f'Existing CSV header mismatch: {path}')
        for line in old[1:]:
            by_run[int(line.split(',', 1)[0])] = line
    for run, row in rows: by_run[run] = row  # Recomputed rows replace their old versions.
    with tempfile.NamedTemporaryFile(mode='w', dir=output, prefix='.compare_livetimes_', delete=False) as f:
        f.write(header+'\n')
        for run in sorted(by_run): f.write(by_run[run]+'\n')
        pending = f.name
    os.replace(pending, path)
    print(f'[merge] {path}: {len(by_run)} total runs ({len(rows)} updated)')
failure_path = output/'compare_livetimes_failed_jobs.csv'
with failure_path.open('w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(('run', 'kin', 'target', 'status', 'log'))
    for run, kin, target, status in failed:
        writer.writerow((run, kin, target, status, str(output/f'compare_livetimes_logs/run{run}.log')))
print(f'[done] {sum(len(v[1]) for v in groups.values())} succeeded; {len(failed)} failed')
if failed: sys.exit(1)
PY

# Plot merged data only after all writers have exited. Generate one multipage PDF per target.
plot_status=0
if (( !no_plots )); then
  plots_dir="${plots_dir:-$output/plots}"
  csv_files=()
  while IFS= read -r kin; do
    csv="$output/compare_livetimes_${kin}.csv"
    [[ -f "$csv" ]] && csv_files+=("$csv")
  done < <(cut -f2 "$tmp/jobs.tsv" | sort -u)
  if ((${#csv_files[@]})); then
    export MPLCONFIGDIR="${MPLCONFIGDIR:-$tmp/matplotlib}"  # Writable font cache in batch shells.
    python3 "$here/plot_compare_livetimes.py" --output-dir "$plots_dir" \
      "${csv_files[@]}" || plot_status=$?
  fi
fi
if ((merge_status)); then exit "$merge_status"; fi
if ((plot_status)); then exit "$plot_status"; fi
