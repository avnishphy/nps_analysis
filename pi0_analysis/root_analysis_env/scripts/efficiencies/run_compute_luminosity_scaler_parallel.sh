#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  run_compute_luminosity_parallel.sh --runs <run tokens...> [options]

Required:
  --runs <tokens...>              Run tokens, e.g. 4398 4399 4400-4410
                                  Also accepts comma-separated input, e.g.
                                  4398,4399,4400-4410

Options:
  --binary <path>                 Path to compute executable
                                  (default: ./compute_luminosity_scaler)
  --root-dir <path>               ROOT files directory passed to compute
  --db <path>                     DB file passed to compute
  --out-csv <path>                Final merged CSV output
                                  (default: livetime_results_parallel.csv)
  --jobs <N>                      Parallel workers (default: nproc)
  --chunk-size <N>                Runs per job (default: auto from runs/jobs)
  --default-ps <value>            Forwarded to compute
  --current-window <uA>           Forwarded to compute
  --current-correction <uA>       Forwarded to compute
  --keep-temp                     Keep temporary chunk CSVs
  -h, --help                      Show this help

Examples:
  ./run_compute_luminosity_parallel.sh \
    --runs 4398 4399 4400-4422 \
    --root-dir /lustre24/expphy/volatile/hallc/nps/singhav/ROOTfiles/root_analysis_env_skim/x60_4b \
    --out-csv livetime_results_recomputed_parallel.csv

  ./run_compute_luminosity_parallel.sh \
    --runs 4398,4399,4400-4422 \
    --jobs 8
EOF
}

expand_run_token() {
  local token="$1"
  if [[ "$token" =~ ^[0-9]+-[0-9]+$ ]]; then
    local lo="${token%-*}"
    local hi="${token#*-}"
    if (( hi < lo )); then
      local tmp="$lo"
      lo="$hi"
      hi="$tmp"
    fi
    local r
    for ((r = lo; r <= hi; ++r)); do
      printf '%s\n' "$r"
    done
  elif [[ "$token" =~ ^[0-9]+$ ]]; then
    printf '%s\n' "$token"
  else
    printf 'Invalid run token: %s\n' "$token" >&2
    return 1
  fi
}

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

binary="$script_dir/compute_luminosity_scaler"
out_csv="$script_dir/livetime_results_parallel.csv"
root_dir=""
db_path=""
default_ps=""
current_window=""
current_correction=""
keep_temp=0

if command -v nproc >/dev/null 2>&1; then
  jobs="$(nproc)"
else
  jobs=1
fi

chunk_size=""
run_tokens=()

while (( $# > 0 )); do
  case "$1" in
    --runs)
      shift
      while (( $# > 0 )) && [[ ! "$1" =~ ^-- ]]; do
        run_tokens+=("$1")
        shift
      done
      continue
      ;;
    --binary)
      binary="$2"
      shift 2
      ;;
    --root-dir)
      root_dir="$2"
      shift 2
      ;;
    --db)
      db_path="$2"
      shift 2
      ;;
    --out-csv)
      out_csv="$2"
      shift 2
      ;;
    --jobs)
      jobs="$2"
      shift 2
      ;;
    --chunk-size)
      chunk_size="$2"
      shift 2
      ;;
    --default-ps)
      default_ps="$2"
      shift 2
      ;;
    --current-window)
      current_window="$2"
      shift 2
      ;;
    --current-correction)
      current_correction="$2"
      shift 2
      ;;
    --keep-temp)
      keep_temp=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      printf 'Unknown argument: %s\n' "$1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

if (( ${#run_tokens[@]} == 0 )); then
  printf 'Error: --runs is required\n' >&2
  usage >&2
  exit 1
fi

if [[ ! -x "$binary" ]]; then
  printf 'Error: executable not found or not executable: %s\n' "$binary" >&2
  exit 1
fi

if [[ ! "$jobs" =~ ^[0-9]+$ ]] || (( jobs < 1 )); then
  printf 'Error: --jobs must be a positive integer\n' >&2
  exit 1
fi

expanded_runs=()
for tok in "${run_tokens[@]}"; do
  IFS=',' read -r -a pieces <<< "$tok"
  for part in "${pieces[@]}"; do
    part="${part#"${part%%[![:space:]]*}"}"
    part="${part%"${part##*[![:space:]]}"}"
    [[ -z "$part" ]] && continue
    while IFS= read -r r; do
      expanded_runs+=("$r")
    done < <(expand_run_token "$part")
  done
done

if (( ${#expanded_runs[@]} == 0 )); then
  printf 'Error: no valid runs after parsing --runs\n' >&2
  exit 1
fi

mapfile -t unique_runs < <(printf '%s\n' "${expanded_runs[@]}" | sort -n | uniq)
num_runs="${#unique_runs[@]}"

if [[ -z "$chunk_size" ]]; then
  chunk_size=$(( (num_runs + jobs - 1) / jobs ))
  (( chunk_size < 1 )) && chunk_size=1
fi

if [[ ! "$chunk_size" =~ ^[0-9]+$ ]] || (( chunk_size < 1 )); then
  printf 'Error: --chunk-size must be a positive integer\n' >&2
  exit 1
fi

tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/luminosity_parallel.XXXXXX")"
cleanup() {
  if (( keep_temp == 0 )) && [[ -d "$tmpdir" ]]; then
    rm -rf "$tmpdir"
  fi
}
trap cleanup EXIT

printf 'Total unique runs: %d\n' "$num_runs"
printf 'Parallel jobs: %d\n' "$jobs"
printf 'Chunk size: %d\n' "$chunk_size"
printf 'Temporary dir: %s\n' "$tmpdir"

chunk_list="$tmpdir/chunks.list"
: > "$chunk_list"

chunk_idx=0
for ((i = 0; i < num_runs; i += chunk_size)); do
  chunk_file="$tmpdir/chunk_${chunk_idx}.runs"
  : > "$chunk_file"
  for ((j = i; j < i + chunk_size && j < num_runs; ++j)); do
    printf '%s\n' "${unique_runs[j]}" >> "$chunk_file"
  done
  printf '%s\n' "$chunk_file" >> "$chunk_list"
  ((chunk_idx += 1))
done

worker="$tmpdir/worker.sh"
cat > "$worker" <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

chunk_file="$1"
chunk_base="$(basename "$chunk_file" .runs)"
out_csv_chunk="$TMP_WORKDIR/${chunk_base}.csv"

mapfile -t runs < "$chunk_file"
if (( ${#runs[@]} == 0 )); then
  printf 'Skipping empty chunk: %s\n' "$chunk_file" >&2
  exit 0
fi

args=(
  "$BINARY"
  --runs "${runs[@]}"
  --out-csv "$out_csv_chunk"
)

if [[ -n "$ROOT_DIR" ]]; then
  args+=(--root-dir "$ROOT_DIR")
fi
if [[ -n "$DB_PATH" ]]; then
  args+=(--db "$DB_PATH")
fi
if [[ -n "$DEFAULT_PS" ]]; then
  args+=(--default-ps "$DEFAULT_PS")
fi
if [[ -n "$CURRENT_WINDOW" ]]; then
  args+=(--current-window "$CURRENT_WINDOW")
fi
if [[ -n "$CURRENT_CORRECTION" ]]; then
  args+=(--current-correction "$CURRENT_CORRECTION")
fi

"${args[@]}"
EOF
chmod +x "$worker"

export TMP_WORKDIR="$tmpdir"
export BINARY="$binary"
export ROOT_DIR="$root_dir"
export DB_PATH="$db_path"
export DEFAULT_PS="$default_ps"
export CURRENT_WINDOW="$current_window"
export CURRENT_CORRECTION="$current_correction"

xargs -P "$jobs" -n 1 "$worker" < "$chunk_list"

mapfile -t csv_chunks < <(ls -1 "$tmpdir"/chunk_*.csv 2>/dev/null | sort)
if (( ${#csv_chunks[@]} == 0 )); then
  printf 'Error: no chunk CSV files were generated\n' >&2
  exit 2
fi

header="$(head -n 1 "${csv_chunks[0]}")"
if [[ -z "$header" ]]; then
  printf 'Error: chunk CSV header is empty in %s\n' "${csv_chunks[0]}" >&2
  exit 2
fi

merge_tmp="$tmpdir/merged_body_unsorted.csv"
: > "$merge_tmp"
for c in "${csv_chunks[@]}"; do
  if [[ ! -s "$c" ]]; then
    continue
  fi
  tail -n +2 "$c" >> "$merge_tmp"
done

if [[ ! -s "$merge_tmp" ]]; then
  printf 'Error: merged CSV body is empty\n' >&2
  exit 2
fi

{
  printf '%s\n' "$header"
  sort -t, -k1,1n "$merge_tmp"
} > "$out_csv"

printf 'Wrote merged CSV: %s\n' "$out_csv"
printf 'Chunks merged: %d\n' "${#csv_chunks[@]}"
if (( keep_temp == 1 )); then
  printf 'Temporary files kept at: %s\n' "$tmpdir"
fi


# KinC_x60_4b runs comma separated:
# 4253,4254,4255,4256,4257,4258,4259,4260,4261,4262,4263,4264,4265,4266,4267,4268,4269,4270,4300,4301,4302,4303,4304,4305,4306,4308,4309,4311,4312,4313,4314,4315,4316,4317,4318,4319,4320,4321,4323,4324,4349,4350,4351,4352,4353,4354,4355,4356,4357,4358,4359,4360,4361,4362,4363,4364,4365,4366,4367,4368,4369,4370,4371,4372,4397,4398,4399,4400,4401,4402,4403,4404,4405,4406,4407,4408,4409,4410,4411,4414,4415,4417,4418,4419,4420,4422,4483,4484,4485,4487,4488,4489,4490,4491,4492,4493,4494,4495,4496,4498,4499,4500,4501,4502,4503,4504,4505,4506,4507,4508,4509,4551,4552,4555,4556,4557,4558,4559,4560,4561,4562,4563,4564,4565,4566,4568
# KinC_x60_3b runs comma separated:
# 4076,4077,4078,4079,4080,4081,4082,4083,4084,4085,4086,4087,4088,4089,4090,4091,4092,4093,4094,4095,4096,4097,4098,4099,4100,4101,4102,4103,4104,4105,4127,4128,4129,4130,4131,4132,4133,4134,4135,4136,4137,4138,4139,4140,4141,4142,4143,4144,4145,4163,4164,4165,4166,4167,4168,4169,4170,4171,4172,4173,4174,4175,4176,4177,4178,4179,4180,4181,4205,4206,4207,4208,4209,4210,4211,4212,4213,4214,4215