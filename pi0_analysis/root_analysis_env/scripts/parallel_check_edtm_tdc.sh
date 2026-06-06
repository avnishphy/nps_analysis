#!/usr/bin/env bash
# parallel_check_edtm_tdc.sh
# Usage: ./parallel_check_edtm_tdc.sh "<run_list>" <root_dir> <output_dir> [max_jobs]
# Example: ./parallel_check_edtm_tdc.sh "4387-4390,4400 4500-4502" /path/to/rootfiles ./pdfs 4

set -e

run_list="$1"
root_dir="$2"
out_dir="$3"
max_jobs="${4:-4}"

if [[ -z "$run_list" || -z "$root_dir" || -z "$out_dir" ]]; then
  echo "Usage: $0 \"<run_list>\" <root_dir> <output_dir> [max_jobs]"
  exit 1
fi

mkdir -p "$out_dir"

# Expand run list (ranges, commas, spaces)
expand_runs() {
  local input="$1"
  local runs=()
  input="${input//,/ }"
  for token in $input; do
    if [[ "$token" =~ ^([0-9]+)-([0-9]+)$ ]]; then
      for ((i=${BASH_REMATCH[1]}; i<=${BASH_REMATCH[2]}; ++i)); do
        runs+=("$i")
      done
    else
      runs+=("$token")
    fi
  done
  # Remove duplicates and sort
  printf "%s\n" "${runs[@]}" | sort -n | uniq
}

runs=( $(expand_runs "$run_list") )
total=${#runs[@]}
completed=0

# Job pool
job_count() { jobs -pr | wc -l; }

for run in "${runs[@]}"; do
  out_pdf="$out_dir/output_${run}.pdf"
  ./check_edtm_tdc "$run" "$root_dir" "$out_pdf" &
  while (( $(job_count) >= max_jobs )); do
    wait -n
    ((completed++))
    echo -ne "\rProgress: $completed/$total runs finished"
  done
done

# Wait for remaining jobs
while (( $(job_count) > 0 )); do
  wait -n
  ((completed++))
  echo -ne "\rProgress: $completed/$total runs finished"
done
echo -e "\nAll runs complete."
