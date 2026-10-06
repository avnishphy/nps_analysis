#!/usr/bin/env bash
# Diagnostic objective validation. No real-data production or extraction.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO"
OUT="${NPS_PI0_OBJECTIVE_OUTPUT:-$REPO/output/pi0_objective_reproduction}"
JOBS="${NPS_PI0_OBJECTIVE_JOBS:-8}"
LIB="$OUT/build/libprimitive_candidate_v4.so"
export OUT LIB
case "${1:-help}" in
  build)
    command -v root-config >/dev/null || { echo 'Load /group/nps/singhav/setup.csh first.' >&2; exit 1; }
    [[ ! -e "$LIB" ]] || { echo "Refusing existing $LIB" >&2; exit 1; }
    mkdir -p "$OUT/build"
    g++ -O3 -std=c++17 -shared -fPIC scripts/pi0_primitive_objective.cpp \
      $(root-config --cflags --libs) -lMinuit2 -o "$LIB.tmp"
    mv "$LIB.tmp" "$LIB" ;;
  screen)
    python3 scripts/diagnose_pi0_fit_bias.py --output "$OUT/old_bias_regression"
    python3 scripts/study_pi0_objectives.py --library "$LIB" --output "$OUT/fixed_screen" \
      --experiments 2000 --fixed-shape --modes 0,1,2 --cases 0,1,2,3,4
    python3 scripts/study_pi0_objectives.py --library "$LIB" --output "$OUT/onoff_all_screen_v4" \
      --experiments 2000 --modes 0,1,2 --cases 0,1,2,3,4,5
    python3 scripts/study_pi0_objectives.py --library "$LIB" --output "$OUT/six_screen_v4" \
      --experiments 2000 --modes 0 --cases 0,1,2,3,4 --six-categories
    python3 scripts/study_pi0_objectives.py --library "$LIB" --output "$OUT/six_low_small_screen_v4" \
      --experiments 4000 --modes 0 --cases 5 --six-categories
    python3 scripts/study_pi0_objectives.py --library "$LIB" --output "$OUT/six_low_small_fixed_v4" \
      --experiments 4000 --modes 0 --cases 5 --six-categories --fixed-shape ;;
  coverage)
    mkdir -p "$OUT/logs"
    for case in 0 1 2 4 5; do for start in 0 100 200 300 400 500; do
      printf '%s %s\n' "$case" "$start"
    done; done | xargs -n2 -P "$JOBS" bash -c '
      python3 scripts/study_pi0_objectives.py --library "$LIB" \
        --output "$OUT/coverage_v4_c${1}_s${2}" --experiments 100 --start "$2" \
        --bootstrap 150 --modes 0 --six-categories --cases "$1" \
        > "$OUT/logs/coverage_c${1}_s${2}.log" 2>&1' _
    python3 scripts/collect_pi0_objective_toys.py --inputs "$OUT"/coverage_v4_c*_s* --output "$OUT/coverage_v4" ;;
  boundary)
    mkdir -p "$OUT/logs"
    for case in 0 1 2 3 4; do for start in 0 200 400; do
      printf '%s %s\n' "$case" "$start"
    done; done | xargs -n2 -P "$JOBS" bash -c '
      python3 scripts/study_pi0_objectives.py --library "$LIB" \
        --output "$OUT/boundary_v4_c${1}_s${2}" --experiments 200 --start "$2" \
        --bootstrap 100 --modes 0 --six-categories --cases "$1" \
        --amplitudes 0,0.05,0.10,0.20,0.40 > "$OUT/logs/boundary_c${1}_s${2}.log" 2>&1' _
    python3 scripts/collect_pi0_objective_toys.py --inputs "$OUT"/boundary_v4_c*_s* --output "$OUT/boundary_v4" ;;
  release-check)
    python3 scripts/check_pi0_objective_release.py --screen "$OUT/six_screen_v4" "$OUT/six_low_small_screen_v4" \
      --coverage "$OUT/coverage_v4" --boundary "$OUT/boundary_v4" --library "$LIB" \
      --output "$OUT/release_gate" --require-ready ;;
  *) echo 'Usage: bash scripts/run_pi0_objective_studies.sh {build|screen|coverage|boundary|release-check}'
     echo 'Load the Hall C ROOT environment. NPS_PI0_OBJECTIVE_OUTPUT selects a new output directory.'
     echo 'NPS_PI0_OBJECTIVE_JOBS controls independent toy workers (default 8).'
     echo 'Real-data regeneration is withheld until estimator validation passes.' ;;
esac
