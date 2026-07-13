#!/usr/bin/env bash

# Usage examples:
#   ./jcache.sh 6874..6880
#   ./jcache.sh 6874-6880 6901 7000..7002
#   ./jcache.sh 6874..6880,6901..6903

# Automatically detect all segment files
# for run in {4224..4568}; do
#     for file in /mss/hallc/c-nps/analysis/pass2/replays/updated/nps_hms_coin_${run}_*_1_-1.root; do
#         if [ -e "$file" ]; then
#             jcacheCMD="jcache get ${file} -x"
#             echo "Executing command: ${jcacheCMD}"
#             eval ${jcacheCMD}
#         else
#             echo "No files found for run ${run}"
#         fi
#     done
# done

expand_runs() {
    local token="$1"
    local part start end run

    # Support comma-separated chunks in each CLI argument.
    IFS=',' read -r -a chunks <<< "$token"
    for part in "${chunks[@]}"; do
        if [[ "$part" =~ ^([0-9]+)\.\.([0-9]+)$ ]] || [[ "$part" =~ ^([0-9]+)-([0-9]+)$ ]]; then
            start="${BASH_REMATCH[1]}"
            end="${BASH_REMATCH[2]}"

            if (( start <= end )); then
                for (( run=start; run<=end; run++ )); do
                    echo "$run"
                done
            else
                for (( run=start; run>=end; run-- )); do
                    echo "$run"
                done
            fi
        elif [[ "$part" =~ ^[0-9]+$ ]]; then
            echo "$part"
        else
            echo "Invalid run token: $part" >&2
            return 1
        fi
    done
}

if (( $# == 0 )); then
    set -- "6874..6880"
fi

for arg in "$@"; do
    while IFS= read -r run; do
        for file in /mss/hallc/c-nps/analysis/pass2/replays/updated/nps_hms_coin_${run}_*_1_-1.root; do
            if [ -e "$file" ]; then
                echo "Executing command: jcache get ${file} -x"
                jcache get "${file}" -x
            else
                echo "No files found for run ${run}"
            fi
        done
    done < <(expand_runs "$arg") || exit 1
done


# Runs for dc calib ranges left earlier:
# 6874-6880,2020-2183,3728-3782,6931-6939,3150-3173,3196-3219,3243-3269,3295-3318,3353-3354,3119-3149,3174-3195,3220-3242,3270-3294,3319-3352,3783-3806,3834-3857,3807-3832,3858-3882,3900-3943,3944-3988

# Runs for dc calib ranges KinC_x60_3b (Josh) and Kinc_x60_4b(Avnish)
# /work/hallc/nps/singhav/nps_analysis/jcache.sh 4253..4270 4300..4324 4349..4372 4397..4411 4413..4422 4483..4509 4551..4568 4076..4105 4127..4145 4163..4181 4203..4215

# Runs for lumonosity_2026 
# 1514-1530 6845-6854 7003-7007

