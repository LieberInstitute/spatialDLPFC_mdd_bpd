#!/usr/bin/env bash
set -euo pipefail

## run all non-MAPK3 target checkpoints with bounded chromosome-independent workers
base=${1:?usage: 11_run_targeted_susie_all.sh BASE [JOBS] [none|eb]}
jobs=${2:-8}
mismatch_mode=${3:-none}
if [[ "$mismatch_mode" != none && "$mismatch_mode" != eb ]]; then
    echo "mismatch mode must be none or eb" >&2
    exit 2
fi
script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
runner="$script_dir/10_run_targeted_susie_one.R"
manifest="$base/target_manifest.tsv"
result_name=results
log_name=logs
if [[ "$mismatch_mode" == eb ]]; then
    result_name=results_eb
    log_name=logs_eb
fi
log_dir="$base/$log_name"
status_dir="$log_dir/status"
queue="$log_dir/target_queue.tsv"
mkdir -p "$log_dir" "$status_dir" "$base/$result_name"

## largest-first ordering balances workers and leaves completed checkpoints untouched
awk -F '\t' '
NR == 1 {
    for (i=1; i<=NF; i++) {
        if ($i == "target_id") target_col=i
        if ($i == "mapk3_pilot_complete") mapk3_col=i
        if ($i == "n_nominal") n_col=i
    }
    next
}
$mapk3_col == "FALSE" {print $target_col "\t" $n_col}
' "$manifest" | sort -t $'\t' -k2,2nr > "$queue"

run_one() {
    local id=$1 base=$2 runner=$3 log_dir=$4 status_dir=$5 mode=$6 result_name status started ended
    result_name=results
    if [[ "$mode" == eb ]]; then result_name=results_eb; fi
    if [[ -s "$base/$result_name/${id}.rds" ]]; then
        printf '%s\t%s\t%s\n' "$id" 0 skipped > "$status_dir/${id}.tsv"
        return 0
    fi
    started=$(date '+%Y-%m-%d %H:%M:%S %z')
    set +e
    env OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4 \
        VECLIB_MAXIMUM_THREADS=4 NUMEXPR_NUM_THREADS=4 \
        /usr/bin/time -v Rscript "$runner" "$base" "$id" "$mode" \
        > "$log_dir/${id}.out" 2> "$log_dir/${id}.err"
    status=$?
    set -e
    ended=$(date '+%Y-%m-%d %H:%M:%S %z')
    printf '%s\t%s\t%s\t%s\n' "$id" "$status" "$started" "$ended" \
        > "$status_dir/${id}.tsv"
    return "$status"
}
export -f run_one

## profiling supports eight workers on this 32-CPU, 503-GiB host
cut -f1 "$queue" | xargs -r -P "$jobs" -I '{}' \
    bash -c 'run_one "$1" "$2" "$3" "$4" "$5" "$6"' _ \
    '{}' "$base" "$runner" "$log_dir" "$status_dir" "$mismatch_mode"
