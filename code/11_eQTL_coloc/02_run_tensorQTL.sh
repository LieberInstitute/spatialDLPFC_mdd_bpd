#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "${script_dir}/../.." && pwd)

in_dir="${repo_root}/processed-data/11_eQTL_coloc/tqtl_in"
out_dir="${repo_root}/processed-data/11_eQTL_coloc/tqtl_out"
manifest="${in_dir}/prep_manifest.csv"
mapping_py="${script_dir}/02a_tensorQTL_cis.py"
py_bin="${script_dir}/.venv/bin/python"

dry_run=0
force=0
start_from=""
only_regex=""
declare -a user_datasets=()

usage() {
  cat <<'EOF'
Usage:
  ./02_run_tensorQTL.sh [options] [dataset_id ...]

Runs tensorQTL (02a_tensorQTL_cis.py) for prepared datasets in:
  processed-data/11_eQTL_coloc/tqtl_in

Options:
  --dry-run            Show detected datasets/commands only; do not run
  --force              Re-run datasets even if <dataset>.gene.map_cis.tab.gz exists
  --start-from <id>    Start at this dataset (sorted order), skipping previous
  --only <regex>       Keep only dataset IDs matching regex
  -h, --help           Show this help

If dataset_id arguments are provided, only those are considered.
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dry-run)
      dry_run=1
      shift
      ;;
    --force)
      force=1
      shift
      ;;
    --start-from)
      [[ $# -ge 2 ]] || { echo "Missing value for --start-from" >&2; exit 1; }
      start_from="$2"
      shift 2
      ;;
    --only)
      [[ $# -ge 2 ]] || { echo "Missing value for --only" >&2; exit 1; }
      only_regex="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    --)
      shift
      while [[ $# -gt 0 ]]; do
        user_datasets+=("$1")
        shift
      done
      ;;
    -*)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 1
      ;;
    *)
      user_datasets+=("$1")
      shift
      ;;
  esac
done

[[ -d "${in_dir}" ]] || { echo "Missing input directory: ${in_dir}" >&2; exit 1; }
[[ -f "${mapping_py}" ]] || { echo "Missing runner script: ${mapping_py}" >&2; exit 1; }
mkdir -p "${out_dir}" "${out_dir}/logs"

if [[ ! -x "${py_bin}" ]]; then
  py_bin="python"
fi

if [[ ${#user_datasets[@]} -gt 0 ]]; then
  mapfile -t dataset_ids < <(printf "%s\n" "${user_datasets[@]}" | sed '/^$/d' | sort -u)
elif [[ -f "${manifest}" ]]; then
  mapfile -t dataset_ids < <(awk -F',' 'NR>1 && $1 != "" {print $1}' "${manifest}" | sort -u)
else
  mapfile -t dataset_ids < <(
    find "${in_dir}" -maxdepth 1 -type f -name '*.gene.expr.bed.gz' -printf '%f\n' |
      sed 's/\.gene\.expr\.bed\.gz$//' | sort -u
  )
fi

if [[ -n "${only_regex}" ]]; then
  mapfile -t dataset_ids < <(printf "%s\n" "${dataset_ids[@]}" | grep -E "${only_regex}" || true)
fi

if [[ ${#dataset_ids[@]} -eq 0 ]]; then
  echo "No dataset IDs found to run." >&2
  exit 1
fi

# Prioritize run order:
#  1) all-donor datasets (no _f/_m suffix)
#  2) female-only datasets (_f)
#  3) male-only datasets (_m)
mapfile -t dataset_ids < <(
  printf "%s\n" "${dataset_ids[@]}" |
    awk '
      function tier(ds) {
        if (ds ~ /_f$/) return 2
        if (ds ~ /_m$/) return 3
        return 1
      }
      { print tier($0) "\t" $0 }
    ' |
    sort -k1,1n -k2,2 |
    cut -f2
)

if [[ -n "${start_from}" ]]; then
  seen=0
  declare -a trimmed=()
  for ds in "${dataset_ids[@]}"; do
    if [[ ${seen} -eq 0 && "${ds}" == "${start_from}" ]]; then
      seen=1
    fi
    if [[ ${seen} -eq 1 ]]; then
      trimmed+=("${ds}")
    fi
  done
  if [[ ${seen} -eq 0 ]]; then
    echo "--start-from dataset not found: ${start_from}" >&2
    exit 1
  fi
  dataset_ids=("${trimmed[@]}")
fi

for ds in "${dataset_ids[@]}"; do
  expr="${in_dir}/${ds}.gene.expr.bed.gz"
  cov="${in_dir}/${ds}.gene.covars.txt"
  [[ -f "${expr}" ]] || { echo "Missing expression BED for ${ds}: ${expr}" >&2; exit 1; }
  [[ -f "${cov}" ]] || { echo "Missing covariates for ${ds}: ${cov}" >&2; exit 1; }
done

echo "Python: ${py_bin}"
echo "Input dir: ${in_dir}"
echo "Output dir: ${out_dir}"
echo "Datasets (${#dataset_ids[@]}): ${dataset_ids[*]}"

if [[ ${dry_run} -eq 1 ]]; then
  echo "Dry run complete."
  exit 0
fi

run_stamp=$(date '+%Y%m%d_%H%M%S')
ran=0
skipped=0

for ds in "${dataset_ids[@]}"; do
  cis_out="${out_dir}/${ds}.gene.map_cis.tab.gz"
  logfile="${out_dir}/logs/${run_stamp}_${ds}.log"

  if [[ ${force} -eq 0 && -f "${cis_out}" ]]; then
    echo "Skipping ${ds} (exists: ${cis_out}); use --force to re-run."
    skipped=$((skipped + 1))
    continue
  fi

  echo "[$(date '+%F %T')] Running ${ds}"
  "${py_bin}" -u "${mapping_py}" \
    --input-dir "${in_dir}" \
    --output-dir "${out_dir}" \
    "${ds}" 2>&1 | sed -u 's/\r/\n/g' | tee "${logfile}"
  ran=$((ran + 1))
done

echo "Completed. Ran: ${ran}; skipped: ${skipped}; total considered: ${#dataset_ids[@]}."
