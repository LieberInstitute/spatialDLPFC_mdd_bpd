#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "${script_dir}/../.." && pwd)

in_dir="${repo_root}/processed-data/11_eQTL_coloc/seurat/tqtl_in"
out_dir="${repo_root}/processed-data/11_eQTL_coloc/seurat/tqtl_out"
mapping_py="${script_dir}/02a_tensorQTL_cis.py"
py_bin="${script_dir}/.venv/bin/python"

dry_run=0
force=0
analysis="seurat"
splits_csv="all"
splits_explicit=0
input_dir_explicit=0
output_dir_explicit=0
start_from=""
only_regex=""
declare -a user_datasets=()

usage() {
  cat <<'EOF'
Usage:
  ./02_run_tensorQTL.sh [options] [dataset_id ...]

Runs tensorQTL (02a_tensorQTL_cis.py) for prepared datasets.

Options:
  --dry-run            Show detected datasets/commands only; do not run
  --force              Re-run datasets even if <dataset>.gene.map_cis.tab.gz exists
  --analysis <name>    Input/output series: seurat or custom_cluster (default: seurat)
  --splits <csv>       Dataset splits to run: all,male,female (default: all)
  --include-sex-splits Shorthand for --splits all,male,female
  --input-dir <path>   Override prepared tensorQTL input directory
  --output-dir <path>  Override tensorQTL output directory
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
    --analysis)
      [[ $# -ge 2 ]] || { echo "Missing value for --analysis" >&2; exit 1; }
      analysis="$2"
      shift 2
      ;;
    --splits)
      [[ $# -ge 2 ]] || { echo "Missing value for --splits" >&2; exit 1; }
      splits_csv="$2"
      splits_explicit=1
      shift 2
      ;;
    --include-sex-splits)
      splits_csv="all,male,female"
      splits_explicit=1
      shift
      ;;
    --input-dir)
      [[ $# -ge 2 ]] || { echo "Missing value for --input-dir" >&2; exit 1; }
      in_dir="$2"
      input_dir_explicit=1
      shift 2
      ;;
    --output-dir)
      [[ $# -ge 2 ]] || { echo "Missing value for --output-dir" >&2; exit 1; }
      out_dir="$2"
      output_dir_explicit=1
      shift 2
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

case "${analysis}" in
  seurat)
    if [[ ${input_dir_explicit} -eq 0 ]]; then
      in_dir="${repo_root}/processed-data/11_eQTL_coloc/seurat/tqtl_in"
    fi
    if [[ ${output_dir_explicit} -eq 0 ]]; then
      out_dir="${repo_root}/processed-data/11_eQTL_coloc/seurat/tqtl_out"
    fi
    ;;
  custom_cluster)
    if [[ ${input_dir_explicit} -eq 0 ]]; then
      in_dir="${repo_root}/processed-data/11_eQTL_coloc/custom_cluster/tqtl_in"
    fi
    if [[ ${output_dir_explicit} -eq 0 ]]; then
      out_dir="${repo_root}/processed-data/11_eQTL_coloc/custom_cluster/tqtl_out"
    fi
    ;;
  *)
    echo "Unsupported --analysis: ${analysis}; expected seurat or custom_cluster" >&2
    exit 1
    ;;
esac

IFS=',' read -r -a selected_splits <<< "${splits_csv}"
declare -A selected_split_map=()
declare -a selected_splits_clean=()
for split in "${selected_splits[@]}"; do
  split=$(echo "${split}" | xargs)
  case "${split}" in
    all|male|female)
      if [[ -z "${selected_split_map[${split}]+x}" ]]; then
        selected_split_map["${split}"]=1
        selected_splits_clean+=("${split}")
      fi
      ;;
    "")
      ;;
    *)
      echo "Unsupported split in --splits: ${split}" >&2
      exit 1
      ;;
  esac
done
if [[ ${#selected_split_map[@]} -eq 0 ]]; then
  echo "--splits must include at least one of: all,male,female" >&2
  exit 1
fi
splits_csv=$(IFS=','; echo "${selected_splits_clean[*]}")

manifest="${in_dir}/prep_manifest.csv"

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

if [[ ${#user_datasets[@]} -eq 0 || ${splits_explicit} -eq 1 ]]; then
  if [[ -f "${manifest}" ]]; then
    mapfile -t split_dataset_ids < <(
      awk -F',' -v splits="${splits_csv}" '
        BEGIN {
          n = split(splits, arr, ",")
          for (i = 1; i <= n; i++) keep[arr[i]] = 1
        }
        NR == 1 {
          for (i = 1; i <= NF; i++) {
            if ($i == "dataset_id") ds_col = i
            if ($i == "split") split_col = i
          }
          next
        }
        ds_col > 0 && split_col > 0 && keep[$split_col] && $ds_col != "" {
          print $ds_col
        }
      ' "${manifest}" | sort -u
    )
    mapfile -t dataset_ids < <(
      awk 'NR == FNR {keep[$0] = 1; next} keep[$0]' \
        <(printf "%s\n" "${split_dataset_ids[@]}") \
        <(printf "%s\n" "${dataset_ids[@]}")
    )
  else
    mapfile -t dataset_ids < <(
      printf "%s\n" "${dataset_ids[@]}" |
        awk -v splits="${splits_csv}" '
          BEGIN {
            n = split(splits, arr, ",")
            for (i = 1; i <= n; i++) keep[arr[i]] = 1
          }
          function ds_split(ds) {
            if (ds ~ /_f$/) return "female"
            if (ds ~ /_m$/) return "male"
            return "all"
          }
          keep[ds_split($0)]
        '
    )
  fi
fi

if [[ -n "${only_regex}" ]]; then
  mapfile -t dataset_ids < <(printf "%s\n" "${dataset_ids[@]}" | grep -E "${only_regex}" || true)
fi

if [[ ${#dataset_ids[@]} -eq 0 ]]; then
  echo "No dataset IDs found to run." >&2
  exit 1
fi

## prioritize run order:
##  1) all-donor datasets (no _f/_m suffix)
##  2) female-only datasets (_f)
##  3) male-only datasets (_m)
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
echo "Analysis: ${analysis}"
echo "Splits: ${splits_csv}"
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
