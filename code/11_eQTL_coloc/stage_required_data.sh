#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<USAGE
Usage: $(basename "$0") [--check-only] [--include-custom-cluster]

Checks required tensorQTL inputs under local processed-data.
If --check-only is omitted, missing files are fetched from JHPCE via rsync over ssh alias 'jt'.
By default, only Seurat all-donor workflow inputs are required.
USAGE
}

check_only=0
include_custom_cluster=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --check-only)
      check_only=1
      shift
      ;;
    --include-custom-cluster)
      include_custom_cluster=1
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown argument: $1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "${script_dir}/../.." && pwd)
jhpce_loc_file="${repo_root}/jhpce-location.txt"

if [[ ! -f "${jhpce_loc_file}" ]]; then
  echo "Missing ${jhpce_loc_file}" >&2
  exit 1
fi

jhpce_root=$(tr -d ' \t\r\n' < "${jhpce_loc_file}")
if [[ -z "${jhpce_root}" ]]; then
  echo "Empty JHPCE root path in ${jhpce_loc_file}" >&2
  exit 1
fi
## normalize to no trailing slash
jhpce_root=${jhpce_root%/}

required_files=(
  "processed-data/06_pseudobulk/Seurat/spe_n119_pseudo-no-lowUMI_sample-seurat-pc30_norm-filt.Rdata"
  "processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv"
  "processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_seurat-pc30_dx-sex_degs-F-test-t-test.csv"
  "processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv"
  "processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_smoothed-k9-1663_dx-sex_degs-F-test-t-test.csv"
  "processed-data/ref/granges.qs2"
)

custom_cluster_files=(
  "processed-data/06_pseudobulk/custom_cluster/spe_n119_pseudo_sample-custom-cluster_norm-filt.Rdata"
  "processed-data/07_dx_DE/layer-adjusted-pc3-age-nspots_custom-cluster_dx-sex_degs-F-test-t-test.csv"
  "processed-data/07_dx_DE/layer-restricted-pc3-age-nspots_custom-cluster_dx-sex_degs-F-test-t-test.csv"
)

if [[ ${include_custom_cluster} -eq 1 ]]; then
  required_files+=("${custom_cluster_files[@]}")
fi

geno_required=(
  "processed-data/00_genotypes/plink2/merged_maf05.pgen"
  "processed-data/00_genotypes/plink2/merged_maf05.psam"
  "processed-data/00_genotypes/plink2/merged_maf05.pvar"
)

missing=0
fetched=0

echo "Repo root: ${repo_root}"
echo "Remote root: ${jhpce_root}"
echo "Include custom cluster: ${include_custom_cluster}"

for rel in "${required_files[@]}"; do
  local_path="${repo_root}/${rel}"
  if [[ -f "${local_path}" ]]; then
    echo "[OK] ${rel}"
    continue
  fi
  echo "[MISS] ${rel}"
  missing=$((missing + 1))
  if [[ ${check_only} -eq 0 ]]; then
    remote_path="${jhpce_root}/${rel}"
    mkdir -p "$(dirname "${local_path}")"
    rsync -av "jt:${remote_path}" "${local_path}"
    if [[ -f "${local_path}" ]]; then
      echo "[FETCHED] ${rel}"
      fetched=$((fetched + 1))
    else
      echo "[ERROR] failed to fetch ${rel}" >&2
      exit 1
    fi
  fi
done

for rel in "${geno_required[@]}"; do
  local_path="${repo_root}/${rel}"
  if [[ ! -f "${local_path}" ]]; then
    echo "[ERROR] missing required genotype file: ${rel}" >&2
    exit 1
  fi
  echo "[OK] ${rel}"
done

echo "Summary: missing=${missing}, fetched=${fetched}, check_only=${check_only}"

if [[ ${check_only} -eq 1 && ${missing} -gt 0 ]]; then
  exit 2
fi
