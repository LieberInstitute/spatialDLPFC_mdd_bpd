#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd "${script_dir}/../../.." && pwd)

prefix="${repo_root}/processed-data/00_genotypes/plink2/merged_maf05"

if [[ ! -f "${prefix}.pgen" || ! -f "${prefix}.psam" || ! -f "${prefix}.pvar" ]]; then
  echo "ERROR: missing genotype inputs for prefix ${prefix}" >&2
  exit 1
fi

echo "Using PLINK prefix: ${prefix}"

# LD prune first, then compute top ancestry PCs on the pruned variant set.
plink2 --pfile "${prefix}" --indep-pairwise 50 5 0.2 --out "${prefix}_ldpruned"
plink2 --pfile "${prefix}" \
  --extract "${prefix}_ldpruned.prune.in" \
  --pca 5 \
  --out "${prefix}_pca"

echo "Wrote: ${prefix}_pca.eigenvec"
