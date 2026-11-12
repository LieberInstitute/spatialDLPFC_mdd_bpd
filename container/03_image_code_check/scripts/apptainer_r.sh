#!/usr/bin/env bash
#SBATCH --job-name=mbv-spe-r
#SBATCH --mem=5G
#SBATCH -o container/outputs/image_code_check/apptainer/%x_%j.log

UTC_DIR=${1:?Usage: sbatch container/03_image_code_check/scripts/apptainer_r.sh <UTC-dir>}
if [[ ! "$UTC_DIR" =~ ^[0-9]{8}T[0-9]{6}Z$ ]]; then
  echo "UTC-dir must look like 20261112T141200Z" >&2
  exit 2
fi

R_SIF="container/outputs/image_build_test/apptainer/$UTC_DIR/r-jhpce/mbv-r-jhpce.sif"

module load apptainer/1.5.3

apptainer exec \
  --bind "$PWD:/work" \
  --pwd /work \
  "$R_SIF" \
  Rscript --vanilla container/03_image_code_check/scripts/spe.R container/outputs/image_code_check/apptainer/ \
  2>&1

status=$?
printf 'exit status: %s\n' "$status"
exit "$status"
