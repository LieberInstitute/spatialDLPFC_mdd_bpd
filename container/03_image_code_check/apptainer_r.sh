#!/usr/bin/env bash
#SBATCH --job-name=mbv-spe-r
#SBATCH --mem=5G
#SBATCH -o container/outputs/image_code_check/apptainer/%x_%j.log

R_SIF="container/outputs/image_build_test/apptainer/20261111T211530Z/r-jhpce/mbv-r-jhpce.sif"

module load apptainer/1.5.3

apptainer exec \
  --bind "$PWD:/work" \
  --pwd /work \
  "$R_SIF" \
  Rscript --vanilla container/03_image_code_check/spe.R container/outputs/image_code_check/apptainer/ \
  2>&1

status=$?
printf 'exit status: %s\n' "$status"
exit "$status"
