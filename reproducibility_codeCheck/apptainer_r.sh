#!/usr/bin/env bash
#SBATCH --job-name=mbv-spe-r
#SBATCH --mem=5G
#SBATCH -o test_results/apptainer/%x_%j.log

R_SIF="reproducibility_outputs/apptainer/20261111T211530Z/r-jhpce/mbv-r-jhpce.sif"

module load apptainer/1.5.3

apptainer exec \
  --bind "$PWD:/work" \
  --pwd /work \
  "$R_SIF" \
  Rscript --vanilla reproducibility_codeCheck/apptainer-spe.R test_results/apptainer/ \
  2>&1

status=$?
printf 'exit status: %s\n' "$status"
exit "$status"
