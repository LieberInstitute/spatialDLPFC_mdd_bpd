#!/usr/bin/env bash
#SBATCH --job-name=mbv-aucell-python
#SBATCH --mem=5G
#SBATCH -o test_results/apptainer/%x_%j.log

PYTHON_SIF="reproducibility_outputs/apptainer/20261111T141200Z/python-jhpce/mbv-python-jhpce.sif"

module load apptainer/1.5.3

apptainer exec \
  --bind "$PWD:/work" \
  --bind "/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/processed-data/09_DEG_GRN:/source-data:ro" \
  --pwd /work \
  "$PYTHON_SIF" \
  python3 -u reproducibility_codeCheck/apptainer-aucell.py \
  2>&1

status=$?
printf 'exit status: %s\n' "$status"
exit "$status"
