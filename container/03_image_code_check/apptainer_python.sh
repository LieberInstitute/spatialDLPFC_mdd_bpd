#!/usr/bin/env bash
#SBATCH --job-name=mbv-aucell-python
#SBATCH --mem=5G
#SBATCH -o container/outputs/image_code_check/apptainer/%x_%j.log

PYTHON_SIF="container/outputs/image_build_test/apptainer/20261111T141200Z/python-jhpce/mbv-python-jhpce.sif"

module load apptainer/1.5.3

apptainer exec \
  --bind "$PWD:/work" \
  --bind "/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/processed-data/09_DEG_GRN:/source-data:ro" \
  --pwd /work \
  "$PYTHON_SIF" \
  python3 -u container/03_image_code_check/apptainer-aucell.py \
  2>&1

status=$?
printf 'exit status: %s\n' "$status"
exit "$status"
