#!/usr/bin/env bash
#SBATCH --job-name=mbv-aucell-python
#SBATCH --mem=5G
#SBATCH -o container/outputs/image_code_check/apptainer/%x_%j.log

UTC_DIR=${1:?Usage: sbatch container/03_image_code_check/scripts/apptainer_python.sh <UTC-dir>}
if [[ ! "$UTC_DIR" =~ ^[0-9]{8}T[0-9]{6}Z$ ]]; then
  echo "UTC-dir must look like 20261112T141200Z" >&2
  exit 2
fi

PYTHON_SIF="container/outputs/image_build_test/apptainer/$UTC_DIR/python-jhpce/mbv-python-jhpce.sif"

module load apptainer/1.5.3

apptainer exec \
  --bind "$PWD:/work" \
  --bind "/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/processed-data/09_DEG_GRN:/source-data:ro" \
  --pwd /work \
  "$PYTHON_SIF" \
  python3 -u container/03_image_code_check/scripts/apptainer-aucell.py \
  2>&1

status=$?
printf 'exit status: %s\n' "$status"
exit "$status"
