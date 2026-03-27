#!/bin/bash

#SBATCH --mem=20G
#SBATCH --job-name=cisTarget_119_DE-input-13162_regulons_filtered_top20percent
#SBATCH -o code/10_SCENIC/logs/%x_%j.log

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node(s): ${SLURM_NODELIST}"
echo "Node memory requested: ${SLURM_MEM_PER_NODE}"
echo "n Tasks: ${SLURM_NTASKS}"

module load conda
module list
source activate pyscenic_bioconda

python3 /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/10_SCENIC/06_cisTarget_enrichment.py

echo "**** Job ends ****"
date
