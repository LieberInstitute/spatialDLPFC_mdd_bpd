#!/bin/bash

#SBATCH --mem=5G
#SBATCH --job-name=plot_aucell_DE-input-13162_modules-DEG-subset-refined_fixed-thold-05_swap-modules
#SBATCH -o code/09_DEG_GRN/logs/%x_%j_%a.log

echo "**** Job starts ****"
date

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

python3 /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/09_DEG_GRN/08_plot_AUCell-modules.py

echo "**** Job ends ****"
date
