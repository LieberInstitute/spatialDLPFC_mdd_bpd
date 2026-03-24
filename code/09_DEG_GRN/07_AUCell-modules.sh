#!/bin/bash

#SBATCH --mem=100G
#SBATCH --job-name=aucell_DE-input-13162_modules-DEG-subset-refined_fixed-thold-05_seurat-label_swap-modules
#SBATCH -o code/09_DEG_GRN/logs/%x_%j_%a.log
#SBATCH --ntasks=10
#SBATCH --array=1-8

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

cluster=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print$3}' code/09_DEG_GRN/array_thresholds.txt)

echo $cluster
echo ""
echo ""

python3 /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/09_DEG_GRN/07_AUCell-modules.py $cluster

echo "**** Job ends ****"
date
