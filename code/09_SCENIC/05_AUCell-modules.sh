#!/bin/bash

#SBATCH --mem=100G
#SBATCH --job-name=aucell_modules_fixed-thold_custom-cluster
#SBATCH -o code/09_SCENIC/logs/%x_%j_%a.log
#SBATCH --ntasks=10
#SBATCH --array=1,4-12

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

nGenes=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print$2}' code/09_SCENIC/array_thresholds.txt)
cluster=$(awk -v Index=$SLURM_ARRAY_TASK_ID '$1==Index {print$3}' code/09_SCENIC/array_thresholds.txt)

echo $cluster
echo "cluster-specific 5% detected genes cutoff: ${nGenes}"
echo ""
echo ""

python3 /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/09_SCENIC/05_AUCell-modules.py $nGenes $cluster

echo "**** Job ends ****"
date
