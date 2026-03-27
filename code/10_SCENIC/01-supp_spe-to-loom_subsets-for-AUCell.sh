#!/bin/bash
#SBATCH --mem=50G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jthom338@jh.edu
#SBATCH --job-name=spe-logcounts_109_to_loom_13162-gene-names_seurat-subsets
#SBATCH --output=/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/10_SCENIC/logs/%x_%j.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node(s): ${SLURM_NODELIST}"
echo "Node memory requested: ${SLURM_MEM_PER_NODE}"
echo "n Tasks: ${SLURM_NTASKS}"

module load conda_R/4.4.x
module list
Rscript /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/10_SCENIC/01-supp_spe-to-loom_subsets-for-AUCell.r

echo "**** Job ends ****"
date
