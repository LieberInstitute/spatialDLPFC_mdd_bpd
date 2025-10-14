#!/bin/bash
#SBATCH --mem=3G
#SBATCH --job-name=plot_LA-LR-overlap_rev-gene-input
#SBATCH --output=/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/07_dx_DE/logs/%x_%j_%a.log
#SBATCH --array=1-6

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

input=$(head -n $SLURM_ARRAY_TASK_ID /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/07_dx_DE/05-supp_LA-LR-overlap_groups.txt | tail -n 1)
echo $input
Rscript /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/07_dx_DE/05_LA-LR-overlap_dx-sex-dir-genes.r $input

echo "**** Job ends ****"
date
