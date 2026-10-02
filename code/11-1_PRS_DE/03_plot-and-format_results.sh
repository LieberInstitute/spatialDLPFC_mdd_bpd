#!/bin/bash
#SBATCH --mem=3G
#SBATCH --job-name=plot_format_layer-adjusted_no-lowUMI
#SBATCH --output=/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/tmp_PRS_DE/logs/%x_%j_%a.log
#SBATCH --array=1-3

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node(s): ${SLURM_NODELIST}"
echo "Node memory requested: ${SLURM_MEM_PER_NODE}"
echo "n Tasks: ${SLURM_NTASKS}"

input=$(head -n $SLURM_ARRAY_TASK_ID /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/tmp_PRS_DE/PRS_list.txt | tail -n 1)
echo $input

module load conda_R/4.4.x
module list

Rscript /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/tmp_PRS_DE/03_plot-and-format_results.r $input

echo "**** Job ends ****"
date
