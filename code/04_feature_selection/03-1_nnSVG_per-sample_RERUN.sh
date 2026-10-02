#!/bin/bash
#SBATCH --ntasks=12
#SBATCH --mem=30G
#SBATCH --time=0-05:00:00
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jthom338@jh.edu
#SBATCH --nodes=1
#SBATCH --job-name=standard-nnSVG_per-sample_min-100_RERUN
#SBATCH --output=/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/04_feature_selection/array_logs/%x_%j_%a.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node(s): ${SLURM_NODELIST}"
echo "Node memory requested: ${SLURM_MEM_PER_NODE}"
echo "n Tasks: ${SLURM_NTASKS}"

input=$(head -n $SLURM_ARRAY_TASK_ID /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/processed-data/04_feature_selection/per-sample_spe_RERUN_list.txt | tail -n 1)
echo $input
module load conda_R/4.4.x
module list
Rscript /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/04_feature_selection/03_nnSVG_per-sample_HDF5.r $input

echo "**** Job ends ****"
date
