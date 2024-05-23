#!/bin/bash
#SBATCH --mem=30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jthom338@jh.edu
#SBATCH --job-name=spatial-plot_outlier-overview_TEST
#SBATCH --error=/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/03_QC/logs/%x_%j.err
#SBATCH --output=/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/03_QC/logs/%x_%j.out

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node(s): ${SLURM_NODELIST}"
echo "n Tasks: ${SLURM_NTASKS}"

date
module load conda_R/devel
module list
Rscript /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/03_QC/02_spatial-plot_outlier-overview_TEST.r

echo "Finished!"
date
