#!/bin/bash
#SBATCH --mem=50G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jthom338@jh.edu
#SBATCH --job-name=test_bindev_HDF5
#SBATCH --output=/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/04_feature_selection/logs/%x_%j.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node(s): ${SLURM_NODELIST}"
echo "Node memory requested: ${SLURM_MEM_PER_NODE}"
echo "n Tasks: ${SLURM_NTASKS}"

module load conda_R/devel
module list
export TMPDIR=$MYSCRATCH
echo "new temp dir: $TMPDIR"
Rscript /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/04_feature_selection/test_bindev_HDF5.r

echo "check temp dir"
ls $MYSCRATCH
echo "**** Job ends ****"
date
