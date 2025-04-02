#!/bin/bash
#SBATCH --mem=10G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jthom338@jh.edu
#SBATCH --job-name=prep-precast_HDF5_H-M-markers
#SBATCH --output=/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/05_clustering/PRECAST/logs/%x_%j.log

set -e
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

# these lines necessary only if performing multiple k
#echo "set stack size to unlimited"
#ulimit -s unlimited
#ulimit -s

Rscript /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/05_clustering/PRECAST/01-alt_prep-PRECAST_H-M_HDF5.r

echo "**** Job ends ****"
date
