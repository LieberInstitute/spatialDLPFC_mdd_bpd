#!/bin/bash
#SBATCH --mem=200GB
#SBATCH --ntasks=12
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jthom338@jh.edu
#SBATCH --job-name=precast_n120_n1051_k9
#SBATCH --output=/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/05_clustering/logs/%x_%j.log

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

Rscript /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/05_clustering/02_PRECAST-n120.r

echo "**** Job ends ****"
date
