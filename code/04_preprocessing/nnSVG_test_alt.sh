#!/bin/bash
#SBATCH --ntasks=12
#SBATCH --mem-per-cpu=80G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=jthom338@jh.edu
#SBATCH --job-name=nnSVG_test_alt
#SBATCH --error=./%x_%j.err
#SBATCH --output=./%x_%j.out

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node(s): ${SLURM_NODELIST}"
echo "Node memory requested: ${SLURM_MEM_PER_NODE}"
echo "n Tasks: ${SLURM_ARRAY_TASK_COUNT}"

date
module load conda_R/devel
module list
Rscript /users/jthompso/nnSVG_test_alt.r

echo "Finished!"
date
