#!/bin/bash

#SBATCH --mem=500G
#SBATCH --job-name=aucell_full
#SBATCH -o code/09_SCENIC/logs/%x_%j.log
#SBATCH --ntasks=24

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
source activate pyscenic

reg_path="processed-data/09_SCENIC/spe-n119_21077_reg.csv"
echo $reg_path
if [ ! -f $reg_path ]; then
    echo "File not found!"
fi

loom_path="processed-data/09_SCENIC/spe-n119_21077-genes.loom"
out_path="processed-data/09_SCENIC/spe-n119_21077_AUCell-output.loom"

pyscenic aucell $loom_path $reg_path \
	--output $out_path \
	--seed 1234 --num_workers 10

#echo $out_path

echo "**** Job ends ****"
date
