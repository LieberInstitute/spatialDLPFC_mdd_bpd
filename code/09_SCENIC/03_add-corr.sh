#!/bin/bash

#SBATCH --mem=300G
#SBATCH --job-name=add-cor_full_DE-input_no-mask
#SBATCH -o code/09_SCENIC/logs/%x_%j.log

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
source activate pyscenic_bioconda

adj_path="processed-data/09_SCENIC/spe-n119_13844_adj.csv"
echo $adj_path
if [ ! -f $adj_path ]; then
    echo "File not found!"
fi

loom_path="processed-data/09_SCENIC/spe-n119_13844-genes.loom"

out_path="processed-data/09_SCENIC/spe-n119_13844_adj_with-corr-no-mask.csv"

pyscenic add_cor \
	--expression_mtx_fname $loom_path \
	--output $out_path \
	$adj_path $loom_path

# --mask_dropouts

echo "**** Job ends ****"
date
