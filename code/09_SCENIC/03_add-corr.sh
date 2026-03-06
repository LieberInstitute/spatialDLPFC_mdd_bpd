#!/bin/bash

#SBATCH --mem=300G
#SBATCH --job-name=add-cor_full-logcounts_DE-input_no-lowUMI
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


adj_path="processed-data/09_SCENIC/spe-n119_13844-no-lowUMI_adj.csv"
echo $adj_path
if [ ! -f $adj_path ]; then
    echo "Adjacency file not found!"
    echo ""
    echo "**** Forced stop ****"
    date
    scancel ${SLURM_JOB_ID}
fi

loom_path="processed-data/09_SCENIC/spe-n119_13844-genes_no-lowUMI_logcounts.loom"
echo $loom_path
if [ ! -f $loom_path ]; then
    echo "Loom file not found!"
    echo ""
    echo "**** Forced stop ****"
    date
    scancel ${SLURM_JOB_ID}
fi

out_path="processed-data/09_SCENIC/spe-n119_13844-no-lowUMI_adj_with-logcounts-corr_duplicate.csv"
echo $out_path
if [ -f $out_path ]; then
    echo "File already exists at designated output path!"
    echo "Overwriting existing file..."
#    echo "Check output path and comment out these lines to force overwrite."
#    echo ""
#    echo "**** Forced stop ****"
#    date
#    scancel ${SLURM_JOB_ID}
fi

pyscenic add_cor \
	--expression_mtx_fname $loom_path \
	--output $out_path \
	--mask_dropouts \
	$adj_path $loom_path

# --mask_dropouts

echo "**** Job ends ****"
date
