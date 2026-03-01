#!/bin/bash

#SBATCH --mem=30G
#SBATCH --job-name=grn_single-slide
#SBATCH -o code/09_SCENIC/logs/%x_%j.log
#SBATCH --ntasks=20

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

f_loom_path_scenic="processed-data/09_SCENIC/single-slide_V13B23-329_17300-genes.loom"
echo $f_loom_path_scenic
if [ ! -f $f_loom_path_scenic ]; then
    echo "File not found!"
fi


f_tfs="raw-data/SCENIC_aux/tf_lists/allTFs_hg38.txt"
out_path="processed-data/09_SCENIC/single-slide_adj.csv"
pyscenic grn $f_loom_path_scenic $f_tfs -o $out_path --num_workers 20 --seed 1234

#echo $out_path

echo "**** Job ends ****"
date
