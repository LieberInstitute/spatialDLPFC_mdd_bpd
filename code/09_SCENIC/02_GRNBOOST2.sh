#!/bin/bash

#SBATCH --mem=100G
#SBATCH --job-name=grn_full_dask-5GB-sparse
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

loom_path="processed-data/09_SCENIC/spe-n119_21077-genes.loom"
echo $loom_path
if [ ! -f $loom_path ]; then
    echo "File not found!"
fi


f_tfs="raw-data/SCENIC_aux/tf_lists/allTFs_hg38.txt"
out_path="processed-data/09_SCENIC/spe-n119_21077_adj.csv"

pyscenic grn $loom_path $f_tfs -o $out_path --num_workers 20 --seed 1234 --sparse

#arboreto_with_multiprocessing.py \
#    $loom_path $f_tfs \
#    --method grnboost2 \
#    --output $out_path \
#    --num_workers 12 \
#    --seed 1234

#echo $out_path

echo "**** Job ends ****"
date
