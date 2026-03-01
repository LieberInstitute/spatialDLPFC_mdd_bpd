#!/bin/bash

#SBATCH --mem=30G
#SBATCH --job-name=ctx_single-sample
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
source activate pyscenic

adj_path="processed-data/09_SCENIC/single-sample_adj.csv"
echo $adj_path
if [ ! -f $adj_path ]; then
    echo "File not found!"
fi

loom_path="processed-data/09_SCENIC/single-sample_V13B23-329-A1_17300-genes.loom"

annot_path='/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/raw-data/SCENIC_aux/genome_annotation/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/raw-data/SCENIC_aux/genome_annotation/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
motif_path="/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/raw-data/SCENIC_aux/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

out_path="processed-data/09_SCENIC/single-sample_reg.csv"

pyscenic ctx $adj_path $annot_path \
	--annotations_fname $motif_path \
	--expression_mtx_fname $loom_path \
	--output $out_path \
	--mask_dropouts --all_modules --num_workers 20

#echo $out_path

echo "**** Job ends ****"
date
