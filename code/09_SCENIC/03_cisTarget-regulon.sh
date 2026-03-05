#!/bin/bash

#SBATCH --mem=500G
#SBATCH --job-name=ctx_full_no-mask
#SBATCH -o code/09_SCENIC/logs/%x_%j.log
#SBATCH --ntasks=10

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

adj_path="processed-data/09_SCENIC/spe-n119_21077_adj.csv"
echo $adj_path
if [ ! -f $adj_path ]; then
    echo "File not found!"
fi

loom_path="processed-data/09_SCENIC/spe-n119_21077-genes.loom"

# motif-based anno
annot_path='/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/raw-data/SCENIC_aux/genome_annotation/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/raw-data/SCENIC_aux/genome_annotation/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
motif_path="/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/raw-data/SCENIC_aux/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

# track based anno
#annot_path='/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/raw-data/SCENIC_aux/genome_annotation/encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__10kb_up_and_down_tss.max.genes_vs_tracks.rankings.feather /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/raw-data/SCENIC_aux/genome_annotation/encode_20190621__ChIP_seq_transcription_factor.hg38__refseq-r80__500bp_up_and_100bp_down_tss.max.genes_vs_tracks.rankings.feather'
#motif_path="/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/raw-data/SCENIC_aux/track2tf/encode_project_20190621__ChIP-seq_transcription_factor.homo_sapiens.hg38.bigwig_signal_pvalue.track_to_tf_in_motif_to_tf_format.tsv"

out_path="processed-data/09_SCENIC/spe-n119_21077_reg_no-mask.csv"

pyscenic ctx $adj_path $annot_path \
	--annotations_fname $motif_path \
	--expression_mtx_fname $loom_path \
	--output $out_path \
	--num_workers 10
#--thresholds 0.5
#--mask_dropouts --all_modules

echo "**** Job ends ****"
date
