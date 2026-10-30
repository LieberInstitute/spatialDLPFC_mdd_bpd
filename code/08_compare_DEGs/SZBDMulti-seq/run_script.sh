#!/bin/bash
#SBATCH --mem=5G
#SBATCH --job-name=plot-format_fixed-raw_age-ncells
#SBATCH --output=/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/08_compare_DEGs/SZBDMulti-seq/logs/%x_%j.log

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

#Rscript /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/08_compare_DEGs/SZBDMulti-seq/01_post-transfer_sn-qc.R
#Rscript /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/08_compare_DEGs/SZBDMulti-seq/02_pseudobulk.R
#Rscript /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/08_compare_DEGs/SZBDMulti-seq/03_norm-QC-PCA.R
#Rscript /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/08_compare_DEGs/SZBDMulti-seq/04_covariate-selection.R
#Rscript /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/08_compare_DEGs/SZBDMulti-seq/05_DE.R
Rscript /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/08_compare_DEGs/SZBDMulti-seq/06_plot-format-DE-results.R

echo "**** Job ends ****"
date
