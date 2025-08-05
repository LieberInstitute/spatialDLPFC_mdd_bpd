#!/bin/bash
#SBATCH --mem=80G
#SBATCH -c 8
#SBATCH -p katun
#SBATCH --job-name=spaceranger_segment
#SBATCH -o logs/spaceranger_segment_%a.txt
#SBATCH -e logs/spaceranger_segment_%a.txt
#SBATCH --array=1-24

#   Segment the full-res image for each Visium capture area

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load SpaceRanger
module load spaceranger/4.0.1

## List current modules for reproducibility
module list

repo_dir=/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd
sample_csv=$repo_dir/raw-data/sample_info/VistoSeg_samples_n24.csv
SAMPLE=$(awk -F',' 'NR>1 {print $1}' "$sample_csv" | sed -n "${SLURM_ARRAY_TASK_ID}p")
IMG_PATH=/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/raw-data/images/${SAMPLE}.tif
OUT_DIR=$repo_dir/processed-data/spaceranger_segment/${SAMPLE}

mkdir -p ${OUT_DIR}

echo "Processing sample ${SAMPLE}"

spaceranger segment \
    --id=${SAMPLE} \
    --tissue-image=${IMG_PATH} \
    --output-dir=${OUT_DIR} \
    --localcores=8 \
    --localmem=80 \
    --disable-ui

echo "**** Job ends ****"
date
