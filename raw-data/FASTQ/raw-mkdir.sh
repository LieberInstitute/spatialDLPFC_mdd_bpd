#!/bin/bash
#SBATCH --array=1-16

SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" raw-mkdir.txt)

mkdir -p ./logs/
mkdir -p /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/raw-data/FASTQ/${SAMPLE}/

mv slurm-*.out ./logs
