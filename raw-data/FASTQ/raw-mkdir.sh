#!/bin/bash
#SBATCH -o logs/mkdir-%a.txt
#SBATCH --array=1-96

SAMPLE=$(awk 'BEGIN {FS="\t"} {print $1}' raw-mkdir.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
LINKMAKE=$(awk 'BEGIN {FS="\t"} {print $2}' raw-mkdir.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
REGPATERN='*fastq.gz'

mkdir -p ./logs/
mkdir -p /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/raw-data/FASTQ/${SAMPLE}/

#ln -s /dcs04/lieber/lcolladotor/rawDataTDSC_LIBD001/raw-data/2024-08-23_Psomagen/${LINKMAKE}/22GHKFLT4/${REGPATERN} /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/raw-data/FASTQ/${SAMPLE}/

#mv slurm-*.out ./logs
