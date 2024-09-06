#!/bin/bash
#SBATCH --mem=80G
#SBATCH -o logs/shift-splitSlide.txt 
#SBATCH --job-name=MBv-splitslide
#SBATCH --array=1
#SBATCH -t 2-00:00:00

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"s
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## load MATLAB
module load matlab/R2023a

## Load toolbox for VistoSeg
toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/VistoSeg/code'
samplelist="01-shifted-splitSlide-list.txt"

#shifted slides - V13B23-311.tif(C1,D1), V13B23-332.tif(A1, B1, C1, D1), V13B23-380.tif(A1)

## Read inputs from 01-shifted-splitSlide-list.txt file
fname=$(awk 'BEGIN {FS="\t"} {print $1}' ${samplelist} | awk "NR==${SLURM_ARRAY_TASK_ID}")
A1=$(awk 'BEGIN {FS="\t"} {print $2}' ${samplelist} | awk "NR==${SLURM_ARRAY_TASK_ID}")
B1=$(awk 'BEGIN {FS="\t"} {print $3}' ${samplelist} | awk "NR==${SLURM_ARRAY_TASK_ID}")
C1=$(awk 'BEGIN {FS="\t"} {print $4}' ${samplelist} | awk "NR==${SLURM_ARRAY_TASK_ID}")
D1=$(awk 'BEGIN {FS="\t"} {print $5}' ${samplelist} | awk "NR==${SLURM_ARRAY_TASK_ID}")

## Run refineVNS function
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide('$fname',$A1,$B1,$C1,$D1)"

echo "**** Job ends ****"
date



