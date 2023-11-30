#!/bin/bash
#SBATCH --mem=80G
#SBATCH -o logs/slurm-o_splitSlide.txt 
#SBATCH -e logs/slurm-e_splitSlide.txt
#SBATCH --job-name=MBv-splitslide
#SBATCH --array=1-4

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
samplelist="splitSlide_list.txt"

## Read inputs from splitSlide_list.txt file
fname=$(awk 'BEGIN {FS="\t"} {print $1}' ${samplelist} | awk "NR==${SLURM_ARRAY_TASK_ID}")
A1=$(awk 'BEGIN {FS="\t"} {print $2}' ${samplelist} | awk "NR==${SLURM_ARRAY_TASK_ID}")
B1=$(awk 'BEGIN {FS="\t"} {print $3}' ${samplelist} | awk "NR==${SLURM_ARRAY_TASK_ID}")
C1=$(awk 'BEGIN {FS="\t"} {print $4}' ${samplelist} | awk "NR==${SLURM_ARRAY_TASK_ID}")
D1=$(awk 'BEGIN {FS="\t"} {print $5}' ${samplelist} | awk "NR==${SLURM_ARRAY_TASK_ID}")

## Run refineVNS function
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), splitSlide('$fname',$A1,$B1,$C1,$D1)"

echo "**** Job ends ****"
date



