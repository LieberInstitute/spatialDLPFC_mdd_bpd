#!/bin/bash
#SBATCH --mem=80G
#SBATCH --job-name=mbv-refineVNS
#SBATCH -o logs/mbv-refineVNS250521o-%a.txt
#SBATCH --array=1-24

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "****"
echo "Sample id: $(cat /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/VistoSeg/code/refineVNS_list.txt | awk '{print $NF}' | awk "NR==${SLURM_ARRAY_TASK_ID}")"
echo "****"

## load MATLAB
module load matlab/R2023b

## Load toolbox for VistoSeg
toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/VistoSeg/code'

## Read inputs from refineVNS_list.txt file
fname=$(awk 'BEGIN {FS="\t"} {print $1}' refineVNS_list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
M=$(awk 'BEGIN {FS="\t"} {print $2}' refineVNS_list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")

## Run refineVNS function
matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), refineVNS('$fname',$M)"

echo "**** Job ends ****"
date



