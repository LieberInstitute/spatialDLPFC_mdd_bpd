#!/bin/bash
#SBATCH --mem=80G
#SBATCH --job-name=mbv-countNuc
#SBATCH -o logs/mbv-countNuc250522o-%a.txt
#SBATCH --array=1-24%4

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "****"
echo "Sample id: $(cat /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/VistoSeg/code/countNuclei_list.txt | awk '{print $NF}' | awk "NR==${SLURM_ARRAY_TASK_ID}")"
echo "****"

## load MATLAB
module load matlab/R2023b

## Load toolbox for VistoSeg
toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/VistoSeg/code'

## Read parameters
mask=$(awk 'BEGIN {FS="\t"} {print $1}' countNuclei_list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
jsonname=$(awk 'BEGIN {FS="\t"} {print $2}' countNuclei_list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
posname=$(awk 'BEGIN {FS="\t"} {print $3}' countNuclei_list.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")


matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), countNuclei('$mask','$jsonname', '$posname')"

echo "**** Job ends ****"
date

