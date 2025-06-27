#!/bin/bash
#SBATCH --mem=80G
#SBATCH --job-name=mbv-VNS
#SBATCH -o logs/mbv-VNS250516o-%a.txt
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
echo "Sample id: $(cat /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/VistoSeg/code/VNS_list.txt | awk '{print $NF}' | awk "NR==${SLURM_ARRAY_TASK_ID}")"
echo "****"

module load matlab/R2023b


toolbox='/dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/VistoSeg/code'
fname=$(cat /dcs04/lieber/marmaypag/spatialDLPFC_mdd_bpd_LIBD4100/spatialDLPFC_mdd_bpd/code/VistoSeg/code/VNS_list.txt | awk '{print $NF}' | awk "NR==${SLURM_ARRAY_TASK_ID}")

matlab -nodesktop -nosplash -nojvm -r "addpath(genpath('$toolbox')), VNS('$fname',5)"

echo "**** Job ends ****"
date



