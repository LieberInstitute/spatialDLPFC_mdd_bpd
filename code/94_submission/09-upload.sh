#!/bin/bash
#SBATCH -p transfer
#SBATCH --mem=5G
#SBATCH --job-name=09-upload
#SBATCH -c 1
#SBATCH -o ../../processed-data/94_submission/nda-sub/logs/09-upload_submit.txt
#SBATCH -e ../../processed-data/94_submission/nda-sub/logs/09-upload_submit.txt

set -e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${SLURMD_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load nda-tools/0.7.0
module list

repo_dir=$(git rev-parse --show-toplevel)
meta_dir=$repo_dir/processed-data/94_submission/nda-sub
data_dir=$repo_dir/processed-data/94_submission/nda-sub/to_upload

#   Try uploading all data except the imaging data, which had mysterious
#   validation issues
 vtcmd \
     -l $data_dir \
     -b \
     -c 5229 \
     -t "LIBD spatial DLPFC NTC/MDD/BPD RNA-seq and imaging submission" \
     -d "All imaging and RNA-seq FASTQ data for the spatialDLPFC NTC/MDD/BPD project" \
     -u rmillerlibd \
     $meta_dir/rna_seq.csv $meta_dir/visium_image.csv $meta_dir/genomics_subject.csv

#   Upload imaging data later, after fixing the validation issues
# vtcmd \
#     -l $data_dir \
#     -b \
#     -c 5229 \
#     -t "LIBD spatial DLPFC imaging submission" \
#     -d "All H&E and IF images for the spatialDLPFC project" \
#     -u nickeagles77 \
#     $meta_dir/visium_image.csv

#   Resubmit metadata files after QA issues and fixes. Had to be done
#   interactively because of yes/no prompts
#vtcmd \
#    -b \
#    -rs 65813 \
#    -u rmillerlibd \
#    $meta_dir/rna_seq.csv $meta_dir/snp_array.csv $meta_dir/genomics_subject.csv

#vtcmd \
#    -b \
#    -rs 65841 \
#    -u rmillerlibd \
#    $meta_dir/visium_image.csv

echo "**** Job ends ****"
date

## This script was made using slurmjobs version 1.2.1
## available from http://research.libd.org/slurmjobs/
