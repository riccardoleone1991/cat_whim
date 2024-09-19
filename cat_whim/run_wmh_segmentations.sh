#!/bin/bash
####    Run WMH Segmentation: Slurm -- Version 1.0.0
# Last edit:  2024-08-14
# Authors:    Leone, Riccardo (RL)
# Notes:      - Script to run WMH segmentations
#             - Release notes:
#                 * Initial release for longitudinal data
# To do:      - 
# Comments:   
# Sources:    
####


# SLURM settings
#SBATCH --job-name="wmh"
#SBATCH --time=2:00:00
#SBATCH --output=/home/leoner/Projects/cat_whim/slurmlogs/outputs/%x-%A-%a.out
#SBATCH --error=/home/leoner/Projects/cat_whim/slurmlogs/errors/%x-%A-%a.err
#SBATCH --nodes=1
#SBATCH --array=1-730
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=HPC-CPUs
#SBATCH --ntasks=1

# Define variables
HOME_DIR=/home/leoner
PROJ_DIR=${HOME_DIR}/Projects/cat_whim
BIDS_DIR=${PROJ_DIR}/data/interim/bids
DER_DIR=${BIDS_DIR}/derivatives
CODE_DIR=${PROJ_DIR}/cat_whim

SUBJECT_N=$( sed -n -E "$((${SLURM_ARRAY_TASK_ID} + 1))s/sub-(\S*)\>.*/\1/gp" ${BIDS_DIR}/participants.tsv )
SUBJECT_ID="sub-"${SUBJECT_N}

echo "Processing ${SUBJECT_ID}"
wmh_script=${CODE_DIR}/wmh_segmentation.py
. ${HOME_DIR}/miniconda3/etc/profile.d/conda.sh

conda activate cat_whim
echo ">>> Activated cat_whim environment"
python ${wmh_script} ${SUBJECT_ID}

