#!/bin/bash
####    Run fMRIPrep: Slurm -- Version 2.1.2
# Last edit:  2024-08-12
# Authors:    Leone, Riccardo (RL)
# Notes:      - Script to run sMRIPrep on bids data from ADNI database
#             - Release notes:
#                 * Initial release
# To do:      - 
# Comments:   
# Sources:    https://github.com/andrewjahn/OpenScience_Scripts/blob/master/fmriprep_Scripted.sh
#             https://andysbrainbook.readthedocs.io/en/latest/OpenScience/OS/fMRIPrep.html
#             https://blog.ronin.cloud/slurm-job-arrays/
#             https://dartbrains.org/content/fmriprep_on_discovery.html
#             https://www.nipreps.org/apps/singularity/
####


# SLURM settings
#SBATCH --job-name="smriprep"
#SBATCH --time=01:00:00
#SBATCH --output=/home/leoner/Projects/cat_whim/slurmlogs/outputs/%x-%A-%a.out
#SBATCH --error=/home/leoner/Projects/cat_whim/slurmlogs/errors/%x-%A-%a.err
#SBATCH --nodes=1
#SBATCH --array=1-2
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=HPC-CPUs
#SBATCH --ntasks=1


# module load singularity

HOME_DIR=/home/leoner
PROJ_DIR=${HOME_DIR}/Projects/cat_whim

DATA_DIR=${PROJ_DIR}/data/interim
BIDS_DIR=${DATA_DIR}/bids
DER_DIR=${BIDS_DIR}/derivatives
DOCS_DIR=${PROJ_DIR}/docs

mkdir -p ${DER_DIR}

subjectsfile=${BIDS_DIR}/participants.tsv
export LICENSE=${DOCS_DIR}/license.txt
export SINGMAGE=${PROJ_DIR}/singims/smriprep-0.16.0.sif

# Parse the participants.tsv file and extract one subject ID from the line corresponding to this SLURM task.
# subject=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $subjectsfile)
subject_n=$( sed -n -E "$((${SLURM_ARRAY_TASK_ID} + 1))s/sub-(\S*)\>.*/\1/gp" ${BIDS_DIR}/participants.tsv )
subject="ADNI"${subject_n}
echo $subject
#find ${fsurf_dir}/$subject/ -name "*IsRunning*" -type f -delete

# Make temporary directories
GLOB_TEMP_DIR=${DATA_DIR}/temp_dir
SUBJ_TEMP_DIR=${GLOB_TEMP_DIR}/${subject}

mkdir -p ${GLOB_TEMP_DIR}
mkdir -p ${SUBJ_TEMP_DIR}

# Mount directories
singularity_cmd="singularity run --cleanenv -B $BIDS_DIR:/data:ro \
    -B $DER_DIR:/out -B $SUBJ_TEMP_DIR:/work \
    -B $LICENSE:/freesurfer_license.txt:ro $SINGMAGE"

# sMRIPrep flags 
cmd="${singularity_cmd} /data /out participant --participant-label $subject \
    --fs-license-file $LICENSE --notrack \
    --longitudinal \
    --omp-nthreads 8 --nthreads 12 \
    --output-spaces MNI152NLin6Asym:res-01 anat \
    --stop-on-first-crash -w /work"

# Run sMRIPrep in Singularity
echo Running task ${SLURM_ARRAY_TASK_ID}
echo Commandline: $cmd
eval $cmd
exitcode=$?

# echo "sub-$subject ${SLURM_ARRAY_TASK_ID} $exitcode" \
#     >> ${SLURM_JOB_NAME}${SLURM_ARRAY_JOB_ID}${subject}.tsv
echo Finished task ${SLURM_ARRAY_TASK_ID} with exit code $exitcode
exit $exitcode
# ------------------------------------------------------------------------ End