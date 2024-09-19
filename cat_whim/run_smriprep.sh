#!/bin/bash
####    Run sMRIPrep: Slurm -- Version 1.0.0
# Last edit:  2024-08-12
# Authors:    Leone, Riccardo (RL)
# Notes:      - Script to run sMRIPrep on bids data from ADNI database / longitudinal
#             - Release notes:
#                 * Initial release for longitudinal data
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
#SBATCH --time=24:00:00
#SBATCH --output=/home/leoner/Projects/cat_whim/slurmlogs/outputs/%x-%A-%a.out
#SBATCH --error=/home/leoner/Projects/cat_whim/slurmlogs/errors/%x-%A-%a.err
#SBATCH --nodes=1
#SBATCH --array=1-730
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=10G
#SBATCH --partition=HPC-CPUs
#SBATCH --ntasks=1

module load singularity

HOME_DIR=/home/leoner
PROJ_DIR=${HOME_DIR}/Projects/cat_whim

DATA_DIR=${PROJ_DIR}/data/interim
BIDS_DIR=${DATA_DIR}/bids
DER_DIR=${BIDS_DIR}/derivatives
FSURF_DIR=${DER_DIR}/freesurfer

DER_DIR_M00=${DER_DIR}/ses-M00
DER_DIR_M01=${DER_DIR}/ses-M01
DER_DIR_M02=${DER_DIR}/ses-M02
DER_DIR_M03=${DER_DIR}/ses-M03

FSURF_DIR_M00=${DER_DIR_M00}/freesurfer
FSURF_DIR_M01=${DER_DIR_M01}/freesurfer
FSURF_DIR_M02=${DER_DIR_M02}/freesurfer
FSURF_DIR_M03=${DER_DIR_M03}/freesurfer

DOCS_DIR=${PROJ_DIR}/docs
SCRIPTS_DIR=${PROJ_DIR}/cat_whim

include_ses_M00=${SCRIPTS_DIR}/include_ses_M00.json
include_ses_M01=${SCRIPTS_DIR}/include_ses_M01.json
include_ses_M02=${SCRIPTS_DIR}/include_ses_M02.json
include_ses_M03=${SCRIPTS_DIR}/include_ses_M03.json

subjectsfile=${BIDS_DIR}/participants.tsv

export LICENSE=${DOCS_DIR}/license.txt
export SINGMAGE=${PROJ_DIR}/singims/smriprep-0.16.0.sif

# Set up directories and variables
setup_dirs() {
    mkdir -p ${FSURF_DIR_M00}
    mkdir -p ${FSURF_DIR_M01}
    mkdir -p ${FSURF_DIR_M02}
    mkdir -p ${FSURF_DIR_M03}

    mkdir -p ${DATA_DIR}/temp_dir
}

#### Function to run sMRIPrep for a given subject and session
run_smriprep() {
    local subject_id=$1
    local session=$2
    local der_dir=$3
    local fsurf_dir=$4
    local include_json=$5

    local temp_dir=${DATA_DIR}/temp_dir/${subject_id}_smriprep
    mkdir -p ${temp_dir}

    local singularity_cmd="singularity run --cleanenv -B $BIDS_DIR:/data \
        -B $der_dir:/out -B $temp_dir:/work \
        -B $LICENSE:/freesurfer_license.txt:ro $SINGMAGE"

    local cmd_cross="${singularity_cmd} /data /out participant --participant-label $subject_id \
        --fs-license-file $LICENSE --notrack \
        --fs-subjects-dir /data/derivatives/${session}/freesurfer \
        --omp-nthreads 8 --nthreads 12 \
        --bids-filter-file ${include_json} \
        --output-spaces MNI152NLin6Asym:res-01 anat \
        --stop-on-first-crash -w /work"

    echo "Moving ${FSURF_DIR}/sub-${subject_id}_${session} to ${fsurf_dir}/sub-${subject_id}"
    mv "${FSURF_DIR}/sub-${subject_id}_${session}" "${fsurf_dir}/sub-${subject_id}"

    echo "Removing any IsRunning instances for session ${session}"
    find "${fsurf_dir}/sub-${subject_id}/" -name "*IsRunning*" -type f -delete

    echo "Running sMRIPrep for session ${session}:"
    echo "Commandline: $cmd_cross"
    eval $cmd_cross
    return $?
}


# Main script
main() {
    setup_dirs

    local subject_id=$(sed -n -E "$((${SLURM_ARRAY_TASK_ID} + 1))s/sub-(\S*)\>.*/\1/gp" ${BIDS_DIR}/participants.tsv)
    echo "Processing subject: $subject_id"

    run_smriprep $subject_id "ses-M00" $DER_DIR_M00 $FSURF_DIR_M00 $include_ses_M00
    local exitcode_m00=$?

    run_smriprep $subject_id "ses-M01" $DER_DIR_M01 $FSURF_DIR_M01 $include_ses_M01
    local exitcode_m01=$?

    run_smriprep $subject_id "ses-M02" $DER_DIR_M02 $FSURF_DIR_M02 $include_ses_M02
    local exitcode_m02=$?

    run_smriprep $subject_id "ses-M03" $DER_DIR_M03 $FSURF_DIR_M03 $include_ses_M03
    local exitcode_m03=$?


    echo "$subject_id M00 $exitcode_m00 M01 $exitcode_m01 M02 $exitcode_m02 M03 $exitcode_m03" \
        >> ${DER_DIR}/smriprep_exit_codes.tsv

    echo "Finished task ${SLURM_ARRAY_TASK_ID} with exit codes: M00 $exitcode_m00, M01 $exitcode_m01, M02 $exitcode_m02, M03 $exitcode_m03"

    exit $((exitcode_m00 + exitcode_m01 + exitcode_m02 + exitcode_m03))
}

# Run the main function
main
