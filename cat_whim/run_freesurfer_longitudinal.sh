#!/bin/bash
####    Run FreeSurfer Longitudinal: Slurm -- Version 1.0.0
# Last edit:  2024-08-14
# Authors:    Leone, Riccardo (RL)
# Notes:      - Script to run FreeSurfer longitudinal pipeline
#               for data preprocessed with sMRIPrep
#             - Release notes:
#                 * Initial release for longitudinal data
# To do:      - 
# Comments:   
# Sources:    
####


# SLURM settings
#SBATCH --job-name="freesurfer_long"
#SBATCH --time=3-12:00:00
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
DATA_DIR=${PROJ_DIR}/data/interim
BIDS_DIR=${DATA_DIR}/bids
DER_DIR=${BIDS_DIR}/derivatives
FS_DIR=${DER_DIR}/freesurfer
DOCS_DIR=${PROJ_DIR}/docs

if [ ! -d "${FS_DIR}" ]; then
  mkdir -p "${FS_DIR}"
fi

SUBJECT_ID=$( sed -n -E "$((${SLURM_ARRAY_TASK_ID} + 1))s/sub-(\S*)\>.*/\1/gp" ${BIDS_DIR}/participants.tsv )

echo "Processing ${SUBJECT_ID}"

# Make temporary directories
GLOB_TEMP_DIR=${DATA_DIR}/temp_dir
SUBJ_TEMP_DIR=${GLOB_TEMP_DIR}/${SUBJECT_ID}

mkdir -p ${GLOB_TEMP_DIR}
mkdir -p ${SUBJ_TEMP_DIR}

export SINGIMAGE=${PROJ_DIR}/singims/freesurfer_7.4.1-202309.sif
export LICENSE=${DOCS_DIR}/license.txt

# Load Singularity in order to run the sMRIPrep container containing FreeSurfer
module load singularity

singularity_cmd="singularity run --cleanenv -B $BIDS_DIR:/data \
    -B $FS_DIR:/out -B $SUBJ_TEMP_DIR:/work \
    -B $LICENSE:/license.txt:ro $SINGIMAGE"

cmd="${singularity_cmd} /data /out participant --participant_label $SUBJECT_ID \
  --license_file /license.txt \
  --skip_bids_validator \
  --measurements area volume thickness thicknessstd meancurv\
  --refine_pial FLAIR \
  --steps cross-sectional template longitudinal"

echo Running task ${SLURM_ARRAY_TASK_ID}
echo Commandline: $cmd
eval $cmd
exitcode=$?
echo Finished task ${SLURM_ARRAY_TASK_ID} with exit code:  $exitcode
rm -r ${SUBJ_TEMP_DIR}
echo "${SUBJECT_ID}  ${exitcode}" >> ${FS_DIR}/exit_codes_fs.tsv
exit $exitcode