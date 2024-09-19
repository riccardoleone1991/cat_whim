#!/bin/bash

# Define paths
HOME_DIR=/home/leoner
PROJ_DIR=${HOME_DIR}/Projects/cat_whim
DATA_DIR=${PROJ_DIR}/data/interim
BIDS_DIR=${DATA_DIR}/bids
DER_DIR=${BIDS_DIR}/derivatives
FS_DIR=${DER_DIR}/freesurfer

SINGULARITY_IMG=${PROJ_DIR}/singims/freesurfer_7.4.1-202309.sif
QDEC_OUTPUT_FILE=${FS_DIR}/qdec_table.dat

# Create the qdec file for easier gathering of stats
echo "fsid fsid-base" > $QDEC_OUTPUT_FILE
for subject in $FS_DIR/sub-*M0{0,1,2,3,4}; do
    # Extract the subject ID and session (assuming directory name format is sub-XX_ses-YY)
    subject=$(basename $subject)
    subject_id=${subject%%.long*}
    base_id="${subject_id%_ses-*}"
    echo "$subject_id $base_id" >> $QDEC_OUTPUT_FILE
done

# Load the Singularity module
module load singularity 

# Execute the FreeSurfer commands within the Singularity container
singularity exec -B $FS_DIR:/data $SINGULARITY_IMG /bin/bash -c "\
    export SUBJECTS_DIR=/data && \
    aparcstats2table --qdec-long ${QDEC_OUTPUT_FILE} --hemi lh --meas thickness --tablefile /data/aparc_lh_table.tsv --skip && \
    aparcstats2table --qdec-long ${QDEC_OUTPUT_FILE} --hemi rh --meas thickness --tablefile /data/aparc_rh_table.tsv --skip"

# asegstats2table --qdec-long ${QDEC_OUTPUT_FILE} --stats aseg.stats --meas thickness --tablefile aseg_table.tsv --skip && \
# cp ${FS_DIR}/aseg_table.tsv ${DATA_DIR}/aseg_table.tsv
cp ${FS_DIR}/aparc_lh_table.tsv ${DATA_DIR}/aparc_lh_table.tsv
cp ${FS_DIR}/aparc_rh_table.tsv ${DATA_DIR}/aparc_rh_table.tsv