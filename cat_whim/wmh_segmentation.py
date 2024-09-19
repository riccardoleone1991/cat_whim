# %%
import sys
import ants
import antspynet
import nibabel as nib

from loguru import logger

from wmh_segmentation_utils import *


def main(subj_name):
    WMH_SUBJ_DIR = WMH_DIR / subj_name

    SUBJ_SES00_DIR = BIDS_DATA_DIR / "derivatives" / "ses-M00" / "smriprep" / subj_name
    SUBJ_SES01_DIR = BIDS_DATA_DIR / "derivatives" / "ses-M01" / "smriprep" / subj_name
    SUBJ_SES02_DIR = BIDS_DATA_DIR / "derivatives" / "ses-M02" / "smriprep" / subj_name
    SUBJ_SES03_DIR = BIDS_DATA_DIR / "derivatives" / "ses-M03" / "smriprep" / subj_name

    list_ses_subj_dir = [SUBJ_SES00_DIR, SUBJ_SES01_DIR, SUBJ_SES02_DIR, SUBJ_SES03_DIR]

    if not Path.exists(WMH_SUBJ_DIR):
        Path.mkdir(WMH_SUBJ_DIR)

    for subj_dir in list_ses_subj_dir:
    
        for ses_dir in subj_dir.iterdir():

            ses_name = ses_dir.name

            if ses_name.startswith("ses-"):
                logger.info(f"Processing session {ses_name}")
                WMH_SUBJ_SES_DIR = WMH_SUBJ_DIR / ses_name

                if not Path.exists(WMH_SUBJ_SES_DIR):
                    Path.mkdir(WMH_SUBJ_SES_DIR)

                # Filenames for loading (strings because ants does not like Path)
                t1_file = str(
                    BIDS_DATA_DIR
                    / subj_name
                    / ses_name
                    / "anat"
                    / f"{subj_name}_{ses_name}_T1w.nii.gz"
                )
                t1_bm_file = str(
                    ses_dir / "anat" / f"{subj_name}_{ses_name}_desc-brain_mask.nii.gz"
                )
                flair_file = str(
                    BIDS_DATA_DIR
                    / subj_name
                    / ses_name
                    / "anat"
                    / f"{subj_name}_{ses_name}_FLAIR.nii.gz"
                )
                t1_mni6_file = str(
                    ses_dir
                    / "anat"
                    / f"{subj_name}_{ses_name}_space-MNI152NLin6Asym_res-01_desc-preproc_T1w.nii.gz"
                )
                t1_mni6_bm_file = str(
                    ses_dir
                    / "anat"
                    / f"{subj_name}_{ses_name}_space-MNI152NLin6Asym_res-01_desc-brain_mask.nii.gz"
                )

                # Filenames for saving in subject and MNI152NLin6Asym space
                wmh_mask_file = str(
                    WMH_SUBJ_SES_DIR / f"{subj_name}_{ses_name}_label-wmh_mask_desc-raw.nii.gz"
                )
                wmh_mask_postproc_file = str(
                    WMH_SUBJ_SES_DIR
                    / f"{subj_name}_{ses_name}_label-wmh_mask_desc-postprocessed.nii.gz"
                )
                wmh_mask_mni6_file = str(
                    WMH_SUBJ_SES_DIR
                    / f"{subj_name}_{ses_name}_space-MNI152NLin6Asym_res-01_label-wmh_mask_desc-raw.nii.gz"
                )
                wmh_mask_mni6_postproc_file = str(
                    WMH_SUBJ_SES_DIR
                    / f"{subj_name}_{ses_name}_space-MNI152NLin6Asym_res-01_label-wmh_mask_desc-postprocessed.nii.gz"
                )

                # Load T1 file and T1 brain mask file in subject space and perform brain extraction
                t1 = ants.image_read(t1_file)
                t1_bm = ants.image_read(t1_bm_file)
                t1_bet = ants.mask_image(t1, t1_bm)

                # Load FLAIR file in subject space and perform brain extraction
                flair = ants.image_read(flair_file)
                flair_prob_bm = antspynet.brain_extraction(flair, modality="flair")
                flair_bm = ants.get_mask(flair_prob_bm, low_thresh=0.5)
                flair_bet = ants.mask_image(flair, flair_bm)

                # Load T1 and brain mask file in MNI6 space and perform brain extraction
                t1_mni6 = ants.image_read(t1_mni6_file)
                t1_mni6_bm = ants.image_read(t1_mni6_bm_file)
                t1_mni6_bet = ants.mask_image(t1_mni6, t1_mni6_bm)

                # Register FLAIR to T1 in subject space
                flair_to_t1_registration = ants.registration(
                    fixed=t1_bet, moving=flair_bet, type_of_transform="SyN"
                )
                flair_bet_reg = flair_to_t1_registration["warpedmovout"]

                # Perform wmh segmentation in subject space
                logger.info("Performing WMH segmentation in subject space...")
                wmh_mask_prob = antspynet.sysu_media_wmh_segmentation(
                    flair=flair_bet_reg, t1=t1_bet, verbose=True
                )
                # Binarize, convert to Nibabel
                wmh_mask_bin = ants.get_mask(wmh_mask_prob, low_thresh=0.9, cleanup=-1)
                # wmh_mask_bin_nib = to_nibabel(wmh_mask_bin)

                # Now we want to register the WMH mask in MNI152NLin6Asym space for later processing with LQT
                # Perform registration of T1 to MNI6 space
                logger.info("Registering T1 to MNI space...")
                registration = ants.registration(
                    fixed=t1_mni6_bet, moving=t1_bet, type_of_transform="SyN"
                )
                # Apply transformation to the WMH mask
                logger.info("Applying transformation to WMH mask")
                wmh_mni6_mask_bin = ants.apply_transforms(
                    moving=wmh_mask_bin,
                    fixed=registration["warpedmovout"],
                    transformlist=registration["fwdtransforms"],
                    interpolator="genericLabel",
                )
                wmh_mni6_mask_bin_nib = to_nibabel(wmh_mni6_mask_bin)
                wmh_mni6_mask_bin_nib_postproc = postproc_wmh_mask_mni6(wmh_mni6_mask_bin_nib)
                wmh_mni6_mask_bin_postproc = from_nibabel(wmh_mni6_mask_bin_nib_postproc)
                logger.info("Registering post-processed WMH mask to subject space")
                # Retransform the postprocessed to subject space
                wmh_mask_postproc_bin = ants.apply_transforms(
                    moving=wmh_mni6_mask_bin_postproc,
                    fixed=t1_bet,
                    transformlist=registration["invtransforms"],
                    interpolator="genericLabel",
                )

                # wmh_mask_bin_postproc_nib = to_nibabel(wmh_mask_postproc_bin)
                logger.info("Saving created files...")
                # Save both the raw and postprocessed maps for comparison in subject and MNI space
                wmh_mask_bin.to_filename(wmh_mask_file)
                wmh_mask_postproc_bin.to_filename(wmh_mask_postproc_file)
                wmh_mni6_mask_bin.to_filename(wmh_mask_mni6_file)
                wmh_mni6_mask_bin_postproc.to_filename(wmh_mask_mni6_postproc_file)
                logger.success(f"Succesfully processed {subj_name}!")


subj_name = sys.argv[1]
if __name__ == "__main__":
    main(subj_name)
    logger.info("Cleaning unused temporary files...")
    clear_tmp_dir(TMP_DIR)
    logger.success(f"WMH segmentation for {subj_name} fully completed!")
