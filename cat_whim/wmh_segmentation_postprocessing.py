#%%
import shutil
import numpy as np
import pandas as pd
import nibabel as nib
import nibabel.imagestats as nibstats

from tqdm import tqdm
from loguru import logger

from wmh_segmentation_utils import WMH_DIR, INTERIM_DATA_DIR, BIDS_DATA_DIR, INCL_EXCL_FLW_DIR


def is_there_sequence(sequence_name):
    try: 
        nib.load(sequence_name)
        return True
    except:
        logger.warning(f"{sequence_name} not found")
        return False


def check_errors_smriprep():

    logger.info("Checking if there were errors in smriprep..")

    list_errors_smriprep = []

    for subj_dir in WMH_DIR.iterdir():
        if subj_dir.is_dir() and subj_dir.name.startswith("sub-"):
            subj = subj_dir.name
            for ses_dir in subj_dir.iterdir():
                ses = ses_dir.name

                if ses == "ses-M00":
                    SMRIPREP_DIR = BIDS_DATA_DIR / "derivatives" / "ses-M00" / "smriprep"
                
                elif ses == "ses-M01":
                    SMRIPREP_DIR = BIDS_DATA_DIR / "derivatives" / "ses-M01" / "smriprep"
                
                elif ses == "ses-M02":
                    SMRIPREP_DIR = BIDS_DATA_DIR / "derivatives" / "ses-M02" / "smriprep"
                
                elif ses == "ses-M03":
                    SMRIPREP_DIR = BIDS_DATA_DIR / "derivatives" / "ses-M03" / "smriprep"
                

                T1_dir = SMRIPREP_DIR / subj / ses / "anat"
                t1_filename = T1_dir / f"{subj}_{ses}_desc-preproc_T1w.nii.gz"
                t1_filename_mni = T1_dir / f"{subj}_{ses}_space-MNI152NLin6Asym_res-01_desc-preproc_T1w.nii.gz"
                mask_filename = ses_dir / f"{subj}_{ses}_label-wmh_mask_desc-postprocessed.nii.gz"
                mask_filename_mni = ses_dir / f"{subj}_{ses}_space-MNI152NLin6Asym_res-01_label-wmh_mask_desc-postprocessed.nii.gz"

                is_seq = is_there_sequence(t1_filename_mni)
                if not is_seq:
                    list_errors_smriprep.append(subj)

    n_errors_smriprep = len(list_errors_smriprep)
    df = pd.read_csv(INCL_EXCL_FLW_DIR / "df_inclusion_exclusion_longitudinal_tau_amyloid.csv")
    df["N errors smriprep"] = n_errors_smriprep
    df.to_csv(INCL_EXCL_FLW_DIR / "df_inclusion_exclusion_longitudinal_tau_amyloid.csv")

    logger.warning(f"There were {n_errors_smriprep} errors in smriprep. I will not process these subjects further")

    return list_errors_smriprep

def postprocess_wmh_segmentations(list_errors_smriprep):
    list_ptid = []
    list_ses = []
    list_wmh_vol = []

    len_subjs = len([subj for subj in WMH_DIR.iterdir() if subj.name.startswith("sub-") and subj.is_dir()])
    logger.info("I will create empty masks for subjects who don't have WMH and calculate the WMH volumes for those that do...")
    with tqdm(total=len_subjs) as pbar:
        for subj_dir in WMH_DIR.iterdir():
            if subj_dir.is_dir() and subj_dir.name.startswith("sub-"):
                subj = subj_dir.name
                #If there were no errors in processing, continue
                if subj not in list_errors_smriprep:
                    for ses_dir in subj_dir.iterdir():
                        ses = ses_dir.name

                        if ses == "ses-M00":
                            SMRIPREP_DIR = BIDS_DATA_DIR / "derivatives" / "ses-M00" / "smriprep"
                        
                        elif ses == "ses-M01":
                            SMRIPREP_DIR = BIDS_DATA_DIR / "derivatives" / "ses-M01" / "smriprep"
                        
                        elif ses == "ses-M02":
                            SMRIPREP_DIR = BIDS_DATA_DIR / "derivatives" / "ses-M02" / "smriprep"
                        
                        elif ses == "ses-M03":
                            SMRIPREP_DIR = BIDS_DATA_DIR / "derivatives" / "ses-M03" / "smriprep"
                        
                        T1_dir = SMRIPREP_DIR / subj / ses / "anat"
                        t1_filename = T1_dir / f"{subj}_{ses}_desc-preproc_T1w.nii.gz"
                        t1_filename_mni = T1_dir / f"{subj}_{ses}_space-MNI152NLin6Asym_res-01_desc-preproc_T1w.nii.gz"

                        mask_filename = ses_dir / f"{subj}_{ses}_label-wmh_mask_desc-postprocessed.nii.gz"
                        mask_filename_mni = ses_dir / f"{subj}_{ses}_space-MNI152NLin6Asym_res-01_label-wmh_mask_desc-postprocessed.nii.gz"
                        
                        try:
                            # Create empty mask in subject space
                            wmh = nib.load(mask_filename)
                            vol = nibstats.mask_volume(wmh)
                            list_ptid.append(subj)
                            list_ses.append(ses)
                            list_wmh_vol.append(vol)
                        except:
                            logger.info(f"There is no mask for {subj}, {ses}, creating an empty mask in subject space...")
                            t1 = nib.load(t1_filename)
                            t1_affine = t1.affine
                            t1_arr = t1.get_fdata()
                            msk = np.zeros_like(t1_arr)
                            msk_img = nib.Nifti1Image(msk, t1_affine)
                            msk_img.to_filename(mask_filename)
                            list_ptid.append(subj)
                            list_ses.append(ses)
                            list_wmh_vol.append(0)
                        try:
                            nib.load(mask_filename_mni)
                        except:
                            # Create empty mask in MNI space
                            logger.info(f"There is no mask for {subj}, {ses}, creating an empty mask in MNI space...")
                            t1_mni = nib.load(t1_filename_mni)
                            t1_mni_affine = t1_mni.affine
                            t1_mni_arr = t1_mni.get_fdata()
                            msk_mni = np.zeros_like(t1_mni_arr)
                            msk_img_mni = nib.Nifti1Image(msk_mni, t1_mni_affine)
                            msk_img_mni.to_filename(mask_filename_mni)
            pbar.update(1)
        df = pd.DataFrame({"PTID": list_ptid, "session": list_ses, "wmh_vol": list_wmh_vol})
        df.to_csv(INTERIM_DATA_DIR / "df_wmh_vols.csv", index=None)


LQT_DIR = BIDS_DATA_DIR / "derivatives" / "LQT"
LQT_MASKS_DIR = LQT_DIR / "lesion_masks"

if not LQT_DIR.exists():
    LQT_DIR.mkdir()

if not LQT_MASKS_DIR.exists():
    LQT_MASKS_DIR.mkdir()


list_errors_smriprep = check_errors_smriprep()
for subj_dir in WMH_DIR.iterdir():
    if subj_dir.is_dir() and subj_dir.name.startswith("sub-"):
        subj = subj_dir.name
        #If there were no errors in processing, continue
        if subj not in list_errors_smriprep:
            for ses_dir in subj_dir.iterdir():
                ses = ses_dir.name
                filename = f"{subj}_{ses}_space-MNI152NLin6Asym_res-01_label-wmh_mask_desc-postprocessed.nii.gz"
                old_name = ses_dir / filename
                new_name = LQT_MASKS_DIR / filename
                shutil.copy(old_name, new_name)

if __name__ == "__main__":
    logger.info("Starting to post-process WMH data.")
    list_err_smriprep = check_errors_smriprep()
    postprocess_wmh_segmentations(list_err_smriprep)
    logger.success("WMH post-processing succesfully completed!")