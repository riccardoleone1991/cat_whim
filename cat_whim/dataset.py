# %%
import csv
import subprocess
import shutil
import pandas as pd
import numpy as np

from pathlib import Path
from loguru import logger
from tqdm import tqdm

from cat_whim.config import (
    RAW_DATA_DIR,
    INTERIM_DATA_DIR,
    BIDS_DATA_DIR,
    SCRIPTS_DIR,
    SINGIMS_DIR,
    INCL_EXCL_FLW_DIR,
)


def download_dataset_adni(
        input_path: Path = RAW_DATA_DIR / "cat_whim_url_list_2024-08-21.csv",
        out_path: Path = INTERIM_DATA_DIR):

    """
    This function downloads the urls from the url_list into cat_whim/data/raw. 
    Then, it copies the folder into cat_whim/data/interim.
    """
    # Read the entire CSV file into a list
    with open(input_path, mode="r") as file:
        reader = csv.reader(file)
        lines = list(reader)

    for i in range(len(lines)):
        # Extract the URL from the current line
        url = lines[i][0]

        # Extract the last part of the URL
        file_name = url.split("/")[-1]
        print(f"Downloading file: {file_name}")

        # Download the file using wget to RAW_DIR
        subprocess.run(["wget", "-P", str(RAW_DATA_DIR), url])

        # Unzip the downloaded file if it's a zip file
        if file_name.endswith(".zip"):
            print(f"unzipping {file_name}")
            subprocess.run(["unzip", "-d", str(RAW_DATA_DIR), str(RAW_DATA_DIR / file_name)])
            logger.info(f"Extracted {file_name} successfully!")
            print(f"removing {file_name} folder")
            Path.unlink(RAW_DATA_DIR / file_name)
    logger.info("Copying downloaded folder to data/interim")
    cmd = ["cp", "-r", f"{RAW_DATA_DIR}/ADNI", f"{out_path}"]
    # Run the command
    subprocess.run(cmd)


def check_download_match(
    input_path: Path = RAW_DATA_DIR / "ADNI",
    input_file: Path = INCL_EXCL_FLW_DIR / "subjects_to_download.csv",
    in_out_file: Path = INCL_EXCL_FLW_DIR / "df_inclusion_exclusion_longitudinal_tau_amyloid.csv"
):
    """
    Check that all data were downloaded correctly based on a list of subjects to download.

    Parameters:
    input_path (Path): Path to a folder where data was downloaded
    input_file (Path): Path to a .csv file containing all subjects to be downloaded
    in_out_file (Path): Path to a .csv file containing inclusion and exclusion numbers
    
    Returns:
    None

    """
    list_downloaded_subjs = []
    # Create the list of subjects you downloaded
    for dir in input_path.iterdir():
        subj = dir.name
        list_downloaded_subjs.append(subj)

    # Load the csv file with the names of the subjects with available tau and amyloid data
    df_subjs_to_download = pd.read_csv(input_file, header=None)
    df_subjs_to_download.columns = ["PTID"]
    # Check that all subjects you had to download are in the downloaded folder
    list_bool_check_download = [
        subj in list_downloaded_subjs for subj in df_subjs_to_download["PTID"]
    ]
    # The sum of the boolean list should be equal to the number of subjects to download
    is_download_complete = sum(list_bool_check_download) == df_subjs_to_download["PTID"].shape[0]
    list_missing = [
        subj for subj in df_subjs_to_download["PTID"] if not subj in list_downloaded_subjs
    ]

    if not is_download_complete:
        logger.warning(
            f"{len(list_missing)} subjects went missing when downloading... {list_missing}. You should check what happened! (e.g., maybe the subject doesn't have any 3D T1 or FLAIR)"
        )
        df_incl_excl = pd.read_csv(in_out_file)
        df_incl_excl["excluded_after_download"] = len(list_missing)
        df_incl_excl.to_csv(in_out_file, index=None)


def load_df_amy_tau():
    return pd.read_csv(INTERIM_DATA_DIR / "df_amyBoth_tau_to_download.csv")


def load_df_adnimerge():
    return pd.read_csv(RAW_DATA_DIR / "ADNIMERGE_18Jul2024.csv", dtype="object")


def load_df_mri():

    df_images = pd.read_csv(RAW_DATA_DIR / "cat_whim_8_21_2024.csv")
    df_images_filtered = df_images[
        df_images["Description"].str.lower().str.contains("t1|rage|spgr|flair", regex=True)
    ].copy()
    df_images_filtered["Acq Date"] = pd.to_datetime(df_images_filtered["Acq Date"])
    df_images_filtered = df_images_filtered.rename(
        columns={"Subject": "PTID", "Acq Date": "MRIDATE"}
    )
    return df_images_filtered[["PTID", "MRIDATE", "Description"]]


def filter_adnimerge_based_on_pet(df_pet, df_adni):
    """
    Filter the adnimerge dataframe to keep only those subjects/timepoints for which we have PET data.

    Parameters:
    df_pet (DataFrame): PET data.
    df_adni (DataFrame): Clinical data from adnimerge.

    Returns:
    DataFrame: Merged dataframe.
    """

    df_adni["EXAMDATE"] = pd.to_datetime(df_adni["EXAMDATE"])
    df_pet["SCANDATE_amy"] = pd.to_datetime(df_pet["SCANDATE_amy"])

    # Select the columns from the original adnimerge
    adni_cols = [col for col in df_adni.columns if col != "VISCODE"]
    # Add the session column that we want to keep for reference between datasets
    adni_cols = adni_cols + ["session"]
    # Merge PET and ADNIMERGE
    merged_df = pd.merge(df_pet, df_adni, on="PTID", suffixes=("_pet", "_adnimerge"))
    # Calculate the difference between the dates in days
    merged_df["date_diff"] = np.abs((merged_df["SCANDATE_amy"] - merged_df["EXAMDATE"]).dt.days)
    # We consider as the same timepoints sessions that were performed less than 6 months apart
    result_df = merged_df[merged_df["date_diff"] <= 180]
    # Recreate the adnimerge only for data we want to include
    df_adni_new = result_df.drop(columns=["VISCODE_adnimerge"])[adni_cols]
    df_adni_new.to_csv(INTERIM_DATA_DIR / "df_adnimerge_filtered_pet.csv", index=None)
    return df_adni_new


def filter_mri_based_on_pet(df_pet, df_mri):
    """
    Filter the MRI dataframe to keep only those subjects/timepoints for which we have PET data.

    Parameters:
    df_pet (DataFrame): PET data.
    df_mri (DataFrame): Downloaded MRI data.

    Returns:
    DataFrame: Merged dataframe.
    """

    df_pet["SCANDATE_amy"] = pd.to_datetime(df_pet["SCANDATE_amy"])
    df_mri["MRIDATE"] = pd.to_datetime(df_mri["MRIDATE"])

    # Select the columns from the original adnimerge
    mri_cols = [col for col in df_mri.columns if col != "VISCODE"]
    # Add the session column that we want to keep for reference between datasets
    mri_cols = mri_cols + ["session"]
    # Merge PET and ADNIMERGE
    merged_df = pd.merge(df_pet, df_mri, on="PTID")
    # Calculate the difference between the dates in days
    merged_df["date_diff"] = np.abs((merged_df["SCANDATE_amy"] - merged_df["MRIDATE"]).dt.days)
    # We consider as the same timepoint, sessions that were performed less than 6 months apart
    result_df = merged_df[merged_df["date_diff"] <= 180]
    # Recreate the adnimerge only for data we want to include
    df_mri_new = result_df[mri_cols]
    return df_mri_new


def select_complete_mri(df_filtered):
    """
    Preprocess dataframe of downloaded images to keep only subjects with all MRI data needed --> T1 and FLAIR.

    Parameters:
    df_filtered (DataFrame): Df containing downloaded images filtered for subjects that also have tau /amyloid pet.

    Returns:
    DataFrame: Filtered df containing subjects/sessions combinations that have at least all the 2 sequences
    """

    def check_description(description_list):
        lower_i = [el.lower() for el in description_list]
        if (
            any("spgr" in el for el in lower_i)
            or any("rage" in el for el in lower_i)
            or any("t1" in el for el in lower_i)
        ):
            if any("flair" in el for el in lower_i):
                return True
            else:
                return False
        else:
            return False

    df_counts = df_filtered.groupby(["PTID", "MRIDATE"]).count()
    # Count those PTID/MRIDATE combinations that have less than 2 sequences
    less_than_two = df_counts[df_counts["Description"] <= 1].index
    # Select only subjects that have 2 or more sequence for each timepoint 
    df_filtered = df_filtered[
        ~df_filtered.set_index(["PTID", "MRIDATE"]).index.isin(less_than_two)
    ]

    result = df_filtered.groupby(["PTID", "MRIDATE"])["Description"].unique()
    check_results = result.apply(check_description)
    df_check_results = check_results.reset_index(name="Contains_T1_and_Flair")
    df_mri_complete = df_check_results[df_check_results["Contains_T1_and_Flair"] == True]
    df_mri_complete_session = (
        pd.merge(
            df_mri_complete,
            df_filtered[["PTID", "MRIDATE", "session", "Description"]],
            on=["PTID", "MRIDATE"],
        )
        .drop_duplicates()
        .reset_index(drop=True)
    )

    return df_mri_complete_session


def perform_checks(df_amy_tau, df_adni_filtered, df_mri_filtered_complete):

    n_unique_amy_tau = df_amy_tau['PTID'].nunique()
    n_unique_adnimerge = df_adni_filtered['PTID'].nunique()
    n_unique_mri = df_mri_filtered_complete['PTID'].nunique()
    
    if n_unique_adnimerge < n_unique_mri:
        n_min = n_unique_adnimerge
        which_df = "adnimerge"

    else:
        n_min = n_unique_mri
        which_df = "mri"

    logger.info(
        f"The number of subjects initially with amyloid and tau at the same timepoint is: {n_unique_amy_tau}"
    )   

    # Check which df has the minimum number of subjects, then save the number of subjects that get excluded
    # due to the fact that they don't have clinical or imaging at the same timepoint as PET.
    logger.info(
        f"That of subjects with complete adnimerge is: {n_unique_adnimerge}"
    )

    logger.info(
        f"That of subjects with complete MRI is: {n_unique_mri}"
    )
    
    # Save the numbers in the incl/excl flowchart
    df_incl_excl = pd.read_csv(
        INCL_EXCL_FLW_DIR / "df_inclusion_exclusion_longitudinal_tau_amyloid.csv", index_col=0
    )

    df_incl_excl[
        "N with amy and tau with two or more timepoints after quality check and same tracer and imaging and clinical data"
    ] = n_min
    df_incl_excl.to_csv(INCL_EXCL_FLW_DIR / "df_inclusion_exclusion_longitudinal_tau_amyloid.csv")
    logger.info(f"The df with less subjects is {which_df}, I will use this to filter the other two for the same subjects...")
    return which_df


def filter_save_dfs(which_df, df_amy_tau, df_adni, df_mri):
    if which_df == "mri":
        # Filter the dfs with more data to have the same subjects as the one with less data.
        df_amy_tau_final = df_amy_tau[df_amy_tau["PTID"].isin(df_mri["PTID"])].copy()
        df_adni_final = df_adni[
            df_adni["PTID"].isin(df_mri["PTID"])
            ].copy()
        df_mri_final = df_mri.copy()
        
    elif which_df == "adnimerge":
        # Filter the dfs with more data to have the same subjects as the one with less data.
        df_amy_tau_final = df_amy_tau[df_amy_tau["PTID"].isin(df_adni["PTID"])].copy()
        df_adni_final = df_adni.copy()
        df_mri_final = df_adni[
            df_adni["PTID"].isin(df_adni["PTID"])
            ].copy()
      
    # Save everything in the interim (there are stil some subjects that could get excluded in the preprocessing)
    df_amy_tau_final.to_csv(INTERIM_DATA_DIR / "df_amy_tau_after_cross_checking.csv", index=None)
    df_adni_final.to_csv(INTERIM_DATA_DIR / "df_adnimerge_after_cross_checking.csv", index=None)
    df_mri_final.to_csv(INTERIM_DATA_DIR / "df_mri_after_cross_checking.csv", index=None)
    return df_amy_tau_final, df_adni_final, df_mri_final

def clean_downloaded_folders(df_mri, input_path: Path = INTERIM_DATA_DIR / "ADNI"):
    """
    Cleans the folder containing downloaded images based on the df containing the sequences we want.
    """
    df_mri_temp = df_mri.copy()
    df_mri_temp["Description"] = df_mri_temp["Description"].str.replace(" ", "_")
    df_mri_temp["MRIDATE"] = df_mri_temp["MRIDATE"].astype(str)
    # Loop through the folder structure
    for subject_dir in input_path.iterdir():
        if subject_dir.is_dir():  # Check if it's a directory (subject)
            subject_name = subject_dir.name
            # If the subject is not in the df, remove the subject folder
            if not subject_name in df_mri["PTID"].values:
                shutil.rmtree(subject_dir)
                logger.info(f"Removed subject: {subject_dir}")
            # If the subject is included, check each sequence/date combination
            else:
                for sequence_dir in subject_dir.iterdir():
                    if sequence_dir.is_dir():  # Check if it's a directory (sequence)
                        sequence_name = sequence_dir.name
                        # If subjects were already converted their folder should start with M0s
                        if not sequence_name.startswith("M0"):
                            for date_dir in sequence_dir.iterdir():
                                if date_dir.is_dir():  # Check if it's a directory (date)
                                    date_name = date_dir.name[:-11]
                                    # Check if the combination of subject, date, and sequence exists in the DataFrame
                                    if not (
                                        (df_mri_temp["PTID"] == subject_name)
                                        & (df_mri_temp["MRIDATE"] == date_name)
                                        & (df_mri_temp["Description"] == sequence_name)
                                    ).any():
                                        # If the combination does not exist, remove the date folder
                                        shutil.rmtree(date_dir)
                                        logger.info(f"Removed sequence: {date_dir}")
                                    else:
                                        session_number = df_mri_temp[
                                            (df_mri_temp["PTID"] == subject_name)
                                            & (df_mri_temp["MRIDATE"] == date_name)
                                            & (df_mri_temp["Description"] == sequence_name)
                                        ]["session"].values[0]
                                        # If the combination exists, reorganize the folder structure
                                        new_path = (
                                            input_path
                                            / subject_name
                                            / f"{session_number}"
                                            / sequence_name
                                        )
                                        new_path.parent.mkdir(parents=True, exist_ok=True)
                                        shutil.move(str(date_dir), str(new_path))
                                        logger.info(f"Moved sequence: {date_dir} to {new_path}")
                            is_empty = not any(sequence_dir.iterdir())
                            if is_empty:
                                Path.rmdir(sequence_dir)
    # Make sure that the number of folder matches the unique subjects in the df
    assert df_mri_temp["PTID"].nunique() == len(
        [subj_dir for subj_dir in (INTERIM_DATA_DIR / "ADNI").iterdir() if subj_dir.is_dir()]
    ), "Something wrong happened when cleaning!"


def run_heudiconv(input_path: Path = INTERIM_DATA_DIR / "ADNI"):
    
    def add_adni_name(name):
        "Function to add ADNI to a string containing the subject name"
        if not "ADNI" in name:
            return f"{name.split('-')[0]}-ADNI{name.split('-')[1]}"
        else:
            return name

    logger.info("Running Heudiconv")
    # Loop through subject directories
    for subject_dir in input_path.iterdir():
        if subject_dir.is_dir():
            subj = subject_dir.name
            # Check for sessions between 00 and max_sessions
            list_ses = [ses_dir.name for ses_dir in subject_dir.iterdir()]
            for ses in list_ses:
                ses_dir = subject_dir / ses
                if ses_dir.is_dir():
                    # Fast enough not to be parallelized
                    # Construct the Singularity command
                    cmd = [
                        "singularity",
                        "run",
                        "-B",
                        f"{input_path}:/raw:ro",
                        "-B",
                        f"{BIDS_DATA_DIR}:/out",
                        "-B",
                        f"{SCRIPTS_DIR}:/scripts",
                        f"{SINGIMS_DIR}/heudiconv_1.0.1.sif",
                        "-d",
                        "/raw/{subject}/{session}/*/*/*.dcm",
                        "-s",
                        subj,
                        "--ses",
                        ses,
                        "-f",
                        "/scripts/heuristic.py",
                        "-c",
                        "dcm2niix",
                        "-b",
                        "--overwrite",
                        "-g",
                        "all",
                        "-o",
                        "/out",
                    ]
                    # Run the command
                    subprocess.run(cmd)
                else:
                    print(f"Session folder {ses_dir} not found for subject {subj}. Skipping...")

    # Check if all subjects were correctly converted into BIDS format and nifti and update the flowchart
    ls_names_bids = [
        subj_dir.name
        for subj_dir in BIDS_DATA_DIR.iterdir()
        if subj_dir.is_dir() and subj_dir.name.startswith("sub-")
    ]
    ls_names = [
        "sub-ADNI" + subj_dir.name.replace("_", "")
        for subj_dir in (INTERIM_DATA_DIR / "ADNI").iterdir()
        if subj_dir.is_dir()
    ]
    n_errors_heudiconv = len([sub for sub in ls_names if not sub in ls_names_bids])
    logger.warning(f"There has been {n_errors_heudiconv} subjects with errors in Heudiconv")
    df_incl_excl = pd.read_csv(
        INCL_EXCL_FLW_DIR / "df_inclusion_exclusion_longitudinal_tau_amyloid.csv", index_col=0
    )
    df_incl_excl["Errors in Heudiconv"] = n_errors_heudiconv
    df_incl_excl.to_csv(INCL_EXCL_FLW_DIR / "df_inclusion_exclusion_longitudinal_tau_amyloid.csv")
    # For some reasons some participants do not get included into the participants.tsv file automatically.
    # Here we check that, if the subject was converted into bids format by Heudiconv, then it should be also in
    # the participants.tsv file...
    df_participants_bids = pd.read_csv(BIDS_DATA_DIR / "participants.tsv", delimiter="\t")
    df_participants_bids["participant_id"] = df_participants_bids["participant_id"].apply(
        add_adni_name
    )
    not_in_participants_tsv = [
        sub for sub in ls_names_bids if sub not in df_participants_bids["participant_id"].values
    ]
    # We create fake age and sex since we are not going to use them from here anyways...
    df_to_add = pd.DataFrame(
        {
            "participant_id": not_in_participants_tsv,
            "age": [0 for _ in not_in_participants_tsv],
            "sex": [0 for _ in not_in_participants_tsv],
            "group": ["control" for _ in not_in_participants_tsv],
        }
    )

    df_new_participants_bids = pd.concat([df_participants_bids, df_to_add]).reset_index(drop=True)
    df_new_participants_bids.to_csv(BIDS_DATA_DIR / "participants.tsv", sep="\t", index=None)


def main():

    logger.info("Processing dataset...")
    logger.info("Downloading data...")
    download_dataset_adni()
    logger.success("Data downloaded correctly!")
    check_download_match()
    logger.info("Merging Amy/tau dfs with adnimerge and imaging.")

    df_adnimerge = load_df_adnimerge()
    df_mri_download = load_df_mri()
    df_amy_tau = load_df_amy_tau()

    df_adni_filtered = filter_adnimerge_based_on_pet(df_amy_tau, df_adnimerge)
    df_mri_filtered = filter_mri_based_on_pet(df_amy_tau, df_mri_download)
    df_mri_filtered_with_both_seqs = select_complete_mri(df_mri_filtered)
    name_of_df_min_subjs = perform_checks(df_amy_tau, df_adni_filtered, df_mri_filtered_with_both_seqs, )
    df_amy_tau_final, df_adni_final, df_mri_final = filter_save_dfs(name_of_df_min_subjs, df_amy_tau, df_adni_filtered, df_mri_filtered_with_both_seqs)
    logger.success("Dataframes processed correctly. Find them in data/interim.")
    logger.info("Cleaning the downloaded data folders")
    clean_downloaded_folders(df_mri_final)
    logger.success("Cleanup completed!")
    logger.info("Converting to BIDS format using Heudiconv...")
    run_heudiconv()
    logger.success("Data converted to BIDS correctly")
    logger.success("Dataset preprocessing complete!")

if __name__ == "__main__":
    is_sure = input("Are you sure you want to download the data and process it? [y/n]")
    if is_sure.lower() in ["y", "yes"]:
        main()
