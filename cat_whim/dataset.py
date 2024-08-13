#%%
import typer
import csv
import subprocess
import shutil
import pandas as pd
import numpy as np

from pathlib import Path
from loguru import logger
from tqdm import tqdm

from cat_whim.config import RAW_DATA_DIR, INTERIM_DATA_DIR, BIDS_DATA_DIR, SCRIPTS_DIR, SINGIMS_DIR, INCL_EXCL_FLW_DIR

app = typer.Typer()


def download_dataset_adni(input_path: Path = RAW_DATA_DIR / "cat_whim_url_list_2024-08-06.csv"):
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
        subprocess.run(['wget', '-P', str(RAW_DATA_DIR), url])

        # Unzip the downloaded file if it's a zip file
        if file_name.endswith(".zip"):
            print(f"unzipping {file_name}")
            subprocess.run(['unzip', '-d', str(RAW_DATA_DIR), str(RAW_DATA_DIR / file_name)])
            logger.info(f"Extracted {file_name} successfully!")
            print(f"removing {file_name} folder")
            Path.unlink(RAW_DATA_DIR / file_name)
    logger.info("Copying downloaded folder to data/interim")
    cmd = ["cp", "-r", f"{RAW_DATA_DIR}/ADNI", f"{INTERIM_DATA_DIR}"]
    # Run the command
    subprocess.run(cmd)



def check_download_match(
    input_path: Path = RAW_DATA_DIR / "ADNI",
    input_file: Path = INCL_EXCL_FLW_DIR / "subjects_to_download.csv",
):

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
        logger.warning(f"Something went wrong! You are missing {len(list_missing)} subjects... {list_missing}. You should check what happened! (e.g., maybe the subject doesn't have T1 or FLAIR")
        df_incl_excl = pd.read_csv(INCL_EXCL_FLW_DIR / "df_inclusion_exclusion_longitudinal_tau_amyloid.csv", index_col=0)
        df_incl_excl["excluded_after_download"] = len(list_missing)
        df_incl_excl.to_csv(INCL_EXCL_FLW_DIR / "df_inclusion_exclusion_longitudinal_tau_amyloid.csv")

def load_df_amy_tau():
    return pd.read_csv(INTERIM_DATA_DIR / "df_amyFBP_tau_final.csv")


def load_df_adnimerge():
    return pd.read_csv(RAW_DATA_DIR / "ADNIMERGE_18Jul2024.csv", dtype="object")


def load_df_mri():

    df_images = pd.read_csv(RAW_DATA_DIR / "cat_whim_8_06_2024.csv")
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
    # Add the visit column that we want to keep for reference between datasets
    adni_cols = adni_cols + ["visit"]
    # Merge PET and ADNIMERGE
    merged_df = pd.merge(df_pet, df_adni, on="PTID", suffixes=("_pet", "_adnimerge"))
    # Calculate the difference between the dates in days
    merged_df["date_diff"] = np.abs((merged_df["SCANDATE_amy"] - merged_df["EXAMDATE"]).dt.days)
    # We consider as the same timepoints visits that were performed less than 6 months apart
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
    # Add the visit column that we want to keep for reference between datasets
    mri_cols = mri_cols + ["visit"]
    # Merge PET and ADNIMERGE
    merged_df = pd.merge(df_pet, df_mri, on="PTID")
    # Calculate the difference between the dates in days
    merged_df["date_diff"] = np.abs((merged_df["SCANDATE_amy"] - merged_df["MRIDATE"]).dt.days)
    # We consider as the same timepoint, visits that were performed less than 6 months apart
    result_df = merged_df[merged_df["date_diff"] <= 180]
    # Recreate the adnimerge only for data we want to include
    df_mri_new = result_df[mri_cols]
    return df_mri_new


def select_complete_mri(df_filtered):
    """
    Preprocess dataframe of downloaded images to keep only subjects with all MRI data needed T1 and FLAIR.
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

    df_counts = df_filtered.groupby(["PTID", "MRIDATE"]).count()
    less_than_two = df_counts[df_counts["Description"] <= 1].index
    df_filtered = df_filtered[
        ~df_filtered.set_index(["PTID", "MRIDATE"]).index.isin(less_than_two)
    ]

    result = df_filtered.groupby(["PTID", "MRIDATE"])["Description"].unique()
    check_results = result.apply(check_description)
    df_check_results = check_results.reset_index(name="Contains_T1_and_Flair")
    df_mri_complete = df_check_results[df_check_results["Contains_T1_and_Flair"] == True]
    df_mri_complete_visit = (
        pd.merge(
            df_mri_complete,
            df_filtered[["PTID", "MRIDATE", "visit", "Description"]],
            on=["PTID", "MRIDATE"],
        )
        .drop_duplicates()
        .reset_index(drop=True)
    )

    return df_mri_complete_visit


def perform_checks(df_mri_filtered_complete, df_adni_filtered, df_amy_tau):

    # Check that all subjects in df_adni_filtered and df_mri_filtered_complete have both visit 0 and 1
    # Find where the timepoints are still 2 for a subject.
    idx_to_keep_mri = np.argwhere(
        df_mri_filtered_complete.groupby("PTID")["visit"].unique().apply(lambda x: len(x) == 2)
    ).squeeze()
    # Find the names of these subjects
    ptid_to_keep_mri = (
        df_mri_filtered_complete.groupby("PTID")["visit"]
        .unique()
        .iloc[idx_to_keep_mri]
        .index.to_list()
    )
    # Keep only these subjects
    df_mri_filtered_complete_to_save = df_mri_filtered_complete[
        df_mri_filtered_complete["PTID"].isin(ptid_to_keep_mri)
    ].copy()

    # Do the same for adnimerge
    idx_to_keep_adni = np.argwhere(
        df_adni_filtered.groupby("PTID")["visit"].unique().apply(lambda x: len(x) == 2)
    ).squeeze()
    ptid_to_keep_adni = (
        df_adni_filtered.groupby("PTID")["visit"].unique().iloc[idx_to_keep_adni].index.to_list()
    ).copy()
    df_adni_complete = df_adni_filtered[df_adni_filtered["PTID"].isin(ptid_to_keep_adni)]

    # Check which df has the minimum number of subjects, then save the number of subjects that get excluded
    # due to the fact that they don't have clinical or imaging at the same timepoint as PET.
    logger.info(
        f"The number of subjects with complete MRI is: {df_mri_filtered_complete_to_save['PTID'].unique().shape[0]}"
    )
    logger.info(
        f"That of subjects with complete adnimerge is: {df_adni_complete['PTID'].unique().shape[0]}"
    )
    logger.info(
        f"That of subjects initially with amyloid and tau: {df_amy_tau['PTID'].unique().shape[0]}"
    )
    # Save the numbers in the incl/excl flowchart
    df_incl_excl = pd.read_csv(
        INCL_EXCL_FLW_DIR / "df_inclusion_exclusion_longitudinal_tau_amyloid.csv", index_col=0
    )
    df_incl_excl[
        "N with amy and tau with two or more timepoints after quality check and same tracer (FBP) and imaging and clinical data"
    ] = (df_mri_filtered_complete_to_save["PTID"].unique().shape[0])
    df_incl_excl.to_csv(
        INCL_EXCL_FLW_DIR / "df_inclusion_exclusion_longitudinal_tau_amyloid.csv", index=None
    )
    logger.info("Filtering amy tau and adnimerge dfs for subjects with available imaging")

    return df_mri_filtered_complete_to_save, df_adni_complete


def filter_save_dfs(df_mri_final, df_adni_filtered_cross, df_amy_tau):
    # Filter the dfs with more data to have the same subjects as the one with less data.
    df_amy_tau_final = df_amy_tau[df_amy_tau["PTID"].isin(df_mri_final["PTID"])].copy()
    df_adni_final = df_adni_filtered_cross[
        df_adni_filtered_cross["PTID"].isin(df_mri_final["PTID"])
    ].copy()

    # Save everything in the interim (there are stil some subjects that could get excluded in the preprocessing)
    df_mri_final.to_csv(INTERIM_DATA_DIR / "df_mri_after_cross_checking.csv", index=None)
    df_amy_tau_final.to_csv(INTERIM_DATA_DIR / "df_amy_tau_after_cross_checking.csv", index=None)
    df_adni_final.to_csv(INTERIM_DATA_DIR / "df_adnimerge_after_cross_checking.csv", index=None)


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
                        # If subjects were already converted their folder should start with M0 
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
                                        visit_number = df_mri_temp[
                                            (df_mri_temp["PTID"] == subject_name)
                                            & (df_mri_temp["MRIDATE"] == date_name)
                                            & (df_mri_temp["Description"] == sequence_name)]["visit"].values[0]
                                        # If the combination exists, reorganize the folder structure
                                        new_path = (
                                            input_path / subject_name / f"M0{visit_number}" / sequence_name
                                        )
                                        new_path.parent.mkdir(parents=True, exist_ok=True)
                                        shutil.move(str(date_dir), str(new_path))
                                        logger.info(f"Moved sequence: {date_dir} to {new_path}")
                            is_empty = not any(sequence_dir.iterdir())
                            if is_empty:
                                Path.rmdir(sequence_dir)


def run_heudiconv(input_path : Path = INTERIM_DATA_DIR / "ADNI"):
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
                        "singularity", "run",
                        "-B", f"{input_path}:/raw:ro",
                        "-B", f"{BIDS_DATA_DIR}:/out",
                        "-B", f"{SCRIPTS_DIR}:/scripts",
                        f"{SINGIMS_DIR}/heudiconv_1.0.1.sif",
                        "-d", "/raw/{subject}/{session}/*/*/*.dcm",
                        "-s", subj,
                        "--ses", ses,
                        "-f", "/scripts/heuristic.py",
                        "-c", "dcm2niix", "-b",
                        "--overwrite",
                        "-g", "all",
                        "-o", "/out"
                    ]
                    # Run the command
                    subprocess.run(cmd)
                else:
                    print(f"Session folder {ses_dir} not found for subject {subj}. Skipping...")


@app.command()
def main():

    logger.info("Processing dataset...")
    logger.info("Downloading data...")
    # download_dataset_adni()
    # check_download_match()
    # logger.success("Data downloaded correctly!")
    
    # logger.info("Merging Amy (FBP)/tau dfs with adnimerge and imaging.")

    # df_adnimerge = load_df_adnimerge()
    # df_mri_download = load_df_mri()
    # df_amy_tau = load_df_amy_tau()

    # df_adni_filtered = filter_adnimerge_based_on_pet(df_amy_tau, df_adnimerge)
    # df_mri_filtered = filter_mri_based_on_pet(df_amy_tau, df_mri_download)
    # df_mri_filtered_with_both_seqs = select_complete_mri(df_mri_filtered)

    # df_mri_final, df_adni_filtered_cross = perform_checks(
    #     df_mri_filtered_with_both_seqs, df_adni_filtered, df_amy_tau
    # )
    # filter_save_dfs(df_mri_final, df_adni_filtered_cross, df_amy_tau)

    # logger.success("Dataframes processed correctly. Find them in data/interim.")
    # logger.info("Cleaning the downloaded data folders")
    # clean_downloaded_folders(df_mri_final)
    # logger.success("Cleanup completed!")
    # logger.info("Converting to BIDS format using Heudiconv...")
    # run_heudiconv()
    # logger.success("Data converted to BIDS correctly")
    logger.success("Dataset preprocessing complete!")

if __name__ == "__main__":
    app()
