######################################################################################
# This script selects subjects with longitudinal tau/amyloid data from ADNI that was
# already preprocessed and saves a .csv files with the names in data/utils. By copying 
# this file into ADNI/Download/Image Collections/Advanced Search/Subject ID.
# REMEMBER that you can't directly compare different tracers and we have 2 different
# tracers for amyloid, so you need to divide into two datasets!
# TODO: Remove the 17Jul2024 before releasing so that people can use different dates
# TODO: Currently we have only one TAU pet tracer but this could not be the case for
# other datasets...
######################################################################################
# %%
import pandas as pd
import numpy as np

from loguru import logger
from pandas.tseries.offsets import MonthEnd

from cat_whim.config import RAW_DATA_DIR, INTERIM_DATA_DIR, INCL_EXCL_FLW_DIR


def load_df_berkeley(date_amy, date_tau):
    """
    This function loads the df amy and tau (partial volume corrected) from UC Berkeley, 
    it selects scans with good quality and with no nans.
    """
    
    # Load already preprocessed data from UC Berkeley
    data_amy = RAW_DATA_DIR / f"UCBERKELEY_AMY_6MM_{date_amy}.csv"
    data_tau = RAW_DATA_DIR / f"UCBERKELEY_TAUPVC_6MM_{date_tau}.csv"

    df_amy = pd.read_csv(data_amy)
    df_tau = pd.read_csv(data_tau)

    df_amy["SCANDATE"] = pd.to_datetime(df_amy["SCANDATE"], errors="coerce")
    n_df_amy_init = df_amy["PTID"].unique().shape[0]

    df_tau["SCANDATE"] = pd.to_datetime(df_tau["SCANDATE"], errors="coerce")
    n_df_tau_init = df_amy["PTID"].unique().shape[0]

    return n_df_amy_init, n_df_tau_init, df_amy, df_tau


def filter_dfs_same_timepoint(df_amy, df_tau):
    """
    This function filters for subject that have tau and amyloid data within 6 months from one another
    """
    # Merge on PTID, because dates of amyloid and tau PET can differ
    df_amy_tau = pd.merge(df_amy, df_tau, on=["PTID", "VISCODE"], suffixes=("_amy", "_tau"))
    # Filter based on the date condition of 6 months interval between the two scans at most
    df_amy_tau_filtered = df_amy_tau[
        (df_amy_tau["SCANDATE_tau"] - MonthEnd(3) <= df_amy_tau["SCANDATE_amy"])
        & (df_amy_tau["SCANDATE_tau"] + MonthEnd(3) >= df_amy_tau["SCANDATE_amy"])
    ].copy()

    # Also sort by date
    df_amy_tau_filtered = df_amy_tau_filtered.sort_values(by=["PTID", "SCANDATE_amy"])
    n_amy_tau_same_tp = df_amy_tau_filtered["PTID"].unique().shape[0]

    return n_amy_tau_same_tp, df_amy_tau_filtered


def get_longitudinal_number(df_amy_tau_filtered):

    # Assign ranks to each SCANDATE_amy within each subject (PTID)
    df_amy_tau_filtered["session"] = df_amy_tau_filtered.groupby("PTID").cumcount()
    df_amy_tau_filtered["session"] = df_amy_tau_filtered["session"].apply(lambda x : "M0" + str(x))

    return df_amy_tau_filtered


def get_dk_labels_from_df_berkeley(df):
    """
    Defines a list of columns from the DK atlas starting from the UC Berkeley csv files
    """

    list_cols_dk_pre = df.columns[
        (df.columns.str.contains("SUVR"))
        & ((df.columns.str.contains("LH")) | (df.columns.str.contains("RH")))
    ].to_list()

    # "CTX_ENTORHINAL_SUVR" gets picked out even though we don't want it, so we remove it
    list_dk_labels = [col for col in list_cols_dk_pre if col != "CTX_ENTORHINAL_SUVR"]

    return list_dk_labels


def qc_berkeley(df):
    """
    Performs QC for Amyloid and Tau PET data and drops subjects with nans in DK atlas SUVR regions
    """
    # Select cortical labels of DK atlas
    list_dk_labels = get_dk_labels_from_df_berkeley(df)
    # Drop bad scans: Full preprocessing quality control is good if 2 (partial if 1, failed if 0)
    # Note that when using TAUPVC df, there is no qc_flag possibly because they processed only
    # good subjects, if you use the "simple" Tau df (not pvc), then change as follows...
    # ... df_good = df[(df["qc_flag"] > 1) & (df["qc_flag_tau"] > 1)]
    df_good = df[df["qc_flag"] > 1]
    # Drop subjects who have nas in the cols of interest
    df_good_no_na = df_good.dropna(subset=list_dk_labels)
    n_good_no_na = df_good_no_na["PTID"].unique().shape[0]
    return n_good_no_na, df_good_no_na

def check_same_tracer(df):
    """
    Check that the tracer for amyloid or tau is the same. We can't compare different amyloid/tau tracers to one another.
    """
    list_tracers = ["TRACER_amy", "TRACER_tau"]
    list_subjs_different_tracer = []
    for tracer in list_tracers:
        # Compute nunique values of tracer for each PTID
        nunique_counts = df.groupby("PTID")[tracer].nunique()
        # Filter out subjects where nunique is not equal to 1
        subjs_different_tracer = [sub for sub in nunique_counts[nunique_counts != 1].index]
        if len(subjs_different_tracer) >= 1:
            for subj in subjs_different_tracer:
                list_subjs_different_tracer.append(subj)
    df_amy_tau_final = df[~df["PTID"].isin(list_subjs_different_tracer)].copy()
    n_subjs_same_tracer = df_amy_tau_final["PTID"].unique().shape[0]
    return n_subjs_same_tracer, df_amy_tau_final


def save_df_incl_excl(
    n_amy_init,
    n_tau_init,
    n_amy_tau_same_tp,
    n_passed_qc,
    n_same_tracer,
    n_amyFBP_final,
    n_amyFBB_final,
):
    """
    This function saves the numbers of subjects that were at the beginning of each step, so that we can recreate an inclusion/exclusion flowchart.
    """


    df_incl_excl = pd.DataFrame(
        {
            "N initial amyloid": [n_amy_init],
            "N initial tau": [n_tau_init],
            "N initial with both amy and tau same timepoint": [n_amy_tau_same_tp],
            "N with amy and tau same timepoint after quality check": [n_passed_qc],
            "N with amy and tau same timepoint after quality check and same tracer in both timepoints": [
                n_same_tracer
            ],
            "N amy and tau with two or more timepoints after quality check and same tracer (FBP)": n_amyFBP_final,
            "N amy and tau with two or more timepoints after quality check and same tracer (FBB)": n_amyFBB_final,
        }
    )

    df_incl_excl.to_csv(INCL_EXCL_FLW_DIR / "df_inclusion_exclusion_longitudinal_tau_amyloid.csv", index=None)
    return df_incl_excl

def main():

    logger.info(
        "Starting to select subjects with available amyloid/tau PET data at the same timepoint"
    )
    # TODO: uncomment before releasing
    # date_amy = input(
    #     "Add AMY PET date in this format: e.g., 17Jul2024 (the date we used)"
    # )  # TODO: Should assert the regex expression for the date...
    # date_tau = input(
    #     "Add TAU PET (partial volume corrected (PVC)) date in this format: e.g., 17Jul2024 (the date we used)"
    # )
    date_amy = "17Jul2024"
    date_tau = "17Jul2024"
    # 1. Load the df with preprocessed amyloid / tau data
    n_amy_init, n_tau_init, df_amy, df_tau = load_df_berkeley(date_amy, date_tau)
    logger.info(f"UCBERKELEY_AMY_6MM_{date_amy}.csv contains {n_amy_init} unique subjects")
    logger.info(f"UCBERKELEY_TAUPVC_6MM_{date_amy}.csv contains {n_tau_init} unique subjects")

    # 2. Make sure that amyloid and tau scans were acquired within 6 months from one another
    n_amy_tau_same_tp, df_amy_tau_filtered = filter_dfs_same_timepoint(df_amy, df_tau)

    # 3. Get the session number 
    df_amy_tau_with_fu = get_longitudinal_number(df_amy_tau_filtered)
    logger.info(f"There are {n_amy_tau_same_tp} subjects with amyloid and tau at the same timepoint")

    # 4. Perform quality check for amy/tau pet scan and check nans
    n_passed_qc, df_amy_tau_with_fu_qc = qc_berkeley(df_amy_tau_with_fu)
    logger.info(
        f"There are {n_passed_qc} subjects with amyloid and tau at the same timepoint for 2 timepoints after quality check"
    )

    # 5. Check that the tracer for the two timepoints for the same subject is the same, otherwise we can't compare
    n_same_tracer, df_amy_tau_with_fu_qc_same_tracer = check_same_tracer(
        df_amy_tau_with_fu_qc
    )
    n_different_tracer = n_passed_qc - n_same_tracer
    logger.warning(f"{n_different_tracer} subjects have different tracers at the baseline and one of the follow-ups. These are not comparable, so we drop these subjects.")
    logger.info(f"{n_same_tracer} subjects passed all preprocessing and have the same tracer in all studies.")

    # 6. Separate the dataframes based on amy tracer (tau is the same in our case)
    df_amyFBP_tau_final = df_amy_tau_with_fu_qc_same_tracer[
        df_amy_tau_with_fu_qc_same_tracer["TRACER_amy"] == "FBP"
    ]
    df_amyFBB_tau_final = df_amy_tau_with_fu_qc_same_tracer[
        df_amy_tau_with_fu_qc_same_tracer["TRACER_amy"] == "FBB"
    ]

    n_amyFBP_final = df_amyFBP_tau_final["PTID"].unique().shape[0]
    n_amyFBB_final = df_amyFBB_tau_final["PTID"].unique().shape[0]

    df_amy_tau_with_fu_qc_same_tracer.to_csv(INTERIM_DATA_DIR / "df_amyBoth_tau_to_download.csv", index=None)
    df_amyFBP_tau_final.to_csv(INTERIM_DATA_DIR / "df_amyFBP_tau_to_download.csv", index=None)
    df_amyFBB_tau_final.to_csv(INTERIM_DATA_DIR / "df_amyFBB_tau_to_download.csv", index=None)
    logger.success("Created complete and separate dfs for subjects who underwent both scans with the same tracer (FBP-PET or FBB-PET)")
    logger.info(f"The FBP-PET df contains {n_amyFBP_final} subects...")
    logger.info(f"The FBB-PET df contains {n_amyFBB_final} subects...")

    # 7. Get the list of subjects to download
    subjs_to_download = df_amy_tau_with_fu_qc_same_tracer["PTID"].unique()
    np.savetxt(str(INCL_EXCL_FLW_DIR / "subjects_to_download.csv"), subjs_to_download, delimiter=",", fmt="%s")
    logger.success(
        f"Created list of subjects to download, find it at {INCL_EXCL_FLW_DIR}/subjects_to_download.csv"
    )

    # 8. Save the number of subjects we included for each step
    save_df_incl_excl(
        n_amy_init,
        n_tau_init,
        n_amy_tau_same_tp,
        n_passed_qc,
        n_same_tracer,
        n_amyFBP_final,
        n_amyFBB_final,
    )
    logger.success(
        f"Created inclusion and exclusion dataframe, find it at {INCL_EXCL_FLW_DIR}/df_inclusion_exclusion_longitudinal_tau_amyloid.csv"
    )

    for subj in subjs_to_download:
        print(subj, ",")
    logger.info("Copy/Paste the above list separated with commas for ease of use.")
    logger.success("UC Berkeley preprocessing completed succesfully! Now you can download the subjects' MRI data from the ADNI website!")

if __name__ == "__main__":
    main()
