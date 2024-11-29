# %%
import pandas as pd

from loguru import logger
from pathlib import Path
from cat_whim.config import (
    INCL_EXCL_FLW_DIR,
    RAW_DATA_DIR,
    INTERIM_DATA_DIR,
    UTILS_DATA_DIR,
    BIDS_DATA_DIR,
    PROCESSED_DATA_DIR,
)

LQT_DIR = BIDS_DATA_DIR / "derivatives" / "LQT"


# Functions
def format_df_aparc(df):

    new_col_names = "CTX_" + df.columns.str.upper()
    df.columns = new_col_names

    if "CTX_LH.APARC.THICKNESS" in df.columns:

        df = df.rename(columns={"CTX_LH.APARC.THICKNESS": "PTID"})

    elif "CTX_RH.APARC.THICKNESS" in df.columns:

        df = df.rename(columns={"CTX_RH.APARC.THICKNESS": "PTID"})

    df["PTID"] = df["PTID"].str.split(".long").str[0]
    df["session"] = df["PTID"].str.split("_ses-").str[1]
    df["PTID"] = df["PTID"].str.split("_ses-").str[0]

    return df


def format_df_adnimerge(df):
    adnimerge_cols_to_keep = [
        "PTID",
        "RID",
        "EXAMDATE",
        "AGE",
        "PTGENDER",
        "RAVLT_forgetting",
        "ADAS13",
        "DX_bl",
        "DX",
        "PTEDUCAT",
        "APOE4",
        "MMSE",
        "session",
        "Years_bl",
        "HMHYPERT",
    ]

    df = df[adnimerge_cols_to_keep].copy()
    df = df.rename(columns={"AGE": "AGE_bl"})
    df["PTID"] = df["PTID"].apply(lambda x: "sub-ADNI" + x.replace("_", ""))
    df["HMHYPERT"] = df["HMHYPERT"].astype(str)
    # NOTE: Adnimerge records the date as baseline in AGE (see: https://groups.google.com/g/adni-data/c/cBJW72gyJzg), while
    # Years_bl are the years from the baseline, so we can get the age of the subject at a particular session by summing them.
    df["AGE"] = df["AGE_bl"] + df["Years_bl"]

    return df


def format_df_disconn(df, dk_regions):

    df = df[["ID"] + dk_regions]
    df.columns = [col + "_disconn" for col in df.columns]
    df = df.rename(columns={"ID_disconn": "PTID"})
    df["session"] = df["PTID"].str.split("_ses-").str[1]
    df["PTID"] = df["PTID"].str.split("_ses-").str[0]
    df["PTID"] = df["PTID"].apply(lambda x: "sub-" + x)

    return df


def format_df_amy_tau(df, dk_regions):

    amy_labels = [reg + "_SUVR_amy" for reg in dk_regions]
    tau_labels = [reg + "_SUVR_tau" for reg in dk_regions]
    df["PTID"] = df["PTID"].apply(lambda x: "sub-ADNI" + x.replace("_", ""))
    filter_labels = ["PTID", "session", "TRACER_amy", "TRACER_tau", "AMYLOID_STATUS", "META_TEMPORAL_SUVR"] + amy_labels + tau_labels
    df = df[filter_labels].copy()

    return df

def check_session_m01_has_m00(df):

    # Filter rows for M00 and M01 sessions
    m00_sessions = df[df["session"] == "M00"].copy()
    m01_sessions = df[df["session"] == "M01"].copy()

    # Get the PTIDs for both M00 and M01 sessions
    m00_ptids = set(m00_sessions["PTID"])
    m01_ptids = set(m01_sessions["PTID"])

    # Find PTIDs that have M01 but don't have M00
    missing_m00 = m01_ptids - m00_ptids

    # Check if all M01 PTIDs also have M00
    if len(missing_m00) == 0:
        logger.success("All PTIDs with session M01 also have session M00.")
    else:
        logger.warning(f"PTIDs with session M01 but missing session M00: {missing_m00}")
        # Make the session M01 a session M00, since it is the same for this study
        for ptid in missing_m00:
            df.loc[df["PTID"] == ptid, "session"] = "M00"

    return df

def calculate_years_m00(group):

    # Check if session 'm00' exists
    if "M00" in group["session"].values:
        # Get the Years_bl for session 'm00'
        years_m00 = group.loc[group["session"] == "M00", "Years_bl"].values[0]
        # Subtract it from the Years_bl of all sessions
        group["Years_m00"] = group["Years_bl"] - years_m00
    else:
        # If 'm00' is not found, assign NaN or handle as needed
        group["Years_m00"] = pd.NA
    
    return group

def update_df_incl_excl(df_data, 
                        common_ptids,
                        ptids_amy_tau,
                        incl_excl_file : Path = INCL_EXCL_FLW_DIR / "df_inclusion_exclusion_longitudinal_tau_amyloid.csv"):

    df_to_calc = df_data[df_data["session"] == "M00"].copy()
    df_incl_excl = pd.read_csv(incl_excl_file, index_col=0)

    # Update the df_inclusion_exclusion using m00 not to duplicate the numbers
    df_incl_excl["N errors in freesurfer/smriprep/wmh segmentation"] = len(ptids_amy_tau) - len(
        common_ptids
    )
    df_incl_excl["N included before covariates"] = len(common_ptids)
    df_incl_excl["N without AGE"] = df_to_calc["AGE"].isna().sum()
    df_incl_excl["N without DX"] = df_to_calc["DX"].isna().sum()
    df_incl_excl["N without MMSE"] = df_to_calc["MMSE"].isna().sum()
    df_incl_excl["N without APOE4"] = df_to_calc["APOE4"].isna().sum()
    df_to_calc.dropna(subset=["MMSE", "DX", "AGE"], inplace=True)
    df_to_calc_complete_apoe = df_to_calc.dropna(subset=["APOE4"])
    df_incl_excl["N included analyses without APOE"] = df_to_calc["PTID"].nunique()
    df_incl_excl["N included analyses with APOE"] = df_to_calc_complete_apoe["PTID"].nunique()

    df_incl_excl.to_csv(incl_excl_file, index=None)

def create_df_long_format(df_m00_m01, df_m01, dk_regions, mean_or_region = "mean"):

    """These are mainly for plotting correlations at different timepoints"""

    if mean_or_region == "mean":

        thick_rate_year = (df_m01["mean_thickness"].values - df_m00_m01["mean_thickness"].values) / (df_m00_m01["mean_thickness"].values * df_m01["Years_m00"].values)
        amy_rate_year = (df_m01["mean_SUVR_amy"].values - df_m00_m01["mean_SUVR_amy"].values) / (df_m00_m01["mean_SUVR_amy"].values* df_m01["Years_m00"].values)
        tau_rate_year = (df_m01["mean_SUVR_tau"].values - df_m00_m01["mean_SUVR_tau"].values) / (df_m00_m01["mean_SUVR_tau"].values* df_m01["Years_m00"].values)
        disconn_rate_year = (df_m01["mean_disconn"].values - df_m00_m01["mean_disconn"].values) / (df_m00_m01["mean_disconn"].values* df_m01["Years_m00"].values)

        df_long = pd.DataFrame({
            "PTID": df_m00_m01["PTID"],
            "SUVR_amy": df_m00_m01["mean_SUVR_amy"].values,
            "SUVR_tau": df_m00_m01["mean_SUVR_tau"].values,
            "disconn": df_m00_m01["mean_disconn"].values,
            "thickness_m00": df_m00_m01["mean_thickness"].values,
            "thickness_m01": df_m01["mean_thickness"].values,
            "thick_rate_year": thick_rate_year,
            "amy_rate_year": amy_rate_year,
            "tau_rate_year": tau_rate_year,
            "disconn_rate_year": disconn_rate_year,
            })

    elif mean_or_region == "region":

        df_long = pd.DataFrame()

        for region_label in dk_regions:

            region_label_amy = region_label + "_SUVR_amy"
            region_label_tau = region_label + "_SUVR_tau"
            region_label_disconn = region_label + "_disconn"
            region_label_thick = region_label + "_THICKNESS"

            thick_rate_year = (df_m01[region_label_thick].values - df_m00_m01[region_label_thick].values) / (df_m00_m01[region_label_thick].values * df_m01["Years_m00"].values)

            df_region = pd.DataFrame({
                "PTID": df_m00_m01["PTID"],
                "region": [region_label for _ in range(df_m00_m01.shape[0])],
                "SUVR_amy": df_m00_m01[region_label_amy].values,
                "SUVR_tau": df_m00_m01[region_label_tau].values,
                "disconn": df_m00_m01[region_label_disconn].values,
                "thickness_m00": df_m00_m01[region_label_thick].values,
                "thickness_m01": df_m01[region_label_thick].values,
                "thick_rate_year": thick_rate_year,
                "Years_m00": df_m01["Years_m00"].values
                })

            df_long = pd.concat([df_long, df_region]).reset_index(drop=True)

    df_long = pd.merge(df_long, df_m00_m01[["PTID", "DX"]], on = ["PTID"])

    return df_long


def exclude_infarcts(df):
    
    df_infarcts = pd.read_csv(UTILS_DATA_DIR / "MRI_INFARCTS_01_29_21_30Oct2024.csv")
    df_inf = df_infarcts[df_infarcts["RID"].isin(df["RID"])]
    df_inf_final = df_inf[df_inf["SIDE"] != "-"].sort_values(by=["RID", "MRI.DATE1"]).drop_duplicates(subset=["RID"], keep="first")
    infarcts_to_remove = []
    for rid in df_inf_final["RID"]:
        date_df = df[df["RID"] == rid]["EXAMDATE"].values[0]
        date_inf = df_inf_final[df_inf_final["RID"] == rid]["MRI.DATE1"].values[0]
        if date_df > date_inf:
            infarcts_to_remove.append(rid)
        else:
            continue
    df = df[~df["RID"].isin(infarcts_to_remove)].copy()
    n_infarcts = len(infarcts_to_remove)
    return df, n_infarcts

def exclude_ad(df):
    n_dementia = df[df["DX"] == "Dementia"]["PTID"].nunique()
    df = df[df["DX"]!= "Dementia"].reset_index(drop=True).copy()
    return df, n_dementia


def main():
    # Load dataframes
    df_adnimerge = pd.read_csv(INTERIM_DATA_DIR / "df_adnimerge_after_cross_checking.csv")
    df_adnimerge = df_adnimerge.drop_duplicates(subset=["PTID", "session"])
    df_htn = pd.read_csv(RAW_DATA_DIR / "MODHACH_09Sep2024.csv")
    df_adnimerge = pd.merge(df_adnimerge, df_htn[["RID", "HMHYPERT"]], on="RID")
    df_amy_tau = pd.read_csv(INTERIM_DATA_DIR / "df_amy_tau_after_cross_checking.csv")
    df_aparc_lh = pd.read_csv(INTERIM_DATA_DIR / "aparc_lh_table.tsv", sep="\t")
    df_aparc_rh = pd.read_csv(INTERIM_DATA_DIR / "aparc_rh_table.tsv", sep="\t")
    df_disconn = pd.read_csv(LQT_DIR / "results" / "dataframes" / "parc_discon.csv", index_col=0)

    # Format dataframes
    df_adnimerge = format_df_adnimerge(df_adnimerge)
    df_aparc_lh = format_df_aparc(df_aparc_lh)
    df_aparc_rh = format_df_aparc(df_aparc_rh)

    # Merge the aparc dfs
    df_aparc_lh_rh = pd.merge(
        df_aparc_lh.drop(columns=["CTX_BRAINSEGVOLNOTVENT", "CTX_ETIV"]),
        df_aparc_rh,
        on=["PTID", "session"],
    )  # Drop CTX_BRAINSEGVOLNOTVENT and CTX_ETIV as these are the same also in df_rh_aparc

    # Get the region labels for the atlas
    dk_regions_lh_rh = [reg.split("_THICKNESS")[0] for reg in df_aparc_lh.columns[1:-4]] + [
        reg.split("_THICKNESS")[0] for reg in df_aparc_rh.columns[1:-4]
    ]

    pd.DataFrame(dk_regions_lh_rh).to_csv(PROCESSED_DATA_DIR / "dk_region_names.csv", index=None)
    # Format remaining dataframes based on names of regions
    df_disconn = format_df_disconn(df_disconn, dk_regions_lh_rh)
    df_amy_tau = format_df_amy_tau(df_amy_tau, dk_regions_lh_rh)

    # Find PTIDs from each DataFrame
    ptids_amy_tau = set(df_amy_tau["PTID"])
    ptids_disconn = set(df_disconn["PTID"])
    ptids_adnimerge = set(df_adnimerge["PTID"])
    ptids_aparc_lh_rh = set(df_aparc_lh_rh["PTID"])

    # Find PTIDs that are in all four DataFrames
    common_ptids = ptids_amy_tau & ptids_disconn & ptids_adnimerge & ptids_aparc_lh_rh

    # Filter for subjects that are in all DataFrames
    df_amy_tau_final = df_amy_tau[df_amy_tau["PTID"].isin(common_ptids)].copy()
    df_disconn_final = df_disconn[df_disconn["PTID"].isin(common_ptids)].copy()
    df_aparc_lh_rh_final = df_aparc_lh_rh[df_aparc_lh_rh["PTID"].isin(common_ptids)].copy()
    df_adnimerge_final = df_adnimerge[df_adnimerge["PTID"].isin(common_ptids)].copy()

    # Create intermediate dataset with all imaging data
    df_amy_tau_disconn = pd.merge(df_amy_tau_final, df_disconn_final, on=["PTID", "session"])
    df_imaging = pd.merge(df_amy_tau_disconn, df_aparc_lh_rh_final, on=["PTID", "session"])

    # Create final dataset with imaging and clinical data
    df_adnimerge_imaging = pd.merge(df_imaging, df_adnimerge_final, on=["PTID", "session"])

    # Check that who has ses-M01 also has ses-M00
    df_adnimerge_imaging = check_session_m01_has_m00(df_adnimerge_imaging)

    # Calculate the years from M00
    df_adnimerge_imaging = (
        df_adnimerge_imaging.groupby("PTID", group_keys=True)
        .apply(calculate_years_m00)
        .drop(columns="Years_bl")
    )
    df_adnimerge_imaging.reset_index(drop=True, inplace=True)

    df_adnimerge_imaging, n_infarcts = exclude_infarcts(df_adnimerge_imaging)
    df_adnimerge_imaging, n_dementia = exclude_ad(df_adnimerge_imaging)

    df_incl_excl = pd.read_csv(INCL_EXCL_FLW_DIR / "df_inclusion_exclusion_longitudinal_tau_amyloid.csv")
    df_incl_excl["Infarcts"] = n_infarcts
    df_incl_excl["Dementia"] = n_dementia
    df_incl_excl.to_csv(INCL_EXCL_FLW_DIR / "df_inclusion_exclusion_longitudinal_tau_amyloid.csv", index=None)

    # Update the inclusion/exclusion DataFrame
    update_df_incl_excl(df_adnimerge_imaging, common_ptids, ptids_amy_tau)
    df_adnimerge_imaging.dropna(subset=["MMSE", "AGE", "DX", "PTEDUCAT"], inplace=True)
    
    # Save what we need
    df_adnimerge_imaging.to_csv(PROCESSED_DATA_DIR / "df_cat_whim_to_analyze.csv", index=None)
    logger.success("Datasets created correctly! You can start analyzing!")

if __name__ == "__main__":
    main()
# %%
