
import numpy as np
import pandas as pd
import statsmodels.api as sm
import seaborn as sns
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

from scipy import stats
from matplotlib.patches import Patch
from scipy.stats import zscore, false_discovery_control


def make_demographics_table(df):
    df_m00 = df[df["session"] == "M00"].copy()
    # Continuous variables
    continuous_vars = ["AGE", "PTEDUCAT", "MMSE", "wmh_vol_log", "RAVLT_forgetting"]
    normality_results = {}
    desc_stats = {}

    for var in continuous_vars:
        # Perform Shapiro-Wilk test for normality
        stat, p_value = stats.shapiro(df_m00[var].dropna())  # Drop NA values for the test
        normality_results[var] = p_value

        # Calculate mean and SD or median and IQR based on normality
        if p_value > 0.05:  # Data is normal
            desc_stats[var] = {
                "Variable": var,
                "Statistic": f"{df_m00[var].mean():.2f} ± {df_m00[var].std():.2f}",
                "Note": "Mean (SD)",
            }
            print(var, "mean (SD)")
        else:  # Data is not normal
            desc_stats[var] = {
                "Variable": var,
                "Statistic": f"{df_m00[var].median():.2f} ({df_m00[var].quantile(0.25):.2f} - {df_m00[var].quantile(0.75):.2f})",
                "Note": "Median (IQR)",
            }

    # Categorical variables
    categorical_vars = ["PTGENDER", "DX", "AMYLOID_STATUS"]
    cat_stats = []

    for var in categorical_vars:
        # Add the variable name as the first row
        cat_stats.append({"Variable": var, "Statistic": "", "Note": "Frequency and Proportion"})

        # Get frequency and proportions
        freq = df_m00[var].value_counts()
        proportions = df_m00[var].value_counts(normalize=True)

        # Add each group as a separate row
        for index in freq.index:
            cat_stats.append(
                {
                    "Variable": f"{var} - {index}",
                    "Statistic": f"{freq[index]} ({proportions[index] * 100:.1f}%)",
                    "Note": "",
                }
            )

    # Combine continuous and categorical tables
    table = pd.DataFrame.from_dict(desc_stats, orient="index")
    cat_table = pd.DataFrame(cat_stats)

    df_fu = df[df['Years_m00'] != 0].copy()
    demo_fu = pd.DataFrame({"Variable": [f"Follow-Up, years (n={df_fu['PTID'].nunique()})"],
                "Statistic": [f"{df_fu['Years_m00'].median():.2f} ({df_fu['Years_m00'].quantile(0.25):.2f} - {df_fu['Years_m00'].quantile(0.75):.2f})"],
                })

    # Combine both tables into one final table
    final_table = pd.concat([table, cat_table, demo_fu], axis=0, ignore_index=True)
    num_subj = df_m00["PTID"].nunique()
    demo_table = final_table[["Variable", "Statistic"]].copy()
    demo_table = demo_table.rename(columns={"Statistic": f"All (n = {num_subj})"})
    return demo_table


def plot_regplot_with_corr(df, x, y, ax=None):
    """
    Plots a regression plot with Seaborn and adds the correlation coefficient
    (Pearson or Spearman) to the title of the plot.

    Parameters:
    - x: array-like or pandas Series
    - y: array-like or pandas Series
    - ax: matplotlib axis object, optional

    Returns:
    - ax: The matplotlib axis with the plot
    """
    # Create a matplotlib axis if not provided
    if ax is None:
        fig, ax = plt.subplots()

    # Test normality of x and y using Shapiro-Wilk test
    _, p_x = stats.shapiro(df[x])
    _, p_y = stats.shapiro(df[y])

    # Choose correlation method: Spearman for non-normal data, Pearson for normal data
    if p_x < 0.05 or p_y < 0.05:
        corr_method = "Spearman"
        corr_coef, p = stats.spearmanr(df[x], df[y])
    else:
        corr_method = "Pearson"
        corr_coef, p = stats.pearsonr(df[x], df[y])

    # Plot regression plot using Seaborn
    sns.regplot(data=df, x=x, y=y, ax=ax, line_kws={"color":"red"}, scatter_kws={"color":"black", "alpha": 0.4, "s":20})

    # Add the correlation coefficient to the title
    ax.set_title(f"{corr_method} r = {corr_coef:.2f}, p = {p:.2f}", fontsize=8)

    return ax

def plot_boxplot_with_comparison(df, ax, x, y, pairs):
    from statannotations.Annotator import Annotator

    # Plot each boxplot and add statistical annotations
    sns.boxplot(data=df, x=x, y=y, ax=ax, palette="tab10", hue=x)
    annotator = Annotator(ax, pairs, data=df, x=x, y=y)
    annotator.configure(test='Mann-Whitney', text_format='star', loc='inside', verbose=2)
    annotator.apply_and_annotate()

def calculate_cross_sectional_associations(df_m00, roi_list, roi_names, wmh_vol, analysis_type="base"):

    df_m00_reg = df_m00.copy()
    
    # Create lists for the beta coefficients for disconnections and wmh...
    list_beta_wmh = []
    list_beta_amy_wmh = []
    list_beta_tau_wmh = []

    # Lower CI ..
    list_wmh_ci_lower = []
    list_amy_wmh_ci_lower = []
    list_tau_wmh_ci_lower = []
    
    # and upper CI
    list_wmh_ci_upper = []
    list_amy_wmh_ci_upper = []
    list_tau_wmh_ci_upper = []

    # Do the same for p-values...
    list_p_wmh = []
    list_p_amy_wmh = []
    list_p_tau_wmh = []
    
    # Number of analyzed subjects
    list_n_wmh = []

    for region in roi_list: 

        thickness = region + "_THICKNESS"
        disconn = region + "_disconn"
        disconn_ratio = region + "_disconn_ratio"
        amy = region + "_SUVR_amy"
        tau = region + "_SUVR_tau"
         
        if analysis_type == "basis":
            formula_wmh = f"{thickness} ~ {wmh_vol} + {amy} + {tau} + AGE + PTGENDER + C(HMHYPERT) + PTEDUCAT + RAVLT_forgetting + CTX_ETIV + TRACER_amy"

        elif analysis_type == "basis_disconn":
            # Only subjs with some 1% or more disconnection on average for that brain region
            df_m00_disconn = df_m00_reg[df_m00_reg[disconn] > 0.01].copy()
            # in the disconnectivity analyses we are interested only in that coefficient, so we change this to keep the same code
            wmh_vol = disconn_ratio
            formula_wmh = f"{thickness} ~ wmh_vol_log + {wmh_vol} + {amy} + {tau} + AGE + PTGENDER + C(HMHYPERT) + PTEDUCAT + RAVLT_forgetting + CTX_ETIV + TRACER_amy"

        elif analysis_type == "basis_zscore":
            df_m00_reg[thickness] = zscore(df_m00_reg[thickness])
            df_m00_reg[wmh_vol] = zscore(df_m00_reg[wmh_vol])
            formula_wmh = f"{thickness} ~ {wmh_vol} + {amy} + {tau} + AGE + PTGENDER + C(HMHYPERT) + PTEDUCAT + RAVLT_forgetting + CTX_ETIV + TRACER_amy"


        elif analysis_type == "sensitivity_less_covariates": #This analysis was requested by the reviewer but is not in the manuscript
            formula_wmh = f"{thickness} ~ {wmh_vol} + {amy} + {tau} + AGE + PTGENDER + CTX_ETIV + TRACER_amy"
        
        elif analysis_type == "global_tau_amy": #This analysis was requested by the reviewer but is not in the manuscript
            formula_wmh = f"{thickness} ~ {wmh_vol} + {amy} + {tau} + AGE + PTGENDER + C(HMHYPERT) + PTEDUCAT + RAVLT_forgetting + CTX_ETIV + TRACER_amy + global_tau + global_amy"
               
        if not "disconn" in wmh_vol:
            lr_wmh = smf.ols(data=df_m00_reg, formula=formula_wmh)
        else:
            lr_wmh = smf.ols(data=df_m00_disconn, formula=formula_wmh)

        res_wmh = lr_wmh.fit()
    
        # Extract parameter estimates and p-values for all-things WMH global
        list_p_wmh.append(res_wmh.pvalues[wmh_vol])
        list_beta_wmh.append(res_wmh.params[wmh_vol])
        list_p_amy_wmh.append(res_wmh.pvalues[amy])
        list_beta_amy_wmh.append(res_wmh.params[amy])
        list_p_tau_wmh.append(res_wmh.pvalues[tau])
        list_beta_tau_wmh.append(res_wmh.params[tau])

        wmh_ci = res_wmh.conf_int().loc[wmh_vol]
        amy_wmh_ci = res_wmh.conf_int().loc[amy]
        tau_wmh_ci = res_wmh.conf_int().loc[tau]

        # Append confidence intervals to the lists
        list_wmh_ci_lower.append(wmh_ci[0])
        list_wmh_ci_upper.append(wmh_ci[1])
        list_amy_wmh_ci_lower.append(amy_wmh_ci[0])
        list_amy_wmh_ci_upper.append(amy_wmh_ci[1])
        list_tau_wmh_ci_lower.append(tau_wmh_ci[0])
        list_tau_wmh_ci_upper.append(tau_wmh_ci[1])

        if not "disconn" in wmh_vol:
            list_n_wmh.append(df_m00_reg["PTID"].nunique())
        else:
            list_n_wmh.append(df_m00_disconn["PTID"].nunique())

    # Convert p-values to arrays and apply FDR correction
    arr_p_wmh = np.array(list_p_wmh)
    arr_p_amy_wmh = np.array(list_p_amy_wmh)
    arr_p_tau_wmh = np.array(list_p_tau_wmh)
    arr_p_wmh_fdr = false_discovery_control(arr_p_wmh)
    arr_p_amy_wmh_fdr = false_discovery_control(arr_p_amy_wmh)
    arr_p_tau_wmh_fdr = false_discovery_control(arr_p_tau_wmh)

    df_res_pretty_wmh = pd.DataFrame({
        "Region": [reg.capitalize()for reg in roi_names],
        "WMH_N": list_n_wmh,
        "WMH_beta": [f"{b:.3f} ({ci_lower:.3f}, {ci_upper:.3f})" for b, ci_lower, ci_upper in zip(list_beta_wmh, list_wmh_ci_lower, list_wmh_ci_upper)],
        "WMH_p_val": [f"{p:.2g}" if p >= 0.0001 else "p < 0.0001" for p in arr_p_wmh],
        "WMH_p_val_fdr": [f"{p:.2g}" if p >= 0.0001 else "p < 0.0001" for p in arr_p_wmh_fdr],
        "Amy_beta": [f"{b:.3f} ({ci_lower:.3f}, {ci_upper:.3f})" for b, ci_lower, ci_upper in zip(list_beta_amy_wmh, list_amy_wmh_ci_lower, list_amy_wmh_ci_upper)],
        "Amy_p_val": [f"{p:.2g}" if p >= 0.0001 else "p < 0.0001" for p in arr_p_amy_wmh],
        "Amy_p_val_fdr": [f"{p:.2g}" if p >= 0.0001 else "p < 0.0001" for p in arr_p_amy_wmh_fdr],
        "Tau_beta": [f"{b:.3f} ({ci_lower:.3f}, {ci_upper:.3f})" for b, ci_lower, ci_upper in zip(list_beta_tau_wmh, list_tau_wmh_ci_lower, list_tau_wmh_ci_upper)],
        "Tau_p_val": [f"{p:.2g}" if p >= 0.0001 else "p < 0.0001" for p in arr_p_tau_wmh],
        "Tau_p_val_fdr": [f"{p:.2g}" if p >= 0.0001 else "p < 0.0001" for p in arr_p_tau_wmh_fdr]
    })

    label = "GB-WMH" if analysis_type != "basis_disconn" else "FB-WMH"
    df_res_pretty_wmh.columns = pd.MultiIndex.from_tuples([

        ("", "region"),
        (label, "N"),
        (label, "Beta (95% CI)"),
        (label, "p-val"),
        (label, "p-val (FDR)"),
        
        ("Amyloid", "Beta (95% CI)"),
        ("Amyloid", "p-val"),
        ("Amyloid", "p-val (FDR)"),

        ("Tau", "Beta (95% CI)"),
        ("Tau", "p-val"),
        ("Tau", "p-val (FDR)"),
    ])

    return df_res_pretty_wmh

def zscore_also_longitudinal(df, biomarker):

    df_m00 = df[df["session"] == "M00"].copy()
    mean = df_m00[biomarker].mean()
    st_dev = df_m00[biomarker].std()

    zscored_biomarker = (df[biomarker] - mean) / st_dev

    return zscored_biomarker

def calculate_longitudinal_associations(df, roi_list, roi_names, wmh_vol, analysis_type="basis"):

    subjs_with_fu = df[df["session"] == "M01"]["PTID"].unique()
    df_long = df[df["PTID"].isin(subjs_with_fu)].copy()
    # Get baseline WMH for each subject
    baseline_wmh = df_long[df_long["session"] == "M00"][["PTID", wmh_vol]].rename(columns={f"{wmh_vol}": f"{wmh_vol}_baseline"})
    baseline_ravlt = df_long[df_long["session"] == "M00"][["PTID", "RAVLT_forgetting"]].rename(columns={"RAVLT_forgetting": f"RAVLT_forgetting_baseline"})
    # Merge back into long dataframe
    df_long = df_long.merge(baseline_wmh, on="PTID", how="left")
    df_long = df_long.merge(baseline_ravlt, on="PTID", how="left")

    list_beta_wmh = []
    list_p_wmh = []
    list_wmh_ci_lower = []
    list_wmh_ci_upper = []

    list_n_wmh = []

    for region in roi_list: 

        disconn = region + "_disconn"
        disconn_ratio = region + "_disconn_ratio"
        thickness = region + "_THICKNESS"
        amy = region + "_SUVR_amy"
        tau = region + "_SUVR_tau"
        
        baseline_amy = df_long[df_long["session"] == "M00"][["PTID", amy]].rename(columns={f"{amy}": f"{amy}_baseline"})
        baseline_tau = df_long[df_long["session"] == "M00"][["PTID", tau]].rename(columns={f"{tau}": f"{tau}_baseline"})
        baseline_disconn_ratio = df_long[df_long["session"] == "M00"][["PTID", disconn_ratio]].rename(columns={f"{disconn_ratio}": f"{disconn_ratio}_baseline"})
        
        df_long = df_long.merge(baseline_amy, on="PTID", how="left")
        df_long = df_long.merge(baseline_tau, on="PTID", how="left")
        df_long = df_long.merge(baseline_disconn_ratio, on="PTID", how="left")

        df_m00 = df[df["session"] == "M00"].copy()

        df_m00_disconn_ptids = df_m00[df_m00[disconn] > 0.01]["PTID"] # subjects with no disconnections at the first session
        df_long_disconn = df_long[df_long["PTID"].isin(df_m00_disconn_ptids)].copy()
        
        if analysis_type == "basis":
            formula_wmh = f"{thickness} ~ {wmh_vol}_baseline*Years_m00 + {amy}_baseline + {tau}_baseline + AGE + PTGENDER + C(HMHYPERT) + RAVLT_forgetting_baseline + PTEDUCAT + CTX_ETIV + TRACER_amy"
        
        elif analysis_type == "basis_disconn":
            # in the disconnectivity analyses we are interested only in that coefficient, so we change this to keep the same code
            wmh_vol = disconn_ratio
            formula_wmh = f"{thickness} ~ wmh_vol_log + {wmh_vol}_baseline*Years_m00 + {amy}_baseline + {tau}_baseline + AGE + PTGENDER + C(HMHYPERT) + RAVLT_forgetting_baseline + PTEDUCAT + CTX_ETIV + TRACER_amy"

        elif analysis_type == "sensitivity_less_covariates": # asked by reviewer, but not included
            formula_wmh = f"{thickness} ~ {wmh_vol}_baseline*Years_m00 + {amy}_baseline + {tau}_baseline + AGE + PTGENDER + CTX_ETIV + TRACER_amy"

        elif analysis_type == "global_tau_amy": # asked by reviewer, but not included
            formula_wmh = f"{thickness} ~ {wmh_vol}_baseline*Years_m00 + {amy}_baseline + {tau}_baseline + AGE + PTGENDER + C(HMHYPERT) + PTEDUCAT + RAVLT_forgetting_baseline + CTX_ETIV + TRACER_amy + global_tau + global_amy"
        
        elif analysis_type == "tau": # asked by reviewer, but not included
            formula_wmh = f"{tau} ~ {wmh_vol}_baseline*Years_m00 + {amy}_baseline + AGE + PTGENDER + C(HMHYPERT) + PTEDUCAT + RAVLT_forgetting_baseline + CTX_ETIV + TRACER_amy"
            
        # Fit the models
        if not "disconn" in analysis_type:
            lr_wmh = smf.mixedlm(data=df_long, formula=formula_wmh, groups=df_long["PTID"])
            res_wmh = lr_wmh.fit()
        else:
            lr_wmh = smf.mixedlm(data=df_long_disconn, formula=formula_wmh, groups=df_long_disconn["PTID"])
            res_wmh = lr_wmh.fit()
        
        list_p_wmh.append(res_wmh.pvalues[f"{wmh_vol}_baseline:Years_m00"])
        list_beta_wmh.append(res_wmh.params[f"{wmh_vol}_baseline:Years_m00"])

        # Extract 95% confidence intervals for each parameter
        wmh_ci = res_wmh.conf_int().loc[f"{wmh_vol}_baseline:Years_m00"]

        list_wmh_ci_lower.append(wmh_ci[0])
        list_wmh_ci_upper.append(wmh_ci[1])
        if not "disconn" in analysis_type:
            list_n_wmh.append(df_long["PTID"].nunique())
        else:
            list_n_wmh.append(df_long_disconn["PTID"].nunique())

    # Convert p-values to arrays and apply FDR correction
    arr_p_wmh = np.array(list_p_wmh)

    # arr_p_disconn_fdr = false_discovery_control(arr_p_disconn)
    arr_p_wmh_fdr = false_discovery_control(arr_p_wmh)

    df_res_pretty_wmh = pd.DataFrame({
        "Region": [reg.capitalize()for reg in roi_names],
        "WMH-global_N": list_n_wmh,
        "WMH-global_beta_wmh": [f"{b:.3f} ({ci_lower:.3f}, {ci_upper:.3f})" for b, ci_lower, ci_upper in zip(list_beta_wmh, list_wmh_ci_lower, list_wmh_ci_upper)],
        "WMH-global_p_val": [f"{p:.2f}" for p in arr_p_wmh],
        "WMH-global_p_val_fdr": [f"{p:.2f}" for p in arr_p_wmh_fdr],
        })
    label = "GB-WMH" if analysis_type != "basis_disconn" else "FB-WMH"
    df_res_pretty_wmh.columns = pd.MultiIndex.from_tuples([
        ("", "region"),
        (label, "N"),
        (label, "Beta (95% CI)"),
        (label, "p-value"),
        (label, "FDR-p"),
    ])

    return df_res_pretty_wmh


# Fix the parser to correctly handle strings with multiple spaces
def parse_beta_ci(series):
    betas = []
    errors = []
    for s in series:
        beta_str = s.split(" ")[0]
        ci_str = s[s.find("(")+1:s.find(")")]
        beta = float(beta_str)
        ci_low, ci_high = map(float, ci_str.split(","))
        err_low = beta - ci_low
        err_high = ci_high - beta
        betas.append(beta)
        errors.append((err_low, err_high))
    return np.array(betas), np.array(errors).T  # transpose to get (2, N)