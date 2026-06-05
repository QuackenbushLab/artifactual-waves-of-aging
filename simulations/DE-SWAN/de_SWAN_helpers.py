import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

import statsmodels.api as sm
import statsmodels.formula.api as smf

from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import os

def run_de_SWAN(window, age_col, exp_data, min_age=40, max_age=65, wilcoxon=False, covariates=None):
    """
    Simple DE-SWAN-style sliding window test.

    For each whole-number midpoint:
        young group: [midpt - window, midpt)
        old group:   [midpt, midpt + window)

    For each molecule:
        expression ~ age_group + optional covariates

    Then BH-adjust p-values across molecules at each midpoint.
    """

    if covariates is None:
        covariates = []
    results = []
    molecules = sorted(exp_data["molecule"].unique())
    half_window = window / 2
    for midpt in range(min_age, max_age):  # whole-number midpoints
        pvals = []

        for mol in molecules:
            exp_data_m = exp_data[exp_data["molecule"] == mol].copy()
            age_m = exp_data_m[age_col]

            young_df = exp_data_m[(age_m >= midpt - half_window) & (age_m < midpt)].copy()
            old_df = exp_data_m[(age_m >= midpt) & (age_m < midpt + half_window)].copy()
            young_df["age_group"] = 0
            old_df["age_group"] = 1

            window_df = pd.concat([young_df, old_df], axis=0).copy()

            try:
                if wilcoxon:
                    print("Running Mann-Whitney U test for molecule", mol, "at midpoint", midpt)
                    window_df = window_df.dropna(subset=["expression", "age_group"])

                    young_vals = window_df.loc[window_df["age_group"] == 0, "expression"].values

                    old_vals = window_df.loc[window_df["age_group"] == 1, "expression"].values

                    _, p = mannwhitneyu(young_vals, old_vals,alternative="two-sided")

                else:
                    window_df = window_df.dropna(
                        subset=["expression", "age_group"] + covariates
                    ) 

                    rhs_terms = ["age_group"] + covariates
                    formula = "expression ~ " + " + ".join(rhs_terms)

                    fit = smf.ols(formula=formula, data=window_df).fit()
                    p = fit.pvalues["age_group"]

            except Exception:
                p = np.nan

            pvals.append(p)

        pvals = np.array(pvals, dtype=float)
        valid = ~np.isnan(pvals)
        pvals_adj = np.full_like(pvals, np.nan, dtype=float)

        if valid.sum() > 0:
            pvals_adj[valid] = multipletests(pvals[valid], method="fdr_bh")[1]

        n_sig = np.sum(pvals_adj < 0.05)

        results.append({"midpoint": midpt, "n_sig": n_sig})

    de_swan_df = pd.DataFrame(results)

    return de_swan_df

def run_all_de_swan(mols_num, exp_df, window=20, age_col="age", min_age=25, max_age=75, wilcoxon_indicator=False):
    # Reshape exp_df to have molecule column and expression column
    exp_df_reshaped = []
    for mol_idx in range(mols_num):
        mol_df = exp_df[["age", f"mol_{mol_idx}"]].copy()
        mol_df["molecule"] = mol_idx
        mol_df.rename(columns={f"mol_{mol_idx}": "expression"}, inplace=True)
        exp_df_reshaped.append(mol_df)

    exp_df_long = pd.concat(exp_df_reshaped, ignore_index=True)

    de_swan_df = run_de_SWAN(
        window=window,
        age_col=age_col,
        exp_data=exp_df_long,
        min_age=min_age,
        max_age=max_age,
        wilcoxon=wilcoxon_indicator
    )

    return de_swan_df

def plot_de_swan_curves(
    de_swan_df,
    path,
    y_lim_lower=0,
    y_lim_upper=None,
    filename="de_swan_results"
):
    os.makedirs(path, exist_ok=True)

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, ax = plt.subplots(figsize=(5.5, 3.2))

    ax.plot(
        de_swan_df["midpoint"],
        de_swan_df["n_sig"],
        marker="o",
        markersize=3,
        linewidth=1.8,
        label="DE-SWAN hits"
    )

    ax.set_xlabel("Age midpoint")
    ax.set_ylabel("Significant molecules")
    ax.set_title("DE-SWAN simulation results")
    ax.set_ylim(y_lim_lower, y_lim_upper)
    ax.tick_params(axis="both", width=0.8, length=3)

    ax.legend(
        frameon=False,
        fontsize=7,
        handlelength=1.6,
        borderaxespad=0.3,
        loc="best"
    )

    fig.tight_layout()

    fig.savefig(
        os.path.join(path, f"{filename}.pdf"),
        bbox_inches="tight"
    )

    fig.savefig(
        os.path.join(path, f"{filename}.png"),
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

def run_de_SWAN(window, age_col, exp_data, min_age=40, max_age=65, wilcoxon=False, covariates=None):
    """
    Simple DE-SWAN-style sliding window test.

    For each whole-number midpoint:
        young group: [midpt - window/2, midpt)
        old group:   [midpt, midpt + window/2)

    For each molecule:
        expression ~ age_group + optional covariates

    Then BH-adjust p-values across molecules at each midpoint.
    """
    
    half_window = window / 2

    if covariates is None:
        covariates = []
    results = []
    molecules = sorted(exp_data["molecule"].unique())

    for midpt in range(min_age, max_age):  # whole-number midpoints
        pvals = []

        for mol in molecules:
            exp_data_m = exp_data[exp_data["molecule"] == mol].copy()
            age_m = exp_data_m[age_col]

            young_df = exp_data_m[(age_m >= midpt - half_window) & (age_m < midpt)].copy()
            old_df = exp_data_m[(age_m >= midpt) & (age_m < midpt + half_window)].copy()
            young_df["age_group"] = 0
            old_df["age_group"] = 1

            window_df = pd.concat([young_df, old_df], axis=0).copy()

            try:
                if wilcoxon:
                    print("Running Mann-Whitney U test for molecule", mol, "at midpoint", midpt)
                    window_df = window_df.dropna(subset=["expression", "age_group"])

                    young_vals = window_df.loc[window_df["age_group"] == 0, "expression"].values

                    old_vals = window_df.loc[window_df["age_group"] == 1, "expression"].values

                    _, p = mannwhitneyu(young_vals, old_vals,alternative="two-sided")

                else:
                    window_df = window_df.dropna(
                        subset=["expression", "age_group"] + covariates
                    ) 

                    rhs_terms = ["age_group"] + covariates
                    formula = "expression ~ " + " + ".join(rhs_terms)

                    fit = smf.ols(formula=formula, data=window_df).fit()
                    p = fit.pvalues["age_group"]

            except Exception:
                p = np.nan

            pvals.append(p)

        pvals = np.array(pvals, dtype=float)
        valid = ~np.isnan(pvals)
        pvals_adj = np.full_like(pvals, np.nan, dtype=float)

        if valid.sum() > 0:
            pvals_adj[valid] = multipletests(pvals[valid], method="fdr_bh")[1]

        n_sig = np.sum(pvals_adj < 0.05)

        results.append({"midpoint": midpt, "n_sig": n_sig})

    de_swan_df = pd.DataFrame(results)

    return de_swan_df

def run_all_de_swan(mols_num, exp_df, window=20, age_col="age", min_age=25, max_age=75, wilcoxon_indicator=False):
    # Reshape exp_df to have molecule column and expression column
    exp_df_reshaped = []
    for mol_idx in range(mols_num):
        mol_df = exp_df[["age", f"mol_{mol_idx}"]].copy()
        mol_df["molecule"] = mol_idx
        mol_df.rename(columns={f"mol_{mol_idx}": "expression"}, inplace=True)
        exp_df_reshaped.append(mol_df)

    exp_df_long = pd.concat(exp_df_reshaped, ignore_index=True)

    de_swan_df = run_de_SWAN(
        window=window,
        age_col=age_col,
        exp_data=exp_df_long,
        min_age=min_age,
        max_age=max_age,
        wilcoxon=wilcoxon_indicator
    )

    return de_swan_df

def plot_de_swan_curves(
    de_swan_df,
    path,
    y_lim_lower=None,
    y_lim_upper=None,
    filename="de_swan_results"
    ):
    os.makedirs(path, exist_ok=True)

    plt.rcParams.update({
        "font.family": "serif",
        "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 7,
        "axes.linewidth": 0.8,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig, ax = plt.subplots(figsize=(5.5, 3.2))

    ax.plot(
        de_swan_df["midpoint"],
        de_swan_df["n_sig"],
        marker="o",
        markersize=3,
        linewidth=1.8,
        label="DE-SWAN hits"
    )

    ax.set_xlabel("Age midpoint")
    ax.set_ylabel("Significant molecules")
    ax.set_title("DE-SWAN simulation results")
    ax.set_ylim(y_lim_lower, y_lim_upper)
    ax.tick_params(axis="both", width=0.8, length=3)

    ax.legend(
        frameon=False,
        fontsize=7,
        handlelength=1.6,
        borderaxespad=0.3,
        loc="best"
    )

    fig.tight_layout()

    fig.savefig(
        os.path.join(path, f"{filename}.pdf"),
        bbox_inches="tight"
    )

    fig.savefig(
        os.path.join(path, f"{filename}.png"),
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)

def generate_exp_data(
    n,
    mols_num,
    true_slope_sd=0.1,
    true_intercept_sd=0.1,
    noise_sd=0.85,
    dist="unif",
    noise_trend=None,
    norm_mean=None,
    norm_sd=None,
    unif_lower=None,
    unif_upper=None,
    mixture_prop=None,
    mixture_mean1=None,
    mixture_sd1=None,
    mixture_mean2=None,
    mixture_sd2=None,
    noise_shape="u",
    flip=False,
    noise_power=2.0,
    min_noise_multiplier=0.5,
    max_noise_multiplier=1.5,
    exp_strength=5.0,
    fixed_x = None,
    outlier_inds=None,
    outlier_multiplier = 3.0
    ):

    if fixed_x is not None:
        X = fixed_x
    else:
        if dist == "unif":
            X = np.random.uniform(unif_lower, unif_upper, size=n)

        elif dist == "gaussian":
            X = np.random.normal(norm_mean, norm_sd, size=n)

        elif dist == "bimodal":
            n1 = int(n * mixture_prop)
            n2 = n - n1

            X1 = np.random.normal(loc=mixture_mean1, scale=mixture_sd1, size=n1)
            X2 = np.random.normal(loc=mixture_mean2, scale=mixture_sd2, size=n2)

            X = np.concatenate([X1, X2])
            np.random.shuffle(X)

        else:
            raise ValueError(f"Unknown dist: {dist}")

    age_scaled = (X - X.min()) / (X.max() - X.min())

    if noise_trend == "increasing":

        exp_scaled = (np.exp(exp_strength * age_scaled) - 1) / (
            np.exp(exp_strength) - 1
        )

        noise_sd_by_age = noise_sd * (
            min_noise_multiplier
            + (max_noise_multiplier - min_noise_multiplier) * exp_scaled
        )

    elif noise_trend == "exp_decreasing":

        age_scaled = (X - X.min()) / (X.max() - X.min())

        exp_scaled = (np.exp(exp_strength * age_scaled) - 1) / (
            np.exp(exp_strength) - 1
        )

        noise_sd_by_age = noise_sd * (
            max_noise_multiplier
            - (max_noise_multiplier - min_noise_multiplier) * exp_scaled
        )
    elif noise_trend == "exp_increasing":

        age_scaled = (X - X.min()) / (X.max() - X.min())

        exp_scaled = (np.exp(exp_strength * age_scaled) - 1) / (
            np.exp(exp_strength) - 1
        )

        noise_sd_by_age = noise_sd * (
            min_noise_multiplier
            + (max_noise_multiplier - min_noise_multiplier) * exp_scaled
        )
    elif noise_trend == "decreasing":
        noise_sd_by_age = noise_sd * (max_noise_multiplier - age_scaled)

    elif noise_trend == "ends":
        age_mid = (X.min() + X.max()) / 2
        age_half_range = (X.max() - X.min()) / 2
        age_dist_from_mid = np.abs(X - age_mid) / age_half_range

        min_multiplier = min_noise_multiplier
        max_multiplier = max_noise_multiplier

        if noise_shape == "v":
            # Linear V-shape
            shape = age_dist_from_mid

        elif noise_shape == "u":
            # Smooth U-shape
            shape = age_dist_from_mid ** noise_power

        elif noise_shape == "laplace_u":
            # Exponential/Laplace-like U-shape:
            # low in middle, high at ends, nonlinear rise toward edges
            shape = 1 - np.exp(-laplace_decay * age_dist_from_mid)
            shape = shape / shape.max()

        else:
            raise ValueError(f"Unknown noise_shape: {noise_shape}")

        if flip:
            shape = 1 - shape

        noise_multiplier = min_multiplier + (
            max_multiplier - min_multiplier
        ) * shape

        noise_sd_by_age = noise_sd * noise_multiplier
    elif noise_trend == "middle":
        age_mid = (X.min() + X.max()) / 2
        age_half_range = (X.max() - X.min()) / 2
        age_dist_from_mid = np.abs(X - age_mid) / age_half_range
        noise_sd_by_age = noise_sd * (1.5 - age_dist_from_mid)

    else:
        noise_sd_by_age = np.repeat(noise_sd, n)

    data = {
        "age": X,
    }

    for mol_idx in range(mols_num):
        true_slope = np.random.normal(loc=0, scale=true_slope_sd)
        true_intercept = np.random.normal(loc=0, scale=true_intercept_sd)

        noise = np.random.normal(loc=0, scale=noise_sd_by_age, size=n)
        y = true_slope * X + true_intercept + noise

        if outlier_inds is not None:
            outlier_mask = np.zeros(n, dtype=bool)
            outlier_mask[outlier_inds] = True

            y[outlier_mask] += np.random.normal(
                loc=0,
                scale=outlier_multiplier * noise_sd,
                size=outlier_mask.sum()
            )

        data[f"mol_{mol_idx}"] = y
    exp_df = pd.DataFrame(data)

    return exp_df

def get_window_realized_noise_variance(
    exp_df,
    window,
    age_col="age",
    min_age=25,
    max_age=75
    ):
    half_window = window / 2
    noise_cols = [
        col for col in exp_df.columns
        if col.startswith("noise_")
        and col not in ["noise_sd_by_age", "noise_var_by_age"]
    ]

    if len(noise_cols) == 0:
        raise ValueError(
            "No realized noise columns found. Make sure exp_df has columns like noise_0, noise_1, ..."
        )

    results = []

    for midpt in range(min_age, max_age):
        window_df = exp_df[
            (exp_df[age_col] >= midpt - half_window) &
            (exp_df[age_col] < midpt + half_window)
        ]

        # All realized noise values for all people x all molecules in this window
        window_noise_values = window_df[noise_cols].to_numpy().ravel()

        # Drop NaNs just in case
        window_noise_values = window_noise_values[~np.isnan(window_noise_values)]

        if len(window_noise_values) > 1:
            realized_noise_var = np.var(window_noise_values, ddof=1)
            realized_noise_sd = np.std(window_noise_values, ddof=1)
        else:
            realized_noise_var = np.nan
            realized_noise_sd = np.nan

        results.append({
            "midpoint": midpt,
            "realized_noise_var": realized_noise_var,
            "realized_noise_sd": realized_noise_sd,
            "n_people": window_df.shape[0],
            "n_noise_values": len(window_noise_values)
        })

    return pd.DataFrame(results)