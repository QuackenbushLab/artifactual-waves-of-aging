import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import glob
import re
from functools import reduce

import pandas as pd
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.append(PROJECT_ROOT)

from de_SWAN_helpers import generate_exp_data, run_all_de_swan
import mild_configs as config

# Generate data
os.makedirs(config.OUTDIR, exist_ok=True)

for sim in range(config.N_SIMS):

    X = np.random.uniform(config.MIN_AGE, config.MAX_AGE, size=config.N)

    n_outliers = config.N_OUTLIERS  # 5% of samples as outliers
    age_center = np.median(X)

    center_order = np.argsort(np.abs(X - age_center))
    outlier_inds = center_order[:n_outliers]

    exp_df = generate_exp_data(
        n=config.N,
        mols_num=config.MOLS_NUM,
        true_slope_sd=config.TRUE_SLOPE_SD,
        true_intercept_sd=config.TRUE_INTERCEPT_SD,
        noise_sd=config.NOISE_SD,
        dist=config.DIST,
        noise_trend=config.NOISE_TREND,
        unif_lower = config.MIN_AGE,
        unif_upper = config.MAX_AGE,
        outlier_multiplier=config.OUTLIER_MULTIPLIER,
        outlier_inds=outlier_inds,
        fixed_x = X,
        )

    # Optional: save example molecule plot only for sim 0
    if sim == 0:
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

        random_mols = np.random.choice(
            [f"mol_{i}" for i in range(config.MOLS_NUM)],
            size=6,
            replace=False
        )

        fig, axes = plt.subplots(
            2, 3,
            figsize=(3.2, 2.15),
            sharex=True,
            sharey=True
        )

        axes = axes.flatten()

        for ax, mol in zip(axes, random_mols):
            ax.scatter(
                exp_df["age"],
                exp_df[mol],
                alpha=0.45,
                s=4,
                linewidths=0
            )

            ax.tick_params(axis="both", width=0.8, length=2.5)

            # No individual axis labels; use shared labels below
            ax.set_xlabel("")
            ax.set_ylabel("")

        # Hide x tick labels on the top row
        for ax in axes[:3]:
            ax.tick_params(labelbottom=False)

        # Hide y tick labels on non-left columns
        for ax in [axes[1], axes[2], axes[4], axes[5]]:
            ax.tick_params(labelleft=False)

        # Shared labels for the full mini-panel
        fig.supxlabel("Age", y=0.02, fontsize=9)
        fig.supylabel("Expression", x=0.02, fontsize=9)

        fig.tight_layout(pad=0.35, w_pad=0.35, h_pad=0.35)

        fig.savefig(
            os.path.join(config.OUTDIR, "random_6_molecules_2x3.pdf"),
            bbox_inches="tight"
        )

    # Run DE-SWAN
    de_swan_df = run_all_de_swan(
        mols_num=config.MOLS_NUM,
        exp_df=exp_df,
        window=config.WINDOW,
        age_col="age",
        min_age=config.MIN_AGE,
        max_age=config.MAX_AGE
    )


    de_swan_df = de_swan_df[["midpoint", "n_sig"]].copy()
    de_swan_df = de_swan_df.rename(columns={"n_sig": f"sim_{sim}_n_sig"})

    out_file = os.path.join(config.OUTDIR, f"de_swan_sim_{sim}_n={config.N}_p={config.MOLS_NUM}.csv")
    de_swan_df.to_csv(out_file, index=False)

    print(f"Finished simulation {sim}. Saved to {out_file}")


##### Combine individual DE-SWAN simulation result files
files = sorted(
    glob.glob(os.path.join(config.OUTDIR, f"de_swan_sim_*_n={config.N}_p={config.MOLS_NUM}.csv"))
)
print(f"Found {len(files)} files")

if len(files) == 0:
    raise ValueError(f"No files found in {config.OUTDIR}")

dfs = []

for file in files:
    df = pd.read_csv(file)
    basename = os.path.basename(file)

    match = re.search(r"de_swan_sim_(\d+)_n=" + str(config.N) + r"_p=" + str(config.MOLS_NUM) + r"\.csv", basename)

    if match is None:
        raise ValueError(f"Could not extract sim number from {basename}")

    sim_id = match.group(1)

    if "n_sig" in df.columns:
        tmp = df[["midpoint", "n_sig"]].copy()
        tmp = tmp.rename(columns={"n_sig": f"sim_{sim_id}_n_sig"})

    else:
        sim_cols = [
            col for col in df.columns
            if col.startswith("sim_") and col.endswith("_n_sig")
        ]

        if len(sim_cols) != 1:
            raise ValueError(
                f"Expected one sim_*_n_sig column in {basename}, found {sim_cols}. "
                f"Columns are: {df.columns.tolist()}"
            )

        tmp = df[["midpoint", sim_cols[0]]].copy()
        tmp = tmp.rename(columns={sim_cols[0]: f"sim_{sim_id}_n_sig"})

    dfs.append(tmp)


master_df = reduce(
    lambda left, right: left.merge(right, on="midpoint", how="outer"),
    dfs
)

master_df = master_df.sort_values("midpoint").reset_index(drop=True)


##### Summarize across simulations
sim_cols = [
    col for col in master_df.columns
    if col.startswith("sim_") and col.endswith("_n_sig")
]

master_df["mean_n_sig"] = master_df[sim_cols].mean(axis=1)
master_df["sd_n_sig"] = master_df[sim_cols].std(axis=1)

combined_rslts = os.path.join(config.OUTDIR, f"de_swan_n={config.N}_p={config.MOLS_NUM}.csv")
master_df.to_csv(combined_rslts, index=False)

