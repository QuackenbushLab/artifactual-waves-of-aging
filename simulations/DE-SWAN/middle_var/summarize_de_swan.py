import os
import sys
import glob
import re
from functools import reduce

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

##### Import config from parent directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MIDDLE_VAR_DIST_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))

import middle_var_configs as config


if MIDDLE_VAR_DIST_DIR not in sys.path:
    sys.path.insert(0, MIDDLE_VAR_DIST_DIR)

##### Paths
results_dir = config.RESULTS_DIR

plot_dir = os.path.join(results_dir, "plots")
os.makedirs(plot_dir, exist_ok=True)

y_lower = 500
y_upper = 1700
x_lower = config.MIN_AGE + 3
x_upper = config.MAX_AGE - 3
fig_width = 5.0
fig_height = 4

##### Publication-style matplotlib settings for LaTeX / Overleaf
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 13,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "axes.linewidth": 0.9,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


plot_dir = os.path.join(results_dir, "plots")
os.makedirs(plot_dir, exist_ok=True)


##### Helper functions

def save_pdf_png(fig, plot_dir, filename):
    fig.savefig(
        os.path.join(plot_dir, f"{filename}.pdf"),
        bbox_inches="tight"
    )

    fig.savefig(
        os.path.join(plot_dir, f"{filename}.png"),
        dpi=300,
        bbox_inches="tight"
    )


##### Combine individual DE-SWAN simulation result files
files = sorted(
    glob.glob(os.path.join(results_dir, f"de_swan_sim_*_n={config.N}_p={config.MOLS_NUM}.csv"))
)
print(f"Found {len(files)} files")

if len(files) == 0:
    raise ValueError(f"No files found in {results_dir}")

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


##### Save combined results before plotting
out_file = os.path.join(results_dir, f"de_swan_all_sims_n={config.N}_p={config.MOLS_NUM}.csv")
master_df.to_csv(out_file, index=False)

print(f"Saved combined results to {out_file}")
print(master_df.head())


##### Plot 1: all DE-SWAN simulation curves with mean
fig, ax = plt.subplots(figsize=(fig_width, fig_height))
# Plotting-only filter: remove first edge artifact but keep master_df unchanged
plot_df = master_df[master_df["midpoint"] > config.MIN_AGE].copy()
for col in sim_cols:
    ax.plot(
        plot_df["midpoint"],
        plot_df[col],
        alpha=0.12,
        linewidth=0.8
    )

ax.plot(
    plot_df["midpoint"],
    plot_df["mean_n_sig"],
    linewidth=2.2,
    label="Mean across simulations"
)

ax.set_xlabel("Age midpoint")
ax.set_xlim(x_lower, x_upper)
ax.set_ylabel("Significant molecules")
ax.set_ylim(y_lower, y_upper)
ax.set_title(f"DE-SWAN curves across {len(files)} simulations")
ax.tick_params(axis="both", width=0.8, length=3)

fig.tight_layout()
save_pdf_png(fig, plot_dir, "de_swan_all_curves_with_mean")
plt.show()


# ##### Plot 2: mean DE-SWAN curve overlaid with uniform age density
x_age = np.linspace(
    plot_df["midpoint"].min(),
    plot_df["midpoint"].max(),
    500
)

if config.NOISE_TREND == "ends":
    age_mid = (x_age.min() + x_age.max()) / 2
    age_half_range = (x_age.max() - x_age.min()) / 2
    age_dist_from_mid = np.abs(x_age - age_mid) / age_half_range

    min_multiplier = config.MIN_NOISE_MULTIPLIER
    max_multiplier = config.MAX_NOISE_MULTIPLIER

    if config.NOISE_SHAPE == "v":
        # Linear V-shape
        shape = age_dist_from_mid
    else:
        raise ValueError(f"Unknown noise_shape: {config.NOISE_SHAPE}")

    if config.FLIP:
        shape = 1 - shape

    noise_multiplier = min_multiplier + (
        max_multiplier - min_multiplier
    ) * shape

    noise_sd_by_age = config.NOISE_SD * noise_multiplier
fig, ax1 = plt.subplots(figsize=(fig_width, fig_height))

ax1.plot(
    plot_df["midpoint"],
    plot_df["mean_n_sig"],
    linewidth=2.2,
    label="Mean DE-SWAN hits"
)

ax1.fill_between(
    plot_df["midpoint"],
    plot_df["mean_n_sig"] - plot_df["sd_n_sig"],
    plot_df["mean_n_sig"] + plot_df["sd_n_sig"],
    alpha=0.18,
    label="Mean ± 1 SD"
)

ax1.set_xlabel("Age midpoint")
ax1.set_xlim(x_lower, x_upper)
ax1.set_ylabel("Significant molecules")
ax1.set_ylim(y_lower, y_upper)
ax1.tick_params(axis="both", width=0.8, length=3)

ax2 = ax1.twinx()

ax2.plot(
    x_age,
    noise_sd_by_age,
    linestyle="--",
    linewidth=1.8,
    label=rf"Noise Trend"
)

ax2.set_ylabel("Noise SD by age")
ax2.tick_params(axis="both", width=0.8, length=3)

ax1.set_title("DE-SWAN hits overlaid with noise trend")

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

fig.tight_layout()
save_pdf_png(fig, plot_dir, "de_swan_mean_with_noise_overlay")
plt.show()

