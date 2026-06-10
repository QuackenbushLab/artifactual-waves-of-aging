import os
import sys
import glob
import re
from functools import reduce

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import norm


##### Import config from parent directory
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
UNIF_DIST_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))

if UNIF_DIST_DIR not in sys.path:
    sys.path.insert(0, UNIF_DIST_DIR)

import bimodal_dist_configs as config

##### Paths
results_dir = config.RESULTS_DIR

plot_dir = os.path.join(results_dir, "plots")
os.makedirs(plot_dir, exist_ok=True)

y_lower = 0
y_upper = 2000
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

##### Helper functions

import numpy as np
from scipy.stats import norm


def truncated_normal_mean(lower, upper, mu, sd):
    """
    E[X | lower <= X < upper] for X ~ Normal(mu, sd).
    """
    a = (lower - mu) / sd
    b = (upper - mu) / sd

    prob = norm.cdf(b) - norm.cdf(a)

    if prob == 0:
        return np.nan

    return mu + sd * (norm.pdf(a) - norm.pdf(b)) / prob


def truncated_bimodal_mean(lower, upper, prop, mu1, sd1, mu2, sd2):
    """
    E[X | lower <= X < upper] for 
    X ~ prop * N(mu1, sd1) + (1 - prop) * N(mu2, sd2).
    """

    # Probability each component puts in the interval
    p1 = norm.cdf(upper, loc=mu1, scale=sd1) - norm.cdf(lower, loc=mu1, scale=sd1)
    p2 = norm.cdf(upper, loc=mu2, scale=sd2) - norm.cdf(lower, loc=mu2, scale=sd2)

    # Mixture probability of falling in the interval
    p_mix = prop * p1 + (1 - prop) * p2

    if p_mix == 0:
        return np.nan

    # Component-specific truncated means
    m1 = truncated_normal_mean(lower, upper, mu1, sd1)
    m2 = truncated_normal_mean(lower, upper, mu2, sd2)

    # Posterior mixture weights inside this interval
    w1 = prop * p1 / p_mix
    w2 = (1 - prop) * p2 / p_mix

    return w1 * m1 + w2 * m2


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
    glob.glob(
        os.path.join(
            results_dir,
            f"de_swan_sim_*_n={config.N}_p={config.MOLS_NUM}.csv"
        )
    )
)

print(
        os.path.join(
            results_dir,
            f"de_swan_sim_*_n={config.N}_p={config.MOLS_NUM}.csv"
        )
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


##### Compute bimodal-distribution diagnostics

# Probability that age falls inside the full DE-SWAN window:
# [midpoint - WINDOW, midpoint + WINDOW)
master_df["bimodal_window_prob"] = (
    config.MIXTURE_PROP1 * (
        norm.cdf(
            master_df["midpoint"] + config.WINDOW / 2,
            loc=config.MIXTURE_MEAN1,
            scale=config.MIXTURE_SD1
        )
        -
        norm.cdf(
            master_df["midpoint"] - config.WINDOW / 2,
            loc=config.MIXTURE_MEAN1,
            scale=config.MIXTURE_SD1
        )
    )
    +
    (1 - config.MIXTURE_PROP1) * (
        norm.cdf(
            master_df["midpoint"] + config.WINDOW / 2,
            loc=config.MIXTURE_MEAN2,
            scale=config.MIXTURE_SD2
        )
        -
        norm.cdf(
            master_df["midpoint"] - config.WINDOW / 2,
            loc=config.MIXTURE_MEAN2,
            scale=config.MIXTURE_SD2
        )
    )
)

p_young = (
    config.MIXTURE_PROP1 * (
        norm.cdf(
            master_df["midpoint"],
            loc=config.MIXTURE_MEAN1,
            scale=config.MIXTURE_SD1
        )
        -
        norm.cdf(
            master_df["midpoint"] - config.WINDOW / 2,
            loc=config.MIXTURE_MEAN1,
            scale=config.MIXTURE_SD1
        )
    )
    +
    (1 - config.MIXTURE_PROP1) * (
        norm.cdf(
            master_df["midpoint"],
            loc=config.MIXTURE_MEAN2,
            scale=config.MIXTURE_SD2
        )
        -
        norm.cdf(
            master_df["midpoint"] - config.WINDOW / 2,
            loc=config.MIXTURE_MEAN2,
            scale=config.MIXTURE_SD2
        )
    )
)

p_old = (
    config.MIXTURE_PROP1 * (
        norm.cdf(
            master_df["midpoint"] + config.WINDOW / 2,
            loc=config.MIXTURE_MEAN1,
            scale=config.MIXTURE_SD1
        )
        -
        norm.cdf(
            master_df["midpoint"],
            loc=config.MIXTURE_MEAN1,
            scale=config.MIXTURE_SD1
        )
    )
    +
    (1 - config.MIXTURE_PROP1) * (
        norm.cdf(
            master_df["midpoint"] + config.WINDOW / 2,
            loc=config.MIXTURE_MEAN2,
            scale=config.MIXTURE_SD2
        )
        -
        norm.cdf(
            master_df["midpoint"],
            loc=config.MIXTURE_MEAN2,
            scale=config.MIXTURE_SD2
        )
    )
)

master_df["expected_n_young"] = config.N * p_young
master_df["expected_n_old"] = config.N * p_old
master_df["expected_n_total"] = (
    master_df["expected_n_young"] + master_df["expected_n_old"]
)

master_df["expected_n_eff"] = (
    2
    * master_df["expected_n_young"]
    * master_df["expected_n_old"]
    / master_df["expected_n_total"]
)

# Expected mean age in each half-window
young_mean_age = []
old_mean_age = []

for midpt in master_df["midpoint"]:
    young_mean_age.append(
        truncated_bimodal_mean(
            lower=midpt - config.WINDOW / 2,
            upper=midpt,
            prop=config.MIXTURE_PROP1,
            mu1=config.MIXTURE_MEAN1,
            sd1=config.MIXTURE_SD1,
            mu2=config.MIXTURE_MEAN2,
            sd2=config.MIXTURE_SD2
        )
    )

    old_mean_age.append(
        truncated_bimodal_mean(
            lower=midpt,
            upper=midpt + config.WINDOW / 2,
            prop=config.MIXTURE_PROP1,
            mu1=config.MIXTURE_MEAN1,
            sd1=config.MIXTURE_SD1,
            mu2=config.MIXTURE_MEAN2,
            sd2=config.MIXTURE_SD2
        )
    )

master_df["expected_young_mean_age"] = young_mean_age
master_df["expected_old_mean_age"] = old_mean_age
master_df["expected_abs_age_diff"] = np.abs(
    master_df["expected_old_mean_age"] - master_df["expected_young_mean_age"]
)


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

# ax.legend(
#     frameon=False,
#     fontsize=7,
#     handlelength=1.6,
#     borderaxespad=0.3,
#     loc="best"
# )

fig.tight_layout()
save_pdf_png(fig, plot_dir, "de_swan_all_curves_with_mean")
plt.show()


##### Plot 2: mean DE-SWAN curve overlaid with bimodal age density
x_age = np.linspace(
    plot_df["midpoint"].min(),
    plot_df["midpoint"].max(),
    500
)

bimodal_density = (
    config.MIXTURE_PROP1 * norm.pdf(
        x_age,
        loc=config.MIXTURE_MEAN1,
        scale=config.MIXTURE_SD1
    )
    +
    (1 - config.MIXTURE_PROP1) * norm.pdf(
        x_age,
        loc=config.MIXTURE_MEAN2,
        scale=config.MIXTURE_SD2
    )
)

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
    bimodal_density,
    linestyle="--",
    linewidth=1.8,
    label=rf"Age density"
)

ax2.set_ylabel("Age density")
ax2.tick_params(axis="both", width=0.8, length=3)

ax1.set_title("DE-SWAN hits overlaid with simulated age density")

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()


fig.tight_layout()
save_pdf_png(fig, plot_dir, "de_swan_mean_with_bimodal_density")
plt.show()


##### Plot 4: mean DE-SWAN curve overlaid with bimodal probability mass in each sliding window
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
    plot_df["midpoint"],
    plot_df["bimodal_window_prob"],
    linestyle="--",
    linewidth=1.8,
    label=rf"$P(a \in [m-{config.WINDOW/2}, m+{config.WINDOW/2}])$"
)

ax2.set_ylabel("Probability in sliding window")
ax2.tick_params(axis="both", width=0.8, length=3)

ax1.set_title("DE-SWAN hits overlaid with age-window probability")

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

fig.tight_layout()
save_pdf_png(fig, plot_dir, "de_swan_mean_with_bimodal_window_probability")
plt.show()


##### Plot 5: mean DE-SWAN curve overlaid with expected sample count per window
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
    plot_df["midpoint"],
    plot_df["expected_n_total"],
    linestyle="--",
    linewidth=1.8,
    label="Expected samples in window"
)

ax2.set_ylabel("Expected sample count per sliding window")
ax2.tick_params(axis="both", width=0.8, length=3)

ax1.set_title("DE-SWAN hits overlaid with expected window sample count")

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()

fig.tight_layout()
save_pdf_png(fig, plot_dir, "de_swan_mean_with_expected_window_count")
plt.show()


##### Plot 6: mean DE-SWAN curve overlaid with expected sample count per window
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
    plot_df["midpoint"],
    np.minimum(
    plot_df["expected_n_young"],
    plot_df["expected_n_old"]),
    linestyle="--",
    linewidth=1.8
)

ax2.set_ylabel("Expected samples in smallest group (young/old)")
ax2.tick_params(axis="both", width=0.8, length=3)

# ax1.set_title("DE-SWAN hits overlaid with expected window sample count")

lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.set_title("DE-SWAN hits with expected smallest group count")

fig.tight_layout()
save_pdf_png(fig, plot_dir, "de_swan_mean_with_smallest_expected_count")
plt.show()


