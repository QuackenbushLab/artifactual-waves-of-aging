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
INC_VAR_DIST_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))

if INC_VAR_DIST_DIR not in sys.path:
    sys.path.insert(0, INC_VAR_DIST_DIR)

import exp_inc_configs as config

##### Paths
results_dir = config.RESULTS_DIR

plot_dir = os.path.join(results_dir, "plots")
os.makedirs(plot_dir, exist_ok=True)

y_lower = 600
y_upper = 1600
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


##### Compute uniform-distribution diagnostics

# Probability that age falls inside the full DE-SWAN window:
# [midpoint - WINDOW, midpoint + WINDOW)
# master_df["uniform_window_prob"] = (
#     uniform.cdf(
#         master_df["midpoint"] + config.WINDOW,
#         loc=config.UNIF_LOWER,
#         scale=config.UNIF_UPPER - config.UNIF_LOWER
#     )
#     -
#     uniform.cdf(
#         master_df["midpoint"] - config.WINDOW,
#         loc=config.UNIF_LOWER,
#         scale=config.UNIF_UPPER - config.UNIF_LOWER
#     )
# )

# # Young-half probability: [midpoint - WINDOW, midpoint)
# p_young = (
#     uniform.cdf(
#         master_df["midpoint"],
#         loc=config.UNIF_LOWER,
#         scale=config.UNIF_UPPER - config.UNIF_LOWER
#     )
#     -
#     uniform.cdf(
#         master_df["midpoint"] - config.WINDOW,
#         loc=config.UNIF_LOWER,
#         scale=config.UNIF_UPPER - config.UNIF_LOWER
#     )
# )

# # Old-half probability: [midpoint, midpoint + WINDOW)
# p_old = (
#     uniform.cdf(
#         master_df["midpoint"] + config.WINDOW,
#         loc=config.UNIF_LOWER,
#         scale=config.UNIF_UPPER - config.UNIF_LOWER
#     )
#     -
#     uniform.cdf(
#         master_df["midpoint"],
#         loc=config.UNIF_LOWER,
#         scale=config.UNIF_UPPER - config.UNIF_LOWER
#     )
# )

# master_df["expected_n_young"] = config.N * p_young
# master_df["expected_n_old"] = config.N * p_old
# master_df["expected_n_total"] = (
#     master_df["expected_n_young"] + master_df["expected_n_old"]
# )

# master_df["expected_n_eff"] = (
#     2
#     * master_df["expected_n_young"]
#     * master_df["expected_n_old"]
#     / master_df["expected_n_total"]
# )

# # Expected mean age in each half-window
# young_mean_age = []
# old_mean_age = []

# for midpt in master_df["midpoint"]:
#     young_mean_age.append(
#         truncated_unif_mean(
#             lower=midpt - config.WINDOW,
#             upper=midpt,
#             a=config.UNIF_LOWER,
#             b=config.UNIF_UPPER
#         )
#     )

#     old_mean_age.append(
#         truncated_unif_mean(
#             lower=midpt,
#             upper=midpt + config.WINDOW,
#             a=config.UNIF_LOWER,
#             b=config.UNIF_UPPER
#         )
#     )

# master_df["expected_young_mean_age"] = young_mean_age
# master_df["expected_old_mean_age"] = old_mean_age
# master_df["expected_abs_age_diff"] = np.abs(
#     master_df["expected_old_mean_age"] - master_df["expected_young_mean_age"]
# )


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

# noise_sd_by_age = config.NOISE_SD * (config.MIN_NOISE_MULTIPLIER + x_age)
age_scaled = (x_age - x_age.min()) / (x_age.max() - x_age.min())

exp_scaled = (np.exp(config.EXP_STRENGTH * age_scaled) - 1) / (
            np.exp(config.EXP_STRENGTH) - 1
        )

noise_sd_by_age = config.NOISE_SD  * (
            config.MIN_NOISE_MULTIPLIER
            + (config.MAX_NOISE_MULTIPLIER - config.MIN_NOISE_MULTIPLIER) * exp_scaled
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
ax.set_xlim(x_lower, x_upper)
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


# ##### Plot 4: mean DE-SWAN curve overlaid with uniform probability mass in each sliding window
# fig, ax1 = plt.subplots(figsize=(6.0, 3.6))

# ax1.plot(
#     plot_df["midpoint"],
#     plot_df["mean_n_sig"],
#     linewidth=2.2,
#     label="Mean DE-SWAN hits"
# )

# ax1.fill_between(
#     plot_df["midpoint"],
#     plot_df["mean_n_sig"] - plot_df["sd_n_sig"],
#     plot_df["mean_n_sig"] + plot_df["sd_n_sig"],
#     alpha=0.18,
#     label="Mean ± 1 SD"
# )

# ax1.set_xlabel("Age midpoint")
# ax1.set_ylabel("Significant molecules")
# ax1.tick_params(axis="both", width=0.8, length=3)

# ax2 = ax1.twinx()

# ax2.plot(
#     plot_df["midpoint"],
#     plot_df["uniform_window_prob"],
#     linestyle="--",
#     linewidth=1.8,
#     label=rf"$P(a \in [m-{config.WINDOW}, m+{config.WINDOW}])$"
# )

# ax2.set_ylabel("Probability in sliding window")
# ax2.tick_params(axis="both", width=0.8, length=3)

# ax1.set_title("DE-SWAN hits overlaid with age-window probability")

# lines1, labels1 = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()

# # ax1.legend(
# #     lines1 + lines2,
# #     labels1 + labels2,
# #     loc="upper right",
# #     frameon=False,
# #     fontsize=7,
# #     handlelength=1.6,
# #     borderaxespad=0.3
# # )

# fig.tight_layout()
# save_pdf_png(fig, plot_dir, "de_swan_mean_with_uniform_window_probability")
# plt.show()


# ##### Plot 5: mean DE-SWAN curve overlaid with expected sample count per window
# fig, ax1 = plt.subplots(figsize=(6.0, 3.6))

# ax1.plot(
#     plot_df["midpoint"],
#     plot_df["mean_n_sig"],
#     linewidth=2.2,
#     label="Mean DE-SWAN hits"
# )

# ax1.fill_between(
#     plot_df["midpoint"],
#     plot_df["mean_n_sig"] - plot_df["sd_n_sig"],
#     plot_df["mean_n_sig"] + plot_df["sd_n_sig"],
#     alpha=0.18,
#     label="Mean ± 1 SD"
# )

# ax1.set_xlabel("Age midpoint")
# ax1.set_ylabel("Significant molecules")
# ax1.tick_params(axis="both", width=0.8, length=3)

# ax2 = ax1.twinx()

# ax2.plot(
#     plot_df["midpoint"],
#     plot_df["expected_n_total"],
#     linestyle="--",
#     linewidth=1.8,
#     label="Expected samples in window"
# )

# ax2.set_ylabel("Expected sample count per sliding window")
# ax2.tick_params(axis="both", width=0.8, length=3)

# ax1.set_title("DE-SWAN hits overlaid with expected window sample count")

# lines1, labels1 = ax1.get_legend_handles_labels()
# lines2, labels2 = ax2.get_legend_handles_labels()

# # ax1.legend(
# #     lines1 + lines2,
# #     labels1 + labels2,
# #     loc="upper right",
# #     frameon=False,
# #     fontsize=7,
# #     handlelength=1.6,
# #     borderaxespad=0.3
# # )

# fig.tight_layout()
# save_pdf_png(fig, plot_dir, "de_swan_mean_with_expected_window_count")
# plt.show()


# # ##### Plot 6: mean DE-SWAN curve overlaid with expected effective sample size
# # fig, ax1 = plt.subplots(figsize=(5.5, 3.2))

# # ax1.plot(
# #     master_df["midpoint"],
# #     master_df["mean_n_sig"],
# #     linewidth=2.2,
# #     label="Mean DE-SWAN hits"
# # )

# # ax1.fill_between(
# #     master_df["midpoint"],
# #     master_df["mean_n_sig"] - master_df["sd_n_sig"],
# #     master_df["mean_n_sig"] + master_df["sd_n_sig"],
# #     alpha=0.18,
# #     label="Mean ± 1 SD"
# # )

# # ax1.set_xlabel("Age midpoint")
# # ax1.set_ylabel("Significant molecules")
# # ax1.tick_params(axis="both", width=0.8, length=3)

# # ax2 = ax1.twinx()

# # ax2.plot(
# #     master_df["midpoint"],
# #     master_df["expected_n_eff"],
# #     linestyle="--",
# #     linewidth=1.8,
# #     label="Expected effective sample size"
# # )

# # ax2.set_ylabel("Expected effective sample size")
# # ax2.tick_params(axis="both", width=0.8, length=3)

# # ax1.set_title("DE-SWAN hits overlaid with effective sample size")

# # lines1, labels1 = ax1.get_legend_handles_labels()
# # lines2, labels2 = ax2.get_legend_handles_labels()

# # ax1.legend(
# #     lines1 + lines2,
# #     labels1 + labels2,
# #     loc="upper right",
# #     frameon=False,
# #     fontsize=7,
# #     handlelength=1.6,
# #     borderaxespad=0.3
# # )

# # fig.tight_layout()
# # save_pdf_png(fig, plot_dir, "de_swan_mean_with_expected_effective_sample_size")
# # plt.show()


# # ##### Plot 7: mean DE-SWAN curve overlaid with expected absolute young-old age difference
# # fig, ax1 = plt.subplots(figsize=(5.5, 3.2))

# # ax1.plot(
# #     master_df["midpoint"],
# #     master_df["mean_n_sig"],
# #     linewidth=2.2,
# #     label="Mean DE-SWAN hits"
# # )

# # ax1.fill_between(
# #     master_df["midpoint"],
# #     master_df["mean_n_sig"] - master_df["sd_n_sig"],
# #     master_df["mean_n_sig"] + master_df["sd_n_sig"],
# #     alpha=0.18,
# #     label="Mean ± 1 SD"
# # )

# # ax1.set_xlabel("Age midpoint")
# # ax1.set_ylabel("Significant molecules")
# # ax1.tick_params(axis="both", width=0.8, length=3)

# # ax2 = ax1.twinx()

# # ax2.plot(
# #     master_df["midpoint"],
# #     master_df["expected_abs_age_diff"],
# #     linestyle="--",
# #     linewidth=1.8,
# #     label="Expected young-old age difference"
# # )

# # ax2.set_ylabel("Expected age difference")
# # ax2.tick_params(axis="both", width=0.8, length=3)

# # ax1.set_title("DE-SWAN hits overlaid with expected age separation")

# # lines1, labels1 = ax1.get_legend_handles_labels()
# # lines2, labels2 = ax2.get_legend_handles_labels()

# # ax1.legend(
# #     lines1 + lines2,
# #     labels1 + labels2,
# #     loc="upper right",
# #     frameon=False,
# #     fontsize=7,
# #     handlelength=1.6,
# #     borderaxespad=0.3
# # )

# # fig.tight_layout()
# # save_pdf_png(fig, plot_dir, "de_swan_mean_with_expected_age_difference")
# # plt.show()


# # ##### Save final summary with mean, SD, and diagnostic columns
# # summary_out_file = os.path.join(results_dir, "de_swan_summary_n=1000_p=500.csv")
# # master_df.to_csv(summary_out_file, index=False)

# # print(f"Saved summary with mean, SD, and diagnostics to {summary_out_file}")
# # print(f"Saved plots to {plot_dir}")