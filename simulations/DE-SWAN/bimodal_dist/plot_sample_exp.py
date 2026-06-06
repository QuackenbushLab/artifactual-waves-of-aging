import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(PROJECT_ROOT)

from de_SWAN_helpers import generate_exp_data
import bimodal_dist_configs as config


def main():
    
    plot_DIR = config.RESULTS_DIR + "/plots"
    os.makedirs(plot_DIR, exist_ok=True)

    # Generate data
    exp_df = generate_exp_data(
        n=config.N,
        mols_num=config.MOLS_NUM,
        true_slope_sd=config.TRUE_SLOPE_SD,
        true_intercept_sd=config.TRUE_INTERCEPT_SD,
        noise_sd=config.NOISE_SD,
        dist=config.DIST,
        noise_trend=config.NOISE_TREND,
        mixture_prop=config.MIXTURE_PROP1,
        mixture_mean1=config.MIXTURE_MEAN1,
        mixture_sd1=config.MIXTURE_SD1,
        mixture_mean2=config.MIXTURE_MEAN2,
        mixture_sd2=config.MIXTURE_SD2
    )

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
        os.path.join(plot_DIR, "random_6_molecules_2x3.pdf"),
        bbox_inches="tight"
    )
    plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()