import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(PROJECT_ROOT)

from de_SWAN_helpers import generate_exp_data, run_all_de_swan
import bimodal_dist_configs as config


def main():
    sim = int(sys.argv[1])

    os.makedirs(config.RESULTS_DIR, exist_ok=True)

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


    # Run DE-SWAN
    de_swan_df = run_all_de_swan(
        mols_num=config.MOLS_NUM,
        exp_df=exp_df,
        window=config.WINDOW,
        age_col="age",
        min_age=config.MIN_AGE,
        max_age=config.MAX_AGE,
    )

    de_swan_df = de_swan_df[["midpoint", "n_sig"]].copy()
    de_swan_df = de_swan_df.rename(columns={"n_sig": f"sim_{sim}_n_sig"})

    out_file = os.path.join(config.RESULTS_DIR, f"de_swan_sim_{sim}_n={config.N}_p={config.MOLS_NUM}.csv")
    de_swan_df.to_csv(out_file, index=False)

    print(f"Finished simulation {sim}. Saved to {out_file}")


if __name__ == "__main__":
    main()