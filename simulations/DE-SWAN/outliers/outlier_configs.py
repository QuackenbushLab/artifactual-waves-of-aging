from datetime import datetime
import os
import sys

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.append(PROJECT_ROOT)

N = 200
MOLS_NUM = 1500
N_SIMS = 50

WINDOW = 20
MIN_AGE = 25
MAX_AGE = 75

TRUE_SLOPE_SD = 0.1
TRUE_INTERCEPT_SD = 0.1
NOISE_SD = 0.85

DIST = "unif"
NOISE_TREND = None

N_OUTLIERS = 5
OUTLIER_MULTIPLIER = 5

date_time_str = datetime.now().strftime("%m-%d-%y_%H-%M")
OUTDIR = os.path.join(PROJECT_ROOT, "outliers", "results",f"results_{date_time_str}")

