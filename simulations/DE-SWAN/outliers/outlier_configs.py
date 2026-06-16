from datetime import datetime
import os
import sys

PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.append(PROJECT_ROOT)

N = 200 # number of samples
MOLS_NUM = 1500 # number of molecules (features)
N_SIMS = 5 # number of simulations

WINDOW = 20 # DE-SWAN window (younger = [midpoint - WINDOW/2, midpoint], older = [midpoint, midpoint + WINDOW/2])
MIN_AGE = 25 # midpoint age at which to start DE-SWAN algorithm
MAX_AGE = 75 # midpoint age at which to stop DE-SWAN algorithm

TRUE_SLOPE_SD = 0.1 # Standard Deviation of true slope values across molecules 
TRUE_INTERCEPT_SD = 0.1 # Standard Deviation of true intercept values across molecules 
NOISE_SD = 0.85 # Standard Deviation of noise term

UNIF_LOWER = 25 # lower bound of uniform distribution for age sampling
UNIF_UPPER = 75 # upper bound of uniform distribution for age sampling

DIST = "unif"  # distribution specification
NOISE_TREND = None

N_OUTLIERS = 5 # Number of samples with (potentially) extreme expression values
OUTLIER_MULTIPLIER = 5 # Multiplicative factor for the NOISE_SD for outliers

date_time_str = datetime.now().strftime("%m-%d-%y_%H-%M")
OUTDIR = os.path.join(PROJECT_ROOT, "outliers", "results",f"results_{date_time_str}")

