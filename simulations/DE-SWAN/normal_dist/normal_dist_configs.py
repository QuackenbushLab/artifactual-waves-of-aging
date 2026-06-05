N = 1000
MOLS_NUM = 2000
N_SIMS = 100 # Note: to run N sims, run this script N times with sim number as an argument 
             # (e.g. python run_one_sim.py 1, python run_one_sim.py 2, etc.)

WINDOW = 20 # 10 years for younger group and 10 years for older group
MIN_AGE = 25 
MAX_AGE = 75

TRUE_SLOPE_SD = 0.1
TRUE_INTERCEPT_SD = 0.1
NOISE_SD = 0.85

DIST = "gaussian" # normal
NORM_MEAN = 50
NORM_SD = 7
NOISE_TREND = None

RESULTS_DIR = "simulations/DE-SWAN/normal_dist/results/results_06-05-26"
