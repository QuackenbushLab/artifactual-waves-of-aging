N = 1000 # number of samples
MOLS_NUM = 2000 # number of molecules (features)
N_SIMS = 100 # Note: to run N_SIMS simulations, run run_one_sim.py script N_SIMS times with sim number as an argument 
             # (e.g. python run_one_sim.py 1, python run_one_sim.py 2, ...,  python run_one_sim.py {N_SIMS})

WINDOW = 20 # DE-SWAN window (younger = [midpoint - WINDOW/2, midpoint], older = [midpoint, midpoint + WINDOW/2])
MIN_AGE = 25 # midpoint age at which to start DE-SWAN algorithm
MAX_AGE = 75 # midpoint age at which to stop DE-SWAN algorithm

TRUE_SLOPE_SD = 0.1 # Standard Deviation of true slope values across molecules 
TRUE_INTERCEPT_SD = 0.1 # Standard Deviation of true intercept values across molecules 
NOISE_SD = 0.85 # Standard Deviation of noise term

UNIF_LOWER = 25 # lower bound of uniform distribution for age sampling
UNIF_UPPER = 75 # upper bound of uniform distribution for age sampling

DIST = "unif"  # distribution specification
NOISE_TREND = "exp_increasing"
EXP_STRENGTH = 10.0 # exponential power for scaling age (see simulations/DE-SWAN/de_SWAN_helpers.py) 
MIN_NOISE_MULTIPLIER = 2.0 # multiplicative factor for noise at minimum variance point 
MAX_NOISE_MULTIPLIER = 8.0# multiplicative factor for noise at maximum variance point 


RESULTS_DIR = "simulations/DE-SWAN/exp_inc/results/results_06-10-26" # folder for result
