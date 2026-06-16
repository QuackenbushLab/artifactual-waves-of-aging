N = 1000 # number of samples
MOLS_NUM = 2000 # number of molecules (features)
N_SIMS = 100 # Note: to run N_SIMS simulations, run run_one_sim.py script N_SIMS times with sim number as an argument 
             # (e.g. python run_one_sim.py 1, python run_one_sim.py 2, etc.)


WINDOW = 20 # DE-SWAN window (younger = [midpoint - WINDOW/2, midpoint], older = [midpoint, midpoint + WINDOW/2])
MIN_AGE = 25 # midpoint age at which to start DE-SWAN algorithm
MAX_AGE = 75 # midpoint age at which to stop DE-SWAN algorithm

TRUE_SLOPE_SD = 0.1 # Standard Deviation of true slope values across molecules 
TRUE_INTERCEPT_SD = 0.1 # Standard Deviation of true intercept values across molecules 
NOISE_SD = 0.85 # Standard Deviation of noise term

DIST = "gaussian" # distribution specification 
NORM_MEAN = 50 # mean of first normal distribution in bimodal mixture
NORM_SD = 7 # standard deviation of first normal distribution in bimodal mixture
NOISE_TREND = None # keep variance constant

RESULTS_DIR = "simulations/DE-SWAN/normal_dist/results/results_06-05-26" # folder for result
