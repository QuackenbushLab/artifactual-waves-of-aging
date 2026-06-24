# Artifactual Waves of Aging

This repository contains code supporting _Artifactual "waves" of molecular aging arise when coupling LOESS with DE-SWAN_. 

## Contents

+ `r-notebooks` contains the scripts for _Artifact 1: LOESS can induce clustered trajectories in null data_ and some components of _Artifact 2: LOESS+DE-SWAN can elicit artifactual waves of aging that are not supported by null simulations and permutation testing_.
+ `ipop-data-analysis` contains the scripts for the remaining components of _Artifact 2_ and for _Artifact 3: DE-SWAN identifies non-linear patterns of aging in linear data due to age-associated differences in statistical power_. 
+ `utils` provides supporting utility functions for `r-notebooks`.
+ `simulations/DE-SWAN` contains code for DE-SWAN simulations shown in section _Artifact 3: DE-SWAN identifies non-linear patterns of aging in linear data due to age-associated differences in statistical power_.

## Getting Started:

Navigate to your project directory
```bash
cd path/to/your/projects
```
Clone the repository
```bash
git clone https://github.com/QuackenbushLab/artifactual-waves-of-aging.git
``` 
Navigate to the repository
```bash
cd artifactual-waves-of-aging
```

### For Python code:

Set up a virtual environment:

Mac:
```bash
python3 -m venv .venv
```
Windows: 
```cmd
python -m venv .venv
```
Activate the environment:

Mac:
```bash
source .venv/bin/activate
```
Windows: 
```cmd
.venv\Scripts\activate.bat
```
Install dependencies:
```bash
pip install -r requirements.txt
```

### For R code

Required packages are listed at the top of each notebook in the `r-notebooks` directory. We are in the process of preparing an environment file to make this process easier. 

### Simulation studies on random noise

RMarkdown files in the `r-notebooks` directory can be used to reproduce the results of the simulations that identify patterns in random noise. This directory also includes the script to plot the literature search results.

### iPOP Data Analysis

Replication of the transcriptomic analysis in _Nonlinear dynamics of multi-omics profiles during human aging_ can be found in the Jupyter Notebooks here (`ipop_data_analysis/transcriptomics_replication.ipynb`). Instructions for downloading data and reproducing the analysis can be found in the notebook itself. When downloaded, data should be stored in `ipop_data_analysis/nonlinear_aging_data`. Due to the amount of time permutation testing can take, previous results have been stored in `ipop_data_analysis/permutation_rslts`. Lastly, `ipop_data_analysis/loess.R` contains code from the original analysis [repository](https://github.com/jaspershen-lab/ipop_aging) to run LOESS-fitting.

We have additionally included the the proteomic pipeline (`ipop_data_analysis/proteomics_replication.ipynb`).

### DE-SWAN Simulations

Code for the DE-SWAN simulations can be found in `simulations/DE-SWAN`. In each folder within `DE-SWAN`, there are different folders containing code that modulates different aspects of the DGP: 

* `uniform_dist`: Sets age distribution to be uniform (Figure 8: Row 1)
* `normal_dist`: Sets age distribution to be normal (Figure 8: Row 2)
* `bimodal_dist`: Sets age distribution to be bimodal (specification 1) (Figure 8: Row 3)
* `bimodal_2_dist`: Sets age distribution to be bimodal (specification 2) (Figure 8: Row 4)
* `middle_var`: Sets variance to increase and then decrease (Figure 10: Row 1)
* `ends_var`: Sets variance to decrease towards the center of the distribution and then increase (Figure 10: Row 2)
* `exp_inc`: Sets variance to increase nonlinearly (Figure 10: Row 3)
* `exp_dec`: Sets variance to increase nonlinearly (Figure 10: Row 3)
* `outliers`: Simulates outliers at the center of the age distribution (Figure S3)

Within each folder, there is a configuration file with parameters for each simulation setting (`XXX_configs.py`), a python executable file for one simulation (`run_one_sim.py`), and two visualization files: `summarize_de_swan.py` and `plot_sample_exp.py`. `summarize_de_swan.py` plots DE-SWAN cuves and additional analysis and `plot_sample_exp.py` shows 6 examples of molecular expression data that are simulated by the parameters. `simulations/DE-SWAN/de_SWAN_helpers.py` contains helper functions used by all simulation settings.

In general, simulations should be run according to the following procedure:
1. Set paramters including the results directory. Recommended format is `"simulations/DE-SWAN/bimodal_2_dist/results/results_{DATE}"`
2. Visualize the data generated with those parameters using 
```bash
python simulations/DE-SWAN/bimodal_2_dist/plot_sample_exp.py
```
3. Run one simulation. From the root folder (`artifactual-waves-of-aging`) and run 
```bash
python simulations/DE-SWAN/normal_dist/run_one_sim.py {SIM_NUM}
```
where `{SIM_NUM}` is some integer value indexing the simulation. 

4. Run multiple simulations, by looping through the script in 3. For large simulations (>15), we recommend using external computing resources. 

5. Visualize results by running 
```bash
python simulations/DE-SWAN/bimodal_2_dist/summarize_de_swan.py
```
Note that, depending on parameters, parameters such as `y_lower`,`y_upper`, `x_lower`, `x_upper`, `fig_width`, and `fig_height` might need to be changed.

In the case of the outliers simulation, all simulations can be run together as they are typically smaller. To run, use 

```bash 
python simulations/DE-SWAN/outliers/run_all_sims.py
```

The outliers simulations have been structured differently. To run the outliers simulations and generate all plots, the only command needed is

```bash
python simulations/DE-SWAN/outliers/run_all_sims.py
```
All outlier simulation results will be stored in a folder with time/date stamp.

  
## Authors

Code in this repository was written by Madeleine Carbonneau (mtcarbonneau) and Kate Shutta (katehoffshutta).

## Citation

Add citation here when preprint is up. 

