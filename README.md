# Artifactual Waves of Aging

This repository contains code supporting _Artifactual "waves" of molecular aging arise when coupling LOESS with DE-SWAN_. 

<!-- ## USE: 
``` 
git clone ___
cd ___
virtualenv venv to create your new environment (called 'venv' here)
source venv/bin/activate to enter the virtual environment
pip install -r requirements.txt 
```  -->

## Contents

+ `r-notebooks` contains the scripts for Section 3: "Simulation example on a single gene", Section 4: "4	Waves of aging appear in null data and depend on the age distribution", and Figures 2, 3, and 4.
+ `ipop-data-analysis` contains the scripts for Section 5: "DE-SWAN alone does not detect 'waves of aging' in the iPOP data", Section 6: "Use of permutation testing to validate the procedure", and Figures 5, 6, and 7.
+ `utils` provides supporting utility functions for `r-notebooks`.

### DE-SWAN Simulations

Code for the DE-SWAN simulations can be found in `simulations/DE-SWAN`. In each folder within `DE-SWAN`, there are different folders containing code that modulates different aspects of the DGP: 

* `uniform_dist`: Sets age distribution to be uniform (Figure XX: Row 1)
* `normal_dist`: Sets age distribution to be normal (Figure XX: Row 2)
* `bimodal_dist`: Sets age distribution to be bimodal (specification 1) (Figure XX: Row 3)
* `bimodal_2_dist`: Sets age distribution to be bimodal (specification 2) (Figure XX: Row 4)
* MC TO-DO: Fill in the rest of these

Within each folder, there is a configuration file with parameters for each simulation setting (`XXX_configs.py`), a python executable file for one simulation (`run_one_sim.py`), and two visualization files: `summarize_de_swan.py` and `plot_sample_exp.py`. `summarize_de_swan.py` plots DE-SWAN cuves and additional analysis and `plot_sample_exp.py` shows 6 examples of molecular expression data that are simulated by the parameters. To run one simulation, navigate to the root folder (`artifactual-waves-of-aging`) and use

```python
simulations/DE-SWAN/normal_dist/run_one_sim.py {SIM_NUM}
```

where `{SIM_NUM}` is some integer value indexing the simulation. To run multiple simulations, calling this script can be looped. For large simulations, we recommend using external computing resources. 

`simulations/DE-SWAN/de_SWAN_helpers.py` contains helper functions used by all simulation settings.


  
## Authors

Code in this repository was written by Madeleine Carbonneau (mtcarbonneau) and Kate Shutta (katehoffshutta).

## Citation

Add citation here when preprint is up. 

