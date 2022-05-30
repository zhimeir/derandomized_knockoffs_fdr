# Paper repository
The repository contains code for reproducing the numerical results in the paper "Derandomized Knockoffs: Leveraging E-values for False Discovery Rate Control".

## Overview
We develop a method for derandomzing the model-X knockoffs procedure 
([Candès et al. 2018](https://candes.su.domains/publications/downloads/MX_Knockoffs.pdf))
with provable false discovery rate (FDR) control. The code in this directory 
contains the code to reproduce simulation results in Section 4 and the Appendix D,
and the data analysis of the HIV dataset from the 
[HIV Drug Resistance Database](https://hivdb.stanford.edu/download/GenoPhenoDatasets/PI_DataSet.txt).

All the scripts are written in [R](https://www.r-project.org/) (version 4.0.3).
To run the scripts in this repository, the following packages need to 
be installed: [tidyverse](https://www.tidyverse.org/), 
[glmnet](https://cran.r-project.org/web/packages/glmnet/index.html),
[knockoff](https://cran.r-project.org/web/packages/knockoff/index.html), 
[gam](https://cran.r-project.org/web/packages/gam/index.html),
[adaptiveKnockoff](https://github.com/zhimeir/adaptiveKnockoffs).

## Folders
- `/simulations` contains the executable R scripts to reproduce the simulation results.
  * `simulation_linear.R` generates (one run of) the results corresponding to the Gaussian linear model simulation setup in Section 4.
  * `simulation_binom.R` generates (one run of) the results corresponding to the logistic model simulation setup in Section 4.
  * `simulation_highdim.R` generates (one run of) the results corresponding to the high-dimensional simulation setup in Appendix D1.
  * `simulation_linear_pfer.R` generates (one run of) the results corresponding to comparison with the PFER version in Appendix D2.
  * `simulation_linear_robust.R` generates (one run of) the results corresponding to simulation with imperfect knowledge of X distribution in Appendix D3. 
  * `simulation_multi.R` generates (one run of) the results corresponding to the multi-environment simulation setup in Appendix D4.
  * `simulation_sideinfo.R` generates (one run of) the results corresponding to setting with side informatino in Appendix D5.
- `/real_data` contains the code to reproduce the results in Section 5.
- `/data` contains the preprocessed HIV datasets.
- `/utils` contains core functions to implement the methods.

## Usage
### A single run
To run `simulation_linear.R`, `simulation_binom.R`, `simulation_highdim.R`, `simulation_multi.R` or `simulation_sideinfo.R`,
one needs to specify the signal amplitude and the random seed. For example,
if we want to run `simulation_binom.R` with random seed 1 and signal amplitude 19,
we can run the following command in terminal:
```{r}
cd simulations
Rscript simulation_binom.R 1 19
```

### Multiple runs



