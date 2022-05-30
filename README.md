# Paper repository
The repository contains code for reproducing the numerical results in the paper "Derandomized Knockoffs: Leveraging E-values for False Discovery Rate Control".

## Overview
We develop a method for derandomzing the model-X knockoffs procedure 
([Candes et al. 2018](https://candes.su.domains/publications/downloads/MX_Knockoffs.pdf))
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
- `/simulations` contains the R scripts to reproduce the simulation results.
  * `simulation_linear.R` reproduces the upper row of Figure 1.
  * `simulation_binom.R` reproduces the lower row of Figure 1.
  * `simulation_highdim.R` reproduces Figure D.1.
  * `simulation_linear_pfer.R` reproduces Figure D.2.
  * `simulation_linear_robust.R` reproduces Figure D.3.
  * `simulation_multi.R` reproduces Figure D.4.
  * `simulation_sideinfo.R` reproduces Figure D.5.
- `/real_data` contains the code to reproduce the results in Section 5.
- `/data` contains the preprocessed HIV datasets.
- `/utils` contains core functions to implement the methods.

## Usage
