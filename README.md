# Paper repository
The repository contains code for reproducing the numerical results in the paper "Derandomized Knockoffs: Leveraging E-values for False Discovery Rate Control".

## Overview
We develop a method for derandomzing the model-X knockoffs procedure 
([Candès et al. 2018](https://candes.su.domains/publications/downloads/MX_Knockoffs.pdf))
with provable false discovery rate (FDR) control. The code in this directory 
reproduces simulation results in Section 4 and Section D of the Appendix,
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
  * `simulation_linear.R` reproduces the simulation results from the Gaussian linear model in Section 4.
  * `simulation_binom.R` reproduces the simulation results from the logistic mode in Section 4.
  * `simulation_highdim.R` reproduces the simulation results from the high-dimensional model in Appendix D1.
  * `simulation_original_low.R` reproduces the comparison with the original knockoffs procedure under a low-dimensional setting 
considered in [Candès et al. 2018](https://candes.su.domains/publications/downloads/MX_Knockoffs.pdf).
  * `simulation_original_high.R` reproduces the comparison with the original knockoffs procedure under a high-dimensional setting 
considered in [Candès et al. 2018](https://candes.su.domains/publications/downloads/MX_Knockoffs.pdf).
  * `simulation_pfer.R` reproduces the comparison with the PFER version in Appendix D2.
  * `simulation_robustness.R` reproduces the simulation with imperfect knowledge of X distribution in Appendix D3. 
  * `simulation_multienv.R` reproduces the results from the multi-environment simulation in Appendix D4.
  * `simulation_sideinfo.R` reproduces the the results from the setting with side informatino in Appendix D5.
- `/real_data` contains the R script `knockoffs.R` to reproduce the real data analysis in Section 5.
- `/data` contains the preprocessed HIV datasets.
- `/utils` contains the core functions to implement the methods.

## Usage
### A single run
To run `simulation_linear.R`, `simulation_binom.R`, `simulation_highdim.R`, 
`simulation_original_low.R`, `simulation_original_high.R`, `simulation_multienv.R` or `simulation_sideinfo.R`,
one needs to specify the signal amplitude and the random seed. For example,
if we want to run `simulation_linear.R` with `seed=1` and signal amplitude `A=4.5`,
we can run the following command in terminal:
```{r}
cd simulations
Rscript simulation_binom.R 1 4.5
```
To run `simulation_pfer.R`, we need to specify the random seed, 
the signal amplitude and the parameter v. For instance, if using 
`seed=1,A=4.5,v=2`, we can execute the following command in the terminal:
```{r}
cd simulations
Rscript simulation_pfer.R 1 4.5 2
```

To run `simulation_robustness.R`, we need to specify the random seed, 
the signal amplitude and the number of unlabelled samples. For instance, if using 
`seed=1,A=4.5,N=600`, we can execute the following command in terminal:
```{r}
cd simulations
Rscript simulation_robustness.R 1 4.5 600
```


### Multiple runs
The simulation results presented in the paper are averaged over multiple runs. 
The bash files in the folder `/batch` call the desired functions in the batch 
mode. For example, if we want to run `simulation_binom.R` with signal amplitudes
`A={16,19,22,25,28}` and `seed={1,2,...,100}`, we can run the following commands
in the terminal:
```{r}
cd batch
bash run_simulation_binom.sh
```
Note that it may take a long time to run all the repetitions locally.
It is recommended to run the batch mode on a server parallelly.


## Acknowledgement
The code to implement the [multi-environment knockoff filter (MEKF)](https://academic.oup.com/biomet/advance-article-abstract/doi/10.1093/biomet/asab055/6415825?redirectedFrom=fulltext&login=false)
is adapted from [https://github.com/lsn235711/MEKF_code](https://github.com/lsn235711/MEKF_code).

