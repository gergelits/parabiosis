## parabiosis
**Probabilistic migration events drive transient tissue residency of lymphocytes during homeostasis**


## Overview

This repository contains Stan and R scripts and data for modelling and analysing lymphocyte migration dynamics, focusing on probabilistic migration events that result in transient tissue residency during homeostasis. It provides Markov Chain Monte Carlo (MCMC) simulation and statistical analysis of lymphocyte migration patterns across different tissues under normal physiological conditions.

## Features

- Stan and R scripts implementing models of lymphocyte migration  
- Data sets used for model validation and analysis  
- Scripts for statistical analysis and visualisation of migration behaviour  


## Usage

1. **Clone the repository**:
```bash
git clone https://github.com/gergelits/parabiosis.git
```
3. **Run the full model suite** (all cell type–tissue combinations) from the command line:
```bash
bash batch_ssub_args_cluster_m0003_i5004.sh
```
5. **Run a selected single model** (one cell type–tissue combination) from within R:
```R
source("batch_args_cluster.r")
```
Customise the script arguments as needed for your selected model.


## License

See the LICENSE file in the repository for license details.
