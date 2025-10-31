## parabiosis
**Probabilistic migration events drive transient tissue residency of lymphocytes during homeostasis**


## Overview

This repository contains Stan and R scripts and data for modelling and analysing lymphocyte migration dynamics, focusing on probabilistic migration events that result in transient tissue residency during homeostasis. It provides Markov Chain Monte Carlo (MCMC) simulation to estimate Markov chain model, simulation of cell behaviour based on the fitted Markov chain models and statistical analysis of lymphocyte migration patterns across different tissues under normal physiological conditions.

## Features

- Stan and R scripts implementing models of lymphocyte migration  
- Data sets used for model validation and analysis  
- R Scripts for statistical analysis and visualisation of migration behaviour  

## System Requirements
- **OS:** Linux (≥ Ubuntu 20.04), macOS (≥ 11), Windows (≥ 11).
- **R:** ≥ 4.5.1  
- **R packages:** `digest`, `expm`, `gridExtra`, `gtools`, `modeest`, `rstan`, `stringr`, `scales`, `magrittr`, `ggplot2`, `dplyr`, `readr`, `tidyr`, `forcats`, `patchwork`, `ggh4x`, `posterior`, `ggpubr`
- **Hardware:** ≥ 16 GB RAM recommended

## Installation Guide

1. **Clone the repository**:
```bash
git clone https://github.com/gergelits/parabiosis.git
cd parabiosis
```
2. **Then, open R and run**
```R
install.packages(c(
  "digest", "expm", "gridExtra", "gtools", "modeest", "rstan",
  "stringr", "scales", "magrittr", "ggplot2", "dplyr", "readr",
  "tidyr", "forcats", "patchwork", "ggh4x", "posterior", "ggpubr"
))
```
**Typical installation time:**  
~20–30 minutes on a normal desktop computer


## Usage
**Run the full model suite** (all cell type–tissue combinations) from the command line:
```bash
bash batch_ssub_args_cluster_m0003_i5004.sh
```

**Demo: Run a selected single model** (one cell type–tissue combination) from within R:
```R
source("batch_args_cluster.r")
```
Customise the script arguments as needed for your selected model.

**Expected output**
- Posterior parameter estimates saved as ```.csv```, MCMC chains saved as ```.rda``` objects in ```Analyzed/``` directory.
- Figures saved in ```Results/``` directory.

**Expected runtime for demo**  
~20 minutes on a "normal" desktop computer for one MCMC model (1 chain and 1300 iterations).
This includes Stan code compilation, MCMC simulation, simulation from the estimated Markov
chain model, and figure generation.

## License

See the LICENSE file in the repository for license details.
