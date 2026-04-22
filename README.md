# Bayesian Stopping Rules for Machine Learning-Assisted Systematic Review Screening

**Author:** Ashley Razo  
**Institution:** Hertie School  
**Date:** April 2026

---

## Overview

This repository contains all code, data, and outputs associated with the thesis:

> *When to Stop: Bayesian Updating for Robust Stopping Rules in 
> Machine Learning-Assisted Systematic Review Screening*

The thesis develops and evaluates the first Bayesian stopping rule for 
machine learning-assisted systematic review screening. The method models 
the screening process as a biased urn process, jointly inferring the total 
number of relevant documents and the ranking bias parameter from the observed 
screening sequence. Stopping is triggered once the posterior probability of 
having achieved a target recall level exceeds a user-specified confidence 
threshold.

The evaluation proceeds in three stages:
1. Controlled Monte Carlo simulations under constant bias
2. Robustness analysis under time-varying bias scenarios
3. Real-world evaluation on ten completed systematic reviews via the PIK 
   computing cluster

---

## Repository Structure

---

## Simulations (`/simulations`)

All simulation code is written in R. The scripts should be run in the 
following order:

| Script | Description |
|--------|-------------|
| `Simulation_Code.R` | Core simulation engine: biased urn generator and Bayesian inference procedure |
| `simulation_monte_carlo.R` | Stage 1: Monte Carlo evaluation across four configurations and nine τ/γ combinations |
| `monte_carlo_heatmap.R` | Generates heatmap visualisations from Monte Carlo results |
| `simulation_with_bias.R` | Stage 2: Four bias scenarios (constant, linear decay, sudden drop, remaining relevance) |
| `simulating_change_in_omega.R` | Sensitivity analysis over decay strength λ with and without conservative adjustment |
| `decay_curve_with_and_without_adjustment.R` | Generates comparison plot of unadjusted vs adjusted method |

### Dependencies

```r
install.packages(c("BiasedUrn", "ggplot2", "dplyr", "parallel"))
```

### Key Parameters

- Baseline configuration: N = 10,000, prevalence = 5%, ω = 3
- Recall target: τ ∈ {0.90, 0.95, 0.99}
- Confidence level: γ ∈ {0.90, 0.95, 0.99}
- Monte Carlo runs: 100 per configuration

---

## Real-World Evaluation (`/pik-cluster`)

| File | Description |
|------|-------------|
| `bayesian_stopping_experiments.ipynb` | Jupyter notebook containing the Python implementation of the Bayesian stopping rule before scaling to the PIK cluster |
| `results/summary_results.csv` | Aggregated summary statistics across all reviews and confidence levels |
| `results/results_review_X.csv` | Per-run results for each of the ten systematic reviews (reviews 3--12) |
| `results/results_analysis.ipynb` | Combined analysis notebook aggregating all per-review CSV files into a single dataset and generating all figures used in the thesis (Figures 6, 7, and 8) |

### Result Files

Each `results_review_X.csv` file contains the following columns:

| Column | Description |
|--------|-------------|
| `review_id` | Review identifier |
| `run_id` | Screening run identifier |
| `N` | Total corpus size |
| `true_K` | True number of relevant documents |
| `confidence` | Confidence level γ |
| `stop_n_bayes` | Stopping point under Bayesian method |
| `true_recall_bayes` | Achieved recall under Bayesian method |
| `stop_n_baseline` | Stopping point under frequentist baseline |
| `true_recall_baseline` | Achieved recall under frequentist baseline |

### Dependencies
pip install numpy pandas matplotlib jupyter

---

## Figures

All figures used in the thesis are available in their respective directories:

**From simulations (`/simulations/figures`):**
- `heatmap_recall.png` — Median recall at stopping across configurations
- `heatmap_worksaved.png` — Median work saved across configurations  
- `heatmap_missrate.png` — Miss rate across configurations
- `plot_for_change_in_omega.png` — Recall under four bias scenarios
- `plot_decay_comparison.png` — Conservative adjustment sensitivity analysis

**From real-world evaluation (`/pik-cluster/results`):**
- `plot_miss_rate.png` — Miss rate vs confidence level (Bayesian vs baseline)
- `plot_work_saved.png` — Work saved vs confidence level (Bayesian vs baseline)
- `plot_additional_effort_hist.png` — Distribution of additional effort at γ = 0.95

---

## Data

The real-world screening data used in Stage 3 is sourced from:

> Callaghan, M. W., & Müller-Hansen, F. (2020). Statistical stopping criteria 
> for automated screening in systematic reviews. *Systematic Reviews*, 9(1), 273.
> https://doi.org/10.1186/s13643-020-01521-4

The dataset is publicly available on Zenodo:
https://zenodo.org/records/18164103


The dataset is not included in this repository due to its size. 
Please download it directly from the Zenodo link above and place 
it in a local `data/` folder in the root of this repository before 
running the PIK cluster notebook.

### Processed Results

The processed outputs from running the Bayesian stopping rule and 
frequentist baseline on the Zenodo dataset are included in this 
repository under `pik-cluster/results/`. These CSV files contain 
all stopping points, achieved recall, and work saved values used 
in the thesis results and are sufficient to reproduce all figures 
and tables without re-running the full pipeline.

---

## Thesis

The final submitted thesis PDF is available in `/thesis/Thesis_Draft.pdf`.

---

## Citation

If you use or build on this work, please cite:

Razo, A. (2026). When to Stop: Bayesian Updating for Robust Stopping Rules
in Machine Learning-Assisted Systematic Review Screening.
Master's Thesis, Hertie School.

---

## Contact

For questions about this repository, please contact Ashley Razo.
