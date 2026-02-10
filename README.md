# Bayesian Ranking Algorithm for Olympic Medal Table

This repository implements a **Bayesian ranking model** for Olympic medal performance.  The method aims to estimate each nation’s *long-run per-capita medal probability*, producing rankings that better balances performance across large and small countries by applying **Bayesian shrinkage**.  

---

## 1. Overview  

Traditional medal tables reward population size. Per-capita tables correct this but exaggerate small-country fluctuations. This Bayesian framework models medal counts as **Poisson processes**.

The model is run using the `targets` pipeline in R for reproducibility and scalability across multiple Olympic Games datasets.  

---

## 2. Core Idea  

Each country’s medal-winning process is modeled through the probability that an individual wins exactly *i* medals, for *i = 1–4*. These probabilities are tied together through hierarchical priors that borrow strength across all countries, enabling shrinkage toward global averages.  

The key estimated quantity is the **posterior distribution of the per-capita medal rate**:  

$$
E(M_c/n_c) = p_{1,c} + 2p_{2,c} + 3p_{3,c} + 4p_{4,c}
$$ 

where $p_{i,c}$ is the probability that an individual from country $c$ wins exactly *i* medals.  

Posterior samples are **ranked within each iteration**, producing **posterior distributions of ranks** from which each country is ranked based on a central measure - here taken as the mean.

---

## 3. Priors and Model Variants  

The algorithm supports several prior specifications, each defined in `R/jags_model.R`:  

| Model | Prior on $p_c$ | Notes |
|--------|------------------|-------|
| **`beta`** | $p_c \sim \text{Beta}(a, b)$ with hyperpriors $a \sim U(0,1)$, $b \sim U(10^4,10^8)$ | Conditional probabilities shared across all countries |
| **`beta-2`** | Adds country-level variation: $X_1[i], X_2[i] \sim \text{Beta}(a_k,b_k)$ | Country specific conditional probability for double medallists |
| **`logit-normal`** | $\text{logit}(p_c) \sim N(\mu, \tau^{-1})$ | heavyier tail prior|
| **`mixture-beta`** | Mixture of multiple Beta components | Should capture multimodal medal-performance distributions - if they exist. |

All models assume the number of medal-winning athletes for a country follows Poisson processes

---

## 4. Pipeline Description  

The project uses the **{targets}** framework for automation. Key pipeline stages:  

| Target | Function | Purpose |
|---------|-----------|----------|
| `medalcounts` | Load CSV medal data | Input data per Olympic Games. See [here](https://github.com/Cormy1/Processing-of-Medal-Data-for-Olympics) for how this might be done |
| `model_datalist` | `make_model_list()` | Prepare model inputs |
| `jags_beta` | `jags_model()` | Select prior model type |
| `mcmc_beta` | `jags_run()` | Run Gibbs sampling in JAGS |
| `process_mcmc_beta` | `tidy_mcmc()` | Combine and format posterior chains |
| `posterior_ranks_beta` | `rank_mcmc()` | Compute posterior rank distributions |
| `probs_beta` | `post_prob_summary()` | Calculate summaries of posterior probability distributions |
| `ranks_beta` | `post_rank_rank()` | Compute posterior mean and credible intervals for distribution of ranks |
| `results.ranking` | `results()` | Merge observed and estimated metrics |
| `export_results` | Write CSV | Save final ranked table |

---

## 5. How to Run  

1. Place raw data in `0_data/0_raw/Medalcounts_<games>.csv`.  (I have included data in this repository for Tokyo, Rio, London, Beijing and Athens)
2. Edit the `games` vector to include desired Olympic editions.  eg. tokyo-2020
3. Run the pipeline:  
   ```r
   tar_make()
   ```
4. Output will be written as 
```
0_data/1_output/BayesianRanking_<games>.csv
```
---
## Dependancies
```
library(targets)
library(tarchetypes)
library(tidyverse)
library(countrycode)
library(rjags)
```
---
**Author’s Note:**
This pipeline is modular, transparent, and hopefully extendable to future Olympic cycles or alternative priors with minimal alteration. For teh unfamiliar, targets automates data–model–output flow, ensuring full reproducibility of Bayesian medal rankings. (To visualise the pipeline run `tar_visnetwork`)
