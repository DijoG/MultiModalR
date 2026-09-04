# MultiModalR 🏔️ 📊 🔬

[![R](https://img.shields.io/badge/R-≥4.0-blue?style=for-the-badge&logo=r)](https://www.r-project.org/)
[![C++](https://img.shields.io/badge/C++-RcppArmadillo-blue?style=for-the-badge&logo=cplusplus)](https://isocpp.org/)
[![Bayesian](https://img.shields.io/badge/Bayesian-MCMC-green?style=for-the-badge)](https://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo)

<img src="https://latex.codecogs.com/svg.latex?\color{white}w_1\mathcal{N}(\mu_1,\sigma_1^2)+w_2\mathcal{N}(\mu_2,\sigma_2^2)+w_3\mathcal{N}(\mu_3,\sigma_3^2)" height="30" align="center">


**MultiModalR** performs Bayesian mixture modeling for multimodal data. It detects subpopulations and assigns probabilistic memberships using two advanced Markov Chain Monte Carlo (MCMC) algorithms implemented in optimized C++:

1. **Metropolis-Hastings within Gibbs Sampler for Gaussian Mixture Models** - Fast and robust
2. **Dirichlet-Multinomial (collapsed Gibbs)** - Slower and rigorously robust

## 🎯 Features
- 🚀 **Dual MCMC algorithms**: Choose between Metropolis-Hastings (speed) or Dirichlet-Multinomial (robustness)
- 🔍 **Enhanced Mode Detection**: Height-aware peak detection with four bandwidth methods
- 📊 **Bayesian Probability Assignment**: Soft assignment with probability estimates
- 🎪 **Subpopulation Detection**: Automatic detection of multimodal components
- ⚡ **Parallel Processing**: Built-in multi-core computation for stratified models
- ✅ **Validation Tools**: Built-in plotting and validation functions
- 🏗️ **Hierarchical Mixture Model**: Unified model with category-specific parameters and shared components
- 🔬 **MCMC Diagnostics**: Effective Sample Size (ESS) and Gelman-Rubin convergence checks
- ⚙️ **Optimized C++ Core**: Blazing fast MCMC sampling with RcppArmadillo

## 📋 Prerequisites

### R Package Dependencies
`Rcpp`, `RcppArmadillo`, `dplyr`, `furrr`, `future`

### System Requirements
- R ≥ 4.0 
- Multiple CPU cores for parallel processing

## 💾 Installation
```r
devtools::install_github("DijoG/MultiModalR")
```

## 🚀 Quick Start Example

### 1. Stratified Mixture Model (fuss_PARALLEL_mcmc)
When to use: Categories are independent and should not share information. Each category gets its own separate mixture model.

```r
library(MultiModalR)

# Load data
df <- MultiModalR::multimodal_dummy

# Run stratified analysis (parallel)
results <- fuss_PARALLEL_mcmc(
  data = df,
  varCLASS = "Category",
  varY = "Value",
  varID = "ID"
)

# View results summary
summary(results)
```
### 2. Hierarchical Mixture Model (fuss_COVARIATE_mcmc)
When to use: Categories share the same underlying components (same K) but may have different means, variances, and weights. Borrows information across categories, making it more stable for small categories.

```r
# Run hierarchical model
result <- fuss_COVARIATE_mcmc(
  data = df,
  varY = "Value",
  varCLASS = "Category",
  K = 3,                    # Auto-detected if NULL
  out_dir = "output",       # Auto-writes CSV files
  n_iter = 10000,
  burnin = 2000
)

# View results
summary(result)
plot(result)
head(result$prob_matrix)

# Check convergence
MultiModalR::check_convergence(result)
```

## ⚙️ Parameters
### Sratified Model (fuss_PARALLEL_mcmc)
```r
MultiModalR::fuss_PARALLEL_mcmc(
  data = df,                  # 📦 -> required
  varCLASS = "Category",      # 🏷️ -> required
  varY = "Value",             # 📈 -> required
  varID = "ID",               # 🆔 -> required
  method = "sj-dpi",          # 📏 /default
  within = 1,                 # 🎯 /default
  maxNGROUP = 5,              # 🔢 /default
  out_dir = ".../output",     # 💾 -> optional 
  n_workers = 3,              # ⚡ /default
  n_iter = NULL,              # 🔄 /default
  burnin = NULL,              # 🔥 /default
  proposal_sd = 0.15,         # 📊 /default
  sj_adjust = 0.5,            # ⚖️ /default
  mcmc_method = "metropolis", # 🧮 /default
  dirichlet_alpha = 2.0       # 🎲 /default
)
```
### Hierarchical Model (fuss_COVARIATE_mcmc)
```r
MultiModalR::fuss_COVARIATE_mcmc(
  data = df,                  # 📦 -> required
  varY = "Value",             # 📈 -> required
  varCLASS = "Category",      # 🏷️ -> required
  varID = "ID",               # 🆔 -> optional
  K = NULL,                   # 🔢 auto-detected
  out_dir = NULL,             # 💾 -> optional (auto-writes CSV)
  n_iter = 10000,             # 🔄 /default
  burnin = 2000,              # 🔥 /default
  proposal_sd = 0.15,         # 📊 /default
  adaptive = TRUE,            # 🧠 /default
  alpha0 = 3.0,               # 📐 /default
  beta0 = 2.0,                # 📐 /default
  alpha_dirichlet = 5.0,      # 🎲 /default
  method = "sj-dpi",          # 📏 /default
  sj_adjust = 0.5,            # ⚖️ /default
  within = 1.0,               # 🎯 /default
  seed = 123                  # 🎲 /default
)
```
## 📚 Detailed Examples

### Data
```r
library(MultiModalR)
library(ggplot2)
library(dplyr)

# Load the built-in dataset
df <- MultiModalR::multimodal_dummy

# View the data structure
head(df)
str(df)
```
### Data Visualization
```r
library(ggplot2)

# Plot 01 ~ subpopulations/subgroups not shown
ggplot(df, aes(x = Value)) +
  geom_density(color = NA, fill = "grey98", adjust = .8) +
  facet_wrap(~Category) +
  theme_dark() +
  labs(title = "Multimodal Data ~ Density", 
       x = "Value", y = "Density") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme(legend.position = "top",
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5))
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/MMR/MMR_001.png" width="550">

```r
# Plot 02 ~ subgroups shown
ggplot(df, aes(x = Value, fill = Subpopulation)) +
  geom_density(alpha = 0.5, color = NA) +
  scale_fill_manual(values = c("firebrick2", "forestgreen", "cyan3"), 
                     name = "Subgroups") +
  facet_wrap(~Category) +
  theme_dark() +
  labs(title = "Multimodal Data ~ Density with Subgroups", 
       x = "Value", y = "Density") +
  scale_y_continuous(expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  theme(legend.position = "top",
        legend.key = element_rect(fill = "transparent", color = NA),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = .5)) +
  guides(fill = guide_legend(override.aes = list(alpha = .6)))
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/MMR/MMR_002.png" width="550">

### Running Stratified Analysis 
```r
# Run stratified model with parallel processing
cores <- 3

results <- fuss_PARALLEL_mcmc(
  data = df,
  varCLASS = "Category",
  varY = "Value",
  varID = "ID",
  out_dir = "D:/test/MMR_strat",
  n_workers = cores
)

# View results
head(results)
```
### Running Hierarchical Analysis
```r
# Run hierarchical model with auto CSV output
result <- fuss_COVARIATE_mcmc(
  data = df,
  varY = "Value",
  varCLASS = "Category",
  K = 3,
  out_dir = "D:/test/MMR_hier",
  n_iter = 10000,
  burnin = 2000,
  proposal_sd = 0.12
)

# View results
summary(result)
plot(result)
head(result$prob_matrix)

# Check convergence
ess <- check_convergence(result)
```
### Output 

Both functions generate:
  - **Data CSV** files: Original data with assigned subgroups and probabilities
  
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/MMR/MMR_csv.png" width="550">

A **Data CSV** file consists of the following fields (maxNGROUP = 5):
  - `y`: Original/observed value
  - `Group`: Original/observed subgroup
  - `Group_1`: Predicted belonging probability to
  - `Group_2`: Predicted belonging probability to
  - `Group_3`: Predicted belonging probability to
  - `Group_4`: Predicted belonging probability to
  - `Group_5`: Predicted belonging probability to
  - `Assigned_Group`: Assigned/predicted subgroup
  - `Min_Assigned`: Minimum value of the assigned/predicted range
  - `Max_Assigned`: Maximum value of the assigned/predicted range
  - `Mean_Assigned`: Mean value of the assigned/predicted range
  - `Mode_Assigned`: Mode of the assigned/predicted range
  - `Main_Class`: Category/main group/class

### Validation Visualization
The `plot_VALIDATION()` function creates a density plot with jittered points colored by assigned group:
```r
# For stratified model
plot_VALIDATION(
  "D:/test/MMR_strat",
  df,
  subpop_col = "Subpopulation",
  value_col = "Value",
  id_col = "ID"
)
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/MMR/MMR_strat.png" width="550">

```r
# For hierarchical model
plot_VALIDATION(
  "D:/test/MMR_hier",
  df,
  subpop_col = "Subpopulation",
  value_col = "Value",
  id_col = "ID"
)
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/MMR/MMR_hier.png" width="550">

**Validation results show accurate subgroup assignments across categories.**

## 📊 Performance Comparison

### Accuracy Metrics
| **Aspect** | **Stratified** | **Hierarchical** | **Winner** |
|:---|:---:|:---:|:---:|
| **Accuracy** | **92.6%** | 90.7% | 🏆 Stratified |
| **Speed** | **12.59s** | 104.11s | 🏆 Stratified |
| **Parallel** | ✅ Yes | ❌ No | 🏆 Stratified |
| **Convergence Diagnostics** | ❌ No | ✅ Yes | 🏆 Hierarchical |
| **Information Borrowing** | ❌ No | ✅ Yes | 🏆 Hierarchical |
| **Best For** | Exploratory, speed | Publication, small categories | — |

### Confusion Matrix

Stratified Model:
```text
         Group 1 Group 2 Group 3
  1        209       9       0
  2         16     206      15
  3          0      10     210
```
Hierarchical Model:
```text
         Group 1 Group 2 Group 3
  1        215      19       0
  2         10     182      10
  3          0      24     215
```
### Per-Category Accuracy
| Category | Stratified | Hierarchical | Difference |
|:---|:---:|:---:|:---:|
| AA | 93.3% | 92.0% | ⬇️ -1.3% |
| BB | 92.0% | 89.3% | ⬇️ -2.7% |
| CC | 82.7% | 78.7% | ⬇️ -4.0% |
| DD | 🏆 **100.0%** | 🏆 **100.0%** | ➖ 0.0% |
| EE | 94.7% | 93.3% | ⬇️ -1.3% |
| FF | 85.3% | 82.7% | ⬇️ -2.7% |
| GG | 🏆 **100.0%** | 97.3% | ⬇️ -2.7% |
| HH | 94.7% | 94.7% | ➖ 0.0% |
| II | 90.7% | 88.0% | ⬇️ -2.7% |
| **Mean** | **🏆 92.6%** | **90.7%** | **⬇️ -1.9%** |

### MCMC Convergence Diagnostics (Hierarchical Model)
```text
Metric	         Value
Min ESS	          141
Mean ESS	       253
Max R_hat	       < 1.1 (converged)
```
## 🎯 When to Use Which Model

| Scenario | Recommended Model | Reason |
|:---|:---|:---|
| Categories are truly independent | Stratified (`fuss_PARALLEL_mcmc`) | No information should be shared |
| Categories share the same components | **Hierarchical** (`fuss_COVARIATE_mcmc`) | Borrows strength across categories |
| Small sample sizes per category | **Hierarchical** (`fuss_COVARIATE_mcmc`) | Stabilizes estimates |
| Quick exploratory analysis | Stratified (`fuss_PARALLEL_mcmc`) | Faster (parallel processing) |
| Publication-quality inference | **Hierarchical** (`fuss_COVARIATE_mcmc`) | More principled, convergence diagnostics |
| Maximum accuracy (synthetic data) | Stratified (`fuss_PARALLEL_mcmc`) | Higher accuracy on this benchmark |


## 📦 Generate Custom Data

You can also generate custom multimodal data with different parameters:
```r
# Generate custom dataset
custom_data <- MultiModalR::create_multimodal_dummy(
  seed = 12,
  n_categories = 6,
  n_per_group = 30,
  n_subgroups = 4
)
```
### Happy multimoda(e)ling! 🏔️️ 📊 🎯

## 📝 Citation

If you use **MultiModalR** in your research, please cite the original paper: 
*after publication*