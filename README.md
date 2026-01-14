# MultiModalR 🏔️ <img src="https://latex.codecogs.com/svg.latex?\color{green}w_1\mathcal{N}(\mu_1,\sigma_1^2)+w_2\mathcal{N}(\mu_2,\sigma_2^2)+w_3\mathcal{N}(\mu_3,\sigma_3^2)" height="30" align="center">

[![R](https://img.shields.io/badge/R-≥4.0-green?style=for-the-badge&logo=r)](https://www.r-project.org/)
[![C++](https://img.shields.io/badge/C++-RcppArmadillo-green?style=for-the-badge&logo=cplusplus)](https://isocpp.org/)


**MultiModalR** performs Bayesian mixture modeling for multimodal data. It detects subpopulations and assigns probabilistic memberships using a fast C++ implementation of Metropolis-Hastings-within-partial-Gibbs sampling.

## Features

- **Fast MCMC**: Speed-optimized Metropolis-Hastings-within-partial-Gibbs implementation
- **Subpopulation Detection**: Automatic detection of multimodal components
- **Bayesian Probability Assignment**: Soft assignment with probability estimates
- **Parallel Processing**: Built-in support for multi-core computation
- **Validation Tools**: Built-in plotting and validation functions
- **Flexible Input**: Works with various categorical data structures

## Installation

### From GitHub
```r
# Install from GitHub
devtools::install_github("DijoG/MultiModalR")

# Or using remotes
remotes::install_github("DijoG/MultiModalR")
```
## Dependencies

`Rcpp`, `dplyr`, `furrr`, `future`

## Quick Start

### Basic Usage
```r
library(MultiModalR)

# Prepare data and run analysis
MultiModalR::fuss_PARALLEL(
  data = df_GROUPS,
  varCLASS = "Category", 
  varY = "Value", 
  varID = "ID",
  method = "dpi", 
  within = 1, 
  maxNGROUP = 5, 
  df_prob = FALSE, 
  out_dir = ".../output", 
  n_workers = 3
)
```
## Detailed Example

### Data Generation

For dummy data creation the **truncnorm** package is needed.

```r
library(tidyverse);library(truncnorm)

# Set seed for reproducibility
set.seed(5)

# Define nine categories with three subpopulations
categories <- rep(LETTERS[1:9], each = 75)
subpopulations <- rep(rep(c("Group 1", "Group 2", "Group 3"), each = 25), times = 9)

# Generate data with single-peaked distributions within each subgroup
values <- c(
  rtruncnorm(25, a = 5, b = 10, mean = 6, sd = 0.4), rtruncnorm(25, a = 5, b = 10, mean = 7.5, sd = 0.4), rtruncnorm(25, a = 5, b = 10, mean = 9, sd = 0.4),
  rtruncnorm(25, a = 5, b = 10, mean = 6.2, sd = 0.5), rtruncnorm(25, a = 5, b = 10, mean = 7.7, sd = 0.5), rtruncnorm(25, a = 5, b = 10, mean = 9.2, sd = 0.5),
  rtruncnorm(25, a = 5, b = 10, mean = 5.8, sd = 0.6), rtruncnorm(25, a = 5, b = 10, mean = 7.4, sd = 0.6), rtruncnorm(25, a = 5, b = 10, mean = 8.9, sd = 0.6),
  rtruncnorm(25, a = 5, b = 10, mean = 6.1, sd = 0.4), rtruncnorm(25, a = 5, b = 10, mean = 7.8, sd = 0.4), rtruncnorm(25, a = 5, b = 10, mean = 9.3, sd = 0.4),
  rtruncnorm(25, a = 5, b = 10, mean = 6.3, sd = 0.5), rtruncnorm(25, a = 5, b = 10, mean = 7.8, sd = 0.5), rtruncnorm(25, a = 5, b = 10, mean = 9.4, sd = 0.5),
  rtruncnorm(25, a = 5, b = 10, mean = 5.9, sd = 0.6), rtruncnorm(25, a = 5, b = 10, mean = 7.5, sd = 0.6), rtruncnorm(25, a = 5, b = 10, mean = 9.2, sd = 0.6),
  rtruncnorm(25, a = 5, b = 10, mean = 6.4, sd = 0.4), rtruncnorm(25, a = 5, b = 10, mean = 7.9, sd = 0.4), rtruncnorm(25, a = 5, b = 10, mean = 9.5, sd = 0.4),
  rtruncnorm(25, a = 5, b = 10, mean = 6.0, sd = 0.5), rtruncnorm(25, a = 5, b = 10, mean = 7.6, sd = 0.5), rtruncnorm(25, a = 5, b = 10, mean = 9.3, sd = 0.5),
  rtruncnorm(25, a = 5, b = 10, mean = 6.2, sd = 0.6), rtruncnorm(25, a = 5, b = 10, mean = 7.8, sd = 0.6), rtruncnorm(25, a = 5, b = 10, mean = 9.6, sd = 0.6)
)

# Create data frame
df <- data.frame(Category = categories, Subpopulation = subpopulations, Value = values) %>%
  mutate(ID = 1:nrow(.))    
```
### Data Visualization
```r
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
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/MMR/MMR_01.png">

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
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/MMR/MMR_02.png">


### Parallel Processing Setup 
```r
# Configure parallel processing
cores <- 3
```
### Running Analysis
```r
library(furrr)

# Analysis with probability output
MultiModalR::fuss_PARALLEL(
  data = df,
  varCLASS = "Category",
  varY = "Value",
  varID = "ID",
  method = "dpi",
  within = 1,
  maxNGROUP = 5,
  out_dir = "D:/MultiModalR/test",  
  n_workers = cores
)
tictoc::toc()
# Processing time: 9.09 sec (3 cores)
```
### Output 

The function generates:
  - **Data CSV** files: Original data with assigned subgroups and probabilities
  
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/MMR/MMR_csv.png">

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

```r
# Validate subgroup assignments
MultiModalR::plot_VALIDATION(
  "D:/MultiModalR/test", 
  df, 
  subpop_col = "Subpopulation", 
  value_col = "Value",
  id_col = "ID")
```
<img align="bottom" src="https://raw.githubusercontent.com/DijoG/storage/main/MMR/MMR_03.png">

**Validation results show accurate subgroup assignment across categories.**

**Happy multimoda(e)ling!** 🏔️