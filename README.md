# Simple Bayesian Models

This repository contains two small Bayesian modelling projects built on the same input data. The projects use the Adult Census Income dataset to model the probability that an observation belongs to the high-income group.

The repository is organised as a reproducible analytical workspace rather than as a collection of course files. Each model directory contains the source code, input data, generated tables or figures, and the final report where available.

## Models

### Bayesian Linear Probability Model

Directory: `bayesian_linear_probability_model/`

This project estimates a Bayesian linear probability model for the binary outcome `income_high`, where `income_high = 1` indicates income above USD 50,000. The model uses a normal-gamma prior structure and derives posterior simulation draws analytically from the conjugate Bayesian regression setup. The explanatory variables include education, age, age squared, weekly working hours, sex, marital status, and a capital-gain indicator.

The model is useful as a transparent baseline: coefficients are directly interpreted as changes in probability, while posterior uncertainty is summarised through simulated draws and HPDI intervals.

### Bayesian Logistic Regression Model

Directory: `bayesian_logistic_regression_model/`

This project estimates a Bayesian logistic regression model for the same binary outcome. The likelihood is Bernoulli-logit and the regression coefficients use Student-t prior distributions elicited on the odds-ratio scale. Posterior inference is performed with Stan through the `rstan` interface using MCMC sampling.

The model is a probabilistically coherent extension of the linear probability specification: fitted probabilities are constrained to `[0, 1]`, odds ratios are reported directly, and posterior predictive checks are used to assess whether the model reproduces the observed event rate.

## Repository Structure

```text
.
|-- bayesian_linear_probability_model/
|   |-- code/
|   |-- data/
|   |-- figures/
|   `-- report/
|-- bayesian_logistic_regression_model/
|   |-- code/
|   |-- data/
|   |-- figures/
|   |-- outputs/
|   |-- report/
|   |-- stan/
|   `-- tables/
|-- LICENSE
`-- README.md
```

## Input Data

Both projects use the Adult Census Income dataset. The same raw input files are included in each model directory under `data/`:

- `adult.data`
- `adult.names`
- `old.adult.names`
- `Index`

The modelling scripts restrict the data to observations from the United States and draw a reproducible sample of 500 observations. Random seeds and MCMC configuration parameters are defined explicitly in the R scripts.

## Running the Code

The code is written in R. Run each script from either the model directory or its `code/` subdirectory.

For the linear probability model:

```r
setwd("bayesian_linear_probability_model")
source("code/bayesian_linear_probability_model.R")
```

For the logistic regression model:

```r
setwd("bayesian_logistic_regression_model")
source("code/bayesian_logit_student_t.R")
```

The logistic model requires a working Stan toolchain. On Windows this usually means installing Rtools compatible with the installed R version. The script can use a local `Makevars.win` file if one is supplied by the user, but the repository does not require such a file for platforms where the Stan toolchain is already configured.

## Technologies

The projects use:

- R
- Stan through `rstan`
- `MASS`
- `coda`
- `posterior`
- `ggplot2`
- `xtable`
- LaTeX for the written report in the logistic regression project

## Educational and Analytical Purpose

The repository illustrates how two Bayesian approaches can be applied to the same binary income-classification problem. The linear probability model provides a simple conjugate Bayesian baseline, while the logistic regression model introduces a more appropriate likelihood, prior elicitation through odds ratios, MCMC estimation, diagnostics, and posterior predictive assessment.
