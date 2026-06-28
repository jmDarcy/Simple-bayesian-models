# Simple Bayesian Models

This repository contains two Bayesian modelling projects based on the Adult Census Income dataset. Both projects analyse the probability that an individual belongs to the high-income group, defined as annual income exceeding USD 50,000.

The repository is organised as a reproducible analytical workspace rather than as a collection of course files. Each model directory contains the source code, input data, generated tables or figures, and the final report where available.

The two models are intentionally built on a comparable empirical setup. The first project uses a Bayesian linear probability model as a transparent conjugate baseline. The second project extends the analysis to a Bayesian logistic regression model, which is more appropriate for a binary response variable.

## Research Question

The central modelling question is:

> Which observable demographic and labour-market characteristics are associated with the probability of earning more than USD 50,000 per year?

The response variable is defined as

$$
Y_i =
\begin{cases}
1, & \text{if individual } i \text{ earns more than USD 50,000 per year}, \\
0, & \text{otherwise}.
\end{cases}
$$

The explanatory variables include:

- education level,
- age,
- squared standardised age,
- weekly working hours,
- sex,
- marital status,
- positive capital-gain indicator.

The empirical motivation follows the human-capital interpretation of income determination. Education and experience-related variables are expected to be positively associated with income, while the squared age term allows the income profile to flatten over the life cycle.

## Models

### 1. Bayesian Linear Probability Model

Directory: `bayesian_linear_probability_model/`

The first project estimates a Bayesian linear probability model:

$$
Y_i = x_i^\top \beta + \varepsilon_i,
$$

where

$$
\mathbb{E}(Y_i \mid x_i) = \mathbb{P}(Y_i = 1 \mid x_i) = x_i^\top \beta.
$$

The model uses the normal likelihood

$$
y \mid \beta, \sigma^2 \sim \mathcal{N}(X\beta, \sigma^2 I_n),
$$

with a conjugate normal-gamma prior structure:

$$
\beta \mid \sigma^2 \sim \mathcal{N}(b_0, \sigma^2 B_0),
\qquad
\tau = \sigma^{-2} \sim \mathrm{Gamma}(a_0, d_0).
$$

This specification allows posterior simulation from closed-form conditional distributions. It is therefore useful as a mathematically transparent Bayesian baseline.

The model is easy to interpret: coefficients approximate changes in the probability of high income, measured in percentage points. For example, a coefficient of `0.13` on standardised education can be read as an approximate 13 percentage point increase in the probability of high income when education increases by one standard deviation, holding other variables constant.

However, the linear probability model has important limitations:

- it may generate fitted values outside the interval `[0, 1]`,
- the error variance is heteroskedastic for binary outcomes,
- the coefficients describe conditional associations, not causal effects.

The report therefore treats the model as an interpretable baseline rather than as the final probabilistic specification.

### 2. Bayesian Logistic Regression Model

Directory: `bayesian_logistic_regression_model/`

The second project estimates a Bayesian logistic regression model for the same binary outcome:

$$
Y_i \sim \mathrm{Bernoulli}(p_i),
$$

The probability of high income is linked to the predictors through the log-odds transformation:

```math
\log\left(\frac{p_i}{1-p_i}\right)
=
x_i^\top \beta.
```

Equivalently, the fitted probability is

```math
p_i
=
\frac{\exp(x_i^\top \beta)}{1+\exp(x_i^\top \beta)}.
```

The model constrains fitted probabilities to the interval `[0, 1]` and provides a natural interpretation through odds ratios:

$$
\mathrm{OR}_j = \exp(\beta_j).
$$

Thus, if `exp(beta_j) = 2`, then a one-unit increase in variable `x_j`, holding other variables fixed, doubles the odds of belonging to the high-income group.

The prior distributions are specified as Student-t distributions:

$$
\beta_j \sim t_\nu(m_j, s_j),
$$

with prior locations elicited on the odds-ratio scale:

$$
m_j = \log(\mathrm{OR}_j).
$$

For the intercept, the prior location is obtained from an assumed baseline probability:

$$
m_0 =
\log\left(\frac{p_0}{1-p_0}\right).
$$

The Student-t prior has heavier tails than a normal prior, which makes it a reasonable compromise between economically motivated prior information and robustness to stronger effects in the data.

Posterior inference is performed using Stan through the `rstan` interface. The model is sampled with NUTS, the No-U-Turn Sampler, an adaptive variant of Hamiltonian Monte Carlo.

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
