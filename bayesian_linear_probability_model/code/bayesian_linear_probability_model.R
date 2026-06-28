# Bayesian linear probability model with a normal-gamma prior.

rm(list=ls())
cat("\014")
if (!is.null(dev.list())) dev.off()

args_file <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(args_file) > 0) {
  script_dir <- dirname(normalizePath(sub("^--file=", "", args_file[1]), winslash = "/"))
  project_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/")
} else {
  current_dir <- normalizePath(getwd(), winslash = "/")
  if (basename(current_dir) == "code") {
    project_dir <- normalizePath(file.path(current_dir, ".."), winslash = "/")
  } else {
    project_dir <- current_dir
  }
}

data_dir <- file.path(project_dir, "data")

library(MASS)

# ============================================================
# Nazwy kolumn
# ============================================================

columns <- c(
  "age",
  "workclass",
  "fnlwgt",
  "education",
  "education_num",
  "marital_status",
  "occupation",
  "relationship",
  "race",
  "sex",
  "capital_gain",
  "capital_loss",
  "hours_per_week",
  "native_country",
  "income"
)

# ============================================================
# Wczytanie zbioru
# ============================================================

adult_data <- read.csv(
  file.path(data_dir, "adult.data"),
  header = FALSE,
  col.names = columns,
  na.strings = "?",
  strip.white = TRUE
)

# ============================================================
# Ograniczenie zbioru do USA
# ============================================================

adult_us <- subset(adult_data, native_country == "United-States")

set.seed(123)  # dla powtarzalności wyników

nrow(adult_us)

# ============================================================
# Losowanie próby (N = 500)
# ============================================================

adult_sample <- adult_us[sample(nrow(adult_us), 500), ]

nrow(adult_sample)

prop.table(table(adult_us$income))
prop.table(table(adult_sample$income))

# ============================================================
# Przygotowanie próby do modelu
# ============================================================

adult_model <- adult_sample

adult_model$income_high <- ifelse(adult_model$income == ">50K", 1, 0)

adult_model$sex_male <- ifelse(adult_model$sex == "Male", 1, 0)

adult_model$married <- ifelse(
  adult_model$marital_status %in% c("Married-civ-spouse", "Married-spouse-absent"),
  1, 0
)

adult_model$capital_gain_ind <- ifelse(adult_model$capital_gain > 0, 1, 0)

# standaryzacja zmiennych ilościowych
adult_model$education_z <- scale(adult_model$education_num)
adult_model$age_z <- scale(adult_model$age)
adult_model$age_z2 <- adult_model$age_z^2
adult_model$hours_z <- scale(adult_model$hours_per_week)

# macierz X
X <- model.matrix(
  income_high ~ education_z + age_z + age_z2 + hours_z +
    sex_male + married + capital_gain_ind,
  data = adult_model
)

y <- adult_model$income_high

n <- length(y)
k <- ncol(X)


b0 <- c(
  0.25,   # intercept: bazowe prawdopodobieństwo wysokiego dochodu
  0.08,   # education_z: edukacja dodatnio wpływa na dochód
  0.05,   # age_z: wiek/doświadczenie dodatnio
  -0.03,   # age_z2: efekt malejący po pewnym wieku
  0.04,   # hours_z: więcej godzin pracy -> wyższy dochód
  0.05,   # sex_male: dodatni efekt dla mężczyzn
  0.12,   # married: osoby w małżeństwie częściej mają >50K
  0.20    # capital_gain_ind: dodatnie zyski kapitałowe silny sygnał dochodu
)

prior_sd <- c(
  0.20,  # intercept
  0.08,  # education
  0.08,  # age
  0.06,  # age^2
  0.08,  # hours
  0.10,  # male
  0.10,  # married
  0.15   # capital gain
)

B0 <- diag(prior_sd^2)

a0 <- 2
d0 <- 1

B0_inv <- solve(B0)

Bn <- solve(B0_inv + t(X) %*% X)

bn <- Bn %*% (B0_inv %*% b0 + t(X) %*% y)

an <- a0 + n / 2

dn <- d0 + 0.5 * (
  t(y) %*% y +
    t(b0) %*% B0_inv %*% b0 -
    t(bn) %*% solve(Bn) %*% bn
)

posterior_mean_beta <- as.vector(bn)

S <- 10000
k <- ncol(X)

# losowanie z posteriora
tau_draws <- rgamma(S, shape = an, rate = dn)   # tau = 1/sigma^2

beta_draws <- matrix(NA, nrow = S, ncol = k)
colnames(beta_draws) <- colnames(X)

for (s in 1:S) {
  beta_draws[s, ] <- mvrnorm(
    n = 1,
    mu = as.vector(bn),
    Sigma = (1 / tau_draws[s]) * Bn
  )
}

# tabela wartości oczekiwanych a posteriori
posterior_table <- data.frame(
  parameter = colnames(X),
  posterior_mean = colMeans(beta_draws),
  posterior_sd = apply(beta_draws, 2, sd)
)

posterior_table

# ============================================================
# Rozkłady brzegowe a posteriori
# ============================================================

par(mfrow = c(3, 3))

green_area <- rgb(24, 121, 104, 160, names = NULL, maxColorValue = 255)

for (j in 1:k) {
  hist(
    beta_draws[, j],
    breaks = 40,
    probability = TRUE,
    main = colnames(X)[j],
    xlab = expression(beta),
    border = "white",
    col = green_area
    )
  lines(density(beta_draws[, j]), lwd = 2)
  abline(v = 0, lty = 2)
}

# ============================================================
# Rozkłady brzegowe a priori i a posteriori parametrów beta
# ============================================================

dt_scaled <- function(x, mean, scale, df) {
  dt((x - mean) / scale, df = df) / scale
}

green_area <- rgb(24, 121, 104, 80, maxColorValue = 255)
green_line <- rgb(24, 121, 104, 255, maxColorValue = 255)

grey_area <- rgb(180, 180, 180, 90, maxColorValue = 255)
grey_line <- rgb(90, 90, 90, 255, maxColorValue = 255)

par(mfrow = c(3, 3))

for (j in 1:k) {
  
  # prior marginalny beta_j
  prior_mean_j <- b0[j]
  prior_scale_j <- sqrt((d0 / a0) * B0[j, j])
  prior_df <- 2 * a0
  
  # posterior marginalny beta_j
  post_mean_j <- as.numeric(bn[j])
  post_scale_j <- sqrt(as.numeric((dn / an) * Bn[j, j]))
  post_df <- 2 * an
  
  # wspólna siatka wartości beta
  beta_min <- min(
    prior_mean_j - 4 * prior_scale_j,
    post_mean_j - 4 * post_scale_j
  )
  
  beta_max <- max(
    prior_mean_j + 4 * prior_scale_j,
    post_mean_j + 4 * post_scale_j
  )
  
  beta_space <- seq(beta_min, beta_max, length.out = 1000)
  
  prior_dens <- dt_scaled(
    beta_space,
    mean = prior_mean_j,
    scale = prior_scale_j,
    df = prior_df
  )
  
  post_dens <- dt_scaled(
    beta_space,
    mean = post_mean_j,
    scale = post_scale_j,
    df = post_df
  )
  
  ymax <- max(prior_dens, post_dens) * 1.15
  
  plot(
    beta_space,
    prior_dens,
    type = "l",
    lwd = 2,
    col = grey_line,
    ylim = c(0, ymax),
    main = colnames(X)[j],
    xlab = expression(beta),
    ylab = "gęstość",
    bty = "n",
    las = 1
  )
  
  polygon(
    c(beta_space, rev(beta_space)),
    c(prior_dens, rep(0, length(beta_space))),
    col = grey_area,
    border = NA
  )
  
  lines(beta_space, prior_dens, lwd = 2, col = grey_line)
  
  polygon(
    c(beta_space, rev(beta_space)),
    c(post_dens, rep(0, length(beta_space))),
    col = green_area,
    border = NA
  )
  
  lines(beta_space, post_dens, lwd = 2, col = green_line)
  
  abline(v = 0, lty = 2)
  abline(v = prior_mean_j, col = grey_line, lwd = 2)
  abline(v = post_mean_j, col = green_line, lwd = 2)
  
  legend(
    "topright",
    legend = c("a priori", "a posteriori"),
    fill = c(grey_area, green_area),
    border = NA,
    bty = "n",
    cex = 0.75
  )
}

# ============================================================
# HPDI
# ============================================================

par(mfrow = c(1, 1))

hpdi <- function(x, prob = 0.95) {
  x <- sort(x)
  n <- length(x)
  m <- floor(prob * n)
  
  widths <- x[(m + 1):n] - x[1:(n - m)]
  i <- which.min(widths)
  
  c(
    lower = x[i],
    upper = x[i + m]
  )
}

hpdi_table <- data.frame(
  parameter = colnames(X),
  posterior_mean = colMeans(beta_draws),
  HPDI_lower = apply(beta_draws, 2, function(x) hpdi(x)[1]),
  HPDI_upper = apply(beta_draws, 2, function(x) hpdi(x)[2]),
  P_beta_gt_0 = apply(beta_draws, 2, function(x) mean(x > 0)),
  P_beta_lt_0 = apply(beta_draws, 2, function(x) mean(x < 0))
)

hpdi_table


# ============================================================
# czynniki Bayesa
# ============================================================

# funkcja gęstości brzegowej t-Studenta
dt_scaled <- function(x, mean, scale, df) {
  dt((x - mean) / scale, df = df) / scale
}

BF_table <- data.frame(
  parameter = colnames(X),
  BF_01 = NA,
  BF_10 = NA
)

for (j in 1:k) {
  
  # prior brzegowy beta_j
  prior_mean_j <- b0[j]
  prior_scale_j <- sqrt((d0 / a0) * B0[j, j])
  prior_df <- 2 * a0
  
  prior_density_0 <- dt_scaled(
    x = 0,
    mean = prior_mean_j,
    scale = prior_scale_j,
    df = prior_df
  )
  
  # posterior brzegowy beta_j
  post_mean_j <- bn[j]
  post_scale_j <- sqrt((dn / an) * Bn[j, j])
  post_df <- 2 * an
  
  post_density_0 <- dt_scaled(
    x = 0,
    mean = post_mean_j,
    scale = post_scale_j,
    df = post_df
  )
  
  # Savage-Dickey:
  # BF_01 = p(beta_j = 0 | y) / p(beta_j = 0)
  BF_01 <- post_density_0 / prior_density_0
  BF_10 <- 1 / BF_01
  
  BF_table$BF_01[j] <- BF_01
  BF_table$BF_10[j] <- BF_10
}

BF_table

# ============================================================
# Predykcja
# ============================================================

i <- 1

x_new <- X[i, , drop = FALSE]
y_true <- y[i]

# posterior predictive dla wartości oczekiwanej E[y|x]
mu_draws <- as.vector(x_new %*% t(beta_draws))

point_forecast <- mean(mu_draws)
interval_mu <- hpdi(mu_draws, prob = 0.95)

point_forecast
interval_mu
y_true

sigma_draws <- sqrt(1 / tau_draws)

y_pred_draws <- rnorm(
  S,
  mean = mu_draws,
  sd = sigma_draws
)

pred_point <- mean(y_pred_draws)
pred_interval <- hpdi(y_pred_draws, prob = 0.95)

pred_point
pred_interval


point_forecast_clipped <- min(max(point_forecast, 0), 1)
interval_mu_clipped <- pmin(pmax(interval_mu, 0), 1)

point_forecast_clipped
interval_mu_clipped

