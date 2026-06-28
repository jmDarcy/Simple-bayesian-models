# EKONOMETRIA BAYESOWSKA - ZADANIE DOMOWE 2

# ============================================================
# 1. Pakiety i ustawienia
# ============================================================

# ustawiamy ścieżki lokalnie dla bieżącej sesji R.
rtools_paths <- c(
  "C:/rtools43/usr/bin",
  "C:/rtools43/x86_64-w64-mingw32.static.posix/bin",
  "/usr/bin",
  "/x86_64-w64-mingw32.static.posix/bin"
)
existing_rtools_paths <- rtools_paths[dir.exists(rtools_paths)]
if (length(existing_rtools_paths) > 0) {
  Sys.setenv(PATH = paste(c(existing_rtools_paths, Sys.getenv("PATH")), collapse = .Platform$path.sep))
}

# Ustalenie katalogu projektu i danych.

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

# Pakiety potrzebne do analizy
library(rstan)
library(coda)
library(ggplot2)

has_bayesplot <- requireNamespace("bayesplot", quietly = TRUE)
if (has_bayesplot) {
  library(bayesplot)
} else {
  message("Pakiet bayesplot nie jest zainstalowany; pomijam wykresy bayesplot.")
}


options(mc.cores = parallel::detectCores())
rstan_options(auto_write = FALSE)

# Parametry
DATA_SEED <- 121732
MCMC_SEED <- 20260609
N_SAMPLE <- 500
MCMC_CHAINS <- 4
MCMC_ITER <- 4000
MCMC_WARMUP <- 1000
MCMC_THIN <- 1
ADAPT_DELTA <- 0.95
MAX_TREEDEPTH <- 12
PRIOR_DF <- 4

DATA_SEED
MCMC_SEED

# ============================================================
# 2. Wczytanie danych Adult Census Income
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

adult_raw <- read.csv(
  file.path(data_dir, "adult.data"),
  header = FALSE,
  col.names = columns,
  na.strings = "?",
  strip.white = TRUE
)

# Oglądamy pierwsze obserwacje i strukturę danych.
head(adult_raw)
str(adult_raw)

#View(adult_raw)

# ============================================================
# 3. Podstawowy podgląd danych
# ============================================================

dim(adult_raw)
summary(adult_raw[, c("age", "education_num", "capital_gain", "capital_loss", "hours_per_week")])

table(adult_raw$income)
prop.table(table(adult_raw$income))

table(adult_raw$native_country, useNA = "ifany")[1:10]

# ============================================================
# 4. Ograniczenie danych do obserwacji z USA
# ============================================================

adult_us <- subset(adult_raw, native_country == "United-States")

dim(adult_us)
prop.table(table(adult_us$income))

# ============================================================
# 5. Losowanie próby 500 obserwacji
# ============================================================

set.seed(DATA_SEED)
adult_sample <- adult_us[sample(seq_len(nrow(adult_us)), N_SAMPLE), ]

dim(adult_sample)
head(adult_sample)
prop.table(table(adult_sample$income))

# View(adult_sample)

# ============================================================
# 6. Konstrukcja zmiennej income_high
# ============================================================

adult_model <- adult_sample
adult_model$income_high <- ifelse(adult_model$income == ">50K", 1L, 0L)

table(adult_model$income, adult_model$income_high)
prop.table(table(adult_model$income_high))

# ============================================================
# 7. Konstrukcja zmiennych objaśniających
# ============================================================

adult_model$sex_male <- ifelse(adult_model$sex == "Male", 1L, 0L)

adult_model$married <- ifelse(
  adult_model$marital_status %in% c("Married-civ-spouse", "Married-spouse-absent"),
  1L,
  0L
)

adult_model$capital_gain_ind <- ifelse(adult_model$capital_gain > 0, 1L, 0L)

head(adult_model[, c("income_high", "sex", "sex_male", "marital_status", "married",
                     "capital_gain", "capital_gain_ind")])

table(adult_model$sex_male)
table(adult_model$married)
table(adult_model$capital_gain_ind)

# ============================================================
# 8. Standaryzacja zmiennych ilościowych
# ============================================================

adult_model$education_z <- as.numeric(scale(adult_model$education_num))
adult_model$age_z <- as.numeric(scale(adult_model$age))
adult_model$age_z2 <- adult_model$age_z^2
adult_model$hours_z <- as.numeric(scale(adult_model$hours_per_week))

head(adult_model[, c("education_num", "education_z", "age", "age_z", "age_z2",
                     "hours_per_week", "hours_z")])

summary(adult_model[, c("education_z", "age_z", "age_z2", "hours_z")])

# ============================================================
# 9. Podstawowa eksploracja próby
# ============================================================

data_summary <- data.frame(
  wielkość = c(
    "obserwacje w pliku adult.data",
    "obserwacje z USA",
    "wylosowana próba",
    "udział income_high = 1 w próbie"
  ),
  wartość = c(
    nrow(adult_raw),
    nrow(adult_us),
    nrow(adult_model),
    round(mean(adult_model$income_high), 3)
  )
)

data_summary

aggregate(income_high ~ sex_male, data = adult_model, mean)
aggregate(income_high ~ married, data = adult_model, mean)
aggregate(income_high ~ capital_gain_ind, data = adult_model, mean)

hist(adult_model$age, breaks = 25, col = rgb(24, 121, 104, 120, maxColorValue = 255),
     border = "white", main = "Wiek w próbie", xlab = "wiek")

hist(adult_model$education_num, breaks = 15, col = rgb(24, 121, 104, 120, maxColorValue = 255),
     border = "white", main = "Liczba lat edukacji", xlab = "education_num")

# ============================================================
# 10. Elicytacja priorów Studenta-t
# ============================================================

# W modelu logitowym beta oznacza zmianę logarytmu szans.
# Dlatego wiedzę a priori wygodnie zapisać przez ilorazy szans:
# prior_location_beta = log(oczekiwany_OR).

parameter_labels <- c(
  "(Intercept)" = "wyraz wolny",
  "education_z" = "edukacja, +1 SD",
  "age_z" = "wiek, +1 SD",
  "age_z2" = "kwadrat wieku",
  "hours_z" = "godziny pracy, +1 SD",
  "sex_male" = "mężczyzna",
  "married" = "w małżeństwie",
  "capital_gain_ind" = "dodatnie zyski kapitałowe"
)

prior_or_location <- c(
  "(Intercept)" = NA,
  "education_z" = 1.60,
  "age_z" = 1.30,
  "age_z2" = 0.85,
  "hours_z" = 1.25,
  "sex_male" = 1.40,
  "married" = 2.00,
  "capital_gain_ind" = 3.00
)

prior_location <- c(
  "(Intercept)" = qlogis(0.25),
  "education_z" = unname(log(prior_or_location["education_z"])),
  "age_z" = unname(log(prior_or_location["age_z"])),
  "age_z2" = unname(log(prior_or_location["age_z2"])),
  "hours_z" = unname(log(prior_or_location["hours_z"])),
  "sex_male" = unname(log(prior_or_location["sex_male"])),
  "married" = unname(log(prior_or_location["married"])),
  "capital_gain_ind" = unname(log(prior_or_location["capital_gain_ind"]))
)

prior_scale <- c(
  "(Intercept)" = 0.90,
  "education_z" = 0.35,
  "age_z" = 0.35,
  "age_z2" = 0.30,
  "hours_z" = 0.30,
  "sex_male" = 0.45,
  "married" = 0.45,
  "capital_gain_ind" = 0.70
)

prior_table <- data.frame(
  parameter = names(prior_location),
  opis = unname(parameter_labels[names(prior_location)]),
  prior_df = PRIOR_DF,
  prior_location_beta = as.numeric(prior_location),
  prior_scale_beta = as.numeric(prior_scale),
  central_OR = as.numeric(prior_or_location[names(prior_location)])
)

prior_table
# View(prior_table)

# ============================================================
# 11. Przygotowanie danych dla Stan
# ============================================================

X <- model.matrix(
  income_high ~ education_z + age_z + age_z2 + hours_z +
    sex_male + married + capital_gain_ind,
  data = adult_model
)

y <- as.integer(adult_model$income_high)

head(X)
summary(y)
colnames(X)

i_pred <- 1L
X_pred <- X[i_pred, , drop = FALSE]
y_true_pred <- y[i_pred]

X_pred
y_true_pred

stan_data <- list(
  N = nrow(X),
  K = ncol(X),
  X = X,
  y = y,
  beta_prior_location = as.numeric(prior_location[colnames(X)]),
  beta_prior_scale = as.numeric(prior_scale[colnames(X)]),
  beta_prior_df = PRIOR_DF,
  N_pred = 1L,
  X_pred = X_pred
)

str(stan_data)

# ============================================================
# 12. Kod modelu Stan
# ============================================================

# Model:
# y_i ~ Bernoulli(p_i)
# logit(p_i) = x_i' beta
# beta_j ~ Student-t(df, location_j, scale_j)

stan_code <- "
data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  int<lower=0, upper=1> y[N];
  vector[K] beta_prior_location;
  vector<lower=0>[K] beta_prior_scale;
  real<lower=1> beta_prior_df;
  int<lower=1> N_pred;
  matrix[N_pred, K] X_pred;
}

parameters {
  vector[K] beta;
}

model {
  beta ~ student_t(beta_prior_df, beta_prior_location, beta_prior_scale);
  y ~ bernoulli_logit(X * beta);
}

generated quantities {
  vector[N] eta;
  vector[N] p_hat;
  vector[N] log_lik;
  int<lower=0, upper=1> y_rep[N];
  vector[N_pred] p_pred;
  int<lower=0, upper=1> y_pred_rep[N_pred];

  eta = X * beta;

  for (n in 1:N) {
    p_hat[n] = inv_logit(eta[n]);
    log_lik[n] = bernoulli_logit_lpmf(y[n] | eta[n]);
    y_rep[n] = bernoulli_rng(p_hat[n]);
  }

  for (m in 1:N_pred) {
    p_pred[m] = inv_logit(dot_product(row(X_pred, m), beta));
    y_pred_rep[m] = bernoulli_rng(p_pred[m]);
  }
}
"

cat(stan_code)

# ============================================================
# 13. Estymacja bayesowskiego modelu logitowego
# ============================================================

# Jeżeli obiekt fit już istnieje w Environment, poniższy fragment go nie nadpisze.

if (!exists("fit")) {
  stan_model_logit <- rstan::stan_model(model_code = stan_code)

  fit <- rstan::sampling(
    object = stan_model_logit,
    data = stan_data,
    seed = MCMC_SEED,
    chains = MCMC_CHAINS,
    iter = MCMC_ITER,
    warmup = MCMC_WARMUP,
    thin = MCMC_THIN,
    control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH)
  )
} else {
  message("Obiekt fit już istnieje w środowisku. Pomijam ponowną estymację.")
}

fit
print(fit)

# ============================================================
# 14. Diagnostyka MCMC
# ============================================================

pars_beta <- paste0("beta[", seq_len(ncol(X)), "]")

fit_summary_raw <- rstan::summary(fit, pars = c(pars_beta, "p_pred"))$summary
fit_summary_raw

diagnostics_table <- data.frame(
  parameter = colnames(X),
  n_eff = fit_summary_raw[pars_beta, "n_eff"],
  Rhat = fit_summary_raw[pars_beta, "Rhat"]
)

diagnostics_table

# Rhat blisko 1 sugeruje zgodność łańcuchów.
# Duże n_eff oznacza wysoką efektywną liczbę próbek.

rstan::traceplot(fit, pars = pars_beta, inc_warmup = FALSE, nrow = 4)

if (has_bayesplot) {
  posterior_array <- as.array(fit, pars = pars_beta)
  bayesplot::mcmc_trace(posterior_array)
  bayesplot::mcmc_dens_overlay(posterior_array)
}

sampler_params <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
divergences <- sum(vapply(sampler_params, function(x) sum(x[, "divergent__"]), numeric(1)))
max_treedepth_hits <- sum(vapply(sampler_params, function(x) sum(x[, "treedepth__"] >= MAX_TREEDEPTH), numeric(1)))
mean_accept_stat <- mean(unlist(lapply(sampler_params, function(x) x[, "accept_stat__"])))

divergences
max_treedepth_hits
mean_accept_stat

# ============================================================
# 15. Podsumowanie rozkładów a posteriori
# ============================================================

draws <- rstan::extract(fit, permuted = TRUE, inc_warmup = FALSE)

beta_draws <- draws$beta
colnames(beta_draws) <- colnames(X)

head(beta_draws)
dim(beta_draws)

hpdi <- function(x, prob = 0.95) {
  x <- sort(as.numeric(x))
  n <- length(x)
  m <- floor(prob * n)
  widths <- x[(m + 1):n] - x[seq_len(n - m)]
  i <- which.min(widths)
  c(lower = x[i], upper = x[i + m])
}

posterior_summary <- data.frame(
  parameter = colnames(X),
  opis = unname(parameter_labels[colnames(X)]),
  posterior_mean = colMeans(beta_draws),
  posterior_median = apply(beta_draws, 2, median),
  posterior_sd = apply(beta_draws, 2, sd),
  HPDI_lower = apply(beta_draws, 2, function(x) hpdi(x)[1]),
  HPDI_upper = apply(beta_draws, 2, function(x) hpdi(x)[2]),
  P_beta_gt_0 = apply(beta_draws, 2, function(x) mean(x > 0)),
  P_beta_lt_0 = apply(beta_draws, 2, function(x) mean(x < 0))
)

posterior_summary
# View(posterior_summary)

# beta jest na skali log-szans.
# Dodatnia beta oznacza wzrost szans income_high = 1.

# ============================================================
# 16. Wizualizacja priorów i posteriorów
# ============================================================

green_area <- rgb(24, 121, 104, 90, maxColorValue = 255)
green_line <- rgb(24, 121, 104, 255, maxColorValue = 255)
grey_area <- rgb(180, 180, 180, 80, maxColorValue = 255)
grey_line <- rgb(90, 90, 90, 255, maxColorValue = 255)

dt_scaled <- function(x, mean, scale, df) {
  dt((x - mean) / scale, df = df) / scale
}

# Histogramy brzegowych rozkładów posteriori.
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))
for (j in seq_len(ncol(X))) {
  hist(
    beta_draws[, j],
    breaks = 40,
    probability = TRUE,
    main = colnames(X)[j],
    xlab = expression(beta),
    border = "white",
    col = green_area
  )
  lines(density(beta_draws[, j]), lwd = 2, col = green_line)
  abline(v = 0, lty = 2)
}
par(mfrow = c(1, 1))

# Porównanie prioru i posterioru dla każdego parametru.
par(mfrow = c(3, 3), mar = c(4, 4, 3, 1))
for (j in seq_len(ncol(X))) {
  par_name <- colnames(X)[j]
  prior_mean_j <- prior_location[par_name]
  prior_scale_j <- prior_scale[par_name]
  post_mean_j <- mean(beta_draws[, j])
  post_sd_j <- sd(beta_draws[, j])

  beta_min <- min(prior_mean_j - 4 * prior_scale_j, post_mean_j - 4 * post_sd_j)
  beta_max <- max(prior_mean_j + 4 * prior_scale_j, post_mean_j + 4 * post_sd_j)
  beta_space <- seq(beta_min, beta_max, length.out = 1000)

  prior_dens <- dt_scaled(beta_space, prior_mean_j, prior_scale_j, PRIOR_DF)
  post_dens <- density(beta_draws[, j], from = beta_min, to = beta_max, n = 1000)
  ymax <- max(prior_dens, post_dens$y) * 1.15

  plot(
    beta_space,
    prior_dens,
    type = "l",
    lwd = 2,
    col = grey_line,
    ylim = c(0, ymax),
    main = par_name,
    xlab = expression(beta),
    ylab = "gęstość",
    bty = "n",
    las = 1
  )
  polygon(c(beta_space, rev(beta_space)), c(prior_dens, rep(0, length(beta_space))),
          col = grey_area, border = NA)
  lines(beta_space, prior_dens, lwd = 2, col = grey_line)
  polygon(c(post_dens$x, rev(post_dens$x)), c(post_dens$y, rep(0, length(post_dens$x))),
          col = green_area, border = NA)
  lines(post_dens$x, post_dens$y, lwd = 2, col = green_line)
  abline(v = 0, lty = 2)
  abline(v = prior_mean_j, col = grey_line, lwd = 2)
  abline(v = post_mean_j, col = green_line, lwd = 2)
  legend("topright", legend = c("a priori", "a posteriori"),
         fill = c(grey_area, green_area), border = NA, bty = "n", cex = 0.75)
}
par(mfrow = c(1, 1))

# ============================================================
# 17. Prawdopodobieństwa P(beta > 0 | y)
# ============================================================

prob_positive_table <- posterior_summary[, c("parameter", "opis", "P_beta_gt_0", "P_beta_lt_0")]
prob_positive_table

ggplot(prob_positive_table, aes(x = P_beta_gt_0, y = reorder(opis, P_beta_gt_0))) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_point(size = 2.4, color = green_line) +
  xlim(0, 1) +
  labs(x = "P(beta > 0 | y)", y = NULL,
       title = "Prawdopodobieństwo dodatniego wpływu") +
  theme_minimal(base_size = 11)

# P(beta > 0 | y) bliskie 1 oznacza silne wsparcie dla dodatniego efektu.
# P(beta > 0 | y) bliskie 0 oznacza silne wsparcie dla efektu ujemnego.

# ============================================================
# 18. Ilorazy szans exp(beta)
# ============================================================

or_draws <- exp(beta_draws)
colnames(or_draws) <- colnames(X)

odds_ratio_table <- data.frame(
  parameter = colnames(X),
  opis = unname(parameter_labels[colnames(X)]),
  OR_mean = colMeans(or_draws),
  OR_median = apply(or_draws, 2, median),
  OR_HPDI_lower = apply(or_draws, 2, function(x) hpdi(x)[1]),
  OR_HPDI_upper = apply(or_draws, 2, function(x) hpdi(x)[2])
)

odds_ratio_table
# View(odds_ratio_table)

# exp(beta) > 1 oznacza wzrost szans income_high = 1.
# exp(beta) < 1 oznacza spadek szans income_high = 1.

odds_ratio_plot_data <- subset(odds_ratio_table, parameter != "(Intercept)")

ggplot(odds_ratio_plot_data, aes(x = OR_median, y = reorder(opis, OR_median))) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = OR_HPDI_lower, xmax = OR_HPDI_upper),
                 height = 0.18, color = green_line, linewidth = 0.7) +
  geom_point(size = 2.4, color = green_line) +
  scale_x_log10() +
  labs(x = "iloraz szans exp(beta), skala logarytmiczna", y = NULL,
       title = "Ilorazy szans i 95% HPDI") +
  theme_minimal(base_size = 11)

# ============================================================
# 19. Predykcja bayesowska dla pierwszej obserwacji
# ============================================================

p_pred_draws <- as.numeric(draws$p_pred[, 1])
y_pred_rep <- as.integer(draws$y_pred_rep[, 1])

prediction_table <- data.frame(
  selected_observation = i_pred,
  observed_y = y_true_pred,
  posterior_mean_probability = mean(p_pred_draws),
  posterior_median_probability = median(p_pred_draws),
  HPDI_lower_probability = hpdi(p_pred_draws)[1],
  HPDI_upper_probability = hpdi(p_pred_draws)[2],
  predictive_probability_y1 = mean(y_pred_rep == 1L)
)

prediction_table

prediction_covariates <- data.frame(
  variable = colnames(X),
  value = as.numeric(X_pred[1, ])
)

prediction_covariates

hist(p_pred_draws, breaks = 40, probability = TRUE,
     col = green_area, border = "white",
     main = "Predykcja prawdopodobieństwa dla pierwszej obserwacji",
     xlab = "P(income > 50K | x, y)")
lines(density(p_pred_draws), col = green_line, lwd = 2)
abline(v = mean(p_pred_draws), col = green_line, lwd = 2)
abline(v = hpdi(p_pred_draws), lty = 2, col = "grey35")

# ============================================================
# 20. Posterior predictive check
# ============================================================

# Posterior predictive check sprawdza, czy model generuje dane podobne
# do zaobserwowanych. Tutaj porównujemy udział income_high = 1.

y_rep_draws <- draws$y_rep
observed_high_income_rate <- mean(y)
rep_high_income_rate <- rowMeans(y_rep_draws)

observed_high_income_rate
summary(rep_high_income_rate)
hpdi(rep_high_income_rate)

posterior_predictive_check_table <- data.frame(
  statistic = c(
    "udział y=1 w próbie",
    "średni udział y=1 w replikacjach",
    "HPDI dolny udział y=1 w replikacjach",
    "HPDI górny udział y=1 w replikacjach"
  ),
  value = c(
    observed_high_income_rate,
    mean(rep_high_income_rate),
    hpdi(rep_high_income_rate)[1],
    hpdi(rep_high_income_rate)[2]
  )
)

posterior_predictive_check_table

hist(rep_high_income_rate, breaks = 35, probability = TRUE,
     col = grey_area, border = "white",
     main = "Posterior predictive check: udział income_high = 1",
     xlab = "udział income_high = 1")
abline(v = observed_high_income_rate, col = green_line, lwd = 2)
legend("topright", legend = c("obserwowany udział"), lwd = 2,
       col = green_line, bty = "n")

# ============================================================
# 21. Krótkie komentarze analityczne
# ============================================================

# Najważniejsze obiekty do obejrzenia:
prior_table
diagnostics_table
posterior_summary
prob_positive_table
odds_ratio_table
prediction_table
posterior_predictive_check_table

# Interpretacja:
# - beta to efekt na skali log-szans.
# - exp(beta) to iloraz szans.
# - P(beta > 0 | y) pokazuje posteriorowe prawdopodobieństwo dodatniego wpływu.
# - HPDI to bayesowski przedział wiarygodności o największej gęstości.
# - posterior predictive check pokazuje, czy model odtwarza prostą cechę danych,
#   tutaj udział obserwacji z income_high = 1.

# Przykładowe odczyty po estymacji:
subset(posterior_summary, P_beta_gt_0 > 0.95)
subset(odds_ratio_table, parameter != "(Intercept)")
prediction_table

