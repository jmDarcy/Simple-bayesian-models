# EKONOMETRIA BAYESOWSKA - ZADANIE DOMOWE 2
# Bayesowski model logitowy z priorami Studenta-t i estymacja MCMC w rstan.

rm(list = ls())
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
stan_dir <- file.path(project_dir, "stan")
figures_dir <- file.path(project_dir, "figures")
tables_dir <- file.path(project_dir, "tables")
outputs_dir <- file.path(project_dir, "outputs")

for (dir_path in c(figures_dir, tables_dir, outputs_dir)) {
  if (!dir.exists(dir_path)) dir.create(dir_path, recursive = TRUE)
}

# Projektowy Makevars.win wskazuje kompilator Rtools dla rstan.
makevars_user <- file.path(project_dir, "Makevars.win")
if (file.exists(makevars_user)) {
  Sys.setenv(R_MAKEVARS_USER = makevars_user)
}

required_packages <- c("rstan", "coda", "posterior", "ggplot2", "xtable", "knitr")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Brak wymaganych pakietow R: ", paste(missing_packages, collapse = ", "))
}

# rstan kompiluje model C++. Na tej maszynie Rtools jest zainstalowany,
# ale nie zawsze jest widoczny w PATH podczas uruchomienia R z PowerShella.
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

library(rstan)
library(coda)
library(posterior)
library(ggplot2)
library(xtable)

# Parametry replikacji.
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

set.seed(DATA_SEED)

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

adult_data <- read.csv(
  file.path(data_dir, "adult.data"),
  header = FALSE,
  col.names = columns,
  na.strings = "?",
  strip.white = TRUE
)

adult_us <- subset(adult_data, native_country == "United-States")
adult_sample <- adult_us[sample(seq_len(nrow(adult_us)), N_SAMPLE), ]

adult_model <- adult_sample
adult_model$income_high <- ifelse(adult_model$income == ">50K", 1L, 0L)
adult_model$sex_male <- ifelse(adult_model$sex == "Male", 1L, 0L)
adult_model$married <- ifelse(
  adult_model$marital_status %in% c("Married-civ-spouse", "Married-spouse-absent"),
  1L,
  0L
)
adult_model$capital_gain_ind <- ifelse(adult_model$capital_gain > 0, 1L, 0L)

adult_model$education_z <- as.numeric(scale(adult_model$education_num))
adult_model$age_z <- as.numeric(scale(adult_model$age))
adult_model$age_z2 <- adult_model$age_z^2
adult_model$hours_z <- as.numeric(scale(adult_model$hours_per_week))

data_summary <- data.frame(
  wielkosc = c(
    "obserwacje w pliku adult.data",
    "obserwacje z USA",
    "wylosowana próba",
    "udział income_high = 1 w próbie"
  ),
  wartosc = c(
    nrow(adult_data),
    nrow(adult_us),
    nrow(adult_model),
    round(mean(adult_model$income_high), 3)
  ),
  stringsAsFactors = FALSE
)
write.csv(data_summary, file.path(tables_dir, "data_summary.csv"), row.names = FALSE)
print(
  xtable(data_summary, digits = c(0, 0, 3)),
  file = file.path(tables_dir, "data_summary.tex"),
  include.rownames = FALSE,
  floating = FALSE
)

X <- model.matrix(
  income_high ~ education_z + age_z + age_z2 + hours_z +
    sex_male + married + capital_gain_ind,
  data = adult_model
)
y <- as.integer(adult_model$income_high)

stopifnot(nrow(X) == N_SAMPLE)
stopifnot(all(y %in% c(0L, 1L)))

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

# Elicytacja priorow na skali ilorazow szans.
# Wartosci centralne OR sa przeksztalcane przez log(OR), poniewaz beta
# w logitcie oznacza zmiane logarytmu szans.
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
  central_OR = as.numeric(prior_or_location[names(prior_location)]),
  stringsAsFactors = FALSE
)

write.csv(prior_table, file.path(tables_dir, "prior_elicitation.csv"), row.names = FALSE)
print(
  xtable(prior_table, digits = c(0, 0, 0, 0, 3, 3, 3)),
  file = file.path(tables_dir, "prior_elicitation.tex"),
  include.rownames = FALSE,
  floating = FALSE
)

i_pred <- 1L
X_pred <- X[i_pred, , drop = FALSE]
y_true_pred <- y[i_pred]

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

options(mc.cores = min(MCMC_CHAINS, parallel::detectCores()))
rstan_options(auto_write = TRUE)

fit_rds_path <- file.path(outputs_dir, "fit_income_logit_student_t.rds")
if (file.exists(fit_rds_path)) {
  fit <- readRDS(fit_rds_path)
} else {
  fit <- stan(
    file = file.path(stan_dir, "income_logit_student_t.stan"),
    data = stan_data,
    seed = MCMC_SEED,
    chains = MCMC_CHAINS,
    iter = MCMC_ITER,
    warmup = MCMC_WARMUP,
    thin = MCMC_THIN,
    control = list(adapt_delta = ADAPT_DELTA, max_treedepth = MAX_TREEDEPTH)
  )
  saveRDS(fit, fit_rds_path)
}
saveRDS(adult_model, file.path(outputs_dir, "adult_model_sample.rds"))

pars_beta <- paste0("beta[", seq_len(ncol(X)), "]")
fit_summary <- summary(fit, pars = c(pars_beta, "p_pred"))$summary
write.csv(fit_summary, file.path(tables_dir, "stan_summary_raw.csv"))

draws <- rstan::extract(fit, permuted = TRUE, inc_warmup = FALSE)
beta_draws <- draws$beta
colnames(beta_draws) <- colnames(X)
or_draws <- exp(beta_draws)
colnames(or_draws) <- colnames(X)

hpdi <- function(x, prob = 0.95) {
  x <- sort(as.numeric(x))
  n <- length(x)
  m <- floor(prob * n)
  widths <- x[(m + 1):n] - x[seq_len(n - m)]
  i <- which.min(widths)
  c(lower = x[i], upper = x[i + m])
}

prior_density_at_zero <- function(j) {
  dt(
    (0 - prior_location[j]) / prior_scale[j],
    df = PRIOR_DF
  ) / prior_scale[j]
}

posterior_density_at_zero <- function(x) {
  dens <- density(x, adjust = 1.2, n = 4096)
  approx(dens$x, dens$y, xout = 0, rule = 2)$y
}

posterior_table <- data.frame(
  parameter = colnames(X),
  opis = unname(parameter_labels[colnames(X)]),
  prior_location = as.numeric(prior_location[colnames(X)]),
  posterior_mean = colMeans(beta_draws),
  posterior_median = apply(beta_draws, 2, median),
  posterior_sd = apply(beta_draws, 2, sd),
  HPDI_lower = apply(beta_draws, 2, function(x) hpdi(x)[1]),
  HPDI_upper = apply(beta_draws, 2, function(x) hpdi(x)[2]),
  P_beta_gt_0 = apply(beta_draws, 2, function(x) mean(x > 0)),
  P_beta_lt_0 = apply(beta_draws, 2, function(x) mean(x < 0)),
  stringsAsFactors = FALSE
)

posterior_table$OR_mean <- colMeans(or_draws)
posterior_table$OR_median <- apply(or_draws, 2, median)
posterior_table$OR_HPDI_lower <- apply(or_draws, 2, function(x) hpdi(x)[1])
posterior_table$OR_HPDI_upper <- apply(or_draws, 2, function(x) hpdi(x)[2])

posterior_table$prior_density_0 <- vapply(colnames(X), prior_density_at_zero, numeric(1))
posterior_table$posterior_density_0 <- apply(beta_draws, 2, posterior_density_at_zero)
posterior_table$BF_01 <- posterior_table$posterior_density_0 / posterior_table$prior_density_0
posterior_table$BF_10 <- 1 / posterior_table$BF_01

write.csv(posterior_table, file.path(tables_dir, "posterior_summary.csv"), row.names = FALSE)
print(
  xtable(
    posterior_table[, c(
      "parameter", "opis", "posterior_mean", "posterior_sd",
      "HPDI_lower", "HPDI_upper", "P_beta_gt_0", "BF_10"
    )],
    digits = c(0, 0, 0, 3, 3, 3, 3, 3, 2)
  ),
  file = file.path(tables_dir, "posterior_summary.tex"),
  include.rownames = FALSE,
  floating = FALSE
)

or_table <- posterior_table[, c(
  "parameter", "opis", "OR_mean", "OR_median", "OR_HPDI_lower", "OR_HPDI_upper"
)]
write.csv(or_table, file.path(tables_dir, "odds_ratio_summary.csv"), row.names = FALSE)
print(
  xtable(or_table, digits = c(0, 0, 0, 3, 3, 3, 3)),
  file = file.path(tables_dir, "odds_ratio_summary.tex"),
  include.rownames = FALSE,
  floating = FALSE
)

p_pred_draws <- as.numeric(draws$p_pred[, 1])
y_pred_rep <- as.integer(draws$y_pred_rep[, 1])
prediction_table <- data.frame(
  selected_observation = i_pred,
  observed_y = y_true_pred,
  posterior_mean_probability = mean(p_pred_draws),
  posterior_median_probability = median(p_pred_draws),
  HPDI_lower_probability = hpdi(p_pred_draws)[1],
  HPDI_upper_probability = hpdi(p_pred_draws)[2],
  predictive_probability_y1 = mean(y_pred_rep == 1L),
  stringsAsFactors = FALSE
)

prediction_covariates <- data.frame(
  variable = colnames(X),
  value = as.numeric(X_pred[1, ]),
  stringsAsFactors = FALSE
)

write.csv(prediction_table, file.path(tables_dir, "prediction_summary.csv"), row.names = FALSE)
write.csv(prediction_covariates, file.path(tables_dir, "prediction_covariates.csv"), row.names = FALSE)
print(
  xtable(prediction_table, digits = c(0, 0, 0, 3, 3, 3, 3, 3)),
  file = file.path(tables_dir, "prediction_summary.tex"),
  include.rownames = FALSE,
  floating = FALSE
)

diagnostic_summary <- summary(fit, pars = pars_beta)$summary
diagnostic_table <- data.frame(
  parameter = colnames(X),
  n_eff = diagnostic_summary[, "n_eff"],
  Rhat = diagnostic_summary[, "Rhat"],
  stringsAsFactors = FALSE
)
write.csv(diagnostic_table, file.path(tables_dir, "mcmc_diagnostics.csv"), row.names = FALSE)
print(
  xtable(diagnostic_table, digits = c(0, 0, 0, 3)),
  file = file.path(tables_dir, "mcmc_diagnostics.tex"),
  include.rownames = FALSE,
  floating = FALSE
)

sampler_params <- get_sampler_params(fit, inc_warmup = FALSE)
divergences <- sum(vapply(sampler_params, function(x) sum(x[, "divergent__"]), numeric(1)))
max_treedepth_hits <- sum(vapply(sampler_params, function(x) sum(x[, "treedepth__"] >= MAX_TREEDEPTH), numeric(1)))
accept_stat_mean <- mean(unlist(lapply(sampler_params, function(x) x[, "accept_stat__"])))

rep_draws <- draws$y_rep
observed_high_income_rate <- mean(y)
rep_high_income_rate <- rowMeans(rep_draws)

ppc_table <- data.frame(
  statistic = c(
    "sredni_udzial_y_1_w_probie",
    "sredni_udzial_y_1_w_replikacjach",
    "HPDI_dolny_udzial_y_1_w_replikacjach",
    "HPDI_gorny_udzial_y_1_w_replikacjach"
  ),
  value = c(
    observed_high_income_rate,
    mean(rep_high_income_rate),
    hpdi(rep_high_income_rate)[1],
    hpdi(rep_high_income_rate)[2]
  ),
  stringsAsFactors = FALSE
)
write.csv(ppc_table, file.path(tables_dir, "posterior_predictive_check.csv"), row.names = FALSE)

replication_info <- data.frame(
  parameter = c(
    "DATA_SEED", "MCMC_SEED", "N_SAMPLE", "MCMC_CHAINS", "MCMC_ITER",
    "MCMC_WARMUP", "MCMC_THIN", "ADAPT_DELTA", "MAX_TREEDEPTH",
    "PRIOR_DF", "divergences", "max_treedepth_hits", "mean_accept_stat"
  ),
  value = c(
    DATA_SEED, MCMC_SEED, N_SAMPLE, MCMC_CHAINS, MCMC_ITER,
    MCMC_WARMUP, MCMC_THIN, ADAPT_DELTA, MAX_TREEDEPTH,
    PRIOR_DF, divergences, max_treedepth_hits, accept_stat_mean
  ),
  stringsAsFactors = FALSE
)
write.csv(replication_info, file.path(tables_dir, "replication_info.csv"), row.names = FALSE)

capture.output(sessionInfo(), file = file.path(outputs_dir, "session_info.txt"))
capture.output(print(fit, pars = c(pars_beta, "p_pred")), file = file.path(outputs_dir, "stan_fit_print.txt"))

green_area <- rgb(24, 121, 104, 90, maxColorValue = 255)
green_line <- rgb(24, 121, 104, 255, maxColorValue = 255)
grey_area <- rgb(180, 180, 180, 80, maxColorValue = 255)
grey_line <- rgb(90, 90, 90, 255, maxColorValue = 255)

pdf(file.path(figures_dir, "posterior_marginals_beta.pdf"), width = 10, height = 8)
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
dev.off()

dt_scaled <- function(x, mean, scale, df) {
  dt((x - mean) / scale, df = df) / scale
}

pdf(file.path(figures_dir, "prior_vs_posterior_beta.pdf"), width = 10, height = 8)
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
    ylab = "gestosc",
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
dev.off()

or_plot_data <- subset(or_table, parameter != "(Intercept)")
or_plot_data$opis <- factor(or_plot_data$opis, levels = rev(or_plot_data$opis))
p_or <- ggplot(or_plot_data, aes(x = OR_median, y = opis)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = OR_HPDI_lower, xmax = OR_HPDI_upper),
                 height = 0.18, color = green_line, linewidth = 0.7) +
  geom_point(size = 2.4, color = green_line) +
  scale_x_log10() +
  labs(x = "iloraz szans exp(beta), skala logarytmiczna", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(figures_dir, "odds_ratios_hpdi.pdf"), p_or, width = 8, height = 4.8)

trace_pdf <- file.path(figures_dir, "traceplot_beta.pdf")
pdf(trace_pdf, width = 10, height = 8)
rstan::traceplot(fit, pars = pars_beta, inc_warmup = FALSE, nrow = 4)
dev.off()

posterior_prob_plot_data <- posterior_table
posterior_prob_plot_data$opis <- factor(posterior_prob_plot_data$opis, levels = rev(posterior_prob_plot_data$opis))
p_prob <- ggplot(posterior_prob_plot_data, aes(x = P_beta_gt_0, y = opis)) +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_point(size = 2.2, color = green_line) +
  xlim(0, 1) +
  labs(x = "P(beta > 0 | y)", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank())
ggsave(file.path(figures_dir, "posterior_probability_positive.pdf"), p_prob, width = 8, height = 4.8)

pred_df <- data.frame(p = p_pred_draws)
p_pred <- ggplot(pred_df, aes(x = p)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40, fill = green_area, color = "white") +
  geom_density(color = green_line, linewidth = 0.9) +
  geom_vline(xintercept = mean(p_pred_draws), color = green_line, linewidth = 0.8) +
  geom_vline(xintercept = hpdi(p_pred_draws), linetype = "dashed", color = "grey35") +
  labs(x = "P(income > 50K | x, y)", y = "gestosc") +
  theme_minimal(base_size = 11)
ggsave(file.path(figures_dir, "prediction_probability_selected_observation.pdf"), p_pred, width = 7, height = 4.6)

ppc_df <- data.frame(rep_high_income_rate = rep_high_income_rate)
p_ppc <- ggplot(ppc_df, aes(x = rep_high_income_rate)) +
  geom_histogram(aes(y = after_stat(density)), bins = 35, fill = grey_area, color = "white") +
  geom_vline(xintercept = observed_high_income_rate, color = green_line, linewidth = 1) +
  labs(x = "udzial obserwacji z income_high = 1", y = "gestosc predykcyjna") +
  theme_minimal(base_size = 11)
ggsave(file.path(figures_dir, "posterior_predictive_check_income_rate.pdf"), p_ppc, width = 7, height = 4.6)

cat("Analiza zakonczona.\n")
cat("PDF-y wykresow zapisano w:", figures_dir, "\n")
cat("Tabele zapisano w:", tables_dir, "\n")
