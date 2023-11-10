# ------------------------------------------------------------------------------
# Kyle Butts & Chris Walters
# 10/2023
# This program simulates data from a school value-added model
# and estimates the model parameters using empirical Bayes methods

library(fixest)
library(data.table)
library(broom)
library(glue)
library(here)
library(tidyverse)
set.seed(1028)

# Change:
setwd(here::here("Labs/Lab-1"))

n_students = 2500
n_schools = 50
sigma_y = 1
sigma_theta = 0.2
df = as.data.table(read_csv("vam_example_data.csv"))

# Estimate value-added models --------------------------------------------------

### No controls ----------------------------------------------------------------

est_unc = feols(
  Y ~ 0 + i(D), df, vcov = "hc1"
)

# Extract value-added estimates
va_unc = as.data.table(broom::tidy(est_unc))
va_unc = va_unc[grepl("D::", va_unc$term), ]
va_unc$D = as.numeric(gsub("D::", "", va_unc$term))

va_unc = va_unc[, .(D, thetahat_unc = estimate, se_unc = std.error)]

### With control for X ---------------------------------------------------------

est_con = feols(
  Y ~ 0 + i(D) + X, df, vcov = "hc1"
)

va_con = as.data.table(broom::tidy(est_con))
va_con = va_con[grepl("D::", va_con$term), ]
va_con$D = as.numeric(gsub("D::", "", va_con$term))

va_con = va_con[, .(D, thetahat_con = estimate, se_con = std.error)]

## Collapse to school-level ----------------------------------------------------

# Get number of students and true value added for each school
school_va = df[, .(N = .N, theta = theta_D[1]), by = "D"]

# Merge in estimates
school_va = merge(school_va, va_unc, by = "D")
school_va = merge(school_va, va_con, by = "D")

## Estimate hyperparameters for each model -------------------------------------

#' @param thetahat Vector of school-level estimates
#' @param se Vector of standard errors of school-level estimates
#' @return List with estimated hyperparameters `muhat` and `sigmahat`
est_hyperparameters = function(thetahat, se) {
  muhat = mean(thetahat)
  c = (thetahat - muhat)^2 - se^2

  c_mean = mean(c)
  # If mean(c) < 0, then the we should set c = 0
  c_mean = max(c_mean, 0)

  sigmahat = sqrt(c_mean)

  return(list(muhat = muhat, sigmahat = sigmahat))
}

hp_unc = est_hyperparameters(school_va$thetahat_unc, school_va$se_unc)
hp_con = est_hyperparameters(school_va$thetahat_con, school_va$se_con)

# Display results
glue("Estimated hyperparameters for model `unc`: mean = {round(hp_unc$muhat, 2)}, SD = {round(hp_unc$sigmahat, 3)}")
glue("Estimated hyperparameters for model `con`: mean = {round(hp_con$muhat, 2)}, SD = {round(hp_con$sigmahat, 3)}")

## Compute linear shrinkage estimates ------------------------------------------

#' @param thetahat Vector of school-level estimates
#' @param se Vector of standard errors of school-level estimates
#' @param muhat Estimated mean of school-level estimates (hyperparameter)
#' @param sigmahat Estimated SD of school-level estimates (hyperparameter)
#' @return Vector of linear shrinkage estimates
linear_shrinkage <- function(thetahat, se, muhat, sigmahat) {
  # Shrinkage factor 
  lambda = sigmahat^2 / (sigmahat^2 + se^2)

  # Linear shrinkage formula
  thetastar = lambda * thetahat + (1 - lambda) * muhat
}

school_va$thetastar_unc = linear_shrinkage(
  school_va$thetahat_unc, school_va$se_unc, 
  hp_unc$muhat, hp_unc$sigmahat
)

school_va$thetastar_con = linear_shrinkage(
  school_va$thetahat_con, school_va$se_con, 
  hp_con$muhat, hp_con$sigmahat
)

## Assess performance of estimators --------------------------------------------

# Compare dispersion of unbiased estimates, estimated prior, shrunk posteriors, and true parameters
SD_true = sd(school_va$theta)

SD_unc_1 = sd(school_va$thetahat_unc)
SD_unc_2 = sd(school_va$thetastar_unc)

glue("
  Model `unc`: 
    SD estimates = {round(SD_unc_1, 3)}, 
    SD prior = {hp_unc$sigmahat}, 
    SD posteriors = {SD_unc_2}, 
    SD truth = {SD_true}
")

SD_con_1 = sd(school_va$thetahat_con)
SD_con_2 = sd(school_va$thetastar_con)
glue("
  Model `con`: 
    SD estimates = {round(SD_con_1, 3)}, 
    SD prior = {hp_con$sigmahat}, 
    SD posteriors = {SD_con_2}, 
    SD truth = {SD_true}
")

# Compare MSE of each estimate
calc_mse = function(est, true) mean((est - true)^2)

MSE_thetahat_unc = calc_mse(school_va$thetahat_unc, school_va$theta)
MSE_thetastar_unc = calc_mse(school_va$thetastar_unc, school_va$theta)
MSE_thetahat_con = calc_mse(school_va$thetahat_con, school_va$theta)
MSE_thetastar_con = calc_mse(school_va$thetastar_con, school_va$theta)

glue("
  MSE thetahat_unc:  {round(MSE_thetahat_unc, 3)}
  MSE thetastar_unc: {round(MSE_thetastar_unc, 3)}
  MSE thetahat_con:  {round(MSE_thetahat_con, 3)}
  MSE thetastar_con: {round(MSE_thetastar_con, 3)}
")


## Two-way plot of truth and estimate ------------------------------------------

# For controlled model, graph true VA against unbiased estimates and linear shrinkage estimates, along w/45-deg line

# Make plot square
extent = max(abs(range(school_va$theta, school_va$thetahat_con, school_va$thetastar_con)))
limits = c(-extent, extent) * 1.1

ggplot(school_va) + 
  geom_point(
    aes(x = thetahat_con, y = theta, color = "thetahat"),
    size = 2
  ) +
  geom_point(
    aes(x = thetastar_con, y = theta, color = "thetastar"),
    size = 2
  ) +
  geom_abline(
    slope = 1, intercept = 0, 
    color = "grey", linetype = "dashed"
  ) +
  scale_color_manual(
    values = c("thetahat" = "navy", "thetastar" = "maroon"), 
    labels = c("thetahat" = "Unbiased Estimates", "thetastar" = "Linear Shrinkage Estimates")
  ) +
  # Make scale of x and y axes the same
  scale_x_continuous(
    limits = limits
  ) + 
  scale_y_continuous(
    limits = limits
  ) +
  labs(
    x = "Estimated value-added",
    y = "True value-added",
    color = NULL
  ) +
  theme_bw(base_size = 20) + 
  theme(
    legend.position = "bottom"
  )

## Histogram of estimates ------------------------------------------------------

width = 0.06
binwidth = 0.06 / 1.5
prior_dist_func = function(x) {
  n_schools * width * dnorm(x, hp_con$muhat, hp_con$sigmahat)
}
hist_limits = c(-4 * sigma_theta, 4 * sigma_theta)

ggplot(school_va) + 
  geom_histogram(
    aes(x = thetahat_con, color = "thetahat", fill = "thetahat"),
    binwidth = binwidth
  ) +
  geom_histogram(
    aes(x = thetastar_con, color = "thetastar", fill = "thetastar"),
    binwidth = binwidth
  ) +
  # Plot prior distribution
  stat_function(
    fun = prior_dist_func, 
    color = "grey20", linewidth = 1
  ) +
  scale_color_manual(
    values = c("thetahat" = "navy", "thetastar" = "maroon"), 
    labels = c("thetahat" = "Unbiased Estimates", "thetastar" = "Linear Shrinkage Estimates")
  ) +
  scale_fill_manual(
    values = c("thetahat" = "#ffffff00", "thetastar" = "#b0306055"), 
    labels = c("thetahat" = "Unbiased Estimates", "thetastar" = "Linear Shrinkage Estimates")
  ) +
  scale_x_continuous(limits = hist_limits) +
  labs(
    x = "Estimated value-added",
    y = "True value-added",
    title = "Plot of value-added estimates and prior distribution",
    color = NULL, fill = NULL
  ) +
  theme_bw(base_size = 20) + 
  theme(
    legend.position = "bottom"
  )



