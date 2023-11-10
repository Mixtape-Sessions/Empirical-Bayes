# ------------------------------------------------------------------------------
# Kyle Butts & Chris Walters
# 10/2023
# This program applies empirical Bayes methods
# to study employment discrimination using data
# from Kline, Rose, and Walters (QJE 2022)

library(fixest)
library(data.table)
library(tidyverse)
library(here)
library(glue)
library(deconvolveR)

# Switches
estimate_gaps = 1
multiple_test = 1

# Load and format data ---------------------------------------------------------
data = read_csv(here("Labs/Lab-2/krw_data.csv"))

# Summarize data
summary(data)

# Overall treatment effects (robust and clustered ses)
est = feols(callback ~ i(white), data)

esttable(
  summary(est, vcov = "hc1"),
  summary(est, cluster = ~job)
)

# Collapse to data set of firm-specific gaps and SEs
firm = data |> 
  # Collapse per job-opening
  summarize(
    ybar_white = mean(callback[white == 1]),
    ybar_black = mean(callback[white == 0]),
    n_j = n(),
    firm = firm[1],
    .by = "job"
  ) |> 
  mutate(thetahat_jf = ybar_white - ybar_black) |>
  ungroup() |>
  # Collapse per firm
  filter(!is.na(thetahat_jf)) |> 
  summarize(
    thetahat_f = mean(thetahat_jf),
    std = sd(thetahat_jf),
    J = n(),
    .by = "firm"
  ) |> 
  mutate(
    s_f = sqrt((1/J) * std^2)
  ) |> 
  select(firm, thetahat_f, s_f)

# Estimate distribution of race gaps -------------------------------------------
if (estimate_gaps == 1) {

  # Estimate variance of mixing distribution
  mu = mean(firm$thetahat_f)
  N_firms = length(firm$thetahat_f)

  v = (firm$thetahat_f - mu)^2 - 
        (((N_firms - 1) / N_firms) * (firm$s_f^2))

  sigma = sqrt(mean(v))
  glue::glue("Estimated SD of firm gaps: {round(sigma, 3)}")

  # Linear shrinkage estimates
	firm$lambda = sigma^2 / (sigma^2 + firm$s_f^2)
	firm$thetastar = 
    firm$lambda * firm$thetahat_f + 
    (1 - firm$lambda) * mu

  # Histogram of estimates vs. linear shrinkage posteriors
  ggplot(firm) + 
    geom_histogram(
      aes(x = thetahat_f, color = "thetahat", fill = "thetahat"),
      size = 1.2
    ) +
    geom_histogram(
      aes(x = thetastar, color = "thetastar", fill = "thetastar"),
      size = 1.2
    ) +
    scale_color_manual(
      values = c("thetahat" = "navy", "thetastar" = "maroon"), 
      labels = c("thetahat" = "Unbiased Estimates", "thetastar" = "Linear Shrinkage Estimates")
    ) +
    scale_fill_manual(
      values = c("thetahat" = "#ffffff00", "thetastar" = "#b0306055"), 
      labels = c("thetahat" = "Unbiased Estimates", "thetastar" = "Linear Shrinkage Estimates")
    ) +
    labs(
      x = "White/Black contact gap",
      y = "Number of firms",
      color = NULL, fill = NULL
    ) +
    theme_bw(base_size = 20) + 
    theme(
      legend.position = "bottom"
    )

  # Nonparametric deconvolution ----
  
    # construct z-scores
    firm$z_f = firm$thetahat_f / firm$s_f

    ## Choose tuning parameters ----
    # Spline order
    p = 5

    # Support (limits and number of points)
    n_pts = 200
    z_low = min(firm$z_f)
    z_high = max(firm$z_f)
    width = (z_high - z_low) / n_pts
    supp_z = seq(z_low, z_high, width)
    
    theta_low = min(z_low * firm$s_f)
    theta_high = max(z_high * firm$s_f)
    width_theta= (theta_high - theta_low) / n_pts
    supp_theta = seq(theta_low, theta_high, width_theta)

    # MLE penalization parameter
    c = 0.1

    ## Deconvolution ----
    decon_results <- deconv(
      tau = supp_z, X = firm$z_f, family = "Normal", c0 = c, pDegree = p
    )
    g_mu <- decon_results$stats[, "g"] / width
    mus <- decon_results$stats[, "theta"]

  # Convert from z-scores to levels of theta, assuming independence of mu and s ----

    # Estimate density of log(s) (this integrates to 1)
    sdens <- stats::density(log(firm$s_f, base = exp(1)), n = 2^10) 

    # Change of variables to go from g(z) and f(s) to g(z * s)
    deltadens <- function(delta) {
      support <- delta / exp(sdens$x)
            
      fz <- sapply(support, function(x) {
        if (x < min(supp_z)) {
          return(0) # assume density is zero outside support for mus
        }
        else if (x > max(supp_z)) {
          return(0) # assume density is zero outside support for mus
        }
        
        g_mu[which.min(abs(supp_z - x))]
      })

      inside <- exp(-sdens$x) * fz * sdens$y
      dens <- sum(inside) / sum(sdens$y)
      return(dens)
    }
    
    g_theta <- sapply(supp_theta, deltadens)  
    decon_est = data.frame(supp_theta, g_theta)

  # Plot deconvolution results ----
  ggplot() + 
    geom_histogram(
      aes(x = thetahat_f, color = "thetahat"),
      data = firm, 
      fill = NA, size = 1.2,
      key_glyph = draw_key_smooth
    ) + 
    geom_line(
      aes(x = supp_theta, y = g_theta, color = "np"),
      data = decon_est, size = 1.5
    ) +    
    labs(
      x = "White/Black contact gap",
      y = "Density",
      color = NULL
    ) +
    scale_color_manual(
      values = c("thetahat" = "navy", "np" = "maroon"), 
      labels = c("thetahat" = "Unbiased Estimates", "np" = "Nonparametric Prior"),
    ) +
    theme_bw(base_size = 20) + 
    theme(
      legend.position = "bottom"
    )

}

# Multiple testing: Classify firms controlling FDR -----------------------------
if (multiple_test == 1) { 

  # Choose cutoff lambda for calculating bound on share of true nulls (pi0)
  lambda = 0.5
  
  # Choose False Discovery Rate control threshold
  FDR = 0.05

  # Compute one-tailed tests of H0: theta = 0 vs. HA: theta > 0
  firm$z_f = firm$thetahat_f / firm$s_f
  firm$p_f = 1 - pnorm(firm$z_f)
        
  # Bound pi_0
  c = (firm$p_f * (firm$p_f > lambda)) / (1 - lambda)
  pi_0 = mean(c)
  glue::glue("Bound on pi_0: {round(pi_0, 3)}")
    
  # Get CDF of p-value
  firm = firm[order(firm$p_f), ]
  firm$F_p = order(firm$p_f) / nrow(firm)
  
  # Compute q-values
  firm$q_f = (firm$p_f * pi_0) / firm$F_p
    
  # Count firms with q-vals below FDR threshold, and find p-val cutoff
  below_FDR = firm$q_f < FDR
  N_firms_below = sum(below_FDR)
  p_cutoff = max(firm[below_FDR, ]$p_f)
					
  # Plot histogram of p-vals with lambda, pi0, and p-val cutoff for FDR control
  pi0 = round(pi_0, 2)
  lambda_loc = lambda + 0.005
  pi0_loc = pi0 + 0.25
  p_loc = p_cutoff + 0.01
  breaks = seq(0, 1, 0.05)

  ggplot() + 
    geom_histogram(
      aes(x = p_f, y = after_stat(density)), 
      data = firm, breaks = breaks, 
      color = "navy", fill = NA, size = 1.2
    ) +
    # Pi_0
    geom_hline(
      yintercept = pi0, 
      color = "maroon", linetype = "dashed"
    ) + 
    annotate(
      "text", label = as.character(glue("pi_0 = {pi0}")),
      x = 0.725, y = pi0_loc, hjust = 0.5, size = 6
    ) +
    # lambda
    geom_vline(
      xintercept = lambda,
      linetype = "dashed"
    ) + 
    annotate(
      "text", label = as.character(glue("lambda = {lambda}")),
      x = lambda_loc, y = 3.5, hjust = 0, size = 6
    ) +
    # p_cutoff
    geom_vline(
      xintercept = p_cutoff,
      linetype = "dashed"
    ) +
    annotate(
      "text", label = as.character(glue("{N_firms_below} firms with q-vals < {FDR}")),
      x = p_loc, y = 5.5, hjust = 0, size = 6
    ) +
    labs(
      x = "P-value from test of no disc. against Black applicants",
      y = "Density"
    ) +
    theme_bw(base_size = 20) 
}




