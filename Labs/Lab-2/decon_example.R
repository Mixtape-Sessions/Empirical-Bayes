# This code is called from within Stata

# Libraries
library(foreign)
library(tidyverse)
library(devtools)
library(ggtext)
library(ggplot2)
library(gtools)
library(deconvolveR)

# Set working directory
setwd("/Users/christopherwalters/Dropbox/teaching/short_courses/mixtape/assignments/lab2")

# Load data
data <- read.csv("disc_estimates.csv", header = F)
id_f <- unlist(data[1])
thetahat_f <- unlist(data[2])
s_f <- unlist(data[3])

# Construct z-scores
z_f <- thetahat_f / s_f

# Choose tuning parameters

# Spline order
p <- 5

# Support (limits and number of points)
z_low <- min(z_f)
z_high <- max(z_f)
n_pts <- 200
width <- (z_high - z_low) / n_pts
supp_z <- seq(z_low, z_high, width)
theta_low <- min(z_low * s_f)
theta_high <- max(z_high * s_f)
width_theta <- (theta_high - theta_low) / n_pts
supp_theta <- seq(theta_low, theta_high, width_theta)

# MLE penalization parameter
c <- 0.05

# Deconvolution
decon_results <- deconv(
  tau = supp_z, X = z_f, family = "Normal",
  c0 = c, pDegree = p
)
g_mu <- decon_results$stats[, "g"] / width
mus <- decon_results$stats[, "theta"]

# Convert from z-scores to levels of theta, assuming independence of mu and s

# Estimate density of log(s)
sdens <- stats::density(log(s_f, base = exp(1)), n = 2^10) # this integrates to 1

# Change of variables to go from g(z) and f(s) to g(z*s)
deltadens <- function(delta) {
  support <- delta / exp(sdens$x)
  fz <- sapply(support, function(x) {
    if (x < min(supp_z)) {
      return(0) # assume density is zero outside support for mus
    } else if (x > max(supp_z)) {
      return(0) # assume density is zero outside support for mus
    }

    g_mu[which.min(abs(supp_z - x))]
  })
  inside <- exp(-sdens$x) * fz * sdens$y
  dens <- sum(inside) / sum(sdens$y)
  return(dens)
}
g_theta <- sapply(supp_theta, deltadens)

#Posterior means
postmean <- function(z) {
  postdens <- dnorm(z-mus)*g_mu
  return(sum(postdens*mus)/sum(postdens))
}
post_mean_mu_f <- sapply(z_f, postmean)
post_mean_theta_f <- post_mean_mu_f*s_f

# Save results
output_g <- data.frame(supp_theta, g_theta)
write.table(output_g, "decon_estimates.csv", sep = ",", col.names = F, row.names = F)
output_postmean <- data.frame(id_f,post_mean_theta_f)
write.table(output_postmean, "postmean_estimates.csv", sep = ",", col.names = F, row.names = F)
