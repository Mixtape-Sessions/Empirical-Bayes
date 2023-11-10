#' @param n_students Number of students
#' @param n_schools Number of schools
#' @param sigma_y Residual SD of outcome
#' @param sigma_x SD of covariate X
#' @param beta Coefficient on X in outcome equation
#' @param va_model Parametric form for VA distribution (normal by default, or lognormal)
#' @param mu_theta Mean of school value-added
#' @param sigma_theta SD of school value-added
#' @param sigma_delta Dispersion of school mean utilities
#' @param sigma_gamma Dispersion of school utility coefficients on student covariate X
#' @return Data.frame consisting of student-level data
simulate_data <- function(
  n_students = 2500, n_schools = 50, sigma_y = 1,
  sigma_x = 1, beta = 1, 
  va_model = "normal", mu_theta = 0, sigma_theta = 0.2,
  sigma_delta = 1, sigma_gamma = 0.5
) {
  ## Draw school-level parameters ----------------------------------------------

  # Student-by-school dataset
  df = data.table::CJ(student_id = 1:n_students, D = 1:n_schools)

  # Value-added
  df[, 
    theta_D := rnorm(1, mu_theta, sigma_theta), 
    by = D
  ]
  if (va_model == "lognormal") {
    df[, theta_D := exp(theta_D)]
  }

  # Utility parameters
  df[, 
    `:=`(
      delta = rnorm(1, 0, sigma_delta), 
      gamma = rnorm(1, 0, sigma_gamma)
    ), 
    by = D
  ]

  ## Generate student-level data -----------------------------------------------

  # Student covariate
  df[, X := rnorm(1, 0, sigma_x), by = student_id]

  # For each student, calculate utility for each school
  df[, Uj := delta + gamma * X - log(-log(runif(.N)))]

  # Draw the highest utility school
  # Order schools by utility
  setorder(df, student_id, -Uj)
  # Select max U for each student i
  df = df[, .SD[1], by = student_id]

  # Generate student outcomes
  df[, Y := theta_D + beta * X + rnorm(.N, 0, sigma_y)]

  ## Clean-up data and return
  df$Uj = df$delta = df$gamma = NULL
  return(df)
}

# Simulate data ----------------------------------------------------------------
library(fixest)
library(data.table)
library(broom)
library(glue)
library(here)
library(tidyverse)
set.seed(1028)

# Change:
setwd(here::here("Labs/Lab-1"))

df = simulate_data()
write_csv(df, "vam_example_data.csv")

