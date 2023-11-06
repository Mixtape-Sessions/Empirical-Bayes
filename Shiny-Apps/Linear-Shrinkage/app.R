## -----------------------------------------------------------------------------
## app.R
## Kyle Butts, CU Boulder Economics
##
## This is a shiny application to build intuition about linear shrinkage models
## -----------------------------------------------------------------------------

library(shiny)
library(shiny.tailwind)
library(fixest)
library(data.table)
library(broom)
library(glue)
library(tidyverse)
set.seed(1028)

# setwd("Shiny-Apps/Linear-Shrinkage")
source("components.R")

ui <- page(
  mid_wrapper(
    shiny.tailwind::use_tailwind(
      css = "style.css",
      tailwindConfig = "tailwind.config.js"
    ),
    shiny::withMathJax(),
    # Header
    gradient_header("Linear Shrinkage Estimators", from = "#4b6cb7", to = "#182848"),
    div(
      class = "mt-6 prose lg:prose-lg xl:prose-xl",
      p("Play around the simulation parameters below to see how the linear shrinkage estimates performs compared to standard fixed-effect estimates."),
    ),
    border_container(
      class = "my-4", 
      # Inputs 
      div(class = "hide_grid my-4 grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 md:gap-x-2 gap-y-4",
        div(
            caps_header("Number of Observations"),
            sliderInput("N_students", "Number of Students:", 2500, min = 100, max = 10000),
            sliderInput("N_schools", "Number of Schools:", 50, min = 10, max = 100),
        ),
        div(
          caps_header("Student-level parameters"),
          sliderInput("sigma_x", "SD of X:", 1, min = 0.01, max = 5),
          sliderInput("beta", "Coefficient on X in outcome equation:", 1, min = -3, max = 3),
          sliderInput("sigma_y", "Residual SD of outcome:", 1, min = 0.1, max = 2),
        ),
        div(
          caps_header("School-level parameters"),
          sliderInput("mu_theta", "Mean of school value-added:", 0, min = -3, max = 3),
          sliderInput("sigma_theta", "SD of school value-added:", 0.2, min = 0.01, max = 1),
          sliderInput("sigma_delta", "Dispersion of school mean utilities:", 1, min = 0.01, max = 5),
          sliderInput("sigma_gamma", "Dispersion of school utility coefficients on student covariate X:", 0.5, min = 0.01, max = 2.5),
        ),
      )
    ),
    border_container(
      class = "my-4",
      div(
        caps_header("Plot of value-added estimators"),
        shiny::plotOutput("va_plot", width = "100%", height = "100%")
      ),
      div(
        class = "my-4",
        caps_header("Histogram of value-added estimators"),
        shiny::plotOutput("va_hist", width = "100%", height = "100%")
      ),
      div(
        class = "my-4",
        caps_header("Performance of value-added estimators"),
        div(class = "grid grid-cols-1 gap-y-4 md:grid-cols-2 md:gap-x-4 ",
          uiOutput("va_sd"), 
          uiOutput("va_mse") 
        )
      ), 
    ),
  ),
)

gen_data <- function(N_students, N_schools, mu_theta, sigma_theta, sigma_delta, sigma_gamma, sigma_x, beta, sigma_y) {
  ## Draw school-level parameters ----------------------------------------------
  # Student-by-school dataset
  df = data.table::CJ(i = 1:N_students, j = 1:N_schools)

  # Value-added
  df[, 
    theta := rnorm(1, mu_theta, sigma_theta), 
    by = "j"
  ]

  # Utility parameters
  df[, 
    `:=`(
      delta = rnorm(1, 0, sigma_delta), 
      gamma = rnorm(1, 0, sigma_gamma)
    ), 
    by = "j"
  ]

  ## Generate student-level data -----------------------------------------------
  # Student covariate
  df[, X := rnorm(1, 0, sigma_x), by = "i"]

  # For each student, calculate utility for each school
  df$Uj = df$delta + df$gamma * df$X - log(-log(runif(nrow(df))))

  # Draw the highest utility school
  # Order schools by utility
  setorderv(df, c("i", "Uj"), c(1, -1))
  # Select max U for each student i
  df = df[, .SD[1], by = "i"]

  # Generate student outcomes
  df$Y = df$theta + beta * df$X + rnorm(nrow(df), 0, sigma_y)

  df
}

est_school_va <- function(df) {
  # Estimate model -------------------------------------------------------------

  ## Value-added models --------------------------------------------------------

  ### No controls --------------------------------------------------------------

  est_unc = feols(
    Y ~ 0 + i(j), data = df, vcov = "hc1"
  )

  # Extract value-added estimates
  va_unc = as.data.table(broom::tidy(est_unc))
  va_unc = va_unc[grepl("j::", va_unc$term), ]
  va_unc$j = as.numeric(gsub("j::", "", va_unc$term))

  va_unc = va_unc[, .(j, thetahat_unc = estimate, se_unc = std.error)]

  ### With control for X -------------------------------------------------------

  est_con = feols(
    Y ~ 0 + i(j) + X, data = df, vcov = "hc1"
  )

  va_con = as.data.table(broom::tidy(est_con))
  va_con = va_con[grepl("j::", va_con$term), ]
  va_con$j = as.numeric(gsub("j::", "", va_con$term))

  va_con = va_con[, .(j, thetahat_con = estimate, se_con = std.error)]

  ## Collapse to school-level --------------------------------------------------

  # Get number of students and true value added for each school
  school_va = df[, .(N = .N, theta = theta[1]), by = j]

  # Merge in estimates
  school_va = merge(school_va, va_unc, by = "j")
  school_va = merge(school_va, va_con, by = "j")
  school_va

  ## Estimate hyperparameters for each model -----------------------------------

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

  ## Compute linear shrinkage estimates ----------------------------------------

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

  school_va
}

# set.seed(1028)
# (df = gen_data(2500, 50, 0, 0.2, 1, 0.5, 1, 1, 1))
# est_school_va(df)

server <- function(input, output, session) {
  set.seed(1028)
  
  df = reactive({ 
    gen_data(
      input$N_students, input$N_schools, input$mu_theta, input$sigma_theta, input$sigma_delta, input$sigma_gamma, input$sigma_x, input$beta, input$sigma_y
    )
  })

  value_added_ests = reactive({
    est_school_va(df())
  })

  output$va_sd <- shiny::renderUI({
    school_va = value_added_ests()

    #' @param thetahat Vector of school-level estimates
    #' @param se Vector of standard errors of school-level estimates
    #' @return List with estimated hyperparameters `muhat` and `sigmahat`
    est_hyperparameters = function(thetahat, se) {
      muhat = mean(thetahat)
      c = (thetahat - muhat)^2 - se^2
      sigmahat = sqrt(mean(c))

      return(list(muhat = muhat, sigmahat = sigmahat))
    }

    hp_unc = est_hyperparameters(school_va$thetahat_unc, school_va$se_unc)
    hp_con = est_hyperparameters(school_va$thetahat_con, school_va$se_con)

    SD_true = sd(school_va$theta)

    SD_unc_1 = sd(school_va$thetahat_unc)
    SD_unc_2 = sd(school_va$thetastar_unc)
    SD_con_1 = sd(school_va$thetahat_con)
    SD_con_2 = sd(school_va$thetastar_con)

    list(
      tags$div(
        tags$h3(class = "font-semibold", "Without covariate"),
        tags$ul(
          class = "list-disc list-inside",
          tags$li("SD estimates", round(SD_unc_1, 3)),
          tags$li("SD prior", hp_unc$sigmahat),
          tags$li("SD posteriors", round(SD_unc_2, 3)),
          tags$li("SD truth", round(SD_true, 3))
        )
      ),
      tags$div(
        tags$h3(class = "font-semibold my-4", "With covariate"),
        tags$ul(
          class = "list-disc list-inside",
          tags$li("SD estimates", round(SD_con_1, 3)),
          tags$li("SD prior", hp_con$sigmahat),
          tags$li("SD posteriors", round(SD_con_2, 3)),
          tags$li("SD truth", round(SD_true, 3))
        )
      )
    )

  })

  output$va_mse <- shiny::renderUI({
    school_va = value_added_ests()

    # Compare MSE of each estimate
    calc_mse = function(est, true) mean((est - true)^2)

    MSE_thetahat_unc = calc_mse(school_va$thetahat_unc, school_va$theta)
    MSE_thetastar_unc = calc_mse(school_va$thetastar_unc, school_va$theta)
    MSE_thetahat_con = calc_mse(school_va$thetahat_con, school_va$theta)
    MSE_thetastar_con = calc_mse(school_va$thetastar_con, school_va$theta)

    tags$div(
      tags$h3(class = "font-semibold", "MSE of value-added estimates"),
      tags$ul(
        class = "list-disc list-inside",
        tags$li("Unbiased estimates: ", round(MSE_thetahat_unc, 3)),
        tags$li("Linear shrinkage estimates: ", round(MSE_thetastar_unc, 3)),
        tags$li("Unbiased estimates: ", round(MSE_thetahat_con, 3)),
        tags$li("Linear shrinkage estimates: ", round(MSE_thetastar_con, 3))
      )
    )
  })

  output$va_plot <- shiny::renderPlot(
    {
      school_va = value_added_ests()

      # For controlled model, graph true VA against unbiased estimates and linear shrinkage estimates, along w/45-deg line
      # Make plot square
      extent = max(abs(range(school_va$theta, school_va$thetahat_con, school_va$thetastar_con)))
      limits = c(-extent, extent) * 1.1

      # Plot
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
    },
    # https://github.com/rstudio/shiny/issues/650
    height = function() {
      session$clientData$output_va_plot_width * 7/16
    }, 
    res = 72,
  )

  output$va_hist <- shiny::renderPlot(
    {
      school_va = value_added_ests()

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
      
      width = 0.06
      binwidth = 0.06 / 1.5
      prior_dist_func = function(x) {
        input$N_schools * width * dnorm(x, hp_con$muhat, hp_con$sigmahat)
      }
      hist_limits = c(-4 * input$sigma_theta, 4 * input$sigma_theta)

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
    },
    # https://github.com/rstudio/shiny/issues/650
    height = function() {
      session$clientData$output_va_hist_width * 7/16
    }, 
    res = 72,
  )



}

shinyApp(ui = ui, server = server)


