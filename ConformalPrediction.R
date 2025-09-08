## ======================================================================
## Distributional Conformal Prediction (DCP) for Markov(p=1)
## ======================================================================

#-----------------------------------------------------------------------
# 1. Package Dependencies and Setup
#-----------------------------------------------------------------------

# Check and install packages if not present
if (!require("VGAM")) install.packages("VGAM")
if (!require("dplyr")) install.packages("dplyr")
if (!require("knitr")) install.packages("knitr")
if (!require("doParallel")) install.packages("doParallel")
if (!require("foreach")) install.packages("foreach")

suppressPackageStartupMessages({
  library(VGAM)
  library(dplyr)
  library(knitr)
  library(doParallel)
  library(foreach)
})

# Set a random seed for reproducibility
set.seed(123)

#-----------------------------------------------------------------------
# 2. Kernel and CDF Helper Functions
#-----------------------------------------------------------------------

# Smoothed Heaviside function using clamped Normal CDF
kernel_lambda <- function(x) {
  clamped_x <- pmax(-2, pmin(2, x))
  pnorm(clamped_x)
}

# Standard Gaussian kernel
kernel_K <- function(u) {
  dnorm(u)
}

# Conditional CDF estimator: F̂(x | y) for Markov(p=1)
estimate_conditional_cdf <- function(x_val, y_cond, h, series, h0 = h^2) {
  m <- length(series)
  if (m < 2) stop("Series is too short for p=1.")
  
  x_train <- series[2:m]
  y_train <- series[1:(m-1)]
  
  weights <- kernel_K((y_cond - y_train) / h)
  sum_weights <- sum(weights)
  
  if (!is.finite(sum_weights) || sum_weights < 1e-12) {
    return(mean(x_train <= x_val))
  }
  
  numerator <- sum(weights * kernel_lambda((x_val - x_train) / h0))
  numerator / sum_weights
}

# Bandwidth selection via Kolmogorov-Smirnov test p-value maximization
select_bandwidth_ks <- function(series, h_grid) {
  m <- length(series)
  if (m < 3) stop("Series is too short for bandwidth selection.")
  
  ks_p_values <- vapply(h_grid, function(h) {
    h0 <- h^2
    v <- numeric(m - 1)
    for (t in 2:m) {
      v[t - 1] <- estimate_conditional_cdf(
        x_val = series[t], 
        y_cond = series[t - 1], 
        h = h, 
        series = series, 
        h0 = h0
      )
    }
    v <- pmax(1e-6, pmin(1 - 1e-6, v)) # Stabilize for ks.test
    suppressWarnings(ks.test(v, "punif")$p.value)
  }, numeric(1))
  
  h_grid[which.max(ks_p_values)]
}

#-----------------------------------------------------------------------
# 4. Main DCP Prediction Interval Function
#-----------------------------------------------------------------------

dcp_prediction_interval <- function(x, h_grid, alpha = 0.05) {
  
  n <- length(x)
  
  # --- Optimized Bandwidth Selection (done once) ---
  h_sel <- select_bandwidth_ks(x, h_grid)
  h0_sel <- h_sel^2
  
  # --- Candidate Grid for the next value ---
  # Heuristic grid centered around a simple forecast
  last_val <- x[n]
  forecast_center <- sin(last_val)
  ytrial <- seq(forecast_center - 4, forecast_center + 4, length.out = 100)
  
  yconfidence <- c()
  
  for (y in ytrial) {
    x_aug <- c(x, y)
    n_aug <- n + 1
    
    # Compute v-statistics for the augmented series using the pre-selected bandwidth
    v_stats <- vapply(2:n_aug, function(t) {
      estimate_conditional_cdf(
        x_val = x_aug[t], 
        y_cond = x_aug[t - 1], 
        h = h_sel, 
        series = x_aug, 
        h0 = h0_sel
      )
    }, numeric(1))
    
    v_stats <- pmax(1e-6, pmin(1 - 1e-6, v_stats))
    
    # The v-statistic for the candidate point y
    v_new <- v_stats[n_aug - 1] 
    
    # Calculate the one-sided p-value based on rank
    p_value <- mean(v_stats >= v_new)
    
    if (p_value > alpha) {
      yconfidence <- c(yconfidence, y)
    }
  }
  
  if (length(yconfidence) > 0) {
    interval <- range(yconfidence)
    list(lower = interval[1], upper = interval[2])
  } else {
    list(lower = NA, upper = NA)
  }
}

#-----------------------------------------------------------------------
# 5. Monte Carlo Simulation Driver (CORRECTED)
#-----------------------------------------------------------------------

run_simulation_parallel <- function(n, error_dist = c("Normal", "Laplace"),
                                    alpha = 0.05, num_sims = 50, burn_in = 500, B = 250) {
  error_dist <- match.arg(error_dist)
  
  # The foreach loop returns a data frame with all raw results
  results_df <- foreach(i = seq_len(num_sims), .combine = rbind) %dopar% {
    
    # --- Generate Data ---
    total_len <- n + burn_in
    eps <- if (error_dist == "Normal") {
      rnorm(total_len)
    } else {
      # Need to load VGAM on the worker nodes
      VGAM::rlaplace(total_len, location = 0, scale = 1 / sqrt(2))
    }
    
    x <- numeric(total_len)
    x[1] <- rnorm(1)
    for (t in 2:total_len) {
      x[t] <- sin(x[t - 1]) + eps[t]
    }
    
    x_train <- x[(burn_in + 1):total_len]
    
    # --- Generate the single true next value to check coverage against ---
    true_error <- if (error_dist == "Normal") {
      rnorm(B)
    } else {
      VGAM::rlaplace(B, location = 0, scale = 1 / sqrt(2))
    }
    x_true_next <- sin(x[total_len]) + true_error
    
    # Heuristic bandwidth grid
    h_grid <- seq(0.1, 1.5, length.out = 50) 
    
    # --- Get Prediction Interval ---
    interval <- try(dcp_prediction_interval(x_train, h_grid, alpha), silent = TRUE)
    
    # --- Record Results ---
    if (inherits(interval, "try-error") || is.na(interval$lower)) {
      data.frame(covered = NA, interval_length = NA)
    } else {
      covered <- mean(x_true_next >= interval$lower & x_true_next <= interval$upper)
      interval_length <- interval$upper - interval$lower
      data.frame(covered = covered, interval_length = interval_length)
    }
  }
  
  # --- Aggregate and Return Summary Statistics ---
  data.frame(
    n = n,
    error_dist = error_dist,
    nominal_coverage = 1 - alpha,
    CVR = mean(results_df$covered, na.rm = TRUE),
    LEN_mean = mean(results_df$interval_length, na.rm = TRUE),
    LEN_sd = sd(results_df$interval_length, na.rm = TRUE)
  )
}

#-----------------------------------------------------------------------
# 6. Main Execution Block
#-----------------------------------------------------------------------

# --- Setup Parallel Backend ---
num_cores <- detectCores() - 1 
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat(sprintf("Registered parallel backend with %d cores.\n", getDoParWorkers()))

# Export necessary functions and packages to the workers
clusterEvalQ(cl, {
  library(VGAM) # Make sure VGAM is available on each worker for Laplace dist
})
clusterExport(cl, c("kernel_lambda", "kernel_K", "estimate_conditional_cdf", 
                    "select_bandwidth_ks", "dcp_prediction_interval"))

# --- Define Simulation Parameters ---
param_grid <- expand.grid(
  n = c(50, 100, 200),
  error_dist = c("Normal", "Laplace"),
  alpha = c(0.05, 0.10),
  stringsAsFactors = FALSE
)


# --- Run All Simulations ---
all_results <- bind_rows(lapply(seq_len(nrow(param_grid)), function(i) {
  run_simulation_parallel(
    n = param_grid$n[i],
    error_dist = param_grid$error_dist[i],
    alpha = param_grid$alpha[i],
    num_sims = 500, # Number of Monte Carlo simulations
    burn_in = 500,
    B = 250
  )
}))

# --- Stop Cluster ---
stopCluster(cl)

# --- Print Final Results ---
cat("\nDCP Simulation Results:\n")
print(kable(all_results %>%
              mutate(
                CVR = round(CVR, 3), 
                LEN_mean = round(LEN_mean, 3),
                LEN_sd = round(LEN_sd, 3)
              ) %>%
              arrange(error_dist, nominal_coverage, n),
            format = "markdown"))