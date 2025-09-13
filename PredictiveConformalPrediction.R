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
if (!require("np")) install.packages("np")

suppressPackageStartupMessages({
  library(VGAM)
  library(dplyr)
  library(knitr)
  library(doParallel)
  library(foreach)
  library(np)
})

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

# Make the train data to fit the local constant estimator for conditional CDF
make_train_xy <- function(x, p) {
  n <- length(x)
  if (n <= p) stop("Time series is too short for the specified order p.")
  x_train <- x[(p + 1):n]
  y_train <- embed(x, p + 1)[, -1, drop = FALSE] # Read the y_train from the first row to the second row and so on
  list(x_train = x_train, y_train = y_train)
}


# Conditional CDF Estimator (Leave-One-Out Version)
estimate_conditional_cdf <- function(x_val, y_cond, h, series, h0 = h^2, exclude_index) {
  m <- length(series)
  p <- 1
  if (m < p + 1) stop("Series is too short for p=1.")
  
  indices_to_keep <- 1:(m-p)
  
  loo_pair_index <- exclude_index - p
  if (loo_pair_index > 0 && loo_pair_index <= m-p) {
    indices_to_keep <- indices_to_keep[-loo_pair_index]
  }
  
  # original training data
  x_train <- series[2:m]
  y_train <- series[1:(m-1)]
  
  # leave-one-out training data
  x_train <- x_train[indices_to_keep]
  y_train <- y_train[indices_to_keep]
  
  weights <- kernel_K((y_cond - y_train) / h)
  sum_weights <- sum(weights)
  
  if (!is.finite(sum_weights) || sum_weights < 1e-12) {
    return(mean(x_train <= x_val))
  }
  
  numerator <- sum(weights * kernel_lambda((x_val - x_train) / h0))
  numerator / sum_weights
}

# Bandwidth selection via Kolmogorov-Smirnov test p-value maximization
select_bandwidth_ks <- function(series, p,  h_grid) {
  m <- length(series)
  if (m < 3) stop("Series is too short for bandwidth selection.")
  
  ks_p_values <- vapply(h_grid, function(h) {
    h0 <- h^2
    v <- numeric(m - p)
    for (t in (p+1):m) {
      v[t - p] <- estimate_conditional_cdf(
        x_val = series[t], 
        y_cond = series[t - p], 
        h = h, 
        series = series, 
        h0 = h0,
        exclude_index = t
      )
    }
    v <- pmax(1e-6, pmin(1 - 1e-6, v)) # Stabilize for ks.test
    suppressWarnings(ks.test(v, "punif")$p.value)
  }, numeric(1))
  
  h_grid[which.max(ks_p_values)]
}


# ANOTHER CHOICE to pick the badnwidth by npcdistbw function


select_bandwidth_cvls <- function(series, p) {
  
  original_training_data <- make_train_xy(series, p)
  x_train <- original_training_data$x_train
  y_train <- original_training_data$y_train
  data <- data.frame(x_train, y_train)
  h_cv_ls = npcdistbw(formula = x_train ~ y_train, data)
  h_x = h_cv_ls$ybw # since we treat x as the Y_{p+1} and y as Y_{p}
  h_y = h_cv_ls$xbw  
  
  return(list(h_x = h_x, h_y = h_y))
}

#-----------------------------------------------------------------------
# 4. Main DCP Prediction Interval Function
#-----------------------------------------------------------------------

dcp_prediction_interval <- function(x, p = 1) {
  
  n <- length(x)
  
  # --- Optimized Bandwidth Selection (done once) ---
  h_cvls <- select_bandwidth_cvls(x, p) #############################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!########################
  h_sel <- h_cvls$h_x
  h0_sel <- h_cvls$h_y
  
  
  # --- Candidate Grid for the next value ---
  # Heuristic grid centered around a simple forecast
  last_val <- x[n]
  ytrial <- seq(-max(abs(x)), max(abs(x)), length.out = 100)
  
  yconfidence_90 <- c()
  yconfidence_95 <- c()
  
  for (y in ytrial) {
    x_aug <- c(x, y)
    n_aug <- n + 1
    
    # Compute v-statistics for the augmented series using the pre-selected bandwidth
    v_stats <- vapply((p+1):n_aug, function(t) {
      estimate_conditional_cdf(
        x_val = x_aug[t], 
        y_cond = x_aug[(t - p):(t-1)], #############################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!######################## 
        # For p = 1, it is fine. If we want to consider order p > 1, it is not appropriate, we should use x_aug[(t - p):(t-1)]
        h = h_sel, 
        series = x_aug, 
        h0 = h0_sel,
        exclude_index = t
      )
    }, numeric(1))
    
    #v_stats <- pmax(1e-6, pmin(1 - 1e-6, v_stats))
    v_stats <- abs(v_stats - 1/2)
    # The v-statistic for the candidate point y
    v_new <- v_stats[n_aug - 1] 
    
    # Calculate the one-sided p-value based on rank
    p_value <- mean(v_stats >= v_new)
    
    if (p_value > 0.05) {
      yconfidence_95 <- c(yconfidence_95, y)
    }
    if (p_value > 0.1) {
      yconfidence_90 <- c(yconfidence_90, y)
    }
  }
  
  if (length(yconfidence_95) > 0) {
    interval_95 <- range(yconfidence_95)
    list_95 = list(lower_95 = interval_95[1], upper_95 = interval_95[2])
  } else {
    list_95 = list(lower_95 = NA, upper_95 = NA)
  }
  
  if (length(yconfidence_90) > 0) {
    interval_90 <- range(yconfidence_90)
    list_90 = list(lower_90 = interval_90[1], upper_90 = interval_90[2])
  } else {
    list_90 = list(lower_90 = NA, upper_90 = NA)
  }
  
  return(c(list_90,list_95))
}

#-----------------------------------------------------------------------
# 5. Monte Carlo Simulation Driver (CORRECTED)
#-----------------------------------------------------------------------

run_simulation_parallel_DCP <- function(n, error_dist = c("Normal", "Laplace"),
                                        num_sims = 50, S = 5000, burn_in = 500, B = 250,  p = 1) {
  #error_dist <- match.arg(error_dist)
  
  # The foreach loop returns a data frame with all raw results
  results_df <- foreach(i = seq_len(num_sims), .combine = rbind,.packages = c("VGAM")) %dopar% {
    
    # --- Generate Data ---
    set.seed(i)
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
      rnorm(S)
    } else {
      VGAM::rlaplace(S, location = 0, scale = 1 / sqrt(2))
    }
    x_true_next <- sin(x[total_len]) + true_error
    
    # --- Get Prediction Interval ---
    interval <- try(dcp_prediction_interval(x_train, p), silent = TRUE)
    
    # --- Record Results ---
    if (inherits(interval, "try-error") || is.na(interval$lower_90)) {
      data.frame(covered = NA, interval_length = NA)
    } else {
      covered_90 <- mean(x_true_next >= interval$lower_90 & x_true_next <= interval$upper_90)
      interval_length_90 <- interval$upper_90 - interval$lower_90
      covered_95 <- mean(x_true_next >= interval$lower_95 & x_true_next <= interval$upper_95)
      interval_length_95 <- interval$upper_95 - interval$lower_95
      data.frame(covered_90 = covered_90, interval_length_90 = interval_length_90,
                 covered_95 = covered_95, interval_length_95 = interval_length_95)
    }
  }
  
  # --- Aggregate and Return Summary Statistics ---
  data.frame(
    n = n,
    error_dist = error_dist,
    CVR_90 = mean(results_df$covered_90, na.rm = TRUE),
    CVR_95 = mean(results_df$covered_95, na.rm = TRUE),
    LEN_mean_90 = mean(results_df$interval_length_90, na.rm = TRUE),
    LEN_mean_95 = mean(results_df$interval_length_95, na.rm = TRUE),
    LEN_sd_90 = sd(results_df$interval_length_90, na.rm = TRUE),
    LEN_sd_95 = sd(results_df$interval_length_95, na.rm = TRUE) # Added standard deviation
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
                    "select_bandwidth_ks", "dcp_prediction_interval", "npcdistbw", "make_train_xy","select_bandwidth_cvls"))

# --- Define Simulation Parameters ---
param_grid <- expand.grid(
  n = c(50, 100, 200),
  error_dist = c("Normal", "Laplace"),
  stringsAsFactors = FALSE
)


# --- Run All Simulations ---
all_results_DCP <- bind_rows(lapply(seq_len(nrow(param_grid)), function(i) {
  print(paste("Run setting n =",param_grid$n[i],"error_dist = ",param_grid$error_dist[i], sep = " "))
  run_simulation_parallel_DCP(
    n = param_grid$n[i],
    error_dist = param_grid$error_dist[i],
    num_sims = 5, # Number of Monte Carlo simulations
    burn_in = 500,
    S = 5000,
    B = 250,
    p = 1
  )
}))

# --- Stop Cluster ---
stopCluster(cl)

# # --- Print Final Results ---
# cat("\nDCP Simulation Results:\n")
# print(kable(all_results %>%
#               mutate(
#                 CVR = round(CVR, 3), 
#                 LEN_mean = round(LEN_mean, 3),
#                 LEN_sd = round(LEN_sd, 3)
#               ) %>%
#               arrange(error_dist, nominal_coverage, n),
#             format = "markdown"))