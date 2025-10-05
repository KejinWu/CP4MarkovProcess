## ======================================================================
## Smoothed Model-free Prediction (SMF) for Markov(p=1)
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
# 2. Kernel Functions
#-----------------------------------------------------------------------

kernel_lambda <- function(x) {
  clamped_x <- pmax(-2, pmin(2, x))
  pnorm(clamped_x)
}

kernel_K <- function(u) {
  dnorm(u)
}

#-----------------------------------------------------------------------
# 3. Helper Functions
#-----------------------------------------------------------------------

draw_consecutive <- function(x, p) {
  n <- length(x)
  if (p > n) stop("p cannot be greater than the time series length.")
  start <- sample(1:(n - p + 1), 1)
  x[start:(start + p - 1)]
}


# x_train is Y_{p+1}
# y_train is Y_{1,2,3..p} and so on.
make_train_xy <- function(x, p) {
  n <- length(x)
  if (n <= p) stop("Time series is too short for the specified order p.")
  x_train <- x[(p + 1):n]
  y_train <- embed(x, p + 1)[, -1, drop = FALSE] # Read the y_train from the first row to the second row and so on
  list(x_train = x_train, y_train = y_train)
}


compute_weights <- function(y_cond, y_train, h) {
  kernel_K((y_cond - y_train)/h)
}

# y_cond is the p continuing predictor of the time series
# u is value of the p+1 index of the time series

estimate_conditional_cdf_MF <- function(u, y_cond, x_train, y_train, h, h0) {
  weights <- compute_weights(y_cond, y_train, h)
  sum_weights <- mean(weights)
  if (sum_weights < 1e-12) return(mean(x_train <= u))
  numerator <- mean(weights * kernel_lambda((u - x_train) / h0))
  numerator / sum_weights
}

inverse_conditional_cdf <- function(q, y_cond, x_train, y_train, h, h0) {
  cdf_minus_q <- function(x_val) {
    estimate_conditional_cdf_MF(x_val, y_cond, x_train, y_train, h, h0) - q # minus q so that later we can solve to get tthe quantile
  }
  rng <- range(x_train)
  pad <- 3 * sd(x_train)
  a <- rng[1] - pad
  b <- rng[2] + pad
  out <- try(uniroot(cdf_minus_q, c(a, b))$root, silent = TRUE)
  if (inherits(out, "try-error")) {
    as.numeric(quantile(x_train, probs = q, type = 1))
  } else {
    out
  }
}

compute_transformed_v <- function(x, p, h, h0) {
  n <- length(x)
  training_data <- make_train_xy(x, p)
  x_train <- training_data$x_train
  y_train <- training_data$y_train
  v <- numeric(n - p)
  for (t in (p + 1):n) {
    v[t - p] <- estimate_conditional_cdf_MF(
      u = x_train[t - p], y_cond = y_train[(t - p), ],
      x_train = x_train, y_train = y_train,
      h = h, h0 = h0
    )
  }
  pmax(1e-6, pmin(1 - 1e-6, v))
}

#-----------------------------------------------------------------------
# 4. Main SMF Bootstrap Function
#-----------------------------------------------------------------------

smf_bootstrap_interval_optimal <- function(x, h, h0, p = 1, B = 250, M = NULL) {
  n <- length(x)
  if (n <= p + 1) stop("Time series is too short.")
  if (is.null(M)) M <- max(p, floor(0.5 * n))
  original_training_data <- make_train_xy(x, p)
  x_train <- original_training_data$x_train
  y_train <- original_training_data$y_train
  v <- compute_transformed_v(x, p, h, h0)
  y_test <- x[n:(n-p+1)]
  beta <- seq(0.1, 0.9, 9)
  x_hat_n_1 <- rep(0, length(beta))
  for (i in 1:length(beta)){
    v_q <- quantile(v, beta[i])
    x_hat_n_1[i] <- vapply(v_q, function(q){
      inverse_conditional_cdf(q, y_test, x_train, y_train, h, h0)
    }, numeric(1))

  # x_hat_n_1 <- mean(vapply(v, function(q) {
  #   inverse_conditional_cdf(q, y_test, x_train, y_train, h, h0)
  # }, numeric(1)))
  # v_q <- quantile(v, beta)
  # x_hat_n_1 <- vapply(v_q, function(q){
  #     inverse_conditional_cdf(q, y_test, x_train, y_train, h, h0)
  #   }, numeric(1))
###
# length optimization depends on the q(1-alpha/2 + gamma ) - q(alpha/2 + gamma) only
# does not depend on the "center" x_hat_n_1
# unless we use similar operations (find quantile) in Bootstrap
    roots <- numeric(B)
    for (b in seq_len(B)) {
      total_len <- (M+1) + (n+1)
      v_star <- sample(v, size = total_len, replace = TRUE)
      x_star <- numeric(total_len)
      x_star[1:p] <- draw_consecutive(x, p)
      for (t in (p + 1): (total_len - 1)) {
        y_cond_star <- x_star[(t - 1):(t - p)]
        x_star[t] <- inverse_conditional_cdf(v_star[t], y_cond_star, x_train, y_train, h, h0)
      }
      y_test <- x[n:(n-p+1)]
      x_star[total_len] <- inverse_conditional_cdf(v_star[total_len], y_test, x_train, y_train, h, h0)

      x_star_train <- x_star[(total_len - n):(total_len - 1)]
      training_data_star <- make_train_xy(x_star_train, p)
      x_star_train <- training_data_star$x_train
      y_star_train <- training_data_star$y_train
    # x_hat_star <- mean(vapply(v_star[(total_len - n + p):(total_len - 1)], function(q){
    #   inverse_conditional_cdf(q, y_test, x_star_train, y_star_train, h, h0)
    # }, numeric(1)))
      v_star_q <- quantile(v_star[(total_len - n + p):(total_len - 1)], beta[i])
      x_hat_star <- vapply(v_star_q, function(q){
        inverse_conditional_cdf(q, y_test, x_star_train, y_star_train, h, h0)
      }, numeric(1))
      roots[b] <- x_star[total_len] - x_hat_star
    }

    grid <- seq(-0.1/2, 0.1/2, 20)
    grid2 <- seq(-0.05/2, 0.05/2, 20)
    len_90 <- matrix(0, nrow = length(beta), ncol = length(grid))
    len_95 <- matrix(0, nrow = length(beta), ncol = length(grid))
    for (j in 1:length(grid)){
      quantiles_90 <- quantile(roots, c(0.1 / 2 + grid[j], 1 - 0.1 / 2 + grid[j]), na.rm = TRUE)
      quantiles_95 <- quantile(roots, c(0.05 / 2 + grid2[j], 1 - 0.05 / 2 + grid2[j]), na.rm = TRUE)
      len_90[i, j] <- quantiles_90[2] - quantiles_90[1]
      len_95[i, j] <- quantiles_95[2] - quantiles_90[1]
    }
    idx_90 <- which(len_90 == min(len_90, na.rm = TRUE), arr.ind = TRUE)
    row_90 <- idx_90[1, 1]
    col_90 <- idx_90[1, 2]

    idx_95 <- which(len_95 == min(len_95, na.rm = TRUE), arr.ind = TRUE)
    row_95 <- idx_95[1, 1]
    col_95 <- idx_95[1, 2]

    quantiles_90 <- quantile(roots, c(0.1 / 2 + grid[col_90], 1 - 0.1 / 2 + grid[which.min(len_90)]), na.rm = TRUE)
    quantiles_95 <- quantile(roots, c(0.05 / 2 + grid2[col_95], 1 - 0.05 / 2 + grid2[which.min(len_95)]), na.rm = TRUE)
    x_hat_n_1_90 <- x_hat_n_1[row_90]
    x_hat_n_1_95 <- x_hat_n_1[row_95]
  }

  list(
    lower_90 = x_hat_n_1_90 + quantiles_90[1],
    upper_90 = x_hat_n_1_90 + quantiles_90[2],
    lower_95 = x_hat_n_1_95 + quantiles_95[1],
    upper_95 = x_hat_n_1_95 + quantiles_95[2]
    # x_tilde = x_hat_n_1,
    # roots = roots
  )
}

#-----------------------------------------------------------------------
# 5. Monte Carlo Simulation Driver
#-----------------------------------------------------------------------

#' @param num_sims the number of simulations we want to replicate
#' @param B the number of bootstrap we take to do predictions
#' @param S the number of future values we generate to evaluate the coverage
#' @param M the burn in for MF prediction approach
#' @param p the order of the time series model


run_simulation_parallel_MF <- function(n, error_dist = c("Normal", "Laplace"),
                                       num_sims = 50, burn_in = 500, S = 5000, B = 250, p = 1, M = NULL) {
  error_dist <- match.arg(error_dist)
  
  # The foreach loop now returns a data frame with all raw results
  results_df <- foreach(i = seq_len(num_sims), .combine = rbind, .packages = c("VGAM")) %dopar% {
    
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
    
    all_x <- x[(burn_in + 1):total_len]
    
    # --- Generate the single true next value to check coverage against ---
    true_error <- if (error_dist == "Normal") {
      rnorm(S)
    } else {
      VGAM::rlaplace(S, location = 0, scale = 1 / sqrt(2))
    }
    x_true_next <- sin(x[total_len]) + true_error
    
    # h <- 1.06 * sd(x_train) * n^(-1/5)
    
    # NOTE can only work for markov(1) process for this moment
    # We use default setting to find the optimal bandwidth for x_train-- h_cv_ls$ybw and y_train--h_cv_ls$xbw
    # Since we treat x_train as Y 
    original_training_data <- make_train_xy(all_x, p)
    x_train <- original_training_data$x_train
    y_train <- original_training_data$y_train
    data <- data.frame(x_train, y_train)
    h_cv_ls = npcdistbw(formula = x_train ~ y_train, data)
    h_x = h_cv_ls$ybw
    h_y = h_cv_ls$xbw   
    # NOTE we define x as y; SO ybw is corresponding with h_x i.e., h
    # Select the bandwidth by cv.ls
    
    result <- try(smf_bootstrap_interval_optimal(all_x, p, h = h_x, h0 = h_y, B = B, M = M), silent = TRUE)
    
    if (inherits(result, "try-error")) {
      data.frame(covered = NA, interval_length = NA)
    } else {
      covered_90 <- mean(x_true_next >= result$lower_90 & x_true_next <= result$upper_90)
      interval_length_90 <- result$upper_90 - result$lower_90
      covered_95 <- mean(x_true_next >= result$lower_95 & x_true_next <= result$upper_95)
      interval_length_95 <- result$upper_95 - result$lower_95
      data.frame(covered_90 = covered_90, interval_length_90 = interval_length_90,
                 covered_95 = covered_95, interval_length_95 = interval_length_95)
    }
  }
  
  # NEW: Calculate summary statistics from the collected data frame
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

num_cores <- detectCores() - 1 
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat(sprintf("Registered parallel backend with %d cores.\n", getDoParWorkers()))

clusterExport(cl, c("smf_bootstrap_interval_optimal", "draw_consecutive",
                    "make_train_xy", "compute_weights",
                    "estimate_conditional_cdf_MF", "inverse_conditional_cdf",
                    "compute_transformed_v", "kernel_lambda", "kernel_K", "npcdistbw"),
              envir = environment())

param_grid <- expand.grid(
  n = c(50, 100, 200), # Expanded for better comparison
  error_dist = c("Normal", "Laplace"),
  alpha = c(0.05, 0.1),
  stringsAsFactors = FALSE
)


# The lapply call remains the same, as the new calculations are inside the function
all_results_MF <- bind_rows(lapply(seq_len(nrow(param_grid)), function(i) {
  print(paste("Run setting n =",param_grid$n[i],"error_dist = ",param_grid$error_dist[i], sep = " "))
  run_simulation_parallel_MF(
    n = param_grid$n[i],
    error_dist = param_grid$error_dist[i],
    num_sims = 5,
    burn_in = 500,
    B = 10,
    S = 10,
    p = 1,
    M = floor(param_grid$n[i] / 2)
  )
}
))

stopCluster(cl)

# cat("\nSimulation Results:\n")
# print(kable(all_results %>%
#               mutate(
#                 CVR = round(CVR, 3), 
#                 LEN_mean = round(LEN_mean, 3),
#                 LEN_sd = round(LEN_sd, 3) # Round the new column
#               ),
#             format = "markdown"))
# 
# 


