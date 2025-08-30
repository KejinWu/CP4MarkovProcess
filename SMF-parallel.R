## ======================================================================
## Conformal / SMF for Markov(p): one-file runnable script
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
# 2. Kernel Functions (保持不变)
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

make_train_xy <- function(x, p) {
  n <- length(x)
  if (n <= p) stop("Time series is too short for the specified order p.")
  x_train <- x[(p + 1):n]
  y_train <- embed(x, p + 1)[, -1, drop = FALSE]
  list(x_train = x_train, y_train = y_train)
}

compute_weights <- function(y_cond, y_train, h) {
  U <- sweep(y_train, 2, y_cond, "-")
  d <- sqrt(rowSums(U^2)) / h
  kernel_K(d)
}

estimate_conditional_cdf <- function(u, y_cond, x_train, y_train, h, h0) {
  weights <- compute_weights(y_cond, y_train, h)
  sum_weights <- mean(weights)
  if (sum_weights < 1e-12) return(mean(x_train <= u))
  numerator <- mean(weights * kernel_lambda((u - x_train) / h0))
  numerator / sum_weights
}

inverse_conditional_cdf <- function(q, y_cond, x_train, y_train, h, h0) {
  cdf_minus_q <- function(x_val) {
    estimate_conditional_cdf(x_val, y_cond, x_train, y_train, h, h0) - q
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
    y_cond <- y_train[t - p, ]
    v[t - p] <- estimate_conditional_cdf(
      u = x[t], y_cond = y_cond,
      x_train = x_train, y_train = y_train,
      h = h, h0 = h0
    )
  }
  pmax(1e-6, pmin(1 - 1e-6, v))
}

#-----------------------------------------------------------------------
# 4. Main SMF Bootstrap Function
#-----------------------------------------------------------------------

smf_bootstrap_interval <- function(x, p = 1, h, B = 250, alpha = 0.05, M = NULL) {
  n <- length(x)
  if (n <= p + 1) stop("Time series is too short.")
  if (is.null(M)) M <- max(p, floor(0.5 * n))
  h0 <- h^2
  original_training_data <- make_train_xy(x, p)
  x_train <- original_training_data$x_train
  y_train <- original_training_data$y_train
  v <- compute_transformed_v(x, p, h, h0)
  y_last <- y_train[n - p, ]
  x_tilde_n1 <- mean(vapply(v, function(q) {
    inverse_conditional_cdf(q, y_last, x_train, y_train, h, h0)
  }, numeric(1)))
  roots <- numeric(B)
  for (b in seq_len(B)) {
    v_star <- sample(v, size = (n - p + 1 + M), replace = TRUE)
    total_len <- n + 1 + M
    x_star <- numeric(total_len)
    x_star[1:p] <- draw_consecutive(x, p)
    for (t in (p + 1):total_len) {
      y_cond_star <- rev(x_star[(t - 1):(t - p)])
      q <- v_star[t - p]
      x_star[t] <- inverse_conditional_cdf(q, y_cond_star, x_train, y_train, h, h0)
    }
    x_star_train <- x_star[(M + 1):(M + n)]
    y_star_n1 <- x_star[M + n + 1]
    training_data_star <- make_train_xy(x_star_train, p)
    x_train_star <- training_data_star$x_train
    y_train_star <- training_data_star$y_train
    v_train_star <- compute_transformed_v(x_star_train, p, h, h0)
    y_last_star <- y_train_star[n - p, ]
    x_tilde_star_n1 <- mean(vapply(v_train_star, function(q) {
      inverse_conditional_cdf(q, y_last_star, x_train_star, y_train_star, h, h0)
    }, numeric(1)))
    roots[b] <- y_star_n1 - x_tilde_star_n1
  }
  quantiles <- quantile(roots, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  list(
    lower = x_tilde_n1 + quantiles[1],
    upper = x_tilde_n1 + quantiles[2],
    x_tilde = x_tilde_n1,
    roots = roots
  )
}

#-----------------------------------------------------------------------
# 5. Monte Carlo Simulation Driver (MODIFIED)
#-----------------------------------------------------------------------

run_simulation_parallel <- function(n, error_dist = c("Normal", "Laplace"),
                                    alpha = 0.05, num_sims = 50, B = 250, p = 1, M = NULL) {
  error_dist <- match.arg(error_dist)
  
  # The foreach loop now returns a data frame with all raw results
  results_df <- foreach(i = seq_len(num_sims), .combine = rbind, .packages = c("VGAM")) %dopar% {
    # Data generation remains the same
    eps <- if (error_dist == "Normal") {
      rnorm(n + 1, 0, 1)
    } else {
      VGAM::rlaplace(n + 1, location = 0, scale = 1 / sqrt(2))
    }
    
    x <- numeric(n + 1)
    x[1] <- rnorm(1)
    for (t in 1:n) {
      x[t + 1] <- sin(x[t]) + eps[t]
    }
    
    x_train <- x[1:n]
    x_true_next <- x[n + 1]
    h <- 1.06 * sd(x_train) * n^(-1/5)
    
    result <- try(smf_bootstrap_interval(x_train, p, h, B, alpha, M), silent = TRUE)
    
    if (inherits(result, "try-error")) {
      data.frame(covered = NA, interval_length = NA)
    } else {
      covered <- (x_true_next >= result$lower && x_true_next <= result$upper)
      interval_length <- result$upper - result$lower
      data.frame(covered = covered, interval_length = interval_length)
    }
  }
  
  # NEW: Calculate summary statistics from the collected data frame
  data.frame(
    n = n,
    error_dist = error_dist,
    nominal_coverage = 1 - alpha,
    CVR = mean(results_df$covered, na.rm = TRUE),
    LEN_mean = mean(results_df$interval_length, na.rm = TRUE),
    LEN_sd = sd(results_df$interval_length, na.rm = TRUE) # Added standard deviation
  )
}

#-----------------------------------------------------------------------
# 6. Main Execution Block (MODIFIED for new column names)
#-----------------------------------------------------------------------

num_cores <- detectCores() - 1 
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat(sprintf("Registered parallel backend with %d cores.\n", getDoParWorkers()))

clusterExport(cl, c("smf_bootstrap_interval", "draw_consecutive",
                    "make_train_xy", "compute_weights",
                    "estimate_conditional_cdf", "inverse_conditional_cdf",
                    "compute_transformed_v", "kernel_lambda", "kernel_K"),
              envir = environment())

param_grid <- expand.grid(
  n = c(50, 100), # Expanded for better comparison
  error_dist = c("Normal", "Laplace"),
  alpha = c(0.05),
  stringsAsFactors = FALSE
)

# The lapply call remains the same, as the new calculations are inside the function
all_results <- bind_rows(lapply(seq_len(nrow(param_grid)), function(i) {
  run_simulation_parallel(
    n = param_grid$n[i],
    error_dist = param_grid$error_dist[i],
    alpha = param_grid$alpha[i],
    num_sims = 500,
    B = 250,
    p = 1,
    M = floor(param_grid$n[i] / 2)
  )
}))

stopCluster(cl)

cat("\nSimulation Results:\n")
print(kable(all_results %>%
              mutate(
                CVR = round(CVR, 3), 
                LEN_mean = round(LEN_mean, 3),
                LEN_sd = round(LEN_sd, 3) # Round the new column
              ),
            format = "markdown"))