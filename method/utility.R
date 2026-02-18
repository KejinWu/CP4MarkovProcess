if (!require("VGAM")) install.packages("VGAM")
if (!require("dplyr")) install.packages("dplyr")
if (!require("np")) install.packages("np")

suppressPackageStartupMessages({
  library(VGAM)
  library(dplyr)
  library(np)
})



#-----------------------------------------------------------------------
# 1. Kernel Functions
#-----------------------------------------------------------------------

kernel_lambda <- function(x) {
  clamped_x <- pmax(-2, pmin(2, x))
  pnorm(clamped_x)
}

kernel_K <- function(u) {
  dnorm(u)
}

#-----------------------------------------------------------------------
# 2. Basic Functions
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
  y_train <- embed(x, p + 1)[, -1, drop = FALSE] # Read the y_train from the first row to the second row and so on
  list(x_train = x_train, y_train = y_train)
}


compute_weights <- function(y_cond, y_train, h) {
  kernel_K((y_cond - y_train)/h)
}

#-----------------------------------------------------------------------
# 3. CDF Functions
#-----------------------------------------------------------------------

estimate_conditional_cdf <- function(x_val, y_cond, series, p, h, h0) {
  n = length(series)
  x_train <- series[(p+1):n]
  y_train <- embed(series,p+1)[, -1, drop = FALSE]
  weights <- compute_weights(y_cond, y_train, h)
  sum_weights <- mean(weights)
  if (sum_weights < 1e-12) return(mean(x_train <= x_val))
  numerator <- mean(weights * kernel_lambda((x_val - x_train) / h0))
  numerator / sum_weights
}

inverse_conditional_cdf <- function(q, y_cond, series, p, h, h0) {
  cdf_minus_q <- function(x_val) {
    estimate_conditional_cdf(x_val, y_cond, series, p, h, h0) - q # minus q so that later we can solve to get tthe quantile
  }
  n = length(series)
  x_train <- series[(p+1):n]
  # y_train <- embed(series,p+1)[, -1, drop = FALSE]
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


estimate_conditional_cdf_loo <- function(x_val, y_cond, p, h, series, h0 = h^2, exclude_index) {
  m <- length(series)
  if (m < p + 1) stop(paste("Series is too short for p =", p))
  indices_to_keep <- 1:(m-p)
  
  loo_pair_index <- exclude_index - p
  if (loo_pair_index > 0 && loo_pair_index <= m-p) {
    indices_to_keep <- indices_to_keep[-loo_pair_index]
  }
  
  # original training data
  x_train <- series[(p+1):m]
  y_train <- embed(series, p + 1)[, -1]
  y_train <- as.matrix(y_train)
  
  # leave-one-out training data
  x_train <- x_train[indices_to_keep]

  y_train <- y_train[indices_to_keep, , drop = FALSE]
  
  
  weights <- kernel_K((y_cond - y_train) / h)
  sum_weights <- sum(weights)
  
  if (!is.finite(sum_weights) || sum_weights < 1e-12) {
    return(mean(x_train <= x_val))
  }
  
  numerator <- sum(weights * kernel_lambda((x_val - x_train) / h0))
  numerator / sum_weights
}



select_bandwidth_cvls <- function(series, p) {
  original_training_data <- make_train_xy(series, p)
  x_train <- original_training_data$x_train
  y_train <- original_training_data$y_train

  h_cv_ls = npcdistbw(xdat = y_train, ydat = x_train)
  h_x = h_cv_ls$ybw
  h_y = h_cv_ls$xbw

  return(list(h_x = h_x, h_y = h_y))
}