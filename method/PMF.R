source("utility.R")

compute_transformed_v_PMF <- function(x, p, h, h0) {
  n <- length(x)
  training_data <- make_train_xy(x, p)
  x_train <- training_data$x_train
  y_train <- training_data$y_train
  v <- numeric(n - p)
  for (t in (p + 1):n) {
    y_cond <- y_train[t - p, ]
    v[t - p] <- estimate_conditional_cdf_loo(
      u = x_train[t-p], y_cond = y_train[(t-p), ], p = p,
      series = x,
      h = h, h0 = h0, exclude_index = t
    )
  }
  pmax(1e-6, pmin(1 - 1e-6, v))
}

PMF <- function(x, p = 1, B = 250, M = NULL, alpha = 0.05) {
  n <- length(x)
  if (n <= p + 1) stop("Time series is too short.")
  if (is.null(M)) M <- max(p, floor(0.5 * n))
  
  
  # --- Optimized Bandwidth Selection (done once) ---
  h_cvls <- select_bandwidth_cvls(x, p) 
  h <- h_cvls$h_x
  h0 <- h_cvls$h_y
  
  original_training_data <- make_train_xy(x, p)
  x_train <- original_training_data$x_train
  y_train <- original_training_data$y_train
  v_P <- compute_transformed_v_PMF(x, p, h, h0) 
  
  y_test <- x[n:(n-p+1)]
  x_hat_n_1 <- mean(vapply(v_P, function(q) {
    inverse_conditional_cdf(q, y_test, x,p, h, h0)
  }, numeric(1)))
  roots <- numeric(B)
  for (b in seq_len(B)) {
    total_len <- (M+1) + (n+1)
    v_star <- sample(v_P, size = total_len, replace = TRUE)
    x_star <- numeric(total_len)
    x_star[1:p] <- draw_consecutive(x, p)
    for (t in (p + 1): (total_len - 1)) {
      y_cond_star <- x_star[(t - 1):(t - p)]
      x_star[t] <- inverse_conditional_cdf(v_star[t], y_cond_star, x,p, h, h0)
    }
    y_test <- x[n:(n-p+1)]
    # v_star_p <- sample(v, size = 1, replace = TRUE)
    x_star[total_len] <- inverse_conditional_cdf(v_star[total_len], y_test, x,p, h, h0)
    
    x_star_train <- x_star[(total_len - n):(total_len - 1)]
    # training_data_star <- make_train_xy(x_star_train, p)
    # x_star_train <- training_data_star$x_train
    # y_star_train <- training_data_star$y_train
    x_hat_star <- mean(vapply(v_star[(total_len - n + p):(total_len - 1)], function(q){
      inverse_conditional_cdf(q, y_test, x_star_train,p , h, h0)
    }, numeric(1)))
    roots[b] <- x_star[total_len] - x_hat_star
  }
  quantiles <- quantile(roots, c(alpha/ 2, 1 - alpha / 2), na.rm = TRUE)
  # quantiles_95 <- quantile(roots, c(0.05 / 2, 1 - 0.05 / 2), na.rm = TRUE)
  
  CI = list(
    lower = x_hat_n_1 + quantiles[1],
    upper = x_hat_n_1 + quantiles[2]
  )
  return(CI)
}