source("utility.R")

compute_transformed_v <- function(x, p, h, h0) {
  n <- length(x)
  training_data <- make_train_xy(x, p)
  x_train <- training_data$x_train
  y_train <- training_data$y_train
  v <- numeric(n - p)
  for (t in (p + 1):n) {
    v[t - p] <- estimate_conditional_cdf(
      x_val = x_train[t - p], y_cond = y_train[(t - p), ], p = p,
      series = x,
      h = h, h0 = h0
    )
  }
  pmax(1e-6, pmin(1 - 1e-6, v))
}



SMF <- function(x, p = 1, B = 250, M = NULL, alpha = 0.05) {
  n <- length(x)
  if (n <= p + 1) stop("Time series is too short.")
  if (is.null(M)) M <- max(p, floor(0.5 * n))
  # original_training_data <- make_train_xy(x, p)
  # x_train <- original_training_data$x_train
  # y_train <- original_training_data$y_train
  
  # --- Optimized Bandwidth Selection (done once) ---
  h_cvls <- select_bandwidth_cvls(x, p)
  h <- h_cvls$h_x
  h0 <- h_cvls$h_y
  
  v <- compute_transformed_v(x, p, h, h0)
  y_test <- x[n:(n-p+1)]
  x_hat_n_1 <- mean(vapply(v, function(q) {
    inverse_conditional_cdf(q, y_test, x, p, h, h0)
  }, numeric(1)))
  roots <- numeric(B)
  for (b in seq_len(B)) {
    total_len <- (M+1) + (n+1)
    v_star <- sample(v, size = total_len, replace = TRUE)
    x_star <- numeric(total_len)
    x_star[1:p] <- draw_consecutive(x, p)
    for (t in (p + 1): (total_len - 1)) {
      y_cond_star <- x_star[(t - 1):(t - p)]
      x_star[t] <- inverse_conditional_cdf(v_star[t], y_cond_star, x,p, h, h0)
    }
    y_test <- x[n:(n-p+1)]
    x_star[total_len] <- inverse_conditional_cdf(v_star[total_len], y_test, x, p,h, h0)
    
    x_star_train <- x_star[(total_len - n):(total_len - 1)]
    # training_data_star <- make_train_xy(x_star_train, p)
    # x_star_train <- training_data_star$x_train
    # y_star_train <- training_data_star$y_train
    x_hat_star <- mean(vapply(v_star[(total_len - n + p):(total_len - 1)], function(q){
      inverse_conditional_cdf(q, y_test, x_star_train, p, h, h0)
    }, numeric(1)))
    roots[b] <- x_star[total_len] - x_hat_star
  }
  quantiles <- quantile(roots, c(alpha/ 2, 1 - alpha / 2), na.rm = TRUE)
  CI = list(
    lower = x_hat_n_1 + quantiles[1],
    upper = x_hat_n_1 + quantiles[2]
  )
  return(CI)
}
