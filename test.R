set.seed(1)
n = 50
burn_in = 100
error_dist = "Normal"
p = 1
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

orinianl_x_train <- x[(burn_in + 1):total_len]

1.06 * sd(x_train) * n^(-1/5)

original_training_data <- make_train_xy(orinianl_x_train, p)
x_train <- original_training_data$x_train
y_train <- original_training_data$y_train
y_train_2 <- y_train
data <- data.frame(x_train, y_train, y_train_2)
test = npcdistbw(formula = x_train ~ y_train + y_train_2, data)
test$xbw
test$ybw





run_simulation_parallel_MF <- function(n, error_dist = c("Normal"),
                                       alpha = 0.05, num_sims = 50, burn_in = 500, S = 5000, B = 250, p = 1, M = 20) {
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
    #print(h_y)
    # NOTE we define x as y; SO ybw is corresponding with h_x i.e., h
    # Select the bandwidth by cv.ls
    
    result <- try(smf_bootstrap_interval(all_x, p, h = h_x, h0 = h_y, B = B, alpha, M = M), silent = TRUE)
    
    if (inherits(result, "try-error")) {
      data.frame(covered = NA, interval_length = NA)
    } else {
      covered <- mean(x_true_next >= result$lower & x_true_next <= result$upper)
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


run_simulation_parallel_DCP <- function(n, error_dist = c("Normal"),
                                        num_sims = 50, S = 5000, burn_in = 500, B = 250,  p = 1) {
  error_dist <- match.arg(error_dist)
  
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
    
    # Heuristic bandwidth grid
    h_grid <- seq(0.1, 1.5, length.out = 50) 
    
    # --- Get Prediction Interval ---
    interval <- try(dcp_prediction_interval(x_train, p, h_grid), silent = TRUE)
    
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

