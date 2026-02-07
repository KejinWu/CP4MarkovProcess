source("utility.R")

PMDCP <- function(x, p = 1) {
  
  n <- length(x)
  
  # --- Optimized Bandwidth Selection (done once) ---
  h_cvls <- select_bandwidth_cvls(x, p) 
  h_sel <- h_cvls$h_x
  h0_sel <- h_cvls$h_y
  
  
  # --- Candidate Grid for the next value ---
  # Heuristic grid centered around a simple forecast
  last_val <- x[n]
  ytrial <- seq(-max(abs(x)), max(abs(x)), length.out = 500)
  
  yconfidence_90 <- c()
  yconfidence_95 <- c()
  
  for (y in ytrial) {
    x_aug <- c(x, y)
    n_aug <- n + 1
    
    # Compute v-statistics for the augmented series using the pre-selected bandwidth
    v_stats <- vapply((p+1):n_aug, function(t) {
      estimate_conditional_cdf_loo(
        x_val = x_aug[t], 
        y_cond = x_aug[(t - p):(t-1)], p = p,
        h = h_sel, 
        series = x_aug, 
        h0 = h0_sel,
        exclude_index = t
      )
    }, numeric(1))
    
    #v_stats <- pmax(1e-6, pmin(1 - 1e-6, v_stats))
    v_stats <- abs(v_stats - 1/2)
    # The v-statistic for the candidate point y
    v_new <- v_stats[length(v_stats)] 
    
    # Calculate the one-sided p-value based on rank
    p_value <- mean(v_stats >= v_new)
    
    if (p_value > alpha) {
      yconf <- c(yconf, y)
    }
  }
  
  if (length(yconf) > 0) {
    range_CI <- range(yconf)
    CI = list(lower = range_CI[1], upper = range_CI[2])
  } else {
    CI = list(lower = NA, upper = NA)
  }
  
  return(CI)
}
