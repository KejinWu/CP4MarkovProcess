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

naive_conformal_prediction <- function(x, p = 1){
  n <- length(x)
  x_train <- x[1:(n-p)]
  y_train <- x[(p+1):n]
  regData <- data.frame(x_train, y_train)
  fit_lm <- lm(y_train ~ x_train, data=regData)
  e_hat <- abs(fit_lm$residuals)
  
  last_val <- x[n]
  print(x[n])
  mu_hat <- predict(fit_lm, newdata = data.frame(x_train = last_val))
  alpha <- 0.05
  LB_95 <- mu_hat - quantile(e_hat, 1- alpha)
  UB_95 <- mu_hat + quantile(e_hat, 1- alpha)
  
  alpha <- 0.1
  LB_90 <- mu_hat - quantile(e_hat, 1- alpha)
  UB_90 <- mu_hat + quantile(e_hat, 1- alpha)
  
  list_95 <- list(lower_95 = LB_95, upper_95 = UB_95)
  list_90 <- list(lower_90 = LB_90, upper_90 = UB_90)
  return(c(list_90,list_95))
}

run_simulation_parallel_naive_CP <- function(n, error_dist = c("Normal","Laplace"),
                                             num_sims = 50, S = 5000, burn_in = 500, B = 250, p = 1){
  results_df <- foreach(i = seq_len(num_sims), .combine = rbind,
                        .packages = c("VGAM","np"),
                        .export = c("naive_conformal_prediction")) %dopar% {
    set.seed(i)
    total_len <- n + burn_in
    eps <- if (error_dist == "Normal") {
      rnorm(total_len)
    } else {
      # Need to load VGAM on the worker nodes
      VGAM::rlaplace(total_len, location = 0, scale = 1 / sqrt(2))
    }
    
    x <- numeric(total_len)
    x[1] <- 0
    for (t in 2:total_len) {
      x[t] <- sin(x[t - 1]) + eps[t]
    }
    
    x_train <- x[(burn_in + 1):total_len]  
    true_error <- if (error_dist == "Normal") {
      rnorm(S)
    } else {
      VGAM::rlaplace(S, location = 0, scale = 1 / sqrt(2))
    }
    x_true_next <- sin(x[total_len]) + true_error
    interval <-naive_conformal_prediction(x_train)
    covered_90 <- mean(x_true_next >= interval$lower_90 & x_true_next <= interval$upper_90)
    interval_length_90 <- interval$upper_90 - interval$lower_90
    covered_95 <- mean(x_true_next >= interval$lower_95 & x_true_next <= interval$upper_95)
    interval_length_95 <- interval$upper_95 - interval$lower_95
    data.frame(covered_90 = covered_90, interval_length_90 = interval_length_90,
               covered_95 = covered_95, interval_length_95 = interval_length_95)
    
  }
  data.frame(
    n = n,
    error_dist = error_dist,
    CVR_90 = mean(results_df$covered_90, na.rm = TRUE),
    CVR_95 = mean(results_df$covered_95, na.rm = TRUE),
    LEN_mean_90 = mean(results_df$interval_length_90, na.rm = TRUE),
    LEN_mean_95 = mean(results_df$interval_length_95, na.rm = TRUE),
    LEN_sd_90 = sd(results_df$interval_length_90, na.rm = TRUE),
    LEN_sd_95 = sd(results_df$interval_length_95, na.rm = TRUE)
  )
}

# --- Setup Parallel Backend ---
num_cores <- detectCores() - 1 
cl <- makeCluster(num_cores)
registerDoParallel(cl)
cat(sprintf("Registered parallel backend with %d cores.\n", getDoParWorkers()))

# Export necessary functions and packages to the workers
clusterEvalQ(cl, {
  library(VGAM) # Make sure VGAM is available on each worker for Laplace dist
})
clusterExport(cl, c("naive_conformal_prediction"))

# --- Define Simulation Parameters ---
param_grid <- expand.grid(
  n = c(50, 100, 200),
  error_dist = c("Normal", "Laplace"),
  stringsAsFactors = FALSE
)


# --- Run All Simulations ---
all_results_naive_CP <- bind_rows(lapply(seq_len(nrow(param_grid)), function(i) {
  print(paste("Run setting n =",param_grid$n[i],"error_dist = ",param_grid$error_dist[i], sep = " "))
  run_simulation_parallel_naive_CP(
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
cat("\nnaive CP Simulation Results:\n")
print(all_results_naive_CP)
# print(kable(all_results_naive_CP %>%
#               mutate(
#                 CVR = round(CVR, 3),
#                 LEN_mean = round(LEN_mean, 3),
#                 LEN_sd = round(LEN_sd, 3)
#               ) %>%
#               arrange(error_dist, nominal_coverage, n),
#             format = "markdown"))
