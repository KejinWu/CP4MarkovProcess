#!/usr/bin/env Rscript

library(readr)
library(parallel)

parse_args <- function(args) {
  out <- list(
    input = "real_data/UNRATE.csv",
    train_ratio = 0.5,
    alphas = c(0.05, 0.1),
    seed = 2026,
    cores = NA_integer_,
    output = "real_data_new_results.csv",
    p_output = "real_data_new_p_estimation.csv"
  )

  for (a in args) {
    if (startsWith(a, "--input=")) {
      out$input <- sub("^--input=", "", a)
    } else if (startsWith(a, "--train-ratio=")) {
      out$train_ratio <- as.numeric(sub("^--train-ratio=", "", a))
    } else if (startsWith(a, "--seed=")) {
      out$seed <- as.integer(sub("^--seed=", "", a))
    } else if (startsWith(a, "--cores=")) {
      out$cores <- as.integer(sub("^--cores=", "", a))
    } else if (startsWith(a, "--output=")) {
      out$output <- sub("^--output=", "", a)
    } else if (startsWith(a, "--p-output=")) {
      out$p_output <- sub("^--p-output=", "", a)
    } else if (startsWith(a, "--alphas=")) {
      alpha_txt <- sub("^--alphas=", "", a)
      out$alphas <- as.numeric(strsplit(alpha_txt, ",", fixed = TRUE)[[1]])
    } else if (a == "--help") {
      cat("Usage:\n")
      cat("  Rscript real_data_new.R [options]\n\n")
      cat("Options:\n")
      cat("  --input=PATH            Input CSV path (default: real_data/UNRATE.csv)\n")
      cat("  --train-ratio=NUM       Train ratio in (0, 1), default: 0.5\n")
      cat("  --alphas=A,B            Alpha list, default: 0.05,0.1\n")
      cat("  --seed=INT              Random seed (default: 2026)\n")
      cat("  --cores=INT             Parallel cores (default: min(tasks, detectCores))\n")
      cat("  --output=PATH           Results table CSV (default: real_data_new_results.csv)\n")
      cat("  --p-output=PATH         P-estimation table CSV (default: real_data_new_p_estimation.csv)\n")
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", a)
    }
  }

  if (!is.finite(out$train_ratio) || out$train_ratio <= 0 || out$train_ratio >= 1) {
    stop("--train-ratio must be in (0, 1).")
  }
  if (length(out$alphas) == 0 || any(!is.finite(out$alphas)) || any(out$alphas <= 0 | out$alphas >= 1)) {
    stop("--alphas must be valid probabilities in (0, 1), e.g. --alphas=0.05,0.1")
  }
  out
}

get_script_dir <- function() {
  args_full <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args_full, value = TRUE)
  if (length(file_arg) == 0L) {
    return(getwd())
  }
  dirname(normalizePath(sub("^--file=", "", file_arg[1])))
}

run_parallel <- function(items, fun, cores = NA_integer_) {
  n_jobs <- length(items)
  if (n_jobs == 0L) {
    return(list())
  }
  if (is.na(cores)) {
    n_cores <- max(1L, min(detectCores(logical = TRUE), n_jobs))
  } else {
    n_cores <- max(1L, min(as.integer(cores), n_jobs))
  }
  if (.Platform$OS.type == "windows" || n_cores == 1L) {
    return(lapply(items, fun))
  }
  mclapply(items, fun, mc.cores = n_cores)
}

list_to_df <- function(x) {
  if (length(x) == 0L) {
    return(data.frame())
  }
  do.call(rbind.data.frame, c(lapply(x, as.data.frame), list(stringsAsFactors = FALSE)))
}

rolling_test <- function(series, train_n, method_name, p, alpha) {
  covered <- logical(0)
  interval_length <- numeric(0)
  n_series <- length(series)

  for (t in train_n:(n_series - 1)) {
    x_train <- series[(t - train_n + 1):t]

    interval <- switch(
      method_name,
      MDCP = MDCP(x = x_train, p = p, alpha = alpha),
      SMF = SMF(x = x_train, p = p, B = 250, M = NULL, alpha = alpha),
      PMF = PMF(x = x_train, p = p, B = 250, M = NULL, alpha = alpha),
      PMDCP = PMDCP(x = x_train, p = p, alpha = alpha),
      stop("Unknown method name: ", method_name)
    )

    lower <- interval$lower
    upper <- interval$upper
    if (!is.finite(lower) || !is.finite(upper)) {
      next
    }

    y_next <- series[t + 1]
    covered <- c(covered, (y_next >= lower && y_next <= upper))
    interval_length <- c(interval_length, upper - lower)
  }

  list(
    rolling_coverage = mean(covered),
    mean_pi_length = mean(interval_length),
    sd_pi_length = sd(interval_length),
    valid_steps = length(interval_length)
  )
}

script_dir <- get_script_dir()
setwd(script_dir)

args <- parse_args(commandArgs(trailingOnly = TRUE))
set.seed(args$seed)

source("method/order_est.R")

UNRATE <- read_csv(args$input, show_col_types = FALSE)
x <- diff(UNRATE$UNRATE)
n <- length(x)
train_n <- floor(args$train_ratio * n)

method_map <- data.frame(
  method_name = c("MDCP", "SMF", "PMF", "PMDCP"),
  order_method = c("conformal", "MF", "PMF", "conformal_predict"),
  stringsAsFactors = FALSE
)

config_grid <- merge(
  method_map,
  data.frame(alpha = args$alphas),
  by = NULL
)

estimate_p_worker <- function(idx) {
  row <- config_grid[idx, ]
  p_hat <- tryCatch(
    order_est(x[1:train_n], method = row$order_method, alpha = row$alpha),
    error = function(e) NA_integer_
  )
  list(
    method = row$method_name,
    order_est_method = row$order_method,
    alpha = row$alpha,
    p_hat = as.integer(p_hat)
  )
}

estimated_p_list <- run_parallel(seq_len(nrow(config_grid)), estimate_p_worker, cores = args$cores)
estimated_p_list <- Filter(function(z) !is.na(z$p_hat), estimated_p_list)

if (length(estimated_p_list) == 0L) {
  stop("No valid p estimate found. Please check order_est or input data.")
}

p_table <- list_to_df(estimated_p_list)
p_table <- p_table[order(p_table$method, p_table$alpha), ]

cat("========================================\n")
cat("P estimation table\n")
cat("train_n:", train_n, "(", round(args$train_ratio, 4), "* n )\n")
print(p_table, row.names = FALSE)

test_worker <- function(item) {
  stats <- rolling_test(
    series = x,
    train_n = train_n,
    method_name = item$method,
    p = item$p_hat,
    alpha = item$alpha
  )
  list(
    method = item$method,
    alpha = item$alpha,
    p_hat = item$p_hat,
    rolling_coverage = stats$rolling_coverage,
    mean_pi_length = stats$mean_pi_length,
    sd_pi_length = stats$sd_pi_length,
    valid_steps = stats$valid_steps
  )
}

results_list <- run_parallel(estimated_p_list, test_worker, cores = args$cores)
results_table <- list_to_df(results_list)
results_table <- results_table[order(results_table$method, results_table$alpha), ]

cat("========================================\n")
cat("Rolling summary table\n")
print(results_table, row.names = FALSE)

write.csv(p_table, args$p_output, row.names = FALSE)
write.csv(results_table, args$output, row.names = FALSE)

cat("========================================\n")
cat("Saved p table to:", args$p_output, "\n")
cat("Saved result table to:", args$output, "\n")
