install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

install_if_missing("changepoint")
install_if_missing("strucchange")

suppressPackageStartupMessages(library(changepoint))
suppressPackageStartupMessages(library(strucchange))

args <- commandArgs(trailingOnly = TRUE)
file_path <- if (length(args) >= 1) args[1] else "/Users/alexanderdai/Github/CP4MarkovProcess/real_data/UNRATE.csv"
minseglen <- if (length(args) >= 2) as.integer(args[2]) else 8

if (!file.exists(file_path)) {
  stop("File not found: ", file_path)
}

d <- read.csv(file_path, stringsAsFactors = FALSE)
if (!all(c("observation_date", "UNRATE") %in% names(d))) {
  stop("CSV must contain columns: observation_date, UNRATE")
}

d$observation_date <- as.Date(d$observation_date)
x <- as.numeric(d$UNRATE)
keep <- is.finite(x)
if (!all(keep)) {
  d <- d[keep, ]
  x <- x[keep]
}

x_diff <- diff(x)
dates_diff <- d$observation_date[-1]

fit_mean <- cpt.mean(x_diff, method = "PELT", penalty = "MBIC", minseglen = minseglen)
fit_var <- cpt.var(x_diff, method = "PELT", penalty = "MBIC", minseglen = minseglen)
fit_meanvar <- cpt.meanvar(x_diff, method = "PELT", penalty = "MBIC", minseglen = minseglen)

map_cp <- function(cp_idx, label) {
  cp_idx <- as.integer(cp_idx)
  cp_idx <- cp_idx[cp_idx > 0 & cp_idx <= length(x_diff)]
  if (length(cp_idx) == 0) {
    return(data.frame(method = character(), idx = integer(), date = as.Date(character()), diff_value = numeric()))
  }
  data.frame(
    method = label,
    idx = cp_idx,
    date = dates_diff[cp_idx],
    diff_value = x_diff[cp_idx]
  )
}

res_mean <- map_cp(cpts(fit_mean), "mean")
res_var <- map_cp(cpts(fit_var), "variance")
res_meanvar <- map_cp(cpts(fit_meanvar), "mean+variance")

bp <- breakpoints(x_diff ~ 1, h = 0.1)
k_hat <- which.min(BIC(bp)) - 1
bp_fit <- breakpoints(bp, breaks = k_hat)
res_bp <- map_cp(bp_fit$breakpoints[is.finite(bp_fit$breakpoints)], "strucchange_BIC")

results <- rbind(res_mean, res_var, res_meanvar, res_bp)

cat("n_raw =", length(x), "\n")
cat("n_diff =", length(x_diff), "\n")
cat("minseglen =", minseglen, "\n\n")

if (nrow(results) == 0) {
  cat("No change points detected by configured methods.\n")
} else {
  print(results)
}

out_csv <- sub("\\.csv$", "_diff_cp_results.csv", file_path)
write.csv(results, out_csv, row.names = FALSE)
cat("\nSaved results to:", out_csv, "\n")
