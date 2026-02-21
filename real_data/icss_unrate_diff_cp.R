install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

install_if_missing("ICSS")
suppressPackageStartupMessages(library(ICSS))

args <- commandArgs(trailingOnly = TRUE)
file_path <- if (length(args) >= 1) args[1] else "/Users/alexanderdai/Github/CP4MarkovProcess/real_data/UNRATE.csv"
demean_flag <- if (length(args) >= 2) as.logical(args[2]) else FALSE

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
n_raw <- length(x)
n_diff <- length(x_diff)

cp_idx <- ICSS::ICSS(x_diff, demean = demean_flag)
cp_idx <- as.integer(cp_idx)
cp_idx <- cp_idx[is.finite(cp_idx) & cp_idx > 0 & cp_idx <= n_diff]

if (length(cp_idx) == 0) {
  results <- data.frame(
    idx = integer(),
    date = as.Date(character()),
    diff_value = numeric(),
    pct_diff = numeric(),
    raw_idx = integer(),
    pct_raw = numeric()
  )
} else {
  results <- data.frame(
    idx = cp_idx,
    date = dates_diff[cp_idx],
    diff_value = x_diff[cp_idx],
    pct_diff = round(cp_idx / n_diff * 100, 2),
    raw_idx = cp_idx + 1L,
    pct_raw = round((cp_idx + 1L) / n_raw * 100, 2)
  )
}

cat("ICSS on diff(UNRATE)\n")
cat("n_raw =", n_raw, "\n")
cat("n_diff =", n_diff, "\n")
cat("demean =", demean_flag, "\n\n")

if (nrow(results) == 0) {
  cat("No variance change points detected by ICSS.\n")
} else {
  print(results)
}

out_csv <- sub("\\.csv$", "_diff_icss_results.csv", file_path)
write.csv(results, out_csv, row.names = FALSE)
cat("\nSaved results to:", out_csv, "\n")
