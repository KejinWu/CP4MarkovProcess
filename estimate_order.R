# install.packages("knnmi")
library(knnmi)

x <- as.numeric(x)
x <- x[!is.na(x)]

make_lags <- function(x, L){
  n <- length(x)
  M <- matrix(NA_real_, n - L, L + 1)
  for (j in 0:L) M[, j + 1] <- x[(L + 1 - j):(n - j)]
  M
}

# CMI: I(X_t ; X_{t-p-1} | X_{t-1..t-p})
cmi_for_p <- function(x, p, k = 3){
  L <- p + 1
  M <- make_lags(x, L)                    # rows = samples
  Xt <- M[, 1]                            # length N
  Z  <- M[, 2:(p + 1), drop = FALSE]      # N x p  (samples x vars)
  W  <- M[, p + 2]                        # length N
  
  # knnmi 要求 Z 是 (vars x samples) = p x N
  cond_mutual_inf(Xt, W, t(Z), k = k)
}

perm_test_cmi <- function(x, p, k = 3, B = 800, seed = 1){
  set.seed(seed)
  L <- p + 1
  M <- make_lags(x, L)
  Xt <- M[, 1]
  Z  <- M[, 2:(p + 1), drop = FALSE]
  W  <- M[, p + 2]
  
  obs <- cond_mutual_inf(Xt, W, t(Z), k = k)
  
  null_vals <- numeric(B)
  for (b in 1:B){
    null_vals[b] <- cond_mutual_inf(Xt, sample(W), t(Z), k = k)
  }
  pval <- (1 + sum(null_vals >= obs)) / (B + 1)
  list(stat = obs, pval = pval)
}

# ---------- main ----------
pmax <- 6          # n=280 建议 6（最多 8，但会更没power）
k <- 3            # 小样本常用 3 或 5
B <- 1000           # 置换次数 >= 500 更稳
alpha <- 0.1
alpha_adj <- alpha / pmax  # Bonferroni

res <- data.frame(p = 1:pmax, CMI = NA_real_, pval = NA_real_)
for (p in 1:pmax){
  out <- perm_test_cmi(x, p, k = k, B = B, seed = 100 + p)
  res$CMI[p] <- out$stat
  res$pval[p] <- out$pval
}
print(res)

# 选阶：最小 p 使得“连续两次”不拒绝（pval > alpha_adj）
good <- res$pval > alpha_adj
p_hat <- NA_integer_
for (p in 1:(pmax - 1)){
  if (good[p] && good[p + 1]) { p_hat <- p; break }
}
p_hat
