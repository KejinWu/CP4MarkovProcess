# --- smooth CDF: clamp to [-2,2] then Phi ---
lambda <- function(x) {
  z <- pmax(-2, pmin(2, x))
  pnorm(z)
}

K <- function(x) dnorm(x)

# 条件CDF估计器  F̂(x | y)  for Markov(p=1)
# series: 训练序列（数值向量），长度 m >= 2
# x_val:  目标值 (标量)    ; y_cond: 条件值 (标量)
D <- function(x_val, y_cond, h, series, h0 = h^2, eps = 1e-12) {
  m <- length(series)
  if (m < 2) stop("series too short for p=1")
  
  # 训练对：x_i = series[2:m] ; y_i = series[1:(m-1)]
  x_train <- series[2:m]
  y_train <- series[1:(m-1)]
  
  w <- K((y_cond - y_train) / h)
  sw <- sum(w)
  if (!is.finite(sw) || sw < eps) return(mean(x_train <= x_val))
  
  num <- sum(w * lambda((x_val - x_train) / h0))
  num / sw
}

# 计算 v_t = F̂(X_t | X_{t-1})，并做 KS 选带宽（p=1）
select <- function(series, h_grid) {
  m <- length(series)
  if (m < 3) stop("series too short")
  
  ks_p <- numeric(length(h_grid))
  for (j in seq_along(h_grid)) {
    h  <- h_grid[j]
    h0 <- h^2
    v  <- numeric(m - 1)
    # t = 2..m  =>  v_{t-1} = F̂( X_t | X_{t-1} )
    for (t in 2:m) {
      v[t - 1] <- D(x_val = series[t], y_cond = series[t - 1], h = h,
                    series = series, h0 = h0)
    }
    v <- pmax(1e-6, pmin(1 - 1e-6, v))  # 稳定
    ks_p[j] <- suppressWarnings(ks.test(v, "punif")$p.value)
  }
  h_grid[ which.max(ks_p) ]
}

set.seed(1)
B <- 500
result  <- numeric(B)  # 覆盖率
result2 <- numeric(B)  # 区间长度

n_window <- 50 # Window size used to do training
burn_out <- 500 # Burn out numeber to remove the effects of initial value effects
h_grid <- seq(0.5 * n_window^(-1/4), n_window^(-1/5), length.out = 400)
alpha <- 0.05

# WE ONLY CONSIDER ONE STEP AHEAD PREDICTION SO FAR
for (step in 1:B) {
  # 生成数据：X_{t+1} = sin(X_t) + e_{t+1}
  total_length <- n_window + burn_out
  e <- rnorm(total_length)
  x <- numeric(total_length); x[1] <- 0
  for (i in 2:total_length) x[i] <- sin(x[i - 1]) + e[i] 
  
  data <- x[(burn_out + 1) : (burn_out + n_window)]  # 训练段，长度 n
  # 测试点用于估计经验覆盖（只是评估用）
  ytest <- rnorm(1000) + sin(data[burn_out + n_window]) # Use 1000 pseduo future values to test the coverage rate
  
  # 网格候选
  ytrial <- seq(-1.25 * max(abs(data)), 1.25 * max(abs(data)), length.out = 100)
  
  # —— DCP: 对每个候选 y 构造“增广数据”并计算 v-rank —— #
  yconfidence <- c()
  for (y in ytrial) {
    data_aug <- c(data, y)          # 增广到长度 n+1
    n_aug    <- length(data_aug)
    
    # 为增广数据选择带宽（也可复用原 h，这里随增广更新）
    h_sel <- select(data_aug, h_grid)
    h0_sel <- h_sel^2
    
    # 计算增广数据上的 v 统计量（t = 2..n+1）
    v <- numeric(n_aug - 1)
    for (t in 2:n_aug) {
      v[t - 1] <- D(x_val = data_aug[t], y_cond = data_aug[t - 1],
                    h = h_sel, series = data_aug, h0 = h0_sel)
    }
    v <- pmax(1e-6, pmin(1 - 1e-6, v))
    
    # DCP 的 p-value（rank-based；可选两端或一端）
    v_new <- v[n_aug - 1]                 # 最后一对 (y | data[n])
    pval  <- mean(v >= v_new)             # 单侧：大于等于
    if (pval > alpha) yconfidence <- c(yconfidence, y)
  }
  
  # 记录覆盖与长度
  if (length(yconfidence) > 0) {
    result[step]  <- mean(ytest <= max(yconfidence) & ytest >= min(yconfidence))
    result2[step] <- max(yconfidence) - min(yconfidence)
  } else {
    result[step]  <- NA_real_
    result2[step] <- NA_real_
  }
}

cat("Coverage (mean):", mean(result,  na.rm = TRUE), "\n")
cat("Length   (mean):", mean(result2, na.rm = TRUE), "\n")
cat("Length     (sd):",   sd(result2, na.rm = TRUE), "\n")






