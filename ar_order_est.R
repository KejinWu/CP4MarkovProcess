x <- as.numeric(x)
x <- x[!is.na(x)]
n <- length(x)

pmax <- 12   # n=280 建议 8~12 足够做 sanity check

BICs <- rep(NA_real_, pmax)
AICs <- rep(NA_real_, pmax)

for (p in 1:pmax) {
  fit <- arima(x, order = c(p, 0, 0), include.mean = TRUE, method = "ML")
  k <- length(fit$coef) + 1   # 参数数：phi(1..p) + mean + sigma^2
  BICs[p] <- -2 * fit$loglik + k * log(n)
  AICs[p] <- -2 * fit$loglik + 2 * k
}

p_bic <- which.min(BICs)
p_aic <- which.min(AICs)

p_bic
BICs

plot(1:pmax, BICs, type = "b", xlab = "p", ylab = "BIC")
abline(v = p_bic, lty = 2)
