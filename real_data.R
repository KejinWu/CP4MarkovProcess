library(readr)
library(latex2exp)
UNRATE <- read_csv("real_data/UNRATE.csv")
date <- UNRATE$observation_date
dat <- UNRATE$UNRATE
x <- diff(dat)

pdf("Unrate.pdf", width = 10, height = 4)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1)) 
plot(date, dat, type = "l", xlab = "Date", ylab = "Unemployment Rate")
plot(date[-1], x, type = "l", xlab = "Date", ylab = TeX("$\\Delta Y_t$"))
dev.off() 

n <- length(x)
# train_n  <- floor(0.7 * n)  
train_n <- floor(0.5 * n)
#train_n <- 83
alpha = 0.05

source("method/order_est.R")
order_est(x[1:train_n], method = "conformal", alpha = alpha)

order_est(x[1:train_n], method = "conformal_predict", alpha = alpha)

order_est(x[1:train_n], method = "PMF", alpha = alpha)

############################################################
##### MDCP ##################################################
############################################################
source("method/MDCP.R")
covered <- logical(0)
length <- c()

mdcp_upper<- c()
mdcp_lower <- c()
p <- 4 # p <- 1
for (t in train_n:(n - 1)) {
  x_train <- x[(t-train_n+1):t]
  interval <- MDCP(x_train, p, alpha = alpha)
  
  lower <- interval$lower
  upper <- interval$upper
 
  mdcp_upper <- c(mdcp_upper, upper)
  mdcp_lower <- c(mdcp_lower, lower)

  y_next <- x[t + 1]
  covered <- c(covered, (y_next >= lower && y_next <= upper))
  length <- c(length, upper - lower)
}

cvr <- mean(covered)
len <- mean(length)
cat("Rolling", (1-alpha)*100, "% coverage:", round(cvr, 4), "\n")
cat("Mean of",(1-alpha)*100,"% PI length:", round(len, 4), "\n")
cat("Sd of",(1-alpha)*100,"% PI length:", round(sd(length), 4),"\n")

############################################################
##### SMF ##################################################
############################################################
source("method/SMF.R")
p = 4
covered_smf <- c()
length_smf <- c()

smf_upper<- c()
smf_lower <- c()

for (t in train_n:(n - 1)) {
  x_train <- x[(t-train_n+1):t]
  B = 250
  M = NULL
  interval <- SMF(x = x_train, p = p, B = B, M = M, alpha = alpha)
  
  lower <- interval$lower
  upper <- interval$upper
  
  smf_upper <- c(smf_upper, upper)
  smf_lower <- c(smf_lower, lower)

  y_next <- x[t + 1]
  
  covered_smf <- c(covered_smf, (y_next >= lower && y_next <= upper))
  length_smf <- c(length_smf, upper - lower)
}

cvr_smf <- mean(covered_smf)
len_smf <- mean(length_smf)

cat("Rolling", (1-alpha)*100, "% coverage:", round(cvr_smf, 4), "\n")
cat("Mean of",(1-alpha)*100,"% PI length:", round(len_smf, 4), "\n")
cat("Sd of",(1-alpha)*100,"% PI length:", round(sd(length_smf), 4),"\n")

############################################################
##### PMF ##################################################
############################################################

covered_pmf <- c()
length_pmf <- c()
pmf_upper<- c()
pmf_lower <- c()
p <- 1
for (t in train_n:(n - 1)) {
  x_train <- x[(t-train_n+1):t]
  
  B = 250
  M = NULL
  interval <- PMF(x = x_train, p = p, B = B, M = M, alpha = alpha)

  lower <- interval$lower
  upper <- interval$upper
  pmf_upper <- c(pmf_upper, upper)
  pmf_lower <- c(pmf_lower, lower)

  y_next <- x[t + 1]
  
  covered_pmf <- c(covered_pmf, (y_next >= lower && y_next <= upper))
  length_pmf <- c(length_pmf, upper-lower)
}


cvr_pmf <- mean(covered_pmf)
len_pmf <- mean(length_pmf)
cat("Rolling", (1-alpha)*100, "% coverage:", round(cvr_pmf, 4), "\n")
cat("Mean of",(1-alpha)*100,"% PI length:", round(len_pmf, 4), "\n")
cat("Sd of",(1-alpha)*100,"% PI length:", round(sd(length_pmf), 4),"\n")


############################################################
##### PMDCP ##################################################
############################################################

source("method/PMDCP.R")
alpha = 0.05
p = 7
covered_pmdcp <- c()
length_pmdcp <- c()
pmdcp_upper<- c()
pmdcp_lower <- c()

for (t in train_n:(n - 1)) {
  x_train <- x[(t-train_n+1):t]
  interval <- PMDCP(x = x_train, p = p, alpha = alpha)
  
  lower <- interval$lower
  upper <- interval$upper
  pmdcp_upper <- c(pmdcp_upper, upper)
  pmdcp_lower <- c(pmdcp_lower, lower)

  y_next <- x[t + 1]
  
  covered_pmdcp <- c(covered_pmdcp, (y_next >= lower && y_next <= upper))
  length_pmdcp <- c(length_pmdcp, upper - lower)
  length_pmdcp <- c(length_pmdcp, upper - lower)
}


cvr_pmdcp <- mean(covered_pmdcp)
len_pmdcp <- mean(length_pmdcp)
cat("Rolling", (1-alpha)*100, "% coverage:", round(cvr_pmdcp, 4), "\n")
cat("Mean of",(1-alpha)*100,"% PI length:", round(len_pmdcp, 4), "\n")
cat("Sd of",(1-alpha)*100,"% PI length:", round(sd(length_pmdcp), 4),"\n")





plot(ts(x[(train_n+1):(n)]), ylim = c(-2,2))
lines(pmdcp_upper_90, col = "red", lty = 2)
lines(pmdcp_lower_90, col = "red", lty = 2)
lines(smf_upper_90, col = "blue", lty = 2)
lines(smf_lower_90, col = "blue", lty = 2)
lines(pmf_upper_90, col = "green", lty = 2)
lines(pmf_lower_90, col = "green", lty = 2)
lines(mdcp_upper_90, col = "yellow", lty = 2)
lines(mdcp_lower_90, col = "yellow", lty = 2)



legend("topleft",
       legend = c("True", "PMDCP-90", "SMF-90", "PMF-90", "MDCP-90"),
       col = c("black", "red", "blue","green","yellow"),
       lty = c(1, 1,1,1,1),
       lwd = c(1, 1,1,1,1),
       cex = 0.5)


pdf("rolling_plot.pdf", width = 8, height = 5)   # 8×5 英寸画布
plot(ts(x[(train_n+1):(n)]), ylim = c(-2,2),
     xlab = "Rolling prediction step",
     ylab = "Value")
lines(pmdcp_upper_90, col = "red", lty = 2)
lines(pmdcp_lower_90, col = "red", lty = 2)
lines(smf_upper_90, col = "blue", lty = 2)
lines(smf_lower_90, col = "blue", lty = 2)
lines(pmf_upper_90, col = "green", lty = 2)
lines(pmf_lower_90, col = "green", lty = 2)
lines(mdcp_upper_90, col = "yellow", lty = 2)
lines(mdcp_lower_90, col = "yellow", lty = 2)



legend("topleft",
       legend = c("True", "PMDCP-90", "SMF-90", "PMF-90", "MDCP-90"),
       col = c("black", "red", "blue","green","yellow"),
       lty = c(1, 1,1,1,1),
       lwd = c(1, 1,1,1,1),
       cex = 0.5)
dev.off()