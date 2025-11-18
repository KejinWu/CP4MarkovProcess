library(readr)
library(latex2exp)
UNRATE <- read_csv("real_data/UNRATE.csv")
date <- UNRATE$observation_date
dat <- UNRATE$UNRATE
x <- diff(dat)
#
pdf("Unrate.pdf", width = 10, height = 4)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1)) 
plot(date, dat, type = "l", xlab = "Date", ylab = "Unemployment Rate")
plot(date[-1], x, type = "l", xlab = "Date", ylab = TeX("$\\Delta Y_t$"))
dev.off() 

n <- length(x)
train_n  <- floor(0.7 * n)  
#train_n <- 83


covered_90 <- logical(0)
covered_95 <- logical(0)
length_90 <- c()
length_95 <- c()

for (t in train_n:(n - 1)) {
  x_train <- x[(t-train_n+1):t]
  interval <- dcp_prediction_interval(x_train)

  lower90 <- interval$lower_90
  upper90 <- interval$upper_90
  lower95 <- interval$lower_95
  upper95 <- interval$upper_95

  y_next <- x[t + 1]
  covered_90 <- c(covered_90, (y_next >= lower90 && y_next <= upper90))
  covered_95 <- c(covered_95, (y_next >= lower95 && y_next <= upper95))
  length_90 <- c(length_90, upper90 - lower90)
  length_95 <- c(length_95, upper95 - lower95)
}

cvr_90 <- mean(covered_90)
cvr_95 <- mean(covered_95)
len_90 <- mean(length_90)
len_95 <- mean(length_95)
cat("Rolling 90% coverage:", round(cvr_90, 4), "\n")
cat("Rolling 95% coverage:", round(cvr_95, 4), "\n")
cat("Mean of 90% PI length:", round(len_90, 4), "\n")
cat("Mean of 95% PI length:", round(len_95, 4), "\n")
cat("Sd of 90% PI length:", round(sd(length_90), 4), "\n")
cat("Sd of 95% PI length:", round(sd(length_95), 4),"\n")


covered_90_smf <- c()
covered_95_smf <- c()
length_90_smf <- c()
length_95_smf <- c()
x_center <- c()
for (t in train_n:(n - 1)) {
  x_train <- x[(t-train_n+1):t]
  
  original_training_data <- make_train_xy(x_train, 1)
  x_train <- original_training_data$x_train
  y_train <- original_training_data$y_train
  data <- data.frame(x_train, y_train)
  h_cv_ls = npcdistbw(formula = x_train ~ y_train, data)
  h_x = h_cv_ls$ybw
  h_y = h_cv_ls$xbw
  
  B = 250
  p = 1
  M = NULL
  interval <- smf_bootstrap_interval(x[(t-train_n+1):t], p, h = h_x, h0 = h_y, B = B, M = M)
  
  lower90 <- interval$lower_90
  upper90 <- interval$upper_90
  lower95 <- interval$lower_95
  upper95 <- interval$upper_95
  x_center <- c(x_center, interval$x_tilde)
  print(interval$x_tilde)
  
  y_next <- x[t + 1]
  
  covered_90_smf <- c(covered_90_smf, (y_next >= lower90 && y_next <= upper90))
  covered_95_smf <- c(covered_95_smf, (y_next >= lower95 && y_next <= upper95))
  length_90_smf <- c(length_90_smf, upper90 - lower90)
  length_95_smf <- c(length_95_smf, upper95 - lower95)
}

cvr_90_smf <- mean(covered_90_smf)
cvr_95_smf <- mean(covered_95_smf)
len_90_smf <- mean(length_90_smf)
len_95_smf <- mean(length_95_smf)

cat("Rolling 90% coverage:", round(cvr_90_smf, 4), "\n")
cat("Rolling 95% coverage:", round(cvr_95_smf, 4), "\n")
cat("Mean of 90% PI length:", round(len_90_smf, 4), "\n")
cat("Mean of 95% PI length:", round(len_95_smf, 4), "\n")
cat("Sd of 90% PI length:", round(sd(length_90_smf), 4), "\n")
cat("Sd of 95% PI length:", round(sd(length_95_smf), 4),"\n")

covered_90_pmf <- c()
covered_95_pmf <- c()
length_90_pmf <- c()
length_95_pmf <- c()
x_center_pmf <- c()
for (t in train_n:(n - 1)) {
  x_train <- x[(t-train_n+1):t]

  original_training_data <- make_train_xy(x_train, 1)
  x_train <- original_training_data$x_train
  y_train <- original_training_data$y_train
  data <- data.frame(x_train, y_train)
  h_cv_ls = npcdistbw(formula = x_train ~ y_train, data)
  h_x = h_cv_ls$ybw
  h_y = h_cv_ls$xbw

  B = 250
  p = 1
  M = NULL
  interval <- pmf_bootstrap_interval(x[(t-train_n+1):t], p, h = h_x, h0 = h_y, B = B, M = M)

  lower90 <- interval$lower_90
  upper90 <- interval$upper_90
  lower95 <- interval$lower_95
  upper95 <- interval$upper_95
  x_center_pmf <- c(x_center_pmf, interval$x_tilde)
  print(interval$x_tilde)

  print(interval$x_tilde)

  y_next <- x[t + 1]

  covered_90_pmf <- c(covered_90_pmf, (y_next >= lower90 && y_next <= upper90))
  covered_95_pmf <- c(covered_95_pmf, (y_next >= lower95 && y_next <= upper95))
  length_90_pmf <- c(length_90_pmf, upper90 - lower90)
  length_95_pmf <- c(length_95_pmf, upper95 - lower95)
}


cvr_90_pmf <- mean(covered_90_pmf)
cvr_95_pmf <- mean(covered_95_pmf)
len_90_pmf <- mean(length_90_pmf)
len_95_pmf <- mean(length_95_pmf)
cat("Rolling 90% coverage:", round(cvr_90_pmf, 4), "\n")
cat("Rolling 95% coverage:", round(cvr_95_pmf, 4), "\n")
cat("Mean of 90% PI length:", round(len_90_pmf, 4), "\n")
cat("Mean of 95% PI length:", round(len_95_pmf, 4), "\n")
cat("Sd of 90% PI length:", round(sd(length_90_pmf), 4), "\n")
cat("Sd of 95% PI length:", round(sd(length_95_pmf), 4),"\n")



covered_90_pmdcp <- c()
covered_95_pmdcp <- c()
length_90_pmdcp <- c()
length_95_pmdcp <- c()
for (t in train_n:(n - 1)) {
  x_train <- x[(t-train_n+1):t]
  interval <- pmdcp_prediction_interval(x_train)
  
  lower90 <- interval$lower_90
  upper90 <- interval$upper_90
  lower95 <- interval$lower_95
  upper95 <- interval$upper_95
  
  y_next <- x[t + 1]
  
  covered_90_pmdcp <- c(covered_90_pmdcp, (y_next >= lower90 && y_next <= upper90))
  covered_95_pmdcp <- c(covered_95_pmdcp, (y_next >= lower95 && y_next <= upper95))
  length_90_pmdcp <- c(length_90_pmdcp, upper90 - lower90)
  length_95_pmdcp <- c(length_95_pmdcp, upper95 - lower95)
}


cvr_90 <- mean(covered_90_pmdcp)
cvr_95 <- mean(covered_95_pmdcp)
len_90_pmdcp <- mean(length_90_pmdcp)
len_95_pmdcp <- mean(length_95_pmdcp)
cat("Rolling 90% coverage:", round(cvr_90, 4), "\n")
cat("Rolling 95% coverage:", round(cvr_95, 4), "\n")
cat("Mean of 90% PI length:", round(len_90_pmdcp, 4), "\n")
cat("Mean of 95% PI length:", round(len_95_pmdcp, 4), "\n")
cat("Sd of 90% PI length:", round(sd(length_90_pmdcp), 4), "\n")
cat("Sd of 95% PI length:", round(sd(length_95_pmdcp), 4),"\n")

