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
train_n  <- floor(0.7 * n)  

covered_90 <- logical(0)
covered_95 <- logical(0)

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
}

cvr_90 <- mean(covered_90)
cvr_95 <- mean(covered_95)

cat("Rolling 90% coverage:", round(cvr_90, 4), "\n")
cat("Rolling 95% coverage:", round(cvr_95, 4), "\n")
cat("Number of test points:", length(covered_90), "\n")


covered_90_smf <- logical(0)
covered_95_smf <- logical(0)

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
  interval <- smf_bootstrap_interval(x_train, p, h = h_x, h0 = h_y, B = B, M = M)
  
  lower90 <- interval$lower_90
  upper90 <- interval$upper_90
  lower95 <- interval$lower_95
  upper95 <- interval$upper_95
  
  y_next <- x[t + 1]
  
  covered_90_smf <- c(covered_90_smf, (y_next >= lower90 && y_next <= upper90))
  covered_95_smf <- c(covered_95_smf, (y_next >= lower95 && y_next <= upper95))
}

cvr_90_smf <- mean(covered_90_smf)
cvr_95_smf <- mean(covered_95_smf)

cat("Rolling 90% coverage:", round(cvr_90_smf, 4), "\n")
cat("Rolling 95% coverage:", round(cvr_95_smf, 4), "\n")


covered_90_pmf <- c()
covered_95_pmf <- c()
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
  interval <- pmf_bootstrap_interval(this_x, p, h = h_x, h0 = h_y, B = B, M = M)
  
  lower90 <- interval$lower_90
  upper90 <- interval$upper_90
  lower95 <- interval$lower_95
  upper95 <- interval$upper_95
  
  y_next <- x[t + 1]
  
  covered_90_pmf <- c(covered_90_pmf, (y_next >= lower90 && y_next <= upper90))
  covered_95_pmf <- c(covered_95_pmf, (y_next >= lower95 && y_next <= upper95))
}


cvr_90 <- mean(covered_90_pmf)
cvr_95 <- mean(covered_95_pmf)

cat("Rolling 90% coverage:", round(cvr_90, 4), "\n")
cat("Rolling 95% coverage:", round(cvr_95, 4), "\n")



covered_90_pmdcp <- c()
covered_95_pmdcp <- c()
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
}


cvr_90 <- mean(covered_90_pmdcp)
cvr_95 <- mean(covered_95_pmdcp)

cat("Rolling 90% coverage:", round(cvr_90, 4), "\n")
cat("Rolling 95% coverage:", round(cvr_95, 4), "\n")
