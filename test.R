set.seed(1)
n = 50
burn_in = 100
error_dist = "Normal"
p = 1
total_len <- n + burn_in
eps <- if (error_dist == "Normal") {
  rnorm(total_len)
} else {
  # Need to load VGAM on the worker nodes
  VGAM::rlaplace(total_len, location = 0, scale = 1 / sqrt(2))
}

x <- numeric(total_len)
x[1] <- rnorm(1)
for (t in 2:total_len) {
  x[t] <- sin(x[t - 1]) + eps[t]
}

orinianl_x_train <- x[(burn_in + 1):total_len]

1.06 * sd(x_train) * n^(-1/5)

original_training_data <- make_train_xy(orinianl_x_train, p)
x_train <- original_training_data$x_train
y_train <- original_training_data$y_train
y_train_2 <- y_train
data <- data.frame(x_train, y_train, y_train_2)
test = npcdistbw(formula = x_train ~ y_train + y_train_2, data)
test$xbw
test$ybw
