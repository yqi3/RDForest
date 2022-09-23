######################################
# Table 1: (7-9) Random Forests
# Simulation: Univariate score using data from Lee (2008)

# We implement random forests using grf on the following platform
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS release 6.5 (Final)
######################################
rm(list=ls())
set.seed(1234)

library(grf)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Lee_Helper_Fns.R")

c_r0 = 0.52
c_r1 = 0.84
c_r2 = -3
c_r3 = 7.99
c_r4 = -9.01
c_r5 = 3.56
c_l0 = 0.48
c_l1 = 1.27
c_l2 = 7.18
c_l3 = 20.21
c_l4 = 21.54
c_l5 = 7.33
c_rz = 0
c_lz = 0
sigma_y = 0.1295

truth <- c_r0 - c_l0
focal <- matrix(0)

#### MC Simulation ####
t <- 1000 # number of Monte Carlo iterations

# Number of observations in each Monte Carlo iteration
# Comment/uncomment the three lines below to reproduce columns (7) - (9) of results in the paper
n <- 1000  # (7)
# n <- 5000  # (8)
# n <- 10000  # (9)

B <- 5000 # number of trees
scale <- 0.4  # for tuning sample fraction
tune.num.trees <- 50 # number of trees in each forest
tune.num.reps <- 100 # number of forests to be built
tune.num.draws <- 1000 # number of random draws, evaluated using the Dice Kriging model
moderate.num.trees <- tune.num.trees * 4
tune.parameters <- c("alpha", "imbalance.penalty", "honesty.fraction", "honesty.prune.leaves")  # in accordance with Helper_Fns
num.params <- length(tune.parameters)
lower_bound_default <- (1+ log(1-0.05) / log(0.05))^(-1)  # d=1, alpha=0.05
result_rf <- c()
coverage_rf <- c()
BLB1 <- c()
BLB0 <- c()

for (i in 1:t){
  print(i)
  # generate simulation data
  X <- 2*rbeta(n,2,4)-1
  epsilon <- rnorm(n, 0, sigma_y)
  Y <- ifelse(X>=0, c_r0+c_r1*X+c_r2*X^2+c_r3*X^3+c_r4*X^4+c_r5*X^5+epsilon,
              c_l0+c_l1*X+c_l2*X^2+c_l3*X^3+c_l4*X^4+c_l5*X^5+epsilon)

  data <- as.data.frame(cbind("Y" = Y, "X"=X))

  ### RF ###
  trt <- data[data$X>=0,]
  ctrl <- data[data$X<0,]
  feature_trt <- as.data.frame(trt[,2])
  feature_ctrl <- as.data.frame(ctrl[,2])

  # -----------------------------------------------------------------------------------
  # Modified tuning procedure based on the source code of grf; see
  # https://github.com/grf-labs/grf/blob/master/r-package/grf/R/tune_forest.R
  # -----------------------------------------------------------------------------------

  ## Parameter tuning for treated
  fit.draws <- get_fit_draws(tune.num.reps)
  draw.parameters <- get_draw_parameters(fit.draws)
  draw.parameters <- cbind(draw.parameters, rep(tune.num.trees, tune.num.reps)) # also put tune.num.trees into draw.parameters for row-wise operation below
  draw.parameters <- cbind(draw.parameters, rep("trt", tune.num.reps)) # also put type into draw.parameters for row-wise operation below

  # 1. Obtain forest errors (average debiased error across observations)
  small.forest.errors <- apply(draw.parameters, 1, get_forest_debiased_error)

  # 2. Fit the 'dice kriging' model to these error estimates.
  kriging.model <- get_kriging_model(fit.draws, small.forest.errors)

  # 3. To determine the optimal parameter values, predict using the kriging model at a large
  # number of random values, then select those that produced the lowest error.
  optimize.draws <- get_fit_draws(tune.num.draws)
  model.surface <- predict(kriging.model, newdata = data.frame(optimize.draws), type = "SK")$mean
  tuned.params <- get_draw_parameters(optimize.draws)

  grid <- cbind(error = c(model.surface), tuned.params)
  small.forest.optimal.draw <- which.min(grid[, "error"])

  # To avoid the possibility of selection bias, re-train a moderately-sized forest
  # at the value chosen by the method above
  retrained.forest.params <- grid[small.forest.optimal.draw, -1]
  retrained.forest.error <- get_forest_debiased_error(c(retrained.forest.params, moderate.num.trees, "trt"))

  # 4. Train a forest with default parameters, and check its predicted error.
  s_default <- scale*ceiling(nrow(trt)^lower_bound_default)/nrow(trt)
  forestT_default <- regression_forest(feature_trt, trt$Y, mtry=1, num.trees=B, sample.fraction=s_default)
  default.forest.error <- mean(forestT_default$debiased.error)

  # After tuning, run forests with tuned/default parameters
  if (default.forest.error < retrained.forest.error) {
    predictionT <- as.numeric(predict(forestT_default, focal))
    BLB1 <- append(BLB1, predict(forestT_default, focal, estimate.variance = TRUE)$variance.estimates)
  } else {
    a <- as.double(retrained.forest.params[1])
    imb <- as.double(retrained.forest.params[2])
    h <- as.double(retrained.forest.params[3])
    prune <- as.double(retrained.forest.params[4])
    lower_bound <- (1+ log(1-a) / log(a))^(-1)  # d=1
    s <- scale*ceiling(nrow(trt)^lower_bound)/nrow(trt)

    forestT <- regression_forest(feature_trt, trt$Y, mtry=1, num.trees=B, sample.fraction=s, honesty.fraction = h, alpha = a,
                                 imbalance.penalty = imb, honesty.prune.leaves = prune)
    predictionT <- as.numeric(predict(forestT, focal))
    BLB1 <- append(BLB1, predict(forestT, focal, estimate.variance = TRUE)$variance.estimates)
  }

  ## Parameter tuning for control
  fit.draws <- get_fit_draws(tune.num.reps)
  draw.parameters <- get_draw_parameters(fit.draws)
  draw.parameters <- cbind(draw.parameters, rep(tune.num.trees, tune.num.reps)) # also put tune.num.trees into draw.parameters for row-wise operation below
  draw.parameters <- cbind(draw.parameters, rep("ctrl", tune.num.reps)) # also put type into draw.parameters for row-wise operation below

  # 1. Obtain forest errors (average debiased error across observations)
  small.forest.errors <- apply(draw.parameters, 1, get_forest_debiased_error)

  # 2. Fit the 'dice kriging' model to these error estimates.
  kriging.model <- get_kriging_model(fit.draws, small.forest.errors)

  # 3. To determine the optimal parameter values, predict using the kriging model at a large
  # number of random values, then select those that produced the lowest error.
  optimize.draws <- get_fit_draws(tune.num.draws)
  model.surface <- predict(kriging.model, newdata = data.frame(optimize.draws), type = "SK")$mean
  tuned.params <- get_draw_parameters(optimize.draws)

  grid <- cbind(error = c(model.surface), tuned.params)
  small.forest.optimal.draw <- which.min(grid[, "error"])

  # To avoid the possibility of selection bias, re-train a moderately-sized forest
  # at the value chosen by the method above
  retrained.forest.params <- grid[small.forest.optimal.draw, -1]
  retrained.forest.error <- get_forest_debiased_error(c(retrained.forest.params, moderate.num.trees, "ctrl"))

  # 4. Train a forest with default parameters, and check its predicted error.
  # This improves our chances of not doing worse than default
  s_default <- scale*ceiling(nrow(ctrl)^lower_bound_default)/nrow(ctrl)
  forestC_default <- regression_forest(feature_ctrl, ctrl$Y, mtry=1, num.trees=B, sample.fraction=s_default)
  default.forest.error <- mean(forestC_default$debiased.error)

  # After tuning, run forests with tuned/default parameters
  if (default.forest.error < retrained.forest.error) {
    predictionC <- as.numeric(predict(forestC_default, focal))
    BLB0 <- append(BLB0, predict(forestC_default, focal, estimate.variance = TRUE)$variance.estimates)
  } else {
    a <- as.double(retrained.forest.params[1])
    imb <- as.double(retrained.forest.params[2])
    h <- as.double(retrained.forest.params[3])
    prune <- as.double(retrained.forest.params[4])
    lower_bound <- (1+ log(1-a) / log(a))^(-1)  # d=1
    s <- scale*ceiling(nrow(ctrl)^lower_bound)/nrow(ctrl)

    forestC <- regression_forest(feature_ctrl, ctrl$Y, mtry=1, num.trees=B, sample.fraction=s, honesty.fraction = h, alpha = a,
                                 imbalance.penalty = imb, honesty.prune.leaves = prune)

    predictionC <- as.numeric(predict(forestC, focal))
    BLB0 <- append(BLB0, predict(forestC, focal, estimate.variance = TRUE)$variance.estimates)
  }

  # get estimate for treatment effect
  prediction <- predictionT - predictionC
  result_rf <- append(result_rf, prediction)

  variance_BLB <- BLB1[i]+BLB0[i]
  se_BLB <- sqrt(variance_BLB)

  if ((prediction-1.96*se_BLB <= truth) & (truth <= prediction+1.96*se_BLB)){
    coverage_rf <- append(coverage_rf, 1)
  }
  else{
    coverage_rf <- append(coverage_rf, 0)
  }

  try({
    print("----------------")
    print(paste0("RF: ", prediction))
    print("----------------")
    print(paste0("Mean of RF:", mean(result_rf, na.rm = TRUE)))
    print("----------------")
    print(paste0("Var of RF: ", var(result_rf, na.rm = TRUE)))
    print(paste0("Var of RF (BLB): ", variance_BLB))
    print("----------------")
    print(paste0("MSE of RF: ", (mean(result_rf,na.rm = TRUE)-truth)^2+var(result_rf,na.rm = TRUE)))
    print("----------------")
    print(paste0("CR of RF: ", mean(coverage_rf, na.rm = TRUE)))
    print("----------------")
  })
}

