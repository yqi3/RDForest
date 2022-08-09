######################################
# Simulation: Univariate score using
# data from Lee (2008)
######################################
library(grf)

set.seed(1234)

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

t <- 1000 
n <- 5000

result_rf <- c() # vector to collect results
coverage_rf_nocov <- c()
BLB1 <- c()
BLB0 <- c()
focal <- matrix(0)


for (i in 1:t) {
  print(i)
  
  X <- 2*rbeta(n,2,4)-1
  epsilon <- rnorm(n, 0, sigma_y)
  Y <- ifelse(X>=0, c_r0+c_r1*X+c_r2*X^2+c_r3*X^3+c_r4*X^4+c_r5*X^5+epsilon, 
              c_l0+c_l1*X+c_l2*X^2+c_l3*X^3+c_l4*X^4+c_l5*X^5+epsilon)
  
  data <- as.data.frame(cbind("Y" = Y, "X"=X))
  
  ############################
  # RF
  ############################
  trt <- data[data$X>=0,]
  ctrl <- data[data$X<0,]
  feature_trt <- as.data.frame(trt[,2])
  feature_ctrl <- as.data.frame(ctrl[,2])
  
  # -------------------------------------------------
  #- same subsample rate (0.5) for treated and control 
  # -------------------------------------------------
  forestT <- regression_forest(feature_trt, trt$Y, mtry=1, num.trees = 5000, tune.parameters = c("alpha", "imbalance.penalty"))
  forestC <- regression_forest(feature_ctrl, ctrl$Y, mtry=1, num.trees = 5000, tune.parameters = c("alpha", "imbalance.penalty"))
  
  # get estimate for tau
  prediction <- as.numeric(predict(forestT, focal)-predict(forestC, focal))
  result_rf <- append(result_rf, prediction)
  
  #- coverage rate for BLB variance estimator under zero covariance
  BLB1 <- append(BLB1, predict(forestT, focal, estimate.variance = TRUE)$variance.estimates)
  BLB0 <- append(BLB0, predict(forestC, focal, estimate.variance = TRUE)$variance.estimates)
  variance_nocov <- BLB1[i]+BLB0[i]
  se_nocov <- sqrt(variance_nocov)
  
  if ((prediction-1.96*se_nocov <= truth) & (truth <= prediction+1.96*se_nocov)){
    coverage_rf_nocov <- append(coverage_rf_nocov, 1)
  }
  else{
    coverage_rf_nocov <- append(coverage_rf_nocov, 0)
  }
  
  try({
    print(paste0("Prediction: ", prediction))
    print(paste0("Mean:", mean(result_rf)))
    print(paste0("Var:", var(result_rf)))
    print(paste0("Var of RF (nocov): ", variance_nocov))
    print(paste0("MSE: ", (mean(result_rf)-truth)^2+var(result_rf)))
    print(paste0("CR of RF (nocov): ", mean(coverage_rf_nocov, na.rm = TRUE)))
    print("----------------")
  })
}