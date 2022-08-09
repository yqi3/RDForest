# ==========================
# Toy Simulations
# RF simulations with BLB variance estimates without covariance, univariate
# ==========================
library(grf)

# DGP
# ==============================
set.seed(1234)
n = 1000 # sample size
t <- 1000 # number of iterations
tau = 0.3 # true treatment effect for the linear DGP
result_univariate <- c() # vector to collect results
coverage_rf_nocov <- c()
BLB1 <- c()
BLB0 <- c()

for (i in 1:t) {
  # running variable(s): X_i ~ uniform[-1,1]^p
  X1 <- runif(n, -1, 1)
  
  # treatment indicator: D_i = 1 if X_1 >= 0
  D1 <- ifelse(X1 >= 0, 1, 0) # univariate running var
  
  # Generate outcome Y as Bernoulli trials with probability 0.2 + tau*D_i
  Y1 <- rbinom(n, 1, 0.2+tau*D1)
  
  # Dataframe for running vars and response only
  mydata <- as.data.frame(cbind("Y1"= Y1, "X1"=X1, "D1"=D1))
  
  cutoff <- matrix(0)
  
  ### Univariate case
  trt <- mydata[mydata$D1==1,]
  ctrl <- mydata[mydata$D1==0,]
  feature_trt <- as.data.frame(trt[,c(2)])
  feature_ctrl <- as.data.frame(ctrl[,c(2)])
  
  # separate forests on each side of the cutoff
  forestT <- regression_forest(feature_trt, trt$Y1, mtry=1, num.trees=5000, tune.parameters = c("sample.fraction", "honesty.fraction", "alpha"))
  forestC <- regression_forest(feature_ctrl, ctrl$Y1, mtry=1, num.trees=5000, tune.parameters = c("sample.fraction", "honesty.fraction", "alpha"))
  
  # get estimated tau at the cutoff
  RF_univariate_pred <- as.numeric(predict(forestT, cutoff)-predict(forestC, cutoff))
  result_univariate <- append(result_univariate, RF_univariate_pred)
  
  #- coverage rate for BLB variance estimator under zero covariance
  BLB1 <- append(BLB1, predict(forestT, cutoff, estimate.variance = TRUE)$variance.estimates)
  BLB0 <- append(BLB0, predict(forestC, cutoff, estimate.variance = TRUE)$variance.estimates)
  variance_nocov <- BLB1[i]+BLB0[i]
  se_nocov <- sqrt(variance_nocov)
  
  if ((RF_univariate_pred-1.96*se_nocov <= tau) & (tau <= RF_univariate_pred+1.96*se_nocov)){
    coverage_rf_nocov <- append(coverage_rf_nocov, 1)
  }
  else{
    coverage_rf_nocov <- append(coverage_rf_nocov, 0)
  }
  
  try({print(i)
    print(paste0("Univariate: ", RF_univariate_pred))
    print(paste0("Mean of univariate:", mean(result_univariate)))
    print(paste0("Var of univariate:", var(result_univariate)))
    print(paste0("Var of RF (nocov): ", variance_nocov))
    print(paste0("MSE of univariate: ", (mean(result_univariate)-tau)^2+var(result_univariate)))
    print(paste0("CR of RF (nocov): ", mean(coverage_rf_nocov, na.rm = TRUE)))
    print("----------------")
  })
}
