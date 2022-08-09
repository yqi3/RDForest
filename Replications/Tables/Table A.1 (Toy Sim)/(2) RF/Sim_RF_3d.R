# ==========================
# Toy Simulations
# RF simulations with BLB variance estimates without covariance, trivariate
# ==========================
library(grf)

# DGP
# ==============================
set.seed(1234)
n = 1000 # sample size
t <- 1000 # number of iterations
tau = 0.3 # true treatment effect for the linear DGP
result_trivariate <- c()
coverage_rf_nocov <- c()
BLB1 <- c()
BLB0 <- c()

for (i in 1:t) {
  # running variable(s): X_i ~ uniform[-1,1]^p
  X1 <- runif(n, -1, 1)
  X2 <- runif(n, -1, 1)
  X3 <- runif(n, -1, 1)
  
  # treatment indicator: D_i = 1 if X_1 >= 0
  D3 <- ifelse(X1 >= 0 & X2 >= 0 & X3 >= 0, 1, 0) # multivariate running var (3)
  
  # Generate outcome Y as Bernoulli trials with probability 0.2 + tau*D_i
  Y3 <- rbinom(n, 1, 0.2+tau*D3)
  
  # Dataframe for running vars and response only
  mydata <- as.data.frame(cbind("Y3" = Y3, "X1"=X1, "X2"=X2, "X3"=X3, "D3"=D3))
  cutoff <- matrix(c(0,0,0), ncol=3)
  
  ### Trivariate case
  trt <- mydata[mydata$D3==1,]
  ctrl <- mydata[mydata$D3==0,]
  feature_trt <- as.data.frame(trt[,c(2,3,4)])
  feature_ctrl <- as.data.frame(ctrl[,c(2,3,4)])
  
  forestT <- regression_forest(feature_trt, trt$Y3, mtry=1, num.trees=5000, tune.parameters = c("sample.fraction", "honesty.fraction", "alpha"))
  forestC <- regression_forest(feature_ctrl, ctrl$Y3, mtry=1, num.trees=5000, tune.parameters = c("sample.fraction", "honesty.fraction", "alpha"))
  
  # get estimated tau at the cutoff
  RF_trivariate_pred <- as.numeric(predict(forestT, cutoff)-predict(forestC, cutoff))
  result_trivariate <- append(result_trivariate, RF_trivariate_pred)
  
  #- coverage rate for BLB variance estimator under zero covariance
  BLB1 <- append(BLB1, predict(forestT, cutoff, estimate.variance = TRUE)$variance.estimates)
  BLB0 <- append(BLB0, predict(forestC, cutoff, estimate.variance = TRUE)$variance.estimates)
  variance_nocov <- BLB1[i]+BLB0[i]
  se_nocov <- sqrt(variance_nocov)
  
  if ((RF_trivariate_pred-1.96*se_nocov <= tau) & (tau <= RF_trivariate_pred+1.96*se_nocov)){
    coverage_rf_nocov <- append(coverage_rf_nocov, 1)
  }
  else{
    coverage_rf_nocov <- append(coverage_rf_nocov, 0)
  }
  
  try({print(i)
    print(paste0("Univariate: ", RF_trivariate_pred))
    print(paste0("Mean of trivariate:", mean(result_trivariate)))
    print(paste0("Var of trivariate:", var(result_trivariate)))
    print(paste0("Var of RF (nocov): ", variance_nocov))
    print(paste0("MSE of trivariate: ", (mean(result_trivariate)-tau)^2+var(result_trivariate)))
    print(paste0("CR of RF (nocov): ", mean(coverage_rf_nocov, na.rm = TRUE)))
    print("----------------")
  })
}
