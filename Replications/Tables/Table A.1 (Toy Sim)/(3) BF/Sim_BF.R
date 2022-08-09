# ==========================
# Toy Simulations
# BF without covariates, 1-3 running vars
# ==========================
library(grf)

# DGP
# ==============================
set.seed(1234)
n = 1000 # sample size
t <- 1000 # number of iterations
tau = 0.3 # true treatment effect for the linear DGP
result_univariate <- c() # vector to collect results
result_bivariate <- c()
result_trivariate <- c()

for (i in 1:t) {
  # running variable(s): X_i ~ uniform[-1,1]^p
  X1 <- runif(n, -1, 1)
  X2 <- runif(n, -1, 1)
  X3 <- runif(n, -1, 1)
  
  # treatment indicator: D_i = 1 if X_1 >= 0
  D1 <- ifelse(X1 >= 0, 1, 0) # univariate running var
  D2 <- ifelse(X1 >= 0 & X2 >= 0, 1, 0) # multivariate running var (2)
  D3 <- ifelse(X1 >= 0 & X2 >= 0 & X3 >= 0, 1, 0) # multivariate running var (3)
  
  # Generate outcome Y as Bernoulli trials with probability 0.2 + tau*D_i
  Y1 <- rbinom(n, 1, 0.2+tau*D1)
  Y2 <- rbinom(n, 1, 0.2+tau*D2)
  Y3 <- rbinom(n, 1, 0.2+tau*D3)
  
  # Dataframe for running vars and response only
  mydata <- as.data.frame(cbind("Y1"= Y1, "Y2" = Y2, "Y3" = Y3, 
                                "X1"=X1, "X2"=X2, "X3"=X3, "D1"=D1, "D2"=D2, "D3"=D3))
  
  ### Univariate case
  trt <- mydata[mydata$D1==1,]
  ctrl <- mydata[mydata$D1==0,]
  feature_trt <- as.data.frame(trt[,4])
  feature_ctrl <- as.data.frame(ctrl[,4])
  cutoff <- matrix(0)
  
  boostT_univariate <- boosted_regression_forest(feature_trt, trt$Y1, num.trees = 5000, mtry=1, boost.steps=1)
  boostC_univariate <- boosted_regression_forest(feature_ctrl, ctrl$Y1, num.trees = 5000, mtry=1, boost.steps=1)
  BF_univariate_pred <- as.numeric(predict(boostT_univariate, cutoff)-predict(boostC_univariate, cutoff))
  result_univariate <- append(result_univariate, BF_univariate_pred)

  ### Bivariate case
  trt <- mydata[mydata$D2==1,]
  ctrl <- mydata[mydata$D2==0,]
  feature_trt <- as.data.frame(trt[,c(4,5)])
  feature_ctrl <- as.data.frame(ctrl[,c(4,5)])
  cutoff <- matrix(c(0,0), ncol=2)
  
  boostT_bivariate <- boosted_regression_forest(feature_trt, trt$Y2, num.trees = 5000, mtry=1, boost.steps=1)
  boostC_bivariate <- boosted_regression_forest(feature_ctrl, ctrl$Y2, num.trees = 5000, mtry=1, boost.steps=1)
  BF_bivariate_pred <- as.numeric(predict(boostT_bivariate, cutoff)-predict(boostC_bivariate, cutoff))
  result_bivariate <- append(result_bivariate, BF_bivariate_pred)
  
  ### Trivariate case
  trt <- mydata[mydata$D3==1,]
  ctrl <- mydata[mydata$D3==0,]
  feature_trt <- as.data.frame(trt[,c(4,5,6)])
  feature_ctrl <- as.data.frame(ctrl[,c(4,5,6)])
  cutoff <- matrix(c(0,0,0), ncol=3)
  
  boostT_trivariate <- boosted_regression_forest(feature_trt, trt$Y3, num.trees = 5000, mtry=1, boost.steps=1)
  boostC_trivariate <- boosted_regression_forest(feature_ctrl, ctrl$Y3, num.trees = 5000, mtry=1, boost.steps=1)
  BF_trivariate_pred <- as.numeric(predict(boostT_trivariate, cutoff)-predict(boostC_trivariate, cutoff))
  result_trivariate <- append(result_trivariate, BF_trivariate_pred)
  
  try({print(i)
    print("----------------")
    print(paste0("Univariate: ", BF_univariate_pred))
    print(paste0("Mean of univariate:", mean(result_univariate)))
    print(paste0("Var of univariate:", var(result_univariate)))
    print(paste0("MSE of univariate: ", (mean(result_univariate)-tau)^2+var(result_univariate)))
    print("----------------")
    print(paste0("Bivariate: ", BF_bivariate_pred))
    print(paste0("Mean of bivariate:", mean(result_bivariate)))
    print(paste0("Var of bivariate:", var(result_bivariate)))
    print(paste0("MSE of bivariate: ", (mean(result_bivariate)-tau)^2+var(result_bivariate)))
    print("----------------")
    print(paste0("Trivariate: ", BF_trivariate_pred))
    print(paste0("Mean of trivariate:", mean(result_trivariate)))
    print(paste0("Var of trivariate:", var(result_trivariate)))
    print(paste0("MSE of trivariate: ", (mean(result_trivariate)-tau)^2+var(result_trivariate)))
    print("----------------")
  })
}