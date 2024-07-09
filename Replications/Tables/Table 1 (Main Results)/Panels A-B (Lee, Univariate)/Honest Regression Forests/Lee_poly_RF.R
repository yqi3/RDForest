######################################
# Table 1, Panel A, Columns 11-15: Honest Regression Forest
# Simulation: Univariate score using data from polynomial DGP

# We implement random forests using grf on the following platform
# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux trixie/sid
######################################
rm(list=ls())
set.seed(1234)

library(grf)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../../../Helper Fns/KT_Helper_Fns.R")

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
truthT <- c_r0
truthC <- c_l0
focal <- matrix(0)

#### MC Simulation ####
t <- 1000 # number of Monte Carlo iterations

# Number of observations in each Monte Carlo iteration \in {1000,5000,10000,50000}
for (n in c(1000,5000,10000,50000,1e5)){
  result_rf <- c()
  bias_rf <- c()
  coverage_rf <- c()
  BLB1 <- c()
  BLB0 <- c()
  result_trt <- c()
  result_ctrl <- c()
  coverage_trt <- c()
  coverage_ctrl <- c()
  
  for (i in 1:t) {
    print(i)
    # generate simulation data
    X <- 2*rbeta(n,2,4)-1
    # epsilon <- rnorm(n, 0, sigma_y*(1/sqrt(abs(X))))
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
    
    forestT <- tuned_RDForest(feature_trt, trt$Y, tune.parameters=c(), type = "trt")
    predictionT <- predict(forestT, focal, estimate.variance = TRUE)
    BLB1 <- append(BLB1, predictionT$variance.estimates)
    
    forestC <- tuned_RDForest(feature_ctrl, ctrl$Y, tune.parameters=c(), type = "ctrl")
    predictionC <- predict(forestC, focal, estimate.variance = TRUE)
    BLB0 <- append(BLB0, predictionC$variance.estimates)
    
    # get estimate for treatment effect
    result_trt <- append(result_trt, predictionT$predictions)
    result_ctrl <- append(result_ctrl, predictionC$predictions)
    prediction <- predictionT$predictions - predictionC$predictions
    result_rf <- append(result_rf, prediction)
    bias_rf <- append(bias_rf, prediction-truth)
    
    variance_BLB <- BLB1[i]+BLB0[i]
    se_BLB <- sqrt(variance_BLB)
    
    if ((predictionT$predictions-1.96*sqrt(BLB1[i]) <= truthT) & (truthT <= predictionT$predictions+1.96*sqrt(BLB1[i]))){
      coverage_trt <- append(coverage_trt, 1)
    } else{
      coverage_trt <- append(coverage_trt, 0)
    }
    
    if ((predictionC$predictions-1.96*sqrt(BLB0[i]) <= truthC) & (truthC <= predictionC$predictions+1.96*sqrt(BLB0[i]))){
      coverage_ctrl <- append(coverage_ctrl, 1)
    } else{
      coverage_ctrl <- append(coverage_ctrl, 0)
    }
    
    if ((prediction-1.96*se_BLB <= truth) & (truth <= prediction+1.96*se_BLB)){
      coverage_rf <- append(coverage_rf, 1)
    } else{
      coverage_rf <- append(coverage_rf, 0)
    }
    
    try({
      print(paste0("n=",n))
      print("----------------")
      print(paste0("RF: ", prediction))
      print(paste0("RF trt: ", predictionT$predictions))
      print(paste0("RF ctrl: ", predictionC$predictions))
      print("----------------")
      print(paste0("Mean of RF:", mean(result_rf, na.rm = TRUE)))
      print(paste0("Bias of RF:", mean(bias_rf, na.rm = TRUE)))
      print(paste0("Mean of RF trt:", mean(result_trt, na.rm = TRUE)))
      print(paste0("Mean of RF ctrl:", mean(result_ctrl, na.rm = TRUE)))
      print("----------------")
      print(paste0("Var of RF: ", var(result_rf, na.rm = TRUE)))
      print(paste0("Var of RF trt: ", var(result_trt, na.rm = TRUE)))
      print(paste0("Var of RF ctrl: ", var(result_ctrl, na.rm = TRUE)))
      print(paste0("Var of RF (BLB): ", variance_BLB))
      print("----------------")
      print(paste0("CR (truth=", truth, "): ", mean(coverage_rf, na.rm = TRUE)))
      print(paste0("CR trt (truthT=", truthT, "): ", mean(coverage_trt, na.rm = TRUE)))
      print(paste0("CR ctrl (truthT=", truthC, "): ", mean(coverage_ctrl, na.rm = TRUE)))
      print("----------------")
    })
  }
  assign(paste0("result",n), as.data.frame(cbind("n"=n, "rf"=result_rf, "rf.bias"=bias_rf, "rf.cr"=coverage_rf, "rf.se"=sqrt(BLB0+BLB1))))
}

result <- rbind(result1000,result5000,result10000,result50000,`result1e+05`)
result$truth <- truth
write.csv(result, "../Results csv/poly-rf-lee.csv", row.names = F)
