######################################
# Figures B.9: Comparisons of different
# choices of the scaling constant using the
# WGAN DGPs, outcome=democratic vote share
######################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(foreign)
library(grf)
set.seed(1234)
source("../../../Helper Fns/KT_Helper_Fns.R")

#### Lee WGAN ####
data_full <- read.csv("../../../Data/Lee2008/generated_Lee.csv")
focal_df <- read.csv("../../../Data/Lee2008/focal_lee.csv")
truthT <- mean(focal_df$Y[focal_df$W==1])
truthC <- mean(focal_df$Y[focal_df$W==0])
truth <- trutT - truthC
focal <- matrix(0)

## MC Simulation
t <- 1000 # number of Monte Carlo iterations
n <- 5000

for (scale in c(0.2, 0.3, 0.4, 0.5)) {
  result_rf <- c()
  bias_rf <- c()
  coverage_rf <- c()
  
  for (i in 1:t) {
    print(i)
    # generate simulation data
    rand <- sample.int(nrow(data_full), n, replace = F)
    datause <- data_full[rand,]
    datause <- datause[,c("X","W","Y")]
    
    ### RF ###
    trt <- datause[datause$W==1,]
    ctrl <- datause[datause$W==0,]
    feature_trt <- as.data.frame(trt$X)
    feature_ctrl <- as.data.frame(ctrl$X)
    
    # -----------------------------------------------------------------------------------
    # Modified tuning procedure based on the source code of grf; see
    # https://github.com/grf-labs/grf/blob/master/r-package/grf/R/tune_forest.R
    # -----------------------------------------------------------------------------------
    
    forestT <- tuned_RDForest(feature_trt, trt$Y, tune.parameters=c(), scale = scale, type = "trt")
    predictionT <- predict(forestT, focal, estimate.variance = TRUE)
    BLB1 <- predictionT$variance.estimates
    
    forestC <- tuned_RDForest(feature_ctrl, ctrl$Y, tune.parameters=c(), scale = scale, type = "ctrl")
    predictionC <- predict(forestC, focal, estimate.variance = TRUE)
    BLB0 <- predictionC$variance.estimates
    
    # get estimate for treatment effect
    prediction <- predictionT$predictions - predictionC$predictions
    result_rf <- append(result_rf, prediction)
    bias_rf <- append(bias_rf, prediction-truth)
    
    variance_BLB <- BLB1+BLB0
    se_BLB <- sqrt(variance_BLB)
    
    if ((prediction-1.96*se_BLB <= truth) & (truth <= prediction+1.96*se_BLB)){
      coverage_rf <- append(coverage_rf, 1)
    } else{
      coverage_rf <- append(coverage_rf, 0)
    }
    
    try({
      print(paste0("scale=",scale))
      print(paste0("n=",n))
      print("----------------")
      print(paste0("RF: ", prediction))
      print(paste0("Mean of RF:", mean(result_rf, na.rm = TRUE)))
      print(paste0("Bias of RF:", mean(bias_rf, na.rm = TRUE)))
      print("----------------")
      print(paste0("Var of RF: ", var(result_rf, na.rm = TRUE)))
      print(paste0("Var of RF (BLB): ", variance_BLB))
      print("----------------")
      print(paste0("CR (truth=", truth, "): ", mean(coverage_rf, na.rm = TRUE)))
      print("----------------")
    })
  }
  
  assign(paste0("result",scale), as.data.frame(cbind("scale"=scale, "rf"=result_rf, "rf.bias"=bias_rf, "rf.bias.sq"=(bias_rf)^2, "rf.cr"=coverage_rf, "rf.se"=sqrt(BLB0+BLB1))))
}

## All default
result_rf <- c()
bias_rf <- c()
coverage_rf <- c()

for (i in 1:t) {
  print(i)
  # generate simulation data
  rand <- sample.int(nrow(data_full), n, replace = F)
  datause <- data_full[rand,]
  datause <- datause[,c("X","W","Y")]
  
  ### RF ###
  trt <- datause[datause$W==1,]
  ctrl <- datause[datause$W==0,]
  feature_trt <- as.data.frame(trt$X)
  feature_ctrl <- as.data.frame(ctrl$X)
  
  # -----------------------------------------------------------------------------------
  # Modified tuning procedure based on the source code of grf; see
  # https://github.com/grf-labs/grf/blob/master/r-package/grf/R/tune_forest.R
  # -----------------------------------------------------------------------------------
  
  forestT <- regression_forest(feature_trt, trt$Y, num.trees = 5000, mtry = 1)
  predictionT <- predict(forestT, focal, estimate.variance = TRUE)
  BLB1 <- predictionT$variance.estimates
  
  forestC <- regression_forest(feature_ctrl, ctrl$Y, num.trees = 5000, mtry = 1)
  predictionC <- predict(forestC, focal, estimate.variance = TRUE)
  BLB0 <- predictionC$variance.estimates
  
  # get estimate for treatment effect
  prediction <- predictionT$predictions - predictionC$predictions
  result_rf <- append(result_rf, prediction)
  bias_rf <- append(bias_rf, prediction-truth)
  
  variance_BLB <- BLB1+BLB0
  se_BLB <- sqrt(variance_BLB)
  
  if ((prediction-1.96*se_BLB <= truth) & (truth <= prediction+1.96*se_BLB)){
    coverage_rf <- append(coverage_rf, 1)
  } else{
    coverage_rf <- append(coverage_rf, 0)
  }
  
  try({
    print("scale=default")
    print(paste0("n=",n))
    print("----------------")
    print(paste0("RF: ", prediction))
    print(paste0("Mean of RF:", mean(result_rf, na.rm = TRUE)))
    print(paste0("Bias of RF:", mean(bias_rf, na.rm = TRUE)))
    print("----------------")
    print(paste0("Var of RF: ", var(result_rf, na.rm = TRUE)))
    print(paste0("Var of RF (BLB): ", variance_BLB))
    print("----------------")
    print(paste0("CR (truth=", truth, "): ", mean(coverage_rf, na.rm = TRUE)))
    print("----------------")
  })
}

assign(paste0("result_default"), as.data.frame(cbind("scale"="default", "rf"=result_rf, "rf.bias"=bias_rf, "rf.cr"=coverage_rf, "rf.se"=sqrt(BLB0+BLB1))))

result <- rbind(result0.2,result0.3,result0.4,result0.5,result_default)
result$truth <- truth
write.csv(result, paste0("../Scale Data/wgan-rf-lee-scale.csv"), row.names = F)
