######################################
# Figures B.8: Comparisons of different
# choices of the scaling constant using the
# polynomial DGPs, outcome=turnout
######################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(foreign)
library(grf)
set.seed(1234)
source("../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Analysis/Local Poly/distance-functions.R")
source("../../../Helper Fns/KT_Helper_Fns.R")

### Read sample data ###
library(haven)
data <- read_dta("../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Data/Voters/Voters_Final.dta")
data <- data[complete.cases(data$latitude, data$longitude, data$treat, data[,c("e2008g")]),]

### Get focal point on the treatment boundary where we estimate treatment effect ###
pointsALL <- read.dbf("../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Data/BorderSegmentPoints_Project.dbf") # boundary data
pointsALL$latitude= pointsALL$POINT_Y
pointsALL$longitude= pointsALL$POINT_X
focal <- matrix(c(pointsALL$latitude[45],pointsALL$longitude[45]),nrow = 1)

X1 <- data$latitude
X2 <- data$longitude
mean_X1 <- mean(X1)
mean_X2 <- mean(X2)
e2008g <- data$e2008g

model.df <- data.frame(Y=e2008g, X1=X1, X1.2=(X1 - mean_X1)^2, X1.3=(X1 - mean_X1)^3,
                       X2=X2, X2.2=(X2 - mean_X2)^2, X2.3=(X2 - mean_X2)^3,
                       X12=(X1 - mean_X1)*(X2 - mean_X2),
                       W=data$treat)

cef <- glm(Y ~ W * ., family=binomial(link='logit'),data=model.df)
sigma_y <- sqrt(summary(cef)$sigma^2)

X1 <- pointsALL$latitude[45]
X2 <- pointsALL$longitude[45]
truth <- predict(cef, data.frame(X1=X1, X1.2=(X1 - mean_X1)^2, X1.3=(X1 - mean_X1)^3,
                                 X2=X2, X2.2=(X2 - mean_X2)^2, X2.3=(X2 - mean_X2)^3,
                                 X12=(X1 - mean_X1)*(X2 - mean_X2), W=1), type = "response") - predict(cef, data.frame(X1=X1, X1.2=(X1 - mean_X1)^2, X1.3=(X1 - mean_X1)^3,
                                                                                                                       X2=X2, X2.2=(X2 - mean_X2)^2, X2.3=(X2 - mean_X2)^3,
                                                                                                                       X12=(X1 - mean_X1)*(X2 - mean_X2), W=0), type = "response")

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
    datause <- genDGP(n, "turnout")
    
    ### RF ###
    trt <- datause[datause$treat==1,]
    ctrl <- datause[datause$treat==0,]
    feature_trt <- as.data.frame(trt[,c(1,2)])
    feature_ctrl <- as.data.frame(ctrl[,c(1,2)])
    
    # -----------------------------------------------------------------------------------
    # Modified tuning procedure based on the source code of grf; see
    # https://github.com/grf-labs/grf/blob/master/r-package/grf/R/tune_forest.R
    # -----------------------------------------------------------------------------------
    
    forestT <- tuned_RDForest(feature_trt, trt[,4], tune.parameters=c(), scale = scale, type = "trt")
    predictionT <- predict(forestT, focal, estimate.variance = TRUE)
    BLB1 <- predictionT$variance.estimates
    
    forestC <- tuned_RDForest(feature_ctrl, ctrl[,4], tune.parameters=c(), scale = scale, type = "ctrl")
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
  
  assign(paste0("result",scale), as.data.frame(cbind("scale"=scale, "rf"=result_rf, "rf.bias"=bias_rf, "rf.cr"=coverage_rf, "rf.se"=sqrt(BLB0+BLB1))))
}

## All default
result_rf <- c()
bias_rf <- c()
coverage_rf <- c()

for (i in 1:t) {
  print(i)
  # generate simulation data
  datause <- genDGP(n, "turnout")
  
  ### RF ###
  trt <- datause[datause$treat==1,]
  ctrl <- datause[datause$treat==0,]
  feature_trt <- as.data.frame(trt[,c(1,2)])
  feature_ctrl <- as.data.frame(ctrl[,c(1,2)])
  
  # -----------------------------------------------------------------------------------
  # Modified tuning procedure based on the source code of grf; see
  # https://github.com/grf-labs/grf/blob/master/r-package/grf/R/tune_forest.R
  # -----------------------------------------------------------------------------------
  
  forestT <- regression_forest(feature_trt, trt[,4], num.trees = 5000, mtry = 1)
  predictionT <- predict(forestT, focal, estimate.variance = TRUE)
  BLB1 <- predictionT$variance.estimates
  
  forestC <- regression_forest(feature_ctrl, ctrl[,4], num.trees = 5000, mtry = 1)
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
write.csv(result, paste0("../Scale Data/poly-rf-T-scale.csv"), row.names = F)
