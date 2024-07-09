######################################
# Table 1, Panel F, Columns 11-15: Honest Regression Forest for Turnout
# Simulation: Bivariate score using WGAN data based on Keele & Titiunik (2015)

# We implement random forests using grf on the following platform
# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux trixie/sid
######################################
rm(list=ls())
set.seed(1234)

library(foreign)
library(grf)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Analysis/Local Poly/distance-functions.R")
source("../../../../Helper Fns/KT_Helper_Fns.R")

### Get focal point on the treatment boundary where we estimate treatment effect ###
pointsALL <- read.dbf("../../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Data/BorderSegmentPoints_Project.dbf") # boundary data
pointsALL$latitude <- pointsALL$POINT_Y
pointsALL$longitude <- pointsALL$POINT_X
focal <- matrix(c(pointsALL$latitude[45],pointsALL$longitude[45]),nrow = 1)

### Read WGAN-generated population data
data_full <- read.csv("../../../../Data/Keele & Titiunik 2015/WGAN_Epoch1k/turnout_cWGAN.csv")
focal_df <- read.csv("../../../../Data/Keele & Titiunik 2015/WGAN_Epoch1k/focal_turnout.csv")

### Population truth
truthT <- mean(focal_df$Y[focal_df$W==1])
truthC <- mean(focal_df$Y[focal_df$W==0])
truth <- truthT - truthC

#### MC Simulation ####
t <- 1000 # Number of Monte Carlo iterations

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
    rand <- sample.int(nrow(data_full), n, replace = F)
    datause <- data_full[rand,]
    datause <- datause[,c("X1","X2","W","Y")]
    
    ### RF ###
    trt <- datause[datause$W==1,]
    ctrl <- datause[datause$W==0,]
    feature_trt <- as.data.frame(trt[,c(1,2)])
    feature_ctrl <- as.data.frame(ctrl[,c(1,2)])
    
    # -----------------------------------------------------------------------------------
    # Modified tuning procedure based on the source code of grf; see
    # https://github.com/grf-labs/grf/blob/master/r-package/grf/R/tune_forest.R
    # -----------------------------------------------------------------------------------
    
    forestT <- tuned_RDForest(feature_trt, trt[,4], tune.parameters=c(), type = "trt")
    predictionT <- predict(forestT, focal, estimate.variance = TRUE)
    BLB1 <- append(BLB1, predictionT$variance.estimates)
    forestC <- tuned_RDForest(feature_ctrl, ctrl[,4], tune.parameters=c(), type = "ctrl")
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
    } else {
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
write.csv(result, "../Results csv/wgan-rf-turnout.csv", row.names = F)
