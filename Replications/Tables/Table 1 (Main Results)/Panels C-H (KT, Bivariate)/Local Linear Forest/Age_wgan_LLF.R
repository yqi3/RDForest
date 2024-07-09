######################################
# Table 1, Panel H, Columns 16-20: Local Linear Forest for Age
# Simulation: Bivariate score using WGAN data based on Keele & Titiunik (2015)

# We implement local linear forests using grf on the following platform
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

### Get focal point on the treatment boundary where we estimate treatment effect ###
pointsALL <- read.dbf("../../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Data/BorderSegmentPoints_Project.dbf") # boundary data
pointsALL$latitude <- pointsALL$POINT_Y
pointsALL$longitude <- pointsALL$POINT_X
focal <- matrix(c(pointsALL$latitude[45],pointsALL$longitude[45]),nrow = 1)

### Read WGAN-generated population data
data_full <- read.csv("../../../../Data/Keele & Titiunik 2015/WGAN_Epoch1k/age_cWGAN.csv")
focal_df <- read.csv("../../../../Data/Keele & Titiunik 2015/WGAN_Epoch1k/focal_age.csv")

### Population truth
truthT <- mean(focal_df$Y[focal_df$W==1])
truthC <- mean(focal_df$Y[focal_df$W==0])
truth <- truthT - truthC

#### MC Simulation ####
t <- 1000 # Number of Monte Carlo iterations

## Set up parameters for local linear forests
B <- 5000 # number of trees
d <- 2
pi_ <- 1
lower_bound_default <- 1-(1+ (d/(1.3*pi_))*log(0.05)/log(1-0.05))^(-1)  # subsample rate for d=2, alpha=0.05; modified for local linear forests with improved rates
buffer <- 1e-30 # buffer to adjust the focal point to be inside the interior of the support

# Number of observations in each Monte Carlo iteration \in {1000,5000,10000,50000}
for (n in c(1000,5000,10000,50000,1e5)){
  result_llf <- c()
  bias_llf <- c()
  coverage_llf <- c()
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
    
    forestT <- ll_regression_forest(feature_trt, trt[,"Y"], 
                                    num.trees=B, mtry=1,
                                    honesty = TRUE, 
                                    enable.ll.split = TRUE,
                                    ll.split.weight.penalty = TRUE,
                                    sample.fraction = 0.4*ceiling(nrow(trt)^lower_bound_default)/nrow(trt))
    
    predT <- predict(forestT, focal-buffer, 
                     estimate.variance = TRUE,
                     ll.weight.penalty = TRUE
                     )
    predictionT <- as.numeric(predT$predictions)
    BLB1 <- append(BLB1, predT$variance.estimates)
      
    forestC <- ll_regression_forest(feature_ctrl, ctrl[,"Y"], 
                                    num.trees=B, mtry=1,
                                    honesty = TRUE, 
                                    enable.ll.split = TRUE,
                                    ll.split.weight.penalty = TRUE,
                                    sample.fraction = 0.4*ceiling(nrow(ctrl)^lower_bound_default)/nrow(ctrl))
    
    predC <- predict(forestC, focal+buffer, 
                     estimate.variance = TRUE,
                     ll.weight.penalty = TRUE
                     )
    predictionC <- as.numeric(predC$predictions)
    BLB0 <- append(BLB0, predC$variance.estimates)
    
    # get estimate for treatment effect
    prediction <- predictionT - predictionC
    result_llf <- append(result_llf, prediction)
    bias_llf <- append(bias_llf, prediction-truth)
    
    result_trt <- append(result_trt, predictionT)
    result_ctrl <- append(result_ctrl, predictionC)
    
    if ((predictionT-1.96*sqrt(BLB1[i]) <= truthT) & (truthT <= predictionT+1.96*sqrt(BLB1[i]))){
      coverage_trt <- append(coverage_trt, 1)
    }
    else{
      coverage_trt <- append(coverage_trt, 0)
    }
    
    if ((predictionC-1.96*sqrt(BLB0[i]) <= truthC) & (truthC <= predictionC+1.96*sqrt(BLB0[i]))){
      coverage_ctrl <- append(coverage_ctrl, 1)
    }
    else{
      coverage_ctrl <- append(coverage_ctrl, 0)
    }
    
    variance_BLB <- BLB1[i]+BLB0[i]
    se_BLB <- sqrt(variance_BLB)
    
    if ((prediction-1.96*se_BLB <= truth) & (truth <= prediction+1.96*se_BLB)){
      coverage_llf <- append(coverage_llf, 1)
    }
    else{
      coverage_llf <- append(coverage_llf, 0)
    }
    
    try({
      print(paste0("n=",n))
      print("----------------")
      print(paste0("LLF: ", prediction))
      print(paste0("LLF trt: ", predictionT))
      print(paste0("LLF ctrl: ", predictionC))
      print("----------------")
      print(paste0("Mean of LLF:", mean(result_llf, na.rm = TRUE)))
      print(paste0("Bias of LLF:", mean(bias_llf, na.rm = TRUE)))
      print(paste0("Mean of LLF trt:", mean(result_trt, na.rm = TRUE)))
      print(paste0("Mean of LLF ctrl:", mean(result_ctrl, na.rm = TRUE)))
      print("----------------")
      print(paste0("Var of LLF: ", var(result_llf, na.rm = TRUE)))
      print(paste0("Var of LLF trt: ", var(result_trt, na.rm = TRUE)))
      print(paste0("Var of LLF ctrl: ", var(result_ctrl, na.rm = TRUE)))
      print(paste0("Var of LLF (BLB): ", variance_BLB))
      print("----------------")
      print(paste0("CR (", "truth= ", truth, "): ", mean(coverage_llf, na.rm = TRUE)))
      print(paste0("CR trt (", "truthT= ", truthT, "): ", mean(coverage_trt, na.rm = TRUE)))
      print(paste0("CR ctrl (", "truthC= ", truthC, "): ", mean(coverage_ctrl, na.rm = TRUE)))
    })
  }
  assign(paste0("result",n), as.data.frame(cbind("n"=n, "llf"=result_llf, "llf.bias"=bias_llf, "llf.cr"=coverage_llf, "llf.se"=sqrt(BLB0+BLB1))))
}

result <- rbind(result1000,result5000,result10000,result50000,`result1e+05`)
result$truth <- truth
write.csv(result, "../Results csv/wgan-llf-age.csv", row.names = F)