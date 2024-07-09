######################################
# Table 1, Panel D, Columns 16-20: Local Linear Forest for Housing Price
# Simulation: Bivariate score using data from polynomial DGP

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
source("KeeleTitiunik2014-PA-replication-files/Analysis/Local Poly/distance-functions.R")
source("KT_Helper_Fns.R")

### Read sample data ###
library(haven)
data <- read_dta("KeeleTitiunik2014-PA-replication-files/Data/Housing/NJ_House_Final.dta")
data <- data[complete.cases(data$latitude, data$longitude, data$treat, data[,c("sqft_price")]),]

### Get focal point on the treatment boundary where we estimate treatment effect ###
pointsALL <- read.dbf("KeeleTitiunik2014-PA-replication-files/Data/BorderSegmentPoints_Project.dbf") # boundary data
pointsALL$latitude= pointsALL$POINT_Y
pointsALL$longitude= pointsALL$POINT_X
focal <- matrix(c(pointsALL$latitude[45],pointsALL$longitude[45]),nrow = 1)

X1 <- data$latitude
X2 <- data$longitude
mean_X1 <- mean(X1)
mean_X2 <- mean(X2)
price <- data$sqft_price

model.df <- data.frame(Y=price, X1=X1, X1.2=(X1 - mean_X1)^2, X1.3=(X1 - mean_X1)^3,
                       X2=X2, X2.2=(X2 - mean_X2)^2, X2.3=(X2 - mean_X2)^3,
                       X12=(X1 - mean_X1)*(X2 - mean_X2),
                       W=data$treat)

cef <- lm(Y ~ W * ., data = model.df)
sigma_y <- sqrt(summary(cef)$sigma^2)

X1 <- pointsALL$latitude[45]
X2 <- pointsALL$longitude[45]
truth <- predict(cef, data.frame(X1=X1, X1.2=(X1 - mean_X1)^2, X1.3=(X1 - mean_X1)^3,
                                 X2=X2, X2.2=(X2 - mean_X2)^2, X2.3=(X2 - mean_X2)^3,
                                 X12=(X1 - mean_X1)*(X2 - mean_X2), 
                                 W=1)) - predict(cef, data.frame(X1=X1, X1.2=(X1 - mean_X1)^2, X1.3=(X1 - mean_X1)^3,
                                                                 X2=X2, X2.2=(X2 - mean_X2)^2, X2.3=(X2 - mean_X2)^3,
                                                                 X12=(X1 - mean_X1)*(X2 - mean_X2), 
                                                                 W=0))

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
  
  for (i in 1:t){
    print(i)
    # generate simulation data
    datause <- genDGP(n, "price")
    
    ### RF ###
    trt <- datause[datause$treat==1,]
    ctrl <- datause[datause$treat==0,]
    feature_trt <- as.data.frame(trt[,c(1,2)])
    feature_ctrl <- as.data.frame(ctrl[,c(1,2)])
    
    forestT <- ll_regression_forest(feature_trt, trt[,4], 
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
    
    forestC <- ll_regression_forest(feature_ctrl, ctrl[,4], 
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
      print(paste0("Mean of LLF:", mean(result_llf, na.rm = TRUE)))
      print(paste0("Bias of LLF:", mean(bias_llf, na.rm = TRUE)))
      print("----------------")
      print(paste0("Var of LLF: ", var(result_llf, na.rm = TRUE)))
      print(paste0("Var of LLF (BLB): ", variance_BLB))
      print("----------------")
      print(paste0("CR (", "truth= ", truth, "): ", mean(coverage_llf, na.rm = TRUE)))
    })
  }
  assign(paste0("result",n), as.data.frame(cbind("n"=n, "llf"=result_llf, "llf.bias"=bias_llf, "llf.cr"=coverage_llf, "llf.se"=sqrt(BLB0+BLB1))))
}

result <- rbind(result1000,result5000,result10000,result50000,`result1e+05`)
result$truth <- truth
write.csv(result, "../Results csv/poly-llf-price.csv", row.names = F)
