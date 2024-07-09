######################################
# Table 1, Panel D, Columns 11-15: Honest Regression Forest for Housing Price
# Simulation: Bivariate score using data from polynomial DGP

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

### Read sample data ###
library(haven)
data <- read_dta("../../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Data/Housing/NJ_House_Final.dta")
data <- data[complete.cases(data$latitude, data$longitude, data$treat, data[,c("sqft_price")]),]

### Get focal point on the treatment boundary where we estimate treatment effect ###
pointsALL <- read.dbf("../../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Data/BorderSegmentPoints_Project.dbf") # boundary data
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

truthT <- predict(cef, data.frame(X1=X1, X1.2=(X1 - mean_X1)^2, X1.3=(X1 - mean_X1)^3,
                                  X2=X2, X2.2=(X2 - mean_X2)^2, X2.3=(X2 - mean_X2)^3,
                                  X12=(X1 - mean_X1)*(X2 - mean_X2), 
                                  W=1))

truthC <- predict(cef, data.frame(X1=X1, X1.2=(X1 - mean_X1)^2, X1.3=(X1 - mean_X1)^3,
                                  X2=X2, X2.2=(X2 - mean_X2)^2, X2.3=(X2 - mean_X2)^3,
                                  X12=(X1 - mean_X1)*(X2 - mean_X2), 
                                  W=0))

truth <- truthT-truthC

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
  
  for (i in 1:t){
    print(i)
    # generate simulation data
    datause <- genDGP(n, "price")
    
    ### RF ###
    trt <- datause[datause$treat==1,]
    ctrl <- datause[datause$treat==0,]
    feature_trt <- as.data.frame(trt[,c(1,2)])
    feature_ctrl <- as.data.frame(ctrl[,c(1,2)])
    
    # -----------------------------------------------------------------------------------
    # Modified tuning procedure based on the source code of grf; see
    # https://github.com/grf-labs/grf/blob/master/r-package/grf/R/tune_forest.R
    # -----------------------------------------------------------------------------------
      
    forestT <- tuned_RDForest(feature_trt, trt[,4], tune.parameters = c(), type = "trt")
    predictionT <- predict(forestT, focal, estimate.variance = TRUE)
    BLB1 <- append(BLB1, predictionT$variance.estimates)
    forestC <- tuned_RDForest(feature_ctrl, ctrl[,4], tune.parameters = c(), type = "ctrl")
    predictionC <- predict(forestC, focal, estimate.variance = TRUE)
    BLB0 <- append(BLB0, predictionC$variance.estimates)
    
    result_trt <- append(result_trt, predictionT$predictions)
    result_ctrl <- append(result_ctrl, predictionC$predictions)
    
    # get estimate for treatment effect
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
    }
    else{
      coverage_rf <- append(coverage_rf, 0)
    }
    
    try({
      print(paste0("n=", n))
      print("----------------")
      print(paste0("RF: ", prediction))
      print(paste0("Mean of RF:", mean(result_rf, na.rm = TRUE)))
      print(paste0("Bias of RF:", mean(bias_rf, na.rm = TRUE)))
      print("----------------")
      print(paste0("Var of RF: ", var(result_rf, na.rm = TRUE)))
      print(paste0("Var of RF (BLB): ", variance_BLB))
      print(paste0("MSE of RF: ", (mean(result_rf,na.rm = TRUE)-truth)^2+var(result_rf,na.rm = TRUE)))
      print("----------------")
      print(paste0("CR of RF: ", mean(coverage_rf, na.rm = TRUE)))
      print(paste0("CR trt (truthT=", truthT, "): ", mean(coverage_trt, na.rm = TRUE)))
      print(paste0("CR ctrl (truthT=", truthC, "): ", mean(coverage_ctrl, na.rm = TRUE)))
      print("----------------")
    })
  }
  
  assign(paste0("result",n), as.data.frame(cbind("n"=n, "rf"=result_rf, "rf.bias"=bias_rf, "rf.cr"=coverage_rf, "rf.se"=sqrt(BLB0+BLB1))))
}

result <- rbind(result1000,result5000,result10000,result50000,`result1e+05`)
result$truth <- truth
write.csv(result, paste0("../Results csv/poly-rf-price.csv"), row.names = F)

