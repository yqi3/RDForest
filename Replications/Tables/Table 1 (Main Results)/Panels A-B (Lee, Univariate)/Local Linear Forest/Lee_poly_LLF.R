######################################
# Table 1, Panel A, Columns 16-20: Local Linear Forest
# Simulation: Univariate score using data from polynomial DGP

# We implement local linear forests using grf on the following platform
# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux trixie/sid
######################################
rm(list=ls())
set.seed(1234)

library(grf)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

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
focal <- matrix(0)

#### MC Simulation ####
t <- 1000 # number of Monte Carlo iterations

## Set up parameters for local linear forests
B <- 5000 # number of trees
d <- 1
pi_ <- 1
lower_bound_default <- 1-(1+ (d/(1.3*pi_))*log(0.05)/log(1-0.05))^(-1)  # subsample rate for d=1, alpha=0.05; modified for local linear forests with improved rates
buffer <- 1e-30 # buffer to adjust the focal point to be inside the interior of the support

# Number of observations in each Monte Carlo iteration \in {1000,5000,10000,50000}
for (n in c(1000,5000,10000,50000,1e5)){
  result_llf <- c()
  bias_llf <- c()
  coverage_llf <- c()
  BLB1 <- c()
  BLB0 <- c()
  
  for (i in 1:t) {
    print(i)
    # generate simulation data
    X <- 2*rbeta(n,2,4)-1
    epsilon <- rnorm(n, 0, sigma_y)
    Y <- ifelse(X>=0, c_r0+c_r1*X+c_r2*X^2+c_r3*X^3+c_r4*X^4+c_r5*X^5+epsilon, 
                c_l0+c_l1*X+c_l2*X^2+c_l3*X^3+c_l4*X^4+c_l5*X^5+epsilon)
    
    data <- as.data.frame(cbind("Y" = Y, "X"=X))
    
    ### llf ###
    trt <- data[data$X>=0,]
    ctrl <- data[data$X<0,]
    feature_trt <- as.data.frame(trt[,2])
    feature_ctrl <- as.data.frame(ctrl[,2])
    
    # -----------------------------------------------------------------------------------
    # Modified tuning procedure based on the source code of grf; see
    # https://github.com/grf-labs/grf/blob/master/r-package/grf/R/tune_forest.R
    # -----------------------------------------------------------------------------------
    
    forestT <- ll_regression_forest(feature_trt, trt$Y, 
                                    num.trees=B, mtry=1,
                                    honesty = TRUE, 
                                    enable.ll.split = TRUE,
                                    ll.split.weight.penalty = TRUE,
                                    sample.fraction = 0.4*ceiling(nrow(trt)^lower_bound_default)/nrow(trt))
    
    predT <- predict(forestT, focal+buffer, 
                     estimate.variance = TRUE,
                     ll.weight.penalty = TRUE
    )
    predictionT <- as.numeric(predT$predictions)
    BLB1 <- append(BLB1, predT$variance.estimates)
    
    forestC <- ll_regression_forest(feature_ctrl, ctrl$Y, 
                                    num.trees=B, mtry=1,
                                    honesty = TRUE, 
                                    enable.ll.split = TRUE,
                                    ll.split.weight.penalty = TRUE,
                                    sample.fraction = 0.4*ceiling(nrow(ctrl)^lower_bound_default)/nrow(ctrl))
    
    predC <- predict(forestC, focal-buffer, 
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
write.csv(result, "../Results csv/poly-llf-lee.csv", row.names = F)
