######################################
# Table 1, Panel A, Columns 1-5: Local Linear Regression
# Simulation: Univariate score using data from polynomial DGP

# We implement local linear regressions using rdrobust on the following platform
# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux trixie/sid
######################################
rm(list=ls())
library(rdrobust)
set.seed(1234)
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

t <- 1000 # Number of Monte Carlo iterations

# Number of observations in each Monte Carlo iteration \in {1000,5000,10000,50000}
for (n in c(1000,5000,10000,50000,1e5)){
  result_rdrobust <- c()
  bias_rdrobust <- c()
  coverage_rdrobust <- c()
  
  for (i in 1:t) {
    print(i)
    
    X <- 2*rbeta(n,2,4)-1
    epsilon <- rnorm(n, 0, sigma_y)
    Y <- ifelse(X>=0, c_r0+c_r1*X+c_r2*X^2+c_r3*X^3+c_r4*X^4+c_r5*X^5+epsilon, 
                c_l0+c_l1*X+c_l2*X^2+c_l3*X^3+c_l4*X^4+c_l5*X^5+epsilon)
    
    data <- as.data.frame(cbind("Y" = Y, "X"=X))
    
    
    ############################
    # rdrobust
    ############################
    RBC <- rdrobust(y=data$Y, x=data$X, all=TRUE)
    result_rdrobust <- append(result_rdrobust, RBC$Estimate[2])
    bias_rdrobust <- append(bias_rdrobust, RBC$Estimate[2]-truth)
    
    if ((RBC$ci[3,1] <= truth) & (truth <= RBC$ci[3,2])) {
      coverage_rdrobust <- append(coverage_rdrobust, 1)
    }
    else {
      coverage_rdrobust <- append(coverage_rdrobust, 0)
    }
    
    try({
      print(paste0("n=",n))
      print(paste0("Prediction: ", RBC$Estimate[2]))
      print(paste0("Mean:", mean(result_rdrobust, na.rm = TRUE)))
      print(paste0("Bias:", mean(bias_rdrobust, na.rm = TRUE)))
      print(paste0("Var: ", var(result_rdrobust, na.rm = TRUE)))
      print(paste0("MSE: ", (mean(result_rdrobust,na.rm = TRUE)-truth)^2+var(result_rdrobust,na.rm = TRUE)))
      print(paste0("CR: ", mean(coverage_rdrobust, na.rm = TRUE)))
    })
    
    assign(paste0("result",n), as.data.frame(cbind("n"=n, "rdrobust"=result_rdrobust, "rdrobust.bias"=bias_rdrobust, "rdrobust.cr"=coverage_rdrobust)))
  }
}

result <- rbind(result1000,result5000,result10000,result50000,`result1e+05`)
result$truth <- truth
write.csv(result, "../Results csv/poly-lee-rdrobust.csv", row.names = F)

