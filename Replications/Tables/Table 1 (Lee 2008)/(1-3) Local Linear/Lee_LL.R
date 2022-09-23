######################################
# Table 1: (1-3) Local Linear
# Simulation: Univariate score using data from Lee (2008)

# We implement local linear regressions using rdrobust on the following platform
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.5
######################################
rm(list=ls())

library(rdrobust)
set.seed(1234)

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

# Number of observations in each Monte Carlo iteration
# Comment/uncomment the three lines below to reproduce columns (1) - (3) of results in the paper
n <- 1000  # (1)
# n <- 5000  # (2)
# n <- 10000  # (3)

result_rdrobust <- c()
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
  
  if ((RBC$ci[3,1] <= truth) & (truth <= RBC$ci[3,2])) {
    coverage_rdrobust <- append(coverage_rdrobust, 1)
  }
  else {
    coverage_rdrobust <- append(coverage_rdrobust, 0)
  }
  
  try({
  print(paste0("rdrobust: ", RBC$Estimate[2]))
  print("----------------")
  print(paste0("Mean of rdrobust:", mean(result_rdrobust, na.rm = TRUE)))
  print("----------------")
  print(paste0("Var of rdrobust: ", var(result_rdrobust, na.rm = TRUE)))
  print("----------------")
  print(paste0("MSE of rdrobust: ", (mean(result_rdrobust,na.rm = TRUE)-truth)^2+var(result_rdrobust,na.rm = TRUE)))
  print("----------------")
  print(paste0("CR of rdrobust: ", mean(coverage_rdrobust, na.rm = TRUE)))
  print("----------------")
  })
}




