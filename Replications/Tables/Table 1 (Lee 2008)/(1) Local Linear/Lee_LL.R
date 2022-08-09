######################################
# Simulation: Univariate score using
# data from Lee (2008)
######################################
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

t <- 1000 
n <- 5000

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




