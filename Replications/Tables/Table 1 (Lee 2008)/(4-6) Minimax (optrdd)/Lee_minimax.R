######################################
# Table 1: (4-6) Minimax-Optimal
# Simulation: Univariate score using data from Lee (2008)

# We implement minimax-optimal estimator using optrdd on the following platform
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.5
######################################
rm(list=ls())
set.seed(1234)

library(optrdd)

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
focal <- 0

#### MC Simulation ####
t <- 1000 # Number of Monte Carlo iterations

# Number of observations in each Monte Carlo iteration
# Comment/uncomment the three lines below to reproduce columns (4) - (6) of results in the paper
n <- 1000  # (4)
# n <- 5000  # (5)
# n <- 10000  # (6)

i <- 1
result_optrdd <- c()
coverage_optrdd <- c()

while (i <= t){
  print(i)
  # generate simulation data
  X <- 2*rbeta(n,2,4)-1
  epsilon <- rnorm(n, 0, sigma_y)
  Y <- ifelse(X>=0, c_r0+c_r1*X+c_r2*X^2+c_r3*X^3+c_r4*X^4+c_r5*X^5+epsilon,
              c_l0+c_l1*X+c_l2*X^2+c_l3*X^3+c_l4*X^4+c_l5*X^5+epsilon)
  W <- as.numeric(X>=0)
  
  # Estimate 2nd derivative bound using 2nd order polynomial
  X2 <-  X^2
  RDF <- data.frame(W=W, X=X, X2=X2, Y=Y)
  m <- lm(Y ~ W * (X + X2), data = RDF)
  coef <- m$coefficients[c(4,6)]
  b <- max(4*coef)  # can modify this multiplier, results in the paper are generated using 4

  # Run optrdd
  last <- length(result_optrdd)
  try({out.pt = optrdd(X=X, Y=Y, W=W, max.second.derivative = b, estimation.point = focal, verbose = FALSE)
  prediction <- out.pt$tau.hat
  result_optrdd <- append(result_optrdd, prediction)})

  if (last == length(result_optrdd)) {
    next
  }

  # Coverage
  if ((prediction-out.pt$tau.plusminus <= truth) & (truth <= prediction+out.pt$tau.plusminus)){
    coverage_optrdd <- append(coverage_optrdd, 1)
  } else{
    coverage_optrdd <- append(coverage_optrdd, 0)
  }
  
  try({print("----------------")
    print(paste0("Prediction: ", prediction))
    print("----------------")
    print(paste0("Mean:", mean(result_optrdd, na.rm = TRUE)))
    print("----------------")
    print(paste0("Var: ", var(result_optrdd, na.rm = TRUE)))
    print("----------------")
    print(paste0("MSE: ", (mean(result_optrdd, na.rm = TRUE)-truth)^2+var(result_optrdd, na.rm = TRUE)))
    print("----------------")
    print(paste0("CR: ", mean(coverage_optrdd, na.rm = TRUE)))
    print("----------------")
  })

  i <- i+1
}
