# ==========================
# Toy Simulations
# GBRT simulations without covariates, univariate
# ==========================
# library(grf)
library(tidymodels)
source("Helper_Fns_Sim.R") # This script only requires GBRT() from the helper functions

# DGP
# ==============================
set.seed(1234)
n = 1000 # sample size
t <- 1000 # number of iterations
tau = 0.3 # true treatment effect for the linear DGP
result_univariate <- c() # vector to collect results

for (i in 1:t) {
  # running variable(s): X_i ~ uniform[-1,1]^p
  X1 <- runif(n, -1, 1)
  
  # treatment indicator: D_i = 1 if X_1 >= 0
  D1 <- ifelse(X1 >= 0, 1, 0) # univariate running var
  
  # Generate outcome Y as Bernoulli trials with probability 0.2 + tau*D_i
  Y1 <- rbinom(n, 1, 0.2+tau*D1)
  
  # Dataframe for running vars and response only
  mydata <- as.data.frame(cbind("Y1"= Y1, "X1"=X1, "D1"=D1))
  
  ### Univariate case
  treat <- mydata[mydata$D1==1,]
  control <- mydata[mydata$D1==0,]
  treat <- data.matrix(cbind(treat$X1,treat$Y1))
  colnames(treat) <- c("X","Y")
  control <- data.matrix(cbind(control$X1,control$Y1))
  colnames(control) <- c("X","Y")
  
  boostT_univariate <- GBRT(treat)
  boostC_univariate <- GBRT(control)
  
  cutoff <- matrix(0)
  colnames(cutoff) <- "X"
  mu1_hat <- as.numeric(fit(boostT_univariate, data=treat) %>% predict(new_data=cutoff))
  mu0_hat <- as.numeric(fit(boostC_univariate, data=control) %>% predict(new_data=cutoff))
  
  # get estimated tau at the cutoff
  RF_univariate_pred <- mu1_hat - mu0_hat
  result_univariate <- append(result_univariate, RF_univariate_pred)
  
  try({print(i)
    print(paste0("Univariate: ", RF_univariate_pred))
    print(paste0("Mean of univariate:", mean(result_univariate)))
    print(paste0("Var of univariate:", var(result_univariate)))
    print(paste0("MSE of univariate: ", (mean(result_univariate)-tau)^2+var(result_univariate)))
    print("----------------")
  })
}