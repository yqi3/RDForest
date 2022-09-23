######################################
# Table A.1: (4) GBRT
# Panel C: Trivariate running variable
# Toy Simulations

# We implement GBRT using tidymodels on the following platform
# R version 4.1.1 (2021-08-10)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: CentOS release 6.5 (Final)
######################################
rm(list=ls())

library(tidymodels)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Toy_Sim_Helper_Fns.R") # This script only requires GBRT() from the helper functions

# DGP
# ==============================
set.seed(1234)
n = 1000 # sample size
t <- 1000 # number of iterations
tau <- 0.3 # true treatment effect for the linear DGP
result_trivariate <- c()

for (i in 1:t) {
  # running variable(s): X_i ~ uniform[-1,1]^p
  X1 <- runif(n, -1, 1)
  X2 <- runif(n, -1, 1)
  X3 <- runif(n, -1, 1)
  
  # treatment indicator: D_i = 1 if X_1 >= 0
  D3 <- ifelse(X1 >= 0 & X2 >= 0 & X3 >= 0, 1, 0) # multivariate running var (3)
  
  # Generate outcome Y as Bernoulli trials with probability 0.2 + tau*D_i
  Y3 <- rbinom(n, 1, 0.2+tau*D3)
  
  # Dataframe for running vars and response only
  mydata <- as.data.frame(cbind("Y3" = Y3, "X1"=X1, "X2"=X2, "X3"=X3, "D3"=D3))
  
  ### Trivariate case
  treat <- mydata[mydata$D3==1,]
  control <- mydata[mydata$D3==0,]
  treat <- data.matrix(cbind(treat$X1,treat$X2,treat$X3,treat$Y3))
  colnames(treat) <- c("X1","X2","X3","Y")
  control <- data.matrix(cbind(control$X1,control$X2,control$X3,control$Y3))
  colnames(control) <- c("X1","X2","X3","Y")
  
  boostT_trivariate <- GBRT(treat)
  boostC_trivariate <- GBRT(control)
  
  cutoff <- matrix(c(0,0,0),ncol=3)
  colnames(cutoff) <- c("X1","X2","X3")
  mu1_hat <- as.numeric(fit(boostT_trivariate, data=treat) %>% predict(new_data=cutoff))
  mu0_hat <- as.numeric(fit(boostC_trivariate, data=control) %>% predict(new_data=cutoff))
  
  # get estimated tau at the cutoff
  RF_trivariate_pred <- mu1_hat - mu0_hat
  result_trivariate <- append(result_trivariate, RF_trivariate_pred)
  
  try({print(i)
    print(paste0("Trivariate: ", RF_trivariate_pred))
    print(paste0("Mean of trivariate:", mean(result_trivariate)))
    print(paste0("Var of trivariate:", var(result_trivariate)))
    print(paste0("MSE of trivariate: ", (mean(result_trivariate)-tau)^2+var(result_trivariate)))
    print("----------------")
  })
}