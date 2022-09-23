######################################
# Table A.1: (4) GBRT
# Panel B: Bivariate running variable
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
result_bivariate <- c()

for (i in 1:t) {
  # running variable(s): X_i ~ uniform[-1,1]^p
  X1 <- runif(n, -1, 1)
  X2 <- runif(n, -1, 1)
  
  # treatment indicator: D_i = 1 if X_1 >= 0
  D2 <- ifelse(X1 >= 0 & X2 >= 0, 1, 0) # multivariate running var (2)
  
  # Generate outcome Y as Bernoulli trials with probability 0.2 + tau*D_i
  Y2 <- rbinom(n, 1, 0.2+tau*D2)
  
  # Dataframe for running vars and response only
  mydata <- as.data.frame(cbind("Y2" = Y2, "X1"=X1, "X2"=X2, "D2"=D2))
  
  ### Bivariate case
  treat <- mydata[mydata$D2==1,]
  control <- mydata[mydata$D2==0,]
  treat <- data.matrix(cbind(treat$X1,treat$X2,treat$Y2))
  colnames(treat) <- c("X1","X2","Y")
  control <- data.matrix(cbind(control$X1,control$X2,control$Y2))
  colnames(control) <- c("X1","X2","Y")
  
  boostT_bivariate <- GBRT(treat)
  boostC_bivariate <- GBRT(control)
  
  cutoff <- matrix(c(0,0),ncol=2)
  colnames(cutoff) <- c("X1","X2")
  mu1_hat <- as.numeric(fit(boostT_bivariate, data=treat) %>% predict(new_data=cutoff))
  mu0_hat <- as.numeric(fit(boostC_bivariate, data=control) %>% predict(new_data=cutoff))
  
  # get estimated tau at the cutoff
  RF_bivariate_pred <- mu1_hat - mu0_hat
  result_bivariate <- append(result_bivariate, RF_bivariate_pred)
  
  try({print(i)
    print(paste0("Bivariate: ", RF_bivariate_pred))
    print(paste0("Mean of bivariate:", mean(result_bivariate)))
    print(paste0("Var of bivariate:", var(result_bivariate)))
    print(paste0("MSE of bivariate: ", (mean(result_bivariate)-tau)^2+var(result_bivariate)))
    print("----------------")
  })
}