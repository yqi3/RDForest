# ==========================
# Toy Simulations
# Rdrobust simulations without covariates, 1-3 running vars
# ==========================
library(rdrobust)

# DGP
# ==============================
set.seed(1234)
n = 1000 # sample size
t <- 1000 # number of iterations
tau = 0.3 # true treatment effect for the linear DGP
result_univariate <- c() # vector to collect results
result_bivariate <- c()
result_trivariate <- c()
coverage_univariate <- c()
coverage_bivariate <- c()
coverage_trivariate <- c()
i <- 1

while (i <= t) {
  # running variable(s): X_i ~ uniform[-1,1]^p
  X1 <- runif(n, -1, 1)
  X2 <- runif(n, -1, 1)
  X3 <- runif(n, -1, 1)
  
  # treatment indicator: D_i = 1 if X_1 >= 0
  D1 <- ifelse(X1 >= 0, 1, 0) # univariate running var
  D2 <- ifelse(X1 >= 0 & X2 >= 0, 1, 0) # multivariate running var (2)
  D3 <- ifelse(X1 >= 0 & X2 >= 0 & X3 >= 0, 1, 0) # multivariate running var (3)
  
  # Generate outcome Y as Bernoulli trials with probability 0.2 + tau*D_i
  Y1 <- rbinom(n, 1, 0.2+tau*D1)
  Y2 <- rbinom(n, 1, 0.2+tau*D2)
  Y3 <- rbinom(n, 1, 0.2+tau*D3)
  
  # Dataframe for running vars and response only
  mydata <- as.data.frame(cbind("Y1"= Y1, "Y2" = Y2, "Y3" = Y3, 
                                "X1"=X1, "X2"=X2, "X3"=X3, "D1"=D1, "D2"=D2, "D3"=D3))
  
  ### Start with the trivariate case, in case we encounter an error and need to rerun the iteration
  last <- length(result_trivariate)
  poolrv <- ifelse(D3==1, sqrt(X1^2+X2^2+X3^2), -sqrt(X1^2+X2^2+X3^2))
  try({RBC_trivariate <- rdrobust(y=mydata$Y3, x=poolrv, kernel="tri", all=TRUE)
  RBC_trivariate_pred <- RBC_trivariate$Estimate[2]
  result_trivariate <- append(result_trivariate, RBC_trivariate_pred)}, silent=TRUE)
  if (last == length(result_trivariate)) {
    next
  }
  
  # Coverage
  if ((RBC_trivariate$ci[3,1] <= tau) & (tau <= RBC_trivariate$ci[3,2])){
    coverage_trivariate <- append(coverage_trivariate, 1)
  }
  else{
    coverage_trivariate <- append(coverage_trivariate, 0)
  }
  
  ### Univariate case
  RBC_univariate <- rdrobust(y=mydata$Y1, x=mydata$X1, kernel="tri", all=TRUE)
  RBC_univariate_pred <- RBC_univariate$Estimate[2]
  result_univariate <- append(result_univariate, RBC_univariate_pred)
  
  # Coverage
  if ((RBC_univariate$ci[3,1] <= tau) & (tau <= RBC_univariate$ci[3,2])){
    coverage_univariate <- append(coverage_univariate, 1)
  }
  else{
    coverage_univariate <- append(coverage_univariate, 0)
  }
  
  ### Bivariate case
  poolrv <- ifelse(D2==1, sqrt(X1^2+X2^2), -sqrt(X1^2+X2^2))
  RBC_bivariate <- rdrobust(y=mydata$Y2, x=poolrv, kernel="tri", all=TRUE)
  RBC_bivariate_pred <- RBC_bivariate$Estimate[2]
  result_bivariate <- append(result_bivariate, RBC_bivariate_pred)
  
  # Coverage
  if ((RBC_bivariate$ci[3,1] <= tau) & (tau <= RBC_bivariate$ci[3,2])){
    coverage_bivariate <- append(coverage_bivariate, 1)
  }
  else{
    coverage_bivariate <- append(coverage_bivariate, 0)
  }
  
  try({print(i)
    print(paste0("Univariate: ", RBC_univariate_pred))
    print(paste0("Mean of univariate:", mean(result_univariate)))
    print(paste0("Var of univariate:", var(result_univariate)))
    print(paste0("MSE of univariate: ", (mean(result_univariate)-tau)^2+var(result_univariate)))
    print(paste0("CR of univariate:", mean(coverage_univariate)))
    print("----------------")
    print(paste0("Bivariate: ", RBC_bivariate_pred))
    print(paste0("Mean of bivariate:", mean(result_bivariate)))
    print(paste0("Var of bivariate:", var(result_bivariate)))
    print(paste0("MSE of bivariate: ", (mean(result_bivariate)-tau)^2+var(result_bivariate)))
    print(paste0("CR of bivariate:", mean(coverage_bivariate)))
    print("----------------")
    print(paste0("Trivariate: ", RBC_trivariate_pred))
    print(paste0("Mean of trivariate:", mean(result_trivariate)))
    print(paste0("Var of trivariate:", var(result_trivariate)))
    print(paste0("MSE of trivariate: ", (mean(result_trivariate)-tau)^2+var(result_trivariate)))
    print(paste0("CR of trivariate:", mean(coverage_trivariate)))
    print("----------------")
  })
  
  i <- i+1
}


