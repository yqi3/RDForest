######################################
# Table 1, Panel B, Columns 6-10: Minimax Optimal
# Simulation: Univariate score using WGAN data based on Lee (2008)

# We implement minimax-optimal estimator using optrdd on the following platform
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.5
######################################
rm(list=ls())
set.seed(1234)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(optrdd)

data_full <- read.csv("../../../../Data/Lee2008/generated_Lee.csv")
focal_df <- read.csv("../../../../Data/Lee2008/focal_lee.csv")
truthT <- mean(focal_df$Y[focal_df$W==1])
truthC <- mean(focal_df$Y[focal_df$W==0])
truth <- truthT - truthC
focal <- 0

#### MC Simulation ####
t <- 1000 # number of Monte Carlo iterations

# Number of observations in each Monte Carlo iteration \in {1000,5000,10000,50000}
for (n in c(1000,5000,10000,50000,1e5)){
  i <- 1
  result_optrdd <- c()
  bias_optrdd <- c()
  coverage_optrdd <- c()

  while (i <= t){
    print(i)
    # generate simulation data
    rand <- sample.int(nrow(data_full), n, replace = F)
    datause <- data_full[rand,]
    X <- datause$X
    Y <- datause$Y
    W <- datause$W

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
    result_optrdd <- append(result_optrdd, prediction)
    bias_optrdd <- append(bias_optrdd, prediction-truth)})
  
    if (last == length(result_optrdd)) {
      next
    }
  
    # Coverage
    if ((prediction-out.pt$tau.plusminus <= truth) & (truth <= prediction+out.pt$tau.plusminus)){
      coverage_optrdd <- append(coverage_optrdd, 1)
    } else{
      coverage_optrdd <- append(coverage_optrdd, 0)
    }
    
    try({
      print(paste0("n=",n))
      print(paste0("Prediction: ", prediction))
      print(paste0("Mean:", mean(result_optrdd, na.rm = TRUE)))
      print(paste0("Bias:", mean(bias_optrdd, na.rm = TRUE)))
      print(paste0("Var: ", var(result_optrdd, na.rm = TRUE)))
      print(paste0("MSE: ", (mean(result_optrdd, na.rm = TRUE)-truth)^2+var(result_optrdd, na.rm = TRUE)))
      print(paste0("CR: ", mean(coverage_optrdd, na.rm = TRUE)))
    })
    
    assign(paste0("result",n), as.data.frame(cbind("n"=n, "optrdd"=result_optrdd, "optrdd.bias"=bias_optrdd, "optrdd.cr"=coverage_optrdd)))
    
    i <- i+1
  }
}

result <- rbind(result1000,result5000,result10000,result50000,`result1e+05`)
result$truth <- truth
write.csv(result, "../Results csv/wgan-lee-optrdd.csv", row.names = F)
