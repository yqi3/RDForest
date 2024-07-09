######################################
# Table 1, Panel B, Columns 1-5: Local Linear Regression
# Simulation: Univariate score using WGAN data based on Lee (2008)

# We implement local linear regressions using rdrobust on the following platform
# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux trixie/sid
######################################
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(rdrobust)
set.seed(1234)

data_full <- read.csv("../../../../Data/Lee2008/generated_Lee.csv")
focal_df <- read.csv("../../../../Data/Lee2008/focal_lee.csv")
truthT <- mean(focal_df$Y[focal_df$W==1])
truthC <- mean(focal_df$Y[focal_df$W==0])
truth <- truthT - truthC
focal <- 0

t <- 1000 # Number of Monte Carlo iterations

# Number of observations in each Monte Carlo iteration in {1000,5000,10000,50000}
for (n in c(1000,5000,10000,50000,1e5)){
  result_rdrobust <- c()
  bias_rdrobust <- c()
  coverage_rdrobust <- c()

  for (i in 1:t) {
    print(i)
    # generate simulation data
    rand <- sample.int(nrow(data_full), n, replace = F)
    data <- data_full[rand,]

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
      print(paste0("rdrobust: ", RBC$Estimate[2]))
      print(paste0("Mean of rdrobust:", mean(result_rdrobust, na.rm = TRUE)))
      print(paste0("Bias of rdrobust:", mean(bias_rdrobust, na.rm = TRUE)))
      print(paste0("Var of rdrobust: ", var(result_rdrobust, na.rm = TRUE)))
      print(paste0("MSE of rdrobust: ", (mean(result_rdrobust,na.rm = TRUE)-truth)^2+var(result_rdrobust,na.rm = TRUE)))
      print(paste0("CR of rdrobust: ", mean(coverage_rdrobust, na.rm = TRUE)))
    })
    
    assign(paste0("result",n), as.data.frame(cbind("n"=n, "rdrobust"=result_rdrobust, "rdrobust.bias"=bias_rdrobust, "rdrobust.cr"=coverage_rdrobust)))
  }
}

result <- rbind(result1000,result5000,result10000,result50000,`result1e+05`)
result$truth <- truth
write.csv(result, "../Results csv/wgan-lee-rdrobust.csv", row.names = F)

