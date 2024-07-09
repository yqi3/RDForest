######################################
# Table 1, Panel D, Columns 6-10: Minimax Optimal Estimator for Housing Price
# Simulation: Bivariate score using data from polynomial DGP

# We implement minimax-optimal estimator using optrdd on the following platform
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.5
######################################
rm(list=ls())
set.seed(1234)

library(foreign)
library(optrdd)
library(glmnet)
library(splines)
library(tidymodels)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Analysis/Local Poly/distance-functions.R")
source("../../../../Helper Fns/KT_Helper_Fns.R")

### Curvature functions ###
compute_curvature = function(xgrid.pred, nA, nB, centers, binw) {
  curvs0 = sapply(centers, function(idx) {
    c((xgrid.pred[idx + 1] + xgrid.pred[idx - 1] - 2 * xgrid.pred[idx]) / binw^2,
      (xgrid.pred[idx + nA] + xgrid.pred[idx - nA] - 2 * xgrid.pred[idx]) / binw^2,
      (xgrid.pred[idx + nA + 1] + xgrid.pred[idx - nA - 1] - 2 * xgrid.pred[idx]) / binw^2 / 2,
      (xgrid.pred[idx + nA - 1] + xgrid.pred[idx - nA + 1] - 2 * xgrid.pred[idx]) / binw^2 / 2)
  })
  curvs = apply(curvs0, 2, function(iii) max(abs(iii)))
  quantile(curvs, 0.95, na.rm=TRUE)
}

get_curvature = function(xx, outcome, ww, binw = 0.1) {
  
  xx = data.frame(xx)
  yy = datause[,outcome]
  
  names(xx) = c("A", "B")
  
  gridA = seq(min(xx$A) - 1.5 * binw, max(xx$A) + 1.5 * binw, by = binw)
  gridB = seq(min(xx$B) - 1.5 * binw, max(xx$B) + 1.5 * binw, by = binw)
  xgrid = data.frame(expand.grid(A=gridA, B = gridB))
  xspl.all = model.matrix(~ 0 + ns(A, df = 7) * ns(B, df = 7),
                          data = rbind(xx, xgrid))
  xspl = xspl.all[1:nrow(xx),]
  
  fit0 = cv.glmnet(xspl[ww==0,], yy[ww==0], alpha = 0)
  fit1 = cv.glmnet(xspl[ww==1,], yy[ww==1], alpha = 0)
  
  xgrid.spl = xspl.all[nrow(xx) + 1:nrow(xgrid),]
  # In the original code of Imbens and Wager (2019), the option s="lambda.1se" is used in the predict() function. We changed this to s="lambda.min" instead. See: https://stackoverflow.com/a/40740273.
  xgrid.pred.0 = predict(fit0, xgrid.spl, s="lambda.min")
  xgrid.pred.1 = predict(fit1, xgrid.spl, s="lambda.min")
  
  nA = length(gridA)
  nB = length(gridB)
  
  bucketA = as.numeric(cut(xx[,1], gridA))
  bucketB = as.numeric(cut(xx[,2], gridB))
  bucket = bucketA + nA * bucketB
  
  c.hat.0 = compute_curvature(xgrid.pred.0, nA, nB, centers = bucket[ww == 0], binw = binw)
  c.hat.1 = compute_curvature(xgrid.pred.1, nA, nB, centers = bucket[ww == 1], binw = binw)
  max(c.hat.0, c.hat.1)
}

### Read sample data ###
library(haven)
data <- read_dta("../../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Data/Housing/NJ_House_Final.dta")
data <- data[complete.cases(data$latitude, data$longitude, data$treat, data[,c("sqft_price")]),]

### Get focal point on the treatment boundary where we estimate treatment effect ###
pointsALL <- read.dbf("../../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Data/BorderSegmentPoints_Project.dbf") # boundary data
pointsALL$latitude= pointsALL$POINT_Y
pointsALL$longitude= pointsALL$POINT_X
focal <- matrix(c(pointsALL$latitude[45],pointsALL$longitude[45]),nrow = 1)
BOUNDARY_PT <- 45

X1 <- data$latitude
X2 <- data$longitude
mean_X1 <- mean(X1)
mean_X2 <- mean(X2)
price <- data$sqft_price

model.df <- data.frame(Y=price, X1=X1, X1.2=(X1 - mean_X1)^2, X1.3=(X1 - mean_X1)^3,
                       X2=X2, X2.2=(X2 - mean_X2)^2, X2.3=(X2 - mean_X2)^3,
                       X12=(X1 - mean_X1)*(X2 - mean_X2),
                       W=data$treat)

cef <- lm(Y ~ W * ., data = model.df)
sigma_y <- sqrt(summary(cef)$sigma^2)

X1 <- pointsALL$latitude[45]
X2 <- pointsALL$longitude[45]
truth <- predict(cef, data.frame(X1=X1, X1.2=(X1 - mean_X1)^2, X1.3=(X1 - mean_X1)^3,
                                 X2=X2, X2.2=(X2 - mean_X2)^2, X2.3=(X2 - mean_X2)^3,
                                 X12=(X1 - mean_X1)*(X2 - mean_X2), 
                                 W=1)) - predict(cef, data.frame(X1=X1, X1.2=(X1 - mean_X1)^2, X1.3=(X1 - mean_X1)^3,
                                                                 X2=X2, X2.2=(X2 - mean_X2)^2, X2.3=(X2 - mean_X2)^3,
                                                                 X12=(X1 - mean_X1)*(X2 - mean_X2), 
                                                                 W=0))

#### MC Simulation ####
t <- 1000 # Number of Monte Carlo iterations

# Number of observations in each Monte Carlo iteration \in {1000,5000,10000,50000}
for (n in c(1000,5000,10000,50000,1e5)){
  i <- 1
  result_optrdd <- c()
  bias_optrdd <- c()
  coverage_optrdd <- c()
  CI_plusminus <- c()
  
  while (i <= t){
    print(i)
    
    # generate simulation data
    outcome <- "price"
    datause <- genDGP(n, "price")
    WV = datause$treat
    XV = (pi * 6371 /180) * cbind(
      A=datause$longitude - mean(datause$longitude),
      B=datause$latitude - mean(datause$latitude))
    
    estimation.point = (pi * 6371 /180) *
      c(pointsALL$longitude[BOUNDARY_PT] - mean(datause$longitude),
        pointsALL$latitude[BOUNDARY_PT] - mean(datause$latitude))
    
    YV = datause[,outcome]
    last <- length(result_optrdd)
    try({
      max.curv = get_curvature(XV, outcome, WV)
      out.pt = optrdd(X=XV, Y=YV, W=WV, max.second.derivative = max.curv, estimation.point = estimation.point, verbose = FALSE)
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
    
    try({print("----------------")
      print(paste0("n=",n))
      print(paste0("Prediction optrdd: ", prediction))
      print("----------------")
      print(paste0("Mean optrdd:", mean(result_optrdd, na.rm = TRUE)))
      print("----------------")
      print(paste0("Bias optrdd:", mean(bias_optrdd, na.rm = TRUE)))
      print("----------------")
      print(paste0("Var optrdd: ", var(result_optrdd, na.rm = TRUE)))
      print("----------------")
      print(paste0("CR (truth=", truth, "): ", mean(coverage_optrdd, na.rm = TRUE)))
      print("----------------")
    })
    
    assign(paste0("result",n), as.data.frame(cbind("n"=n, "optrdd"=result_optrdd, "optrdd.bias"=bias_optrdd, "optrdd.cr"=coverage_optrdd, "CI_plusminus"=CI_plusminus)))
    
    i <- i+1
  }
}

result <- rbind(result1000,result5000,result10000,result50000,`result1e+05`)
result$truth <- truth
write.csv(result, "../Results csv/poly-optrdd-price.csv", row.names = F)
