######################################
# Simulation: Bivariate score using
# data from Keele & Titiunik (2015)
######################################
set.seed(1234)

library(foreign)
library(grf)
library(tidymodels)
library(rdrobust)
library(splines)
library(fields)
source("KeeleTitiunik2014-PA-replication-files/Analysis/Local Poly/distance-functions.R")
source("Helper_Fns.R")

### Read sample data ###
library(haven)
data <- read_dta("KeeleTitiunik2014-PA-replication-files/Data/Voters/Voters_Final.dta")
data <- data[complete.cases(data$latitude, data$longitude, data$treat, data[,c("age")]),]

### Get focal point on the treatment boundary where we estimate treatment effect ###
pointsALL <- read.dbf("KeeleTitiunik2014-PA-replication-files/Data/BorderSegmentPoints_Project.dbf") # boundary data
pointsALL$latitude= pointsALL$POINT_Y
pointsALL$longitude= pointsALL$POINT_X
focal_rdrobust <- matrix(c(pointsALL$latitude[45],pointsALL$longitude[45]),nrow = 1)

X1 <- data$latitude
X2 <- data$longitude
age <- data$age

model.df <- data.frame(Y=age, X1=X1, X1.2=(X1 - mean(X1))^2, X1.3=(X1 - mean(X1))^3,
                       X2=X2, X2.2=(X2 - mean(X2))^2, X2.3=(X2 - mean(X2))^3,
                       X12=(X1 - mean(X1))*(X2 - mean(X2)),
                       W=data$treat)

cef <- lm(Y ~ W * ., data = model.df)
sigma_y <- sqrt(summary(cef)$sigma^2)

X1 <- pointsALL$latitude[45]
X2 <- pointsALL$longitude[45]
truth <- predict(cef, data.frame(X1=X1, X1.2=(X1 - mean(X1))^2, X1.3=(X1 - mean(X1))^3,
                                 X2=X2, X2.2=(X2 - mean(X2))^2, X2.3=(X2 - mean(X2))^3,
                                 X12=(X1 - mean(X1))*(X2 - mean(X2)), 
                                 W=1)) - predict(cef, data.frame(X1=X1, X1.2=(X1 - mean(X1))^2, X1.3=(X1 - mean(X1))^3,
                                                                 X2=X2, X2.2=(X2 - mean(X2))^2, X2.3=(X2 - mean(X2))^3,
                                                                 X12=(X1 - mean(X1))*(X2 - mean(X2)), 
                                                                 W=0))

#### MC Simulation ####
t <- 1000 # number of iterations
n <- 5000 # sample size
result_rf <- c()
coverage_rf_nocov <- c()
BLB1 <- c()
BLB0 <- c()
focal <- matrix(0)

for (i in 1:t){
  # generate simulation data
  datause <- genDGP(n, "age")
  
  #### rdrobust following KT's chordal distance calculation ####
  
  ## Get covariates: latitude and longitude
  x1 <- datause$latitude
  x2 <- datause$longitude
  
  # Get outcome
  y <- datause$age
  treat <- datause$treat
  
  ## Set up treated data
  indx = (treat == 1)
  x1.Tr <- x1[indx]
  x2.Tr <- x2[indx]
  y.Tr  <- y[indx]
  
  ## Set up control data
  indx = (treat == 0)
  x1.Co <- x1[indx]
  x2.Co <- x2[indx]
  y.Co  <- y[indx]
  
  ####################################
  #### Calculate score in treated area
  ####################################
  
  ## Calculate chordal distance between each (x1,x2) latitude-longitude pair in the data and the point of the boundary (b1,b2) where we are evaluating
  b1 = focal_rdrobust[1] # lat
  b2 = focal_rdrobust[2] # long
  outdis = disvec(lat1=b1, lon1=b2, lat2 = x1.Tr, lon2 = x2.Tr)
  disTr = as.vector(outdis$chord)
  
  ####################################
  #### Calculate score in control area
  ####################################
  ## Calculate chordal distance between each (x1,x2) latitude-longitude pair in the data and the point of the boundary (b1,b2) where we are evaluating
  outdis = disvec(lat1=b1, lon1=b2, lat2 = x1.Co, lon2 = x2.Co)
  disCo= as.vector(outdis$chord)
  
  score = c(disTr, -disCo)
  y = c(y.Tr, y.Co)
  
  ### RF ###
  datause$index <- 1:nrow(datause)
  trt <- datause[datause$treat==1,]
  trt$score <- disTr
  trt$y <- y.Tr
  ctrl <- datause[datause$treat==0,]
  ctrl$score <- -disCo
  ctrl$y <- y.Co
  
  # separate forests on each side of the cutoff
  forestT <- regression_forest(data.frame(trt$score), trt$y, mtry=1, num.trees=5000, tune.parameters = c("sample.fraction", "honesty.fraction", "alpha"))
  forestC <- regression_forest(data.frame(ctrl$score), ctrl$y, mtry=1, num.trees=5000, tune.parameters = c("sample.fraction", "honesty.fraction", "alpha"))

  # get estimate for treatment effect
  prediction <- as.numeric(predict(forestT, focal)-predict(forestC, focal))
  result_rf <- append(result_rf, prediction)
  
  #- coverage rate for BLB variance estimator under zero covariance
  BLB1 <- append(BLB1, predict(forestT, focal, estimate.variance = TRUE)$variance.estimates)
  BLB0 <- append(BLB0, predict(forestC, focal, estimate.variance = TRUE)$variance.estimates)
  variance_nocov <- BLB1[i]+BLB0[i]
  se_nocov <- sqrt(variance_nocov)
  
  if ((prediction-1.96*se_nocov <= truth) & (truth <= prediction+1.96*se_nocov)){
    coverage_rf_nocov <- append(coverage_rf_nocov, 1)
  }
  else{
    coverage_rf_nocov <- append(coverage_rf_nocov, 0)
  }
  
  print(i)
  try({print(paste0("RF: ", prediction))
    print("----------------")
    print(paste0("Mean of RF:", mean(result_rf, na.rm = TRUE)))
    print("----------------")
    print(paste0("Var of RF: ", var(result_rf, na.rm = TRUE)))
    print(paste0("Var of RF (nocov): ", variance_nocov))
    print("----------------")
    print(paste0("MSE of RF: ", (mean(result_rf,na.rm = TRUE)-truth)^2+var(result_rf,na.rm = TRUE)))
    print("----------------")
    print(paste0("CR of RF (nocov): ", mean(coverage_rf_nocov, na.rm = TRUE)))
    print("----------------")
  })
}
