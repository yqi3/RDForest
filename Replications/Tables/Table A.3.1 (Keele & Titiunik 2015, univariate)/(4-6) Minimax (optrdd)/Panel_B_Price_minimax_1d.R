######################################
# Table A.3.1: (4-6) Minimax-Optimal
# Panel B: Housing Price
# Simulation: Bivariate score transformed to univariate using data from 
# Keele & Titiunik (2015)

# We implement minimax-optimal estimator using optrdd on the following platform
# R version 4.1.3 (2022-03-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.5
######################################
rm(list=ls())
set.seed(1234)

library(foreign)
library(fields)
library(optrdd)
library(glmnet)
library(splines)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("KeeleTitiunik2014-PA-replication-files/Analysis/Local Poly/distance-functions.R")
source("KT_Helper_Fns.R")

### Read sample data ###
library(haven)
data <- read_dta("KeeleTitiunik2014-PA-replication-files/Data/Housing/NJ_House_Final.dta")
data <- data[complete.cases(data$latitude, data$longitude, data$treat, data[,c("sqft_price")]),]

### Get focal point on the treatment boundary where we estimate treatment effect ###
pointsALL <- read.dbf("KeeleTitiunik2014-PA-replication-files/Data/BorderSegmentPoints_Project.dbf") # boundary data
pointsALL$latitude= pointsALL$POINT_Y
pointsALL$longitude= pointsALL$POINT_X
focal_rdrobust <- matrix(c(pointsALL$latitude[45],pointsALL$longitude[45]),nrow = 1)

X1 <- data$latitude
X2 <- data$longitude
price <- data$sqft_price

# Price
model.df <- data.frame(Y=price, X1=X1, X1.2=(X1 - mean(X1))^2, X1.3=(X1 - mean(X1))^3,
                       X2=X2, X2.2=(X2 - mean(X2))^2, X2.3=(X2 - mean(X2))^3,
                       X12=(X1 - mean(X1))*(X2 - mean(X2)),
                       W=data$treat)

cef <- lm(Y ~ W * ., data = model.df)
sigma_y <- sqrt(summary(cef)$sigma^2)

X1 <- pointsALL$latitude[45]
X2 <- pointsALL$longitude[45]
truth <- predict(cef, data.frame(X1=X1, X1.2=(X1 - mean(X1))^2, X1.3=(X1 - mean(X1))^3,
                                 X2=X2, X2.2=(X2 - mean(X2))^2, X2.3=(X2 - mean(X2))^3,
                                 X12=(X1 - mean(X1))*(X2 - mean(X2)), W=1), type = "response") - predict(cef, data.frame(X1=X1, X1.2=(X1 - mean(X1))^2, X1.3=(X1 - mean(X1))^3,
                                                                                                                         X2=X2, X2.2=(X2 - mean(X2))^2, X2.3=(X2 - mean(X2))^3,
                                                                                                                         X12=(X1 - mean(X1))*(X2 - mean(X2)), W=0), type = "response")

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
focal <- 0

while (i <= t){
  print(i)
  
  # generate simulation data
  datause <- genDGP(n, "price")
  
  #### rdrobust following KT's chordal distance calculation ####
  
  ## Get covariates: latitude and longitude
  x1 <- datause$latitude
  x2 <- datause$longitude
  
  treat <- datause$treat
  
  ## Set up treated data
  indx = (treat == 1)
  x1.Tr <- x1[indx]
  x2.Tr <- x2[indx]
  
  ## Set up control data
  indx = (treat == 0)
  x1.Co <- x1[indx]
  x2.Co <- x2[indx]
  
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
  
  ## Combine chordal distance with outcome ##
  trt <- datause[datause$treat==1,]
  trt$X <- disTr
  ctrl <- datause[datause$treat==0,]
  ctrl$X <- -disCo
  model_data <- rbind(trt, ctrl)
  
  W <- model_data$treat
  X <- model_data$X
  X2 <-  X^2
  Y <- model_data[,4]
  RDF = data.frame(W=W, X=X, X2=X2, Y=Y)
  
  ## Guess second derivative bound using a second order polynomial ##
  m <- lm(Y ~ W * (X + X2), data = RDF)
  coef <- m$coefficients[c(4,6)]
  b <- max(4*coef)  # can modify this multiplier, results in the paper are generated using 4
  
  ## Run optrdd ##
  last <- length(result_optrdd)
  try({
    out.pt = optrdd(X=X, Y=Y, W=W, max.second.derivative = b, estimation.point = focal, verbose = FALSE)
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