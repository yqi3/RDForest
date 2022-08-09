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
data <- read_dta("KeeleTitiunik2014-PA-replication-files/Data/Housing/NJ_House_Final.dta")
data <- data[complete.cases(data$latitude, data$longitude, data$treat, data[,c("sqft_price")]),]

### Get focal point on the treatment boundary where we estimate treatment effect ###
pointsALL <- read.dbf("KeeleTitiunik2014-PA-replication-files/Data/BorderSegmentPoints_Project.dbf") # boundary data
pointsALL$latitude= pointsALL$POINT_Y
pointsALL$longitude= pointsALL$POINT_X
focal <- matrix(c(pointsALL$latitude[45],pointsALL$longitude[45]),nrow = 1)

X1 <- data$latitude
X2 <- data$longitude
price <- data$sqft_price

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
                                 X12=(X1 - mean(X1))*(X2 - mean(X2)), 
                                 W=1)) - predict(cef, data.frame(X1=X1, X1.2=(X1 - mean(X1))^2, X1.3=(X1 - mean(X1))^3,
                                                                 X2=X2, X2.2=(X2 - mean(X2))^2, X2.3=(X2 - mean(X2))^3,
                                                                 X12=(X1 - mean(X1))*(X2 - mean(X2)), W=0))

#### MC Simulation ####
t <- 1000 # number of iterations
n <- 5000 # sample size
result_rdrobust <- c() # vector to collect results
coverage_rdrobust <- c()

for (i in 1:t){
  # generate simulation data
  datause <- genDGP(n, "price")
  
  #### rdrobust following KT's chordal distance calculation ####

  ## Get covariates: latitude and longitude
  x1 <- datause$latitude
  x2 <- datause$longitude

  # Get outcome
  y <- datause$price
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
  b1 = focal[1] # lat
  b2 = focal[2] # long
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

  ### rdrobust ###
  RBC <- rdrobust(y=y, x=score, c=0)
  result_rdrobust <- append(result_rdrobust, RBC$Estimate[2])

  if ( (RBC$ci[3,1] <= truth) & (truth <= RBC$ci[3,2])){
    coverage_rdrobust <- append(coverage_rdrobust, 1)
  }
  else{
    coverage_rdrobust <- append(coverage_rdrobust, 0)
  }
  
  print(i)
  try({print("----------------")
    print(paste0("rdrobust: ", RBC$Estimate[2]))
    print("----------------")
    print(paste0("Mean of rdrobust:", mean(result_rdrobust, na.rm = TRUE)))
    print("----------------")
    print(paste0("Var of rdrobust: ", var(result_rdrobust, na.rm = TRUE)))
    print("----------------")
    print(paste0("MSE of rdrobust: ", (mean(result_rdrobust,na.rm = TRUE)-truth)^2+var(result_rdrobust,na.rm = TRUE)))
    print("----------------")
    print(paste0("CR of rdrobust: ", mean(coverage_rdrobust, na.rm = TRUE)))
  })
}

