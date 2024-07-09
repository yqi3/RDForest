######################################
# Table 1, Panel C, Columns 1-5: Local Linear Regression for Turnout
# Simulation: Bivariate score using data from polynomial DGP

# We implement local linear regressions using rdrobust on the following platform
# R version 4.3.2 (2023-10-31)
# Platform: x86_64-pc-linux-gnu (64-bit)
# Running under: Debian GNU/Linux trixie/sid
######################################
rm(list=ls())
set.seed(1234)

library(foreign)
library(tidymodels)
library(rdrobust)
library(splines)
library(fields)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Analysis/Local Poly/distance-functions.R")
source("../../../../Helper Fns/KT_Helper_Fns.R")

### Read sample data ###
library(haven)
data <- read_dta("../../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Data/Voters/Voters_Final.dta")
data <- data[complete.cases(data$latitude, data$longitude, data$treat, data[,c("e2008g")]),]

### Get focal point on the treatment boundary where we estimate treatment effect ###
pointsALL <- read.dbf("../../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Data/BorderSegmentPoints_Project.dbf") # boundary data
pointsALL$latitude= pointsALL$POINT_Y
pointsALL$longitude= pointsALL$POINT_X
focal <- matrix(c(pointsALL$latitude[45],pointsALL$longitude[45]),nrow = 1)

X1 <- data$latitude
X2 <- data$longitude
mean_X1 <- mean(X1)
mean_X2 <- mean(X2)
e2008g <- data$e2008g

## Turnout
model.df <- data.frame(Y=e2008g, X1=X1, X1.2=(X1 - mean_X1)^2, X1.3=(X1 - mean_X1)^3,
                       X2=X2, X2.2=(X2 - mean_X2)^2, X2.3=(X2 - mean_X2)^3,
                       X12=(X1 - mean_X1)*(X2 - mean_X2),
                       W=data$treat)

cef <- glm(Y ~ W * ., family=binomial(link='logit'),data=model.df)

sigma_y <- sqrt(summary(cef)$sigma^2)

X1 <- pointsALL$latitude[45]
X2 <- pointsALL$longitude[45]
truth <- predict(cef, data.frame(X1=X1, X1.2=(X1 - mean_X1)^2, X1.3=(X1 - mean_X1)^3,
                                 X2=X2, X2.2=(X2 - mean_X2)^2, X2.3=(X2 - mean_X2)^3,
                                 X12=(X1 - mean_X1)*(X2 - mean_X2), W=1), type = "response") - predict(cef, data.frame(X1=X1, X1.2=(X1 - mean_X1)^2, X1.3=(X1 - mean_X1)^3,
                                                                                                                       X2=X2, X2.2=(X2 - mean_X2)^2, X2.3=(X2 - mean_X2)^3,
                                                                                                                       X12=(X1 - mean_X1)*(X2 - mean_X2), W=0), type = "response")

#### MC Simulation ####
t <- 1000 # Number of Monte Carlo iterations

# Number of observations in each Monte Carlo iteration \in {1000,5000,10000,50000}
for (n in c(1000,5000,10000,50000,1e5)){
  result_rdrobust <- c() # vector to collect results
  bias_rdrobust <- c()
  se_rdrobust <- c()
  coverage_rdrobust <- c()
  
  for (i in 1:t) {
    # generate simulation data
    datause <- genDGP(n, "turnout")
    
    #### rdrobust following KT's chordal distance calculation ####
  
    ## Get covariates: latitude and longitude
    x1 <- datause$latitude
    x2 <- datause$longitude
  
    # Get outcome
    y <- datause$e2008g
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
    se_rdrobust <- append(se_rdrobust, RBC$Estimate[4])
    bias_rdrobust <- append(bias_rdrobust, RBC$Estimate[2]-truth)
  
    if ( (RBC$ci[3,1] <= truth) & (truth <= RBC$ci[3,2])){
      coverage_rdrobust <- append(coverage_rdrobust, 1)
    }
    else{
      coverage_rdrobust <- append(coverage_rdrobust, 0)
    }
    
    print(i)
    try({print("----------------")
      print(paste0("n=", n))
      print(paste0("rdrobust: ", RBC$Estimate[2]))
      print(paste0("Mean of rdrobust:", mean(result_rdrobust, na.rm = TRUE)))
      print(paste0("Bias of rdrobust:", mean(bias_rdrobust, na.rm = TRUE)))
      print(paste0("Var of rdrobust: ", var(result_rdrobust, na.rm = TRUE)))
      print(paste0("MSE of rdrobust: ", (mean(result_rdrobust,na.rm = TRUE)-truth)^2+var(result_rdrobust,na.rm = TRUE)))
      print(paste0("CR of rdrobust: ", mean(coverage_rdrobust, na.rm = TRUE)))
      print("----------------")
    })
  }
  
  assign(paste0("result",n), as.data.frame(cbind("n"=n, "rdrobust"=result_rdrobust, "rdrobust.bias"=bias_rdrobust, "rdrobust.cr"=coverage_rdrobust, "rdrobust.se"=se_rdrobust)))
}

result <- rbind(result1000,result5000,result10000,result50000,`result1e+05`)
result$truth <- truth
write.csv(result, paste0("../Results csv/poly-rdrobust-turnout.csv"), row.names = F)
