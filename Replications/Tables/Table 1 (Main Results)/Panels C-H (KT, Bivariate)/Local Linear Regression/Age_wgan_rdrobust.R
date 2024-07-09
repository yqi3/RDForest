######################################
# Table 1, Panel H, Columns 1-5: Local Linear Regression for Age
# Simulation: Bivariate score using WGAN data based on Keele & Titiunik (2015)

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
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("../../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Analysis/Local Poly/distance-functions.R")
source("../../../../Helper Fns/KT_Helper_Fns.R")

### Get focal point on the treatment boundary where we estimate treatment effect ###
pointsALL <- read.dbf("../../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Data/BorderSegmentPoints_Project.dbf") # boundary data
pointsALL$latitude <- pointsALL$POINT_Y
pointsALL$longitude <- pointsALL$POINT_X
focal <- matrix(c(pointsALL$latitude[45],pointsALL$longitude[45]),nrow = 1)

### Read WGAN-generated population data
data_full <- read.csv("../../../../Data/Keele & Titiunik 2015/WGAN_Epoch1k/age_cWGAN.csv")
focal_df <- read.csv("../../../../Data/Keele & Titiunik 2015/WGAN_Epoch1k/focal_age.csv")

### Population truth
truthT <- mean(focal_df$Y[focal_df$W==1])
truthC <- mean(focal_df$Y[focal_df$W==0])
truth <- truthT - truthC

#### MC Simulation ####
t <- 1000 # Number of Monte Carlo iterations

# Number of observations in each Monte Carlo iteration \in {1000,5000,10000,50000}
for (n in c(1000,5000,10000,50000,1e5)){
  result_rdrobust <- c() # vector to collect results
  bias_rdrobust <- c()
  se_rdrobust <- c()
  coverage_rdrobust <- c()
  count <- length(result_rdrobust)

  while (count < t) {
    # generate simulation data
    rand <- sample.int(nrow(data_full), n, replace = F)
    datause <- data_full[rand,]
    datause <- datause[,c("X1","X2","W","Y")]
    
    #### rdrobust following KT's chordal distance calculation ####
    
    ## Get covariates: latitude and longitude
    x1 <- datause$X1
    x2 <- datause$X2
    
    # Get outcome
    y <- datause$Y
    treat <- datause$W
    
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
    try({RBC <- rdrobust(y=y, x=score, c=0)
    result_rdrobust <- append(result_rdrobust, RBC$Estimate[2])
    se_rdrobust <- append(se_rdrobust, RBC$Estimate[4])
    bias_rdrobust <- append(bias_rdrobust, RBC$Estimate[2]-truth)})
    
    if (count == length(result_rdrobust)) { next }
    
    if ( (RBC$ci[3,1] <= truth) & (truth <= RBC$ci[3,2])){
      coverage_rdrobust <- append(coverage_rdrobust, 1)
    }
    else{
      coverage_rdrobust <- append(coverage_rdrobust, 0)
    }
    
    print(count+1)
    try({print("----------------")
      print(paste0("n=", n))
      print("----------------")
      print(paste0("rdrobust: ", RBC$Estimate[2]))
      print("----------------")
      print(paste0("Mean of rdrobust:", mean(result_rdrobust, na.rm = TRUE)))
      print("----------------")
      print(paste0("Bias of rdrobust:", mean(bias_rdrobust, na.rm = TRUE)))
      print("----------------")
      print(paste0("Var of rdrobust: ", var(result_rdrobust, na.rm = TRUE)))
      print("----------------")
      print(paste0("MSE of rdrobust: ", (mean(result_rdrobust,na.rm = TRUE)-truth)^2+var(result_rdrobust,na.rm = TRUE)))
      print("----------------")
      print(paste0("CR of rdrobust: ", mean(coverage_rdrobust, na.rm = TRUE)))
      print("----------------")
    })
    
    count <- count + 1
  }
  
  assign(paste0("result",n), as.data.frame(cbind("n"=n, "rdrobust"=result_rdrobust, "rdrobust.bias"=bias_rdrobust, "rdrobust.cr"=coverage_rdrobust, "rdrobust.se"=se_rdrobust)))
}

result <- rbind(result1000,result5000,result10000,result50000,`result1e+05`)
result$truth <- truth
write.csv(result, "../Results csv/wgan-rdrobust-age.csv", row.names = F)