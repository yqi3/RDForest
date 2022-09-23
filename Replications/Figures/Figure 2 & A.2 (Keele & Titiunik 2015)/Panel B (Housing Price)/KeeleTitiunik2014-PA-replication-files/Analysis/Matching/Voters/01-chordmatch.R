#######################################
#
# First matching analysis in Full_Spatial_Data.dta
#
# 2/28/2011
#
######################################
rm(list = ls())
library(foreign)
library(Matching)

setwd("~/Documents/Geographic Discontinuities/Media Effects/Analysis")
source("./distance-functions.R")
set.seed(984566)

setwd("~/Documents/Geographic Discontinuities/Media Effects/Data/NJ/Voter File")
data <- read.dta("./Voters_Final.dta")
data <- data[complete.cases(data$treat),]
dim(data)
sum(data$latitud!=0 & data$longitud!=0)
data <- data[data$latitud!=0 & data$longitud!=0,]
dim(data)

# create unique ID to be to identify the observations in order
data$id <- 1:nrow(data)

dat <- data[,c("latitude", "longitude", "treat", "id")]

# free memory
rm(data)
gc()

# Test balance before matching ==> see 02-balance-pre-matching.Rout

# Match on chord distance
Ntr = sum(dat$treat == 1)
Nco = sum(dat$treat == 0)
N = length(dat$treat)
# careful! longitud is negative, but the distance functions we use need all coordinates to be positive, take abs()
lat.tr = abs(dat[dat$treat==1,"latitude"])
lat.co = abs(dat[dat$treat==0,"latitude"])
lon.tr = abs(dat[dat$treat==1,"longitude"])
lon.co = abs(dat[dat$treat==0,"longitude"])
indx.tr = dat[dat$treat==1,"id"]
indx.co = dat[dat$treat==0,"id"]
indx.tr.relative = 1:Nco
  
  colnames <- c("i","indx.tr", "indx.co", "indx.co.rel", "distance")
  results <- matrix(NA, nrow=Nco, ncol = length(colnames), dimnames = list(NULL, colnames))
  
  for(i in 1:Nco) {
    t0 = proc.time()[3]
    dist = as.vector(disvec(lat.co[i], lon.co[i], lat.tr, lon.tr)$chord)
    mindist = min(dist)
    id.tr = indx.tr[dist==mindist]
    id.tr.rel = indx.tr.relative[dist==mindist]
    n = length(id.tr)  
    if(n > 1) {
      k = sample(1:n, 1)
      id.tr = id.tr[k]   # break ties randomly
      id.tr.rel = id.tr.rel[k]
    }
    ##results <- rbind(results, cbind(indx.tr = i, indx.co = id.co, distance = mindist))
    results[i,] <- c(i,indx.co[i],id.tr,id.tr.rel,mindist)
    t1 = proc.time()[3]
    cat("Found match for treatment", i, "of", Nco, " in ", t1-t0, "seconds -- Matched control is", id.tr, "\n")
}

setwd("~/Documents/Geographic Discontinuities/Media Effects/Analysis/NJ/Matching/Voters/Results")
save(results, file="chordmatch_fulldata.RData")
