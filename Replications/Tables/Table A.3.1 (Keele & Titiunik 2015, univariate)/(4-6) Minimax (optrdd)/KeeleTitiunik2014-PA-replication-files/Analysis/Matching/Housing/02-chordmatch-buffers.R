#######################################
#
# Matching on chord distance within each buffer band
# using Buffer_Spatial_Data.dta
#
# 2/28/2011
#
######################################
rm(list = ls())
options(width=300)
library(foreign)
library(Matching)

setwd("~/Documents/Geographic Discontinuities/Media Effects/Analysis")
source("./distance-functions.R")
set.seed(856246)

setwd("~/Documents/Geographic Discontinuities/Media Effects/Data/NJ/Housing")
data <- read.dta("./NJ_House_Final.dta")
data <- data[complete.cases(data$buffer1000),]
dim(data)
sum(data$latitud!=0 & data$longitud!=0)
data <- data[data$latitud!=0 & data$longitud!=0,]
dim(data)

buffer.val <- c(1000, 500)
nbuff <- length(buffer.val)

# create unique ID to be to identify the observations in order
data$id <- 1:nrow(data)

buffer <- cbind(data[,"buffer1000"], data[,"buffer500"])
dat <- data[,c("latitude", "longitude", "treat", "id")]  

#free memory
rm(data)
gc()

## "pre-allocate" an empty list of length ncol(buffer)
finalres <- vector("list", nbuff)
paste("buff", buffer.val, collapse=",", sep="")
names(finalres) <- c("buff1000","buff500")

for(j in 1:nbuff) {
  ## Match on chord distance
  Ntr = sum(dat$treat == 1 & buffer[,j] == 1)
  Nco = sum(dat$treat == 0 & buffer[,j] == 1)
  N = Ntr + Nco
  
  ## careful! longitud is negative, but the distance functions we use need all coordinates to be positive, take abs()
  lat.tr = abs(dat[dat$treat==1 & buffer[,j] == 1,"latitude"])
  lat.co = abs(dat[dat$treat==0 & buffer[,j] == 1,"latitude"])
  lon.tr = abs(dat[dat$treat==1 & buffer[,j] == 1,"longitude"])
  lon.co = abs(dat[dat$treat==0 & buffer[,j] == 1,"longitude"])
  indx.tr = dat[dat$treat==1 & buffer[,j] == 1,"id"]
  indx.co = dat[dat$treat==0 & buffer[,j] == 1,"id"]
  indx.co.relative = 1:Nco
  
  colnames <- c("i","indx.tr", "indx.co", "indx.co.rel", "distance")
  results <- matrix(NA, nrow=Ntr, ncol = length(colnames), dimnames = list(NULL, colnames))
  for(i in 1:Ntr) {
    t0 = proc.time()[3]
    dist = as.vector(disvec(lat.tr[i], lon.tr[i], lat.co, lon.co)$chord)
    mindist = min(dist)
    id.co = indx.co[dist==mindist]
    id.co.rel = indx.co.relative[dist==mindist]
    n = length(id.co)  
    if(n > 1) {
      k = sample(1:n, 1)
      id.co = id.co[k]   # break ties randomly
      id.co.rel = id.co.rel[k]
    }
    ##results <- rbind(results, cbind(indx.tr = i, indx.co = id.co, distance = mindist))
    results[i,] <- c(i,indx.tr[i],id.co,id.co.rel,mindist)
    t1 = proc.time()[3]
    cat("Buffer", buffer.val[j], ": Found match for treatment", i, "of", Ntr, " in ", t1-t0, "seconds -- Matched control is", id.co, "\n")
  }
  finalres[[j]] = results
}

setwd("~/Documents/Geographic Discontinuities/Media Effects/Analysis/NJ/Matching/Housing")
save(finalres, file="chordmatch_bufferdata.RData")
