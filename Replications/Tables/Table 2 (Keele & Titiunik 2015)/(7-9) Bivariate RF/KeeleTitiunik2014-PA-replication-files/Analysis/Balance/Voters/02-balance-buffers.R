#######################################
#
# Test balance for different buffer bands using Buffer_Spatial_Data.dta
#
# 2/28/2011
#
######################################
rm(list = ls())
options(width=300)
library(foreign)
library(Matching)
library(RITools)

setwd("~/Documents/Geographic Discontinuities/Media Effects/Data/NJ/Voter File")
data <- read.dta("./Voters_Final.dta")
data <- data[complete.cases(data$treat),]
dim(data)
sum(data$latitud!=0 & data$longitud!=0)
data <- data[data$latitud!=0 & data$longitud!=0,]
dim(data)

treat <- data[,"treat"]
buffer <- cbind(data[,"buffer1000"],data[,"buffer900"], data[,"buffer800"], data[,"buffer700"], data[,"buffer600"], data[,"buffer500"], data[,"buffer400"],data[,"buffer300"],data[,"buffer200"],data[,"buffer100"])
buffer.val <- c(1000, 900, 800, 700, 600, 500, 400, 300, 200, 100)

# covariates we would like to have balance on 
B <- data[,c("black", "asian", "hisp", "ind", "dem", "female", "educ", "age",  "income", "sfh", "sparents", "poverty", "college", "homeown", "bluecollar", "prof", "whitecollar", "unemp")] 

buf.dat <- cbind(treat,B)

rm(data)
gc()

##MatchBalance(treat ~ as.matrix(B), ks = TRUE, nboots=0)

## "pre-allocate" an empty list of length ncol(buffer)
finalres <- vector("list", ncol(buffer))
paste("buff", buffer.val, collapse=",", sep="")
names(finalres) <- c("buff1000","buff900","buff800","buff700","buff600","buff500","buff400","buff300","buff200","buff100")

for(j in 1:ncol(buffer)) {
  cat("Calculating balance testst for buffer", buffer.val[j], "\n")  
  Nt = sum(treat==1 & buffer[,j] == 1)
  Nc = sum(treat==0 & buffer[,j] == 1)
  
  colnames = c("Var" , "MeanT", "MeanC", "Pval ttest", "Var ratio", "Pval KS", "QQ med diff")
  resbal = matrix(NA, nrow=ncol(B)*2 + 2, ncol = 7, dimnames = list(NULL, colnames))
  
  buf.subset <- buf.dat[buffer[,j]==1,]
  buf.subset <- buf.subset[complete.cases(buf.subset$treat),]
  global.test <- xBalance(treat ~ black + asian + hisp + ind + dem + female + educ + age + income + sfh + sparents + poverty + college + homeown + bluecollar + prof + whitecollar,   data=buf.subset, report=c("chisquare.test"))

  for(i in 1:ncol(B)) {
    
    cat("Calculating balance tests for variable", i, "of", ncol(B), "\n")
    indx = (!is.na(B[,i]) & treat==1 & buffer[,j] == 1)
    xtr  = B[indx,i]
    indx = (!is.na(B[,i]) & treat==0 & buffer[,j] == 1)
    xco  = B[indx,i]
    
    bal<-balanceUV(xtr,xco, ks=TRUE,nboots = 1000, paired=FALSE, match=FALSE)
    
    sigmat<-sqrt(bal$var.Tr)      
    sigmac<-sqrt(bal$var.Co)
    #sigmaboth<-sqrt(var(B[,i])*(1/Nt+1/Nc))
    
    resbal[2*(i-1)+1,1] <- names(B)[i]
    resbal[2*(i-1)+1,2] <- round(bal$mean.Tr,digits=3)
    resbal[2*(i-1)+2,2] <- paste("(",round(sigmat/sqrt(Nt),digits=3),")",sep="")
    resbal[2*(i-1)+1,3] <- round(bal$mean.Co,digits=3)
    resbal[2*(i-1)+2,3] <- paste("(",round(sigmac/sqrt(Nc),digits=3),")",sep="")   # st err of mean.Co
    resbal[2*(i-1)+1,4] <- round(bal$p.value,digits=3)
    resbal[2*(i-1)+1,5] <- round(bal$var.ratio,digits=3)
    resbal[2*(i-1)+1,6] <- round(bal$ks$ks.boot.pvalue,digits=3)  # bootstrapped KS  pvalue
    resbal[2*(i-1)+1,7] <- round(bal$qqsummary$mediandiff,digits=3)  # bootstrapped KS  pvalue
  }
  resbal[ncol(B)*2 + 2, 1] <- "Sample size (T,C)"
  resbal[ncol(B)*2 + 2, 2] <- Nt
  resbal[ncol(B)*2 + 2, 3] <- Nc
  resbal[ncol(B)*2 + 1, 1] <- "Global Balance Test"
  resbal[ncol(B)*2 + 1, 2] <- as.numeric(global.test$overall[1])
  resbal[ncol(B)*2 + 1, 3] <- as.numeric(global.test$overall[3])
  
  cat("Printing results for buffer", buffer.val[j], "\n")
  print(resbal)
      
  #finalres = c(finalres, list(resbal))
  finalres[[j]] = resbal
}

setwd("~/Documents/Geographic Discontinuities/Media Effects/Analysis/NJ/Balance/Voters/Results") 
save(finalres, file="voters_balance_buffers.RData")


