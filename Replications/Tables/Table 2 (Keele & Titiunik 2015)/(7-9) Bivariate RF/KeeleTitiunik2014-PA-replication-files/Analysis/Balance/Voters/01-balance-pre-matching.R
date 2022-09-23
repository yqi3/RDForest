#######################################
#
# First matching analysis in Full_Spatial_Data.dta
# Test balance before matching
#
# 2/28/2011
#
######################################
rm(list = ls())
options(width=300)
library(foreign)
library(Matching)
library(RItools)

setwd("~/Documents/Geographic Discontinuities/Media Effects/Data/NJ/Voter File")
data <- read.dta("./Voters_Final.dta")
data <- data[complete.cases(data$treat),]
dim(data)
sum(data$latitud!=0 & data$longitud!=0)
data <- data[data$latitud!=0 & data$longitud!=0,]
dim(data)

#Global Balance Test
data.bal <- data[complete.cases(data$treat), ]
global.test <- xBalance(treat ~ black + asian + hisp + ind + dem + female + educ + age + income + sfh + sparents + poverty + college + homeown + bluecollar + prof + whitecollar, data=data.bal, report=c("chisquare.test"))

treat <- data[,"treat"]

# covariates we would like to have balance on 
B <- data[,c("black", "asian", "hisp", "ind", "dem", "female", "educ", "age",  "income", "sfh", "sparents", "poverty", "college", "homeown", "bluecollar", "prof", "whitecollar", "unemp")] 
rm(data)
gc()

##MatchBalance(treat ~ as.matrix(B), ks = TRUE, nboots=0)
        
Nt = sum(treat==1)
Nc = sum(treat==0)

colnames = c("Var" , "MeanT", "MeanC", "Pval ttest", "Var ratio", "Pval KS", "QQ med diff")
resbal = matrix(NA, nrow=ncol(B)*2 + 2, ncol = 7, dimnames = list(NULL, colnames))
for(i in 1:ncol(B)) {

  #cat("Calculating balance tests for variable", i, "of", ncol(B), "\n")
  indx = (!is.na(B[,i]) & treat==1)
  xtr  = B[indx,i]
  indx = (!is.na(B[,i]) & treat==0)
  xco  = B[indx,i]
  
  trt <- c(rep(1, length(xtr)), rep(0,length(xco)))
  var <- c(xtr, xco)
  
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
resbal[ncol(B)*2 + 1, 1] <- "Sample size (T,C)"
resbal[ncol(B)*2 + 1, 2] <- Nt
resbal[ncol(B)*2 + 1, 3] <- Nc
resbal[ncol(B)*2 + 2, 1] <- "Global Balance Test"
resbal[ncol(B)*2 + 2, 2] <- as.numeric(global.test$overall[1])
resbal[ncol(B)*2 + 2, 3] <- as.numeric(global.test$overall[3])

print(resbal)
  
setwd("~/Documents/Geographic Discontinuities/Media Effects/Analysis/NJ/Balance/Voters/Results") 
save(resbal, file="prematch_balance_fulldata.RData")
