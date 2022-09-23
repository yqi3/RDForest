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

setwd("~/Documents/Geographic Discontinuities/Media Effects/Data/NJ/Housing")
data <- read.dta("./NJ_House_Final.dta")
dim(data)
sum(data$latitud!=0 & data$longitud!=0)
data <- data[data$latitud!=0 & data$longitud!=0,]
dim(data)

treat <- data[,"treat"]

# covariates we would like to have balance on 
B <- data[,c("salepr", "sqft_price", "lotarea", "bldglivingarea")] 
rm(data)
gc()

##MatchBalance(treat ~ as.matrix(B), ks = TRUE, nboots=0)
        
Nt = sum(treat==1)
Nc = sum(treat==0)

colnames = c("Var" , "MeanT", "MeanC", "Pval ttest", "Var ratio", "Pval KS", "QQ med diff")
resbal = matrix(NA, nrow=ncol(B)*2 + 1, ncol = 7, dimnames = list(NULL, colnames))
for(i in 1:ncol(B)) {

  #cat("Calculating balance tests for variable", i, "of", ncol(B), "\n")
  indx = (!is.na(B[,i]) & treat==1)
  xtr  = B[indx,i]
  indx = (!is.na(B[,i]) & treat==0)
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
resbal[ncol(B)*2 + 1, 1] <- "Sample size (T,C)"
resbal[ncol(B)*2 + 1, 2] <- Nt
resbal[ncol(B)*2 + 1, 3] <- Nc

print(resbal)

setwd("~/Documents/Geographic Discontinuities/Media Effects/Analysis/NJ/Balance/Housing/Results")  
save(resbal, file="prematch_balance_house.RData")
