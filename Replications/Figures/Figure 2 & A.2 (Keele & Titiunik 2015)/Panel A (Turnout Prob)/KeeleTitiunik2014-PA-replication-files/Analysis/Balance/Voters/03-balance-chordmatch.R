#######################################
#
# Test post-matching balance for full dataset
# matching was done in file 01-analysis.Rout
# using Full_Spatial_Data.dta
#
# 2/28/2011
#
######################################
rm(list = ls())
options(width=300)
library(foreign)
library(Matching)
library(RItools)

## Load matched dataset
setwd("~/Documents/Geographic Discontinuities/Media Effects/Analysis/NJ/Matching/Voters/Results")
load(file="chordmatch_fulldata.RData") # object: finalres

## Open original, unmatched dataset
setwd("~/Documents/Geographic Discontinuities/Media Effects/Data/NJ/Voter File")
data <- read.dta("./Voters_Final.dta")
data <- data[complete.cases(data$treat),]
dim(data)
sum(data$latitud!=0 & data$longitud!=0)
data <- data[data$latitud!=0 & data$longitud!=0,]
dim(data)

## keep full matched dataset
indx.tr <- results[,"indx.tr"]
indx.co <- results[,"indx.co"]

data.tr <- data[indx.tr, ]
data.co <- data[indx.co, ]

# covariates we would like to have balance on 
Btr <- data.tr[,c("black", "asian", "hisp", "ind", "dem", "female", "educ", "age",  "income", "sfh", "sparents", "poverty", "college", "homeown", "bluecollar", "prof", "whitecollar", "unemp")] 

Bco <- data.co[,c("black", "asian", "hisp", "ind", "dem", "female", "educ", "age",  "income", "sfh", "sparents", "poverty", "college", "homeown", "bluecollar", "prof", "whitecollar", "unemp")] 


#Global Balance Test
match.data <- rbind(data.tr,data.co)
global.test <- xBalance(treat ~ black + asian + hisp + ind + dem + female + educ + age + income + sfh + sparents + poverty + college + homeown + bluecollar + prof + whitecollar, data=match.data, report=c("chisquare.test"))

rm(data)
rm(data.tr)
rm(data.co)
gc()

##MatchBalance(treat ~ as.matrix(B), ks = TRUE, nboots=0)
 
  Nt = nrow(Btr)
  Nc = nrow(Bco)
  
  colnames = c("Var" , "MeanT", "MeanC", "Pval ttest", "Var ratio", "Pval KS", "QQ med diff")
  resbal = matrix(NA, nrow=ncol(Btr)*2 + 2, ncol = 7, dimnames = list(NULL, colnames))
  for(i in 1:ncol(Btr)) {
    
    #cat("Calculating balance tests for variable", i, "of", ncol(Btr), "\n")
    indx = (!is.na(Btr[,i]) & !is.na(Bco[,i]))
    Nmiss = sum(!indx)            
    xtr  = Btr[indx,i]
    xco  = Bco[indx,i]
    
    bal<-balanceUV(xtr,xco, ks=TRUE,nboots = 1000, paired=TRUE, match=FALSE)
    
    sigmat<-sqrt(bal$var.Tr)      
    sigmac<-sqrt(bal$var.Co)
    #sigmaboth<-sqrt(var(c(Btr[,i],Bco[,i]))*(1/Nt+1/Nc))
    
    resbal[2*(i-1)+1,1] <- names(Btr)[i]
    resbal[2*(i-1)+1,2] <- round(bal$mean.Tr,digits=3)
    resbal[2*(i-1)+2,2] <- paste("(",round(sigmat/sqrt(Nt),digits=3),")",sep="")
    resbal[2*(i-1)+1,3] <- round(bal$mean.Co,digits=3)
    resbal[2*(i-1)+2,3] <- paste("(",round(sigmac/sqrt(Nc),digits=3),")",sep="")   # st err of mean.Co
    resbal[2*(i-1)+1,4] <- round(bal$p.value,digits=3)
    resbal[2*(i-1)+1,5] <- round(bal$var.ratio,digits=3)
    resbal[2*(i-1)+1,6] <- round(bal$ks$ks.boot.pvalue,digits=3)  # bootstrapped KS  pvalue
    resbal[2*(i-1)+1,7] <- round(bal$qqsummary$mediandiff,digits=3)  # bootstrapped KS  pvalue
  }
  
  resbal[ncol(Btr)*2 + 1, 1] <- "Sample size (T,C,miss)"
  resbal[ncol(Btr)*2 + 1, 2] <- Nt
  resbal[ncol(Btr)*2 + 1, 3] <- Nc
  resbal[ncol(Btr)*2 + 1, 4] <- Nmiss
  resbal[ncol(Btr)*2 + 2, 1] <- "Global Balance Test"
  resbal[ncol(Btr)*2 + 2, 2] <- as.numeric(global.test$overall[1])
  resbal[ncol(Btr)*2 + 2, 3] <- as.numeric(global.test$overall[3])
  
  print(resbal)
    
setwd("~/Documents/Geographic Discontinuities/Media Effects/Analysis/NJ/Balance/Voters/Results") 
save(resbal, file="balance_postmatch_fulldata.RData")
