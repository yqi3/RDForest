######################################
# Table 1, Panels A-B: Organization of Results
# Simulation: Univariate score using data from WGAN and polynomial DGPs
######################################
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Polynomial ####
### Read results raw data ###
rdrobust_lee <- read.csv("Results csv/poly-lee-rdrobust.csv")
optrdd_lee <- read.csv("Results csv/poly-lee-optrdd.csv")
rf_lee <- read.csv("Results csv/poly-rf-lee.csv")
llf_lee <- read.csv("Results csv/poly-llf-lee.csv")

### Organize ###
## rdrobust
rdrobust_results_lee <- matrix(0, ncol = 5, nrow = 3)
colnames(rdrobust_results_lee) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(rdrobust_results_lee) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  rdrobust_results_lee["bias",toString(n)] <- mean(rdrobust_lee$rdrobust[rdrobust_lee$n==n])-rdrobust_lee$truth[1]
  rdrobust_results_lee["variance",toString(n)] <- var(rdrobust_lee$rdrobust[rdrobust_lee$n==n])
  rdrobust_results_lee["CR",toString(n)] <- mean(rdrobust_lee$rdrobust.cr[rdrobust_lee$n==n])
}

rdrobust_results_lee <- round(rdrobust_results_lee, 4)

## optrdd
optrdd_results_lee <- matrix(0, ncol = 5, nrow = 3)
colnames(optrdd_results_lee) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(optrdd_results_lee) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  optrdd_results_lee["bias",toString(n)] <- mean(optrdd_lee$optrdd[optrdd_lee$n==n])-optrdd_lee$truth[1]
  optrdd_results_lee["variance",toString(n)] <- var(optrdd_lee$optrdd[optrdd_lee$n==n])
  optrdd_results_lee["CR",toString(n)] <- mean(optrdd_lee$optrdd.cr[optrdd_lee$n==n])
}

optrdd_results_lee <- round(optrdd_results_lee, 4)

## RF
rf_results_lee <- matrix(0, ncol = 5, nrow = 3)
colnames(rf_results_lee) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(rf_results_lee) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  rf_results_lee["bias",toString(n)] <- mean(rf_lee$rf[rf_lee$n==n])-rf_lee$truth[1]
  rf_results_lee["variance",toString(n)] <- var(rf_lee$rf[rf_lee$n==n])
  rf_results_lee["CR",toString(n)] <- mean(rf_lee$rf.cr[rf_lee$n==n])
}

rf_results_lee <- round(rf_results_lee, 4)

## LlF
llf_results_lee <- matrix(0, ncol = 5, nrow = 3)
colnames(llf_results_lee) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(llf_results_lee) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  llf_results_lee["bias",toString(n)] <- mean(llf_lee$llf[llf_lee$n==n])-llf_lee$truth[1]
  llf_results_lee["variance",toString(n)] <- var(llf_lee$llf[llf_lee$n==n])
  llf_results_lee["CR",toString(n)] <- mean(llf_lee$llf.cr[llf_lee$n==n])
}

llf_results_lee <- round(llf_results_lee, 4)

### Write results as .csv
write.csv(cbind(rdrobust_results_lee, optrdd_results_lee, rf_results_lee, llf_results_lee), "Panel_A_Organized.csv")


#### WGAN ####
### Read results raw data ###
rdrobust_lee <- read.csv("Results csv/wgan-lee-rdrobust.csv")
optrdd_lee <- read.csv("Results csv/wgan-lee-optrdd.csv")
rf_lee <- read.csv("Results csv/wgan-rf-lee.csv")
llf_lee <- read.csv("Results csv/wgan-llf-lee.csv")

### Organize ###
## rdrobust
rdrobust_results_lee <- matrix(0, ncol = 5, nrow = 3)
colnames(rdrobust_results_lee) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(rdrobust_results_lee) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  rdrobust_results_lee["bias",toString(n)] <- mean(rdrobust_lee$rdrobust[rdrobust_lee$n==n])-rdrobust_lee$truth[1]
  rdrobust_results_lee["variance",toString(n)] <- var(rdrobust_lee$rdrobust[rdrobust_lee$n==n])
  rdrobust_results_lee["CR",toString(n)] <- mean(rdrobust_lee$rdrobust.cr[rdrobust_lee$n==n])
}

rdrobust_results_lee <- round(rdrobust_results_lee, 4)

## optrdd
optrdd_results_lee <- matrix(0, ncol = 5, nrow = 3)
colnames(optrdd_results_lee) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(optrdd_results_lee) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  optrdd_results_lee["bias",toString(n)] <- mean(optrdd_lee$optrdd[optrdd_lee$n==n])-optrdd_lee$truth[1]
  optrdd_results_lee["variance",toString(n)] <- var(optrdd_lee$optrdd[optrdd_lee$n==n])
  optrdd_results_lee["CR",toString(n)] <- mean(optrdd_lee$optrdd.cr[optrdd_lee$n==n])
}

optrdd_results_lee <- round(optrdd_results_lee, 4)

## RF
rf_results_lee <- matrix(0, ncol = 5, nrow = 3)
colnames(rf_results_lee) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(rf_results_lee) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  rf_results_lee["bias",toString(n)] <- mean(rf_lee$rf[rf_lee$n==n])-rf_lee$truth[1]
  rf_results_lee["variance",toString(n)] <- var(rf_lee$rf[rf_lee$n==n])
  rf_results_lee["CR",toString(n)] <- mean(rf_lee$rf.cr[rf_lee$n==n])
}

rf_results_lee <- round(rf_results_lee, 4)

## LlF
llf_results_lee <- matrix(0, ncol = 5, nrow = 3)
colnames(llf_results_lee) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(llf_results_lee) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  llf_results_lee["bias",toString(n)] <- mean(llf_lee$llf[llf_lee$n==n])-llf_lee$truth[1]
  llf_results_lee["variance",toString(n)] <- var(llf_lee$llf[llf_lee$n==n])
  llf_results_lee["CR",toString(n)] <- mean(llf_lee$llf.cr[llf_lee$n==n])
}

llf_results_lee <- round(llf_results_lee, 4)

### Write results as .csv
write.csv(cbind(rdrobust_results_lee, optrdd_results_lee, rf_results_lee, llf_results_lee), "Panel_B_Organized.csv")
