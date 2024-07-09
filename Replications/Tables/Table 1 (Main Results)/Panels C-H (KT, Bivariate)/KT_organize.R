######################################
# Table 1, Panels C-H: Organization of Results
# Simulation: Bivariate score using data from WGAN and polynomial DGPs
######################################
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Polynomial ####
### Read results raw data ###
rdrobust_T <- read.csv("Results csv/poly-rdrobust-turnout.csv")
rdrobust_P <- read.csv("Results csv/poly-rdrobust-price.csv")
rdrobust_A <- read.csv("Results csv/poly-rdrobust-age.csv")

optrdd_T <- read.csv("Results csv/poly-optrdd-turnout.csv")
optrdd_P <- read.csv("Results csv/poly-optrdd-price.csv")
optrdd_A <- read.csv("Results csv/poly-optrdd-age.csv")

rf_T <- read.csv("Results csv/poly-rf-turnout.csv")
rf_P <- read.csv("Results csv/poly-rf-price.csv")
rf_A <- read.csv("Results csv/poly-rf-age.csv")

llf_T <- read.csv("Results csv/poly-llf-turnout.csv")
llf_P <- read.csv("Results csv/poly-llf-price.csv")
llf_A <- read.csv("Results csv/poly-llf-age.csv")

### Organize ###
## rdrobust
rdrobust_results_T <- matrix(0, ncol = 5, nrow = 3)
colnames(rdrobust_results_T) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(rdrobust_results_T) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  rdrobust_results_T["bias",toString(n)] <- mean(rdrobust_T$rdrobust[rdrobust_T$n==n])-rdrobust_T$truth[1]
  rdrobust_results_T["variance",toString(n)] <- var(rdrobust_T$rdrobust[rdrobust_T$n==n])
  rdrobust_results_T["CR",toString(n)] <- mean(rdrobust_T$rdrobust.cr[rdrobust_T$n==n])
}

rdrobust_results_P <- matrix(0, ncol = 5, nrow = 3)
colnames(rdrobust_results_P) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(rdrobust_results_P) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  rdrobust_results_P["bias",toString(n)] <- mean(rdrobust_P$rdrobust[rdrobust_P$n==n])-rdrobust_P$truth[1]
  rdrobust_results_P["variance",toString(n)] <- var(rdrobust_P$rdrobust[rdrobust_P$n==n])
  rdrobust_results_P["CR",toString(n)] <- mean(rdrobust_P$rdrobust.cr[rdrobust_P$n==n])
}

rdrobust_results_A <- matrix(0, ncol = 5, nrow = 3)
colnames(rdrobust_results_A) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(rdrobust_results_A) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  rdrobust_results_A["bias",toString(n)] <- mean(rdrobust_A$rdrobust[rdrobust_A$n==n])-rdrobust_A$truth[1]
  rdrobust_results_A["variance",toString(n)] <- var(rdrobust_A$rdrobust[rdrobust_A$n==n])
  rdrobust_results_A["CR",toString(n)] <- mean(rdrobust_A$rdrobust.cr[rdrobust_A$n==n])
}

rdrobust_results_T <- round(rdrobust_results_T, 4)
rdrobust_results_P <- round(rdrobust_results_P, 4)
rdrobust_results_A <- round(rdrobust_results_A, 4)

## optrdd
optrdd_results_T <- matrix(0, ncol = 5, nrow = 3)
colnames(optrdd_results_T) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(optrdd_results_T) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  optrdd_results_T["bias",toString(n)] <- mean(optrdd_T$optrdd[optrdd_T$n==n])-optrdd_T$truth[1]
  optrdd_results_T["variance",toString(n)] <- var(optrdd_T$optrdd[optrdd_T$n==n])
  optrdd_results_T["CR",toString(n)] <- mean(optrdd_T$optrdd.cr[optrdd_T$n==n])
}

optrdd_results_P <- matrix(0, ncol = 5, nrow = 3)
colnames(optrdd_results_P) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(optrdd_results_P) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  optrdd_results_P["bias",toString(n)] <- mean(optrdd_P$optrdd[optrdd_P$n==n])-optrdd_P$truth[1]
  optrdd_results_P["variance",toString(n)] <- var(optrdd_P$optrdd[optrdd_P$n==n])
  optrdd_results_P["CR",toString(n)] <- mean(optrdd_P$optrdd.cr[optrdd_P$n==n])
}

optrdd_results_A <- matrix(0, ncol = 5, nrow = 3)
colnames(optrdd_results_A) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(optrdd_results_A) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  optrdd_results_A["bias",toString(n)] <- mean(optrdd_A$optrdd[optrdd_A$n==n])-optrdd_A$truth[1]
  optrdd_results_A["variance",toString(n)] <- var(optrdd_A$optrdd[optrdd_A$n==n])
  optrdd_results_A["CR",toString(n)] <- mean(optrdd_A$optrdd.cr[optrdd_A$n==n])
}

optrdd_results_T <- round(optrdd_results_T, 4)
optrdd_results_P <- round(optrdd_results_P, 4)
optrdd_results_A <- round(optrdd_results_A, 4)

## RF
rf_results_T <- matrix(0, ncol = 5, nrow = 3)
colnames(rf_results_T) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(rf_results_T) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  rf_results_T["bias",toString(n)] <- mean(rf_T$rf[rf_T$n==n])-rf_T$truth[1]
  rf_results_T["variance",toString(n)] <- var(rf_T$rf[rf_T$n==n])
  rf_results_T["CR",toString(n)] <- mean(rf_T$rf.cr[rf_T$n==n])
}

rf_results_P <- matrix(0, ncol = 5, nrow = 3)
colnames(rf_results_P) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(rf_results_P) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  rf_results_P["bias",toString(n)] <- mean(rf_P$rf[rf_P$n==n])-rf_P$truth[1]
  rf_results_P["variance",toString(n)] <- var(rf_P$rf[rf_P$n==n])
  rf_results_P["CR",toString(n)] <- mean(rf_P$rf.cr[rf_P$n==n])
}

rf_results_A <- matrix(0, ncol = 5, nrow = 3)
colnames(rf_results_A) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(rf_results_A) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  rf_results_A["bias",toString(n)] <- mean(rf_A$rf[rf_A$n==n])-rf_A$truth[1]
  rf_results_A["variance",toString(n)] <- var(rf_A$rf[rf_A$n==n])
  rf_results_A["CR",toString(n)] <- mean(rf_A$rf.cr[rf_A$n==n])
}

rf_results_T <- round(rf_results_T, 4)
rf_results_P <- round(rf_results_P, 4)
rf_results_A <- round(rf_results_A, 4)

## LLF
llf_results_T <- matrix(0, ncol = 5, nrow = 3)
colnames(llf_results_T) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(llf_results_T) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  llf_results_T["bias",toString(n)] <- mean(llf_T$llf[llf_T$n==n])-llf_T$truth[1]
  llf_results_T["variance",toString(n)] <- var(llf_T$llf[llf_T$n==n])
  llf_results_T["CR",toString(n)] <- mean(llf_T$llf.cr[llf_T$n==n])
}

llf_results_P <- matrix(0, ncol = 5, nrow = 3)
colnames(llf_results_P) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(llf_results_P) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  llf_results_P["bias",toString(n)] <- mean(llf_P$llf[llf_P$n==n])-llf_P$truth[1]
  llf_results_P["variance",toString(n)] <- var(llf_P$llf[llf_P$n==n])
  llf_results_P["CR",toString(n)] <- mean(llf_P$llf.cr[llf_P$n==n])
}

llf_results_A <- matrix(0, ncol = 5, nrow = 3)
colnames(llf_results_A) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(llf_results_A) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  llf_results_A["bias",toString(n)] <- mean(llf_A$llf[llf_A$n==n])-llf_A$truth[1]
  llf_results_A["variance",toString(n)] <- var(llf_A$llf[llf_A$n==n])
  llf_results_A["CR",toString(n)] <- mean(llf_A$llf.cr[llf_A$n==n])
}

llf_results_T <- round(llf_results_T, 4)
llf_results_P <- round(llf_results_P, 4)
llf_results_A <- round(llf_results_A, 4)

### Write results as .csv
write.csv(rbind(cbind(rdrobust_results_T,optrdd_results_T,rf_results_T,llf_results_T),
                cbind(rdrobust_results_P,optrdd_results_P,rf_results_P,llf_results_P),
                cbind(rdrobust_results_A,optrdd_results_A,rf_results_A,llf_results_A)), "Panels_C-E_Organized.csv")


#### WGAN ####
### Read results raw data ###
rdrobust_T <- read.csv("Results csv/wgan-rdrobust-turnout.csv")
rdrobust_P <- read.csv("Results csv/wgan-rdrobust-price.csv")
rdrobust_A <- read.csv("Results csv/wgan-rdrobust-age.csv")

optrdd_T <- read.csv("Results csv/wgan-optrdd-turnout.csv")
optrdd_P <- read.csv("Results csv/wgan-optrdd-price.csv")
optrdd_A <- read.csv("Results csv/wgan-optrdd-age.csv")

rf_T <- read.csv("Results csv/wgan-rf-turnout.csv")
rf_P <- read.csv("Results csv/wgan-rf-price.csv")
rf_A <- read.csv("Results csv/wgan-rf-age.csv")

llf_T <- read.csv("Results csv/wgan-llf-turnout.csv")
llf_P <- read.csv("Results csv/wgan-llf-price.csv")
llf_A <- read.csv("Results csv/wgan-llf-age.csv")

### Organize ###
## rdrobust
rdrobust_results_T <- matrix(0, ncol = 5, nrow = 3)
colnames(rdrobust_results_T) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(rdrobust_results_T) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  rdrobust_results_T["bias",toString(n)] <- mean(rdrobust_T$rdrobust[rdrobust_T$n==n])-rdrobust_T$truth[1]
  rdrobust_results_T["variance",toString(n)] <- var(rdrobust_T$rdrobust[rdrobust_T$n==n])
  rdrobust_results_T["CR",toString(n)] <- mean(rdrobust_T$rdrobust.cr[rdrobust_T$n==n])
}

rdrobust_results_P <- matrix(0, ncol = 5, nrow = 3)
colnames(rdrobust_results_P) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(rdrobust_results_P) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  rdrobust_results_P["bias",toString(n)] <- mean(rdrobust_P$rdrobust[rdrobust_P$n==n])-rdrobust_P$truth[1]
  rdrobust_results_P["variance",toString(n)] <- var(rdrobust_P$rdrobust[rdrobust_P$n==n])
  rdrobust_results_P["CR",toString(n)] <- mean(rdrobust_P$rdrobust.cr[rdrobust_P$n==n])
}

rdrobust_results_A <- matrix(0, ncol = 5, nrow = 3)
colnames(rdrobust_results_A) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(rdrobust_results_A) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  rdrobust_results_A["bias",toString(n)] <- mean(rdrobust_A$rdrobust[rdrobust_A$n==n])-rdrobust_A$truth[1]
  rdrobust_results_A["variance",toString(n)] <- var(rdrobust_A$rdrobust[rdrobust_A$n==n])
  rdrobust_results_A["CR",toString(n)] <- mean(rdrobust_A$rdrobust.cr[rdrobust_A$n==n])
}

rdrobust_results_T <- round(rdrobust_results_T, 4)
rdrobust_results_P <- round(rdrobust_results_P, 4)
rdrobust_results_A <- round(rdrobust_results_A, 4)

## optrdd
optrdd_results_T <- matrix(0, ncol = 5, nrow = 3)
colnames(optrdd_results_T) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(optrdd_results_T) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  optrdd_results_T["bias",toString(n)] <- mean(optrdd_T$optrdd[optrdd_T$n==n])-optrdd_T$truth[1]
  optrdd_results_T["variance",toString(n)] <- var(optrdd_T$optrdd[optrdd_T$n==n])
  optrdd_results_T["CR",toString(n)] <- mean(optrdd_T$optrdd.cr[optrdd_T$n==n])
}

optrdd_results_P <- matrix(0, ncol = 5, nrow = 3)
colnames(optrdd_results_P) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(optrdd_results_P) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  optrdd_results_P["bias",toString(n)] <- mean(optrdd_P$optrdd[optrdd_P$n==n])-optrdd_P$truth[1]
  optrdd_results_P["variance",toString(n)] <- var(optrdd_P$optrdd[optrdd_P$n==n])
  optrdd_results_P["CR",toString(n)] <- mean(optrdd_P$optrdd.cr[optrdd_P$n==n])
}

optrdd_results_A <- matrix(0, ncol = 5, nrow = 3)
colnames(optrdd_results_A) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(optrdd_results_A) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  optrdd_results_A["bias",toString(n)] <- mean(optrdd_A$optrdd[optrdd_A$n==n])-optrdd_A$truth[1]
  optrdd_results_A["variance",toString(n)] <- var(optrdd_A$optrdd[optrdd_A$n==n])
  optrdd_results_A["CR",toString(n)] <- mean(optrdd_A$optrdd.cr[optrdd_A$n==n])
}

optrdd_results_T <- round(optrdd_results_T, 4)
optrdd_results_P <- round(optrdd_results_P, 4)
optrdd_results_A <- round(optrdd_results_A, 4)

## RF
rf_results_T <- matrix(0, ncol = 5, nrow = 3)
colnames(rf_results_T) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(rf_results_T) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  rf_results_T["bias",toString(n)] <- mean(rf_T$rf[rf_T$n==n])-rf_T$truth[1]
  rf_results_T["variance",toString(n)] <- var(rf_T$rf[rf_T$n==n])
  rf_results_T["CR",toString(n)] <- mean(rf_T$rf.cr[rf_T$n==n])
}

rf_results_P <- matrix(0, ncol = 5, nrow = 3)
colnames(rf_results_P) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(rf_results_P) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  rf_results_P["bias",toString(n)] <- mean(rf_P$rf[rf_P$n==n])-rf_P$truth[1]
  rf_results_P["variance",toString(n)] <- var(rf_P$rf[rf_P$n==n])
  rf_results_P["CR",toString(n)] <- mean(rf_P$rf.cr[rf_P$n==n])
}

rf_results_A <- matrix(0, ncol = 5, nrow = 3)
colnames(rf_results_A) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(rf_results_A) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  rf_results_A["bias",toString(n)] <- mean(rf_A$rf[rf_A$n==n])-rf_A$truth[1]
  rf_results_A["variance",toString(n)] <- var(rf_A$rf[rf_A$n==n])
  rf_results_A["CR",toString(n)] <- mean(rf_A$rf.cr[rf_A$n==n])
}

rf_results_T <- round(rf_results_T, 4)
rf_results_P <- round(rf_results_P, 4)
rf_results_A <- round(rf_results_A, 4)

## LLF
llf_results_T <- matrix(0, ncol = 5, nrow = 3)
colnames(llf_results_T) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(llf_results_T) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  llf_results_T["bias",toString(n)] <- mean(llf_T$llf[llf_T$n==n])-llf_T$truth[1]
  llf_results_T["variance",toString(n)] <- var(llf_T$llf[llf_T$n==n])
  llf_results_T["CR",toString(n)] <- mean(llf_T$llf.cr[llf_T$n==n])
}

llf_results_P <- matrix(0, ncol = 5, nrow = 3)
colnames(llf_results_P) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(llf_results_P) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  llf_results_P["bias",toString(n)] <- mean(llf_P$llf[llf_P$n==n])-llf_P$truth[1]
  llf_results_P["variance",toString(n)] <- var(llf_P$llf[llf_P$n==n])
  llf_results_P["CR",toString(n)] <- mean(llf_P$llf.cr[llf_P$n==n])
}

llf_results_A <- matrix(0, ncol = 5, nrow = 3)
colnames(llf_results_A) <- c("1000", "5000", "10000", "50000", "1e+05")
rownames(llf_results_A) <- c("bias", "variance", "CR")
for (n in c(1000,5000,10000,50000,1e5)) {
  llf_results_A["bias",toString(n)] <- mean(llf_A$llf[llf_A$n==n])-llf_A$truth[1]
  llf_results_A["variance",toString(n)] <- var(llf_A$llf[llf_A$n==n])
  llf_results_A["CR",toString(n)] <- mean(llf_A$llf.cr[llf_A$n==n])
}

llf_results_T <- round(llf_results_T, 4)
llf_results_P <- round(llf_results_P, 4)
llf_results_A <- round(llf_results_A, 4)

### Write results as .csv
write.csv(rbind(cbind(rdrobust_results_T,optrdd_results_T,rf_results_T,llf_results_T),
                cbind(rdrobust_results_P,optrdd_results_P,rf_results_P,llf_results_P),
                cbind(rdrobust_results_A,optrdd_results_A,rf_results_A,llf_results_A)), "Panels_F-H_Organized.csv")
