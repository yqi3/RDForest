######################################
# Figure 2: Univariate score using
# data from Lee (2008)
######################################
rm(list=ls())
set.seed(12345)

library(foreign)
library(rdrobust)
library(ggplot2)
library(gridExtra)
library(ggpubr)

# Lee (2008) data obtained from http://economics.mit.edu/faculty/angrist/data1/mhe
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
d <- read.dta("../../Data/Lee2008/table_two_final.dta")

ptrue <- rdplot(y = d$demsharenext, x = d$difdemshare, binselect = "esmv",
                title = "Panel A: True Data",
                y.label = "Democratic vote share in election at time t+1",
                x.label = "Democratic victory margin at time t", 
                col.dots = "deepskyblue3", col.lines="black")

# Fitted 5-th order polynomial based on data from Lee (2008). 
# Coefficients are obtained by Calonico, Cattaneo, and Titiunik (2014), 
# Supplement Section S.3.1.2.: 
# https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA--Supplemental.pdf. 
# Imbens and Kalyanaraman (2011) obtained the same coefficients in their simulation section.

c_r0 = 0.52
c_r1 = 0.84   
c_r2 = -3 
c_r3 = 7.99 
c_r4 = -9.01 
c_r5 = 3.56
c_l0 = 0.48 
c_l1 = 1.27 
c_l2 = 7.18 
c_l3 = 20.21 
c_l4 = 21.54 
c_l5 = 7.33
c_rz = 0
c_lz = 0
sigma_y = 0.1295

#### Generate polynomial simulation data ####
n <- length(d$difdemshare) # sample size
X <- 2*rbeta(n,2,4)-1  # generate running variable
Y <- ifelse(X>=0, c_r0+c_r1*X+c_r2*X^2+c_r3*X^3+c_r4*X^4+c_r5*X^5, c_l0+c_l1*X+c_l2*X^2+c_l3*X^3+c_l4*X^4+c_l5*X^5)+rnorm(n,sd=sigma_y)
simdata <- as.data.frame(cbind("Y"=Y, "X"=X))

psim_poly <- rdplot(y =Y, x = X, binselect = "esmv",
                    title = "Panel B: Simulation Data (Polynomial)",
                    y.label = "Democratic vote share in election at time t+1",
                    x.label = "Democratic victory margin at time t", 
                    col.dots = "deepskyblue3", col.lines="black", x.lim = c(-1,1), y.lim=c(0,1))

#### Generate WGAN simulation data ####
wgan_d <- read.csv("../../Data/Lee2008/generated_Lee.csv")
rand <- sample(nrow(wgan_d), n, replace = F)
wgan_sample <- wgan_d[rand,]

psim_wgan <- rdplot(y =wgan_sample$Y, x = wgan_sample$X, binselect = "esmv",
                    title = "Panel C: Simulation Data (WGAN)",
                    y.label = "Democratic vote share in election at time t+1",
                    x.label = "Democratic victory margin at time t", 
                    col.dots = "deepskyblue3", col.lines="black", x.lim = c(-1,1), y.lim=c(0,1))

fig <- arrangeGrob(ptrue$rdplot, psim_poly$rdplot, psim_wgan$rdplot, ncol=3, top = text_grob("Univariate Design: Lee (2008)", face = "bold", size = 20))
ggsave("Lee_vis.png", fig, width = 16, height = 5)
