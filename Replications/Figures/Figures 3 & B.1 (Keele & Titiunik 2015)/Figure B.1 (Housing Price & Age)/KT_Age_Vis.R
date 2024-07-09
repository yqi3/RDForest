######################################
# Bivariate score using
# data from Keele & Titiunik (2015)
# Figure A.2 Panel C: Age
######################################
rm(list=ls())

library(haven)
library(foreign)
library(tidyverse)
library(dplyr)
library(gridExtra)
library(ggplot2)
library(ggpubr)

set.seed(888)

### Read sample data ###
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
data <- read_dta("../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Data/Voters/Voters_Final.dta")
data <- data[complete.cases(data$latitude, data$longitude, data$treat, data[,"age"]),]
pointsALL <- read.dbf("../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Data/BorderSegmentPoints_Project.dbf") # boundary data
pointsALL$latitude= pointsALL$POINT_Y
pointsALL$longitude= pointsALL$POINT_X
source("../../../Helper Fns/KT_Helper_Fns.R")

X1 <- data$latitude
X2 <- data$longitude
age <- data$age
mean_X1 <- mean(X1)
mean_X2 <- mean(X2)
age <- data$age

model.df <- data.frame(Y=age, X1=X1, X1.2=(X1 - mean_X1)^2, X1.3=(X1 - mean_X1)^3,
                       X2=X2, X2.2=(X2 - mean_X2)^2, X2.3=(X2 - mean_X2)^3,
                       X12=(X1 - mean_X1)*(X2 - mean_X2),
                       W=data$treat)

cef <- lm(Y ~ W * ., data = model.df)
sigma_y <- sqrt(summary(cef)$sigma^2)

pred_Y <- predict(cef, model.df)
pred_plot <- data.frame(pred_Y, X1, X2)

# For better visualization
age[age>=80]=80
age[age<=25]=25
df_plot <- data.frame(age, X1, X2)

p1 <- ggplot() +
  geom_point(data = df_plot, aes(x=X2, y=X1, color = age), size=0.8) +
  geom_path(aes(pointsALL$longitude, pointsALL$latitude), size=0.8) +
  scale_color_gradient2(low = 'dodgerblue3', midpoint=50, mid="lavender",high = 'firebrick3', limits=c(20, 80)) +
  labs(title="Actual Age Distribution",
       x ="Longitude", y = "Latitude", color = "Age") + 
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
  geom_point(data = pointsALL[c(30, 45, 60),], aes(x=longitude, y=latitude), shape = 15, color = "black", size = 2, show.legend = FALSE) +
  geom_point(data = pointsALL[45,], aes(x=longitude, y=latitude), shape = 1, color = "green", size = 4, stroke=1.5, show.legend = FALSE) +
  xlim(min(pointsALL$longitude), max(pointsALL$longitude)) +
  ylim(min(pointsALL$latitude), max(pointsALL$latitude))

### Generate polynomial simulation data ###
datause <- genDGP(2500, "age")

# For better visualization
datause[datause$age>=80,]$age = 80
datause[datause$age<=20,]$age = 20

p_poly <- ggplot() +
  geom_point(data = datause, aes(x=longitude, y=latitude, color = age), size=0.8) +
  geom_path(aes(pointsALL$longitude, pointsALL$latitude), size=0.8) +
  scale_color_gradient2(low = 'dodgerblue3', midpoint=50, mid="lavender",high = 'firebrick3', limits=c(20, 80)) +
  labs(title="Simulated Age (Polynomial)",
       x ="Longitude", y = "Latitude", color = "Age") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
  geom_point(data = pointsALL[c(30, 45, 60),], aes(x=longitude, y=latitude), shape = 15, color = "black", size = 2, show.legend = FALSE) +
  geom_point(data = pointsALL[45,], aes(x=longitude, y=latitude), shape = 1, color = "green", size = 4, stroke=1.5, show.legend = FALSE) +
  xlim(min(pointsALL$longitude), max(pointsALL$longitude)) +
  ylim(min(pointsALL$latitude), max(pointsALL$latitude))

### Generate WGAN simulation data ###
wgan_d <- read.csv("../../../Data/Keele & Titiunik 2015/WGAN_Epoch1k/age_cWGAN.csv")
rand <- sample(nrow(wgan_d), 2500, replace = F)
wgan_sample <- wgan_d[rand,]

# For better visualization
wgan_sample[wgan_sample$Y>=80,]$Y = 80
wgan_sample[wgan_sample$Y<=20,]$Y = 20

p_wgan <- ggplot() +
  geom_point(data = wgan_sample, aes(x=X2, y=X1, color = Y), size=0.8) +
  geom_path(aes(pointsALL$longitude, pointsALL$latitude), size=0.8) +
  scale_color_gradient2(low = 'dodgerblue3', midpoint=50, mid="lavender",high = 'firebrick3', limits=c(20, 80)) +
  labs(title="Simulated Age (WGAN)",
       x ="Longitude", y = "Latitude", color = "Age") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
  geom_point(data = pointsALL[c(30, 45, 60),], aes(x=longitude, y=latitude), shape = 15, color = "black", size = 2, show.legend = FALSE) +
  geom_point(data = pointsALL[45,], aes(x=longitude, y=latitude), shape = 1, color = "green", size = 4, stroke=1.5, show.legend = FALSE) +
  xlim(min(pointsALL$longitude), max(pointsALL$longitude)) +
  ylim(min(pointsALL$latitude), max(pointsALL$latitude))

fig <- arrangeGrob(p1, p_poly, p_wgan, ncol = 3, top = text_grob("Bivariate Design: Keele & Titiunik (2015)", face = "bold", size = 20))
ggsave("KT_age_vis.png", fig, width = 16, height = 6)

