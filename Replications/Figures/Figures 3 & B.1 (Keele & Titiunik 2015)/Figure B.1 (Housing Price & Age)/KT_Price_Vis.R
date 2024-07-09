######################################
# Bivariate score using
# data from Keele & Titiunik (2015)
# Figure B.1: Housing Price
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
data <- read_dta("../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Data/Housing/NJ_House_Final.dta")
data <- data[complete.cases(data$latitude, data$longitude, data$treat, data[,"sqft_price"]),]
pointsALL <- read.dbf("../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Data/BorderSegmentPoints_Project.dbf") # boundary data
pointsALL$latitude= pointsALL$POINT_Y
pointsALL$longitude= pointsALL$POINT_X
source("../../../Helper Fns/KT_Helper_Fns.R")

X1 <- data$latitude
X2 <- data$longitude
price <- data$sqft_price
mean_X1 <- mean(X1)
mean_X2 <- mean(X2)
price <- data$sqft_price

model.df <- data.frame(Y=price, X1=X1, X1.2=(X1 - mean_X1)^2, X1.3=(X1 - mean_X1)^3,
                       X2=X2, X2.2=(X2 - mean_X2)^2, X2.3=(X2 - mean_X2)^3,
                       X12=(X1 - mean_X1)*(X2 - mean_X2),
                       W=data$treat)

cef <- lm(Y ~ W * ., data = model.df)
sigma_y <- sqrt(summary(cef)$sigma^2)

# For better visualization
pred_Y <- predict(cef, model.df)
pred_Y[pred_Y>=260]=260
pred_Y[pred_Y<=180]=180
pred_plot <- data.frame(pred_Y, X1, X2)

price[price>=260]=260
price[price<=180]=180
df_plot <- data.frame(price, X1, X2)

p1 <- ggplot() +
  geom_point(data = df_plot, aes(x=X2, y=X1, color = price), size=0.8) +
  geom_path(aes(pointsALL$longitude, pointsALL$latitude), size=0.8) +
  scale_color_gradient2(low = 'dodgerblue3', midpoint=220, mid="lavender",high = 'firebrick3', limits=c(180, 260)) +
  labs(title="Actual Housing Price Distribution",
       x ="Longitude", y = "Latitude", color = "Price") + 
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
  geom_point(data = pointsALL[c(30, 45, 60),], aes(x=longitude, y=latitude), shape = 15, color = "black", size = 2, show.legend = FALSE) +
  geom_point(data = pointsALL[45,], aes(x=longitude, y=latitude), shape = 1, color = "green", size = 4, stroke=1.5, show.legend = FALSE) +
  xlim(min(pointsALL$longitude), max(pointsALL$longitude)) +
  ylim(min(pointsALL$latitude), max(pointsALL$latitude))

### Generate polynomial simulation data ###
datause <- genDGP(2500, "price")

# For better visualization
datause[datause$price>=260,]$price = 260
datause[datause$price<=180,]$price = 180

p_poly <- ggplot() +
  geom_point(data = datause, aes(x=longitude, y=latitude, color = price), size=0.8) +
  geom_path(aes(pointsALL$longitude, pointsALL$latitude), size=0.8) +
  scale_color_gradient2(low = 'dodgerblue3', midpoint=220, mid="lavender",high = 'firebrick3', limits=c(180, 260)) +
  labs(title="Simulated Housing Price (Polynomial)",
       x ="Longitude", y = "Latitude", color = "Price") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
  geom_point(data = pointsALL[c(30, 45, 60),], aes(x=longitude, y=latitude), shape = 15, color = "black", size = 2, show.legend = FALSE) +
  geom_point(data = pointsALL[45,], aes(x=longitude, y=latitude), shape = 1, color = "green", size = 4, stroke=1.5, show.legend = FALSE) +
  xlim(min(pointsALL$longitude), max(pointsALL$longitude)) +
  ylim(min(pointsALL$latitude), max(pointsALL$latitude))

### Generate WGAN simulation data ###
wgan_d <- read.csv("../../../Data/Keele & Titiunik 2015/WGAN_Epoch1k/price_cWGAN.csv")
rand <- sample(nrow(wgan_d), 2500, replace = F)
wgan_sample <- wgan_d[rand,]

# For better visualization
wgan_sample[wgan_sample$Y>=260,]$Y = 260
wgan_sample[wgan_sample$Y<=180,]$Y = 180

p_wgan <- ggplot() +
  geom_point(data = wgan_sample, aes(x=X2, y=X1, color = Y), size=0.8) +
  geom_path(aes(pointsALL$longitude, pointsALL$latitude), size=0.8) +
  scale_color_gradient2(low = 'dodgerblue3', midpoint=220, mid="lavender",high = 'firebrick3', limits=c(180, 260)) +
  labs(title="Simulated Price (WGAN)",
       x ="Longitude", y = "Latitude", color = "Price") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
  geom_point(data = pointsALL[c(30, 45, 60),], aes(x=longitude, y=latitude), shape = 15, color = "black", size = 2, show.legend = FALSE) +
  geom_point(data = pointsALL[45,], aes(x=longitude, y=latitude), shape = 1, color = "green", size = 4, stroke=1.5, show.legend = FALSE) +
  xlim(min(pointsALL$longitude), max(pointsALL$longitude)) +
  ylim(min(pointsALL$latitude), max(pointsALL$latitude))

fig <- arrangeGrob(p1, p_poly, p_wgan, ncol = 3, top = text_grob("Bivariate Design: Keele & Titiunik (2015)", face = "bold", size = 20))
ggsave("KT_price_vis.png", fig, width = 16, height = 6)
