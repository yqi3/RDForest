######################################
# Figure 3: Bivariate score using
# data from Keele & Titiunik (2015)
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
data <- data[complete.cases(data$latitude, data$longitude, data$treat, data[,"e2008g"]),]
pointsALL <- read.dbf("../../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Data/BorderSegmentPoints_Project.dbf") # boundary data
pointsALL$latitude= pointsALL$POINT_Y
pointsALL$longitude= pointsALL$POINT_X
source("../../../Helper Fns/KT_Helper_Fns.R")

X1 <- data$latitude
X2 <- data$longitude
mean_X1 <- mean(X1)
mean_X2 <- mean(X2)
e2008g <- data$e2008g

model.df <- data.frame(Y=e2008g, X1=X1, X1.2=(X1 - mean_X1)^2, X1.3=(X1 - mean_X1)^3,
                       X2=X2, X2.2=(X2 - mean_X2)^2, X2.3=(X2 - mean_X2)^3,
                       X12=(X1 - mean_X1)*(X2 - mean_X2),
                       W=data$treat)

cef <- glm(Y ~ W * ., family=binomial(link='logit'),data=model.df)

sigma_y <- sqrt(summary(cef)$sigma^2)

df_plot <- data.frame(e2008g, X1, X2)

df_plot$e2008g <- as.factor(df_plot$e2008g)
p1 <- ggplot() + 
  geom_point(data = df_plot, aes(x=X2, y=X1, color=e2008g), size=0.8) +
  geom_path(aes(pointsALL$longitude, pointsALL$latitude), linewidth=0.8) +
  scale_color_manual(values=c('dodgerblue3', 'firebrick3')) +
  labs(title="Actual Turnout Distribution",
       x ="Longitude", y = "Latitude", color = "Turnout") + 
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
  geom_point(data = pointsALL[c(30, 45, 60),], aes(x=longitude, y=latitude), shape = 15, color = "black", size = 2, show.legend = FALSE) +
  geom_point(data = pointsALL[45,], aes(x=longitude, y=latitude), shape = 1, color = "green", size = 4, stroke=1.5, show.legend = FALSE) +
  xlim(min(pointsALL$longitude), max(pointsALL$longitude)) +
  ylim(min(pointsALL$latitude), max(pointsALL$latitude))

### Generate polynomial simulation data ###
datause <- genDGP(2500, "turnout")
datause$e2008g <- as.factor(datause$e2008g)

p_poly <- ggplot() +
  geom_point(data = datause, aes(x=longitude, y=latitude, color = e2008g), size=0.8) +
  geom_path(aes(pointsALL$longitude, pointsALL$latitude), linewidth=0.8) +
  scale_color_manual(values=c('dodgerblue3', 'firebrick3')) +
  labs(title="Simulated Turnout (Polynomial)",
       x ="Longitude", y = "Latitude", color = "Turnout") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
  geom_point(data = pointsALL[c(30, 45, 60),], aes(x=longitude, y=latitude), shape = 15, color = "black", size = 2, show.legend = FALSE) +
  geom_point(data = pointsALL[45,], aes(x=longitude, y=latitude), shape = 1, color = "green", size = 4, stroke=1.5, show.legend = FALSE) +
  xlim(min(pointsALL$longitude), max(pointsALL$longitude)) +
  ylim(min(pointsALL$latitude), max(pointsALL$latitude))

### Generate WGAN simulation data ###
wgan_d <- read.csv("../../../Data/Keele & Titiunik 2015/WGAN_Epoch1k/turnout_cWGAN.csv")
rand <- sample(nrow(wgan_d), 2500, replace = F)
wgan_sample <- wgan_d[rand,]
wgan_sample$Y <- as.factor(wgan_sample$Y)

p_wgan <- ggplot() +
  geom_point(data = wgan_sample, aes(x=X2, y=X1, color = Y), size=0.8) +
  geom_path(aes(pointsALL$longitude, pointsALL$latitude), linewidth=0.8) +
  scale_color_manual(values=c('dodgerblue3', 'firebrick3')) +
  labs(title="Simulated Turnout (WGAN)",
       x ="Longitude", y = "Latitude", color = "Turnout") +
  theme_bw(base_size = 15) +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") + 
  geom_point(data = pointsALL[c(30, 45, 60),], aes(x=longitude, y=latitude), shape = 15, color = "black", size = 2, show.legend = FALSE) +
  geom_point(data = pointsALL[45,], aes(x=longitude, y=latitude), shape = 1, color = "green", size = 4, stroke=1.5, show.legend = FALSE) +
  xlim(min(pointsALL$longitude), max(pointsALL$longitude)) +
  ylim(min(pointsALL$latitude), max(pointsALL$latitude))

fig <- arrangeGrob(p1, p_poly, p_wgan, ncol = 3, top = text_grob("Bivariate Design: Keele & Titiunik (2015)", face = "bold", size = 20))
ggsave("KT_turnout_vis.png", fig, width = 16, height = 6)

