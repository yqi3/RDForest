######################################
# Figure 2: Bivariate score using
# data from Keele & Titiunik (2015)
######################################
library(haven)
library(foreign)
library(tidyverse)
library(dplyr)

set.seed(888)

### Read sample data ###
data <- read_dta("KeeleTitiunik2014-PA-replication-files/Data/Housing/NJ_House_Final.dta")
data <- data[complete.cases(data$latitude, data$longitude, data$treat, data[,"sqft_price"]),]
points <- read.dbf("KeeleTitiunik2014-PA-replication-files/Data/BorderSegmentPoints_Project.dbf") # boundary data
points$latitude= points$POINT_Y
points$longitude= points$POINT_X

X1 <- data$latitude
X2 <- data$longitude
price <- data$sqft_price

model.df <- data.frame(Y=price, X1=X1, X1.2=(X1 - mean(X1))^2, X1.3=(X1 - mean(X1))^3,
                       X2=X2, X2.2=(X2 - mean(X2))^2, X2.3=(X2 - mean(X2))^3,
                       X12=(X1 - mean(X1))*(X2 - mean(X2)),
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
  geom_path(aes(points$longitude, points$latitude), size=0.8) +
  scale_color_gradient2(low = 'dodgerblue3', midpoint=220, mid="lavender",high = 'firebrick3', limits=c(180, 260)) +
  labs(title="Actual Housing Price Distribution",
       x ="Longitude", y = "Latitude", color = "Price") + 
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(data = points[c(30, 45, 60),], aes(x=longitude, y=latitude), shape = 15, color = "black", size = 2, show.legend = FALSE) +
  geom_point(data = points[45,], aes(x=longitude, y=latitude), shape = 1, color = "green", size = 4, stroke=1.5, show.legend = FALSE)


p2 <- ggplot() +
  geom_point(data = pred_plot, aes(x=X2, y=X1, color = pred_Y), size=0.8) +
  geom_path(aes(points$longitude, points$latitude), size=0.8) +
  scale_color_gradient2(low = 'dodgerblue3', midpoint=220, mid="lavender",high = 'firebrick3', limits=c(180, 260)) +
  labs(title="Predicted Housing Price Distribution",
       x ="Longitude", y = "Latitude", color = "Price") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))+ 
  geom_point(data = points[c(30, 45, 60),], aes(x=longitude, y=latitude), shape = 15, color = "black", size = 2, show.legend = FALSE) +
  geom_point(data = points[45,], aes(x=longitude, y=latitude), shape = 1, color = "green", size = 4, stroke=1.5, show.legend = FALSE)


genDGP <- function(n){

  # sample from the original data, without restricting range  
  rand <- sample_n(data, n, replace=TRUE)
  latitude <- rand$latitude+rnorm(n, 0, sd=0.01)
  longitude <- rand$longitude+rnorm(n, 0, sd=0.01)
  
  simdata <- as.data.frame(cbind("latitude"=latitude, "longitude"=longitude))
  simdata$treat <- NA
  
  for (i in 1:n) {# Identify "region" based on longitude
    closest_x <- which.min(abs(points[,"longitude"] - simdata$longitude[i]))
    
    # Region 1, 3, 5, 7: 1-26, 44-49, 66-73, 79-81
    # Longitude as a function of latitude
    if ((closest_x <= 26) | (closest_x >= 44 & closest_x <= 49) | 
        (closest_x >= 66 & closest_x <= 73) | (closest_x >= 79 & closest_x <= 81)) {
      # Find closest point (index) in terms of latitude in this region
      blocks <- list(c(1:26), c(44:49), c(66:73), c(79:81))
      for (j in 1:length(blocks)) {
        if (closest_x %in% blocks[[j]]) {region <- blocks[[j]]}}
      index1 <- which.min(abs(points[region,"latitude"] - simdata$latitude[i])) + region[1]-1

      if (points$latitude[index1] < simdata$latitude[i]){ 
        # simulation point is higher
        # index should decrease by 1 so latitude is bigger, fit a line between the 2 indices
        
        if (index1 == 1) { # Boundary situation, simply compare longitudes
          if (simdata$longitude[i] >= points$longitude[1]) {
            simdata$treat[i] <- 0}
          else {simdata$treat[i] <- 1}
        }
        
        else { # Not boundary, can subtract 1 to find the other index and fit a line
          index2 <- index1-1
          
          # Fitting a line (note: this is equivalent to "x as a function of y")
          slope_original <- (points$latitude[index1] - points$latitude[index2])/(points$longitude[index1] - points$longitude[index2])
          slope <- 1/slope_original
          intercept <- -(points$latitude[index1] - points$longitude[index1]*slope_original)/slope_original
          
          if (simdata$longitude[i] > simdata$latitude[i]*slope+intercept){
            # Simulation point is above the fitted line, not treated
            simdata$treat[i] <- 0
          }
          else{simdata$treat[i] <- 1}}
      }
      
      else {# simulation point is lower
        # index should increase by 1 so latitude is smaller, fit a line between the 2 indices
        index2 <- index1+1    # No boundary situation here
        
        # Fitting a line (note: this is equivalent to "x as a function of y")  
        slope_original <- (points$latitude[index1] - points$latitude[index2])/(points$longitude[index1] - points$longitude[index2])
        slope <- 1/slope_original
        intercept <- -(points$latitude[index1] - points$longitude[index1]*slope_original)/slope_original
        
        if (simdata$longitude[i] > simdata$latitude[i]*slope+intercept){
          simdata$treat[i] <- 0}
        else{simdata$treat[i] <- 1}}
    }
    
    # Region 2, 4, 6, 8: 26-44, 49-66, 73-79, 81-89
    # Latitude as a function of longitude
    if ((closest_x >= 26 & closest_x <= 44) | (closest_x >= 49 & closest_x <= 66) | 
        (closest_x >= 73 & closest_x <= 79) | (closest_x > 81 & closest_x <= 89)) {
      index1 <- closest_x
      
      if (points$longitude[index1] > simdata$longitude[i]){
        index2 <- index1-1
        slope <- (points$latitude[index1] - points$latitude[index2])/(points$longitude[index1] - points$longitude[index2])
        intercept <- points$latitude[index1] - points$longitude[index1]*slope
        
        if (simdata$latitude[i] > simdata$longitude[i]*slope+intercept) {
          simdata$treat[i] <- 0}
        else {simdata$treat[i] <- 1}}
      
      else {
        if (index1 == 89) { # Boundary situation here
          if (points$latitude[index1] < simdata$latitude[i]) {
            simdata$treat[i] <- 0}
          else {simdata$treat[i] <- 1}}
        
        else {
          index2 <- index1+1
          slope <- (points$latitude[index1] - points$latitude[index2])/(points$longitude[index1] - points$longitude[index2])
          intercept <- points$latitude[index1] - points$longitude[index1]*slope
          
          if (simdata$latitude[i] > simdata$longitude[i]*slope+intercept) {
            simdata$treat[i] <- 0}
          else {simdata$treat[i] <- 1}}
      }
    }
  }
  
  X1 <- simdata$latitude
  X2 <- simdata$longitude
  
  model.df <- data.frame(X1=X1, X1.2=(X1 - mean(X1))^2, X1.3=(X1 - mean(X1))^3,
                         X2=X2, X2.2=(X2 - mean(X2))^2, X2.3=(X2 - mean(X2))^3,
                         X12=(X1 - mean(X1))*(X2 - mean(X2)), 
                         W=simdata$treat)
  simdata$price <- predict(cef, model.df)+rnorm(n,mean=0,sd=sigma_y)
  return(simdata)
}

datause <- genDGP(2500)

# For better visualization
datause[datause$price>=260,]$price = 260
datause[datause$price<=180,]$price = 180

p3 <- ggplot() +
  geom_point(data = datause, aes(x=longitude, y=latitude, color = price), size=0.8) +
  geom_path(aes(points$longitude, points$latitude), size=0.8) +
  scale_color_gradient2(low = 'dodgerblue3', midpoint=220, mid="lavender",high = 'firebrick3', limits=c(180, 260)) +
  labs(title="Simulated Housing Price Distribution",
       x ="Longitude", y = "Latitude", color = "Price") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(data = points[c(30, 45, 60),], aes(x=longitude, y=latitude), shape = 15, color = "black", size = 2, show.legend = FALSE) +
  geom_point(data = points[45,], aes(x=longitude, y=latitude), shape = 1, color = "green", size = 4, stroke=1.5, show.legend = FALSE)

grid.arrange(p1, p2, p3, nrow = 1, ncol = 3)