######################################
# Figure B.7: Surface plots near the 
# boundary point for the polynomial DGP 
# and WGAN DGP, respectively (for 
# outcome=price).
######################################
rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  # set directory to the location of the current script
library(plotly)
source("../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Analysis/Local Poly/distance-functions.R")

#### Helper function for splitting trt & ctrl based on coordinates ####
correct <- function(X1, X2) {
  simdata <- as.data.frame(cbind("latitude"=X1, "longitude"=X2))
  simdata$treat <- NA

  for (i in 1:nrow(simdata)) { # Identify "region" based on longitude
    if (i %% 5e3 == 0) { # show progress
      print(paste0("Progress: ", 100*i/nrow(simdata), "%"))
    }

    closest_x <- which.min(abs(points[,"longitude"] - simdata$longitude[i]))

    # Region 1, 3, 5, 7: 1-26, 43-49, 66-73, 78-81
    # Longitude as a function of latitude
    if ((closest_x <= 26) | (closest_x >= 43 & closest_x <= 49) |
        (closest_x >= 66 & closest_x <= 73) | (closest_x >= 78 & closest_x <= 81)) {
      # Find closest point (index) in terms of latitude in this region
      blocks <- list(c(1:26), c(43:49), c(66:73), c(78:81))
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
        else{simdata$treat[i] <- 1}
      }
    }

    # Region 2, 4, 6, 8: 26-43, 49-66, 73-78, 81-89
    # Latitude as a function of longitude
    if ((closest_x >= 26 & closest_x <= 43) | (closest_x >= 49 & closest_x <= 66) |
        (closest_x >= 73 & closest_x <= 78) | (closest_x > 81 & closest_x <= 89)) {
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

  return(simdata$treat)
}

#### Polynomial DGP ####
## Gather data for adding a scatterplot of original price data
data_sample <- read.csv("../../Data/Keele & Titiunik 2015/Observed Samples/data_price.csv")
sample_X1 <- as.matrix(data_sample$X1, ncol=1)
sample_X2 <- as.matrix(data_sample$X2, ncol=1)
sample_Y <- as.matrix(data_sample$Y, ncol=1)
sample_W <- as.matrix(data_sample$W, ncol=1)

### Get focal point on the treatment boundary where we estimate treatment effect ###
points <- foreign::read.dbf("../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Data/BorderSegmentPoints_Project.dbf") # boundary data
points$latitude <- points$POINT_Y
points$longitude <- points$POINT_X
focal <- matrix(c(points$latitude[45],points$longitude[45]),nrow = 1)

# Polynomial DGP for price
data <- haven::read_dta("../../Data/Keele & Titiunik 2015/KeeleTitiunik2014-PA-replication-files/Data/Housing/NJ_House_Final.dta")
data <- data[complete.cases(data$latitude, data$longitude, data$treat, data[,c("sqft_price")]),]

X1 <- data$latitude
X2 <- data$longitude
mean_X1 <- mean(X1)
mean_X2 <- mean(X2)
price <- data$sqft_price

model.df <- data.frame(Y=price,
                       X1=X1,
                       X1.2=(X1 - mean_X1)^2,
                       X1.3=(X1 - mean_X1)^3,
                       X2=X2,
                       X2.2=(X2 - mean_X2)^2,
                       X2.3=(X2 - mean_X2)^3,
                       X12=(X1 - mean_X1)*(X2 - mean_X2),
                       W=data$treat)

cef <- lm(Y ~ W * ., data = model.df)

# First create a grid
latitude_plot <- seq(min(points$latitude), max(points$latitude), length.out=100)

longitude_plot <- seq(min(points$longitude), max(points$longitude), length.out=100)

# Then fill in E[Y|X] for each X in the neighborhood
signal_poly_T <- matrix(data=NA,
                      nrow=length(latitude_plot),
                      ncol=length(longitude_plot))

for (i in 1:nrow(signal_poly_T)){
  for (j in 1:ncol(signal_poly_T)){
    lat_temp <- latitude_plot[i]
    long_temp <- longitude_plot[j]
    
    if (correct(lat_temp, long_temp) == 1) {
      new_X <- data.frame(X1=lat_temp,
                          X1.2=(lat_temp - mean_X1)^2,
                          X1.3=(lat_temp - mean_X1)^3,
                          X2=long_temp,
                          X2.2=(long_temp - mean_X2)^2,
                          X2.3=(long_temp - mean_X2)^3,
                          X12=(lat_temp - mean_X1)*(long_temp - mean_X2),
                          W=1)
      
      signal_poly_T[i,j] <- predict(cef, new_X)
    } else {
      signal_poly_T[i,j] <- NA
    }
  }
}

signal_poly_C <- matrix(data=NA,
                        nrow=length(latitude_plot),
                        ncol=length(longitude_plot))

for (i in 1:nrow(signal_poly_C)){
  for (j in 1:ncol(signal_poly_C)){
    lat_temp <- latitude_plot[i]
    long_temp <- longitude_plot[j]
    
    if (correct(lat_temp, long_temp) == 0) {
      new_X <- data.frame(X1=lat_temp,
                          X1.2=(lat_temp - mean_X1)^2,
                          X1.3=(lat_temp - mean_X1)^3,
                          X2=long_temp,
                          X2.2=(long_temp - mean_X2)^2,
                          X2.3=(long_temp - mean_X2)^3,
                          X12=(lat_temp - mean_X1)*(long_temp - mean_X2),
                          W=0)
      
      signal_poly_C[i,j] <- predict(cef, new_X)
    } else {
      signal_poly_C[i,j] <- NA
    }
  }
}

focal_T <- predict(cef, data.frame(X1=focal[1],
                                   X1.2=(focal[1] - mean_X1)^2,
                                   X1.3=(focal[1] - mean_X1)^3,
                                   X2=focal[2],
                                   X2.2=(focal[2] - mean_X2)^2,
                                   X2.3=(focal[2] - mean_X2)^3,
                                   X12=(focal[1] - mean_X1)*(focal[2] - mean_X2),
                                   W=1))

focal_C <- predict(cef, data.frame(X1=focal[1],
                                   X1.2=(focal[1] - mean_X1)^2,
                                   X1.3=(focal[1] - mean_X1)^3,
                                   X2=focal[2],
                                   X2.2=(focal[2] - mean_X2)^2,
                                   X2.3=(focal[2] - mean_X2)^3,
                                   X12=(focal[1] - mean_X1)*(focal[2] - mean_X2),
                                   W=0))

# figure for the surface/signal/cef of the polynomial DGP
surface_poly <- plot_ly(x = longitude_plot, y = latitude_plot, z = signal_poly_T, scene='scene1', showscale=F) %>% 
  add_surface(cmin = 180, cmax = 295, colorscale = list(c(0,1),c("skyblue","pink")), opacity = 0.8) %>% 
  add_surface(signal_poly_C, cmin = 180, cmax = 295, colorscale = list(c(0,1),c("skyblue","pink")), opacity = 0.8, showscale = F, showscale=F) %>%
  add_trace(x=focal[2], y=focal[1], z=focal_C, 
            type = 'scatter3d', mode = 'markers', 
            name = 'focal_ctrl',
            marker = list(size = 5, color='black'), showlegend = F) %>%
  add_trace(x=focal[2], y=focal[1], z=focal_T,
            type = 'scatter3d', mode = 'markers',
            name = 'focal_trt',
            marker = list(size = 5, color='red'), showlegend = F) %>%
  add_trace(x=sample_X2[sample_W==1], # longitude
            y=sample_X1[sample_W==1], # latitude
            z=sample_Y[sample_W==1],
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 1.5, color='firebrick', opacity=0.5),
            showlegend=FALSE) %>%
  add_trace(x=sample_X2[sample_W==0], # longitude
            y=sample_X1[sample_W==0], # latitude
            z=sample_Y[sample_W==0],
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 1.5, color='black'),
            showlegend=FALSE, opacity=0.5)

# surface_poly

####  ----------------WGAN--------------- ####
### 1kw epochs
### for the WGAN DGP, for each running variable in the neighborhood latitude_plot x longitude_plot, we need to draw multiple Y (price) to simulate E[Y|X] since there's no closed-form formula.

#### Full grid with both trt & ctrl ####
# Window is the same as that of the boundary
latitude_plot <- seq(min(points$latitude), max(points$latitude), length.out=100)
latitude_plot <- round(latitude_plot, 5)

longitude_plot <- seq(min(points$longitude), max(points$longitude), length.out=100)
longitude_plot <- round(longitude_plot, 5)

focal <- round(focal,5)

### All combinations of latitudes & longitudes in grid, for WGAN input & no need to rerun 
# X_plot <- data.frame(cbind("X1"=rep(latitude_plot,
#                                     each=length(latitude_plot)),
#                            "X2"=rep(longitude_plot,
#                                     times=length(longitude_plot))))
# X_plot$W <- correct(X_plot$X1, X_plot$X2)
# 
# # Add focal to both trt & ctrl DGPs
# X_plot[nrow(X_plot)+1,] <- c(round(focal,5), 1)
# X_plot[nrow(X_plot)+1,] <- c(round(focal,5), 0)
# 
# 
# # Repeat each row in X_plot 2000 times to approximate E[Y|X]
# X_plot <- X_plot[rep(1:nrow(X_plot), each = 2000), ]
# write.csv(X_plot, "surface_plot_wgan_X_full.csv", row.names = F)


### For the WGAN DGP, for each running variable in the neighborhood latitude_plot x longitude_plot, we need to draw multiple Y (price) to simulate E[Y|X] since there's no closed-form formula.

### Match & compute conditional mean for every grid point, no need to rerun as the matched data is saved and uploaded as csv files for convenience
# d <- read.csv("../../Data/Keele & Titiunik 2015/WGAN_Epoch1k/plot_price.csv")
# dT <- d
# dT$Y[dT$W==0] <- NA
# dC <- d
# dC$Y[dC$W==1] <- NA
# 
# ## Trt
# signal_wgan_T <- matrix(data=NA,
#                         nrow=length(latitude_plot),
#                         ncol=length(longitude_plot))
# colnames(signal_wgan_T) <- longitude_plot
# rownames(signal_wgan_T) <- latitude_plot
# 
# for (i in 1:nrow(signal_wgan_T)){
#   print(i)
#   for (j in 1:ncol(signal_wgan_T)){
#     lat_temp <- latitude_plot[i]  # each row corresponds to 1 latitude
#     long_temp <- longitude_plot[j]  # each col corresponds to 1 longitude
# 
#     signal_wgan_T[i,j] <- mean(dT[dT$X1==lat_temp & dT$X2==long_temp,]$Y, na.rm=T)
#     print(signal_wgan_T[i,j])
#   }
# }
# 
# focal_T <- mean(dT[dT$X1==focal[1,1] & dT$X2==focal[1,2],]$Y, na.rm=T)
# 
# ## Ctrl
# signal_wgan_C <- matrix(data=NA,
#                         nrow=length(latitude_plot),
#                         ncol=length(longitude_plot))
# colnames(signal_wgan_C) <- longitude_plot
# rownames(signal_wgan_C) <- latitude_plot
# 
# for (i in 1:nrow(signal_wgan_C)){
#   print(i)
#   for (j in 1:ncol(signal_wgan_C)){
#     lat_temp <- latitude_plot[i]  # each row corresponds to 1 latitude
#     long_temp <- longitude_plot[j]  # each col corresponds to 1 longitude
# 
#     signal_wgan_C[i,j] <- mean(dC[dC$X1==lat_temp & dC$X2==long_temp,]$Y, na.rm=T)
#     print(signal_wgan_C[i,j])
#   }
# }

signal_wgan_T <- read.csv("Plot Data/signal_wgan_T_P.csv")
signal_wgan_C <- read.csv("Plot Data/signal_wgan_C_P.csv")
signal_wgan_T <- as.matrix(signal_wgan_T)
signal_wgan_C <- as.matrix(signal_wgan_C)
signal_wgan_T <- signal_wgan_T[,-1]
signal_wgan_C <- signal_wgan_C[,-1]

focal_T <- 237.22789435
focal_C <- 234.212192465

#### Plot surface & scatterplot ####
surface_wgan1k <- plot_ly(x = longitude_plot, y = latitude_plot, z = signal_wgan_T, scene='scene2', showscale=F) %>% 
  add_surface(cmin = 180, cmax = 295, colorscale = list(c(0,1),c("skyblue","pink")), opacity = 0.8, showscale = T) %>% 
  add_surface(signal_wgan_C, cmin = 180, cmax = 295, colorscale = list(c(0,1),c("skyblue","pink")), opacity = 0.8, showscale=F) %>%
  add_trace(x=focal[2], y=focal[1], z=focal_C, 
            type = 'scatter3d', mode = 'markers', 
            name = 'focal_ctrl',
            marker = list(size = 5, color='black'), showlegend = T) %>%
  add_trace(x=focal[2], y=focal[1], z=focal_T,
            type = 'scatter3d', mode = 'markers',
            name = 'focal_trt',
            marker = list(size = 5, color='red'), showlegend = T) %>%
  add_trace(x=sample_X2[sample_W==1], # longitude
            y=sample_X1[sample_W==1], # latitude
            z=sample_Y[sample_W==1],
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 1.5, color='firebrick'), opacity=0.5, showlegend = F) %>%
  add_trace(x=sample_X2[sample_W==0], # longitude
            y=sample_X1[sample_W==0], # latitude
            z=sample_Y[sample_W==0],
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 1.5, color='black'), opacity=0.5, showlegend = F)

### Organize plots in main text
fig_main <- subplot(surface_poly, surface_wgan1k)
fig_main <- fig_main %>% layout(scene = list(text="Polynomial DGP", xaxis=list(title="Longitude"), 
                                             yaxis=list(title="Latitude"), 
                                             zaxis=list(title="E[price | X]",
                                                        range = list(180, 300)),
                                             aspectmode='cube'),
                                scene2 = list(text="WGAN DGP, 1k Epochs", xaxis=list(title="Longitude"), 
                                              yaxis=list(title="Latitude"), 
                                              zaxis=list(title="E[price | X]",
                                                         range = list(180, 300)),
                                              aspectmode='cube'))
fig_main



### 2k epochs
signal_wgan_T <- read.csv("Plot Data/signal_wgan_T_P_2k.csv")
signal_wgan_C <- read.csv("Plot Data/signal_wgan_C_P_2k.csv")
signal_wgan_T <- as.matrix(signal_wgan_T)
signal_wgan_C <- as.matrix(signal_wgan_C)
signal_wgan_T <- signal_wgan_T[,-1]
signal_wgan_C <- signal_wgan_C[,-1]

focal_T <- 216.72418724
focal_C <- 226.513909585

#### Plot surface & scatterplot ####
surface_wgan2k <- plot_ly(x = longitude_plot, y = latitude_plot, z = signal_wgan_T, scene = "scene1") %>% 
  add_surface(cmin = 180, cmax = 295, colorscale = list(c(0,1),c("skyblue","pink")), opacity = 0.8, showscale = FALSE) %>% 
  add_surface(signal_wgan_C, cmin = 180, cmax = 295, colorscale = list(c(0,1),c("skyblue","pink")), opacity = 0.8, showscale = FALSE) %>%
  add_trace(x=focal[2], y=focal[1], z=focal_C, 
            type = 'scatter3d', mode = 'markers', 
            name = 'focal_ctrl',
            marker = list(size = 5, color='black'),
            showlegend=FALSE) %>%
  add_trace(x=focal[2], y=focal[1], z=focal_T,
            type = 'scatter3d', mode = 'markers',
            name = 'focal_trt',
            marker = list(size = 5, color='red'),
            showlegend=FALSE) %>%
  add_trace(x=sample_X2[sample_W==1], # longitude
            y=sample_X1[sample_W==1], # latitude
            z=sample_Y[sample_W==1],
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 1.5, color='firebrick'),
            showlegend=FALSE, opacity=0.5) %>%
  add_trace(x=sample_X2[sample_W==0], # longitude
            y=sample_X1[sample_W==0], # latitude
            z=sample_Y[sample_W==0],
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 1.5, color='black'),
            showlegend=FALSE, opacity=0.5) 

# surface_wgan2k

### 3kw epochs
signal_wgan_T <- read.csv("Plot Data/signal_wgan_T_P_3k.csv")
signal_wgan_C <- read.csv("Plot Data/signal_wgan_C_P_3k.csv")
signal_wgan_T <- as.matrix(signal_wgan_T)
signal_wgan_C <- as.matrix(signal_wgan_C)
signal_wgan_T <- signal_wgan_T[,-1]
signal_wgan_C <- signal_wgan_C[,-1]

focal_T <- 204.341663355
focal_C <- 232.608909035

#### Plot surface & scatterplot ####
surface_wgan3k <- plot_ly(x = longitude_plot, y = latitude_plot, z = signal_wgan_T, scene='scene2', showscale = F) %>% 
  add_surface(cmin = 180, cmax = 295, colorscale = list(c(0,1),c("skyblue","pink")), opacity = 0.8, showscale = FALSE) %>% 
  add_surface(signal_wgan_C, cmin = 180, cmax = 295, colorscale = list(c(0,1),c("skyblue","pink")), opacity = 0.8) %>%
  add_trace(x=focal[2], y=focal[1], z=focal_C, 
            type = 'scatter3d', mode = 'markers', 
            name = 'focal_ctrl',
            marker = list(size = 5, color='black'), showlegend = F) %>%
  add_trace(x=focal[2], y=focal[1], z=focal_T,
            type = 'scatter3d', mode = 'markers',
            name = 'focal_trt',
            marker = list(size = 5, color='red'), showlegend = F) %>%
  add_trace(x=sample_X2[sample_W==1], # longitude
            y=sample_X1[sample_W==1], # latitude
            z=sample_Y[sample_W==1],
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 1.5, color='firebrick'),
            showlegend=FALSE, opacity=0.5) %>%
  add_trace(x=sample_X2[sample_W==0], # longitude
            y=sample_X1[sample_W==0], # latitude
            z=sample_Y[sample_W==0],
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 1.5, color='black'),
            showlegend=FALSE, opacity=0.5)

# surface_wgan3k

### 4k epochs
signal_wgan_T <- read.csv("Plot Data/signal_wgan_T_P_4k.csv")
signal_wgan_C <- read.csv("Plot Data/signal_wgan_C_P_4k.csv")
signal_wgan_T <- as.matrix(signal_wgan_T)
signal_wgan_C <- as.matrix(signal_wgan_C)
signal_wgan_T <- signal_wgan_T[,-1]
signal_wgan_C <- signal_wgan_C[,-1]

focal_T <- 204.676572845
focal_C <- 235.351111965

#### Plot surface & scatterplot ####
surface_wgan4k <- plot_ly(x = longitude_plot, y = latitude_plot, z = signal_wgan_T, scene = "scene3") %>% 
  add_surface(cmin = 180, cmax = 295, colorscale = list(c(0,1),c("skyblue","pink")), opacity = 0.8, showscale = FALSE) %>% 
  add_surface(signal_wgan_C, cmin = 180, cmax = 295, colorscale = list(c(0,1),c("skyblue","pink")), opacity = 0.8) %>%
  add_trace(x=focal[2], y=focal[1], z=focal_C, 
            type = 'scatter3d', mode = 'markers', 
            name = 'focal_ctrl',
            marker = list(size = 5, color='black')) %>%
  add_trace(x=focal[2], y=focal[1], z=focal_T,
            type = 'scatter3d', mode = 'markers',
            name = 'focal_trt',
            marker = list(size = 5, color='red')) %>%
  add_trace(x=sample_X2[sample_W==1], # longitude
            y=sample_X1[sample_W==1], # latitude
            z=sample_Y[sample_W==1],
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 1.5, color='firebrick'),
            showlegend=FALSE, opacity=0.5) %>%
  add_trace(x=sample_X2[sample_W==0], # longitude
            y=sample_X1[sample_W==0], # latitude
            z=sample_Y[sample_W==0],
            type = 'scatter3d', mode = 'markers',
            marker = list(size = 1.5, color='black'),
            showlegend=FALSE, opacity=0.5) 

# surface_wgan4k

### Organize plots in the appendix
fig_appendix <- subplot(surface_wgan2k, surface_wgan3k, surface_wgan4k)
fig_appendix <- fig_appendix %>% layout(scene = list(text="WGAN DGP, 2k Epochs", xaxis=list(title="Longitude"), 
                                    yaxis=list(title="Latitude"), 
                                    zaxis=list(title="E[price | X]",
                                               range = list(180, 300)),
                                    aspectmode='cube'),
                      scene2 = list(text="WGAN DGP, 3k Epochs", xaxis=list(title="Longitude"), 
                                    yaxis=list(title="Latitude"), 
                                    zaxis=list(title="E[price | X]",
                                               range = list(180, 300)),
                                    aspectmode='cube'),
                      scene3 = list(text="WGAN DGP, 4k Epochs", xaxis=list(title="Longitude"), 
                                    yaxis=list(title="Latitude"), 
                                    zaxis=list(title="E[price | X]",
                                               range = list(180, 300)),
                                    aspectmode='cube'))
fig_appendix
