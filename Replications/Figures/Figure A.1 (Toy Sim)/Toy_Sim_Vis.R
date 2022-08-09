######################################
# Figure A.1: A Toy Simulation
######################################
library(plotly)
set.seed(6) # for reproducible results

# 'genDGP' is a function that generates a sample from our specified simulation GDP
# Arguments:
#           size: sample size
#           tau: treatment effect in the linear case; default is 0.3
#           beta: coeffcients of running var X on the covariates Z
#           delta: coefficients of covaraites Z on the outcome Y
genDGP <- function(n, tau=0.3, 
                   beta=matrix(c(-0.14677601, 2.82716002,
                                 1.92389438, 0.15973524, 
                                 0.08320535, -0.70270603,
                                 1.12296211, 0.19537758,
                                 1.29922865, -1.84057169,
                                 1.21369502, 1.08040631,
                                 1.33298827, -0.70907341,
                                 0.82408303, 1.8264555,
                                 -1.1041983, 2.3728503, 
                                 1.8979988, -0.9914020,
                                 -1.2240557, 1.1042126,
                                 -0.5441380, 0.5972449, 
                                 -1.5626238, -1.2025887,
                                 -0.2091987, -1.8234970,
                                 -0.2299455, 1.2882160), ncol=15), 
                   delta=c(1.215260987, -2.784184991, 
                           2.176880660, 3.014243530, 
                           -1.256691145, 5.063790854, 
                           -0.980040716, 0.009260695, 
                           -0.001064875, -0.001994080, 
                           0.005007084, -0.003165896, 
                           -0.002410976, -0.002257493,  
                           0.002370995)){
  
  # ==== Dataframe for plotting ====
  # Running variables X_i
  X1 <- runif(n, -1, 1)
  X2 <- runif(n, -1, 1)
  X3 <- runif(n, -1, 1)
  
  # Treatment indicator D_i
  D1 <- ifelse(X1 >= 0, 1, 0) # univariate running var
  D2 <- ifelse(X1 >= 0 & X2 >= 0, 1, 0) # multivariate running var
  D3 <- ifelse(X1 >= 0 & X2 >= 0 & X3 >= 0, 1, 0) # multivariate running var
  
  # Covariates Z
  Z1 <- matrix(runif(n*15)+rep(beta[1,],n)*rep(X1, each=15), ncol=15)
  Z2 <- matrix(runif(n*15)+rep(beta[1,],n)*rep(X1, each=15)+rep(beta[2,],n)*rep(X2, each=15), ncol=15)
  
  # Without covariates, outcome Y ~ Bernoulli(0.2 + tau*D_i)
  Y1 <- rbinom(n, 1, 0.2+tau*D1)
  Y2 <- rbinom(n, 1, 0.2+tau*D2)
  Y3 <- rbinom(n, 1, 0.2+tau*D3)
  
  # With covariates, outcome Y ~ Bernoulli(sigmoid(linear function of Zs and Xs))
  Y1_sigmoid <- rbinom(n,1,1/(1+exp(-(Z1 %*% matrix(delta, ncol=1) + D1 - 0.05*X1))))
  Y2_sigmoid <- rbinom(n,1,1/(1+exp(-(Z2 %*% matrix(delta, ncol=1) + D2 - 0.05*X1 + 0.02*X2))))
  
  # Big dataframe for all cases
  colnames(Z1) <- c("Z1_1","Z1_2","Z1_3","Z1_4","Z1_5","Z1_6","Z1_7","Z1_8","Z1_9","Z1_10","Z1_11","Z1_12","Z1_13","Z1_14","Z1_15")
  colnames(Z2) <- c("Z2_1","Z2_2","Z2_3","Z2_4","Z2_5","Z2_6","Z2_7","Z2_8","Z2_9","Z2_10","Z2_11","Z2_12","Z2_13","Z2_14","Z2_15")
  mydata <- as.data.frame(cbind("Y1"= Y1, "Y2" = Y2, "Y3" = Y3, 
                                "Y1_sigmoid" = Y1_sigmoid, 
                                "Y2_sigmoid" = Y2_sigmoid, 
                                "X1"=X1, "X2"=X2, "X3"=X3, 
                                "D1"=D1, "D2"=D2, "D3"=D3, Z1, Z2))
  return(mydata)
}


# ============== Univariate RV ================
mydata1 <- genDGP(500000)
uni_nocov <-RDData(mydata1[c("Y1", "X1")], cutoff = 0)
p1 <- plot_RDscatter(uni_nocov, avg=1000, xlab = "Running variable", ylab = "Outcome") + theme_bw() + geom_point(size=1, color="deepskyblue")

# ============== Bivariate RV ================
breaks2 <- seq(-1,1,by=0.1)
ag2 <- aggregate(mydata1$Y2, FUN=mean,by=list(X1=cut(mydata1$X1, breaks2, include.lowest=T), X2=cut(mydata1$X2, breaks2, include.lowest=T)))
ag2$midX1 <- as.numeric(sapply(strsplit(as.character(ag2$X1), "\\(|\\)|\\[|\\]|\\,"), "[[", 2))+0.05
ag2$midX2 <- as.numeric(sapply(strsplit(as.character(ag2$X2), "\\(|\\)|\\[|\\]|\\,"), "[[", 2))+0.05
p2 <- ggplot(ag2, aes(midX1,midX2))+geom_raster(aes(fill=x))+scale_fill_gradient2(
  low = "deepskyblue", mid="white", high = "red", midpoint=mean(ag2$x), name="Outcome") +
  theme_bw()+xlab("Running variable 1")+ylab("Running variable 2")

grid.arrange(p1, p2, nrow = 1, ncol = 2)

# ============== Trivariate RV ================
breaks3 <- seq(-1,1,by=0.1)
ag3 <- aggregate(mydata1$Y3, FUN=mean,by=list(X1=cut(mydata1$X1, breaks3, include.lowest=T), 
                                              X2=cut(mydata1$X2, breaks3, include.lowest=T), 
                                              X3=cut(mydata1$X3, breaks3, include.lowest=T)))
ag3$midX1 <- as.numeric(sapply(strsplit(as.character(ag3$X1), "\\(|\\)|\\[|\\]|\\,"), "[[", 2))+0.05
ag3$midX2 <- as.numeric(sapply(strsplit(as.character(ag3$X2), "\\(|\\)|\\[|\\]|\\,"), "[[", 2))+0.05
ag3$midX3 <- as.numeric(sapply(strsplit(as.character(ag3$X3), "\\(|\\)|\\[|\\]|\\,"), "[[", 2))+0.05

colnames(ag3)[4] <- "Outcome"
axx <- list(title = "Running variable 1", titlefont = list(size = 10))
axy <- list(title = "Running variable 2", titlefont = list(size = 10))
axz <- list(title = "Running variable 3", titlefont = list(size = 10))
p3 <- plot_ly(ag3, x = ~midX1, y = ~midX2, z = ~midX3) %>%
  add_markers(color = ~Outcome, colors = c("deepskyblue", "white", "red")) %>%
  layout(scene = list(xaxis=axx,yaxis=axy,zaxis=axz))
p3
