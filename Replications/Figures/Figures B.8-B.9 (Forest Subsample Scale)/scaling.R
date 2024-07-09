######################################
# Figures B.8-B.9: Comparisons of different
# choices of the scaling constant using the
# polynomial & WGAN DGPs
######################################

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)
library(gridExtra)
library(directlabels)

poly_lee <- read.csv("Scale Data/poly-rf-lee-scale.csv")
poly_lee$scale[poly_lee$scale=="default"] <- 0.1

poly_T <- read.csv("Scale Data/poly-rf-T-scale.csv")
poly_T$scale[poly_T$scale=="default"] <- 0.1

poly_P <- read.csv("Scale Data/poly-rf-P-scale.csv")
poly_P$scale[poly_P$scale=="default"] <- 0.1

poly_A <- read.csv("Scale Data/poly-rf-A-scale.csv")
poly_A$scale[poly_A$scale=="default"] <- 0.1

wgan_lee <- read.csv("Scale Data/wgan-rf-lee-scale.csv")
wgan_lee$scale[wgan_lee$scale=="default"] <- 0.1

wgan_T <- read.csv("Scale Data/wgan-rf-T-scale.csv")
wgan_T$scale[wgan_T$scale=="default"] <- 0.1

wgan_P <- read.csv("Scale Data/wgan-rf-P-scale.csv")
wgan_P$scale[wgan_P$scale=="default"] <- 0.1

wgan_A <- read.csv("Scale Data/wgan-rf-A-scale.csv")
wgan_A$scale[wgan_A$scale=="default"] <- 0.1

#### Create plots ####
### Polynomial
df_poly_lee <- data.frame("Scale"=c(rep(0.1,4), rep(0.2,4), rep(0.3,4), rep(0.4,4), rep(0.5,4)), "Name"=rep(c("Bias^2 ","Variance ","MSE "," 95% CR"),5), "Value"=rep(NA,20))

for (scale in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
  df_poly_lee$Value[df_poly_lee$Scale==scale] <- c(mean(poly_lee$rf.bias[poly_lee$scale==scale])^2, var(poly_lee$rf[poly_lee$scale==scale]), mean(poly_lee$rf.bias[poly_lee$scale==scale])^2+var(poly_lee$rf[poly_lee$scale==scale]), (mean(poly_lee$rf.cr[poly_lee$scale==scale])-0.73)/100)
}

p_poly_lee <- ggplot(data=df_poly_lee[df_poly_lee$Name!=" 95% CR" & df_poly_lee$Name!="Variance ",], aes(x=Scale, y=Value, group=Name)) +
  geom_line(aes(color=Name), size=1) +
  geom_point(aes(color=Name), cex=2) +
  geom_dl(aes(label = Name, color=Name), method = list(dl.combine("first.points"), cex = 1)) +
  scale_y_continuous(name = "Squared bias, variance, MSE", sec.axis = sec_axis(~.*100+0.73, name="95% CR")) +
  scale_color_discrete(guide = "none") + 
  scale_x_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5), labels=c("Default","0.2","0.3","0.4","0.5"), limits = c(0.01,0.58)) +
  theme_bw(base_size = 15) +
  geom_line(data=df_poly_lee[df_poly_lee$Name=="Variance ",],
            aes(color=Name), size=1) +
  geom_point(aes(color=Name), data=df_poly_lee[df_poly_lee$Name=="Variance ",], cex=2) +
  geom_dl(data=df_poly_lee[df_poly_lee$Name=="Variance ",],
          aes(label = Name, color = Name),
          method=list(dl.combine("first.points"), cex = 1, vjust=1.2)) +
  geom_line(data=df_poly_lee[df_poly_lee$Name==" 95% CR",], 
            aes(color=Name),
            linetype='dashed', 
            col = 'black', size=1)+
  geom_point(data=df_poly_lee[df_poly_lee$Name==" 95% CR",], cex=2)+
  geom_dl(data=df_poly_lee[df_poly_lee$Name==" 95% CR",],
          aes(label = Name), 
          method=list(dl.combine("last.points"), color="black", cex = 1)) +
  ggtitle("Lee (2008): Polynomial")

p_poly_lee


df_poly_T <- data.frame("Scale"=c(rep(0.1,4), rep(0.2,4), rep(0.3,4), rep(0.4,4), rep(0.5,4)), "Name"=rep(c("Bias^2 ","Variance ","MSE "," 95% CR"),5), "Value"=rep(NA,20))

for (scale in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
  df_poly_T$Value[df_poly_T$Scale==scale] <- c(mean(poly_T$rf.bias[poly_T$scale==scale])^2, var(poly_T$rf[poly_T$scale==scale]), mean(poly_T$rf.bias[poly_T$scale==scale])^2+var(poly_T$rf[poly_T$scale==scale]), (mean(poly_T$rf.cr[poly_T$scale==scale])-0.85)/12)
}

p_poly_T <- ggplot(data=df_poly_T[df_poly_T$Name!=" 95% CR" & df_poly_T$Name!="Variance ",], aes(x=Scale, y=Value, group=Name)) +
  geom_line(aes(color=Name), size=1) +
  geom_point(aes(color=Name), cex=2) +
  geom_dl(aes(label = Name, color=Name), method = list(dl.combine("first.points"), cex = 1)) +
  scale_y_continuous(name = "Squared bias, variance, MSE", sec.axis = sec_axis(~.*12+0.85, name="95% CR")) +
  scale_color_discrete(guide = "none") + 
  scale_x_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5), labels=c("Default","0.2","0.3","0.4","0.5"), limits = c(0.01,0.58)) +
  theme_bw(base_size = 15) +
  geom_line(data=df_poly_T[df_poly_T$Name=="Variance ",],
            aes(color=Name), size=1) +
  geom_point(aes(color=Name), data=df_poly_T[df_poly_T$Name=="Variance ",], cex=2) +
  geom_dl(data=df_poly_T[df_poly_T$Name=="Variance ",],
          aes(label = Name, color = Name),
          method=list(dl.combine("first.points"), cex = 1, vjust=1.2)) +
  geom_line(data=df_poly_T[df_poly_T$Name==" 95% CR",], 
            aes(color=Name),
            linetype='dashed', 
            col = 'black', size=1)+
  geom_point(data=df_poly_T[df_poly_T$Name==" 95% CR",], cex=2)+
  geom_dl(data=df_poly_T[df_poly_T$Name==" 95% CR",],
          aes(label = Name), 
          method=list(dl.combine("last.points"), color="black", cex = 1)) +
  ggtitle("Keele & Titiunik (2015): Polynomial Turnout")
p_poly_T


df_poly_P <- data.frame("Scale"=c(rep(0.1,4), rep(0.2,4), rep(0.3,4), rep(0.4,4), rep(0.5,4)), "Name"=rep(c("Bias^2 ","Variance ","MSE "," 95% CR"),5), "Value"=rep(NA,20))

for (scale in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
  df_poly_P$Value[df_poly_P$Scale==scale] <- c(mean(poly_P$rf.bias[poly_P$scale==scale])^2, var(poly_P$rf[poly_P$scale==scale]), mean(poly_P$rf.bias[poly_P$scale==scale])^2+var(poly_P$rf[poly_P$scale==scale]), mean(poly_P$rf.cr[poly_P$scale==scale])*320-240)
}

p_poly_P <- ggplot(data=df_poly_P[df_poly_P$Name!=" 95% CR",], aes(x=Scale, y=Value, group=Name)) +
  geom_line(aes(color=Name), size=1) +
  geom_point(aes(color=Name), cex=2) +
  geom_dl(aes(label = Name, color=Name), method = list(dl.combine("first.points"), cex = 1)) +
  scale_y_continuous(name = "Squared bias, variance, MSE", sec.axis = sec_axis(~(.+240)/320, name="95% CR")) +
  scale_color_discrete(guide = "none") + 
  scale_x_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5), labels=c("Default","0.2","0.3","0.4","0.5"), limits = c(0.02,0.58)) +
  theme_bw(base_size = 15) +
  geom_line(data=df_poly_P[df_poly_P$Name==" 95% CR",], 
            aes(color=Name),
            linetype='dashed', 
            col = 'black', size=1)+
  geom_point(data=df_poly_P[df_poly_P$Name==" 95% CR",], cex=2)+
  geom_dl(data=df_poly_P[df_poly_P$Name==" 95% CR",],
          aes(label = Name), 
          method=list(dl.combine("last.points"), color="black", cex = 1)) +
  ggtitle("Keele & Titiunik (2015): Polynomial Price")
p_poly_P


df_poly_A <- data.frame("Scale"=c(rep(0.1,4), rep(0.2,4), rep(0.3,4), rep(0.4,4), rep(0.5,4)), "Name"=rep(c("Bias^2 ","Variance ","MSE "," 95% CR"),5), "Value"=rep(NA,20))

for (scale in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
  df_poly_A$Value[df_poly_A$Scale==scale] <- c(mean(poly_A$rf.bias[poly_A$scale==scale])^2, var(poly_A$rf[poly_A$scale==scale]), mean(poly_A$rf.bias[poly_A$scale==scale])^2+var(poly_A$rf[poly_A$scale==scale]), (mean(poly_A$rf.cr[poly_A$scale==scale])-0.916)/0.004)
}

p_poly_A <- ggplot(data=df_poly_A[df_poly_A$Name!=" 95% CR" & df_poly_A$Name!="Variance ",], aes(x=Scale, y=Value, group=Name)) +
  geom_line(aes(color=Name), size=1) +
  geom_point(aes(color=Name), cex=2) +
  geom_dl(aes(label = Name, color=Name), method = list(dl.combine("first.points"), cex = 1)) +
  scale_y_continuous(name = "Squared bias, variance, MSE", sec.axis = sec_axis(~.*0.004+0.916, name="95% CR")) +
  scale_color_discrete(guide = "none") + 
  scale_x_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5), labels=c("Default","0.2","0.3","0.4","0.5"), limits = c(0.02,0.58)) +
  theme_bw(base_size = 15) +
  geom_line(data=df_poly_A[df_poly_A$Name=="Variance ",],
            aes(color=Name), size=1) +
  geom_point(aes(color=Name), data=df_poly_A[df_poly_A$Name=="Variance ",], cex=2) +
  geom_dl(data=df_poly_A[df_poly_A$Name=="Variance ",],
          aes(label = Name, color = Name),
          method=list(dl.combine("first.points"), cex = 1, vjust=1.2)) +
  geom_line(data=df_poly_A[df_poly_A$Name==" 95% CR",], 
            aes(color=Name),
            linetype='dashed', 
            col = 'black', size=1)+
  geom_point(data=df_poly_A[df_poly_A$Name==" 95% CR",], cex=2)+
  geom_dl(data=df_poly_A[df_poly_A$Name==" 95% CR",],
          aes(label = Name), 
          method=list(dl.combine("last.points"), color="black", cex = 1)) +
  ggtitle("Keele & Titiunik (2015): Polynomial Age")
p_poly_A


### WGAN
df_wgan_lee <- data.frame("Scale"=c(rep(0.1,4), rep(0.2,4), rep(0.3,4), rep(0.4,4), rep(0.5,4)), "Name"=rep(c("Bias^2 ","Variance ","MSE "," 95% CR"),5), "Value"=rep(NA,20))

for (scale in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
  df_wgan_lee$Value[df_wgan_lee$Scale==scale] <- c(mean(wgan_lee$rf.bias[wgan_lee$scale==scale])^2, var(wgan_lee$rf[wgan_lee$scale==scale]), mean(wgan_lee$rf.bias[wgan_lee$scale==scale])^2+var(wgan_lee$rf[wgan_lee$scale==scale]), (mean(wgan_lee$rf.cr[wgan_lee$scale==scale])-0.85)/100)
}

p_wgan_lee <- ggplot(data=df_wgan_lee[df_wgan_lee$Name!=" 95% CR" & df_wgan_lee$Name!="Variance ",], aes(x=Scale, y=Value, group=Name)) +
  geom_line(aes(color=Name), size=1) +
  geom_point(aes(color=Name), cex=2) +
  geom_dl(aes(label = Name, color=Name), method = list(dl.combine("first.points"), cex = 1)) +
  scale_y_continuous(name = "Squared bias, variance, MSE", sec.axis = sec_axis(~.*100+0.85, name="95% CR")) +
  scale_color_discrete(guide = "none") + 
  scale_x_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5), labels=c("Default","0.2","0.3","0.4","0.5"), limits = c(0.02,0.58)) +
  theme_bw(base_size = 15) +
  geom_line(data=df_wgan_lee[df_wgan_lee$Name=="Variance ",],
            aes(color=Name), size=1) +
  geom_point(aes(color=Name), data=df_wgan_lee[df_wgan_lee$Name=="Variance ",], cex=2) +
  geom_dl(data=df_wgan_lee[df_wgan_lee$Name=="Variance ",],
          aes(label = Name, color = Name),
          method=list(dl.combine("first.points"), cex = 1, vjust=1.2)) +
  geom_line(data=df_wgan_lee[df_wgan_lee$Name==" 95% CR",], 
            aes(color=Name),
            linetype='dashed', 
            col = 'black', size=1)+
  geom_point(data=df_wgan_lee[df_wgan_lee$Name==" 95% CR",], cex=2)+
  geom_dl(data=df_wgan_lee[df_wgan_lee$Name==" 95% CR",],
          aes(label = Name), 
          method=list(dl.combine("last.points"), color="black", cex = 1)) +
  ggtitle("Lee (2008): WGAN")

p_wgan_lee


df_wgan_T <- data.frame("Scale"=c(rep(0.1,4), rep(0.2,4), rep(0.3,4), rep(0.4,4), rep(0.5,4)), "Name"=rep(c("Bias^2 ","Variance ","MSE "," 95% CR"),5), "Value"=rep(NA,20))

for (scale in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
  df_wgan_T$Value[df_wgan_T$Scale==scale] <- c(mean(wgan_T$rf.bias[wgan_T$scale==scale])^2, var(wgan_T$rf[wgan_T$scale==scale]), mean(wgan_T$rf.bias[wgan_T$scale==scale])^2+var(wgan_T$rf[wgan_T$scale==scale]), (mean(wgan_T$rf.cr[wgan_T$scale==scale])-0.82)/12)
}

p_wgan_T <- ggplot(data=df_wgan_T[df_wgan_T$Name!=" 95% CR" & df_wgan_T$Name!="Variance ",], aes(x=Scale, y=Value, group=Name)) +
  geom_line(aes(color=Name), size=1) +
  geom_point(aes(color=Name), cex=2) +
  geom_dl(aes(label = Name, color=Name), method = list(dl.combine("first.points"), cex = 1)) +
  scale_y_continuous(name = "Squared bias, variance, MSE", sec.axis = sec_axis(~.*12+0.82, name="95% CR")) +
  scale_color_discrete(guide = "none") + 
  scale_x_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5), labels=c("Default","0.2","0.3","0.4","0.5"), limits = c(0.02,0.58)) +
  theme_bw(base_size = 15) +
  geom_line(data=df_wgan_T[df_wgan_T$Name=="Variance ",],
            aes(color=Name), size=1) +
  geom_point(aes(color=Name), data=df_wgan_T[df_wgan_T$Name=="Variance ",], cex=2) +
  geom_dl(data=df_wgan_T[df_wgan_T$Name=="Variance ",],
          aes(label = Name, color = Name),
          method=list(dl.combine("first.points"), cex = 1)) +
  geom_line(data=df_wgan_T[df_wgan_T$Name==" 95% CR",], 
            aes(color=Name),
            linetype='dashed', 
            col = 'black', size=1)+
  geom_point(data=df_wgan_T[df_wgan_T$Name==" 95% CR",], cex=2)+
  geom_dl(data=df_wgan_T[df_wgan_T$Name==" 95% CR",],
          aes(label = Name), 
          method=list(dl.combine("last.points"), color="black", cex = 1)) +
  ggtitle("Keele & Titiunik (2015): WGAN Turnout")
p_wgan_T


df_wgan_P <- data.frame("Scale"=c(rep(0.1,4), rep(0.2,4), rep(0.3,4), rep(0.4,4), rep(0.5,4)), "Name"=rep(c("Bias^2 ","Variance ","MSE "," 95% CR"),5), "Value"=rep(NA,20))

for (scale in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
  df_wgan_P$Value[df_wgan_P$Scale==scale] <- c(mean(wgan_P$rf.bias[wgan_P$scale==scale])^2, var(wgan_P$rf[wgan_P$scale==scale]), mean(wgan_P$rf.bias[wgan_P$scale==scale])^2+var(wgan_P$rf[wgan_P$scale==scale]), mean(wgan_P$rf.cr[wgan_P$scale==scale])*320-250)
}

p_wgan_P <- ggplot(data=df_wgan_P[df_wgan_P$Name!=" 95% CR",], aes(x=Scale, y=Value, group=Name)) +
  geom_line(aes(color=Name), size=1) +
  geom_point(aes(color=Name), cex=2) +
  geom_dl(aes(label = Name, color=Name), method = list(dl.combine("first.points"), cex = 1)) +
  scale_y_continuous(name = "Squared bias, variance, MSE", sec.axis = sec_axis(~(.+250)/320, name="95% CR")) +
  scale_color_discrete(guide = "none") + 
  scale_x_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5), labels=c("Default","0.2","0.3","0.4","0.5"), limits = c(0.02,0.58)) +
  theme_bw(base_size = 15) +
  geom_line(data=df_wgan_P[df_wgan_P$Name==" 95% CR",], 
            aes(color=Name),
            linetype='dashed', 
            col = 'black', size=1)+
  geom_point(data=df_wgan_P[df_wgan_P$Name==" 95% CR",], cex=2)+
  geom_dl(data=df_wgan_P[df_wgan_P$Name==" 95% CR",],
          aes(label = Name), 
          method=list(dl.combine("last.points"), color="black", cex = 1)) +
  ggtitle("Keele & Titiunik (2015): WGAN Price")
p_wgan_P


df_wgan_A <- data.frame("Scale"=c(rep(0.1,4), rep(0.2,4), rep(0.3,4), rep(0.4,4), rep(0.5,4)), "Name"=rep(c("Bias^2 ","Variance ","MSE "," 95% CR"),5), "Value"=rep(NA,20))

for (scale in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
  df_wgan_A$Value[df_wgan_A$Scale==scale] <- c(mean(wgan_A$rf.bias[wgan_A$scale==scale])^2, var(wgan_A$rf[wgan_A$scale==scale]), mean(wgan_A$rf.bias[wgan_A$scale==scale])^2+var(wgan_A$rf[wgan_A$scale==scale]), (mean(wgan_A$rf.cr[wgan_A$scale==scale])-0.916)/0.004)
}

p_wgan_A <- ggplot(data=df_wgan_A[df_wgan_A$Name!=" 95% CR" & df_wgan_A$Name!="Variance ",], aes(x=Scale, y=Value, group=Name)) +
  geom_line(aes(color=Name), size=1) +
  geom_point(aes(color=Name), cex=2) +
  geom_dl(aes(label = Name, color=Name), method = list(dl.combine("first.points"), cex = 1)) +
  scale_y_continuous(name = "Squared bias, variance, MSE", sec.axis = sec_axis(~.*0.004+0.916, name="95% CR")) +
  scale_color_discrete(guide = "none") + 
  scale_x_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5), labels=c("Default","0.2","0.3","0.4","0.5"), limits = c(0.02,0.58)) +
  theme_bw(base_size = 15) +
  geom_line(data=df_wgan_A[df_wgan_A$Name=="Variance ",],
            aes(color=Name), size=1) +
  geom_point(aes(color=Name), data=df_wgan_A[df_wgan_A$Name=="Variance ",], cex=2) +
  geom_dl(data=df_wgan_A[df_wgan_A$Name=="Variance ",],
          aes(label = Name, color = Name),
          method=list(dl.combine("first.points"), cex = 1, vjust=1.2)) +
  geom_line(data=df_wgan_A[df_wgan_A$Name==" 95% CR",], 
            aes(color=Name),
            linetype='dashed', 
            col = 'black', size=1)+
  geom_point(data=df_wgan_A[df_wgan_A$Name==" 95% CR",], cex=2)+
  geom_dl(data=df_wgan_A[df_wgan_A$Name==" 95% CR",],
          aes(label = Name), 
          method=list(dl.combine("last.points"), color="black", cex = 1)) +
  ggtitle("Keele & Titiunik (2015): WGAN Age")
p_wgan_A

fig <- arrangeGrob(p_poly_lee, p_poly_T, p_poly_P, p_poly_A, ncol = 2)
ggsave("scaling_poly.png", fig, width = 12, height = 8)
fig <- arrangeGrob(p_wgan_lee, p_wgan_T, p_wgan_P, p_wgan_A, ncol = 2)
ggsave("scaling_wgan.png", fig, width = 12, height = 8)
