######################################
# Figure A.4: Parameter tuning for 
# forests with varying scaling 
# constants for sample fraction
######################################

library(ggplot2)

# "FigureA.4.csv" collects results from running grf with different scaling constants for sample fraction using the tuning procedure described in Appendix A.4
results <- read.csv("FigureA.4.csv")
results$scale <- ifelse(results$scale==0, 0.1, results$scale) # note scale=0.1 means no tuning, i.e., all default parameters are used
cr <- results[results$type=="cr",]
other <- results[results$type!="cr",]

## Panel A: Lee (2008)
lee <- ggplot(other[other$model=="lee",], aes(x = scale, y = result)) +
  geom_line(aes(colour = type)) + 
  geom_point(aes(colour = type))+
  # geom_text(aes(label=result))+
  scale_y_continuous(
    name="Squared bias, variance, MSE",
    sec.axis = sec_axis(trans = ~.*100+0.73,name="95% CR"))+
  geom_line(aes(y = (result-0.73)/100),
              data = cr[cr$model=="lee",], linetype = 2)+
  geom_point(aes(y = (result-0.73)/100),
             data = cr[cr$model=="lee",])+
  # geom_text(aes(label=result,y = (result-0.73)/100),data = cr[cr$model=="lee",])+
  theme_bw()+theme(legend.position="none")+
  scale_x_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5), labels=c("Default","0.2","0.3","0.4","0.5"))+xlab("Scale")+ggtitle("Panel A: Lee (2008)")


## Panel B: Keele & Titiunik (2015): Turnout
turnout <- ggplot(other[other$model=="turnout",], aes(x = scale, y = result)) +
  geom_line(aes(colour = type)) + 
  geom_point(aes(colour = type))+
  scale_y_continuous(
    name="Squared bias, variance, MSE",
    sec.axis = sec_axis(trans = ~.*8+0.85,name="95% CR"))+
  geom_line(aes(y = (result-0.85)/8),
            data = cr[cr$model=="turnout",], linetype = 2)+
  geom_point(aes(y = (result-0.85)/8),
             data = cr[cr$model=="turnout",])+
  theme_bw()+theme(legend.position="none")+
  scale_x_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5), labels=c("Default","0.2","0.3","0.4","0.5"))+xlab("Scale")+ggtitle("Panel B: Keele & Titiunik (2015): Turnout")

## Panel C: Keele & Titiunik (2015): Price
price <- ggplot(other[other$model=="price",], aes(x = scale, y = result)) +
  geom_line(aes(colour = type)) + 
  geom_point(aes(colour = type))+
  scale_y_continuous(
    name="Squared bias, variance, MSE",
    sec.axis = sec_axis(trans = ~.*0.002+0.84,name="95% CR"))+
  geom_line(aes(y = (result-0.84)/0.002),
            data = cr[cr$model=="price",], linetype = 2)+
  geom_point(aes(y = (result-0.84)/0.002),
             data = cr[cr$model=="price",])+
  theme_bw()+theme(legend.position="none")+
  scale_x_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5), labels=c("Default","0.2","0.3","0.4","0.5"))+xlab("Scale")+ggtitle("Panel C: Keele & Titiunik (2015): Price")

## Panel D: Keele & Titiunik (2015): Age
age <- ggplot(other[other$model=="age",], aes(x = scale, y = result)) +
  geom_line(aes(colour = type)) + 
  geom_point(aes(colour = type))+
  scale_y_continuous(
    name="Squared bias, variance, MSE",
    sec.axis = sec_axis(trans = ~.*0.004+0.916,name="95% CR"))+
  geom_line(aes(y = (result-0.916)/0.004),
            data = cr[cr$model=="age",], linetype = 2)+
  geom_point(aes(y = (result-0.916)/0.004),
             data = cr[cr$model=="age",])+
  theme_bw()+theme(legend.position="none")+
  scale_x_continuous(breaks=c(0.1,0.2,0.3,0.4,0.5), labels=c("Default","0.2","0.3","0.4","0.5"))+xlab("Scale")+ggtitle("Panel D: Keele & Titiunik (2015): Age")

ggarrange(lee, turnout, price, age, ncol = 2, nrow = 2)
