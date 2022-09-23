library(haven)
d <- read_dta("Desktop/Grad/Papers/Liu/Lee/Lee2008/group_final.dta")

d$diff <- d$difdemshare
d$diff_2 <- (d$diff)^2
d$diff_3 <- (d$diff)^3
d$diff_4 <- (d$diff)^4
d$diff_5 <- (d$diff)^5

m1_left <- lm(mdemsharenext ~ diff + diff_2 + diff_3 
        + diff_4 + diff_5, data = d[d$diff < 0 & d$diff <= 0.99 & d$diff >= -0.99,])
summary(m1_left)

m1_right <- lm(mdemsharenext ~ diff + diff_2 + diff_3 
             + diff_4 + diff_5, data = d[d$diff >= 0 & d$diff <= 0.99 & d$diff >= -0.99,])
summary(m1_right)


m2_left <- lm(mdemsharenext ~ diff + diff_2 + diff_3 
             + diff_4 + diff_5 + mdemshareprev, data = d[d$diff < 0 & d$diff <= 0.99 & d$diff >= -0.99,])
summary(m2_left)

m2_right <- lm(mdemsharenext ~ diff + diff_2 + diff_3 
              + diff_4 + diff_5 + mdemshareprev, data = d[d$diff >= 0 & d$diff <= 0.99 & d$diff >= -0.99,])
summary(m2_right)

m2z_left <- lm(mdemshareprev ~ diff + diff_2 + diff_3 
              + diff_4 + diff_5, data = d[d$diff < 0 & d$diff <= 0.99 & d$diff >= -0.99,])
summary(m2z_left)

m2z_right <- lm(mdemshareprev ~ diff + diff_2 + diff_3 
               + diff_4 + diff_5, data = d[d$diff >= 0 & d$diff <= 0.99 & d$diff >= -0.99,])
summary(m2z_right)
