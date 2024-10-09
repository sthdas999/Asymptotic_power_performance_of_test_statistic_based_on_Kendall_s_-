library(MASS)
library(kedd)
library(mvtnorm)
library(Matrix)
library(truncnorm)
library(npmlda)
library(stats)
library(EnvStats)
library(lattice)
library(pracma)


## Ex - 7 ##

setwd('C:/Users/hp/Desktop/R-programs-for-gen-semipar/Paper_for_ADAMAS_CONF/Ex - 7,8')
G = read.csv(file = 'Simulated_powers_Tn_r=2,3,4,5,10-Ex7_Robinson-ADAMAS_c=(5.92, -3.78, -10.66, 8.89, -5.45, 9.65, 8.35, -7.89).conf.csv')
mu = G[,1]
Power.Vn.order.2 = G[,2]
Power.Vn.order.3 = G[,3]
Power.Vn.order.4 = G[,4]
Power.Vn.order.5 = G[,5]
Power.Vn.order.10 = G[,6]


par(mar = c(4.5, 4.5, 3, 0.21), xpd = TRUE)
plot(mu, Power.Vn.order.2, xlim=c(mu[1],mu[length(mu)]), ylim=c(0,1), type="l", pch=3, col="red", xlab=expression(mu), ylab=expression("Power"))
# Add a line
lines(mu, Power.Vn.order.3, xlim=c(mu[1],mu[length(mu)]), ylim=c(0,1), pch=3, col="black", type="l")
# Add a line
lines(mu, Power.Vn.order.4, xlim=c(mu[1],mu[length(mu)]), ylim=c(0,1), pch=3, col="blue", type="l")
# Add a line
lines(mu, Power.Vn.order.5, xlim=c(mu[1],mu[length(mu)]), ylim=c(0,1), pch=3, col="yellow", type="l")
# Add a line
lines(mu, Power.Vn.order.10, xlim=c(mu[1],mu[length(mu)]), ylim=c(0,1), pch=3, col="green", type="l")
# Add a legend
legend("bottomright",
       c(expression(paste('r=2'),paste('r=3'),paste('r=4'),paste('r=5'),paste('r=10'))),
       fill=c("red","black","blue","yellow","green"))

