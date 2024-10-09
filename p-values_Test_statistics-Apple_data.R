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
library(boot)

N = 160
setwd('C:/Users/hp/Downloads/Apple_quality data')
apple.dat = read.csv('apple_quality.csv')
apple.dat = head(apple.dat, N)

size = apple.dat$Size
weight = apple.dat$Weight
sweetness = apple.dat$Sweetness
crunchiness = apple.dat$Crunchiness
juiciness = apple.dat$Juiciness
ripeness = apple.dat$Ripeness
quality = apple.dat$Quality  ## qualitative ##

acidity = apple.dat$Acidity  ## response ##

plot(size,acidity,type = 'p')  ## nonparametric relationship ##
plot(crunchiness,acidity,type = 'p')  ## nonparametric relationship ##
plot(weight,acidity,type = 'p')  ## nonparametric relationship ##
plot(sweetness,acidity,type = 'p')  ## nonparametric relationship ##

plot(ripeness,acidity,type = 'p')  ## approximately decreasing relationship ##
plot(juiciness,acidity,type = 'p')  ## approximately increasing relationship ##


######################################################################################

x1 = juiciness
x2 = ripeness

w1 = size
w2 = crunchiness
w3 = weight
w4 = sweetness

Y = acidity

N = length(acidity)

h.x1 = 0.9*min(sd(x1),(IQR(x1)/1.34))*N^(-1/5)
h.x2 = 0.9*min(sd(x2),(IQR(x2)/1.34))*N^(-1/5)

h.w1 = 0.9*min(sd(w1),(IQR(w1)/1.34))*N^(-1/5)
h.w2 = 0.9*min(sd(w2),(IQR(w2)/1.34))*N^(-1/5)
h.w3 = 0.9*min(sd(w3),(IQR(w3)/1.34))*N^(-1/5)
h.w4 = 0.9*min(sd(w4),(IQR(w4)/1.34))*N^(-1/5)

h.Y = 0.9*min(sd(Y),(IQR(Y)/1.34))*N^(-1/5)

#####################################################################################

Kern = function(f) dnorm(f)

w1.star <- c()
x1.star <- c()
x2.star<- c()
w2.star <- c()
w3.star<- c()
w4.star<- c()
Y.star<- c() 
for(i in 1:N)
{
  w1.star[i] = w1[rank(w1)==i]
  x1.star[i]<- x1[which(w1==w1.star[i])]
  x2.star[i]<- x2[which(w1==w1.star[i])]
  w2.star[i]<- w2[which(w1==w1.star[i])]
  w3.star[i]<- w3[which(w1==w1.star[i])]
  w4.star[i]<- w4[which(w1==w1.star[i])]
  Y.star[i]<- Y[which(w1==w1.star[i])]
}

weight.w1.w2.w3.w4 = function(s1,s2,s3,s4)
{
  wt = c()
  for(i in 1:N)
  {
    wt[i] = (1/h.w1)*Kern((s1-w1[i])/h.w1)*(1/h.w2)*Kern((s2-w2[i])/h.w2)*(1/h.w3)*Kern((s3-w3[i])/h.w3)*Kern((s4-w4[i])/h.w4)
  }
  return(wt)
}

density.w1w2w3w4 = c()
for(j in 1:N)
{
  density.w1w2w3w4[j] = mean(weight.w1.w2.w3.w4(w1.star[j],w2.star[j],w3.star[j],w4.star[j]))
}

num.g_Y.W1w2w3w4 = c()
num.g_x1.W1w2w3w4 = c()
num.g_x2.W1w2w3w4 = c()
for(j in 1:N)
{
  num.g_Y.W1w2w3w4[j] = mean(weight.w1.w2.w3.w4(w1.star[j],w2.star[j],w3.star[j],w4.star[j])*Y.star)
  num.g_x1.W1w2w3w4[j] = mean(weight.w1.w2.w3.w4(w1.star[j],w2.star[j],w3.star[j],w4.star[j])*x1.star)
  num.g_x2.W1w2w3w4[j] = mean(weight.w1.w2.w3.w4(w1.star[j],w2.star[j],w3.star[j],w4.star[j])*x2.star)
}
g.hat_Y.W1w2w3w4 = num.g_Y.W1w2w3w4/density.w1w2w3w4
e.hat_Y.W1w2w3w4 = Y.star-g.hat_Y.W1w2w3w4

g.hat_x1.W1w2w3w4 = num.g_x1.W1w2w3w4/density.w1w2w3w4
g.hat_x2.W1w2w3w4 = num.g_x2.W1w2w3w4/density.w1w2w3w4

e.hat_x1.W1w2w3w4 = x1 - g.hat_x1.W1w2w3w4
e.hat_x2.W1w2w3w4 = x2 - g.hat_x2.W1w2w3w4

e.hat_x.matrix <- matrix(nrow = N, ncol = 2)
for(column in 1:2){
  e.hat_x.matrix[, 1] <- e.hat_x1.W1w2w3w4
  e.hat_x.matrix[, 2] <- e.hat_x2.W1w2w3w4
}

beta.hat = solve(t(e.hat_x.matrix) %*% e.hat_x.matrix) %*% t(e.hat_x.matrix) %*% e.hat_Y.W1w2w3w4
beta1.hat = beta.hat[1,]
beta1.hat
beta2.hat = beta.hat[2,]
beta2.hat

Y.star.dash = Y.star-beta1.hat*x1.star-beta2.hat*x2.star
num.g_Y.dash.W1w2w3w4 = c()
for(j in 1:N)
{
  num.g_Y.dash.W1w2w3w4[j] = mean(weight.w1.w2.w3.w4(w1.star[j],w2.star[j],w3.star[j],w4.star[j])*Y.star.dash)
}

m.hat.w1.w2.w3.w4 = num.g_Y.dash.W1w2w3w4/density.w1w2w3w4

#####################################################################################

Y.star.hat = beta1.hat*x1.star+beta2.hat*x2.star+m.hat.w1.w2.w3.w4
df_Y.star.Y.star.hat = cbind.data.frame(Y.star,Y.star.hat)
df_Y.star.Y.star.hat
Tn = function(u, v)
{
  temp = 0
  for(i in 1:(n-1))
  {
    for(j in (i+1):n)
    {
      temp = temp + sign((u[i]-u[j])*(v[i]-v[j]))
    }
  }
  return(temp/choose(n, 2))
}

n = N-10
B = 1000
p_value.Tn = function(r)
{
  Y.star.hat.0 = diff(Y.star.hat,r)[1:n]
  Y.star.0 = diff(Y.star,r)[1:n]
  
  Tn.boot = boot(data = c(df_Y.star.Y.star.hat[,1],df_Y.star.Y.star.hat[,2]),Tn,R = 300)
  Tn.B = Tn.boot$t
  
  Y.star.boot = vector("list",B)
  Y.star.hat.boot = vector("list",B)
  for(j in 1:B)
  {
    Y.star.boot[[j]] = sample(Y.star.0,n)
    Y.star.hat.boot[[j]] = sample(Y.star.hat.0,n)   
  }
  
  Tn.0 = c()
  for(j in 1:B) {Tn.0[j] = Tn(Y.star.boot[[j]],Y.star.hat.boot[[j]])}
  
  Tn.critical = quantile(Tn.B, probs = 0.95)
  pvl = mean(Tn.0>Tn.critical)
  return(pvl)
}

p_value.Tn(2)
p_value.Tn(3)
p_value.Tn(4)
p_value.Tn(5)
p_value.Tn(10)

