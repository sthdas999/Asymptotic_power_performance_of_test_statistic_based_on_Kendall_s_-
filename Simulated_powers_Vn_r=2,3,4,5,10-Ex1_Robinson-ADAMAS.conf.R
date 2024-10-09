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

################## H0 #####################

N = 1200

Dat.XW = c(0.18,-0.06,0.22,-0.13,-0.06,0.14,-0.28,0.19,0.22,-0.28,0.2,0.17,-0.13,0.19,0.17,0.25)
cov.matrix.XW = matrix(Dat.XW, nrow = 4, byrow = T)
determinant(cov.matrix.XW, logarithm = F)$modulus[1]

X.vector = rmvnorm(N, mean = rep(0,4), sigma = cov.matrix.XW)

x1 = X.vector[,1]
x2 = X.vector[,2]
w1 = X.vector[,3]
w2 = X.vector[,4]

eps = rnorm(N, 0, sqrt(0.015))

m = function(t,u) 0.45*t*u-0.25*t^2*u+u^3
m.w1.w2 = m(w1,w2)

Y = rt(N,2,0)

h.x1 = 0.9*min(sd(x1),(IQR(x1)/1.34))*N^(-1/5)
h.x2 = 0.9*min(sd(x2),(IQR(x2)/1.34))*N^(-1/5)
h.w1 = 0.9*min(sd(w1),(IQR(w1)/1.34))*N^(-1/5)
h.w2 = 0.9*min(sd(w2),(IQR(w2)/1.34))*N^(-1/5)
h.Y = 0.9*min(sd(Y),(IQR(Y)/1.34))*N^(-1/5)

Kern = function(f) dnorm(f)

w1.star <- c()
x1.star <- c()
x2.star<- c() 
w2.star <- c()
Y.star<- c() 
for(i in 1:N)
{
  w1.star[i] = w1[rank(w1)==i]
  x1.star[i]<- x1[which(w1==w1.star[i])]
  x2.star[i]<- x2[which(w1==w1.star[i])]
  w2.star[i]<- w2[which(w1==w1.star[i])]
  Y.star[i]<- Y[which(w1==w1.star[i])]
}

weight.w1.w2 = function(s1,s2)
{
  wt = c()
  for(i in 1:N)
  {
    wt[i] = (1/h.w1)*Kern((s1-w1[i])/h.w1)*(1/h.w2)*Kern((s2-w2[i])/h.w2)
  }
  return(wt)
}

density.w1w2 = c()
for(j in 1:N)
{
  density.w1w2[j] = mean(weight.w1.w2(w1.star[j],w2.star[j]))
}

num.g_Y.W1w2 = c()
num.g_x1.W1w2 = c()
num.g_x2.W1w2 = c()
for(j in 1:N)
{
  num.g_Y.W1w2[j] = mean(weight.w1.w2(w1.star[j],w2.star[j])*Y.star)
  num.g_x1.W1w2[j] = mean(weight.w1.w2(w1.star[j],w2.star[j])*x1.star)
  num.g_x2.W1w2[j] = mean(weight.w1.w2(w1.star[j],w2.star[j])*x2.star)
}
g.hat_Y.W1w2 = num.g_Y.W1w2/density.w1w2
e.hat_Y.W1w2 = Y.star-g.hat_Y.W1w2

g.hat_x1.W1w2 = num.g_x1.W1w2/density.w1w2
g.hat_x2.W1w2 = num.g_x2.W1w2/density.w1w2

e.hat_x1.W1w2 = x1 - g.hat_x1.W1w2
e.hat_x2.W1w2 = x2 - g.hat_x2.W1w2

e.hat_x.matrix <- matrix(nrow = N, ncol = 2)
for(column in 1:2){
  e.hat_x.matrix[, 1] <- e.hat_x1.W1w2
  e.hat_x.matrix[, 2] <- e.hat_x2.W1w2
}

beta.hat = solve(t(e.hat_x.matrix) %*% e.hat_x.matrix) %*% t(e.hat_x.matrix) %*% e.hat_Y.W1w2
beta1.hat = beta.hat[1,]
beta2.hat = beta.hat[2,]
cbind(beta1.hat,beta2.hat)

Y.star.dash = Y.star-beta1.hat*x1.star-beta2.hat*x2.star

num.g_Y.dash.W1w2 = c()
for(j in 1:N)
{
  num.g_Y.dash.W1w2[j] = mean(weight.w1.w2(w1.star[j],w2.star[j])*Y.star.dash)
}

m.hat.w1.w2 = num.g_Y.dash.W1w2/density.w1w2
Y.star.hat = beta1.hat*x1.star+beta2.hat*x2.star+m.hat.w1.w2

df_Y.star.Y.star.hat = cbind.data.frame(Y.star,Y.star.hat)
df_Y.star.Y.star.hat

## r = 2, m = 2 ##
r = 3
m = 2
n = 1000
B = 400

Y.star.hat.0 = diff(Y.star.hat,r)[1:n]
Y.star.0 = diff(Y.star,r)[1:n]

Y.star.hat.0.boot1 = remp(B,Y.star.hat.0)
Y.star.hat.0.boot2 = remp(B,Y.star.hat.0)

Y.star.0.boot1 = remp(B,Y.star.0)
Y.star.0.boot2 = remp(B,Y.star.0)

Func = function(p,q,r) {sign(p-q)*sign(p-r)}
Y.star.hat.0.val = function(a)
{
  V = c()
  for(j in 1:B)
  {
    V[j] = Func(a,Y.star.hat.0.boot1[j],Y.star.hat.0.boot2[j])
  }
  return(V)
}
Y.star.hat.0.simulated = c()
for(k in 1:n)
{
  Y.star.hat.0.simulated[k] = mean(Y.star.hat.0.val(Y.star.hat.0[k]))
}

################################################################################

h.Y.star.hat.r = 0.9*min(sd(Y.star.hat.0),(IQR(Y.star.hat.0)/1.34))*n^(-1/5)
h.Y.star.r = 0.9*min(sd(Y.star.0),(IQR(Y.star.0)/1.34))*n^(-1/5)

f0.hat = function(u,v)
{
  vec =  (1/h.Y.star.hat.r)*(1/h.Y.star.r)*mean(Kern((u-Y.star.hat.0)/h.Y.star.hat.r))*mean(Kern((v-Y.star.0)/h.Y.star.r))
  return(vec)
}

################## H1 #####################

c1 = -0.13
c2 = -2.1
c3 = 5.6
c4 = -0.95
c.vec = c(c1,c2,c3,c4)
cond.var.eps = c()
for(i in 1:N)
{
  cond.var.eps[i] = 0.015*abs(1+t(c.vec)%*%X.vector[i,])
}

eps1 = c()
for(i in 1:N)
{
  eps1[i] = rnorm(1,0,sqrt(cond.var.eps[i]))
}

Y1 = beta1.hat*x1+beta2.hat*x2+m.w1.w2+eps1
h.Y1 = 0.9*min(sd(Y1),(IQR(Y1)/1.34))*N^(-1/5)
Y1.star<- c() 
for(i in 1:N)
{
  Y1.star[i]<- Y1[which(w1==w1.star[i])]
}

Y1.star.dash = Y1.star-beta1.hat*x1.star-beta2.hat*x2.star

num.g_Y1.dash.W1w2 = c()
for(j in 1:N)
{
  num.g_Y1.dash.W1w2[j] = mean(weight.w1.w2(w1.star[j],w2.star[j])*Y1.star.dash)
}

m.hat1.w1.w2 = num.g_Y.dash.W1w2/density.w1w2
Y1.star.hat = beta1.hat*x1.star+beta2.hat*x2.star+m.hat1.w1.w2

df_Y1.star.Y1.star.hat = cbind.data.frame(Y1.star,Y1.star.hat)

Y.star.hat.1 = diff(Y1.star.hat,r)[1:n]
Y.star.1 = diff(Y1.star,r)[1:n]

Y.star.hat.1.boot = remp(B,Y.star.hat.1)
Y.star.1.boot = remp(B,Y.star.1)

h.Y1.star.hat.r = 0.9*min(sd(Y.star.hat.1),(IQR(Y.star.hat.1)/1.34))*n^(-1/5)
h.Y1.star.r = 0.9*min(sd(Y.star.1),(IQR(Y.star.1)/1.34))*n^(-1/5)

f.hat = function(u,v)
{
  vec1 = (1/h.Y1.star.hat.r)*(1/h.Y1.star.r)*mean(Kern((u-Y.star.hat.1)/h.Y1.star.hat.r))*mean(Kern((v-Y.star.1)/h.Y1.star.r))
  return(vec1)
}

h_function = function(a,b,c,d) {sign((a-b)*(c-d))}
h.kernel = c()
ratio = c()
for(j in 1:B)
{
  h.kernel[j] = h_function(Y.star.hat.0.boot1[[j]],Y.star.hat.0.boot2[[j]],Y.star.0.boot1[[j]],Y.star.0.boot2[[j]])
  ratio[j] = f.hat(Y.star.hat.0.boot1[[j]],Y.star.0.boot1[[j]])/f0.hat(Y.star.hat.0.boot1[[j]],Y.star.0.boot1[[j]])
}

Lambda.val = h.kernel*ratio
Lambda.sim = mean(Lambda.val)
Lambda.sim
variance.Vn = m^2*(1/9)
alpha = 0.05
crit.val = qnorm(1-alpha,mean = 0, sd = sqrt(variance.Vn))

Pow = function(Y) { 1-pnorm((crit.val-Y*Lambda.sim)/sqrt(variance.Vn)) }
power.Vn = c()
mu = 0:round(sqrt(n),0)
for(i in 1:length(mu)) {power.Vn[i] = Pow(mu[i])}
Power.Vn = round(power.Vn,4)
Power.Vn.table = cbind(mu, Power.Vn)
Power.Vn.table

setwd('C:/Users/hp/Desktop/R-programs-for-gen-semipar/Paper_for_ADAMAS_CONF')
write.csv(Power.Vn.table, file = 'Simulated_powers_Vn_r=3-Ex1_Robinson.csv')
plot(mu, Power.Vn, xlim=c(mu[1],mu[length(mu)]), ylim=c(0,1), type="l", pch=10, col="black", xlab=expression(mu), ylab=expression("Power (Vn)"))
