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

r1 = c(0.18,-0.06,0.32,-0.18,-0.24)
r2 = c(-0.06,0.14,0.25,0.15,-0.18)
r3 = c(0.32,0.25,0.2,0.17,-0.22)
r4 = c(-0.18,0.15,0.17,0.25,0.11)
r5 = c(-0.24,-0.18,-0.22,0.11,0.27)

Dat.XW = c(r1,r2,r3,r4,r5)
cov.matrix.XW = matrix(Dat.XW, nrow = 5, byrow = T)
determinant(cov.matrix.XW, logarithm = F)$modulus[1]

X.vector = rmvnorm(N, mean = rep(0,5), sigma = cov.matrix.XW)

x1 = X.vector[,1]
x2 = X.vector[,2]
w1 = X.vector[,3]
w2 = X.vector[,4]
w3 = X.vector[,5]

eps = rnorm(N, 0, sqrt(0.015))

m = function(t,u,v) 0.45*t*u-0.25*t^2*v+v^3
m.w1.w2.w3 = m(w1,w2,w3)

Y = rt(N,2,0)

h.x1 = 0.9*min(sd(x1),(IQR(x1)/1.34))*N^(-1/5)
h.x2 = 0.9*min(sd(x2),(IQR(x2)/1.34))*N^(-1/5)
h.w1 = 0.9*min(sd(w1),(IQR(w1)/1.34))*N^(-1/5)
h.w2 = 0.9*min(sd(w2),(IQR(w2)/1.34))*N^(-1/5)
h.w3 = 0.9*min(sd(w3),(IQR(w3)/1.34))*N^(-1/5)
h.Y = 0.9*min(sd(Y),(IQR(Y)/1.34))*N^(-1/5)

Kern = function(f) dnorm(f)

w1.star <- c()
x1.star <- c()
x2.star<- c() 
w2.star<- c() 
w3.star <- c()
Y.star<- c() 
for(i in 1:N)
{
  w1.star[i] = w1[rank(w1)==i]
  x1.star[i]<- x1[which(w1==w1.star[i])]
  x2.star[i]<- x2[which(w1==w1.star[i])]
  w2.star[i]<- w2[which(w1==w1.star[i])]
  w3.star[i]<- w3[which(w1==w1.star[i])]
  Y.star[i]<- Y[which(w1==w1.star[i])]
}

####################################################################################

weight.w1.w2.w3 = function(s1,s2,s3)
{
  wt = c()
  for(i in 1:N)
  {
    wt[i] = (1/h.w1)*Kern((s1-w1[i])/h.w1)*(1/h.w2)*Kern((s2-w2[i])/h.w2)*(1/h.w3)*Kern((s3-w3[i])/h.w3)
  }
  return(wt)
}

density.w1w2w3 = c()
for(j in 1:N)
{
  density.w1w2w3[j] = mean(weight.w1.w2.w3(w1.star[j],w2.star[j],w3.star[j]))
}

num.g_Y.W1w2w3 = c()
num.g_x1.W1w2w3 = c()
num.g_x2.W1w2w3 = c()
for(j in 1:N)
{
  num.g_Y.W1w2w3[j] = mean(weight.w1.w2.w3(w1.star[j],w2.star[j],w3.star[j])*Y.star)
  num.g_x1.W1w2w3[j] = mean(weight.w1.w2.w3(w1.star[j],w2.star[j],w3.star[j])*x1.star)
  num.g_x2.W1w2w3[j] = mean(weight.w1.w2.w3(w1.star[j],w2.star[j],w3.star[j])*x2.star)
}
g.hat_Y.W1w2w3 = num.g_Y.W1w2w3/density.w1w2w3
e.hat_Y.W1w2w3 = Y.star-g.hat_Y.W1w2w3

g.hat_x1.W1w2w3 = num.g_x1.W1w2w3/density.w1w2w3
g.hat_x2.W1w2w3 = num.g_x2.W1w2w3/density.w1w2w3

e.hat_x1.W1w2w3 = x1 - g.hat_x1.W1w2w3
e.hat_x2.W1w2w3 = x2 - g.hat_x2.W1w2w3

e.hat_x.matrix <- matrix(nrow = N, ncol = 2)
for(column in 1:2){
  e.hat_x.matrix[, 1] <- e.hat_x1.W1w2w3
  e.hat_x.matrix[, 2] <- e.hat_x2.W1w2w3
}

beta.hat = solve(t(e.hat_x.matrix) %*% e.hat_x.matrix) %*% t(e.hat_x.matrix) %*% e.hat_Y.W1w2w3
beta1.hat = beta.hat[1,]
beta2.hat = beta.hat[2,]
cbind(beta1.hat,beta2.hat)

Y.star.dash = Y.star-beta1.hat*x1.star-beta2.hat*x2.star

num.g_Y.dash.W1w2w3 = c()
for(j in 1:N)
{
  num.g_Y.dash.W1w2w3[j] = mean(weight.w1.w2.w3(w1.star[j],w2.star[j],w3.star[j])*Y.star.dash)
}

m.hat.w1.w2.w3 = num.g_Y.dash.W1w2w3/density.w1w2w3
Y.star.hat = beta1.hat*x1.star+beta2.hat*x2.star+m.hat.w1.w2.w3

df_Y.star.Y.star.hat = cbind.data.frame(Y.star,Y.star.hat)
df_Y.star.Y.star.hat

## r = 2, m = 2 ##
r = 4
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

c1 = -0.77
c2 = 0.34
c3 = 0.72
c4 = -0.56
c5 = 0.65
c.vec = c(c1,c2,c3,c4,c5)
S = c()
df = c()
for(i in 1:N)
{
  df[i] = ceiling(1/(abs(1+t(c.vec)%*%X.vector[i,])))
}
eps1 = c()
for(i in 1:N)
{
  eps1[i] = (rchisq(1,df[i],0)-df[i])/sqrt(2*df[i])
}

Y1 = beta1.hat*x1+beta2.hat*x2+m.w1.w2.w3+eps1
h.Y1 = 0.9*min(sd(Y1),(IQR(Y1)/1.34))*N^(-1/5)
Y1.star<- c() 
for(i in 1:N)
{
  Y1.star[i]<- Y1[which(w1==w1.star[i])]
}

Y1.star.dash = Y1.star-beta1.hat*x1.star-beta2.hat*x2.star

num.g_Y1.dash.w1w2w3 = c()
for(j in 1:N)
{
  num.g_Y1.dash.w1w2w3[j] = mean(weight.w1.w2.w3(w1.star[j],w2.star[j],w3.star[j])*Y1.star.dash)
}

m.hat1.w1.w2.w3 = num.g_Y.dash.W1w2w3/density.w1w2w3
Y1.star.hat = beta1.hat*x1.star+beta2.hat*x2.star+m.hat1.w1.w2.w3

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
Power.Tn = c()
mu = 0:round(sqrt(n),0)
for(i in 1:length(mu)) {Power.Tn[i] = Pow(mu[i])}
Power.Tn = round(Power.Tn,4)
Power.Tn.table = cbind(mu, Power.Tn)
Power.Tn.table

setwd('C:/Users/hp/Desktop/R-programs-for-gen-semipar/Paper_for_ADAMAS_CONF')
write.csv(Power.Tn.table, file = 'Simulated_powers_Tn_r=4-Ex6_Robinson.csv')
plot(mu, Power.Tn, xlim=c(mu[1],mu[length(mu)]), ylim=c(0,1), type="l", pch=10, col="black", xlab=expression(mu), ylab=expression("Power (Vn)"))
