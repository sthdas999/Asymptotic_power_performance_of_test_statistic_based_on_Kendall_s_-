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

####  Tau ####

## r = 2, m = 2 ##
r = 2
m = 2
n = 100
B = 500

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

c1 = -1.5
c2 = -1.7
c3 = 1.2
c4 = 1.3
c.vec = c(c1,c2,c3,c4)
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
variance.tau = m^2*(1/9)
alpha = 0.05
crit.val = qnorm(1-alpha,mean = 0, sd = sqrt(variance.tau))

Pow = function(Y) { 1-pnorm((crit.val-Y*Lambda.sim)/sqrt(variance.tau)) }
power.tau = c()
mu = c(0:30)
for(i in 1:length(mu)) {power.tau[i] = Pow(mu[i])}
power.tau
###################################################################################################

####  Tau* ####

n1 = 5
dt.x<- Y.star.hat.0
dt.y<- Y.star.0
data.pts.x1<- unname(quantile(dt.x, probs = ((1:n1)/n1)))
data.pts.y1<- unname(quantile(dt.y, probs = ((1:n1)/n1)))
points.all<- cbind(data.pts.x1, data.pts.y1)
kernel.T2<-function(a,b)  ## kernel function of Tn2 ##
{
  x<-0
  for(i in 1:(n-3))
  {
    for(j in (i+1):(n-2))
    {
      for(k in (j+1):(n-1))
      {
        for(l in (k+1):n)
        {
          x<-x+(sign(abs(a-dt.x[j])+abs(dt.x[k]-dt.x[l])-abs(a-dt.x[k])-abs(dt.x[j]-dt.x[l]))*sign(abs(b-dt.y[j])+abs(dt.y[k]-dt.y[l])-abs(b-dt.y[k])-abs(dt.y[j]-dt.y[l])))
        }
      }
    }
  }
  return(x/choose(n,4))
}
L<- matrix(0,n1,n1)
{
  for(p in 1:n1)
  {
    for(q in 1:n1)
    {
      L[p,q]<- kernel.T2(data.pts.x1[p],data.pts.y1[q])
    }
  }
}
lmat2<- (L+t(L))/2  ## L and lmat2 have same eigenvalues and eigenvectors. ##
E<- eigen(lmat2) ## eigenvalue decomposition of a real symmetric matrix A: A=VDV^T ##
lambdas2<- E$values ## real eigenvalues ##
V2<- E$vector ## real eigenvectors ##
V2.t<- t(V2)
h1<- ratio-1
b2.vals<- vector("list", n1)
for(k in 1:n1)
{
  for(i in 1:n1)
  {
    b2.vals[[k]][i]<- h1[1:n1][i]*V2[i,k]*V2.t[k,i]
  }
}
b2.vals1<- c()  ## computation of a_k ##
for(k in 1:n1)
{
  b2.vals1[k]<- mean(b2.vals[[k]])
}
Z.val<- vector("list", B)
for(j in 1:B)
{
  Z.val[[j]]<- rnorm(n1)
}
T2.null.vals<- vector("list", B)
for(d in 1:B)
{
  for(k in 1:n1)
  {
    T2.null.vals[[d]][k]<- lambdas2[k]*((Z.val[[d]][k])^2-1)
  }
}
T20<- c()
for(d in 1:B)
{
  T20[d]<- sum(T2.null.vals[[d]])
}
alpha = 0.05
T2.quantile<- quantile(T20, probs = 1-alpha)
T2.alt.vals<- vector("list", length(mu))
for(m in 1:length(mu))
{
  for(d in 1:B)
  {
    T2.alt.vals[[m]][d]<- sum(lambdas2*((Z.val[[d]]+mu[m]*b2.vals1)^2-1))
  }
}
power.tau.star<- c()   ## Powers of T.n2 for different values of gamma ##
for(i in 1:length(mu))
{
  power.tau.star[i]<- mean(T2.alt.vals[[i]]>T2.quantile)  
}
power.tau.star


#### dCov ####

kernel.dCov<-function(a,b)  ## kernel function of Tn3 ##
{
  x<-0
  for(i in 1:(n-3))
  {
    for(j in (i+1):(n-2))
    {
      for(k in (j+1):(n-1))
      {
        for(l in (k+1):n)
        {
          x<-x+(1/4)*((abs(a-dt.x[j])+abs(dt.x[k]-dt.x[l])-abs(a-dt.x[k])-abs(dt.x[j]-dt.x[l]))*(abs(b-dt.y[j])+abs(dt.y[k]-dt.y[l])-abs(b-dt.y[k])-abs(dt.y[j]-dt.y[l])))
        }
      }
    }
  }
  return(x/choose(n,4))
}

L1<- matrix(0,n1,n1)  ## Matrix with elements being the values of kernel function of Tn3 approximated at the bivariate marginal quantiles of X and second order difference of Y ##
{
  for(p in 1:n1)
  {
    for(q in 1:n1)
    {
      L1[p,q]<- kernel.dCov(data.pts.x1[p],data.pts.y1[q])
    }
  }
}

lmat3<- (L1+t(L1))/2  ## Symmeytrization of L1 as lmat3, since L1 and lmat3 have same eigenvalues and eigenvectors. ##

E1<- eigen(lmat3) ## eigenvalue decomposition of real symmetric matrix E1 ##

lambdas3<- E1$values ## real eigenvalues of lmat3 ##

V3<- E1$vector ## real eigenvectors of lmat3 ##

V3.t<- t(V3) ## transpose of V3 ##

## computation of the coefficeints a.star_k, k=1,...,n ##
a.star.k<- vector("list", n1)
for(k in 1:n1)
{
  for(i in 1:n1)
  {
    a.star.k[[k]][i]<- h1[1:n1][i]*V3[i,k]*V3.t[k,i]
  }
}
a.star.k1<- c()
for(k in 1:n1)
{
  a.star.k1[k]<- mean(a.star.k[[k]])
}

## Generation of B i.i.d. standard normal variates ##
Z.val<- vector("list", B)
for(j in 1:B)
{
  Z.val[[j]]<- rnorm(n1)
}

## Computation of B values of Tn3 under H0 ##
dCov.null.vals<- vector("list", B)
for(d in 1:B)
{
  for(k in 1:n1)
  {
    dCov.null.vals[[d]][k]<- lambdas3[k]*((Z.val[[d]][k])^2-1)
  }
}
dCov0<- c()
for(d in 1:B)
{
  dCov0[d]<- sum(dCov.null.vals[[d]])
}

dCov.quantile<- quantile(dCov0, probs = 1-alpha)  ## 5% critical value of Tn3 ##

## Computation of values of Tn3 under Hn ##
dCov.alt.vals<- vector("list", length(mu))  
for(m in 1:length(mu))
{
  for(d in 1:B)
  {
    dCov.alt.vals[[m]][d]<- sum(lambdas3*((Z.val[[d]]+mu[m]*a.star.k1)^2-1))
  }
}

Power_of_dCov<- c()  ## Powers of Tn,3 for different values of gamma ##
for(i in 1:length(mu))
{
  Power_of_dCov[i]<- mean(dCov.alt.vals[[i]]>dCov.quantile)
}

Power_of_dCov

cbind.data.frame(power.tau,power.tau.star,Power_of_dCov)
