library(comprehenr)
library(matlib)
library(MASS)
library(glmnet)
library(R.utils)
library(R.matlab)

p=487
OLSE <- function(X, y,D) {
  X<-X[,2:p]
  res <- cv.glmnet(X, y, family="gaussian",intercept=(1==1),weights=D,alpha=1)
  lambda_min = res$lambda.min
  res <- glmnet(X, y, family = "gaussian",intercept=(1==1), lambda = lambda_min, alpha=1, weights=D)
  beta = as.numeric(coef(res))
  return(beta)
}

estimate_p <-function(n,K){
  return(n/K)
}

estimate_sigma <-function(X,y,D,n,theta_est) {
  a = y - X%*%theta_est
  b = sqrt((t(a)%*%D%*%a)/n)
  return(b)
}

estimate_m<-function(X,y,theta,var,n,p){
  mx <- 100
  t <- X%*%theta
  m <- rep(0,n)
  d <- rep(0,n)
  rng <- 1:mx
  p_arr <- (1-p)^(rng-1)
  rng_rt <- sqrt(rng)
  y <- y - t
  for(i in 1:n){
    dist = dnorm(y[i], 0, rng_rt*rep(var,mx))
    m[i] = sum(p_arr*rng*dist)/sum(p_arr*dist)
    d[i] = sum((p_arr*dist)/rng)/sum(p_arr*dist)
  }
  return(list("m"=m,"d"=d))
}

n = 180
crvalue = 3
theta = c(c(1.5,5,3,2,-0.5), rep(0,p-5))
sigma=1
simu=1
g_p=0.50
iters=200
burnin=40
fres = matrix(0,nrow=simu,ncol=p+2)

for(kk in 1:simu){
  print(kk)
  set.seed(kk)
  
  m <-readMat('VesselData.mat')
  X<-matrix(nrow=n,ncol=487)
  X[,2:487]<-m$Xdata[,15:500]
  X[,1]<-rep(1,n)
  y<-m$Ydata[,11]
  
  d = 1/(rgeom(n,0.5)+1)
  D = diag(d)
  theta_hat <- OLSE(X,y,d)
  p_hat <- 0.5
  print(theta_hat)
  sigma_hat <- estimate_sigma(X,y,D,n,theta_hat)
  print(sigma_hat)
  estt = matrix(0,nrow=iters-burnin, ncol=p+2)
  for(i in 1:iters){
    res <- estimate_m(X,y,theta_hat,sigma_hat,n,p_hat)
    m <- res$m 
    d <- res$d
    for(j in 1:n){
      D[j,j] = d[j]
    }
    theta_hat <- OLSE(X,y,d)
    p_hat     <- estimate_p(n, sum(m))
    sigma_hat <- estimate_sigma(X,y,D,n,theta_hat)
    print(c(p_hat, sigma_hat, theta_hat[1:5]))
    if(i>burnin){
      estt[i-burnin,] = c(p_hat, sigma_hat, theta_hat)
      
    }
  }
  fres[kk, ] = apply(estt,2,mean)
}
print(fres)
rat <- rep(0,simu)
for(i in 1:simu) {
  rat[i]=fres[i,2]/sqrt(fres[i,1])
}
print(apply(fres, 2, mean))
print(apply(fres, 2, sd))

print(mean(rat))
print(sd(rat))
print(1/sqrt(g_p))