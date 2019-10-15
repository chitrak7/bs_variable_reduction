library(comprehenr)
library(matlib)
library(MASS)
library(glmnet)
library(R.utils)


OLSE <- function(X, y,D,m,mu,n) {
  y = y - m*rep(mu,n)
  a = t(X)%*%D%*%X
  b = inv(a)%*%t(X)%*%D%*%y
  return(b)
}

estimate_p <-function(n,K){
  return(n/K)
}

estimate_mu <- function(y, X, theta, m) {
  a = y - X%*%theta 
  return(sum(a)/sum(m))
}

estimate_sigma <-function(X,y,D,n,theta_est,mu) {
  a = y - X%*%theta_est -m*rep(mu,n)
  b = sqrt((t(a)%*%D%*%a)/n)
  return(b)
}

estimate_m<-function(X,y,theta,var,n,p,mu){
  mx <- 100
  t <- X%*%theta
  m <- rep(0,n)
  d <- rep(0,n)
  rng <- 1:mx
  p_arr <- (1-p)^(rng-1)
  rng_rt <- sqrt(rng)
  y <- y - t
  for(i in 1:n){
    dist = dnorm(y[i], rng*rep(mu,mx), rng_rt*rep(var,mx))
    m[i] = sum(p_arr*rng*dist)/sum(p_arr*dist)
    d[i] = sum((p_arr*dist)/rng)/sum(p_arr*dist)
  }
  return(list("m"=m,"d"=d))
}

p=5
n = 100
crvalue = 3
theta = c(c(15,5,3,2,-0.5), rep(0,p-5))
sigma=1
simu=10
g_p=0.25
iters=1000
burnin=200
fres = matrix(0,nrow=simu,ncol=p+2)

for(kk in 1:simu){
  print(kk)
  set.seed(kk)
  m =rgeom(n,g_p) +1
  em = sum(m)/n
  X=matrix(nrow = p, ncol = n)      #generate X
  X[1,]=rep(1,n)
  for(i in 2:p)
  {                                   
    # X[i+1,]=rnorm(n,i,2)    #Normalize X so that E(X)=0 and comptation for the intercept is stable;
    X[i,]=rnorm(n,0,2)
  }
  X <- t(X)
  t<-X%*%theta
  y <- rep(1,n)
  for(i in 1:n){
    y[i] = t[i] + rnorm(1,0, sqrt(m[i])*sigma)
  }
  X = X[,2:p]
  D = diag(n)
  m=rep(1,n)
  mu_hat=10
  theta_hat <- OLSE(X,y,D,m,mu_hat,n)
  print(theta_hat)
  p_hat <- 0.5
  sigma_hat <- estimate_sigma(X,y,D,n,theta_hat,mu_hat)
  print(sigma_hat)
  estt = matrix(0,nrow=iters-burnin, ncol=p+2)
  for(i in 1:iters){
    res <- estimate_m(X,y,theta_hat,sigma_hat,n,p_hat,mu_hat)
    m <- res$m 
    d <- res$d
    for(j in 1:n){
      D[j,j] = d[j]
    }
    theta_hat <- OLSE(X,y,D,m,mu_hat,n)
    p_hat     <- estimate_p(n, sum(m))
    mu_hat    <- estimate_mu(y,X,theta_hat,m)
    sigma_hat <- estimate_sigma(X,y,D,n,theta_hat,mu_hat)
    if(i<10)
    print(c(p_hat, sigma_hat, mu_hat/p_hat,theta_hat))
    if(i>burnin){
      estt[i-burnin,] = c(p_hat, sigma_hat, mu_hat/p_hat,theta_hat)
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