# args<-commandArgs(trailingOnly = TRUE)
# args <- as.numeric(args)
# mu=args[1]
# sigma=args[2]
# g_p=args[3]
# file_str=sprintf("symGSN_mu_em_glmnet_mu=%f_p=%f_sigma=%f.txt",mu,g_p,sigma)
# obj_str=sprintf("symGSN_mu_em_glmnet_mu=%f_p=%f_sigma=%f.rds",mu,g_p,sigma)
# mes = sprintf("symGSN with em and glmnet for mu=%f, p=%f, sigma=%f",mu,g_p,sigma)
# sink(file_str)
# print(mes)
mu=10
sigma=1
g_p=0.25
library(comprehenr)
library(matlib)
library(MASS)
library(glmnet)
library(R.utils)

p=200
OLSE <- function(X, y,D,m,mu,n) {
  y = y - m*rep(mu,n)
  res <- cv.glmnet(X, y, family="gaussian",intercept=(1==0),weights=D,alpha=1)
  lambda_min = res$lambda.min
  res <- glmnet(X, y, family = "gaussian",intercept=(1==0), lambda = lambda_min, alpha=1, weights=D)
  beta = as.numeric(coef(res))
  return(beta[1:p+1])
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
n = 100
crvalue = 3
theta = c(c(1.5,5,3,2,-0.5), rep(0,p-5))

simu=100
iters=1000
burnin=200
fres = matrix(0,nrow=simu,ncol=p+3)

for(kk in 1:simu){
  print(kk)
  set.seed(kk)
  m =rgeom(n,g_p) +1
  em = sum(m)/n
  X=matrix(nrow = p, ncol = n)      #generate X
  for(i in 1:p)
  {                                   
    # X[i+1,]=rnorm(n,i,2)    #Normalize X so that E(X)=0 and comptation for the intercept is stable;
    X[i,]=rnorm(n,0,2)
  }
  X <- t(X)
  t<-X%*%theta
  y <- rep(1,n)
  for(i in 1:n){
    y[i] = t[i] + rnorm(1,mu*m[i], sqrt(m[i])*sigma)
  }
  
  p_hat <- 0.5
  #m=rgeom(n,p_hat)+1
  d = 1/m
  D = diag(d)
  mu_hat=1
  theta_hat <- OLSE(X,y,d,m,mu_hat,n)
  sigma_hat <- estimate_sigma(X,y,D,n,theta_hat,mu_hat)
  print(sigma_hat)
  estt = matrix(0,nrow=iters-burnin, ncol=p+3)
  for(i in 1:iters){
    res <- estimate_m(X,y,theta_hat,sigma_hat,n,p_hat,mu_hat)
    m <- res$m 
    d <- res$d
    for(j in 1:n){
      D[j,j] = d[j]
    }
    theta_hat <- OLSE(X,y,d,m,mu_hat,n)
    p_hat     <- estimate_p(n, sum(m))
    mu_hat    <- estimate_mu(y,X,theta_hat,m)
    sigma_hat <- estimate_sigma(X,y,D,n,theta_hat,mu_hat)
    if(i>burnin){
      estt[i-burnin,] = c(p_hat, sigma_hat, mu_hat,theta_hat)
    }
  }
  fres[kk, ] = apply(estt,2,mean)
  print(fres[kk,1:10])
}
fmean <- apply(fres, 2, mean)
fsd  <- apply(fres, 2, sd)
f0 <- c(g_p, sigma, mu, theta)
fbias <- fmean - f0
print("Mean")
print(fmean)
print("Standard Deviation")
print(fsd)
print("Bias")
print(fbias)
print("Mean Squared Error")
for(i in 1:simu){
  fres[i,] = fres[i,] - f0
}
fres = fres^2
mse = apply(fres, 2, mean)
print(mse)
saveRDS(fres, obj_str)