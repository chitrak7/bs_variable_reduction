library(comprehenr)
library(matlib)
library(MASS)
library(glmnet)
library(R.utils)

p=200

OLSE1 <- function(X, y) {
  res <- cv.glmnet(X, y, family="gaussian",intercept=(1==1),alpha=1)
  lambda_min = res$lambda.min
  res <- glmnet(X, y, family = "gaussian",intercept=(1==1), lambda = lambda_min, alpha=1)
  print(lambda_min)
  print(as.numeric(coef(res)))
  return(res)
}

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

log_likelihood <- function(X,y,p,mu,n){
  m = rgeom(n,p)+1
  d = 1/m
  D = diag(d)
  theta = OLSE(X,y,d,m,mu,n)
  sigma = estimate_sigma(X,y,D,n,theta,mu)
  res <- estimate_m(X,y,theta,sigma,n,p,mu)
  m = res$m
  D = diag(res$d)
  sm = n*log(p) + sum(m-rep(1,n))*log(1-p) - n*log(sigma) - 0.5*sum(log(m))
  y = y - X%*%theta - rep(mu,n)*m
  sm = sm - ((t(y)%*%D%*%y))/(2*(sigma^2))
  return(sm)
}

initial_estimate<-function(X,y,n){
  res <- OLSE1(X,y)
  theta_full <- as.numeric(coef(res))
  p_r = seq(0.05,0.91,by=0.05)
  mu_r = rep(theta_full[1],18)*p_r
  prd = predict.glmnet(res,X)
  ll = rep(0,18)
  for(i in 1:18){
    ll[i] = log_likelihood(X,y,p_r[i],mu_r[i],n)
    print(i)
  }
  print(ll)
  i <- which.max(ll)g_p
  p = p_r[i]
  mu = mu_r[i]
  m = rgeom(n,p)+1
  d = 1/m
  D = diag(d)
  theta = OLSE(X,y,d,m,mu,n)
  sigma = estimate_sigma(X,y,D,n,theta,mu)
  res <- estimate_m(X,y,theta,sigma,n,p,mu)
  m <- res$m
  d <- res$d
  D <- diag(d)
  theta <- OLSE(X,y,d,m,mu,n)
  sigma = estimate_sigma(X,y,D,n,theta,mu)
  
  print(theta)
  print(mu)
  print(p)
  print(sigma)
  print(m)
  x <- list("mu"=mu, "theta"=theta, "sigma"=sigma, "p"=p)
  #return(x)
}

initial_estimate2<-function(X,y,mu,n) {
  D=rep(1,n)
  m=rep(0,n)
  theta <- OLSE(X,y,D,m,mu,n)
  z<- y - X%*%theta
  print(z)
  m = z/mu
  D = 1/m
  theta_hat <- OLSE(X,y,D,m,mu,n)
  sigma_hat <- estimate_sigma(X,y,D,n,theta_hat,mu)
  p_hat <- estimate_p(n, sum(m))
  x <- list("theta"=theta_hat, "sigma"=sigma_hat, "p"=p_r[i])
  return(x)
}

n = 100
crvalue = 3
theta = c(c(1.5,5,3,2,-0.5), rep(0,p-5))
mu=10
sigma=1
simu=10
g_p=0.25
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
  
  sample(1:5,10,replace=T)
  D = diag(1/m)
  mu_hat=10
  theta_hat <- OLSE(X,y,m,m,mu_hat,n)
  p_hat <- 0.5
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

file_str=sprintf("symGSN_mu_em_glmnet_mu=%f_p=%f_sigma=%f.rds",mu,g_p,sigma)
saveRDS(fres, file_str)