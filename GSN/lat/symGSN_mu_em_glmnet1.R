library(comprehenr)
library(matlib)
library(MASS)
library(glmnet)
library(R.utils)

p<-200

OLSE1 <- function(X, y) {
  res <- cv.glmnet(X, y, family="gaussian",intercept=(1==1),alpha=1)
  lambda_min <- res$lambda.min
  res <- glmnet(X, y, family = "gaussian",intercept=(1==1), lambda = lambda_min, alpha=1)
  # res <- glmnet(X, y, family = "gaussian",intercept=(1==1), lambda = 0, alpha=1)
  return(res)
}

OLSE2 <- function(X, y) {
  res <- cv.glmnet(X, y, family="gaussian",intercept=(1==0),alpha=0)
  lambda_min <- res$lambda.min
  res <- glmnet(X, y, family = "gaussian",intercept=(1==0), lambda = lambda_min, alpha=0)
  beta <- as.numeric(coef(res))
  return(beta[1:p+1])
}

OLSE <- function(X, y,d,m,mu,n) {
  print(d)
  y <- y*sqrt(d)
  X <- X*sqrt(d)
  y <- y - rep(mu,n)*(1/sqrt(d))
  res <- cv.glmnet(X, y, family="gaussian",intercept=(1==0),alpha=1)
  lambda_min <- res$lambda.min
  res <- glmnet(X, y, family = "gaussian",intercept=(1==0), lambda = lambda_min, alpha=1)
  # res <- glmnet(X, y, family = "gaussian",intercept=(1==0), lambda = 0, alpha=1)
  beta <- as.numeric(coef(res))
  return(beta[1:p+1])
}

estimate_p <-function(n,K){
  return(n/K)
}

estimate_mu <- function(y, X, theta, m) {
  a <- y - X%*%theta 
  return(sum(a)/sum(m))
}

estimate_sigma <-function(X,y,d,m,n,theta_est,mu) {
  a <- y - X%*%theta_est 
  res <- sum(a*a*d) - 2*sum(a*rep(mu,n)) + sum(m*rep(mu^2, n)) 
  return(sqrt(res/n))
}

estimate_m_conf <- function(X,y,theta,var,n,p,mu,m) {
  t <- X%*%theta
  y <- y - t
  return(dnorm(y,m*mu, sqrt(m)*rep(var,n))*(1-p)^(m-1))
}

estimate_m_mm<-function(X,y,theta,var,n,p,mu){
  mx<-50
  t <- X%*%theta
  m <- rep(0,n)
  rng <- 1:mx
  rng_sqrt <- sqrt(rng)
  p_arr <- (1-p)^(rng-1)
  y <- y - t
  for(i in 1:n){
    dist <- dnorm(y[i], rep(mu,mx)*rng, rng_sqrt*rep(var,mx))
    #print(dist*p_arr)
    m[i] <- which.max(p_arr*dist)
  }
  d <- rep(1,n)/m
  return(list("m"=m,"d"=d))
}


estimate_m<-function(X,y,theta,var,n,p,mu){
  mx <- 100
  t <- X%*%theta
  m <- rep(0,n)
  d <- rep(0,n)
  l <- rep(0,n)
  rng <- 1:mx
  p_arr <- (1-p)^(rng-1)
  rng_rt <- sqrt(rng)
  y <- y - t
  for(i in 1:n){
    dist <- dnorm(y[i], rng*rep(mu,mx), rng_rt*rep(var,mx))
    m[i] <- sum(p_arr*rng*dist)/sum(p_arr*dist)
    d[i] <- sum((p_arr*dist)/rng)/sum(p_arr*dist)
    l[i] <- sum((p_arr*dist)*log(rng))/sum(p_arr*dist)
  }
  return(list("m"=m,"d"=d,"l"=l))
}

ll <- function(X,y,theta,mu,sigma,p,m,d,l,n) {
  sm <- n*log(p) + (sum(m)-n)*log(1-p) - n*log(sigma) - 0.5*sum(l)
  a <- y - X%*%theta
  aa <- sum(a*a*d) - 2*sum(a*rep(mu,n)) + sum(m*rep(mu^2, n)) 
  sm <- sm - aa/(2*(sigma^2))
  return(sm)
}

log_likelihood <- function(X,y,theta,p,mu,n){
  sigma <- 1
  print(y)
  res <- estimate_m(X,y,theta,sigma,n,p,mu)
  m <- res$m
  d <- res$d
  for(i in 1:50){
    theta <- OLSE(X,y,d,m,mu,n)
    sigma <- estimate_sigma(X,y,d,m,n,theta,mu)
    res <- estimate_m(X,y,theta,sigma,n,p,mu)
    m <- res$m
    d <- res$d
    l<-res$l
  }
  #m <- res$m
  #d <- res$d
  return(ll(X,y,theta,mu,sigma,p,m,d,l,n))
}

initial_estimate<-function(X,y,n){
  res <- OLSE1(X,y)
  theta_full <- as.numeric(coef(res))
  z <- y-X%*%theta_full[1:p+1]
  p_r <- seq(0.05,0.9,by=0.05)
  mu_r <- rep(theta_full[1],18)*p_r
  ll <- rep(0,18)
  for(i in 1:18){
    ll[i] <- log_likelihood(X,y,theta_full[1:p+1],p_r[i],mu_r[i],n)
    print(ll[i])
  }
  #print(ll)
  i <- which.max(ll)
  p <- p_r[i]
  mu <- mu_r[i]
  m <- rgeom(n,p)+1
  d <- 1/m
  theta <- OLSE(X,y,d,m,mu,n)
  sigma <- estimate_sigma(X,y,d,m,n,theta,mu)
  res <- estimate_m(X,y,theta,sigma,n,p,mu)
  m <- res$m
  d <- res$d
  theta <- OLSE(X,y,d,m,mu,n)
  sigma <- estimate_sigma(X,y,d,m,n,theta,mu)
  x <- list("mu"=mu, "theta"=theta, "sigma"=sigma, "p"=p, "m"=m, "z"=z, "ll"=ll)
  return(x)
}

initial_estimate2<-function(X,y,mu,n,mact) {
  D=rep(1,n)
  m=rep(0,n)
  theta <- OLSE(X,y,D,m,mu,n)
  print(theta)
  z<- (y - X%*%theta)
  min=0
  pos=1
  D= diag(D)
  sigma=estimate_sigma(X,y,D,m,n,theta,mu)
  mg = rep(0,n)
  for(i in 1:n){
    pdf = dnorm(y[i],rep(mu,n)*(1:20),rep(sigma,n)*sqrt((1:20)))
    pos= which.max(pdf)
    mg[i]=pos
  }
  m = z/mu
  comp = matrix(nrow=n, ncol=5)
  comp[,1]=mact
  comp[,2]=m
  comp[,3]=m-mact
  comp[,4]=mg
  comp[,5]=mg-mact
  hist(comp[,5])
  print(comp)
  print(min(comp[,5]))
  print(max(comp[,5]))
  print(min(comp[,3]))
  print(max(comp[,3]))
  print(1/mean(mg))
  print(sum(comp[,5]*comp[,5]))
  print(skewness(y))
  m <- mg
  d = 1/m
  theta_hat <- OLSE(X,y,d,m,mu,n)
  sigma_hat <- estimate_sigma(X,y,d,m,n,theta_hat,mu)
  p_hat <- estimate_p(n, sum(m))
  x <- list("theta"=theta_hat, "sigma"=sigma_hat, "p"=p_hat)
  return(x)
}

n <- 100
crvalue <- 3
theta <- c(c(1.5,5,3,2,-0.5), rep(0,p-5))
mu<-10
sigma<-1
simu<-10
g_p<-0.3
iters<-1000
burnin<-200
fres <- matrix(0,nrow=simu,ncol=p+3)

for(kk in 1:simu){
  print(kk)
  set.seed(kk)
  m <-rgeom(n,g_p) +1
  X<-matrix(nrow = p, ncol = n)      #generate X
  for(i in 1:p)
  {                                   
    # X[i+1,]=rnorm(n,i,2)    #Normalize X so that E(X)=0 and comptation for the intercept is stable;
    X[i,]<-rnorm(n,0,2)
  }
  X <- t(X)
  t<-X%*%theta
  y <- rep(1,n)
  for(i in 1:n){
    y[i] <- t[i] + rnorm(1,mu*m[i], sqrt(m[i])*sigma)
  }
  
  res <- initial_estimate(X,y,n)
  theta_hat <- res$theta
  mu_hat <- res$mu
  sigma_hat <- res$sigma
  p_hat <- res$p
  estt <- matrix(0,nrow=iters-burnin, ncol=p+3)
  for(i in 1:iters){
    print(c(p_hat, sigma_hat, mu_hat,theta_hat[1:5]))
    res <- estimate_m(X,y,theta_hat,sigma_hat,n,p_hat,mu_hat)
    m <- res$m 
    d <- res$d
    # if(i<100) {
    #   mp <- estimate_m_conf(X,y,theta_hat,sigma_hat,n,p_hat, mu_hat,n)
    #   idx <- order(mp)[51:100]
    #   Xp <- X[idx, ]
    #   yp <- y[idx]
    #   m1 <- round(m[idx])
    #   d1 <- 1/m1
    #   theta_hat <- OLSE(Xp,yp,d1,m1,mu_hat,50)
    # }else {
    # }
    
    theta_hat <- OLSE(X,y,d,m,mu_hat,n)
    p_hat     <- estimate_p(n, sum(m))
    mu_hat    <- estimate_mu(y,X,theta_hat,m)
    sigma_hat <- estimate_sigma(X,y,d,m,n,theta_hat,mu_hat)
    if(i>burnin){
      estt[i-burnin,] <- c(p_hat, sigma_hat, mu_hat,theta_hat)
    }
  }
  fres[kk, ] <- apply(estt,2,mean)
  print(fres[kk,1:10])
}
print(fres)
print(apply(fres, 2, mean))
print(apply(fres, 2, sd))

file_str<-sprintf("symGSN_mu_em_glmnet_mu=%f_p=%f_sigma=%f.rds",mu,g_p,sigma)
saveRDS(fres, file_str)
# 
# res<-initial_estimate(X,y,n)
# plot(res$ll)