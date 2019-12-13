
p<-200
n <- 100
crvalue <- 3
theta <- c(c(1.5,5,3,2,-0.5), rep(0,p-5))
args<-commandArgs(trailingOnly = TRUE)
args <- as.numeric(args)
args<-c(1,1,0.5,1)
mu=args[1]
sigma=args[2]
g_p=args[3]
alpha_s=args[4]
simu<-1000
iters<-200
burnin<-40

mes = sprintf("symGSN with em and glmnet for mu=%f, p=%f, sigma=%f, alpha=%f",mu,g_p,sigma,alpha_s)
print(mes)

library(matlib)
library(MASS)
library(glmnet)
library(doParallel)
library(moments)

fres <- matrix(0,nrow=simu,ncol=p+3+n)
# E-step for m, 1/m, log m
estimate_m<-function(X,y,theta,sd,n,p,mu){
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
    dist <- dnorm(y[i], rng*rep(mu,mx), rng_rt*rep(sd,mx))
    m[i] <- sum(p_arr*rng*dist)/sum(p_arr*dist)
    d[i] <- sum((p_arr*dist)/rng)/sum(p_arr*dist)
    l[i] <- sum((p_arr*dist)*log(rng))/sum(p_arr*dist)
  }
  return(list("m"=m,"d"=d,"l"=l))
}
# M-step for p
estimate_p <-function(n,K){
  return(n/K)
}
# M-step for mu
estimate_mu <- function(y, X, theta, m) {
  a <- y - X%*%theta 
  return(sum(a)/sum(m))
}
# M-step for sigma
estimate_sigma <-function(X,y,d,m,n,theta_est,mu) {
  a <- y - X%*%theta_est 
  res <- sum(a*a*d) - 2*sum(a*rep(mu,n)) + sum(m*rep(mu^2, n)) 
  return(sqrt(res/n))
}
# M-step for theta
OLSE <- function(X, y,d,m,mu,n) {
  y <- y*sqrt(d)
  X <- X*sqrt(d)
  y <- y - rep(mu,n)*(1/sqrt(d))
  res <- cv.glmnet(X, y, family="gaussian",intercept=(1==0),alpha=1)
  lambda_min <- res$lambda.min
  res <- glmnet(X, y, family = "gaussian",intercept=(1==0), lambda = lambda_min, alpha=alpha_s)
  beta <- as.numeric(coef(res))
  return(beta[1:p+1])
}
#OLSE with intercept to get beta for mu/p
OLSE1 <- function(X, y) {
  res <- cv.glmnet(X, y, family="gaussian",intercept=(1==1),alpha=1)
  lambda_min <- res$lambda.min
  res <- glmnet(X, y, family = "gaussian",intercept=(1==1), lambda = lambda_min, alpha=alpha_s)
  return(res)
}
# Calculates expected log-likelihood
ll <- function(X,y,theta,mu,sigma,p,m,d,l,n) {
  sm <- n*log(p) + (sum(m)-n)*log(1-p) - n*log(sigma) - 0.5*sum(l)
  a <- y - X%*%theta
  aa <- sum(a*a*d) - 2*sum(a*rep(mu,n)) + sum(m*rep(mu^2, n)) 
  sm <- sm - aa/(2*(sigma^2))
  return(sm)
}
# Calculates log-likelihood
ll2 <- function(X,y,theta,mu,sigma,p,m,n) {
  sm <- n*log(p) + (sum(m)-n)*log(1-p) - n*log(sigma) - 0.5*sum(log(m))
  a <- y - X%*%theta - rep(mu,n)*m
  aa <- (a^2)/m
  sm <- sm - sum(aa)/(2*(sigma^2))
  return(sm)
}

log_likelihood <- function(X,y,theta,p,mu,n){
  sigma <- sqrt(sum((y-X%*%theta)^2)*(p/n))
  #print(sigma)
  z<- (y - X%*%theta)
  pos=0
  mg = rep(0,n)
  for(i in 1:n){
    pdf = dnorm(z[i], mu*(1:20),sigma*sqrt((1:20)))
    #print(pdf)
    pos= which.max(pdf)
    mg[i]=pos
  }
  #print(mg)
  return(ll2(X,y,theta,mu,sigma,p,mg,n))
}

initial_estimate<-function(X,y,n){
  res <- OLSE1(X,y)
  theta_full <- as.numeric(coef(res))
  #print("ok1")
  z <- y-X%*%theta_full[1:p+1]
  #print("ok2")
  p_r <- seq(0.06,0.9,by=0.01)
  mu_r <- rep(theta_full[1],85)*p_r
  ll <- rep(0,85)
  for(i in 1:85){
    ll[i] <- log_likelihood(X,y,theta_full[1:p+1],p_r[i],mu_r[i],n)
  }
  #print(ll)
  #print("ok3")
  i <- which.max(ll)
  theta <- theta_full[1:p+1]
  p <- p_r[i]
  mu <- mu_r[i]
  sigma <- sqrt(sum(z^2)*(p/n))
  x <- list("mu"=mu, "theta"=theta, "sigma"=sigma, "p"=p)
  return(x)
}

rsimu <- function(kk){
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
  m_hat <- res$m
  estt <- matrix(0,nrow=iters-burnin, ncol=p+3)
  for(i in 1:iters){
    res <- estimate_m(X,y,theta_hat,sigma_hat,n,p_hat,mu_hat)
    m <- res$m 
    d <- res$d
    l <- res$l
    theta_hat <- OLSE(X,y,d,m,mu_hat,n)
    p_hat     <- estimate_p(n, sum(m))
    mu_hat    <- estimate_mu(y,X,theta_hat,m)
    sigma_hat <- estimate_sigma(X,y,d,m,n,theta_hat,mu_hat)
    if(i>burnin){
      estt[i-burnin,] <- c(p_hat, sigma_hat, mu_hat,theta_hat)
    }
  }
  Pzero=apply(estt,2,nnzero)/(iters-burnin)
  zvalue=(Pzero-0.5)*(iters-burnin)/0.5
  pvalue=pnorm(zvalue)
  adj.pv=p.adjust(pvalue, "BH")    #Control false discovery rate;
  estt[,which(adj.pv<=0.05)] <- 0
  f<-apply(estt,2,sum)/(iters-burnin)
  return(f)
}

nc <- detectCores()-1
cl <- makeCluster(nc)
st <- Sys.time()
fres <- foreach(kk=1:simu, .packages = c('glmnet','matlib','MASS')) %dopar% rsimu(kk)
print(Sys.time()-st)
stopCluster(cl)
fres <- matrix(unlist(fres), nrow=simu, byrow = TRUE)
print(apply(fres, 2, mean))
print(apply(fres, 2, sd))
file_str<-sprintf("symGSN_mu_em_glmnet_mu=%f_p=%f_sigma=%f.rds",mu,g_p,sigma)
saveRDS(fres, file_str)