library(comprehenr)
library(matlib)
library(MASS)
library(glmnet)
library(R.utils)


p=200
OLSE <- function(X, y,D) {
  X = X[,2:p]
  res <- cv.glmnet(X, y, family="gaussian",intercept=(1==1),weights=D,alpha=1)
  lambda_min = res$lambda.min
  res <- glmnet(X, y, family = "gaussian",intercept=(1==1), lambda = 0.1, alpha=1, weights=D)
  beta = as.numeric(coef(res))
  return(beta)
}

estimate_p<-function(n,K){
  return(n/K)
}

estimate_sigma<-function(X,y,D,n,theta_est) {
  a = y - X%*%theta_est
  b = sqrt((t(a)%*%D%*%a)/n)
  return(b)
}

estimate_m<-function(X,y,theta,var,n,p){
  t <- X%*%theta
  m <- rep(0,n)
  rng <- 1:50
  rng_sqrt <- sqrt(rng)
  p_arr <- (1-p)^(rng-1)
  y <- y - t
  for(i in 1:n){
    dist = dnorm(y[i], 0, rng_sqrt*rep(var,50))
    #print(dist*p_arr)
    m[i] = which.max(p_arr*dist)
  }
  d <- rep(1,50)/m
  return(list("m"=m,"d"=d))
}

estimate_m_exp<-function(X,y,theta,var,n,p){
  t <- X%*%theta
  m <- rep(0,n)
  rng <- 1:50
  rng_sqrt <- sqrt(rng)
  p_arr <- (1-p)^(rng-1)
  y <- y - t
  for(i in 1:n){
    dist = dnorm(y[i], 0, rng_sqrt*rep(var,50))
    #print(dist*p_arr)
    m[i] = which.max(p_arr*dist)
  }
  
  d <- rep(1,50)/m
  return(list("m"=m,"d"=d))
}




SGD <- function(X, y,n) {
  return(list("theta"=theta_hat,"sigma"=sigma_hat,"p"=p_))
}


sample_z <- function(x,y,theta_h,alpha,c){
  sset <- c()
  for(zz in 1:50){
    i=40
    while(i){
      i = i-1;
      v <- rnorm(1,0,0.5*alpha)
      z <- x%*%theta_h + log(1+2*v^2+2*v*sqrt(1+v^2))
      if(!is.na(z[1,1])&&z[1,1]>c){
        sset <- c(sset, c(z[1,1]))
        i=0
      }
    }
  }
  if(length(sset)==0){
    return(Inf)
  }
  return(sum(sset)/length(sset))
}

n = 100
theta = c(c(1.5,5,3,2,-0.5),rep(0,p-5))
sigma=1
simu=10
g_p=0.5
iters=100
burnin=20

fres = matrix(0,nrow=simu,ncol=p+2)
for(kk in 1:simu){
  print(kk)
  set.seed(kk)
  
  m =rgeom(n,g_p) + 1
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
  
  
  D = diag(n)
  theta_hat <- OLSE(X,y,rep(1,n))
  p_hat <- 0.5
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
    p_hat <- estimate_p(n, sum(m))
    sigma_hat <- estimate_sigma(X,y,D,n,theta_hat)
    print(c(p_hat, sigma_hat, theta_hat[1:5]))
    if(i>burnin){
      estt[i-burnin,] = c(p_hat, sigma_hat, theta_hat)
    }
  }
  #print(estimate_m(X,y,theta_hat,sigma_hat,n,p
  Pzero=apply(estt,2,nnzero)/(iters-burnin)
  zvalue=(Pzero-0.5)*(iters-burnin)/0.5
  pvalue=pnorm(zvalue)
  adj.pv=p.adjust(pvalue, "BH")    #Control false discovery rate;
  estt[,which(adj.pv<=0.05)] <- 0
  printf("final values in %dth iteration\n",kk)
  fres[kk, ] = apply(estt,2,sum)/(iters-burnin)
  print(fres[kk,])
}

print(fres)
fres <- matrix(unlist(fres), nrow=simu, ncol=p+1)
print("Mean")
print(apply(fres,2,mean))
print("Variance")
print(apply(fres,2,sd))
print("Over Estimation")
print(apply(fres, 1, function(c)sum(c!=0))-6)
fmean <- apply(fres,2,mean)

