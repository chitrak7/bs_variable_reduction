library(comprehenr)
library(matlib)
library(MASS)
library(glmnet)
library(R.utils)


OLSE <- function(X, y,M,D) {
  a = t(X)%*%M%*%D%*%M%*%X
  b = inv(a)%*%t(X)%*%M%*%D%*%y
  return(b)
}

estimate_p <-function(n,K){
  return(n/K)
}

estimate_sigma <-function(X,y,M,D,n,theta_est) {
  a = y - M%*%X%*%theta_est
  b = sqrt((t(a)%*%D%*%a)/n)
  return(b)
}

estimate_m<-function(X,y,theta,var,n,p){
  t <- X%*%theta
  m <- rep(0,n)
  d <- rep(0,n)
  p_arr <- (1-p)^(m-1)
  ttt <- 1:100
  tttsq <- sqrt(ttt)
  tttcb <- ttt*tttsq
  for(i in 1:n){
    dist = dnorm(y[i], ttt*t[i], tttsq*rep(var,n))
    m[i] = sum(p_arr*tttsq*dist)/sum((p_arr*dist)/tttsq)
    d[i] = sum((p_arr*dist)/tttcb)/sum((p_arr*dist)/tttsq)
  }
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

p=80
n = 100
crvalue = 3
theta = c(c(1.5,5,3,2,-0.5),rep(0,p-5))
sigma=1
simu=100
g_p=0.10
fres = matrix(0,nrow=simu,ncol=206)
iters=1000
burnin=200

fres = matrix(0,nrow=simu,ncol=p+2)
for(kk in 1:simu){
  print(kk)
  set.seed(kk)
  m =rgeom(n,g_p) +1
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
    y[i] = rnorm(1,m[i]*t[i], m[i]*sigma)
  }
  
  
  
  M = diag(n)
  D = diag(n)
  theta_hat <- OLSE(X,y,M,D)
  p_hat <- 0.5
  sigma_hat <- estimate_sigma(X,y,M,D,n,theta_hat)
  print(sigma_hat)
  estt = matrix(0,nrow=iters-burnin, ncol=p+2)
  for(i in 1:iters){
    res <- estimate_m(X,y,theta_hat,sigma_hat,n,p_hat)
    m <- res$m 
    d <- res$d
    for(j in 1:n){
      M[j,j] = m[j]
      D[j,j] = d[j]
    }
    theta_hat <- OLSE(X,y,M,D)
    p_hat <- estimate_p(n, sum(m))
    sigma_hat <- estimate_sigma(X,y,M,D,n,theta_hat)
    if(i>burnin){
      estt[i-burnin,] = c(p_hat, sigma_hat, theta_hat)
    }
  }
  fres[kk, ] = apply(estt,2,sum)/(iters-burnin)
}
print(fres)
print(apply(fres, 2, mean))
print(apply(fres, 2, sd))
