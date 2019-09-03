library(comprehenr)
library(matlib)
library(MASS)
library(glmnet)
library(R.utils)

p = 205
OLSE <- function(X, y) {
  X <-X[,2:205] 
  r <- seq(0,1,by=0.02)
  min_CVE <- 10000000
  lambda_min <- 1
  alpha_min <- 1
  
  res <- cv.glmnet(X, y, family="gaussian",penalty.factor=c(0,rep(1,p-1)))
  lambda_min = res$lambda.min
  res <- glmnet(X, y, family = "gaussian", lambda = lambda_min, alpha=alpha_min,penalty.factor=c(0,rep(1,p-1)))
  return(as.numeric(coef(res)))
}

erf <- function(x) {
  return((pnorm(x*sqrt(2))-0.5)*2)
}

alpha_mle  <- function(y,z,n){
  return(sqrt((4/n)*sum(to_vec(for(i in 1:n) sinh(0.5*(y[i]-z[i]))^2))))
}

C <- function(x) {
  return(2 + 4/(x^2) - ((sqrt((2*pi)/(x^2))))*(1-erf(sqrt(2/(x^2))))*exp(2/(x^2)))
}

R <- function(y,z,n,alpha){
  return(to_vec(for(i in 1:n)(4/alpha^2)*(cosh(0.5*(y[i]-z[i]))*sinh(0.5*(y[i]-z[i]))) 
                - (sinh(0.5*(y[i]-z[i]))/cosh(0.5*(y[i]-z[i])))))
}

sum_z_sq <- function(y,z,n,alphah) {
  sum = sum(to_vec(for(i in 1:n)((4/alphah^2)*(sinh(0.5*(y[i]-z[i]))^2))))/n
  return(sum+1)
}

SGD <- function(X, y, k, theta_hat, alpha_hat) {
  theta_hat <- OLSE(X,y)
  alpha_hat <- alpha_mle(y, X%*%theta_hat, length(y))
  
  for(i in 1:k){
    y_t <- X%*%theta_hat + (2/C(alpha_hat))*R(y,X%*%theta_hat, length(y),alpha_hat)
    theta_hat <- OLSE(X,y_t)
    ttt <- sum_z_sq(y,X%*%theta_hat, length(y),alpha_hat)
    if(ttt>=4 && alpha_hat>=3){
      print("broken")
      break
    }
    alpha_hat <- (ttt*alpha_hat)/2
  }
  return(list("theta"=theta_hat,"alpha"=alpha_hat))
}

initial_estimate <-function(X,y,delta) {
  idx_obs <- to_vec(for(i in 1:n) if(delta[i]==0) i)
  X_obs <- X[idx_obs,]
  y_obs <- y[idx_obs]
  return(SGD(X_obs, y_obs,300,0,0))
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

n = 50
alpha = 0.5
crvalue = 10
theta = c(c(1.5,5,3,2,-0.5), seq(0,0,length.out = 200))
fres = matrix(nrow=0,ncol=206)
simu=10
for(kk in 1:simu){
  set.seed(kk)
  X=matrix(nrow = p, ncol = n)      #generate X
  X[1,]=rep(1,n)
  for(i in 2:p)
  {                                   
    # X[i+1,]=rnorm(n,i,2)    #Normalize X so that E(X)=0 and comptation for the intercept is stable;
    X[i,]=rnorm(n,0,2)
  }
  X <- t(X)
  t <- exp(X%*%theta)
  v <- rnorm(n,0,0.5*alpha)
  y <- log(to_vec(for(i in 1:n) (t[i]*(1+2*(v[i]^2)+2*v[i]*sqrt(1+v[i]^2)))))
  c <- crvalue + sqrt(2)*rnorm(n, 0, 1)
  delta <- seq(0,0,length.out = n)

  for (i in 1:n) {
    if(y[i]>=c[i]){
      y[i]=c[i]
      delta[i]=1
    }
  }
  res <- initial_estimate(X,y,delta)
  
  
  for (i in 1:300) {
    idx <- c()
    for(j in 1:n){
      if(delta[j]==1){
        sample<-sample_z(X[j,],y[j],res[["theta"]],res[["alpha"]],c[j])
        if(!is.infinite(sample)){
          y[j] = sample
          idx <-c(idx,c(j))
        }
      } else {
        idx <- c(idx,c(j))
      }
    }
    X_obs <- X[idx,]
    y_obs <- y[idx]
    res <- SGD(X_obs,y_obs,10,res[["theta"]],res[["alpha"]])
  }
  printf("final values in %dth iteration\n",kk)
  print(res[["theta"]])
  print(res[["alpha"]])
  fres = rbind(fres, c(c(res["alpha"]), c(res[["theta"]])))
}
print(fres)
fres <- matrix(unlist(fres), nrow=simu, ncol=p+1)
print("Mean")
print(apply(fres,2,mean))
print("Variance")
print(apply(fres,2,sd))
fmean <- apply(fres,2,mean)
original <- c(c(alpha), theta)
print("Bias")
print((fmean-original)/original)
