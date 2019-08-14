library(comprehenr)
library(matlib)
library(MASS)
library(glmnet)

OLSE <- function(X, y) {
  r <- seq(0,1,by=0.02)
  min_CVE <- 10000000
  lambda_min <- 1
  alpha_min <- 1
  for (alpha in r) {
    res <- cv.glmnet(X, y, type.measure="deviance", alpha=alpha, intercept=TRUE)
    if(min(res[["cvm"]])<=min_CVE){
      min_CVE = min(res[["cvm"]])
      lambda_min =  res[["lambda"]][which.min(res[["cvm"]])]
      alpha_min <- alpha
    }
    #print(list(lambda_min, alpha_min, min_CVE, min(res[["cvm"]])))
  }
  res <- glmnet(X, y, family = "gaussian", lambda = lambda_min, alpha=alpha_min,intercept = TRUE)
  res1 <- res[["beta"]]
  res1[1] = res[["a0"]]
  return(res1)
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
  sum = sum(to_vec(for(i in 1:n)((4/alphah^2)*(sinh(0.5*(y[i]-z[i]))*sinh(0.5*(y[i]-z[i]))))))/n
  return(sum+1)
}

SGD <- function(X, y) {
  theta_hat <- OLSE(X,y)
  alpha_hat <- alpha_mle(y, X%*%theta_hat, length(y))
  print(alpha_hat)
  for(i in 1:10){
    y_t <- X%*%theta_hat + (2/C(alpha_hat))*R(y,X%*%theta_hat, length(y),alpha_hat)
    theta_hat <- OLSE(X,y_t)
    alpha_hat <- (alpha_hat*sum_z_sq(y,X%*%theta_hat, length(y),alpha_hat))/2
    print(alpha_hat)
  }
  return(list("theta"=theta_hat,"alpha"=alpha_hat))
}

n = 100
alpha = 0.5
p = 204
crvalue = 1
theta = c(c(1.5,5,3,2,-0.5), seq(0,0,length.out = 200))
set.seed(124)

X = matrix(nrow = 0, ncol = 205)
for(i in 1:n){
  X <- rbind(X,c(c(1),rnorm(p,0,2)))
}

initial_estimate <-function(X,y,delta) {
  idx_obs <- to_vec(for(i in 1:n) if(delta[i]==0) i)
  X_obs <- X[idx_obs,]
  y_obs <- y[idx_obs]
  return(SGD(X_obs, y_obs))
}

t <- exp(X%*%theta)
v <- rnorm(n,0,(0.25/(alpha^2)))
y <- log(to_vec(for(i in 1:n) (t[i]*(1+2*(v[i]^2)+2*v[i]*sqrt(1+v[i]^2)))))
c <- exp(crvalue + sqrt(2)*rnorm(n, 0, 1))
delta <- seq(0,0,length.out = n)

for (i in 1:n) {
  if(y[i]>=c[i]){
    y[i]=c[i]
    delta[i]=1
  }
}

res <- initial_estimate(X,y,delta)

sample_z <- function(x,y,theta,alpha,c){
  while(1){
    v <- rnorm(1,0,(0.25/(alpha^2)))
    y <- (x%*%theta)[1][1] + log(1 + 2*v^2 + 2*v*sqrt(1+v^2))
    if(y>=c){
      print(y,c)
      return(y)
    }
  }
}
print(n)
for (i in 1:1) {
  for(j in 1:n){
    if(delta[j]==1){
      y[j]<-sample_z(X[j,],y[j],res[["theta"]],res[["alpha"]],c[j])
    }
  }
  
  print(y)
  res <- SGD(X,y)
}
