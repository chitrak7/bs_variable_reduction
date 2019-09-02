library(comprehenr)
library(matlib)
library(MASS)
library(glmnet)

p = 215
OLSE <- function(X, y) {
  X <-X[,2:215] 
  r <- seq(0,1,by=0.02)
  min_CVE <- 10000000
  lambda_min <- 1
  alpha_min <- 1
  
  res <- cv.glmnet(X, y, family="gaussian",penalty.factor=c(0,rep(1,p)))
  lambda_min = res$lambda.min
  res <- glmnet(X, y, family = "gaussian", lambda = lambda_min, alpha=alpha_min,penalty.factor=c(0,rep(1,p)))
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
  if(alpha_hat==0) {
    theta_hat <- OLSE(X,y)
    alpha_hat <- alpha_mle(y, X%*%theta_hat, length(y))
  }
  for(i in 1:k){
    y_t <- X%*%theta_hat + (2/C(alpha_hat))*R(y,X%*%theta_hat, length(y),alpha_hat)
    alpha_hat <- (alpha_hat*sum_z_sq(y,X%*%theta_hat, length(y),alpha_hat))/2
    theta_hat <- OLSE(X,y_t)
    print(alpha_hat)
  }
  return(list("theta"=theta_hat,"alpha"=alpha_hat))
}

n = 50
alpha = 0.5
crvalue = 3
theta = c(c(1,1.5,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3), seq(0,0,length.out = 199))
set.seed(1)


X = matrix(nrow = 0, ncol = 215)
for(i in 1:n){
  z <- rnorm(1,0,1)
  a <- rnorm(4,z,0.01)
  z <- rnorm(1,0,1)
  b <- rnorm(5,z,0.01)
  z <- rnorm(1,0,1)
  c <- rnorm(5,z,0.01)
  X <- rbind(X,c(c(1),rnorm(1,0,0.1),a,b,c,rnorm(199,0,2)))
}

initial_estimate <-function(X,y,delta) {
  idx_obs <- to_vec(for(i in 1:n) if(delta[i]==0) i)
  X_obs <- X[idx_obs,]
  y_obs <- y[idx_obs]
  return(SGD(X_obs, y_obs,300,0,0))
}

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
print(res)
sample_z <- function(x,y,theta_h,alpha,c){
  sset <- c()
  for(zz in 1:50){
    i=40
    while(i){
      i = i-1;
      v <- rnorm(1,0,0.5*alpha)
      z <- x%*%theta_h + log(1+2*v^2+2*v*sqrt(1+v^2))
      if(z[1,1]>c){
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

print(n)

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
  print("samples rejected")
  print(n-length(idx))
  X_obs <- X[idx,]
  y_obs <- y[idx]
  print("Resampling and recomputing")
  res <- SGD(X_obs,y_obs,50,res[["theta"]],res[["alpha"]])
  print(res)
}
print("final values")
print(res[["theta"]])
print(res[["alpha"]])
