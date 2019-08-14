BS_CDF <- function(x,alpha,beta) {
  
  return(pnorm((sqrt(x/beta) - sqrt(beta/x))/alpha))
}

BS_PDF <- function(x, alpha, beta) {
  z <- sqrt(x/beta)
  return(dnorm((z^2 -1)/(z*alpha),0,alpha^2)*((z^2)+1)/(2*z*x))
}

z <- BS_PDF(1,1,1)*sqrt(2*3.14)


x <- seq(-10,10,by=0.1)
y <- pnorm(x, 0, 1)
plot(x,y)

