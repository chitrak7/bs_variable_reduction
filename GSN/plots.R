library(ggplot2)
p=0.1

x <- seq(-20,20,by=0.01)
y1 <- dnorm(x,0,1/sqrt(p)) 

rng <- 1:100
rng_rt <- sqrt(rng)
p_arr = (1-p)^(rng-1)
y2 <- rep(0,length(x))
for (i in 1:length(x)) {
  y2[i] = p*sum(p_arr*dnorm(x[i],0,rng_rt))
}
plot(x,y2,type="l",col="red")
lines(x,y1,col="green")