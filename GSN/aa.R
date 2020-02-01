library(moments)
library(ggplot2)


# results <- data.frame(
#   mu <- c(),
#   sigma <- c(),
#   p <- c(),
#   n <- c(),
#   alpha <- c(),
#   mu_hat <- c(),
#   p_hat <- c(),
#   sigma_hat <- c(),
#   l1ymean <- c(),
#   l1ymax<- c(),
#   l1ysd <- c(),
#   l1mean <- c(),
#   l1sd <- c(),
#   l1max <- c(),
#   msemean <- c(),
#   msemax <- c(),
#   msesd <- c(),
#   mseymean <- c(),
#   mseymax<- c(),
#   mseysd <- c()
# )

results <- readRDS("results.rds")


    p<-200
    theta <- c(c(1.5,5,3,2,-0.5), rep(0,p-5))
    args<-commandArgs(trailingOnly = TRUE)
    args <- as.numeric(args)
    args<-c(2,0.4,0.6,0.8,100)
    mu=args[1]
    sigma=args[2]
    g_p=args[3]
    alpha_s=args[4]
    n<-args[5]
    simu<-1000
    if(jj==0.0) {
      simu<-500
    }
    if(ii==0.2 && jj==0.2) {
      simu <- 500
    }
    iters<-200
    file_str<-sprintf("symGSN_mu_em_glmnet_mu=%f_p=%f_sigma=%f_alpha=%f_n=%f.rds",mu,g_p,sigma,alpha_s,n)
    title_str<-sprintf("density plot for mu=%f p=%f sigma=%f alpha=%f n=%f",mu,g_p,sigma,alpha_s,n)
    f=readRDS(file_str)
    burnin<-40
    msey <- rep(0,simu)
    mse <- rep(0,simu)
    l1 <- rep(0,simu)
    l1y <- rep(0,simu)
    ttt <- sample(1:simu, 10)
    par(mfrow=c(2,5))
    
    p_hat <- mean(f[,1])
    sigma_hat <- mean(f[,2])
    mu_hat <- mean(f[,3])
    jj=0
    for(kk in ttt) {
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
      yhat <- X%*%f[kk,1:p+3]
      img_str <- sprintf("symGSN_mu_em_glmnet_residuals_mu=%f_p=%f_sigma=%f_alpha=%f_n=%f_%d.jpg",mu,g_p,sigma,alpha_s,n,jj)
      jpeg(img_str)
      aa <- data.frame(x=y-X%*%theta)
      bb <- data.frame(x=y-yhat)
      gg<- ggplot() + geom_density(aes(x=x,color="red"),aa) + geom_density(aes(x=x,color="green"),bb) +
        scale_color_discrete(name="plot for residualts", labels=c("actual residuals", "estimated residuals"))
      plot(gg)
      dev.off()
      img_str <- sprintf("symGSN_mu_em_glmnet_estimates_mu=%f_p=%f_sigma=%f_alpha=%f_n=%f_%d.jpg",mu,g_p,sigma,alpha_s,n,jj)
      jpeg(img_str)
      bb <- data.frame(x=y, y=yhat + mu_hat/p_hat)
      gg<- ggplot() + geom_point(aes(x=x,y=y),aa)
      plot(gg)
      dev.off()
      msey[kk] <- sum((y-yhat)**2)/n
      mse[kk] <- sum((theta-f[kk,1:p+3])**2)/p
      l1[kk] <- sum(abs(theta-f[kk,1:p+3]))/p
      l1y[kk] <- sum(abs(y-yhat))/n
      jj = jj + 1
    }
    #
    gsnpts <- function(x_arr,mu,sd,p) {
      y <- rep(0,length(x_arr))
      for(i in 1:length(x_arr)) {
        x <- x_arr[i]
        mx <- 100
        rng <- 1:mx
        rng_rt <- sqrt(rng)
        p_arr <- ((1-p)^(rng-1))*p
        val  <- dnorm(x,rep(mu,mx)*rng,rng_rt*rep(sd,mx))
        y[i] <- sum(p_arr*val)
      }
      return(y)
    }
    result <- data.frame(
      mu <- c(mu),
      sigma <- c(sigma),
      p <- c(g_p),
      n <- c(n),
      alpha <- c(alpha_s),
      mu_hat <- c(mu_hat),
      p_hat <- c(p_hat),
      sigma_hat <- c(sigma_hat),
      l1ymean <- c(mean(l1y)),
      l1ymax<- c(max(l1y)),
      l1ysd <- c(sd(l1y)),
      l1mean <- c(mean(l1)),
      l1sd <- c(sd(l1)),
      l1max <- c(max(l1)),
      msemean <- c(mean(mse)),
      msemax <- c(max(mse)),
      msesd <- c(sd(mse)),
      mseymean <- c(mean(msey)),
      mseymax<- c(max(msey)),
      mseysd <- c(sd(msey))
    )
    results <- rbind(results,result)
    x <- seq(-2,20,length.out = 1000)
    y <- gsnpts(x,mu,sigma,g_p)
    y1 <- gsnpts(x,mu_hat,sigma_hat,p_hat)
    
    aa <- data.frame(x=x,y=y,lab="actual pdf")
    
    bb <- data.frame(x=x,y=y1,lab="estimated distribution")
    cc <- rbind(aa,bb)
    
    img_str<-sprintf("symGSN_mu_em_glmnet_mu=%f_p=%f_sigma=%f_alpha=%f_n=%f.jpg",mu,g_p,sigma,alpha_s,n)
    jpeg(img_str)
    gg <- ggplot() + 
      geom_line(data=aa, aes(x=x,y=y,color="red")) + 
      geom_line(data=bb, aes(x=x,y=y,color="green")) + 
      scale_color_discrete(name = "PDF for residuals", labels = c("actual distribution", "estimated distribution")) +
      labs(title = title_str)
    plot(gg)
    dev.off()
    saveRDS(results,"results.rds")