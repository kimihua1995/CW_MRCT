# source for 2 covariates and 3 regions
# non-PH for treatment arms

# Generate Regional Data
gen.dat.2x3r <- function(n,r,a,lambda,gamma){
  a <- matrix(a,ncol = 1)
  X <- matrix(0,nrow=1,ncol = 2)
  A <- R <- NULL
  nR <- length(n)
  for (j in 1:nR){
    while (nrow(X) <= sum(n[1:j])){
        X1 <- runif(1)
        X2 <- rnorm(1,1,1)
      pR <- 1/(1+exp(-c(1,X1^0.5,X2^2) %*% r[,j]))
      R.temp <- rbinom(1,1,pR)
      if (R.temp==1) {X <- rbind(X,c(X1,X2))}
    }
    A <- c(A, rbinom(n[j],1,0.5))
    R <- c(R, rep(j,n[j]))
  }
  X <- X[-1,]
  colnames(X) <- c("X1","X2")
  dR <- matrix(0,nrow = sum(n),ncol = nR) # dummy variable for R
  for (j in 1:nR){
    dR[R==j, j] <- 1
  }
  colnames(dR) <- paste0("R",1:nR)
  
  
  df <- data.frame(A=A,X,R=R,dR)
  MM <- model.matrix(~ -1+A+R2+R3+X1+X2+
                       A*R2+A*R3+A*X1+A*X2+R2*X1+R2*X2+R3*X1+R3*X2,df)
  MM <- MM[,-1]
  # log normal
  #epsilon <- rnorm(n,mean=0,sd=0.1)
  #time_e <- c(exp(epsilon - MM %*% a))
  # Exp
  u <- runif(sum(n))
  time_e <- (-log(u)*exp(-MM %*% a)/lambda[A+1])^(1/gamma[A+1])
  time_c <- runif(sum(n),1,10)
  time <- pmin(time_e,time_c)
  status <- as.numeric(time_e < time_c)
  df$Y <- time
  df$status <- status
  
  
  return(as.data.frame(df %>% group_by(R) %>% arrange(Y, .by_group = T)))
}





MRCT <- function(seed,S=1000,n,a,r,tau,lambda,gamma){
  mu.list <- mu.sd.list <- p.list <- NULL
  nR <- ncol(r)
  set.seed(seed)
  for (s in 1:S){
    df <- gen.dat.2x3r(n,r,a,lambda,gamma)
    
    # Naive Method
    fit.naive <- Naive.Est(nR,df,tau)
    res.naive <- region.diff(fit.naive,nR)
    N.mu <- res.naive$mu
    N.sd <- res.naive$sd
    N.p <- region.cons.test(fit.naive,nR)
    
    # IPW
    fit.ipw <- IPW.Est(nR,2,df,tau)
    res.ipw <- region.diff(fit.ipw,nR)
    IPW.mu <- res.ipw$mu
    IPW.sd <- res.ipw$sd
    IPW.p <- region.cons.test(fit.ipw,nR)
    
    # Entropy Weighting Lagrange Multiplier
    fit.ew <- CW.Est(nR,2,c("C","C"),df,tau)
    res.ew <- region.diff(fit.ew,nR)
    EW.mu <- res.ew$mu
    EW.sd <- res.ew$sd
    EW.p <- region.cons.test(fit.ew,nR)
    
    mu.list <- rbind(mu.list, c(N.mu, IPW.mu, EW.mu))
    mu.sd.list <- rbind(mu.sd.list, c(N.sd, IPW.sd, EW.sd))
    p.list <- rbind(p.list, c(N.p, IPW.p, EW.p))
  
  }
  
  region.diff.pair <- NULL
  for (j in 1:(nR-1)) {
    region.diff.pair <- c(region.diff.pair, paste0(((j+1):nR),j))
  }
  colnames(mu.list) <- colnames(mu.sd.list) <- 
    c(paste0(rep(c("N.mu1","N.mu0","N.delta"),nR),"_",rep(1:nR,each=3)),
      paste0("N.delta_",region.diff.pair),
      paste0(rep(c("IPW.mu1","IPW.mu0","IPW.delta"),nR),"_",rep(1:nR,each=3)),
      paste0("IPW.delta_",region.diff.pair),
      paste0(rep(c("EW.mu1","EW.mu0","EW.delta"),nR),"_",rep(1:nR,each=3)),
      paste0("EW.delta_",region.diff.pair))
  
  colnames(p.list) <- c("Naive","IPW","EW")
  return(list(mu = as.data.frame(mu.list),
              sd = as.data.frame(mu.sd.list),
              p = as.data.frame(p.list)))
}



true.rmst <- function(tau,n,a,r,lambda,gamma){
  nR <- ncol(r)
  fxrj <- function(x1,x2,j){
    num <- function(x,y){dnorm(y,1,1)/(1+exp(-r[1,j]-x^0.5*r[2,j]-y^2*r[3,j]))}
    dnm <- pracma::integral2(fun = num, xmin = 0, xmax = 1, ymin = -5, ymax = 7)$Q
    return(num(x1,x2)/dnm)
  }
  
  fxr <- function(x1,x2){
    mix_dist <- 0
    for (j in 1:nR){
      mix_dist <- mix_dist + n[j]/sum(n)*fxrj(x1,x2,j)
    }
    return(mix_dist)
  }
  
  # region, trt
  fun_11 <- function(x,y,z){
    fxr(x,y)*exp(-lambda[2]*(z^gamma[2])*exp((a[3]+a[7])*x+(a[4]+a[8])*y))
  }
  fun_10 <- function(x,y,z){
    fxr(x,y)*exp(-lambda[1]*(z^gamma[1])*exp(a[3]*x+a[4]*y))
  }
  fun_21 <- function(x,y,z){
    fxr(x,y)*exp(-lambda[2]*(z^gamma[2])*exp(a[1]+a[5]+(a[3]+a[7]+a[9])*x+(a[4]+a[8]+a[10])*y))
  }
  fun_20 <- function(x,y,z){
    fxr(x,y)*exp(-lambda[1]*(z^gamma[1])*exp(a[1]+(a[3]+a[9])*x+(a[4]+a[10])*y))
  }
  fun_31 <- function(x,y,z){
    fxr(x,y)*exp(-lambda[2]*(z^gamma[2])*exp(a[2]+a[6]+(a[3]+a[7]+a[11])*x+(a[4]+a[8]+a[12])*y))
  }
  fun_30 <- function(x,y,z){
    fxr(x,y)*exp(-lambda[1]*(z^gamma[1])*exp(a[2]+(a[3]+a[11])*x+(a[4]+a[12])*y))
  }

  xmin <- 0; xmax <- 1
  ymin <- -5; ymax <- 7

  mu11 <- pracma::integral3(fun_11, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, 
                            zmin=0, zmax=tau, reltol = 1e-5)
  mu10 <- pracma::integral3(fun_10, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, 
                            zmin=0, zmax=tau, reltol = 1e-5)
  mu21 <- pracma::integral3(fun_21, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, 
                            zmin=0, zmax=tau, reltol = 1e-5)
  mu20 <- pracma::integral3(fun_20, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, 
                            zmin=0, zmax=tau, reltol = 1e-5)
  mu31 <- pracma::integral3(fun_31, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, 
                            zmin=0, zmax=tau, reltol = 1e-5)
  mu30 <- pracma::integral3(fun_30, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, 
                            zmin=0, zmax=tau, reltol = 1e-5)
  
  true.mu <- c(mu11,mu10,mu11-mu10,
               mu21,mu20,mu21-mu20,
               mu31,mu30,mu31-mu30)
  
  
  for (j in 1:(nR-1)){
    true.mu <- c(true.mu, true.mu[3*((j+1):nR)] - true.mu[3*j])
  }
  
  true.mu <- as.data.frame(matrix(true.mu,nrow=1))
  region.diff.pair <- NULL
  for (j in 1:(nR-1)) {
    region.diff.pair <- c(region.diff.pair, paste0(((j+1):nR),j))
  }
  colnames(true.mu) <-  c(paste0(rep(c("mu1","mu0","delta"),nR),"_",rep(1:nR,each=3)),
                          paste0("delta_",region.diff.pair))
  return(true.mu)
}



res.output <- function(res,n,a,r,tau,lambda,gamma){
  res.mu <- res$mu; res.sd <- res$sd
  
  res.mu.m <- colMeans(res.mu,na.rm = T); res.sd.m <- colMeans(res.sd,na.rm = T)
  res.mu.esd <- apply(res.mu,2,sd); res.sd.esd <- apply(res.sd,2,sd)
  
  res.l <- res.mu - 1.96 * res.sd
  res.h <- res.mu + 1.96 * res.sd
  
  cover.p <- function(res.l,res.h,true.delta){
    cover <- (res.l <= c(true.delta))*(res.h >= c(true.delta))
    cp <- sum(cover == 1,na.rm = T)/length(cover)
    round(cp*100,1)
  }
  
  
  rmst <- true.rmst(tau,n,a,r,lambda,gamma)
  
  output <- list()
  output[[1]] <- output[[2]] <- output[[3]] <- matrix(NA,nrow=4,ncol = 15)
  for (i in 1:(ncol(res.mu)/3/3)){
    for (j in 0:2){
      output[[j+1]][i,] <- round(unlist(c(
        res.mu.m[ncol(res.mu)*j/3 + (3*i-2):(3*i)],
        res.mu.m[ncol(res.mu)*j/3 + (3*i-2):(3*i)] - rmst[(3*i-2):(3*i)],
        res.sd.m[ncol(res.mu)*j/3 + (3*i-2):(3*i)],
        res.mu.esd[ncol(res.mu)*j/3 + (3*i-2):(3*i)],
        cover.p(res.l[,ncol(res.mu)*j/3 + 3*i-2],res.h[,ncol(res.mu)*j/3 + 3*i-2],rmst[3*i-2]),
        cover.p(res.l[,ncol(res.mu)*j/3 + 3*i-1],res.h[,ncol(res.mu)*j/3 + 3*i-1],rmst[3*i-1]),
        cover.p(res.l[,ncol(res.mu)*j/3 + 3*i],res.h[,ncol(res.mu)*j/3 + 3*i],rmst[3*i])
      )),3)
    }
  }
  
  tab <- NULL
  for (i in 1:(ncol(res.mu)/3/3)){
    for (j in 1:3){
      tab <- rbind(tab, output[[j]][i,])
    }
  }
  
  rownames(tab) <- rep(c("Naive","IPW","CW"),4)
  colnames(tab) <- c(paste0("RMST",1:3),
                     paste0("Bias",1:3),
                     paste0("SD",1:3),
                     paste0("SE",1:3),
                     paste0("CP",1:3))
  return(as.data.frame(tab))
}








# simulation example'
source("source_MRCT.R")
a3 = c(0.3,0.5,0.5,0.5,0.3,0.5,-1.5,-1,0.3,0.3,0.5,0.5)
r3 = matrix(c(0.5,1,-1,0.5,1,-2,0.5,-1,1),nrow=3,ncol=3)
lambda = c(0.5,0.3)
gamma = c(1,0.3)
n = c(200,200,200)
res <- MRCT(seed=12321,S=100,n=n,a=a3,r=r3,tau=5,lambda=lambda,gamma=gamma)
output <- res.output(res,n=n,a=a3,r=r3,tau=5,lambda=lambda,gamma=gamma)


