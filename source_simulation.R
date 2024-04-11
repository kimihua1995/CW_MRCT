# Generate Regional Data
## Log-linear sampling
gen.dat1 <- function(n,r,a,lambda,gamma){
  a <- matrix(a,ncol = 1)
  X <- matrix(NA,nrow=1,ncol = 2)
  A <- R <- rho <- NULL
  nR <- length(n)
  for (j in 1:nR){
    while (nrow(X) <= sum(n[1:j])){
      X1 <- runif(1,0,1)
      X2 <- rnorm(1,1,1)
      pR <- exp(c(1,X1,X2) %*% r[,j])
      R.temp <- rbinom(1,1,min(pR,1))
      if (R.temp==1) {
        X <- rbind(X,c(X1,X2))
        rho <- c(rho, pR)
      }
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
  
  
  df <- data.frame(A=A,X,R=R,dR,rho)
  MM <- model.matrix(~ -1+A+R2+R3+X1+X2+
                       A*R2+A*R3+A*X1+A*X2+
                       R2*X1+R2*X2+R3*X1+R3*X2,df)
  MM <- MM[,-1]
  u <- runif(sum(n))
  time_e <- (-log(u)*exp(- MM %*% a)/lambda[A+1])^(1/gamma[A+1])
  time_c <- rexp(sum(n))/0.1
  time <- pmin(time_e,time_c)
  status <- as.numeric(time_e < time_c)
  df$Y <- time
  df$status <- status
  
  
  return(as.data.frame(df %>% group_by(R) %>% arrange(Y, .by_group = T)))
}


## Logistic-nonlinear sampling
gen.dat2 <- function(n,r,a,lambda,gamma){
  a <- matrix(a,ncol = 1)
  X <- matrix(NA,nrow=1,ncol = 2)
  A <- R <- rho <- NULL
  nR <- length(n)
  for (j in 1:nR){
    while (nrow(X) <= sum(n[1:j])){
      X1 <- runif(1,0,1)
      X2 <- rnorm(1,1,1)
      pR <- 1/(1+exp(-c(1,X1*X2,exp(X2/10)) %*% r[,j]))
      R.temp <- rbinom(1,1,pR)
      if (R.temp==1) {
        X <- rbind(X,c(X1,X2))
        rho <- c(rho, pR)
        }
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
  
  
  df <- data.frame(A=A,X,R=R,dR,rho)
  MM <- model.matrix(~ -1+A+R2+R3+X1+X2+
                       A*R2+A*R3+A*X1+A*X2+
                       R2*X1+R2*X2+R3*X1+R3*X2,df)
  MM <- MM[,-1]
  u <- runif(sum(n))
  time_e <- (-log(u)*exp(- MM %*% a)/lambda[A+1])^(1/gamma[A+1])
  time_c <- rexp(sum(n))/0.1
  time <- pmin(time_e,time_c)
  status <- as.numeric(time_e < time_c)
  df$Y <- time
  df$status <- status
  
  
  return(as.data.frame(df %>% group_by(R) %>% arrange(Y, .by_group = T)))
}





## true values of the RMST and RMST difference
true.rmst.fix <- function(tau,nR,a,lambda,gamma){
  fxr <- function(x1,x2){
    dnorm(x2,1,1)
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
  
  
  mu11 <- pracma::integral3(fun_11, xmin=0, xmax=1, ymin=-4, ymax=6, 
                            zmin=0, zmax=tau, reltol = 1e-6)
  mu10 <- pracma::integral3(fun_10, xmin=0, xmax=1, ymin=-4, ymax=6, 
                            zmin=0, zmax=tau, reltol = 1e-6)
  mu21 <- pracma::integral3(fun_21, xmin=0, xmax=1, ymin=-4, ymax=6, 
                            zmin=0, zmax=tau, reltol = 1e-6)
  mu20 <- pracma::integral3(fun_20, xmin=0, xmax=1, ymin=-4, ymax=6, 
                            zmin=0, zmax=tau, reltol = 1e-6)
  mu31 <- pracma::integral3(fun_31, xmin=0, xmax=1, ymin=-4, ymax=6, 
                            zmin=0, zmax=tau, reltol = 1e-6)
  mu30 <- pracma::integral3(fun_30, xmin=0, xmax=1, ymin=-4, ymax=6, 
                            zmin=0, zmax=tau, reltol = 1e-6)
  
  true.mu <- c(mu11,mu10,mu11-mu10,
               mu21,mu20,mu21-mu20,
               mu31,mu30,mu31-mu30)
  true.mu <- as.data.frame(matrix(true.mu,nrow=1))
  colnames(true.mu) <-  paste0(rep(c("mu1","mu0","delta"),nR),"_",rep(1:nR,each=3))
                        
  
  return(true.mu)
}





region.diff <- function(res,nR){
  res.mu <- res$mu
  res.sd <- res$sd
  res.V <- res.sd^2
  for (r in 1:(nR-1)){
    res.mu <- c(res.mu, res.mu[3*((r+1):nR)] - res.mu[3*r])
    res.sd <- c(res.sd, sqrt(res.V[3*((r+1):nR)] + res.V[3*r]))
  }
  return(list(mu = res.mu, sd = res.sd))
}


## test for regional consistency
region.cons.test <- function(fit, nR){
  rmst_diff <- matrix(fit$mu[3*(1:nR)], ncol = 1)
  rmst_diff_sd <- fit$sd[3*(1:nR)]
  
  V <- diag(rmst_diff_sd^2)
  C <- matrix(0, nrow = nR-1, ncol = nR)
  C[,1] <- -1
  for (i in 1:(nR-1)) {C[i,i+1] <- 1}
  
  Q <- t(C %*% rmst_diff) %*% solve(C %*% V %*% t(C)) %*% (C %*% rmst_diff)
  
  #p <- 1 - pchisq(Q, nR-1)
  
  return(as.vector(Q))
}



## code for output all proposed estimators given weight
gen.est <- function(df,nR,p,tau){
  fit.km <- W.KM.Est(df,nR,p,tau)
  KM.mu <- fit.km$mu
  KM.sd <- fit.km$sd
  
  fit.gf <- W.GF.Est(df,nR,p,tau)
  GF.mu <- fit.gf$mu
  GF.sd <- fit.gf$sd
  
  fit.gf.mis <- W.GF.Est.mis(df,nR,p,tau)
  GFmis.mu <- fit.gf.mis$mu
  GFmis.sd <- fit.gf.mis$sd
  
  skip_to_next <- FALSE
  tryCatch(fit.hj <- W.HJ.Est(df,nR,p,tau),
           error = function(e){skip_to_next <<- TRUE})
  if (skip_to_next) {
    HJ.mu <- rep(NA,3*nR)
    HJ.sd <- rep(NA,3*nR)
  }else{
    HJ.mu <- fit.hj$mu
    HJ.sd <- fit.hj$sd
  }
  
  skip_to_next <- FALSE
  tryCatch(fit.ag <- W.AG.Est(df,nR,p,tau),
           error = function(e){skip_to_next <<- TRUE})
  if (skip_to_next) {
    AG.mu <- rep(NA,3*nR)
    AG.sd <- rep(NA,3*nR)
  }else{
    AG.mu <- fit.ag$mu
    AG.sd <- fit.ag$sd
  }
  
  skip_to_next <- FALSE
  tryCatch(fit.ag.mis <- W.AG.Est.mis(df,nR,p,tau),
           error = function(e){skip_to_next <<- TRUE})
  if (skip_to_next) {
    AGmis.mu <- rep(NA,3*nR)
    AGmis.sd <- rep(NA,3*nR)
  }else{
    AGmis.mu <- fit.ag.mis$mu
    AGmis.sd <- fit.ag.mis$sd
  }
  
  
  mu.list <- c(KM.mu, GF.mu, GFmis.mu, HJ.mu, AG.mu, AGmis.mu)
  sd.list <- c(KM.sd, GF.sd, GFmis.sd, HJ.sd, AG.sd, AGmis.sd)
  
  return(list(mu = mu.list, sd = sd.list))
  
}




## code for simulation 
sim.MRCT.Est <- function(seed,S=1000,n,a,r,tau,lambda,gamma,g,gen){
  nR <- ncol(r)
  mu.Naive.list <- sd.Naive.list <- matrix(NA, nrow = S, ncol = nR*3)
  mu.IPW1.list <- sd.IPW1.list <- matrix(NA, nrow = S, ncol = nR*6*3)  # true PS
  mu.IPW2.list <- sd.IPW2.list <- matrix(NA, nrow = S, ncol = nR*6*3)  # estimated PS
  mu.EW.list <- sd.EW.list <- matrix(NA, nrow = S, ncol = nR*6*3)  # calibration weight
  
  
  colnames(mu.Naive.list) <- colnames(sd.Naive.list) <- 
    paste0(rep(c("N.mu1","N.mu0","N.delta"),nR),"_",rep(1:nR,each=3))
  colnames(mu.IPW1.list) <- colnames(sd.IPW1.list) <-
    c(paste0(rep(c("IPW1.KM.mu1","IPW1.KM.mu0","IPW1.KM.delta"),nR),"_",rep(1:nR,each=3)),
      paste0(rep(c("IPW1.GF.mu1","IPW1.GF.mu0","IPW1.GF.delta"),nR),"_",rep(1:nR,each=3)),
      paste0(rep(c("IPW1.GFmis.mu1","IPW1.GFmis.mu0","IPW1.GFmis.delta"),nR),"_",rep(1:nR,each=3)),
      paste0(rep(c("IPW1.HJ.mu1","IPW1.HJ.mu0","IPW1.HJ.delta"),nR),"_",rep(1:nR,each=3)),
      paste0(rep(c("IPW1.AG.mu1","IPW1.AG.mu0","IPW1.AG.delta"),nR),"_",rep(1:nR,each=3)),
      paste0(rep(c("IPW1.AGmis.mu1","IPW1.AGmis.mu0","IPW1.AGmis.delta"),nR),"_",rep(1:nR,each=3)))
  colnames(mu.IPW2.list) <- colnames(sd.IPW2.list) <-
    c(paste0(rep(c("IPW2.KM.mu1","IPW2.KM.mu0","IPW2.KM.delta"),nR),"_",rep(1:nR,each=3)),
      paste0(rep(c("IPW2.GF.mu1","IPW2.GF.mu0","IPW2.GF.delta"),nR),"_",rep(1:nR,each=3)),
      paste0(rep(c("IPW2.GFmis.mu1","IPW2.GFmis.mu0","IPW2.GFmis.delta"),nR),"_",rep(1:nR,each=3)),
      paste0(rep(c("IPW2.HJ.mu1","IPW2.HJ.mu0","IPW2.HJ.delta"),nR),"_",rep(1:nR,each=3)),
      paste0(rep(c("IPW2.AG.mu1","IPW2.AG.mu0","IPW2.AG.delta"),nR),"_",rep(1:nR,each=3)),
      paste0(rep(c("IPW2.AGmis.mu1","IPW2.AGmis.mu0","IPW2.AGmis.delta"),nR),"_",rep(1:nR,each=3)))
  colnames(mu.EW.list) <- colnames(sd.EW.list) <-
    c(paste0(rep(c("EW.KM.mu1","EW.KM.mu0","EW.KM.delta"),nR),"_",rep(1:nR,each=3)),
      paste0(rep(c("EW.GF.mu1","EW.GF.mu0","EW.GF.delta"),nR),"_",rep(1:nR,each=3)),
      paste0(rep(c("EW.GFmis.mu1","EW.GFmis.mu0","EW.GFmis.delta"),nR),"_",rep(1:nR,each=3)),
      paste0(rep(c("EW.HJ.mu1","EW.HJ.mu0","EW.HJ.delta"),nR),"_",rep(1:nR,each=3)),
      paste0(rep(c("EW.AG.mu1","EW.AG.mu0","EW.AG.delta"),nR),"_",rep(1:nR,each=3)),
      paste0(rep(c("EW.AGmis.mu1","EW.AGmis.mu0","EW.AGmis.delta"),nR),"_",rep(1:nR,each=3)))

  
  set.seed(seed)
  for (s in 1:S){
    if (!s%%100) print(s)
    if (gen == 1){
      df <- gen.dat1(n,r,a,lambda,gamma)
    }else if (gen == 2){
      df <- gen.dat2(n,r,a,lambda,gamma)
    }
    
    #1) Naive Method
    p1 <- NULL
    for (i in 1:nR){
      p1 <- c(p1, rep(1,n[i])/n[i])
    }
    fit.naive <- W.KM.Est(df,nR,p1,tau)
    mu.Naive.list[s,] <- fit.naive$mu
    sd.Naive.list[s,] <- fit.naive$sd
    
    
    ############################################
    #2) True PS
    p2 <- 1/df$rho
    
    Est.IPW1 <- gen.est(df,nR,p2,tau)
    mu.IPW1.list[s,] <- Est.IPW1$mu
    sd.IPW1.list[s,] <- Est.IPW1$sd



    #############################################
    #3) Estimated PS
    fit.mnps <- mnps(factor(R) ~ X1 + X2, data = df,
                     estimand = "ATE",
                     verbose = FALSE,
                     stop.method = "es.mean",
                     n.trees = 3000)
    p3 <- get.weights(fit.mnps, stop.method = "es.mean", estimand = "ATE")
    
    Est.IPW2 <- gen.est(df,nR,p3,tau)
    mu.IPW2.list[s,] <- Est.IPW2$mu
    sd.IPW2.list[s,] <- Est.IPW2$sd
    
    
    
    
    #4) Entropy Weighting Lagrange Multiplier
    skip_to_next <- FALSE
    tryCatch(p4 <- cw2(df,nR,g),
             error = function(e){skip_to_next <<- TRUE})
    if (skip_to_next) {p4 <- p1}
    
    Est.EW <- gen.est(df,nR,p4,tau)
    mu.EW.list[s,] <- Est.EW$mu
    sd.EW.list[s,] <- Est.EW$sd
  
    
  }
  

  return(list(mu.N = as.data.frame(mu.Naive.list),
              sd.N = as.data.frame(sd.Naive.list),
              mu.IPW1 = as.data.frame(mu.IPW1.list),
              sd.IPW1 = as.data.frame(sd.IPW1.list),
              mu.IPW2 = as.data.frame(mu.IPW2.list),
              sd.IPW2 = as.data.frame(sd.IPW2.list),
              mu.EW = as.data.frame(mu.EW.list),
              sd.EW = as.data.frame(sd.EW.list)))
}




# check the results
res.check <- function(res_mu,res_sd,rmst_true){
  rmst <- colMeans(res_mu, na.rm = T); l <- length(rmst)
  
  dat_bias <- matrix(rmst,nrow = 6,ncol = l/6,byrow = T) - 
     matrix(as.numeric(rmst_true), nrow = 6, ncol = l/6, byrow = T)
  
  res.l <- res_mu - 1.96 * res_sd
  res.h <- res_mu + 1.96 * res_sd
  
  cover.p <- function(res.l,res.h,true.delta){
    cover <- (res.l <= c(true.delta))*(res.h >= c(true.delta))
    cp <- sum(cover == 1,na.rm = T)/length(cover)
    return(cp*100)
  }
  
  CP_KM <- CP_GF <- CP_GFmis <- CP_HJ <- CP_AG <- CP_AGmis <-  NULL
  for (i in 1:length(rmst_true)){
    CP_KM <- c(CP_KM, cover.p(res.l[,i],res.h[,i],rmst_true[i]))
    CP_GF <- c(CP_GF, cover.p(res.l[,9+i],res.h[,9+i],rmst_true[i]))
    CP_GFmis <- c(CP_GFmis, cover.p(res.l[,18+i],res.h[,18+i],rmst_true[i]))
    CP_HJ <- c(CP_HJ, cover.p(res.l[,27+i],res.h[,27+i],rmst_true[i]))
    CP_AG <- c(CP_AG, cover.p(res.l[,36+i],res.h[,36+i],rmst_true[i]))
    CP_AGmis <- c(CP_AGmis, cover.p(res.l[,45+i],res.h[,45+i],rmst_true[i]))
  }
  dat_CP <- rbind(CP_KM, CP_GF, CP_GFmis, CP_HJ, CP_AG, CP_AGmis)
  dat_bias <- as.data.frame(dat_bias)
  dat_CP <- as.data.frame(dat_CP)
  rownames(dat_CP) <- rownames(dat_bias) <- 
    c("KM","GF","GFmis","HJ", "AG","AGmis")
  nR=3
  colnames(dat_CP) <- colnames(dat_bias) <- 
    paste0(rep(c("mu1","mu0","delta"),nR),"_",rep(1:nR,each=3))
  
  return(list(bias = dat_bias, CP = dat_CP))
}




# absolute standardized mean difference
SMD <- function(X,R,mu,sd){
  nR <- length(table(R))
  SMD_x <- rep(NA, nR)
  
  for (i in 1:nR){
      Xi <- X[R==i]
      mui <- mean(Xi)
      sdi <- sd(Xi)
      SMD_x[i] <- abs(mui - mu)/sqrt((sdi^2 + sd^2)/2)
  }
  
  return(SMD_x)
}




plot.data.gather <- function(mu, names, labels){
  names_delta <- paste0(names,".delta_")
  names_mu1 <- paste0(names,".mu1_")
  names_mu0 <- paste0(names,".mu0_")
  l <- length(labels)
  
  delta1 <- mu[,paste0(names_delta,1)] %>% tidyr::gather(Method,RMST,1:l) %>%
    mutate(Method = factor(Method, levels = paste0(names_delta,1),
                           labels = labels))
  delta2 <- mu[,paste0(names_delta,2)] %>% tidyr::gather(Method,RMST,1:l) %>%
    mutate(Method = factor(Method, levels = paste0(names_delta,2),
                           labels = labels))
  delta3 <- mu[,paste0(names_delta,3)] %>% tidyr::gather(Method,RMST,1:l) %>%
    mutate(Method = factor(Method, levels = paste0(names_delta,3),
                           labels = labels))
  mu1_1 <- mu[,paste0(names_mu1,1)] %>% tidyr::gather(Method,RMST,1:l) %>%
    mutate(Method = factor(Method, levels = paste0(names_mu1,1),
                           labels = labels))
  mu1_2 <- mu[,paste0(names_mu1,2)] %>% tidyr::gather(Method,RMST,1:l) %>%
    mutate(Method = factor(Method, levels = paste0(names_mu1,2),
                           labels = labels))
  mu1_3 <- mu[,paste0(names_mu1,3)] %>% tidyr::gather(Method,RMST,1:l) %>%
    mutate(Method = factor(Method, levels = paste0(names_mu1,3),
                           labels = labels))
  mu0_1 <- mu[,paste0(names_mu0,1)] %>% tidyr::gather(Method,RMST,1:l) %>%
    mutate(Method = factor(Method, levels = paste0(names_mu0,1),
                           labels = labels))
  mu0_2 <- mu[,paste0(names_mu0,2)] %>% tidyr::gather(Method,RMST,1:l) %>%
    mutate(Method = factor(Method, levels = paste0(names_mu0,2),
                           labels = labels))
  mu0_3 <- mu[,paste0(names_mu0,3)] %>% tidyr::gather(Method,RMST,1:l) %>%
    mutate(Method = factor(Method, levels = paste0(names_mu0,3),
                           labels = labels))
  
  
  return(list(delta1 = delta1, delta2 = delta2, delta3 = delta3,
              mu1_1 = mu1_1, mu1_2 = mu1_2, mu1_3 = mu1_3,
              mu0_1 = mu0_1, mu0_2 = mu0_2, mu0_3 = mu0_3))
}




plot.data <- function(res){
  mu.IPW1 <- cbind(res$mu.N, res$mu.IPW1)
  mu.IPW2 <- cbind(res$mu.N, res$mu.IPW2)
  mu.EW <- cbind(res$mu.N, res$mu.EW)
  
  
  names.IPW1 <- c("N","IPW1.KM","IPW1.GF","IPW1.GFmis","IPW1.HJ","IPW1.AG","IPW1.AGmis")
  names.IPW2 <- c("N","IPW2.KM","IPW2.GF","IPW2.GFmis","IPW2.HJ","IPW2.AG","IPW2.AGmis")
  names.EW <-  c("N","EW.KM","EW.GF","EW.GFmis","EW.HJ","EW.AG","EW.AGmis")
  labels <- c("Naive","KM","GF","GFmis","HJ","AG","AGmis")
              

  dat1 <- plot.data.gather(mu.IPW1, names.IPW1, labels)
  dat2 <- plot.data.gather(mu.IPW2, names.IPW2, labels)
  dat3 <- plot.data.gather(mu.EW, names.EW, labels)
  
  
  return(list(IPW1 = dat1, IPW2 = dat2, CW = dat3))
}







