# Generate Regional Data
## for log-linear sampling setting
gen.dat1 <- function(n,r,a,lambda,gamma){
  a <- matrix(a,ncol = 1)
  X <- matrix(NA,nrow=1,ncol = 2)
  A <- R <- NULL
  nR <- length(n)
  for (j in 1:nR){
    while (nrow(X) <= sum(n[1:j])){
      X1 <- runif(1,0,1)
      X2 <- rnorm(1,1,1)
      pR <- exp(c(1,X1,X2) %*% r[,j]); pR
      pR <- min(pR,1)
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
                       A*R2+A*R3+A*X1+A*X2+
                       R2*X1+R2*X2+R3*X1+R3*X2,df)
  MM <- MM[,-1]
  u <- runif(sum(n))
  time_e <- (-log(u)*exp(- MM %*% a)/lambda[A+1])^(1/gamma[A+1])
  time_c <- rexp(sum(n))/0.15
  time <- pmin(time_e,time_c)
  status <- as.numeric(time_e < time_c)
  df$Y <- time
  df$status <- status
  
  return(as.data.frame(df %>% group_by(R) %>% arrange(Y, .by_group = T)))
}


## for logistic sampling setting
gen.dat2 <- function(n,r,a,lambda,gamma){
  a <- matrix(a,ncol = 1)
  X <- matrix(NA,nrow=1,ncol = 2)
  A <- R <- NULL
  nR <- length(n)
  for (j in 1:nR){
    while (nrow(X) <= sum(n[1:j])){
      X1 <- runif(1,0,1)
      X2 <- rnorm(1,1,1)
      pR <- 1/(1+exp(-c(1,X1*X2,exp(X2/10)) %*% r[,j]))
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
                       A*R2+A*R3+A*X1+A*X2+
                       R2*X1+R2*X2+R3*X1+R3*X2,df)
  MM <- MM[,-1]
  
  u <- runif(sum(n))
  time_e <- (-log(u)*exp(- MM %*% a)/lambda[A+1])^(1/gamma[A+1])
  time_c <- rexp(sum(n))/0.15
  time <- pmin(time_e,time_c)
  status <- as.numeric(time_e < time_c)
  df$Y <- time
  df$status <- status
  
  
  return(as.data.frame(df %>% group_by(R) %>% arrange(Y, .by_group = T)))
}





my_akm_rmst <- function(time, status, weight=NULL, tau=NULL){
  
  data <- data.frame(time, status, weight)
  
  #--- AKM ---
  # Based on 'adjusted.KM' function from {IPWsurvival} package
  # Author: F. Le Borgne and Y. Foucher
  tj <- c(0,sort(unique(data$time[data$status==1])))
  dj <- sapply(tj, function(x){sum(data$weight[data$time==x & data$status==1])})
  yj <- sapply(tj, function(x){sum(data$weight[data$time>=x])})
  st <- cumprod(1-(dj/yj))
  m <- sapply(tj, function(x){sum((data$weight[data$time>=x])^2)})
  mj <- ((yj^2)/m)
  #ft <- data.frame(time=tj, n_risk=yj, n_event=dj, survival=st, variable=i, m=mj)
  ft <- data.frame(tj, yj, dj, st, mj)
  
  #--- RMST ---
  # Based on 'rmst1 function' from {survRM2} package
  # Author: Hajime Uno, Lu Tian, Angel Cronin, Chakib Battioui, Miki Horiguchi
  rtime <- ft$tj<=tau
  tj_r <- sort(c(ft$tj[rtime],tau))
  st_r <- ft$st[rtime]
  yj_r <- ft$yj[rtime]
  dj_r <- ft$dj[rtime]
  time_diff <- diff(c(0, tj_r))
  areas <- time_diff * c(1, st_r)
  rmst <- sum(areas)
  
  #--- Variance ---
  mj_r <- ft$mj[rtime]
  var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(mj_r *(yj_r - dj_r)))
  #var_r <- ifelse((yj_r-dj_r)==0, 0, dj_r /(yj_r *(yj_r - dj_r)))
  var_r <- c(var_r,0)
  rmst_var <- sum(cumsum(rev(areas[-1]))^2 * rev(var_r)[-1])
  
  return(data.frame(mu = rmst, V = rmst_var))     
}



W.Est <- function(df,nR,p,tau){
  EW.mu <- EW.sd <- NULL
  for (r in 1:nR){
    dfr <- df[df$R==r,]
    pr <- p[df$R == r]
    mur1 <- my_akm_rmst(dfr$Y[dfr$A==1], dfr$status[dfr$A==1],pr[dfr$A==1],tau)
    mur0 <- my_akm_rmst(dfr$Y[dfr$A==0], dfr$status[dfr$A==0],pr[dfr$A==0],tau)
    deltar <- mur1$mu - mur0$mu
    EW.mu <- c(EW.mu,mur1$mu,mur0$mu,deltar)
    deltar.sd <- sqrt(mur1$V + mur0$V)
    EW.sd <- c(EW.sd,sqrt(mur1$V),sqrt(mur0$V),deltar.sd)
  }
  return(list(mu = EW.mu, sd = EW.sd))
}




# EW Estimator
EW.p.M <- function(df,nR,f.list,iX,nX,M){
  require(Rsolnp)
  nc <- length(iX)
  varX <- paste0("X",1:nX)
  fn <- function(x){
    sum(x*log(x))
  }
  constraints <- c(M,1)
  p <- NULL
  for (r in 1:nR){
    dfr <- df[df$R == r,]
    eqn <- function(x){
      c <- NULL
      for (i in 1:nc){
        c <- c(c, sum(x*f.list[[i]](dfr[,varX[iX[i]]])))
      }
      c <- c(c,sum(x))
      return(c)
    }
    nr <- nrow(dfr)
    x0 <- rep(1/nr,nr)
    sol <- solnp(pars = x0, fun = fn, eqfun = eqn, eqB = constraints, 
                 LB = rep(0,nr), UB = rep(1,nr), control = list(trace=0, tol=1e-6))
    p <- c(p, sol$pars)
  }
  
  return(p)
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


region.cons.test <- function(fit, nR){
  rmst_diff <- matrix(fit$mu[3*(1:nR)], ncol = 1)
  rmst_diff_sd <- fit$sd[3*(1:nR)]
  
  V <- diag(rmst_diff_sd^2)
  C <- matrix(0, nrow = nR-1, ncol = nR)
  C[,1] <- -1
  for (i in 1:(nR-1)) {C[i,i+1] <- 1}
  
  Q <- t(C %*% rmst_diff) %*% solve(C %*% V %*% t(C)) %*% (C %*% rmst_diff)
  
  p <- 1 - pchisq(Q, nR-1)
  
  return(c(as.vector(Q),p))
}



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




sim.MRCT.Est <- function(seed,S=1000,n,a,r,tau,lambda,gamma,M,setting){
  nR <- ncol(r)
  mu.list <- mu.sd.list <- matrix(NA, nrow = S, ncol = (nR*3+3)*4)
  region.diff.pair <- NULL
  for (j in 1:(nR-1)) {
    region.diff.pair <- c(region.diff.pair, paste0(((j+1):nR),j))
  }
  colnames(mu.list) <- colnames(mu.sd.list) <- 
    c(paste0(rep(c("N.mu1","N.mu0","N.delta"),nR),"_",rep(1:nR,each=3)),
      paste0("N.delta_",region.diff.pair),
      paste0(rep(c("IPW.mu1","IPW.mu0","IPW.delta"),nR),"_",rep(1:nR,each=3)),
      paste0("IPW.delta_",region.diff.pair),
      paste0(rep(c("EW1.mu1","EW1.mu0","EW1.delta"),nR),"_",rep(1:nR,each=3)),
      paste0("EW1.delta_",region.diff.pair),
      paste0(rep(c("EW2.mu1","EW2.mu0","EW2.delta"),nR),"_",rep(1:nR,each=3)),
      paste0("EW2.delta_",region.diff.pair))
  
  
  set.seed(seed)
  for (s in 1:S){
    
    if (setting == 1){ # log-linear sampling setting
      df <- gen.dat1(n,r,a,lambda,gamma)
    }else if (setting == 2){ # logistic sampling setting
      df <- gen.dat2(n,r,a,lambda,gamma)
    }
    
    #1) Naive Method
    p1 <- NULL
    for (i in 1:nR){
      p1 <- c(p1, rep(1,n[i])/n[i])
    }
    fit.naive <- W.Est(df,nR,p1,tau)
    res.naive <- region.diff(fit.naive,nR)
    N.mu <- res.naive$mu
    N.sd <- res.naive$sd

    
    
    #2) IPSW Method
    fit.mnps <- mnps(factor(R) ~ X1 + X2, data = df,
                     estimand = "ATE",
                     verbose = FALSE,
                     stop.method = "es.mean",
                     n.trees = 3000)
    p2 <- get.weights(fit.mnps, stop.method = "es.mean", estimand = "ATE")
    fit.ipw <- W.Est(df,nR,p2,tau)
    res.ipw <- region.diff(fit.ipw,nR)
    IPW.mu <- res.ipw$mu
    IPW.sd <- res.ipw$sd
    
    
    
    #3) Calibration Weighting
    nX <- 2
    # 1st moment only
    f.list1 <- list()
    f.list1[[1]] <- f.list1[[2]] <- function(x) {x}
    iX1 <- c(1,2)
    
    p3 <- EW.p.M(df,nR,f.list1,iX1,nX,M[[1]])
    fit.ew1 <- W.Est(df,nR,p3,tau)
    res.ew1 <- region.diff(fit.ew1,nR)
    EW1.mu <- res.ew1$mu
    EW1.sd <- res.ew1$sd
    
    
    
    # 1st and 2nd moments
    f.list2 <- list()
    f.list2[[1]] <- f.list2[[2]] <- function(x) {x}
    f.list2[[3]] <- f.list2[[4]] <- function(x) {x^2}
    iX2 <- c(1,2,1,2)
    
    p4 <- EW.p.M(df,nR,f.list2,iX2,nX,M[[2]])
    fit.ew2 <- W.Est(df,nR,p4,tau)
    res.ew2 <- region.diff(fit.ew2,nR)
    EW2.mu <- res.ew2$mu
    EW2.sd <- res.ew2$sd

    
    mu.list[s,] <- c(N.mu, IPW.mu, EW1.mu, EW2.mu)
    mu.sd.list[s,] <- c(N.sd, IPW.sd, EW1.sd, EW2.sd)

    

  }
  
  
  return(list(mu = as.data.frame(mu.list),
              sd = as.data.frame(mu.sd.list)))
}







plot.data <- function(res){
  mu <- res$mu
  names_delta <- paste0(c("N","IPW","EW1","EW2"),".delta_")
  names_mu1 <- paste0(c("N","IPW","EW1","EW2"),".mu1_")
  names_mu0 <- paste0(c("N","IPW","EW1","EW2"),".mu0_")
  delta1 <- mu[,paste0(names_delta,1)] %>% tidyr::gather(Method,RMST,1:4) %>%
    mutate(Method = factor(Method, levels = paste0(names_delta,1),
                           labels = c("Naive","IPSW","CW1","CW2")))
  delta2 <- mu[,paste0(names_delta,2)] %>% tidyr::gather(Method,RMST,1:4) %>%
    mutate(Method = factor(Method, levels = paste0(names_delta,2),
                           labels = c("Naive","IPSW","CW1","CW2")))
  delta3 <- mu[,paste0(names_delta,3)] %>% tidyr::gather(Method,RMST,1:4) %>%
    mutate(Method = factor(Method, levels = paste0(names_delta,3),
                           labels = c("Naive","IPSW","CW1","CW2")))
  mu1_1 <- mu[,paste0(names_mu1,1)] %>% tidyr::gather(Method,RMST,1:4) %>%
    mutate(Method = factor(Method, levels = paste0(names_mu1,1),
                           labels = c("Naive","IPSW","CW1","CW2")))
  mu1_2 <- mu[,paste0(names_mu1,2)] %>% tidyr::gather(Method,RMST,1:4) %>%
    mutate(Method = factor(Method, levels = paste0(names_mu1,2),
                           labels = c("Naive","IPSW","CW1","CW2")))
  mu1_3 <- mu[,paste0(names_mu1,3)] %>% tidyr::gather(Method,RMST,1:4) %>%
    mutate(Method = factor(Method, levels = paste0(names_mu1,3),
                           labels = c("Naive","IPSW","CW1","CW2")))
  mu0_1 <- mu[,paste0(names_mu0,1)] %>% tidyr::gather(Method,RMST,1:4) %>%
    mutate(Method = factor(Method, levels = paste0(names_mu0,1),
                           labels = c("Naive","IPSW","CW1","CW2")))
  mu0_2 <- mu[,paste0(names_mu0,2)] %>% tidyr::gather(Method,RMST,1:4) %>%
    mutate(Method = factor(Method, levels = paste0(names_mu0,2),
                           labels = c("Naive","IPSW","CW1","CW2")))
  mu0_3 <- mu[,paste0(names_mu0,3)] %>% tidyr::gather(Method,RMST,1:4) %>%
    mutate(Method = factor(Method, levels = paste0(names_mu0,3),
                           labels = c("Naive","IPSW","CW1","CW2")))
  return(list(delta1 = delta1, delta2 = delta2, delta3 = delta3,
              mu1_1 = mu1_1, mu1_2 = mu1_2, mu1_3 = mu1_3,
              mu0_1 = mu0_1, mu0_2 = mu0_2, mu0_3 = mu0_3))
}







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

