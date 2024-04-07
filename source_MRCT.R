# Generate Regional Data
gen.dat1 <- function(n,r,a,lambda,gamma){
  # seed = 12321
  # n = 100
  # r = matrix(c(1,-1,-1,1,1,-1), nrow = 2, ncol = 3)
  # a = c(0.1,-0.1)*1:13
  # lambda = 0.3
  # gamma = 1
  # tcen = 5
  #set.seed(seed)
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
  # log normal
  #epsilon <- rnorm(n,mean=0,sd=0.1)
  #time_e <- c(exp(epsilon - MM %*% a))
  # Exp
  u <- runif(sum(n))
  time_e <- (-log(u)*exp(- MM %*% a)/lambda[A+1])^(1/gamma[A+1])
  #time_c <- runif(sum(n),1,tcen)
  time_c <- rexp(sum(n))/0.1
  #time_c <- 10
  time <- pmin(time_e,time_c)
  status <- as.numeric(time_e < time_c)
  df$Y <- time
  df$status <- status
  
  
  return(as.data.frame(df %>% group_by(R) %>% arrange(Y, .by_group = T)))
}



gen.dat2 <- function(n,r,a,lambda,gamma){
  # seed = 12321
  # n = 100
  # r = matrix(c(1,-1,-1,1,1,-1), nrow = 2, ncol = 3)
  # a = c(0.1,-0.1)*1:13
  # lambda = 0.3
  # gamma = 1
  # tcen = 5
  #set.seed(seed)
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
  # log normal
  #epsilon <- rnorm(n,mean=0,sd=0.1)
  #time_e <- c(exp(epsilon - MM %*% a))
  # Exp
  u <- runif(sum(n))
  time_e <- (-log(u)*exp(- MM %*% a)/lambda[A+1])^(1/gamma[A+1])
  #time_c <- runif(sum(n),1,tcen)
  time_c <- rexp(sum(n))/0.1
  #time_c <- 10
  time <- pmin(time_e,time_c)
  status <- as.numeric(time_e < time_c)
  df$Y <- time
  df$status <- status
  
  
  return(as.data.frame(df %>% group_by(R) %>% arrange(Y, .by_group = T)))
}



# EW weights
cw1 <- function(df,nR,g){
  p <- NULL
  for (r in 1:nR){
    dfr <- df[df$R == r,]
    fnr <- function(x){
      y1 <- sum(exp(x[1]*dfr$X1 + x[2]*dfr$X2) * 
                  (dfr$X1 - g[1]))
      y2 <- sum(exp(x[1]*dfr$X1 + x[2]*dfr$X2) * 
                  (dfr$X2 - g[2]))
      return(c(y1, y2)) 
    }
    sol_lambdar <- nleqslv(c(1,1), fnr,
                           control = list(btol = 1e-5))
    lambdar <- sol_lambdar$x
    pr <- exp(lambdar[1]*dfr$X1 + lambdar[2]*dfr$X2)
    pr <- pr/sum(pr)
    p <- c(p, pr)
  }
  return(p)
}


cw2 <- function(df,nR,g){
  p <- NULL
  for (r in 1:nR){
    dfr <- df[df$R == r,]
    fnr <- function(x){
      y1 <- sum(exp(x[1]*dfr$X1 + x[2]*dfr$X2 + 
                    x[3]*dfr$X1^2 + x[4]*dfr$X2^2) * 
                  (dfr$X1 - g[1]))
      y2 <- sum(exp(x[1]*dfr$X1 + x[2]*dfr$X2 + 
                  x[3]*dfr$X1^2 + x[4]*dfr$X2^2) * 
                  (dfr$X2 - g[2]))
      y3 <- sum(exp(x[1]*dfr$X1 + x[2]*dfr$X2 + 
                      x[3]*dfr$X1^2 + x[4]*dfr$X2^2) * 
                  (dfr$X1^2 - g[3]))
      y4 <- sum(exp(x[1]*dfr$X1 + x[2]*dfr$X2 + 
                      x[3]*dfr$X1^2 + x[4]*dfr$X2^2) * 
                  (dfr$X2^2 - g[4]))
      return(c(y1, y2, y3, y4)) 
    }
    sol_lambdar <- nleqslv(c(0.5,0.5,0.1,0.1), fnr,
                           control = list(btol = 1e-5))
    lambdar <- sol_lambdar$x
    pr <- exp(lambdar[1]*dfr$X1 + lambdar[2]*dfr$X2 +
                lambdar[3]*dfr$X1^2 + lambdar[4]*dfr$X2^2)
    pr <- pr/sum(pr)
    p <- c(p, pr)
  }
  return(p)
}





# weighted Kaplan Meier estimator
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


W.KM.Est <- function(df,nR,p,tau){
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



# weighted G-formula estimator
fit.rmst.reg <- function(df, tau){
  cov <- data.frame(A=df$A, X1=df$X1, X2=df$X2,
                    AX1=df$A*df$X1, AX2=df$A*df$X2)
  rmst_fit <- my.rmst2reg(y = df$Y,
                          delta = df$status,
                          x = cov,
                          arm = df$A,
                          tau = tau)
  return(rmst_fit)
}


fit.rmst.reg.mis <- function(df, tau){
  cov <- data.frame(A=df$A, X1=df$X1, AX1=df$A*df$X1)
  rmst_fit <- my.rmst2reg(y = df$Y,
                          delta = df$status,
                          x = cov,
                          arm = df$A,
                          tau = tau)
  return(rmst_fit)
}



W.GF.Est <- function(df, nR, p, tau){
  EW.mu <- EW.sd <- NULL
  for (r in 1:nR){
    dfr <- df[df$R==r,]
    pr <- p[df$R == r]
    fitr <- fit.rmst.reg(dfr, tau)
    gammar <- coef(fitr)
    vgammar <- vcov(fitr)
    mur1 <- sum(pr*(gammar[1]+gammar[2]+(gammar[3]+gammar[5])*dfr$X1+
                   (gammar[4]+gammar[6])*dfr$X2))/sum(pr)
    mur0 <- sum(pr*(gammar[1]+gammar[3]*dfr$X1+gammar[4]*dfr$X2))/sum(pr)
    deltar <- mur1 - mur0
    Jr1 <- matrix(c(sum(pr),sum(pr),sum(pr*dfr$X1),sum(pr*dfr$X2),
                    sum(pr*dfr$X1),sum(pr*dfr$X2)), ncol = 1)/sum(pr)
    Jr0 <- matrix(c(sum(pr),sum(pr*dfr$X1),sum(pr*dfr$X2)), ncol = 1)/sum(pr)
    deltar.J <- matrix(c(sum(pr),sum(pr*dfr$X1),sum(pr*dfr$X2)), ncol = 1)/sum(pr)
    sdr1 <- as.numeric(sqrt(t(Jr1) %*% vgammar %*% Jr1))
    sdr0 <- as.numeric(sqrt(t(Jr0) %*% vgammar[c(1,3,4),c(1,3,4)] %*% Jr0))
    deltar.sd <- as.numeric(sqrt(t(deltar.J) %*% vgammar[c(2,5,6),c(2,5,6)] %*% deltar.J))
    
    EW.mu <- c(EW.mu,mur1,mur0,deltar)
    EW.sd <- c(EW.sd,sdr1,sdr0,deltar.sd)
  }
  
  return(list(mu = EW.mu, sd = EW.sd))
}




W.GF.Est.mis <- function(df, nR, p, tau){
  EW.mu <- EW.sd <- NULL
  for (r in 1:nR){
    dfr <- df[df$R==r,]
    pr <- p[df$R == r]
    fitr <- fit.rmst.reg.mis(dfr, tau)
    gammar <- coef(fitr)
    vgammar <- vcov(fitr)
    mur1 <- sum(pr*(gammar[1]+gammar[2]+(gammar[3]+gammar[4])*dfr$X1))/sum(pr)
    mur0 <- sum(pr*(gammar[1]+gammar[3]*dfr$X1))/sum(pr)
    deltar <- mur1 - mur0
    Jr1 <- matrix(c(sum(pr),sum(pr),sum(pr*dfr$X1),sum(pr*dfr$X1)), ncol = 1)/sum(pr)
    Jr0 <- matrix(c(sum(pr),sum(pr*dfr$X1)), ncol = 1)/sum(pr)
    deltar.J <- matrix(c(sum(pr),sum(pr*dfr$X1)), ncol = 1)/sum(pr)
    sdr1 <- as.numeric(sqrt(t(Jr1) %*% vgammar %*% Jr1))
    sdr0 <- as.numeric(sqrt(t(Jr0) %*% vgammar[c(1,3),c(1,3)] %*% Jr0))
    deltar.sd <- as.numeric(sqrt(t(deltar.J) %*% vgammar[c(2,4),c(2,4)] %*% deltar.J))
    
    EW.mu <- c(EW.mu,mur1,mur0,deltar)
    EW.sd <- c(EW.sd,sdr1,sdr0,deltar.sd)
  }
  
  return(list(mu = EW.mu, sd = EW.sd))
}


# weighted Hajek estimator
W.HJ.Est <- function(df, nR, p, tau){
  m_fun <- function(data){
    p <- data$p
    A <- data$A
    Y <- data$Y
    w <- data$w
    function(theta){
      c(p*A*w*(Y-theta[1]),
        p*(1-A)*w*(Y-theta[2]))
    }
  }
  
  EW.mu <- EW.sd <- NULL
  for (r in 1:nR){
    dfr <- df[df$R == r,]
    pr <- p[df$R == r]
    Ar <- dfr$A
    yr <- pmin(dfr$Y, tau)
    dr <- dfr$status
    dr[yr==tau]=1

    dr1=dr[Ar==1]; dr0=dr[Ar==0]
    yr1=yr[Ar==1]; yr0=yr[Ar==0]
    
    fitr1=my.func_surv(yr1, 1-dr1)
    fitr0=my.func_surv(yr0, 1-dr0)
    
    wr1=dr1/rep(pmax(fitr1$surv,0.001), table(yr1))
    wr0=dr0/rep(pmax(fitr0$surv,0.001), table(yr0))
    
    
    dat_modr <- data.frame(p = c(pr[Ar==1], pr[Ar==0]),
                           w = c(wr1, wr0),
                           Y = c(yr1, yr0),
                           A = c(Ar[Ar==1], Ar[Ar==0]))
    
    
    #mur1 <- sum(dat_modr$A*dat_modr$w*dat_modr$p*dat_modr$Y)/
    #        sum(dat_modr$A*dat_modr$w*dat_modr$p)
    #mur0 <- sum((1-dat_modr$A)*dat_modr$w*dat_modr$p*dat_modr$Y)/
    #        sum((1-dat_modr$A)*dat_modr$w*dat_modr$p)
    #deltar <- mur1 - mur0
    
    resr <- geex::m_estimate(
      estFUN = m_fun,
      data = dat_modr,
      root_control = setup_root_control(start = c(0.5,0.5))
    )
    
    coefr <- coef(resr)
    mur1 <- coefr[1]
    mur0 <- coefr[2]
    deltar <- mur1 - mur0
    
    vcovr <- vcov(resr)
    sdr1 <- sqrt(vcovr[1,1])
    sdr0 <- sqrt(vcovr[2,2])
    c1 = matrix(c(1,-1), nrow = 1, ncol = 2)
    deltar.sd <- as.numeric(sqrt(c1%*%vcovr%*%t(c1)))
    
    EW.mu <- c(EW.mu,mur1,mur0,deltar)
    EW.sd <- c(EW.sd,sdr1,sdr0,deltar.sd)
  }
  
  return(list(mu = EW.mu, sd = EW.sd))
}






# weighted Augmented estimator
W.AG.Est <- function(df, nR, p, tau){
  m_fun <- function(data){
    p <- data$p
    A <- data$A
    Y <- data$Y
    w <- data$w
    mu1 <- data$mu1
    mu0 <- data$mu0
    function(theta){
      c(p*A*w*(Y-mu1-theta[1]),
        p*(1-A)*w*(Y-mu0-theta[2]),
        p*(mu1-theta[3]),
        p*(mu0-theta[4]))
    }
  }
  
  

  EW.mu <- EW.sd <- NULL
  for (r in 1:nR){
    dfr <- df[df$R == r,]
    pr <- p[df$R == r]
    Ar <- dfr$A
    yr <- pmin(dfr$Y, tau)
    dr <- dfr$status
    dr[yr==tau]=1
    
    dr1=dr[Ar==1]; dr0=dr[Ar==0]
    yr1=yr[Ar==1]; yr0=yr[Ar==0]
    
    fitr1=my.func_surv(yr1, 1-dr1)
    fitr0=my.func_surv(yr0, 1-dr0)
    
    wr1=dr1/rep(pmax(fitr1$surv,0.001), table(yr1))
    wr0=dr0/rep(pmax(fitr0$surv,0.001), table(yr0))
    
    
    fitr <- fit.rmst.reg(dfr, tau)
    gammar <- coef(fitr)
    vgammar <- vcov(fitr)
    mr1 <-  gammar[1]+gammar[2]+(gammar[3]+gammar[5])*dfr$X1+
            (gammar[4]+gammar[6])*dfr$X2
    mr0 <-  gammar[1]+gammar[3]*dfr$X1+gammar[4]*dfr$X2
    
    dat_modr <- data.frame(p = c(pr[Ar==1], pr[Ar==0]),
                           w = c(wr1, wr0),
                           Y = c(yr1, yr0),
                           A = c(Ar[Ar==1], Ar[Ar==0]),
                           mu1 = c(mr1[Ar==1], mr1[Ar==0]),
                           mu0 = c(mr0[Ar==1], mr0[Ar==0]))
    
    #mur1 <- sum(dat_modr$A*dat_modr$w*dat_modr$p*(dat_modr$Y-dat_modr$mu1))/
    #  sum(dat_modr$A*dat_modr$w*dat_modr$p) +
    #    sum(dat_modr$p*dat_modr$mu1)/sum(dat_modr$p)
    #mur0 <- sum((1-dat_modr$A)*dat_modr$w*dat_modr$p*(dat_modr$Y-dat_modr$mu0))/
    #  sum((1-dat_modr$A)*dat_modr$w*dat_modr$p) +
    #    sum(dat_modr$p*dat_modr$mu0)/sum(dat_modr$p)
    #deltar <- mur1 - mur0
    
    resr <- geex::m_estimate(
      estFUN = m_fun,
      data = dat_modr,
      root_control = setup_root_control(start = c(0.5,0.5,0.5,0.5))
    )
    
    coefr <- coef(resr)
    mur1 <- coefr[1] + coefr[3]
    mur0 <- coefr[2] + coefr[4]
    deltar <- mur1 - mur0
    
    vcovr <- vcov(resr)
    c1 <- matrix(c(1,1), nrow = 1, ncol = 2)
    c2 <- matrix(c(1,-1,1,-1), nrow = 1, ncol = 4)
    sdr1 <- as.numeric(sqrt(c1%*%vcovr[c(1,3),c(1,3)]%*%t(c1)))
    sdr0 <- as.numeric(sqrt(c1%*%vcovr[c(2,4),c(2,4)]%*%t(c1)))
    deltar.sd <- as.numeric(sqrt(c2%*%vcovr%*%t(c2)))
    
    EW.mu <- c(EW.mu,mur1,mur0,deltar)
    EW.sd <- c(EW.sd,sdr1,sdr0,deltar.sd)
  }
  
  return(list(mu = EW.mu, sd = EW.sd))
}




W.AG.Est.mis <- function(df, nR, p, tau){
  m_fun <- function(data){
    p <- data$p
    A <- data$A
    Y <- data$Y
    w <- data$w
    mu1 <- data$mu1
    mu0 <- data$mu0
    function(theta){
      c(p*A*w*(Y-mu1-theta[1]),
        p*(1-A)*w*(Y-mu0-theta[2]),
        p*(mu1-theta[3]),
        p*(mu0-theta[4]))
    }
  }
  
  
  
  EW.mu <- EW.sd <- NULL
  for (r in 1:nR){
    dfr <- df[df$R == r,]
    pr <- p[df$R == r]
    Ar <- dfr$A
    yr <- pmin(dfr$Y, tau)
    dr <- dfr$status
    dr[yr==tau]=1
    
    dr1=dr[Ar==1]; dr0=dr[Ar==0]
    yr1=yr[Ar==1]; yr0=yr[Ar==0]
    
    fitr1=my.func_surv(yr1, 1-dr1)
    fitr0=my.func_surv(yr0, 1-dr0)
    
    wr1=dr1/rep(pmax(fitr1$surv,0.001), table(yr1))
    wr0=dr0/rep(pmax(fitr0$surv,0.001), table(yr0))
    
    
    fitr <- fit.rmst.reg.mis(dfr, tau)
    gammar <- coef(fitr)
    vgammar <- vcov(fitr)
    mr1 <-  gammar[1]+gammar[2]+(gammar[3]+gammar[4])*dfr$X1
    mr0 <-  gammar[1]+gammar[3]*dfr$X1
    
    dat_modr <- data.frame(p = c(pr[Ar==1], pr[Ar==0]),
                           w = c(wr1, wr0),
                           Y = c(yr1, yr0),
                           A = c(Ar[Ar==1], Ar[Ar==0]),
                           mu1 = c(mr1[Ar==1], mr1[Ar==0]),
                           mu0 = c(mr0[Ar==1], mr0[Ar==0]))
    
    #mur1 <- sum(dat_modr$A*dat_modr$w*dat_modr$p*(dat_modr$Y-dat_modr$mu1))/
    #  sum(dat_modr$A*dat_modr$w*dat_modr$p) +
    #    sum(dat_modr$p*dat_modr$mu1)/sum(dat_modr$p)
    #mur0 <- sum((1-dat_modr$A)*dat_modr$w*dat_modr$p*(dat_modr$Y-dat_modr$mu0))/
    #  sum((1-dat_modr$A)*dat_modr$w*dat_modr$p) +
    #    sum(dat_modr$p*dat_modr$mu0)/sum(dat_modr$p)
    #deltar <- mur1 - mur0
    
    resr <- geex::m_estimate(
      estFUN = m_fun,
      data = dat_modr,
      root_control = setup_root_control(start = c(0.5,0.5,0.5,0.5))
    )
    
    coefr <- coef(resr)
    mur1 <- coefr[1] + coefr[3]
    mur0 <- coefr[2] + coefr[4]
    deltar <- mur1 - mur0
    
    vcovr <- vcov(resr)
    c1 <- matrix(c(1,1), nrow = 1, ncol = 2)
    c2 <- matrix(c(1,-1,1,-1), nrow = 1, ncol = 4)
    sdr1 <- as.numeric(sqrt(c1%*%vcovr[c(1,3),c(1,3)]%*%t(c1)))
    sdr0 <- as.numeric(sqrt(c1%*%vcovr[c(2,4),c(2,4)]%*%t(c1)))
    deltar.sd <- as.numeric(sqrt(c2%*%vcovr%*%t(c2)))
    
    EW.mu <- c(EW.mu,mur1,mur0,deltar)
    EW.sd <- c(EW.sd,sdr1,sdr0,deltar.sd)
  }
  
  return(list(mu = EW.mu, sd = EW.sd))
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
  true.mu <- as.data.frame(matrix(true.mu,nrow=1))
  colnames(true.mu) <-  paste0(rep(c("mu1","mu0","delta"),nR),"_",rep(1:nR,each=3))
                          
  
  #for (j in 1:(nR-1)){
  #  true.mu <- c(true.mu, true.mu[3*((j+1):nR)] - true.mu[3*j])
  #}
  
  #true.mu <- as.data.frame(matrix(true.mu,nrow=1))
  #region.diff.pair <- NULL
  #for (j in 1:(nR-1)) {
  #  region.diff.pair <- c(region.diff.pair, paste0(((j+1):nR),j))
  #}
  #colnames(true.mu) <-  c(paste0(rep(c("mu1","mu0","delta"),nR),"_",rep(1:nR,each=3)),
  #                        paste0("delta_",region.diff.pair))
  
  
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
    #gen_dat_success <- FALSE
    #while (!gen_dat_success){
    #  skip_to_next <- FALSE
    #tryCatch(df <- gen.dat(n,r,a,lambda,gamma),
    #         error = function(e){skip_to_next <<- TRUE})
    #if (skip_to_next) {next} else {gen_dat_success <- TRUE} 
    #}
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
    #p3 <- p3/sum(p3)
    
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















