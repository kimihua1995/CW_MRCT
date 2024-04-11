# calculating calibration weights
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




#################################
#1) Functions for weighted Kaplan Meier estimator
#################################
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

# Weighted KM Estimator
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


#################################
#2) Functions for weighted G-formula estimator
#################################
my.func_surv <- function(y, d){
  #--input--
  #y=time
  #d=status
  
  #--
  id=order(y)
  y=y[id]
  d=d[id]
  
  #--
  t_idx = unique(c(0,y))
  ny = length(y)
  
  #--
  Y = N = C = S = H = D = E = rep(0,length(t_idx))
  
  #i=1
  Y[1] = ny
  N[1] = 0
  C[1] = 0
  S[1] = 1
  H[1] = 0
  D[1] = 0
  E[1] = 0
  
  #i>=2
  for(i in 2:length(t_idx)){
    Y[i] = Y[i-1] - N[i-1] - C[i-1]
    N[i] = ifelse(sum(y==t_idx[i] & d==1)>0, sum(y==t_idx[i] & d==1), 0)
    C[i] = ifelse(sum(y==t_idx[i] & d==0)>0, sum(y==t_idx[i] & d==0), 0)
    
    if(Y[i]<0){Y[i] = 0}
    
    S[i] = ifelse(Y[i]==0, S[i-1], S[i-1]*(1-(N[i]/Y[i])))
    H[i] = ifelse(Y[i]*(Y[i]-N[i])==0, 0, N[i]/(Y[i]*(Y[i]-N[i])))
    
    if(S[i]<0){S[i] = 0}
    
    D[i] = sum(H[2:i])
    E[i] = sqrt((S[i]**2)*D[i])
    
    if(is.na(S[i])){S[i] = 0}
    if(is.na(E[i])){E[i] = 0}
  }
  
  #--output--
  out           = as.data.frame(cbind(t_idx, Y, N, C, S, E))
  colnames(out) = c("t_idx", "n_risk", "n_event", "n_censor", "surv", "se")
  
  #--to match the output of survfit--
  out2 = out[t_idx!=0,]
  
  #--
  Z2 = list()
  Z2$out      = out2
  Z2$t_idx    = out2[,"t_idx"]
  Z2$n_risk   = out2[,"n_risk"]
  Z2$n_event  = out2[,"n_event"]
  Z2$n_censor = out2[,"n_censor"]
  Z2$surv     = out2[,"surv"]
  Z2$se       = out2[,"se"]
  
  return(Z2)
}



my.rmst2reg=function(y, delta, arm, x, tau, w=rep(1,length(y))){

    n=length(y)
    x=as.matrix(cbind(1, x))
    p=length(x[1,])
    
    y0=pmin(y, tau)
    d0=delta
    d0[y0==tau]=1
    
    d10=d0[arm==1]
    d00=d0[arm==0]
    y10=y0[arm==1]
    y00=y0[arm==0]
    x1=x[arm==1,]
    x0=x[arm==0,]
    n1=length(d10)
    n0=length(d00)
    
    id1=order(y10)
    y10=y10[id1]
    d10=d10[id1]
    x1=x1[id1,]
    
    id0=order(y00)
    y00=y00[id0]
    d00=d00[id0]
    x0=x0[id0,]
    
    fitc1=my.func_surv(y10, 1-d10)
    fitc0=my.func_surv(y00, 1-d00)
    
    weights1=d10/rep(pmax(fitc1$surv,0.001), table(y10))
    weights0=d00/rep(pmax(fitc0$surv,0.001), table(y00))
    
    w1=w[arm==1]
    w0=w[arm==0]
    w1=w1[id1]
    w0=w0[id0]
    weights=c(weights1, weights0)*c(w1,w0)
    

    fitt=lm(c(y10,y00)~ rbind(x1, x0)-1, weights=weights)
    
    return(fitt)
}

# IPCW RMST regression (correct model)
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

# IPCW RMST regression (mis-specified model)
fit.rmst.reg.mis <- function(df, tau){
  cov <- data.frame(A=df$A, X1=df$X1, AX1=df$A*df$X1)
  rmst_fit <- my.rmst2reg(y = df$Y,
                          delta = df$status,
                          x = cov,
                          arm = df$A,
                          tau = tau)
  return(rmst_fit)
}


# weighted GF estimator (correct outcome model)
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



# weighted GF estimator (mis-specified outcome model)
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


######################################
#3) Functions for weighted Hajek estimator
######################################
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






######################################
#3) Functions for weighted Augmented estimator
######################################
# weighted AG estimator (correct outcome model)
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



# weighted AG estimator (mis-specified outcome model)
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












