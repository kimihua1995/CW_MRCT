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




# Naive Estimator
Naive.Est <- function(nR,df,tau){
  N.mu <- N.sd <- NULL
  for (r in 1:nR){
    dfr <- df[df$R==r,]
    nr <- nrow(dfr)
    pr <- rep(1,nr)/nr
    mur1 <- my_akm_rmst(dfr$Y[dfr$A==1], dfr$status[dfr$A==1],pr[dfr$A==1],tau)
    mur0 <- my_akm_rmst(dfr$Y[dfr$A==0], dfr$status[dfr$A==0],pr[dfr$A==0],tau)
    deltar <- mur1$mu - mur0$mu
    #lratior <- log(mur1$mu / mur0$mu)
    N.mu <- c(N.mu,mur1$mu,mur0$mu,deltar)
    deltar.sd <- sqrt(mur1$V + mur0$V)
    #lratior.sd <- sqrt(mur1$V/(mur1$mu)^2 + mur0$V/(mur0$mu)^2)
    N.sd <- c(N.sd,sqrt(mur1$V),sqrt(mur0$V),deltar.sd)
  }
  
  return(list(mu = N.mu, sd = N.sd, p = rep(1,nrow(df))))
}


# IPW Estimator
IPW.Est <- function(nR,nX,df,tau){
  require(nnet)
  IPW.mu <- IPW.sd <- NULL
  varX <- paste0("X",1:nX)
  name_formula <- paste0("R ~ ",paste(varX, collapse = " + "))
  fit <- multinom(formula=name_formula, data = df, trace = F)
  if (nR > 2){
    fit.values <- fit$fitted.values
  }else if (nR == 2){
    fit.values <- cbind(fit$fitted.values, 1 - fit$fitted.values)
  }
  fit.values[fit.values < 0.001] <- 0.001
  fit.values[fit.values > 0.999] <- 0.999
  p <- 1/fit.values
  for (r in 1:nR){
    dfr <- df[df$R==r,]
    pr <- p[df$R==r,r]
    mur1 <- my_akm_rmst(dfr$Y[dfr$A==1], dfr$status[dfr$A==1],pr[dfr$A==1],tau)
    mur0 <- my_akm_rmst(dfr$Y[dfr$A==0], dfr$status[dfr$A==0],pr[dfr$A==0],tau)
    deltar <- mur1$mu - mur0$mu
    #lratior <- log(mur1$mu / mur0$mu)
    IPW.mu <- c(IPW.mu,mur1$mu,mur0$mu,deltar)
    deltar.sd <- sqrt(mur1$V + mur0$V)
    #lratior.sd <- sqrt(mur1$V/(mur1$mu)^2 + mur0$V/(mur0$mu)^2)
    IPW.sd <- c(IPW.sd,sqrt(mur1$V),sqrt(mur0$V),deltar.sd)
  }
  
  return(list(mu = IPW.mu, sd = IPW.sd, p = p))
}


# CW Estimator
CW.Est <- function(nR,nX,X,df,tau){
  require(Rsolnp)
  fn <- function(p){
    sum(p*log(p))
  }
  
  M1 <- M2 <- NULL
  varX <- paste0("X",1:nX)
  for (i in 1:nX){
    M1 <- c(M1, mean(df[,varX[i]]))
    M2 <- c(M2, mean(df[,varX[i]]^2))
  }  
  
  eqn <- function(p){
    c <- NULL
    for (r in 1:nR){
      dfr <- df[df$R==r,]
      cx <- NULL
      for (i in 1:nX){
        if (X[i] == "C"){
          cr1 <- sum(p[df$R==r]*dfr[,varX[i]]) - M1[i]
          cr2 <- sum(p[df$R==r]*dfr[,varX[i]]^2) - M2[i]
          cx <- c(cx, cr1, cr2)
        }else if (X[i] == "B"){
          cr1 <- sum(p[df$R==r]*dfr[,varX[i]]) - M1[i]
          cx <- c(cx, cr1)
        }
      }
      cr3 <- sum(p[df$R==r]) - 1
      c <- c(c,cx,cr3)
    }
    return(c)
  }
  
  constraints <- rep(0,(2*nX - sum(X == "B") + 1)*nR)
  n <- nrow(df)
  p <- rep(1,n)/n
  sol <- solnp(pars = p, fun = fn, eqfun = eqn, eqB = constraints, 
               LB = rep(0,n), UB = rep(1,n), control = list(trace=0, tol=1e-5))
  p <- sol$pars
  EW.mu <- EW.sd <- NULL
  for (r in 1:nR){
    dfr <- df[df$R==r,]
    pr <- p[df$R == r]
    mur1 <- my_akm_rmst(dfr$Y[dfr$A==1], dfr$status[dfr$A==1],pr[dfr$A==1],tau)
    mur0 <- my_akm_rmst(dfr$Y[dfr$A==0], dfr$status[dfr$A==0],pr[dfr$A==0],tau)
    deltar <- mur1$mu - mur0$mu
    #lratior <- log(mur1$mu / mur0$mu)
    EW.mu <- c(EW.mu,mur1$mu,mur0$mu,deltar)
    deltar.sd <- sqrt(mur1$V + mur0$V)
    #lratior.sd <- sqrt(mur1$V/(mur1$mu)^2 + mur0$V/(mur0$mu)^2)
    EW.sd <- c(EW.sd,sqrt(mur1$V),sqrt(mur0$V),deltar.sd)
  }
  
  return(list(mu = EW.mu, sd = EW.sd, p = p))
}





# test for consistency
region.cons.test <- function(fit, nR){
  rmst_diff <- matrix(fit$mu[3*(1:nR)], ncol = 1)
  rmst_diff_sd <- fit$sd[3*(1:nR)]
  
  V <- diag(rmst_diff_sd^2)
  C <- matrix(0, nrow = nR-1, ncol = nR)
  C[,1] <- -1
  for (i in 1:(nR-1)) {C[i,i+1] <- 1}
  
  Q <- t(C %*% rmst_diff) %*% solve(C %*% V %*% t(C)) %*% (C %*% rmst_diff)
  
  p <- 1 - pchisq(Q, nR-1)
  
  return(as.vector(p))
}




# calculate regional difference: delta i - delta j
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
