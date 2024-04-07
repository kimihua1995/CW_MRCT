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

