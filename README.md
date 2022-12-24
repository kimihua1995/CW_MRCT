# Inference of treatment effect and its regional modifiers using the restricted mean survival time in clinical trials conducted in multiple regions

This project contains the R codes for simulation study for MRCT analysis using RMST by calibration weighting method

source_MRCT.R
```
This code provides functions including 1) weighted RMST (SD), 2) 3 methods, Naive, IPW, and CW, that we compared in the simulation study, 3) regional consistency test for treatment effect.

1) function "my_akm_rmst()", calculating weighted RMST
Arguments:
time: time-to-event
status: censoring indicator
weight: weights
tau: time horizon used for RMST

Values:
a data frame including:
mu = weighted RMST
V = variance of weighted RMST


2) function "Naive.Est()", calculating unweighted RMST
Arguments:
nR: number of regions
df: a data frame including:
  R = region indicator from 1 to nR
  A = treatment indicator where 1 denotes experimental treatment group and 0 denotes control group
  Y = time-to-event
  status = censoring indicator
tau: time horizon used for RMST

Values:
a list including:
mu = regional weighted RMST for each treatment group and regional weighted RMST difference
sd = the standard deviation of mu
p = weights (all 1 for Naive method)



3) function "IPW.Est()", calculating weighted RMST by IPW method
Arguments:
nR: number of regions
nX: number of covariates
df: a data frame including:
  R = region indicator from 1 to nR
  A = treatment indicator where 1 denotes experimental treatment group and 0 denotes control group
  Y = time-to-event
  status = censoring indicator
  X1-Xp = covariates named by from "X1" to "Xp", where p=nX
tau: time horizon used for RMST

Values:
a list including:
mu = regional weighted RMST for each treatment group and regional weighted RMST difference
sd = the standard deviation of mu
p = weights (all 1 for Naive method)




4) function "CW.Est()", calculating weighted RMST by CW method
Arguments:
nR: number of regions
nX: number of covariates
X: type of covariates, "C" for continuous variables and "B" for binary variables
df: a data frame including:
  R = region indicator from 1 to nR
  A = treatment indicator where 1 denotes experimental treatment group and 0 denotes control group
  Y = time-to-event
  status = censoring indicator
  X1-Xp = covariates named by from "X1" to "Xp", where p=nX
tau: time horizon used for RMST

Values:
a list including:
mu = regional weighted RMST for each treatment group and regional weighted RMST difference
sd = the standard deviation of mu
p = weights (all 1 for Naive method)



5) function "region.cons.test()", regional consistency test for treatment effect
Arguments:
fit: results from Naive.Est(), IPW.Est(), or CW.Est()
nR: number of regions
```
