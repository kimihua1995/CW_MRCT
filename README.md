# Inference of treatment effect and its regional modifiers using the restricted mean survival time in multi-regional clinical trials
This project contains the R codes for simulation study in the paper

source_estimator.R
```
This code provides functions for generating weighted estimator of region-specific average RMST difference

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


2) function "W.Est()", calculating the weighted adjusted region-specific RMST
Arguments:
df: a data frame including:
  R = region indicator from 1 to nR
  A = treatment indicator where 1 denotes experimental treatment group and 0 denotes control group
  Y = time-to-event
  status = censoring indicator
nR: number of regions
p: weight
tau: time horizon used for RMST

Values:
a list including:
mu = regional weighted RMST for each treatment group and region-specific weighted RMST difference
sd = the standard deviation of mu




3) function "EW.p.M()", calculating the calibration weights
Arguments:
df: a data frame including:
  R = region indicator from 1 to nR
  A = treatment indicator where 1 denotes experimental treatment group and 0 denotes control group
  Y = time-to-event
  status = censoring indicator
  X1-Xp = covariates named by from "X1" to "Xp", where p=nX
nR: number of regions
f.list: a list of functions of "g()" in Equation (3.4)
iX: vector of denoting the which covariate be used in f.list
nX: number of covariates
M: values of \Tilte(g) in Equation (3.4)

Values:
p = calibration weights


4) function "region.diff()", calculating the weighted adjusted region-specific RMST difference
Arguments:
res: restuls from W.Est()
nR: number of regions


5) function "region.cons.test()", regional consistency test for treatment effect
Arguments:
res: results from W.Est()
nR: number of regions
```


simulation.R
```
This code provides functions for simulation studies including: 1) MRCT data generation, 2) Simulation comparing 3 methods, 3) True RMST, and 4) result output

1) function "gen.dat1()", generating MRCT data under log-linear sampling setting
Arguments:
n: sample sizes in each region, vector of length nR
r: parameters used in the sampling score model (Equation 5.1), nR-by-3 matrix
a: parameters used in the hazard function (Equation 5.2), vector of length 12
lambda: scale parameters for baseline hazard for two treatment groups, vector of length 2
gamma: shape parameters for baseline hazard for two treatment groups, vector of length 2

Values:
a data frame including:
A = treatment indicator
X = covariates
R = region indicator
time = time-to-event
status = censoring indicator


2) function "gen.dat2()", generating MRCT data under logistic sampling setting
Arguments:
n: sample sizes in each region, vector of length nR
r: parameters used in the sampling score model (Equation 5.3), nR-by-3 matrix
a: parameters used in the hazard function (Equation 5.2), vector of length 12
lambda: scale parameters for baseline hazard for two treatment groups, vector of length 2
gamma: shape parameters for baseline hazard for two treatment groups, vector of length 2

Values:
a data frame including:
A = treatment indicator
X = covariates
R = region indicator
time = time-to-event
status = censoring indicator


3) function "true.rmst()", function for calculating true RMST
Arguments:
tau: time horizon for RMST
nR: number of regions
a: parameters used in the hazard function, vector of length 12
lambda: scale parameters for baseline hazard for two treatment groups, vector of length 2
gamma: shape parameters for baseline hazard for two treatment groups, vector of length 2

Values:
true region-specific RMST and RMST difference



4) function "sim.MRCT.Est()", simulation comparing 3 methods
Arguments:
seed: seed
S: iterations
n: sample sizes in each region, vector of length nR
r: parameters used in the sampling score model, nR-by-3 matrix
a: parameters used in the time-to-event outcome model, vector of length 12
lambda: scale parameters for baseline hazard for two treatment groups, vector of length 2
gamma: shape parameters for baseline hazard for two treatment groups, vector of length 2
tau: time horizon for RMST
M: values of \Tilte(g) in Equation (3.4)
setting: 1 for log-linear sampling and 2 for logistic sampling

Values:
a list including:
mu = mean of weighted region-specific RMST and RMST difference
sd = sd of weighted region-specific RMST and RMST difference










```
