# SGD inference simulations
# 
#   The file is written for a single node. The results reported in LLSS (2021) 
# is conducted at the Compute Canada Graham cluster over 100 nodes. The 1000 
# iterations are divided into 100 nodes (10 iterations each) with seed number 
# [15673 -- 15772]. 
#
#
#
# First date  : 2021-05-28
# Last update : 2021-05-29
#

rm(list=ls())

# Requires packages and libraries
source('lib-llss-sgd.R')
library(MASS)
library(mvtnorm) 
library(tictoc)
library(xtable)


tic()

# Command Line Inputs
# [1] array task id
# [2] d: 05, 20
# [3] alpha: 0.505, 0.667
# [4] lr: 0.5, 1
# [5] model: 01-04
# [6] simple.output: TURE / FALSE
# [7] burn: 1 (no burn-in), 2 (1% of n burn-in)

args = commandArgs(trailingOnly = TRUE)   
job.id = args[1]

seed = as.numeric(job.id)+15672
set.seed(seed)

################################################################################
# d: dimension of regressors
# n: sample size
# Sigma: Covariance of regressors
# sigma: standard deviation of the error term eps
# c: parameter for recursive batch-mean method
# alpha: parameter for the learning rate
# lr: parameter for the learning rate ( learning rate(t) = lr * t^(-alpha) )
# itenumber: number of replication (itenumber >=2 )
# q: the significance level (nominal rate = 1- q)
# bt_0: beta_0, the tre parameter value
# burn: parameter for burn-in observations. 
# modeltype: simulation design
#   "Linear": y = x'beta_0 + eps
################################################################################
d = as.numeric(args[2])
n = 1e05
Sigma = diag(d)
sigma = 1
c = 1
alpha = as.numeric(args[3])
lr = as.numeric(args[4])
itenumber = 10
q = 0.05
bt_0 = seq(0,1,length.out=d)
modeltype = "Linear"
model = args[5]
simple.output = args[6]
burn = as.numeric(args[7])

result.dic = paste0("d-",args[2],"/d-",args[2],"-",model,"/")

file.name = paste0(modeltype,'_d',d,'_m',model,'_',job.id)



# Data generation for all iterations
data_x = array(NA, dim=c(itenumber, n, d))
data_y = array(NA, dim=c(itenumber, n))

for (ite in (1:itenumber)) {
  data_x[ite, , ] = rmvnorm(n, sigma = Sigma)
  data_y[ite, ] = data_x[ite, , ] %*% bt_0 + rnorm(n, mean=0, sd=sigma)
}


r.rs = exp_cor_withdata(n, data_x, data_y, bt_0, c, alpha, lr, itenumber, q, burn, batchtype = "Random", modeltype, simple.output)
r.bm = exp_cor_withdata(n, data_x, data_y, bt_0, c, alpha, lr, itenumber, q, burn, batchtype = "Full", modeltype, simple.output)
r.pi = exp_cor_withdata(n, data_x, data_y, bt_0, c, alpha, lr, itenumber, q, burn, batchtype = "Plugin", modeltype, simple.output)

################################################################################
# Output:
#   r.rs: results for the random scaling method
#   r.bm: results for the batch-mean method
#   r.pi: results for the plug-in method
################################################################################
save(r.rs, r.bm, r.pi, file = paste0(result.dic,file.name,'.RData'))




toc()
