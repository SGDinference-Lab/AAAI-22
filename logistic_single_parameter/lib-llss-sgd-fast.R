# Libraries
#
# Last Update: 2021-08-05
#
# 2021-08-05  The critical value is changed from KVB (2000) to AP (1997). 6.811 -> 6.747



# Required libraries
library(MASS)
library(mvtnorm)
library(tictoc)
# DGP of the Simulations by CLTZ in AOS (2020)
#   model: "linear" or "logit"
#   Sigma: variance of regressor: "Identity", "Teoplitz", or "Equi_Corr"
#   sigma: standard deviation of error term
#   n: sample size
#   iter: number of iterations (replications)
#   d: dimension of regressors
#   bt.star: (d x 1) dimension, true parameter values
#   r: correlation coefficient of Sigma. Default = 0


random.scaling.update = function(t.obs, bar.bt_t, a.old, b.old, c.old){

  a.new = a.old + t.obs^2 * bar.bt_t %*% t(bar.bt_t)
  b.new = b.old + t.obs^2 * bar.bt_t
  c.new = c.old + t.obs^2

  V_t = ( a.new - b.new %*% t(bar.bt_t) - bar.bt_t %*% t(b.new) + c.new * bar.bt_t %*% t(bar.bt_t) ) / (t.obs^2)

  return(list(a.new=a.new, b.new=b.new, c.new=c.new, V_t=V_t))
}




#---------------------------------------------------------------------------------------------------
#
# The codes below are imported and modified from Zhu et al. (2021)
#
#---------------------------------------------------------------------------------------------------

########################## Function for experiments ##########################

#runs one trial of the simulation for inference of mean vector (online BM)
#and saves a sequence of CI coverage/CI length/the loss of covariance matrix estimator
#-------------------------------------------------------------------
# Additional Comments by Youngki (2021-05-25)
#   n: sample size ('obs' in the code is the observation index 'i' or 't')
#   Sigma: (d x d) covariance matrix for covariate a
#   sigma: a scalar standard error for the error term e
#       note that linear model defined as b = a'*x_true + e
#   c: constant for ak (notation in the paper: a_m)
#   alpha: parameter in the learning rate eta_i = eta_0 * i^(-alpha)
#          also used in a_m
#   lr: parameter in the learning rate eta_i. eta_0 part above
#   itenumber: number of iterations (200 replications written in the paper)
#   d: dimension of covariate a
#   q: the significance level (e.g. 0.05 for 5%)
#   x_true: the true parameter value (note that x is parameter in their notation)
#   burn: number of burn-in observations to set the initial value
#   batchtype: "Full" or "NOL" (non-overlapping)
#
#I think this code is not correct.
#   - Xbar_old is not properly defined.
#   - Xbar_old is not recursively updated.
#-------------------------------------------------------------------


# with generated data (faster)
# Main Code that we will use
# data_a and data_b is 3-dim array (ite, n, d) and (ite, n, 1), repspectively.

exp_cor_withdata = function(n, data_a, data_b, x_true, c, alpha, lr, itenumber, q, burn = 1,
                            batchtype = "Full", modeltype = "Linear", simple.output=TRUE)
  {
  start_time = Sys.time()

  # Dimension of regressors
  d = length(x_true)
  # Set a
  ak = as.integer(c*((1:100))**(2/(1-alpha)))
  ak[1] = 1


  # Initialize the output variables
  bar_cilength = array(0, dim = c(n, itenumber, d))
  coverrate = array(0, dim = c(n, itenumber,d))
  X_hat_lossbar = matrix(0, nrow = n, ncol = d)
  X_bar_hat = matrix(0, nrow = n, ncol = d)  #added to check the convergence Youngki
  S_hat_bar = matrix(0, nrow = n, ncol = d)
  runt_ave =  rep(0,n)
  runt = rep(0,itenumber)
  X = matrix(0, nrow = itenumber, ncol = d, byrow = TRUE)
  Xbar_old = matrix(0, nrow = itenumber, ncol = d, byrow = TRUE)

  # n-th results
  X_hat_losslast= matrix(0, nrow = itenumber, ncol = d)
  S_last = matrix(0, nrow = itenumber, ncol = d)
  cilengthlast = matrix(0, nrow = itenumber, ncol = d)
  iscoverlast = matrix(0, nrow = itenumber, ncol = d)


  if (modeltype == "Linear") {
    if (burn > 1) {
      for (obs in 1:(burn-1)){
        lrnew = lr*obs**(-alpha)
        for (ite in 1:itenumber){
          grad = grad_comp(X[ite,], t(data_a[ite, obs, ]), data_b[ite, obs], modeltype)
          X[ite,] = X[ite,] - lrnew*grad
        }
      }
    }


    if (batchtype == "Plugin"){
      A_hat = as.list(rep(0, itenumber))
      H_hat = as.list(rep(0, itenumber))
    }else{
      k_start = min(which(ak>=burn))  # burn-in obs may not be on the threshold.
      k_old = rep(k_start-1, itenumber)
      v_old = rep(0, itenumber)
      q_old = rep(0, itenumber)
      P_old = matrix(0, nrow = itenumber, ncol = d)
      V_old = matrix(0, nrow = itenumber, ncol = d)
      W_old = matrix(0, nrow = itenumber, ncol = d)
      l_old = rep(0, itenumber)

      # parameters for Random Scaling method updates
      a.old = array(0, dim = c(itenumber, d, d))
      b.old = array(0, dim = c(itenumber, d, 1))
      c.old = rep(0, itenumber)
    }
    pb <- txtProgressBar(min = burn, max = n, style = 3)
    for (obs in burn:n){
      #cat("Processing........",obs," out of ",n,"\n")
      setTxtProgressBar(pb, obs)  #progress bar
      lrnew = lr*obs**(-alpha)
      for (ite in 1:itenumber){
        t1 = Sys.time()
        grad = grad_comp(X[ite,], t(data_a[ite, obs, ]), data_b[ite, obs], modeltype)
        X[ite,] = X[ite,] - lrnew*grad
        Xbar_old[ite,] = (Xbar_old[ite,]*(obs - burn) + X[ite,])/(obs - burn+1)

        if (batchtype == "NOL"){
          old = recursiveupdate_cor_NOL(obs, l_old[ite], k_old[ite], ak, v_old[ite], q_old[ite],
                                        P_old[ite,], V_old[ite,], W_old[ite,], Xbar_old[ite,], X[ite,])
          l_old[ite]= old$l
          k_old[ite] = old$k
          v_old[ite] = old$v
          q_old[ite] = old$q
          P_old[ite,] = old$P
          V_old[ite,] = old$V
          W_old[ite,] = old$W
          S_hat = old$S
        }
        if (batchtype == "Full"){
          old = recursiveupdate_cor_F(obs, l_old[ite], k_old[ite], ak, v_old[ite], q_old[ite],
                                      P_old[ite,], V_old[ite,], W_old[ite,], Xbar_old[ite,], X[ite,])
          l_old[ite]= old$l
          k_old[ite] = old$k
          v_old[ite] = old$v
          q_old[ite] = old$q
          P_old[ite,] = old$P
          V_old[ite,] = old$V
          W_old[ite,] = old$W
          S_hat = old$S
        }
        if (batchtype == "Random"){
          rs = random.scaling.update(obs, Xbar_old[ite, 1], a.old[ite, 1, 1], b.old[ite, , ], c.old[ite])
          a.old[ite, 1, 1] = rs$a.new
          b.old[ite, 1, 1] = rs$b.new
          c.old[ite] = rs$c.new
          S_hat = diag(rs$V_t) # Should fix this later for our Design.
        }
        if(batchtype == "Plugin"){
          heis = heis_comp(X[ite, ], data_a[ite,obs, ], data_b[ite, ], modeltype)
          #cat(heis, '\n')
          A_hat[[ite]] = (A_hat[[ite]]*(obs-burn) + heis)/(obs-burn+1)
          r = eigen(A_hat[[ite]])
          #cat(r$vectors, '\n')
          #cat(diag(1/pmax(r$values, 0.01)), '\n')
          #cat(t(r$vectors), '\n')
          inv_A_hat = r$vectors%*%diag(1/pmax(r$values, 0.01))%*%t(r$vectors)
          H_hat[[ite]] =  (H_hat[[ite]]*(obs-burn) + grad%*%t(grad))/(obs-burn+1)
          S_hat = diag(inv_A_hat%*%H_hat[[ite]]%*%inv_A_hat)
        }
        t2 = Sys.time()
        runt[ite] = as.numeric(difftime(t2,t1,units="secs") )

        if (batchtype == "Random"){
          #critical.value = 6.811 # From Kiefer et al. (2000) Table 2. 97.5%
          critical.value = 6.747 # From Abadir and Paruolo (1997) Table 1. 97.5%
          X_hat_losslast[ite,] = (Xbar_old[ite,] - x_true)*sqrt(obs)
          cilength_cor = critical.value*sqrt(S_hat)
          cilengthlast[ite,] = 2*(cilength_cor/sqrt(obs))
          iscoverlast[ite, ] = as.numeric(abs(X_hat_losslast[ite,]) < cilength_cor)
          S_last[ite,] = S_hat
        } else{
          X_hat_losslast[ite,] = (Xbar_old[ite,] - x_true)*sqrt(obs)
          cilength_cor = abs(qnorm(q/2))*sqrt(S_hat)
          cilengthlast[ite,] = 2*(cilength_cor/sqrt(obs))
          iscoverlast[ite, ] = as.numeric(abs(X_hat_losslast[ite,]) < cilength_cor)
          S_last[ite,] = S_hat
        }

      }

      runt_ave[obs] = runt_ave[max(1, obs-1)] + mean(runt)
      bar_cilength[obs, , ] = cilengthlast
      coverrate[obs, , ] = iscoverlast
      X_hat_lossbar[obs,] = colMeans(X_hat_losslast[,])
      X_bar_hat[obs,] = Xbar_old[1,]  # Added to check the convergence (Youngki)
      S_hat_bar[obs,] = colMeans(S_last[,])
    }
  }

  # End of Linear model
  #------------------------------------------------------------------------------------------------------------

  if (modeltype == "logistic") {
    if (burn > 1) {
      for (obs in 1:(burn-1)){
        lrnew = lr*obs**(-alpha)
        for (ite in 1:itenumber){
          grad = grad_comp(X[ite,], t(data_a[ite, obs, ]), data_b[ite, obs], modeltype)
          X[ite,] = X[ite,] - lrnew*grad
        }
      }
    }


    if (batchtype == "Plugin"){
      A_hat = as.list(rep(0, itenumber))
      H_hat = as.list(rep(0, itenumber))
    }else{
      k_start = min(which(ak>=burn))  # burn-in obs may not be on the threshold.
      k_old = rep(k_start-1, itenumber)
      v_old = rep(0, itenumber)
      q_old = rep(0, itenumber)
      P_old = matrix(0, nrow = itenumber, ncol = d)
      V_old = matrix(0, nrow = itenumber, ncol = d)
      W_old = matrix(0, nrow = itenumber, ncol = d)
      l_old = rep(0, itenumber)

      # parameters for Random Scaling method updates
      a.old = array(0, dim = c(itenumber, d, d))
      b.old = array(0, dim = c(itenumber, d, 1))
      c.old = rep(0, itenumber)
    }
    pb <- txtProgressBar(min = burn, max = n, style = 3)
    for (obs in burn:n){
      #cat("Processing........",obs," out of ",n,"\n")
      setTxtProgressBar(pb, obs)  #progress bar
      lrnew = lr*obs**(-alpha)
      for (ite in 1:itenumber){
        t1 = Sys.time()
        grad = grad_comp(X[ite,], t(data_a[ite, obs, ]), data_b[ite, obs], modeltype)
        X[ite,] = X[ite,] - lrnew*grad
        Xbar_old[ite,] = (Xbar_old[ite,]*(obs - burn) + X[ite,])/(obs - burn+1)

        if (batchtype == "NOL"){
          old = recursiveupdate_cor_NOL(obs, l_old[ite], k_old[ite], ak, v_old[ite], q_old[ite],
                                        P_old[ite,], V_old[ite,], W_old[ite,], Xbar_old[ite,], X[ite,])
          l_old[ite]= old$l
          k_old[ite] = old$k
          v_old[ite] = old$v
          q_old[ite] = old$q
          P_old[ite,] = old$P
          V_old[ite,] = old$V
          W_old[ite,] = old$W
          S_hat = old$S
        }
        if (batchtype == "Full"){
          old = recursiveupdate_cor_F(obs, l_old[ite], k_old[ite], ak, v_old[ite], q_old[ite],
                                      P_old[ite,], V_old[ite,], W_old[ite,], Xbar_old[ite,], X[ite,])
          l_old[ite]= old$l
          k_old[ite] = old$k
          v_old[ite] = old$v
          q_old[ite] = old$q
          P_old[ite,] = old$P
          V_old[ite,] = old$V
          W_old[ite,] = old$W
          S_hat = old$S
        }
        if (batchtype == "Random"){
          rs = random.scaling.update(obs, Xbar_old[ite, 1], a.old[ite, 1, 1], b.old[ite, 1, 1], c.old[ite])
          a.old[ite, 1, 1] = rs$a.new
          b.old[ite, 1, 1] = rs$b.new
          c.old[ite] = rs$c.new
          S_hat = diag(rs$V_t) # Should fix this later for our Design.
        }
        if(batchtype == "Plugin"){
          heis = heis_comp(X[ite, ], data_a[ite,obs, ], data_b[ite, ], modeltype)
          #cat(heis, '\n')
          A_hat[[ite]] = (A_hat[[ite]]*(obs-burn) + heis)/(obs-burn+1)
          r = eigen(A_hat[[ite]])
          #cat(r$vectors, '\n')
          #cat(diag(1/pmax(r$values, 0.01)), '\n')
          #cat(t(r$vectors), '\n')
          inv_A_hat = r$vectors%*%diag(1/pmax(r$values, 0.01))%*%t(r$vectors)
          H_hat[[ite]] =  (H_hat[[ite]]*(obs-burn) + grad%*%t(grad))/(obs-burn+1)
          S_hat = diag(inv_A_hat%*%H_hat[[ite]]%*%inv_A_hat)
        }
        t2 = Sys.time()
        runt[ite] = as.numeric(difftime(t2,t1,units="secs") )

        if (batchtype == "Random"){
          #critical.value = 6.811 # From Kiefer et al. (2000) Table 2. 97.5%
          critical.value = 6.747 # From Abadir and Paruolo (1997) Table 1. 97.5%
          X_hat_losslast[ite,] = (Xbar_old[ite,] - x_true)*sqrt(obs)
          cilength_cor = critical.value*sqrt(S_hat)
          cilengthlast[ite,] = 2*(cilength_cor/sqrt(obs))
          iscoverlast[ite, ] = as.numeric(abs(X_hat_losslast[ite,]) < cilength_cor)
          S_last[ite,] = S_hat
        } else{
          X_hat_losslast[ite,] = (Xbar_old[ite,] - x_true)*sqrt(obs)
          cilength_cor = abs(qnorm(q/2))*sqrt(S_hat)
          cilengthlast[ite,] = 2*(cilength_cor/sqrt(obs))
          iscoverlast[ite, ] = as.numeric(abs(X_hat_losslast[ite,]) < cilength_cor)
          S_last[ite,] = S_hat
        }

      }

      runt_ave[obs] = runt_ave[max(1, obs-1)] + mean(runt)
      bar_cilength[obs, , ] = cilengthlast
      coverrate[obs, , ] = iscoverlast
      X_hat_lossbar[obs,] = colMeans(X_hat_losslast[,])
      X_bar_hat[obs,] = Xbar_old[1,]  # Added to check the convergence (Youngki)
      S_hat_bar[obs,] = colMeans(S_last[,])
    }
  }

  # End of logistic regression
  #----------------------------------------------------------------------------------------------------



  end_time = Sys.time()
  runt_all = end_time - start_time

  # Simplified version
  coverrate.1    = apply(coverrate[,,1], 1, mean)
  bar_cilength.1 = apply(bar_cilength[,,1], 1, mean)

  if (simple.output) {
    val = list( covrate_ave = coverrate.1,
                cilengthbar = bar_cilength.1,
                runtime_cum = runt_ave
    )
  } else {
    val = list( covrate_ave = coverrate,  cilengthbar = bar_cilength,
                Xbarloss_ave = X_hat_lossbar,   S_hat_bar = S_hat_bar,
                Sigma_hat_last = S_last, Xbarloss_last = X_hat_losslast,
                runtime_all = runt_all, runtime_cum = runt_ave,xtrue = x_true, c = c, alpha = alpha, lr = lr,
                X_bar_hat = X_bar_hat)
  }

  # Output description by Youngki (2021-05-25)
  #   cilengthbar: (n x 1) ci length averaged over iterations (and different coefficients) at each step.
  #   covrate_ave: (n x 1) coverage rates averaged over "iterations and coefficients" at each step.
  #   S_hat_bar:   (n x d) average diag(Sima) over iterations
  #   cilength_last:  (itenumber x d) final step (at n) results for all iterations and coefficients
  #   covrate_last:   (itenumber x d) final step (at n) results. Either it includes the true value or not. Binary (0 / 1) results.
  #   Sigma_hat_last: (itenumber x d) final step (at n) results. Diagonal elements only.
}


########################## Functions used in experiments for statistical inference  ##########################
#function to compute the gradient for linear/logistic model
grad_comp = function(x_old, a_new, b_new, type = "Linear"){
  a_new = matrix(a_new, d, 1)
  x_old = matrix(x_old, d, 1)
  if (type == "Linear"){
    grad = a_new *c(t(a_new)%*%x_old - b_new)
  }else if (type =="logistic"){
    #grad = a_new*(1/c(1+exp(-t(a_new)%*%x_old))-b_new)
    grad = -c((2*b_new-1)*a_new) / c(1+exp( (2*b_new-1) * (t(a_new)%*%x_old)) )
  }
  return (grad)
}

#function to compute the Hessian for linear/logistic model
heis_comp = function(x_old, a_new, b_new, type = "Linear"){
  a_new = matrix(a_new, d, 1)
  x_old = matrix(x_old, d, 1)
  if (type == "Linear"){
    heis = a_new %*% t(a_new)
  }else if (type =="logistic"){
    heis =  a_new%*%t(a_new)/c(1+exp(t(a_new)%*%x_old))/c(1+exp(-t(a_new)%*%x_old))
  }
  return (heis)
}


#only update diagonals of covariance matrix estimate (full overlapping version)
recursiveupdate_cor_F = function(n_new, l_old, k_old, a, v_old, q_old, P_old, V_old, W_old, Xbar_new, X_new){
  if(a[k_old + 1]==n_new){
    k_new = k_old + 1
    l_new = 1
    W_new = X_new
  }
  else{
    k_new = k_old
    l_new = l_old + 1
    W_new = X_new + W_old
  }
  v_new = v_old + l_new
  q_new = q_old + l_new**2
  V_new = V_old + W_new**2
  P_new = P_old + l_new*W_new
  Sp = V_new - 2*P_new*Xbar_new + q_new * Xbar_new**2
  S = Sp/v_new
  val = list(k = k_new, l = l_new, v = v_new, q = q_new, P = P_new, V = V_new,
             W = W_new, S = S)
}


