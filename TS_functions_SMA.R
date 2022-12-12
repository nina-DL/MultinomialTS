
# MC estimate of the allocation probabilities ---------------------------------------------------

rho0_Norm = function(M=1000, means, sds){
  # this function provides MC estimates of TS arm1 allocation probability in a 2-arm case
  # M:: number of MC simulations
  # means:: 2-dim vector with posterior means of the two arms at time t
  # sds:: 2-dim vector with posterior standard deviations of the two arms at time t
  mu1 = rnorm(M, means[1], sds[1])
  mu2 = rnorm(M, means[2], sds[2])
  rho_hat = mean(mu1>mu2)
  
  if(rho_hat==1) rho_hat=0.999
  if(rho_hat==0) rho_hat=0.001
  
  return(rho_hat)
}

rho0_Multi = function(M=1000, params1, params2, y_categories = y_categories){
  mu1 = apply(rdirichlet(M, params1), 1, function(p) sum(p*c(y_categories)))
  mu2 = apply(rdirichlet(M, params2), 1, function(p) sum(p*c(y_categories)))
  rho_hat = mean(mu1>mu2)
  
  if(rho_hat==1) rho_hat=0.999
  if(rho_hat==0) rho_hat=0.001
  
  return(rho_hat)
}


rho0_Binary = function(M=1000, params1, params2){
  mu1 = rbeta(M, params1[1], params1[2])
  mu2 = rbeta(M, params2[1], params2[2])
  rho_hat = mean(mu1>mu2)
  
  if(rho_hat==1) rho_hat=0.999
  if(rho_hat==0) rho_hat=0.001
  
  return(rho_hat)
}



# Multinomial TS ---------------------------------------------------

#The higher value of alpha_i, the greater the value of X_i
#and the greater amount of the total "mass" is assigned to it 
#choice of prior params: https://stats.stackexchange.com/questions/483247/what-is-a-non-informative-choice-of-parameters-for-a-dirichlet-distribution
# alpha = 1 or alpha = 0
#draws <- rdirichlet(200, c(1,1,1))
#bivt.contour(draws)

# use this computing rho probs and saving post parameter at each t
TS_MultiDir_full = function(N, n_arms = 2, y_categories = 1:7,
                            prior_params = list(arm0 = rep(1,7), arm1 = rep(1,7)), 
                            true_params = list(p0=MTurk_p0, p1=MTurk_p1),
                            compute_rho = T){
  
  params_dim = length(y_categories) 
  dirichlet_params = list(arm0 = matrix(NA, nrow=N, ncol=params_dim),
                          arm1 = matrix(NA, nrow=N, ncol=params_dim))
  
  dirichlet_params$arm0[1,] = prior_params$arm0
  dirichlet_params$arm1[1,] = prior_params$arm1
  
  p0_hat = matrix(NA, ncol=params_dim, nrow=N)
  p1_hat = matrix(NA, ncol=params_dim, nrow=N)
  y = c()
  a_hat = c()
  rho0 = c()
  for (i in 1:N) { 
    
    if(compute_rho==T){
      rho0[i] = rho0_Multi(params1 = dirichlet_params$arm0[i,], params2 = dirichlet_params$arm1[i,], 
                     y_categories = y_categories)
    }
    #sample from prior/posterior and save it
    p0_hat[i,] = rdirichlet(1, dirichlet_params$arm0[i,])
    p1_hat[i,] = rdirichlet(1, dirichlet_params$arm1[i,])
    
    #choose estimated best arm
    mean0 = sum(p0_hat[i,]*c(y_categories))
    mean1 = sum(p1_hat[i,]*c(y_categories))
    a_hat[i] = (mean1 > mean0)*1
    
    x = rmultinom(1, size = 1, prob = (true_params$p0*(a_hat[i] == 0) + true_params$p1*(a_hat[i] == 1)))
    
    if(i < N){
      dirichlet_params$arm0[i+1,] = dirichlet_params$arm0[i,] + x*(a_hat[i] == 0)
      dirichlet_params$arm1[i+1,] = dirichlet_params$arm1[i,] + x*(a_hat[i] == 1)
    }
    
    y_tmp = c(y_categories)*x
    y[i] = y_tmp[y_tmp!=0]
    #y0[i] = ifelse(a_hat[i] == 0, y[i], "NA")
    #y1[i] = ifelse(a_hat[i] == 1, y[i], "NA")
    
    
  }
  
  return(list(dirichlet_params = dirichlet_params, p0_hat = p0_hat, p1_hat = p1_hat, 
              #y = suppressWarnings(cbind(as.numeric(y0), as.numeric(y1))),
              y = suppressWarnings((as.numeric(y))),
              #y0 = suppressWarnings(as.numeric(y0)), 
              #y1=suppressWarnings(as.numeric(y1)), 
              a_hat = a_hat, all_prop = table(factor(a_hat, levels = 0:1))/N,
              rho0 = rho0))
}

# computationally & memory lighter version: use this for MTS with aug support
TS_MultiDir = function(N, n_arms = 2, y_categories = 1:7,
                       prior_params = list(arm0 = rep(1,7), arm1 = rep(1,7)), 
                       true_params = list(p0=MTurk_p0, p1=MTurk_p1),
                       compute_rho = F){
  
  params_dim = length(y_categories) 
  dirichlet_params = list(arm0 = prior_params$arm0,
                          arm1 = prior_params$arm1)
  
  dirichlet_params$arm0 = prior_params$arm0
  dirichlet_params$arm1 = prior_params$arm1
  
  #p0_hat = matrix(NA, ncol=params_dim, nrow=N)
  #p1_hat = matrix(NA, ncol=params_dim, nrow=N)
  y = c()
  a_hat = c()
  rho0 = c()
  for (i in 1:N) { 
    
    if(compute_rho==T){
      rho0[i] = rho0_Multi(params1 = dirichlet_params$arm0, params2 = dirichlet_params$arm1, 
                           y_categories = y_categories)
    }
    #sample from prior/posterior, but do not save it in memory
    p0_hat = rdirichlet(1, dirichlet_params$arm0)
    p1_hat = rdirichlet(1, dirichlet_params$arm1)
    
    #choose estimated best arm
    mean0 = sum(p0_hat*c(y_categories))
    mean1 = sum(p1_hat*c(y_categories))
    a_hat[i] = (mean1 > mean0)*1
    
    x = rmultinom(1, size = 1, prob = (true_params$p0*(a_hat[i] == 0) + true_params$p1*(a_hat[i] == 1)))
    
    if(i < N){
      dirichlet_params$arm0 = dirichlet_params$arm0 + x*(a_hat[i] == 0)
      dirichlet_params$arm1 = dirichlet_params$arm1 + x*(a_hat[i] == 1)
    }
    
    y_tmp = c(y_categories)*x
    y[i] = y_tmp[y_tmp!=0]
    
  }
  
  return(list(#p0_hat = p0_hat, p1_hat = p1_hat, 
              #y = suppressWarnings(cbind(as.numeric(y0), as.numeric(y1))),
              y = suppressWarnings((as.numeric(y))),
              #y0 = suppressWarnings(as.numeric(y0)), 
              #y1=suppressWarnings(as.numeric(y1)), 
              a_hat = a_hat, all_prop = table(factor(a_hat, levels = 0:1))/N,
              rho0 = rho0))
}

# use this with real Mturk data
TS_MultiDir_real = function(N=110, n_arms = 2, y_categories = 1:7,
                       prior_params = list(arm0 = rep(1,7), arm1 = rep(1,7)), 
                       true_data = list(arm1 = MTurk_data[MTurk_data$arm == "Arm 1",c(22)],
                                        arm2 = MTurk_data[MTurk_data$arm == "Arm 2",c(22)]),
                       compute_rho = F){
  
  params_dim = length(y_categories) 
  dirichlet_params = list(arm0 = prior_params$arm0,
                          arm1 = prior_params$arm1)
  
  dirichlet_params$arm0 = prior_params$arm0
  dirichlet_params$arm1 = prior_params$arm1
  
  #p0_hat = matrix(NA, ncol=params_dim, nrow=N)
  #p1_hat = matrix(NA, ncol=params_dim, nrow=N)
  y = c()
  a_hat = c()
  rho0 = c()
  for (i in 1:N) { 
    
    if(compute_rho==T){
      rho0[i] = rho0_Multi(params1 = dirichlet_params$arm0, params2 = dirichlet_params$arm1, 
                           y_categories = y_categories)
    }
    #sample from prior/posterior
    p0_hat = rdirichlet(1, dirichlet_params$arm0)
    p1_hat = rdirichlet(1, dirichlet_params$arm1)
    
    #choose estimated best arm
    mean0 = sum(p0_hat*c(y_categories))
    mean1 = sum(p1_hat*c(y_categories))
    a_hat[i] = (mean1 > mean0)*1
    
    #x = rmultinom(1, size = 1, prob = (true_params$p0*(a_hat[i] == 0) + true_params$p1*(a_hat[i] == 1)))
    idx = ifelse(a_hat[i] == 0, sample(na.exclude(true_data$arm1), 1), sample(na.exclude(true_data$arm2), 1))
    x = rep(0, params_dim)
    x[idx] = 1
    
    if(i < N){
      dirichlet_params$arm0 = dirichlet_params$arm0 + x*(a_hat[i] == 0)
      dirichlet_params$arm1 = dirichlet_params$arm1 + x*(a_hat[i] == 1)
    }
    
    y_tmp = c(y_categories)*x
    y[i] = y_tmp[y_tmp!=0]
    #y0[i] = ifelse(a_hat[i] == 0, y[i], "NA")
    #y1[i] = ifelse(a_hat[i] == 1, y[i], "NA")
    
    
  }
  
  return(list(#p0_hat = p0_hat, p1_hat = p1_hat, 
    #y = suppressWarnings(cbind(as.numeric(y0), as.numeric(y1))),
    y = suppressWarnings((as.numeric(y))),
    #y0 = suppressWarnings(as.numeric(y0)), 
    #y1=suppressWarnings(as.numeric(y1)), 
    a_hat = a_hat, all_prop = table(factor(a_hat, levels = 0:1))/N,
    rho0 = rho0))
}


# Normal TS ---------------------------------------------------

#TS_NormNorm (sigma known)
TS_NormNorm = function(N, n_arms = 2, normal_priors = list(mu=c(4,4), sigma2=c(100,100)), y_categories = 1:7,
                       true_params = list(p0=MTurk_p0, p1=MTurk_p1), true_vars = MTurk_vars,
                       compute_rho = F){
  
  params_dim = length(y_categories)
  
  normal_params = list(mu=matrix(NA, nrow=N, ncol=length(normal_priors)),
                       sigma2=matrix(NA, nrow=N, ncol=length(normal_priors)))
  normal_params$mu[1,] = normal_priors$mu
  normal_params$sigma2[1,] = normal_priors$sigma2
  
  mu_hat = matrix(NA, ncol=n_arms, nrow=N)
  #y = matrix(NA, ncol=n_arms, nrow=N)
  y = c()
  a_hat = rho0 = c()
  for (i in 1:N) { 
    
    if(compute_rho==T){
      rho0[i] = rho0_Norm(means = normal_params$mu[i,], sds = sqrt(normal_params$sigma2[i,])) 
    }
    
    mu_hat[i,] = rnorm(n_arms, mean = normal_params$mu[i,], sd = sqrt(normal_params$sigma2[i,]))
    a_hat[i] = (mu_hat[i,1] < mu_hat[i,2])*1
    
    x = rmultinom(1, size = 1, prob = (true_params$p0*(a_hat[i] == 0) + true_params$p1*(a_hat[i] == 1)))
    #y_tmp = sum(c(y_categories)*x)
    
    y[i] = sum(c(y_categories)*x)
    #y[i,] = ifelse(rep(a_hat[i] == 0, n_arms), c(y_tmp, NA), c(NA, y_tmp))
    n=c(sum(a_hat==0), sum(a_hat==1))
    #y_mean = apply(y, 2, mean, na.rm=T)
    y_mean = c(mean(y[a_hat==0]), mean(y[a_hat==1]))
    y_mean[is.na(y_mean)] = 0
    
    if(i < N){
      normal_params$sigma2[i+1,] = 1/(1/normal_priors$sigma2 + n/true_vars)
      normal_params$mu[i+1,] = normal_params$sigma2[i+1,]*(normal_priors$mu/normal_priors$sigma2+y_mean*n/(true_vars))
    }
    
  }
  
  return(list(normal_params = normal_params, mu_hat = mu_hat, y = y, 
              a_hat = a_hat, all_prop = table(factor(a_hat, levels = 0:1))/N, rho0 = rho0))
  
}


# Binary TS ---------------------------------------------------------------

#TS_BetaBern
TS_BetaBern = function(N, n_arms = 2, beta_priors = list(alpha=c(1,1), beta=c(1,1)), y_categories=1:7,
                       true_params = list(p0=MTurk_p0, p1=MTurk_p1), cutoff_value,
                       compute_rho = F){
  
  params_dim = length(y_categories)
  
  beta_params = list(alpha=matrix(NA, nrow=N, ncol=length(beta_priors)),
                     beta=matrix(NA, nrow=N, ncol=length(beta_priors)))
  beta_params$alpha[1,] = beta_priors$alpha
  beta_params$beta[1,] = beta_priors$beta
  
  p_hat = matrix(NA, ncol=n_arms, nrow=N)
  #y = y_mn = matrix(NA, ncol=n_arms, nrow=N)
  y = y_mn = c()
  a_hat = rho0 = c()
  for (i in 1:N) {
    
    if(compute_rho==T){
      rho0[i] = rho0_Binary(params1 = beta_params$alpha[i,], params2 = beta_params$beta[i,])
    }
    
    p_hat[i,] = rbeta(n_arms, beta_params$alpha[i,], beta_params$beta[i,])
    a_hat[i] = (p_hat[i,1] < p_hat[i,2])*1
    
    x = rmultinom(1, size = 1, prob = (true_params$p0*(a_hat[i] == 0) + true_params$p1*(a_hat[i] == 1)))
    #y_tmp = sum(c(y_categories)*x)
    y[i] = sum(c(y_categories)*x)
    
    #y_mn[i,] = ifelse(rep(a_hat[i] == 0, n_arms), c(y_tmp, NA), c(NA, y_tmp))
    
    y_binary = ifelse(y[i]<cutoff_value, 0, 1)
    y_binary = ifelse(rep(a_hat[i] == 0, n_arms), c(y_binary, NA), c(NA, y_binary))
    
    if(i < N){
      beta_params$alpha[i+1,] = colSums(rbind(beta_params$alpha[i,], y_binary), na.rm = T)
      beta_params$beta[i+1,] = colSums(rbind(beta_params$beta[i,], 1-y_binary), na.rm = T)
    }
    
  }
  
  return(list(beta_params = beta_params, p_hat = p_hat, y = y, a_hat = a_hat, 
              all_prop = table(factor(a_hat, levels = 0:1))/N, rho0 = rho0))
  
}


# Average Reward and Prop Optimal action -------------------

# For a single dataset
get_reward = function(my_data, to_t=myN, opt_arm = 1){
  # Note: this function is for the two-arm case (n_arms = 2) a normal rewards
  
  # Function inputs:
  # my_data:: db with the trial data containing:
  # variable "Arm": selected arm 
  # to_t:: number of round t for computing reward
  # opt_arm:: arm that is optimal (1 or 2)
  
  # Function output: a 4-dim vector with
  # N1, N2:: number of times arm 1 and arm 2 were selected
  # reward:: cumulative reward
  # prop_opt_arm:: proportion of times the optimal arm was selected
  
  #cum_reward = sum(my_data$y[1:to_t,], na.rm = T)
  cum_reward = sum(my_data$y[1:to_t], na.rm = T)
  
  N1 = sum(my_data$a_hat[1:to_t] == 0)
  N2 = sum(my_data$a_hat[1:to_t] == 1)
  
  if(opt_arm == 1){
    prop_opt_arm = N1/(N1+N2)
  } else {
    prop_opt_arm = N2/(N1+N2)
  }
  
  return(c(N1 = N1, N2 = N2, cum_reward = cum_reward, prop_opt_arm = prop_opt_arm))
}

# Average (across simulations)
mean_rew_optarm = function(dataMN = TSmn_H1_mturk, dataN = TSn_H1_mturk, dataB = TSb_H1_mturk, dataB1 = TSb1_H1_mturk, rounds = c(10,20,30,40,seq(50,myN,50)), mu = MTurk_means, M = Nsim){
  # Note: this function is for the two-arm case (n_arms = 2) a normal rewards
  
  # Function inputs:
  # data_Pitest, data_BOLS:: db with the trial data by batch simulated according to TS for Pi-test and TS for BOLS based Z-test, containing:
  # variable "Batch": the batch (i.e., time or update) number
  # variable "Arm": selected arm 
  # variable "AlloProb1": the allocation probabilities of arm 1 in each batch
  # batches: batches for which reward performances have to be computed
  
  # Function output: a 2-dim list with Avg_Reward & Avg_Prop_opt_arm, each being a dataframe with columns
  # Batch (or time) number
  # Avg_Reward or Avg_Prop_opt_arm for standard TS (Pi-test)
  # Avg_Reward or Avg_Prop_opt_arm for tuned TS (BOLS-test)
  Avg_Reward = c()
  Avg_Reward_SE = c()
  Avg_Regret = c()
  Avg_Prop_opt_arm = c()
  for(k in rounds){
    print(paste0("Computing Reward Data for Round ", k))
    Rew_TSmn = pblapply(dataMN, get_reward, to_t = k)
    Rew_TSmn = data.frame(matrix(unlist(Rew_TSmn), nrow=length(Rew_TSmn), byrow=TRUE))
    
    Rew_TSn = pblapply(dataN, get_reward, to_t = k)
    Rew_TSn = data.frame(matrix(unlist(Rew_TSn), nrow=length(Rew_TSn), byrow=TRUE))
    
    Rew_TSb = pblapply(dataB, get_reward, to_t = k)
    Rew_TSb = data.frame(matrix(unlist(Rew_TSb), nrow=length(Rew_TSb), byrow=TRUE))
    
    Rew_TSb1 = pblapply(dataB1, get_reward, to_t = k)
    Rew_TSb1 = data.frame(matrix(unlist(Rew_TSb1), nrow=length(Rew_TSb1), byrow=TRUE))
    
    names(Rew_TSmn) = names(Rew_TSb) = names(Rew_TSb1) = names(Rew_TSn) = c("N1", "N2", "Reward", "Arm_Alloc")
    
    Rew_TSmn_mean = apply(Rew_TSmn, 2, mean) #get avg (across sim)
    Rew_TSmn_se = sd(Rew_TSmn$Reward)/sqrt(M)
    Rew_TSn_mean = apply(Rew_TSn, 2, mean) #get avg (across sim)
    Rew_TSn_se = sd(Rew_TSn$Reward)/sqrt(M)
    Rew_TSb_mean = apply(Rew_TSb, 2, mean) #get avg (across sim)
    Rew_TSb_se = sd(Rew_TSb$Reward)/sqrt(M)
    Rew_TSb1_mean = apply(Rew_TSb1, 2, mean) #get avg (across sim)
    Rew_TSb1_se = sd(Rew_TSb1$Reward)/sqrt(M)
    
    Avg_Reward = rbind(Avg_Reward, c(Round = k, MultiN = Rew_TSmn_mean[3],
                                     Normal = Rew_TSn_mean[3],
                                     Binary = Rew_TSb_mean[3],
                                     Binary1 = Rew_TSb1_mean[3]))
    
    Avg_Regret = rbind(Avg_Regret, c(Round = k, MultiN = k*mu[1] - (Rew_TSmn_mean[1]*mu[1] + Rew_TSmn_mean[2]*mu[2]),
                                     Normal = k*mu[1] - (Rew_TSn_mean[1]*mu[1] + Rew_TSn_mean[2]*mu[2]),
                                     Binary = k*mu[1] - (Rew_TSb_mean[1]*mu[1] + Rew_TSb_mean[2]*mu[2]),
                                     Binary1 = k*mu[1] - (Rew_TSb1_mean[1]*mu[1] + Rew_TSb1_mean[2]*mu[2])))
    #Reward SE = Regret SE
    Avg_Reward_SE = rbind(Avg_Reward_SE, c(Round = k, MultiN = Rew_TSmn_se,
                                           Normal = Rew_TSn_se,
                                           Binary = Rew_TSb_se,
                                           Binary1 = Rew_TSb1_se))
    
    Avg_Prop_opt_arm = rbind(Avg_Prop_opt_arm, c(Round = k, MultiN = Rew_TSmn_mean[4],
                                                 Normal = Rew_TSn_mean[4],
                                                 Binary = Rew_TSb_mean[4],
                                                 Binary1 = Rew_TSb1_mean[4]))
    
  }
  
  return(list(Avg_Reward = Avg_Reward, Avg_Reward_SE = Avg_Reward_SE, Avg_Regret = Avg_Regret, Avg_Prop_opt_arm = Avg_Prop_opt_arm))
}

