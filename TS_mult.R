# Load functions and MTurk data ---------------------------------------------------
library(Compositional)
library(MCMCpack)
library(pbapply)
library(ggplot2)
library(dplyr)
library(ggridges)

source("TS_functions_SMA.R")

# H1 based on MTurk data ---------------------------------------------------
load("/Users/blackmamba/Desktop/TS_multinomial/MTurk Data/MTurk_data.RData")

# |-- H1. Multinomial TS ---------------------------------------------------
Nsim = 1e4; myN=1000

set.seed(12345)
start_time <- Sys.time()
TSmn_H1_mturk = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                   prior_params = list(arm0 = rep(3,7), arm1 = rep(3,7)), compute_rho = F)))
end_time <- Sys.time()
time_TSmnH1_mturk = end_time - start_time

#looking at average values
TSmn_H1_alloc_mturk <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
hist(TSmn_H1_alloc_mturk[,1], xlim=c(0,1), breaks = 150, freq = F, xlab = "Emprirical Allocation (Multinomial)", main = "Empirical alloc (H1: Mturk 5.81 vs 5.08)")
mean(TSmn_H1_alloc_mturk[,1]); mean(TSmn_H1_alloc_mturk[,2])
quantile(TSmn_H1_alloc_mturk[,1], probs = c(0.05, 0.1, 0.9, 0.95))
mtext(paste("Empirical Alloc Best Arm = ", round(mean(TSmn_H1_alloc_mturk[,1]), 3)))

TSmn_H1_y0_mturk <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H1_y0_mturk[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical distribution of Y (H1: Mturk 5.81 vs 5.08)")
abline(v=mean(TSmn_H1_y0_mturk[,1]), col="red")
mtext(paste("bias = ", round(MTurk_means[1] - mean(TSmn_H1_y0_mturk[,1]), 3)))

TSmn_H1_y1_mturk <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H1_y1_mturk[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical distribution of Y (H1: Mturk 5.81 vs 5.08)")
abline(v=mean(TSmn_H1_y1_mturk[,1], na.rm=T), col="red")
mtext(paste("bias = ", round(MTurk_means[2] - mean(TSmn_H1_y1_mturk[,1], na.rm=T), 3)))

#looking at a single trajectory
table(TSmn_H1_mturk[[10]]$a_hat)
tail(TSmn_H1_mturk[[10]]$p0_hat); tail(TSmn_H1_mturk[[10]]$p1_hat)
plot(cumsum(TSmn_H1_mturk[[10]]$a_hat==0)/1:myN, type="l", ylim=c(0,1), col = "red")
lines(cumsum(TSmn_H1_mturk[[10]]$a_hat==1)/1:myN, col="blue")
#mean(TSmn_H1_mturk[[10]]$y[x$a_hat==1], na.rm = T); mean(TSmn_H1_mturk[[10]]$y[x$a_hat==1], na.rm = T)


# |-- H1. Normal TS under H1 ---------------------------------------------------

library(gtools)
comb = combinations(n = 7, r = 2, v = 1:7, repeats.allowed = F)

mycov0 = mycov1 = c()
for(i in 1:dim(comb)[1]){
  mycov0[i] = prod(comb[i,])*prod(MTurk_p0[comb[i,]])
  mycov1[i] = prod(comb[i,])*prod(MTurk_p1[comb[i,]])
}

var0 = sum(c(1:7)^2*MTurk_p0*(1-MTurk_p0)) - 2*sum(mycov0)
var1 = sum(c(1:7)^2*MTurk_p1*(1-MTurk_p1)) - 2*sum(mycov1)

set.seed(12345)
start_time <- Sys.time()
TSn_H1_mturk = pbreplicate(Nsim, list(TS_NormNorm(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), true_vars = c(var0, var1),
                                                  normal_priors = list(mu=c(4,4), sigma2=c(10,10)))))
end_time <- Sys.time()
time_TSnH1_mturk = end_time - start_time

#looking at average values
TSn_H1_alloc <- data.frame(matrix(unlist(lapply(TSn_H1_mturk, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
hist(TSn_H1_alloc[,1], breaks = 150, xlim=c(0,1), freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H1: Mturk 5.81 vs 5.08)")
mean(TSn_H1_alloc[,1]); mean(TSn_H1_alloc[,2])
quantile(TSn_H1_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
mtext(paste("Empirical Alloc Best Arm = ", round(mean(TSn_H1_alloc[,1]), 3)))

TSn_H1_y0_mturk <- data.frame(matrix(unlist(lapply(TSn_H1_mturk, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSn_H1_y0_mturk[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical distribution of Y (H1: Mturk 5.81 vs 5.08)")
abline(v=mean(TSn_H1_y0_mturk[,1]), col="red")
mtext(paste("bias = ", round(MTurk_means[1] - mean(TSn_H1_y0_mturk[,1]), 3)))

TSn_H1_y1_mturk <- data.frame(matrix(unlist(lapply(TSn_H1_mturk, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSn_H1_y1_mturk[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical distribution of Y (H1: Mturk 5.81 vs 5.08)")
abline(v=mean(TSn_H1_y1_mturk[,1], na.rm=T), col="red")
mtext(paste("bias = ", round(MTurk_means[2] - mean(TSn_H1_y1_mturk[,1], na.rm=T), 3)))

#looking at a single trajectory
table(TSn_H1_mturk[[1]]$a_hat)
tail(TSn_H1_mturk[[1]]$mu_hat)
plot(cumsum(TSn_H1_mturk[[1]]$a_hat==0)/1:myN, type="l", ylim=c(0,1), col = "red")
lines(cumsum(TSn_H1_mturk[[1]]$a_hat==1)/1:myN, col="blue")
mean(TSn_H1_mturk[[1]]$y[,1], na.rm = T); mean(TSn_H1_mturk[[1]]$y[,2], na.rm = T)

# |-- H1. Run Sim: Binary TS under H1 (cutoff 6) ---------------------------------------------------

# my_true_p = c(mean(MTurk_data$reward[MTurk_data$arm=="Arm 1"]==7, na.rm = T),
#               mean(MTurk_data$reward[MTurk_data$arm=="Arm 2"]==7, na.rm = T))

set.seed(12345)
start_time <- Sys.time()
TSb_H1_mturk = pbreplicate(Nsim, list(TS_BetaBern(N=myN, n_arms = 2, y_categories = 1:7,
                                                  beta_priors = list(alpha=c(10,10), beta=c(10,10)), 
                                                  true_params = list(p0 = MTurk_p0, p1 = MTurk_p1),
                                                  cutoff_value = 6)))
end_time <- Sys.time()
time_TSbH1_mturk = end_time - start_time

#looking at average values
TSb_H1_alloc <- data.frame(matrix(unlist(lapply(TSb_H1_mturk, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
hist(TSb_H1_alloc[,1], breaks = 150, xlim=c(0,1), freq = F, xlab = "Emprirical Allocation (Binary model)", main = "Empirical allocation (H1: Mturk 5.81 vs 5.08)")
mean(TSb_H1_alloc[,1]); mean(TSb_H1_alloc[,2])
quantile(TSb_H1_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
mtext(paste("Empirical Alloc Best Arm = ", round(mean(TSb_H1_alloc[,1]), 3)))

TSb_H1_y0_mturk <- data.frame(matrix(unlist(lapply(TSb_H1_mturk, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSb_H1_y0_mturk[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical distribution of Y (H1: Mturk 0.24 vs 0.19)")
abline(v=mean(TSb_H1_y0_mturk[,1]), col="red")
mtext(paste("bias = ", round(MTurk_means[1] - mean(TSb_H1_y0_mturk[,1]), 3)))

TSb_H1_y1_mturk <- data.frame(matrix(unlist(lapply(TSb_H1_mturk, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSb_H1_y1_mturk[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical distribution of Y (H1: Mturk 0.24 vs 0.19)")
abline(v=mean(TSb_H1_y1_mturk[,1], na.rm=T), col="red")
mtext(paste("bias = ", round(MTurk_means[2] - mean(TSb_H1_y1_mturk[,1], na.rm=T), 3)))

#looking at a single trajectory
table(TSb_H1_mturk[[1]]$a_hat)
tail(TSb_H1_mturk[[1]]$p_hat)
plot(cumsum(TSb_H1_mturk[[1]]$a_hat==0)/1:myN, type="l", ylim=c(0,1), col = "red")
lines(cumsum(TSb_H1_mturk[[1]]$a_hat==1)/1:myN, col="blue")
mean(TSb_H1_mturk[[1]]$y[,1], na.rm = T); mean(TSb_H1_mturk[[1]]$y[,2], na.rm = T)

# |-- H1. Run Sim: Binary TS under H1 (cutoff 7) ---------------------------------------------------


set.seed(12345)
start_time <- Sys.time()
TSb1_H1_mturk = pbreplicate(Nsim, list(TS_BetaBern(N=myN, n_arms = 2, 
                                                  beta_priors = list(alpha=c(10,10), beta=c(10,10)), 
                                                  true_params = list(p0 = MTurk_p0, p1 = MTurk_p1),
                                                  cutoff_value = 7)))
end_time <- Sys.time()
time_TSbH1_mturk = end_time - start_time

#looking at average values
TSb1_H1_alloc <- data.frame(matrix(unlist(lapply(TSb1_H1_mturk, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
hist(TSb1_H1_alloc[,1], breaks = 150, xlim=c(0,1), freq = F, xlab = "Emprirical Allocation (Binary model)", main = "Empirical allocation (H1: Mturk 5.81 vs 5.08)")
mean(TSb1_H1_alloc[,1]); mean(TSb1_H1_alloc[,2])
quantile(TSb1_H1_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
mtext(paste("Empirical Alloc Best Arm = ", round(mean(TSb1_H1_alloc[,1]), 3)))

TSb1_H1_y0_mturk <- data.frame(matrix(unlist(lapply(TSb1_H1_mturk, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSb1_H1_y0_mturk[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical distribution of Y (H1: Mturk 0.24 vs 0.19)")
abline(v=mean(TSb1_H1_y0_mturk[,1]), col="red")
mtext(paste("bias = ", round(MTurk_means[1] - mean(TSb1_H1_y0_mturk[,1]), 3)))

TSb1_H1_y1_mturk <- data.frame(matrix(unlist(lapply(TSb1_H1_mturk, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSb1_H1_y1_mturk[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical distribution of Y (H1: Mturk 0.24 vs 0.19)")
abline(v=mean(TSb1_H1_y1_mturk[,1], na.rm=T), col="red")
mtext(paste("bias = ", round(MTurk_means[2] - mean(TSb1_H1_y1_mturk[,1], na.rm=T), 3)))

#looking at a single trajectory
table(TSb1_H1_mturk[[1]]$a_hat)
tail(TSb1_H1_mturk[[1]]$p_hat)
plot(cumsum(TSb1_H1_mturk[[1]]$a_hat==0)/1:myN, type="l", ylim=c(0,1), col = "red")
lines(cumsum(TSb1_H1_mturk[[1]]$a_hat==1)/1:myN, col="blue")
mean(TSb1_H1_mturk[[1]]$y[,1], na.rm = T); mean(TSb1_H1_mturk[[1]]$y[,2], na.rm = T)

# |-- H1. Computing performance data: Reward and Prop Optimal action -------------------

Rew_data = mean_rew_optarm()

for(i in 1:4){
  Rew_data[[i]] = rbind(rep(0, 5), Rew_data[[i]])
}

Rew_data$Avg_Prop_opt_arm[1,] = c(0, 0.5, 0.5, 0.5, 0.5)

save(Rew_data, file = "Rewards1e4_SMA2.RData")


# H0 with symmetric distribution --------------------------------------------------- 

# |-- H0. Multinomial TS ---------------------------------------------------
Nsim = 1e4; myN=1000

# symmetric (~normal) data
bin_norm = rbinom(1e5, 6, 0.5) + 1
plot(table(bin_norm), type = "h", lwd = 2, main = "True distr (H0: p0=p1; sym mu=4)", xlab="Y (multinom)", ylab = "")
mean(bin_norm)

true_p0_sym = round(table(bin_norm)/length(bin_norm),2)
mean_p0_sym = sum(true_p0_sym*1:7)

set.seed(12345)
start_time <- Sys.time()
TSmn_H0sym = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), compute_rho = F)))
end_time <- Sys.time()
time_TSmnH0_mturk = end_time - start_time

#looking at average values
TSmn_H0sym_alloc <- data.frame(matrix(unlist(lapply(TSmn_H0sym, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSmn_H0sym_alloc[,1], breaks = 150, freq = F, xlab = "Ordered Multinomial", ylab = "Emprirical Allocation Arm 0", main = "Empirical alloc (H0: p0=p1 sym mu=4)")
mean(TSmn_H0sym_alloc[,1]); mean(TSmn_H0sym_alloc[,2])
TSmn_H0sym_quant = quantile(TSmn_H0sym_alloc[,1], probs = 0.5)
abline(v = TSmn_H0sym_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSmn_H0sym_quant[1], at = TSmn_H0sym_quant[1], side = 1, col = "red")
abline(v = TSmn_H0sym_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSmn_H0sym_quant[4], at = TSmn_H0sym_quant[4], side = 1, col = "red")

c(mean(TSmn_H0sym_alloc[,1]<0.05)*100, 
  mean(TSmn_H0sym_alloc[,1]<0.1)*100, 
  mean(TSmn_H0sym_alloc[,1]<=0.55 & TSmn_H0sym_alloc[,1]>=0.45)*100, 
  mean(TSmn_H0sym_alloc[,1]<=0.6 & TSmn_H0sym_alloc[,1]>=0.4)*100, 
  mean(TSmn_H0sym_alloc[,1]>0.9)*100, 
  mean(TSmn_H0sym_alloc[,1]>0.95)*100)

TSmn_H0sym_y0 <- data.frame(matrix(unlist(lapply(TSmn_H0sym, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H0sym_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical dist (H0: p0=p1 sym mu=4)")
abline(v=mean(TSmn_H0sym_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_sym - mean(TSmn_H0sym_y0[,1]), 3)))

TSmn_H0sym_y1 <- data.frame(matrix(unlist(lapply(TSmn_H0sym, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H0sym_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical dist (H0: p0=p1 sym mu=4)")
abline(v=mean(TSmn_H0sym_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_sym - mean(TSmn_H0sym_y1[,1]), 3)))

# |-- H0. Normal TS ---------------------------------------------------

library(gtools)
comb = combinations(n = 7, r = 2, v = 1:7, repeats.allowed = F)

mycov_H0 = c()
for(i in 1:dim(comb)[1]){
  mycov_H0[i] = prod(comb[i,])*prod(true_p0_sym[comb[i,]])
}

var_H0 = sum(c(1:7)^2*true_p0_sym*(1-true_p0_sym)) - 2*sum(mycov_H0)

set.seed(12345)
start_time <- Sys.time()
TSn_H0sym = pbreplicate(Nsim, list(TS_NormNorm(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                               true_vars = c(var_H0, var_H0), compute_rho = F)))
end_time <- Sys.time()
time_TSn_H0sym = end_time - start_time


#looking at average values
TSn_H0sym_alloc <- data.frame(matrix(unlist(lapply(TSn_H0sym, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSn_H0sym_alloc[,1], breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")
mean(TSn_H0sym_alloc[,1]); mean(TSn_H0sym_alloc[,2])
TSn_H0sym_quant = quantile(TSn_H0sym_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
abline(v = TSn_H0sym_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSn_H0sym_quant[1], at = TSn_H0sym_quant[1], side = 1, col = "red")
abline(v = TSn_H0sym_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSn_H0sym_quant[4], at = TSn_H0sym_quant[4], side = 1, col = "red")

#check allocation compared to MultinomialTS
hist(TSn_H0sym_alloc[,1], col = "blue", breaks = 150, freq = F, xlab = " ", ylab = "", main = "Empirical allocation of Arm 1")
hist(TSmn_H0sym_alloc[,1],col = "red", add = T, breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")


TSn_H0sym_y0 <- data.frame(matrix(unlist(lapply(TSn_H0sym, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSn_H0sym_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical distribution of Y (H0: mu0=mu1=4)")
abline(v=mean(TSn_H0sym_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_sym - mean(TSn_H0sym_y0[,1]), 4)))

TSn_H0sym_y1 <- data.frame(matrix(unlist(lapply(TSn_H0sym, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSn_H0sym_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical distribution of Y (H0: mu0=mu1=4)")
abline(v=mean(TSn_H0sym_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_sym - mean(TSn_H0sym_y1[,1]), 4)))


# |-- H0. Binary TS (cutoff 6) ---------------------------------------------------

set.seed(12345)
start_time <- Sys.time()
TSb_H0_sym = pbreplicate(Nsim, list(TS_BetaBern(N=myN, n_arms = 2, y_categories = 1:7,
                                            beta_priors = list(alpha=c(1,1), beta=c(1,1)), 
                                            true_params = list(p0 = true_p0_sym, p1 = true_p0_sym),
                                            cutoff_value = 6, compute_rho = F)))
end_time <- Sys.time()
time_TSb_H0_sym = end_time - start_time

#looking at average values
TSb_H0sym_alloc <- data.frame(matrix(unlist(lapply(TSb_H0_sym, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSb_H0sym_alloc[,1], breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")
mean(TSb_H0sym_alloc[,1]); mean(TSb_H0sym_alloc[,2])
TSb_H0sym_quant = quantile(TSb_H0sym_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
abline(v = TSb_H0sym_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSb_H0sym_quant[1], at = TSb_H0sym_quant[1], side = 1, col = "red")
abline(v = TSb_H0sym_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSb_H0sym_quant[4], at = TSb_H0sym_quant[4], side = 1, col = "red")

#check allocation compared to MultinomialTS
hist(TSb_H0sym_alloc[,1], col = "blue", breaks = 150, freq = F, xlab = " ", ylab = "", main = "Empirical allocation of Arm 1")
hist(TSmn_H0sym_alloc[,1], col = "red", add = T, breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")
hist(TSn_H0sym_alloc[,1], col = "black", add = T, breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")


TSb_H0sym_y0 <- data.frame(matrix(unlist(lapply(TSb_H0_sym, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSb_H0sym_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical distribution of Y (H0: mu0=mu1=4)")
abline(v=mean(TSb_H0sym_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_sym - mean(TSb_H0sym_y0[,1]), 4)))

TSb_H0sym_y1 <- data.frame(matrix(unlist(lapply(TSb_H0_sym, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSb_H0sym_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical distribution of Y (H0: mu0=mu1=4)")
abline(v=mean(TSb_H0sym_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_sym - mean(TSb_H0sym_y1[,1]), 4)))

# |-- H0. Binary TS (cutoff 7) ---------------------------------------------------

set.seed(12345)
start_time <- Sys.time()
TSb1_H0_sym = pbreplicate(Nsim, list(TS_BetaBern(N=myN, n_arms = 2, y_categories = 1:7,
                                                beta_priors = list(alpha=c(1,1), beta=c(1,1)), 
                                                true_params = list(p0 = true_p0_sym, p1 = true_p0_sym),
                                                cutoff_value = 7, compute_rho = F)))
end_time <- Sys.time()
time_TSb1_H0_sym = end_time - start_time

#looking at average values
TSb1_H0sym_alloc <- data.frame(matrix(unlist(lapply(TSb1_H0_sym, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSb1_H0sym_alloc[,1], breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")
mean(TSb1_H0sym_alloc[,1]); mean(TSb1_H0sym_alloc[,2])
TSb1_H0sym_quant = quantile(TSb1_H0sym_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
abline(v = TSb1_H0sym_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSb1_H0sym_quant[1], at = TSb1_H0sym_quant[1], side = 1, col = "red")
abline(v = TSb1_H0sym_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSb1_H0sym_quant[4], at = TSb1_H0sym_quant[4], side = 1, col = "red")

#check allocation compared to MultinomialTS
hist(TSn_H0sym_alloc[,1], xlim=c(0,1), ylim=c(0,2), col = "black", breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")
hist(TSmn_H0sym_alloc[,1], col = "red", add = T, breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")
hist(TSb_H0sym_alloc[,1], col = "green", add = T, breaks = 150, freq = F, xlab = " ", ylab = "", main = "Empirical allocation of Arm 1")
hist(TSb1_H0sym_alloc[,1], xlim=c(0,1), add = T, col = "blue", breaks = 150, freq = F, xlab = " ", ylab = "", main = "Empirical allocation of Arm 1")


TSb1_H0sym_y0 <- data.frame(matrix(unlist(lapply(TSb1_H0_sym, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSb1_H0sym_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical distribution of Y (H0: mu0=mu1=4)")
abline(v=mean(TSb1_H0sym_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_sym - mean(TSb1_H0sym_y0[,1]), 4)))

TSb1_H0sym_y1 <- data.frame(matrix(unlist(lapply(TSb1_H0_sym, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSb1_H0sym_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical distribution of Y (H0: mu0=mu1=4)")
abline(v=mean(TSb1_H0sym_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_sym - mean(TSb1_H0sym_y1[,1]), 4)))


# |-- H0: Compute and Plot performance data: Prop Optimal action -------------------

# Build dataset with different distributionsc("Multinomial TS", "Normal TS")
data_alloc_H0sym <- data.frame(
  type = c( rep("Multinomial TS", Nsim), rep("Normal TS", Nsim), 
            rep("Binary TS (cutoff = 6)", Nsim), rep("Binary TS (cutoff = 7)", Nsim), rep("Oracle",160) ),
  value = c(TSmn_H0sym_alloc[,1], TSn_H0sym_alloc[,1], 
            TSb_H0sym_alloc[,1], TSb1_H0sym_alloc[,1], rep(0.5, 160))
)

# Represent it
p_mn <- data_alloc_H0sym[data_alloc_H0sym$type=="Multinomial TS",] %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram(bins = 100, color="#F8766D", alpha=0.7, position = 'identity') +
  scale_fill_manual(values=c("#F8766D")) +
  #ggtitle("Empirical allocation of Arm 1") +
  #scale_alpha_manual(name = "category", values = c(0.9, 0.2, 1)) +
  theme_ipsum() +
  #theme(plot.title = element_text(size=15)) +
  labs(fill="") + ylab(" ") +  ylim(0,200) + xlim(0,1) +
  xlab(expression(T[1]*"/"*T)) +
  theme(legend.position=c(0.5, 0.95),
        legend.text=element_text(size=10)) +
  theme(axis.title.x = element_text(hjust = 0.5, face = "bold")) +
  theme(plot.margin=unit(c(0.2,0,0,0), "cm"))

p_n <- data_alloc_H0sym[data_alloc_H0sym$type=="Normal TS",] %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram(bins = 100, color="#D16DBE", alpha=0.7, position = 'identity') +
  scale_fill_manual(values=c("#D16DBE")) +
  #ggtitle("Empirical allocation of Arm 1") +
  #scale_alpha_manual(name = "category", values = c(0.9, 0.2, 1)) +
  theme_ipsum() +
  #theme(plot.title = element_text(size=15)) +
  labs(fill="") + ylab(" ") +  ylim(0,200) + xlim(0,1) +
  xlab(expression(T[1]*"/"*T)) +
  theme(legend.position=c(0.5, 0.95),
        legend.text=element_text(size=10)) +
  theme(axis.title.x = element_text(hjust = 0.5, face = "bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(plot.margin=unit(c(0.2,0,0,0), "cm"))

p_b <- data_alloc_H0sym[data_alloc_H0sym$type=="Binary TS (cutoff = 6)",] %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram(bins = 100, color="#9590FF", alpha=0.7, position = 'identity') +
  scale_fill_manual(values=c("#9590FF")) +
  #ggtitle("Empirical allocation of Arm 1") +
  #scale_alpha_manual(name = "category", values = c(0.9, 0.2, 1)) +
  theme_ipsum() +
  #theme(plot.title = element_text(size=15)) +
  labs(fill="") + ylab(" ") +  ylim(0,200) + xlim(0,1) +
  xlab(expression(T[1]*"/"*T)) +
  theme(legend.position=c(0.6, 0.95),
        legend.text=element_text(size=10)) +
  theme(axis.title.x = element_text(hjust = 0.5, face = "bold")) +
  theme(plot.margin=unit(c(0.2,0,0,0), "cm"))

p_b1 <- data_alloc_H0sym[data_alloc_H0sym$type=="Binary TS (cutoff = 7)",] %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram(bins = 100, color="#529EFF", alpha=0.7, position = 'identity') +
  scale_fill_manual(values=c("#529EFF")) +
  #ggtitle("Empirical allocation of Arm 1") +
  #scale_alpha_manual(name = "category", values = c(0.9, 0.2, 1)) +
  theme_ipsum() +
  #theme(plot.title = element_text(size=15)) +
  labs(fill="") + ylab(" ") +  ylim(0,200) + xlim(0,1) +
  xlab(expression(T[1]*"/"*T)) +
  theme(legend.position=c(0.6, 0.95),
        legend.text=element_text(size=10)) +
  theme(axis.title.x = element_text(hjust = 0.5, face = "bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(plot.margin=unit(c(0.2,0,0,0), "cm"))

ggarrange(p_mn, p_n, p_b, p_b1, ncol=2, nrow=2)

# H0. Right / Positively skewed distribution ---------------------------------------------------

true_p0_rs = c(0.2, 0.3, 0.15, 0.15, 0.1, 0.05, 0.05)
plot(true_p0_rs, type = "h", lwd = 2, main = "True distr (H0: p0=p1; skew mu=3)", xlab="Y (multinomial)", ylab = "")
mean_p0_rs = sum(true_p0_rs*1:7)

# |-- H0. Multinomial TS ---------------------------------------------------


set.seed(12345)
TSmn_H0rs = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_rs, p1 = true_p0_rs), compute_rho = F)))

#looking at average values
TSmn_H0rs_alloc <- data.frame(matrix(unlist(lapply(TSmn_H0rs, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSmn_H0rs_alloc[,1], breaks = 150, freq = F, xlab = "Ordered Multinomial", ylab = "Emprirical Allocation Arm 0", main = "Empirical alloc (H0: p0=p1 rs mu=4)")
mean(TSmn_H0rs_alloc[,1]); mean(TSmn_H0rs_alloc[,2])
TSmn_H0rs_quant = quantile(TSmn_H0rs_alloc[,1], probs = c(0.05, 0.1, 0.4, 0.6, 0.9, 0.95))
abline(v = TSmn_H0rs_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSmn_H0rs_quant[1], at = TSmn_H0rs_quant[1], side = 1, col = "red")
abline(v = TSmn_H0rs_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSmn_H0rs_quant[4], at = TSmn_H0rs_quant[4], side = 1, col = "red")

c(mean(TSmn_H0rs_alloc[,1]<0.05)*100, 
  mean(TSmn_H0rs_alloc[,1]<0.1)*100, 
  mean(TSmn_H0rs_alloc[,1]<=0.55 & TSmn_H0rs_alloc[,1]>=0.45)*100, 
  mean(TSmn_H0rs_alloc[,1]<=0.6 & TSmn_H0rs_alloc[,1]>=0.4)*100, 
  mean(TSmn_H0rs_alloc[,1]>0.9)*100, 
  mean(TSmn_H0rs_alloc[,1]>0.95)*100)

TSmn_H0rs_y0 <- data.frame(matrix(unlist(lapply(TSmn_H0rs, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H0rs_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical dist (H0: p0=p1 rs mu=4)")
abline(v=mean(TSmn_H0rs_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_rs - mean(TSmn_H0rs_y0[,1]), 3)))

TSmn_H0rs_y1 <- data.frame(matrix(unlist(lapply(TSmn_H0rs, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H0rs_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical dist (H0: p0=p1 rs mu=4)")
abline(v=mean(TSmn_H0rs_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_rs - mean(TSmn_H0rs_y1[,1]), 3)))

# |-- H0. Normal TS ---------------------------------------------------

mycov_H0 = c()
for(i in 1:dim(comb)[1]){
  mycov_H0[i] = prod(comb[i,])*prod(true_p0_rs[comb[i,]])
}
var_H0 = sum(c(1:7)^2*true_p0_rs*(1-true_p0_rs)) - 2*sum(mycov_H0)

set.seed(12345)
TSn_H0rs = pbreplicate(Nsim, list(TS_NormNorm(N=myN, true_params = list(p0 = true_p0_rs, p1 = true_p0_rs), 
                                               true_vars = c(var_H0, var_H0), compute_rho = F)))


#looking at average values
TSn_H0rs_alloc <- data.frame(matrix(unlist(lapply(TSn_H0rs, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSn_H0rs_alloc[,1], breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")
mean(TSn_H0rs_alloc[,1]); mean(TSn_H0rs_alloc[,2])
TSn_H0rs_quant = quantile(TSn_H0rs_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
abline(v = TSn_H0rs_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSn_H0rs_quant[1], at = TSn_H0rs_quant[1], side = 1, col = "red")
abline(v = TSn_H0rs_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSn_H0rs_quant[4], at = TSn_H0rs_quant[4], side = 1, col = "red")

#check allocation compared to MultinomialTS
hist(TSn_H0rs_alloc[,1], col = "blue", breaks = 150, freq = F, xlab = " ", ylab = "", main = "Empirical allocation of Arm 1")
hist(TSmn_H0rs_alloc[,1],col = "red", add = T, breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")


TSn_H0rs_y0 <- data.frame(matrix(unlist(lapply(TSn_H0rs, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSn_H0rs_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical distribution of Y (H0: mu0=mu1=4)")
abline(v=mean(TSn_H0rs_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_rs - mean(TSn_H0rs_y0[,1]), 4)))

TSn_H0rs_y1 <- data.frame(matrix(unlist(lapply(TSn_H0rs, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSn_H0rs_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical distribution of Y (H0: mu0=mu1=4)")
abline(v=mean(TSn_H0rs_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_rs - mean(TSn_H0rs_y1[,1]), 4)))


# |-- H0. Binary TS (cutoff 6) ---------------------------------------------------

set.seed(12345)
TSb_H0_rs = pbreplicate(Nsim, list(TS_BetaBern(N=myN, n_arms = 2, y_categories = 1:7,
                                                beta_priors = list(alpha=c(1,1), beta=c(1,1)), 
                                                true_params = list(p0 = true_p0_rs, p1 = true_p0_rs),
                                                cutoff_value = 6, compute_rho = F)))


#looking at average values
TSb_H0rs_alloc <- data.frame(matrix(unlist(lapply(TSb_H0_rs, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSb_H0rs_alloc[,1], breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")
mean(TSb_H0rs_alloc[,1]); mean(TSb_H0rs_alloc[,2])
TSb_H0rs_quant = quantile(TSb_H0rs_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
abline(v = TSb_H0rs_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSb_H0rs_quant[1], at = TSb_H0rs_quant[1], side = 1, col = "red")
abline(v = TSb_H0rs_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSb_H0rs_quant[4], at = TSb_H0rs_quant[4], side = 1, col = "red")

#check allocation compared to MultinomialTS
hist(TSb_H0rs_alloc[,1], col = "blue", breaks = 150, freq = F, xlab = " ", ylab = "", main = "Empirical allocation of Arm 1")
hist(TSmn_H0rs_alloc[,1], col = "red", add = T, breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")
hist(TSn_H0rs_alloc[,1], col = "black", add = T, breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")


TSb_H0rs_y0 <- data.frame(matrix(unlist(lapply(TSb_H0_rs, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSb_H0rs_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical distribution of Y (H0: mu0=mu1=4)")
abline(v=mean(TSb_H0rs_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_rs - mean(TSb_H0rs_y0[,1]), 4)))

TSb_H0rs_y1 <- data.frame(matrix(unlist(lapply(TSb_H0_rs, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSb_H0rs_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical distribution of Y (H0: mu0=mu1=4)")
abline(v=mean(TSb_H0rs_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_rs - mean(TSb_H0rs_y1[,1]), 4)))

# |-- H0. Binary TS (cutoff 7) ---------------------------------------------------

set.seed(12345)
TSb1_H0_rs = pbreplicate(Nsim, list(TS_BetaBern(N=myN, n_arms = 2, y_categories = 1:7,
                                                 beta_priors = list(alpha=c(1,1), beta=c(1,1)), 
                                                 true_params = list(p0 = true_p0_rs, p1 = true_p0_rs),
                                                 cutoff_value = 7, compute_rho = F)))

#looking at average values
TSb1_H0rs_alloc <- data.frame(matrix(unlist(lapply(TSb1_H0_rs, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSb1_H0rs_alloc[,1], breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")
mean(TSb1_H0rs_alloc[,1]); mean(TSb1_H0rs_alloc[,2])
TSb1_H0rs_quant = quantile(TSb1_H0rs_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
abline(v = TSb1_H0rs_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSb1_H0rs_quant[1], at = TSb1_H0rs_quant[1], side = 1, col = "red")
abline(v = TSb1_H0rs_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSb1_H0rs_quant[4], at = TSb1_H0rs_quant[4], side = 1, col = "red")

#check allocation compared to MultinomialTS
hist(TSn_H0rs_alloc[,1], xlim=c(0,1), ylim=c(0,2), col = "black", breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")
hist(TSmn_H0rs_alloc[,1], col = "red", add = T, breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")
hist(TSb_H0rs_alloc[,1], col = "green", add = T, breaks = 150, freq = F, xlab = " ", ylab = "", main = "Empirical allocation of Arm 1")
hist(TSb1_H0rs_alloc[,1], xlim=c(0,1), add = T, col = "blue", breaks = 150, freq = F, xlab = " ", ylab = "", main = "Empirical allocation of Arm 1")


TSb1_H0rs_y0 <- data.frame(matrix(unlist(lapply(TSb1_H0_rs, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSb1_H0rs_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical distribution of Y (H0: mu0=mu1=4)")
abline(v=mean(TSb1_H0rs_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_rs - mean(TSb1_H0rs_y0[,1]), 4)))

TSb1_H0rs_y1 <- data.frame(matrix(unlist(lapply(TSb1_H0_rs, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSb1_H0rs_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical distribution of Y (H0: mu0=mu1=4)")
abline(v=mean(TSb1_H0rs_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_rs - mean(TSb1_H0rs_y1[,1]), 4)))


# |-- H0: Compute and Plot performance data: Prop Optimal action -------------------


data_alloc_H0rs <- data.frame(
  type = c( rep("Multinomial TS", Nsim), rep("Normal TS", Nsim), 
            rep("Binary TS (cutoff = 6)", Nsim), rep("Binary TS (cutoff = 7)", Nsim), rep("Oracle",160) ),
  value = c(TSmn_H0rs_alloc[,1], TSn_H0rs_alloc[,1], 
            TSb_H0rs_alloc[,1], TSb1_H0rs_alloc[,1], rep(0.5, 160))
)

library(hrbrthemes)
# Represent it
p_mn <- data_alloc_H0rs[data_alloc_H0rs$type=="Multinomial TS",] %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram(bins = 100, color="#F8766D", alpha=0.7, position = 'identity') +
  scale_fill_manual(values=c("#F8766D")) +
  #ggtitle("Empirical allocation of Arm 1") +
  #scale_alpha_manual(name = "category", values = c(0.9, 0.2, 1)) +
  theme_ipsum() +
  #theme(plot.title = element_text(size=15)) +
  labs(fill="") + ylab(" ") +  ylim(0,200) + xlim(0,1) +
  xlab(expression(T[1]*"/"*T)) +
  theme(legend.position=c(0.5, 0.95),
        legend.text=element_text(size=10)) +
  theme(axis.title.x = element_text(hjust = 0.5, face = "bold")) +
  theme(plot.margin=unit(c(0.2,0,0,0), "cm"))

p_n <- data_alloc_H0rs[data_alloc_H0rs$type=="Normal TS",] %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram(bins = 100, color="#D16DBE", alpha=0.7, position = 'identity') +
  scale_fill_manual(values=c("#D16DBE")) +
  #ggtitle("Empirical allocation of Arm 1") +
  #scale_alpha_manual(name = "category", values = c(0.9, 0.2, 1)) +
  theme_ipsum() +
  #theme(plot.title = element_text(size=15)) +
  labs(fill="") + ylab(" ") +  ylim(0,200) + xlim(0,1) +
  xlab(expression(T[1]*"/"*T)) +
  theme(legend.position=c(0.5, 0.95),
        legend.text=element_text(size=10)) +
  theme(axis.title.x = element_text(hjust = 0.5, face = "bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(plot.margin=unit(c(0.2,0,0,0), "cm"))

p_b <- data_alloc_H0rs[data_alloc_H0rs$type=="Binary TS (cutoff = 6)",] %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram(bins = 100, color="#9590FF", alpha=0.7, position = 'identity') +
  scale_fill_manual(values=c("#9590FF")) +
  #ggtitle("Empirical allocation of Arm 1") +
  #scale_alpha_manual(name = "category", values = c(0.9, 0.2, 1)) +
  theme_ipsum() +
  #theme(plot.title = element_text(size=15)) +
  labs(fill="") + ylab(" ") +  ylim(0,200) + xlim(0,1) +
  xlab(expression(T[1]*"/"*T)) +
  theme(legend.position=c(0.6, 0.95),
        legend.text=element_text(size=10)) +
  theme(axis.title.x = element_text(hjust = 0.5, face = "bold")) +
  theme(plot.margin=unit(c(0.2,0,0,0), "cm"))

p_b1 <- data_alloc_H0rs[data_alloc_H0rs$type=="Binary TS (cutoff = 7)",] %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram(bins = 100, color="#529EFF", alpha=0.7, position = 'identity') +
  scale_fill_manual(values=c("#529EFF")) +
  #ggtitle("Empirical allocation of Arm 1") +
  #scale_alpha_manual(name = "category", values = c(0.9, 0.2, 1)) +
  theme_ipsum() +
  #theme(plot.title = element_text(size=15)) +
  labs(fill="") + ylab(" ") +  ylim(0,200) + xlim(0,1) +
  xlab(expression(T[1]*"/"*T)) +
  theme(legend.position=c(0.6, 0.95),
        legend.text=element_text(size=10)) +
  theme(axis.title.x = element_text(hjust = 0.5, face = "bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(plot.margin=unit(c(0.2,0,0,0), "cm"))

ggarrange(p_mn, p_n, p_b, p_b1, ncol=2, nrow=2)

# H0. Left / Negatively skewed distribution ---------------------------------------------------

true_p0_ls = c(0.05, 0.05, 0.1, 0.15, 0.15, 0.3, 0.2)
plot(true_p0_ls, type = "h", lwd = 2, main = "True distr (H0: p0=p1; skew mu=5)", xlab="Y (multinomial)", ylab = "")
mean_p0_ls = sum(true_p0_ls*1:7)

# |-- H0. Multinomial TS ---------------------------------------------------

set.seed(12345)
TSmn_H0ls = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_ls, p1 = true_p0_ls), compute_rho = F)))

#looking at average values
TSmn_H0ls_alloc <- data.frame(matrix(unlist(lapply(TSmn_H0ls, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSmn_H0ls_alloc[,1], breaks = 150, freq = F, xlab = "Ordered Multinomial", ylab = "Emprirical Allocation Arm 0", main = "Empirical alloc (H0: p0=p1 ls mu=4)")
mean(TSmn_H0ls_alloc[,1]); mean(TSmn_H0ls_alloc[,2])
TSmn_H0ls_quant = quantile(TSmn_H0ls_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
abline(v = TSmn_H0ls_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSmn_H0ls_quant[1], at = TSmn_H0ls_quant[1], side = 1, col = "red")
abline(v = TSmn_H0ls_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSmn_H0ls_quant[4], at = TSmn_H0ls_quant[4], side = 1, col = "red")

c(mean(TSmn_H0ls_alloc[,1]<0.05)*100, 
  mean(TSmn_H0ls_alloc[,1]<0.1)*100, 
  mean(TSmn_H0ls_alloc[,1]<=0.55 & TSmn_H0ls_alloc[,1]>=0.45)*100, 
  mean(TSmn_H0ls_alloc[,1]<=0.6 & TSmn_H0ls_alloc[,1]>=0.4)*100, 
  mean(TSmn_H0ls_alloc[,1]>0.9)*100, 
  mean(TSmn_H0ls_alloc[,1]>0.95)*100)

TSmn_H0ls_y0 <- data.frame(matrix(unlist(lapply(TSmn_H0ls, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H0ls_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical dist (H0: p0=p1 ls mu=4)")
abline(v=mean(TSmn_H0ls_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_ls - mean(TSmn_H0ls_y0[,1]), 3)))

TSmn_H0ls_y1 <- data.frame(matrix(unlist(lapply(TSmn_H0ls, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H0ls_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical dist (H0: p0=p1 ls mu=4)")
abline(v=mean(TSmn_H0ls_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_ls - mean(TSmn_H0ls_y1[,1]), 3)))

# |-- H0. Normal TS ---------------------------------------------------

mycov_H0 = c()
for(i in 1:dim(comb)[1]){
  mycov_H0[i] = prod(comb[i,])*prod(true_p0_ls[comb[i,]])
}
var_H0 = sum(c(1:7)^2*true_p0_ls*(1-true_p0_ls)) - 2*sum(mycov_H0)

set.seed(12345)
TSn_H0ls = pbreplicate(Nsim, list(TS_NormNorm(N=myN, true_params = list(p0 = true_p0_ls, p1 = true_p0_ls), 
                                              true_vars = c(var_H0, var_H0), compute_rho = F)))


#looking at average values
TSn_H0ls_alloc <- data.frame(matrix(unlist(lapply(TSn_H0ls, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSn_H0ls_alloc[,1], breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")
mean(TSn_H0ls_alloc[,1]); mean(TSn_H0ls_alloc[,2])
TSn_H0ls_quant = quantile(TSn_H0ls_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
abline(v = TSn_H0ls_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSn_H0ls_quant[1], at = TSn_H0ls_quant[1], side = 1, col = "red")
abline(v = TSn_H0ls_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSn_H0ls_quant[4], at = TSn_H0ls_quant[4], side = 1, col = "red")

#check allocation compared to MultinomialTS
hist(TSn_H0ls_alloc[,1], col = "blue", breaks = 150, freq = F, xlab = " ", ylab = "", main = "Empirical allocation of Arm 1")
hist(TSmn_H0ls_alloc[,1],col = "red", add = T, breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")


TSn_H0ls_y0 <- data.frame(matrix(unlist(lapply(TSn_H0ls, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSn_H0ls_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical distribution of Y (H0: mu0=mu1=4)")
abline(v=mean(TSn_H0ls_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_ls - mean(TSn_H0ls_y0[,1]), 4)))

TSn_H0ls_y1 <- data.frame(matrix(unlist(lapply(TSn_H0ls, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSn_H0ls_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical distribution of Y (H0: mu0=mu1=4)")
abline(v=mean(TSn_H0ls_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_ls - mean(TSn_H0ls_y1[,1]), 4)))


# |-- H0. Binary TS (cutoff 6) ---------------------------------------------------

set.seed(12345)
TSb_H0_ls = pbreplicate(Nsim, list(TS_BetaBern(N=myN, n_arms = 2, y_categories = 1:7,
                                               beta_priors = list(alpha=c(1,1), beta=c(1,1)), 
                                               true_params = list(p0 = true_p0_ls, p1 = true_p0_ls),
                                               cutoff_value = 6, compute_rho = F)))


#looking at average values
TSb_H0ls_alloc <- data.frame(matrix(unlist(lapply(TSb_H0_ls, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSb_H0ls_alloc[,1], breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")
mean(TSb_H0ls_alloc[,1]); mean(TSb_H0ls_alloc[,2])
TSb_H0ls_quant = quantile(TSb_H0ls_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
abline(v = TSb_H0ls_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSb_H0ls_quant[1], at = TSb_H0ls_quant[1], side = 1, col = "red")
abline(v = TSb_H0ls_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSb_H0ls_quant[4], at = TSb_H0ls_quant[4], side = 1, col = "red")

#check allocation compared to MultinomialTS
hist(TSb_H0ls_alloc[,1], col = "blue", breaks = 150, freq = F, xlab = " ", ylab = "", main = "Empirical allocation of Arm 1")
hist(TSmn_H0ls_alloc[,1], col = "red", add = T, breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")
hist(TSn_H0ls_alloc[,1], col = "black", add = T, breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")


TSb_H0ls_y0 <- data.frame(matrix(unlist(lapply(TSb_H0_ls, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSb_H0ls_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical distribution of Y (H0: mu0=mu1=4)")
abline(v=mean(TSb_H0ls_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_ls - mean(TSb_H0ls_y0[,1]), 4)))

TSb_H0ls_y1 <- data.frame(matrix(unlist(lapply(TSb_H0_ls, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSb_H0ls_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical distribution of Y (H0: mu0=mu1=4)")
abline(v=mean(TSb_H0ls_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_ls - mean(TSb_H0ls_y1[,1]), 4)))

# |-- H0. Binary TS (cutoff 7) ---------------------------------------------------

set.seed(12345)
TSb1_H0_ls = pbreplicate(Nsim, list(TS_BetaBern(N=myN, n_arms = 2, y_categories = 1:7,
                                                beta_priors = list(alpha=c(1,1), beta=c(1,1)), 
                                                true_params = list(p0 = true_p0_ls, p1 = true_p0_ls),
                                                cutoff_value = 7, compute_rho = F)))

#looking at average values
TSb1_H0ls_alloc <- data.frame(matrix(unlist(lapply(TSb1_H0_ls, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSb1_H0ls_alloc[,1], breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")
mean(TSb1_H0ls_alloc[,1]); mean(TSb1_H0ls_alloc[,2])
TSb1_H0ls_quant = quantile(TSb1_H0ls_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
abline(v = TSb1_H0ls_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSb1_H0ls_quant[1], at = TSb1_H0ls_quant[1], side = 1, col = "red")
abline(v = TSb1_H0ls_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSb1_H0ls_quant[4], at = TSb1_H0ls_quant[4], side = 1, col = "red")

#check allocation compared to MultinomialTS
hist(TSn_H0ls_alloc[,1], xlim=c(0,1), ylim=c(0,2), col = "black", breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")
hist(TSmn_H0ls_alloc[,1], col = "red", add = T, breaks = 150, freq = F, xlab = "Emprirical Allocation (Normal model)", main = "Empirical allocation (H0: mu0=mu1=4)")
hist(TSb_H0ls_alloc[,1], col = "green", add = T, breaks = 150, freq = F, xlab = " ", ylab = "", main = "Empirical allocation of Arm 1")
hist(TSb1_H0ls_alloc[,1], xlim=c(0,1), add = T, col = "blue", breaks = 150, freq = F, xlab = " ", ylab = "", main = "Empirical allocation of Arm 1")


TSb1_H0ls_y0 <- data.frame(matrix(unlist(lapply(TSb1_H0_ls, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSb1_H0ls_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical distribution of Y (H0: mu0=mu1=4)")
abline(v=mean(TSb1_H0ls_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_ls - mean(TSb1_H0ls_y0[,1]), 4)))

TSb1_H0ls_y1 <- data.frame(matrix(unlist(lapply(TSb1_H0_ls, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSb1_H0ls_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical distribution of Y (H0: mu0=mu1=4)")
abline(v=mean(TSb1_H0ls_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_ls - mean(TSb1_H0ls_y1[,1]), 4)))


# |-- H0: Compute and Plot performance data: Prop Optimal action -------------------


data_alloc_H0ls <- data.frame(
  type = c( rep("Multinomial TS", Nsim), rep("Normal TS", Nsim), 
            rep("Binary TS (cutoff = 6)", Nsim), rep("Binary TS (cutoff = 7)", Nsim), rep("Oracle",160) ),
  value = c(TSmn_H0ls_alloc[,1], TSn_H0ls_alloc[,1], 
            TSb_H0ls_alloc[,1], TSb1_H0ls_alloc[,1], rep(0.5, 160))
)

library(hrbrthemes)
# Represent it
p_mn<- data_alloc_H0ls[data_alloc_H0ls$type=="Multinomial TS",] %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram(bins = 100, color="#F8766D", alpha=0.7, position = 'identity') +
  scale_fill_manual(values=c("#F8766D")) +
  #ggtitle("Empirical allocation of Arm 1") +
  #scale_alpha_manual(name = "category", values = c(0.9, 0.2, 1)) +
  theme_ipsum() +
  #theme(plot.title = element_text(size=15)) +
  labs(fill="") + ylab(" ") +  ylim(0,200) + xlim(0,1) +
  xlab(expression(T[1]*"/"*T)) +
  theme(legend.position=c(0.5, 0.95),
        legend.text=element_text(size=10)) +
  theme(axis.title.x = element_text(hjust = 0.5, face = "bold")) +
  theme(plot.margin=unit(c(0.2,0,0,0), "cm"))

p_n <- data_alloc_H0ls[data_alloc_H0ls$type=="Normal TS",] %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram(bins = 100, color="#D16DBE", alpha=0.7, position = 'identity') +
  scale_fill_manual(values=c("#D16DBE")) +
  #ggtitle("Empirical allocation of Arm 1") +
  #scale_alpha_manual(name = "category", values = c(0.9, 0.2, 1)) +
  theme_ipsum() +
  #theme(plot.title = element_text(size=15)) +
  labs(fill="") + ylab(" ") +  ylim(0,200) + xlim(0,1) +
  xlab(expression(T[1]*"/"*T)) +
  theme(legend.position=c(0.5, 0.95),
        legend.text=element_text(size=10)) +
  theme(axis.title.x = element_text(hjust = 0.5, face = "bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(plot.margin=unit(c(0.2,0,0,0), "cm"))

p_b <- data_alloc_H0ls[data_alloc_H0ls$type=="Binary TS (cutoff = 6)",] %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram(bins = 100, color="#9590FF", alpha=0.7, position = 'identity') +
  scale_fill_manual(values=c("#9590FF")) +
  #ggtitle("Empirical allocation of Arm 1") +
  #scale_alpha_manual(name = "category", values = c(0.9, 0.2, 1)) +
  theme_ipsum() +
  #theme(plot.title = element_text(size=15)) +
  labs(fill="") + ylab(" ") +  ylim(0,200) + xlim(0,1) +
  xlab(expression(T[1]*"/"*T)) +
  theme(legend.position=c(0.6, 0.95),
        legend.text=element_text(size=10)) +
  theme(axis.title.x = element_text(hjust = 0.5, face = "bold")) +
  theme(plot.margin=unit(c(0.2,0,0,0), "cm"))

p_b1 <- data_alloc_H0ls[data_alloc_H0ls$type=="Binary TS (cutoff = 7)",] %>%
  ggplot( aes(x=value, fill=type)) +
  geom_histogram(bins = 100, color="#529EFF", alpha=0.7, position = 'identity') +
  scale_fill_manual(values=c("#529EFF")) +
  #ggtitle("Empirical allocation of Arm 1") +
  #scale_alpha_manual(name = "category", values = c(0.9, 0.2, 1)) +
  theme_ipsum() +
  #theme(plot.title = element_text(size=15)) +
  labs(fill="") + ylab(" ") +  ylim(0,200) + xlim(0,1) +
  xlab(expression(T[1]*"/"*T)) +
  theme(legend.position=c(0.6, 0.95),
        legend.text=element_text(size=10)) +
  theme(axis.title.x = element_text(hjust = 0.5, face = "bold"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  theme(plot.margin=unit(c(0.2,0,0,0), "cm"))

ggarrange(p_mn, p_n, p_b, p_b1, ncol=2, nrow=2)

p_mn1 = p_mn1 + ggtitle("Scenario 3: right-skewness") + 
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  theme(legend.position="none")
p_mn2 = p_mn2 + ggtitle("Scenario 4: left-skewness") + 
  theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold")) +
  theme(legend.position="none")
ggarrange(p_mn1, p_mn2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")


# H0. Modified Multinomial TS with augmented support -------------------------------------------------

# |-- H0. s = 1 ---------------------------------------------------
s=1
true_p0_1s = c(true_p0_ls, rep(0,s))
plot(true_p0_1s, type = "h", lwd = 2, main = "True distr (H0: p0=p1; skew mu=5)", xlab="Y (multinomial)", ylab = "")
mean_p0_1s = sum(true_p0_ls*1:length(true_p0_3s))


set.seed(12345)
TSmn_H01s = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_1s, p1 = true_p0_1s), 
                                               compute_rho = F,
                                               y_categories = 1:length(true_p0_1s),
                                               prior_params = list(arm0 = rep(1,length(true_p0_1s)), 
                                                                   arm1 = rep(1,length(true_p0_1s))))))

#looking at average values
TSmn_H01s_alloc <- data.frame(matrix(unlist(lapply(TSmn_H01s, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSmn_H01s_alloc[,1], breaks = 150, freq = F, xlab = "Ordered Multinomial", ylab = "Emprirical Allocation Arm 0", main = "Empirical alloc (H0: p0=p1 1s mu=4)")
mean(TSmn_H01s_alloc[,1]); mean(TSmn_H01s_alloc[,2])
TSmn_H01s_quant = quantile(TSmn_H01s_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
abline(v = TSmn_H01s_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSmn_H01s_quant[1], at = TSmn_H01s_quant[1], side = 1, col = "red")
abline(v = TSmn_H01s_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSmn_H01s_quant[4], at = TSmn_H01s_quant[4], side = 1, col = "red")

c(mean(TSmn_H01s_alloc[,1]<0.05)*100, 
  mean(TSmn_H01s_alloc[,1]<0.1)*100, 
  mean(TSmn_H01s_alloc[,1]<=0.55 & TSmn_H01s_alloc[,1]>=0.45)*100, 
  mean(TSmn_H01s_alloc[,1]<=0.6 & TSmn_H01s_alloc[,1]>=0.4)*100, 
  mean(TSmn_H01s_alloc[,1]>0.9)*100, 
  mean(TSmn_H01s_alloc[,1]>0.95)*100)

TSmn_H01s_y0 <- data.frame(matrix(unlist(lapply(TSmn_H01s, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H01s_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical dist (H0: p0=p1 1s mu=4)")
abline(v=mean(TSmn_H01s_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_1s - mean(TSmn_H01s_y0[,1]), 3)))

TSmn_H01s_y1 <- data.frame(matrix(unlist(lapply(TSmn_H01s, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H01s_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical dist (H0: p0=p1 1s mu=4)")
abline(v=mean(TSmn_H01s_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_1s - mean(TSmn_H01s_y1[,1]), 3)))


# |-- H0. s = 3 ---------------------------------------------------
s=3
true_p0_3s = c(true_p0_ls, rep(0,s))
plot(true_p0_3s, type = "h", lwd = 2, main = "True distr (H0: p0=p1; skew mu=5)", xlab="Y (multinomial)", ylab = "")
mean_p0_3s = sum(true_p0_3s*1:length(true_p0_3s))


set.seed(12345)
TSmn_H03s = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_3s, p1 = true_p0_3s), 
                                               compute_rho = F,
                                               y_categories = 1:length(true_p0_3s),
                                               prior_params = list(arm0 = rep(1,length(true_p0_3s)), arm1 = rep(1,length(true_p0_3s))))))

#looking at average values
TSmn_H03s_alloc <- data.frame(matrix(unlist(lapply(TSmn_H03s, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSmn_H03s_alloc[,1], breaks = 150, freq = F, xlab = "Ordered Multinomial", ylab = "Emprirical Allocation Arm 0", main = "Empirical alloc (H0: p0=p1 3s mu=4)")
mean(TSmn_H03s_alloc[,1]); mean(TSmn_H03s_alloc[,2])
TSmn_H03s_quant = quantile(TSmn_H03s_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
abline(v = TSmn_H03s_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSmn_H03s_quant[1], at = TSmn_H03s_quant[1], side = 1, col = "red")
abline(v = TSmn_H03s_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSmn_H03s_quant[4], at = TSmn_H03s_quant[4], side = 1, col = "red")

c(mean(TSmn_H03s_alloc[,1]<0.05)*100, 
  mean(TSmn_H03s_alloc[,1]<0.1)*100, 
  mean(TSmn_H03s_alloc[,1]<=0.55 & TSmn_H03s_alloc[,1]>=0.45)*100, 
  mean(TSmn_H03s_alloc[,1]<=0.6 & TSmn_H03s_alloc[,1]>=0.4)*100, 
  mean(TSmn_H03s_alloc[,1]>0.9)*100, 
  mean(TSmn_H03s_alloc[,1]>0.95)*100)

TSmn_H03s_y0 <- data.frame(matrix(unlist(lapply(TSmn_H03s, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H03s_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical dist (H0: p0=p1 3s mu=4)")
abline(v=mean(TSmn_H03s_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_3s - mean(TSmn_H03s_y0[,1]), 3)))

TSmn_H03s_y1 <- data.frame(matrix(unlist(lapply(TSmn_H03s, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H03s_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical dist (H0: p0=p1 3s mu=4)")
abline(v=mean(TSmn_H03s_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_3s - mean(TSmn_H03s_y1[,1]), 3)))


# |-- H0. s = 5 ---------------------------------------------------
s=5
true_p0_5s = c(true_p0_ls, rep(0,s))
plot(true_p0_5s, type = "h", lwd = 2, main = "True distr (H0: p0=p1; skew mu=5)", xlab="Y (multinomial)", ylab = "")
mean_p0_5s = sum(true_p0_5s*1:length(true_p0_5s))


set.seed(12345)
TSmn_H05s = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_5s, p1 = true_p0_5s), 
                                               compute_rho = F,
                                               y_categories = 1:length(true_p0_5s),
                                               prior_params = list(arm0 = rep(1,length(true_p0_5s)), arm1 = rep(1,length(true_p0_5s))))))

#looking at average values
TSmn_H05s_alloc <- data.frame(matrix(unlist(lapply(TSmn_H05s, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSmn_H05s_alloc[,1], breaks = 150, freq = F, xlab = "Ordered Multinomial", ylab = "Emprirical Allocation Arm 0", main = "Empirical alloc (H0: p0=p1 5s mu=4)")
mean(TSmn_H05s_alloc[,1]); mean(TSmn_H05s_alloc[,2])
TSmn_H05s_quant = quantile(TSmn_H05s_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
abline(v = TSmn_H05s_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSmn_H05s_quant[1], at = TSmn_H05s_quant[1], side = 1, col = "red")
abline(v = TSmn_H05s_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSmn_H05s_quant[4], at = TSmn_H05s_quant[4], side = 1, col = "red")

c(mean(TSmn_H05s_alloc[,1]<0.05)*100, 
  mean(TSmn_H05s_alloc[,1]<0.1)*100, 
  mean(TSmn_H05s_alloc[,1]<=0.55 & TSmn_H05s_alloc[,1]>=0.45)*100, 
  mean(TSmn_H05s_alloc[,1]<=0.6 & TSmn_H05s_alloc[,1]>=0.4)*100, 
  mean(TSmn_H05s_alloc[,1]>0.9)*100, 
  mean(TSmn_H05s_alloc[,1]>0.95)*100)

TSmn_H05s_y0 <- data.frame(matrix(unlist(lapply(TSmn_H05s, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H05s_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical dist (H0: p0=p1 5s mu=4)")
abline(v=mean(TSmn_H05s_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_5s - mean(TSmn_H05s_y0[,1]), 3)))

TSmn_H05s_y1 <- data.frame(matrix(unlist(lapply(TSmn_H05s, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H05s_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical dist (H0: p0=p1 5s mu=4)")
abline(v=mean(TSmn_H05s_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_5s - mean(TSmn_H05s_y1[,1]), 3)))

# |-- H0. s = 10 ---------------------------------------------------
s=10
true_p0_10s = c(true_p0_ls, rep(0,s))
plot(true_p0_10s, type = "h", lwd = 2, main = "True distr (H0: p0=p1; skew mu=5)", xlab="Y (multinomial)", ylab = "")
mean_p0_10s = sum(true_p0_10s*1:length(true_p0_10s))


set.seed(12345)
TSmn_H010s = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_10s, p1 = true_p0_10s), 
                                               compute_rho = F,
                                               y_categories = 1:length(true_p0_10s),
                                               prior_params = list(arm0 = rep(1,length(true_p0_10s)), arm1 = rep(1,length(true_p0_10s))))))

#looking at average values
TSmn_H010s_alloc <- data.frame(matrix(unlist(lapply(TSmn_H010s, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSmn_H010s_alloc[,1], breaks = 150, freq = F, xlab = "Ordered Multinomial", ylab = "Emprirical Allocation Arm 0", main = "Empirical alloc (H0: p0=p1 10s mu=4)")
mean(TSmn_H010s_alloc[,1]); mean(TSmn_H010s_alloc[,2])
TSmn_H010s_quant = quantile(TSmn_H010s_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
abline(v = TSmn_H010s_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSmn_H010s_quant[1], at = TSmn_H010s_quant[1], side = 1, col = "red")
abline(v = TSmn_H010s_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSmn_H010s_quant[4], at = TSmn_H010s_quant[4], side = 1, col = "red")

c(mean(TSmn_H010s_alloc[,1]<0.05)*100, 
  mean(TSmn_H010s_alloc[,1]<0.1)*100, 
  mean(TSmn_H010s_alloc[,1]<=0.55 & TSmn_H010s_alloc[,1]>=0.45)*100, 
  mean(TSmn_H010s_alloc[,1]<=0.6 & TSmn_H010s_alloc[,1]>=0.4)*100, 
  mean(TSmn_H010s_alloc[,1]>0.9)*100, 
  mean(TSmn_H010s_alloc[,1]>0.95)*100)

TSmn_H010s_y0 <- data.frame(matrix(unlist(lapply(TSmn_H010s, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H010s_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical dist (H0: p0=p1 10s mu=4)")
abline(v=mean(TSmn_H010s_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_10s - mean(TSmn_H010s_y0[,1]), 3)))

TSmn_H010s_y1 <- data.frame(matrix(unlist(lapply(TSmn_H010s, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H010s_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical dist (H0: p0=p1 10s mu=4)")
abline(v=mean(TSmn_H010s_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_10s - mean(TSmn_H010s_y1[,1]), 3)))

# |-- H0. s = 20 ---------------------------------------------------
s=20
true_p0_20s = c(true_p0_ls, rep(0,s))
plot(true_p0_20s, type = "h", lwd = 2, main = "True distr (H0: p0=p1; skew mu=5)", xlab="Y (multinomial)", ylab = "")
mean_p0_20s = sum(true_p0_20s*1:length(true_p0_20s))


set.seed(12345)
TSmn_H020s = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_20s, p1 = true_p0_20s), 
                                                compute_rho = F,
                                                y_categories = 1:length(true_p0_20s),
                                                prior_params = list(arm0 = rep(1,length(true_p0_20s)), arm1 = rep(1,length(true_p0_20s))))))

#looking at average values
TSmn_H020s_alloc <- data.frame(matrix(unlist(lapply(TSmn_H020s, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSmn_H020s_alloc[,1], breaks = 150, freq = F, xlab = "Ordered Multinomial", ylab = "Emprirical Allocation Arm 0", main = "Empirical alloc (H0: p0=p1 20s mu=4)")
mean(TSmn_H020s_alloc[,1]); mean(TSmn_H020s_alloc[,2])
TSmn_H020s_quant = quantile(TSmn_H020s_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
abline(v = TSmn_H020s_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSmn_H020s_quant[1], at = TSmn_H020s_quant[1], side = 1, col = "red")
abline(v = TSmn_H020s_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSmn_H020s_quant[4], at = TSmn_H020s_quant[4], side = 1, col = "red")

c(mean(TSmn_H020s_alloc[,1]<0.05)*100, 
  mean(TSmn_H020s_alloc[,1]<0.1)*100, 
  mean(TSmn_H020s_alloc[,1]<=0.55 & TSmn_H020s_alloc[,1]>=0.45)*100, 
  mean(TSmn_H020s_alloc[,1]<=0.6 & TSmn_H020s_alloc[,1]>=0.4)*100, 
  mean(TSmn_H020s_alloc[,1]>0.9)*100, 
  mean(TSmn_H020s_alloc[,1]>0.95)*100)

TSmn_H020s_y0 <- data.frame(matrix(unlist(lapply(TSmn_H020s, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H020s_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical dist (H0: p0=p1 20s mu=4)")
abline(v=mean(TSmn_H020s_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_20s - mean(TSmn_H020s_y0[,1]), 3)))

TSmn_H020s_y1 <- data.frame(matrix(unlist(lapply(TSmn_H020s, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H020s_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical dist (H0: p0=p1 20s mu=4)")
abline(v=mean(TSmn_H020s_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_20s - mean(TSmn_H020s_y1[,1]), 3)))

# |-- H0. s = 50 ---------------------------------------------------
s=50
true_p0_50s = c(true_p0_ls, rep(0,s))
plot(true_p0_50s, type = "h", lwd = 2, main = "True distr (H0: p0=p1; skew mu=5)", xlab="Y (multinomial)", ylab = "")
mean_p0_50s = sum(true_p0_50s*1:length(true_p0_50s))


set.seed(12345)
TSmn_H050s = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_50s, p1 = true_p0_50s), 
                                                compute_rho = F,
                                                y_categories = 1:length(true_p0_50s),
                                                prior_params = list(arm0 = rep(1,length(true_p0_50s)), arm1 = rep(1,length(true_p0_50s))))))

#looking at average values
TSmn_H050s_alloc <- data.frame(matrix(unlist(lapply(TSmn_H050s, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSmn_H050s_alloc[,1], breaks = 150, freq = F, xlab = "Ordered Multinomial", ylab = "Emprirical Allocation Arm 0", main = "Empirical alloc (H0: p0=p1 50s mu=4)")
mean(TSmn_H050s_alloc[,1]); mean(TSmn_H050s_alloc[,2])
TSmn_H050s_quant = quantile(TSmn_H050s_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
abline(v = TSmn_H050s_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSmn_H050s_quant[1], at = TSmn_H050s_quant[1], side = 1, col = "red")
abline(v = TSmn_H050s_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSmn_H050s_quant[4], at = TSmn_H050s_quant[4], side = 1, col = "red")

c(mean(TSmn_H050s_alloc[,1]<0.05)*100, 
  mean(TSmn_H050s_alloc[,1]<0.1)*100, 
  mean(TSmn_H050s_alloc[,1]<=0.55 & TSmn_H050s_alloc[,1]>=0.45)*100, 
  mean(TSmn_H050s_alloc[,1]<=0.6 & TSmn_H050s_alloc[,1]>=0.4)*100, 
  mean(TSmn_H050s_alloc[,1]<=0.51 & TSmn_H050s_alloc[,1]>=0.49)*100, 
  mean(TSmn_H050s_alloc[,1]>0.9)*100, 
  mean(TSmn_H050s_alloc[,1]>0.95)*100)

TSmn_H050s_y0 <- data.frame(matrix(unlist(lapply(TSmn_H050s, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H050s_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical dist (H0: p0=p1 50s mu=4)")
abline(v=mean(TSmn_H050s_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_50s - mean(TSmn_H050s_y0[,1]), 3)))

TSmn_H050s_y1 <- data.frame(matrix(unlist(lapply(TSmn_H050s, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H050s_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical dist (H0: p0=p1 50s mu=4)")
abline(v=mean(TSmn_H050s_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_50s - mean(TSmn_H050s_y1[,1]), 3)))


# |-- H0. s = 100 ---------------------------------------------------
s=100
true_p0_100s = c(true_p0_ls, rep(0,s))
plot(true_p0_100s, type = "h", lwd = 2, main = "True distr (H0: p0=p1; skew mu=5)", xlab="Y (multinomial)", ylab = "")
mean_p0_100s = sum(true_p0_100s*1:length(true_p0_100s))


set.seed(12345)
TSmn_H0100s = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_100s, p1 = true_p0_100s), 
                                                compute_rho = F,
                                                y_categories = 1:length(true_p0_100s),
                                                prior_params = list(arm0 = rep(1,length(true_p0_100s)), arm1 = rep(1,length(true_p0_100s))))))

#looking at average values
TSmn_H0100s_alloc <- data.frame(matrix(unlist(lapply(TSmn_H0100s, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
hist(TSmn_H0100s_alloc[,1], breaks = 150, freq = F, xlab = "Ordered Multinomial", ylab = "Emprirical Allocation Arm 0", main = "Empirical alloc (H0: p0=p1 100s mu=4)")
mean(TSmn_H0100s_alloc[,1]); mean(TSmn_H0100s_alloc[,2])
TSmn_H0100s_quant = quantile(TSmn_H0100s_alloc[,1], probs = c(0.05, 0.1, 0.9, 0.95))
abline(v = TSmn_H0100s_quant[1], col="red")
mtext("5%", at = 0.05, col = "red")
mtext(TSmn_H0100s_quant[1], at = TSmn_H0100s_quant[1], side = 1, col = "red")
abline(v = TSmn_H0100s_quant[4], col="red")
mtext("95%", at = 0.95, col = "red")
mtext(TSmn_H0100s_quant[4], at = TSmn_H0100s_quant[4], side = 1, col = "red")

c(mean(TSmn_H0100s_alloc[,1]<0.05)*100, 
  mean(TSmn_H0100s_alloc[,1]<0.1)*100, 
  mean(TSmn_H0100s_alloc[,1]<=0.55 & TSmn_H0100s_alloc[,1]>=0.45)*100, 
  mean(TSmn_H0100s_alloc[,1]<=0.6 & TSmn_H0100s_alloc[,1]>=0.4)*100, 
  mean(TSmn_H0100s_alloc[,1]>0.9)*100, 
  mean(TSmn_H0100s_alloc[,1]>0.95)*100)

TSmn_H0100s_y0 <- data.frame(matrix(unlist(lapply(TSmn_H0100s, function(x) mean(x$y[x$a_hat==0], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H0100s_y0[,1], breaks = 100, freq = F, xlab = "Y[0]", main = "Empirical dist (H0: p0=p1 100s mu=4)")
abline(v=mean(TSmn_H0100s_y0[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_100s - mean(TSmn_H0100s_y0[,1]), 3)))

TSmn_H0100s_y1 <- data.frame(matrix(unlist(lapply(TSmn_H0100s, function(x) mean(x$y[x$a_hat==1], na.rm=T))), nrow=Nsim, byrow=TRUE))
hist(TSmn_H0100s_y1[,1], breaks = 100, freq = F, xlab = "Y[1]", main = "Empirical dist (H0: p0=p1 100s mu=4)")
abline(v=mean(TSmn_H0100s_y1[,1]), col="red")
mtext(paste("bias = ", round(mean_p0_100s - mean(TSmn_H0100s_y1[,1]), 3)))


# |-- Make ridge plot for different s -------------------

library(ggplot2)
library(ggridges)

data_alloc_s <- data.frame(
  type = c( rep("Multinomial-TS (s = 0)", Nsim), rep("s = 1", Nsim), 
            rep("s = 5", Nsim), rep("s = 10", Nsim), 
            rep("s = 20", Nsim), rep("s = 50", Nsim), 
            rep("s = 100", Nsim)
            ),
  value = c(TSmn_H0ls_alloc[,1], TSmn_H01s_alloc[,1], 
            TSmn_H05s_alloc[,1], TSmn_H010s_alloc[,1], 
            TSmn_H020s_alloc[,1], TSmn_H050s_alloc[,1],
            TSmn_H0100s_alloc[,1]
            ))

data_alloc_s$type <- factor(data_alloc_s$type, levels = c("Multinomial-TS (s = 0)", "s = 1",
                                                          "s = 3", "s = 5", 
                                                          "s = 10", "s = 20", 
                                                          "s = 50", "s = 100"
                                                          ), ordered = TRUE)

ggplot(data_alloc_s, aes(x = value, y = type,
               fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,
                      scale = 3, size = 0.3) +
  scale_fill_gradient(low = "white", high = "plum3",
                      name = "Tail prob.") + # To avoid cut off
  theme_minimal() +
  #theme(plot.title = element_text(size=15)) +
  labs(fill="") + ylab(" ") +  xlim(0,1) +
  xlab(expression(T[1]*"/"*T)) +
  theme(legend.position=c(.9,.85),
        legend.text=element_text(size=10)) 

# H0. Sensitivity to prior (Multinomial TS) -------------------

# |-- alpha = 1 ---------------------------------------------------

set.seed(12345)
TSmn_H0sym = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), compute_rho = F)))

#looking at average values
TSmn_H0sym_alloc <- data.frame(matrix(unlist(lapply(TSmn_H0sym, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))


# |-- alpha = 5 ---------------------------------------------------
set.seed(12345)
TSmn_H0sym_alpha5 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                       prior_params = list(arm0 = rep(5,7), arm1 = rep(5,7)), compute_rho = F)))

#looking at average values
TSmn_H0sym_alpha5_alloc <- data.frame(matrix(unlist(lapply(TSmn_H0sym_alpha5, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))
# |-- alpha = 10 ---------------------------------------------------
set.seed(12345)
TSmn_H0sym_alpha10 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                       prior_params = list(arm0 = rep(10,7), arm1 = rep(10,7)), compute_rho = F)))

#looking at average values
TSmn_H0sym_alpha10_alloc <- data.frame(matrix(unlist(lapply(TSmn_H0sym_alpha10, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))

# |-- alpha = 30 ---------------------------------------------------
set.seed(12345)
TSmn_H0sym_alpha30 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                        prior_params = list(arm0 = rep(30,7), arm1 = rep(30,7)), compute_rho = F)))

#looking at average values
TSmn_H0sym_alpha30_alloc <- data.frame(matrix(unlist(lapply(TSmn_H0sym_alpha30, function(x) x$all_prop)), nrow=Nsim, ncol=2, byrow=TRUE))


# |-- Make ridge plot for different priors -------------------

library(ggplot2)
library(ggridges)

data_alloc_s <- data.frame(
  type = c( rep("alpha = 1", Nsim), rep("alpha = 5", Nsim), 
            rep("alpha = 10", Nsim),
            rep("alpha = 30", Nsim)
  ),
  value = c(TSmn_H0sym_alloc[,1], TSmn_H0sym_alpha3_alloc[,1], 
            TSmn_H0sym_alpha5_alloc[,1],
            TSmn_H0sym_alpha10_alloc[,1]
  ))

data_alloc_s$type <- factor(data_alloc_s$type, levels = c("alpha = 1", "alpha = 5",
                                                          "alpha = 10", "alpha = 30"
), ordered = TRUE)

ggplot(data_alloc_s, aes(x = value, y = type,
                         fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE,
                      scale = 3, size = 0.3) +
  scale_fill_gradient(low = "white", high = "lightblue4",
                      name = "Tail prob.") + # To avoid cut off
  theme_minimal() +
  #theme(plot.title = element_text(size=15)) +
  labs(fill="") + ylab(" ") +  xlim(0,1) +
  xlab(expression(T[1]*"/"*T)) +
  theme(legend.position=c(.9,.85),
        legend.text=element_text(size=10)) +
  scale_y_discrete(breaks = c("alpha = 1", "alpha = 5", "alpha = 10", "alpha = 30"), labels = c(bquote(alpha == 1), bquote(alpha == 5), bquote(alpha == 10), bquote(alpha == 30)))
  


# Pseudo-ROC curve ------

Nsim = 1e3; myN=1000

set.seed(12345)
TSmn_H1_mturk_alpha1 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                   prior_params = list(arm0 = rep(1,7), arm1 = rep(1,7)), compute_rho = F)))

set.seed(12345)
TSmn_H1_mturk_alpha3 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                          prior_params = list(arm0 = rep(3,7), arm1 = rep(3,7)), compute_rho = F)))

set.seed(12345)
TSmn_H1_mturk_alpha5 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                          prior_params = list(arm0 = rep(5,7), arm1 = rep(5,7)), compute_rho = F)))

set.seed(12345)
TSmn_H1_mturk_alpha10 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                          prior_params = list(arm0 = rep(10,7), arm1 = rep(10,7)), compute_rho = F)))

set.seed(12345)
TSmn_H1_mturk_alpha20 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                          prior_params = list(arm0 = rep(20,7), arm1 = rep(20,7)), compute_rho = F)))

set.seed(12345)
TSmn_H1_mturk_alpha30 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                          prior_params = list(arm0 = rep(30,7), arm1 = rep(30,7)), compute_rho = F)))

set.seed(12345)
TSmn_H1_mturk_alpha40 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                          prior_params = list(arm0 = rep(40,7), arm1 = rep(40,7)), compute_rho = F)))

set.seed(12345)
TSmn_H1_mturk_alpha50 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                          prior_params = list(arm0 = rep(50,7), arm1 = rep(50,7)), compute_rho = F)))

set.seed(12345)
TSmn_H1_mturk_alpha100 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                            prior_params = list(arm0 = rep(100,7), arm1 = rep(100,7)), compute_rho = F)))

set.seed(12345)
TSmn_H1_mturk_alpha150 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                            prior_params = list(arm0 = rep(150,7), arm1 = rep(150,7)), compute_rho = F)))

set.seed(12345)
TSmn_H1_mturk_alpha200 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                            prior_params = list(arm0 = rep(200,7), arm1 = rep(200,7)), compute_rho = F)))

set.seed(12345)
TSmn_H1_mturk_alpha250 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                            prior_params = list(arm0 = rep(250,7), arm1 = rep(250,7)), compute_rho = F)))

set.seed(12345)
TSmn_H1_mturk_alpha300 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                            prior_params = list(arm0 = rep(300,7), arm1 = rep(300,7)), compute_rho = F)))

set.seed(12345)
TSmn_H1_mturk_alpha400 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                            prior_params = list(arm0 = rep(400,7), arm1 = rep(400,7)), compute_rho = F)))

set.seed(12345)
TSmn_H1_mturk_alpha500 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                            prior_params = list(arm0 = rep(500,7), arm1 = rep(500,7)), compute_rho = F)))

set.seed(12345)
TSmn_H1_mturk_alpha700 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                            prior_params = list(arm0 = rep(700,7), arm1 = rep(700,7)), compute_rho = F)))

set.seed(12345)
TSmn_H1_mturk_alpha1000 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                             prior_params = list(arm0 = rep(1000,7), arm1 = rep(1000,7)), compute_rho = F)))

set.seed(12345)
TSmn_H1_mturk_alpha2000 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                             prior_params = list(arm0 = rep(2000,7), arm1 = rep(2000,7)), compute_rho = F)))

set.seed(12345)
TSmn_H1_mturk_alpha3000 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = MTurk_p0, p1 = MTurk_p1), 
                                                             prior_params = list(arm0 = rep(3000,7), arm1 = rep(3000,7)), compute_rho = F)))


#looking at average values
TSmn_H1_alloc_alpha1 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha1, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H1_alloc_alpha3 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha3, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H1_alloc_alpha5 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha5, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H1_alloc_alpha10 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha10, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H1_alloc_alpha20 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha20, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H1_alloc_alpha30 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha30, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H1_alloc_alpha40 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha40, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H1_alloc_alpha50 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha50, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H1_alloc_alpha100 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha100, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H1_alloc_alpha150 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha150, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H1_alloc_alpha200 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha200, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H1_alloc_alpha250 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha250, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H1_alloc_alpha300 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha300, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H1_alloc_alpha400 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha400, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H1_alloc_alpha500 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha500, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H1_alloc_alpha700 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha700, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H1_alloc_alpha1000 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha1000, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H1_alloc_alpha2000 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha2000, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H1_alloc_alpha3000 <- data.frame(matrix(unlist(lapply(TSmn_H1_mturk_alpha3000, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))


H1_y = c(mean(TSmn_H1_alloc_alpha1[,1]>0.90)*100,
         mean(TSmn_H1_alloc_alpha3[,1]>0.90)*100,
         mean(TSmn_H1_alloc_alpha5[,1]>0.90)*100,
         mean(TSmn_H1_alloc_alpha10[,1]>0.90)*100,
         
         # mean(TSmn_H1_alloc_alpha20[,1]>0.90)*100,
         # mean(TSmn_H1_alloc_alpha30[,1]>0.90)*100,
         # mean(TSmn_H1_alloc_alpha40[,1]>0.90)*100,
         # mean(TSmn_H1_alloc_alpha50[,1]>0.90)*100,
         # mean(TSmn_H1_alloc_alpha100[,1]>0.90)*100,
         
         mean(TSmn_H1_alloc_alpha150[,1]>0.90)*100,
         mean(TSmn_H1_alloc_alpha200[,1]>0.90)*100,
         mean(TSmn_H1_alloc_alpha250[,1]>0.90)*100,
         mean(TSmn_H1_alloc_alpha300[,1]>0.90)*100,
         mean(TSmn_H1_alloc_alpha400[,1]>0.90)*100,
         mean(TSmn_H1_alloc_alpha500[,1]>0.90)*100,
         mean(TSmn_H1_alloc_alpha700[,1]>0.90)*100,
         mean(TSmn_H1_alloc_alpha1000[,1]>0.90)*100,
         mean(TSmn_H1_alloc_alpha2000[,1]>0.90)*100,
         mean(TSmn_H1_alloc_alpha3000[,1]>0.90)*100)
  

# H0  

set.seed(12345)
TSmn_H0_mturk_alpha1 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                          prior_params = list(arm0 = rep(1,7), arm1 = rep(1,7)), compute_rho = F)))

set.seed(12345)
TSmn_H0_mturk_alpha3 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                          prior_params = list(arm0 = rep(3,7), arm1 = rep(3,7)), compute_rho = F)))

set.seed(12345)
TSmn_H0_mturk_alpha5 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                          prior_params = list(arm0 = rep(5,7), arm1 = rep(5,7)), compute_rho = F)))

set.seed(12345)
TSmn_H0_mturk_alpha10 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                           prior_params = list(arm0 = rep(10,7), arm1 = rep(10,7)), compute_rho = F)))

set.seed(12345)
TSmn_H0_mturk_alpha20 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                           prior_params = list(arm0 = rep(20,7), arm1 = rep(20,7)), compute_rho = F)))

set.seed(12345)
TSmn_H0_mturk_alpha30 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                           prior_params = list(arm0 = rep(30,7), arm1 = rep(30,7)), compute_rho = F)))

set.seed(12345)
TSmn_H0_mturk_alpha40 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                           prior_params = list(arm0 = rep(40,7), arm1 = rep(40,7)), compute_rho = F)))

set.seed(12345)
TSmn_H0_mturk_alpha50 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                           prior_params = list(arm0 = rep(50,7), arm1 = rep(50,7)), compute_rho = F)))

set.seed(12345)
TSmn_H0_mturk_alpha100 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                           prior_params = list(arm0 = rep(100,7), arm1 = rep(100,7)), compute_rho = F)))

set.seed(12345)
TSmn_H0_mturk_alpha150 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                           prior_params = list(arm0 = rep(150,7), arm1 = rep(150,7)), compute_rho = F)))

set.seed(12345)
TSmn_H0_mturk_alpha200 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                           prior_params = list(arm0 = rep(200,7), arm1 = rep(200,7)), compute_rho = F)))

set.seed(12345)
TSmn_H0_mturk_alpha250 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                            prior_params = list(arm0 = rep(250,7), arm1 = rep(250,7)), compute_rho = F)))

set.seed(12345)
TSmn_H0_mturk_alpha300 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                            prior_params = list(arm0 = rep(300,7), arm1 = rep(300,7)), compute_rho = F)))

set.seed(12345)
TSmn_H0_mturk_alpha400 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                            prior_params = list(arm0 = rep(400,7), arm1 = rep(400,7)), compute_rho = F)))

set.seed(12345)
TSmn_H0_mturk_alpha500 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                            prior_params = list(arm0 = rep(500,7), arm1 = rep(500,7)), compute_rho = F)))

set.seed(12345)
TSmn_H0_mturk_alpha700 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                            prior_params = list(arm0 = rep(700,7), arm1 = rep(700,7)), compute_rho = F)))

set.seed(12345)
TSmn_H0_mturk_alpha1000 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                            prior_params = list(arm0 = rep(1000,7), arm1 = rep(1000,7)), compute_rho = F)))

set.seed(12345)
TSmn_H0_mturk_alpha2000 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                             prior_params = list(arm0 = rep(2000,7), arm1 = rep(2000,7)), compute_rho = F)))

set.seed(12345)
TSmn_H0_mturk_alpha3000 = pbreplicate(Nsim, list(TS_MultiDir(N=myN, true_params = list(p0 = true_p0_sym, p1 = true_p0_sym), 
                                                             prior_params = list(arm0 = rep(3000,7), arm1 = rep(3000,7)), compute_rho = F)))


#looking at average values
TSmn_H0_alloc_alpha1 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha1, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H0_alloc_alpha3 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha3, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H0_alloc_alpha5 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha5, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H0_alloc_alpha10 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha10, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H0_alloc_alpha20 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha20, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H0_alloc_alpha30 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha30, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H0_alloc_alpha40 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha40, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H0_alloc_alpha50 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha50, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H0_alloc_alpha100 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha100, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H0_alloc_alpha150 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha150, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H0_alloc_alpha200 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha200, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H0_alloc_alpha250 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha250, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H0_alloc_alpha300 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha300, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H0_alloc_alpha400 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha400, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H0_alloc_alpha500 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha500, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H0_alloc_alpha700 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha700, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H0_alloc_alpha1000 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha1000, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H0_alloc_alpha2000 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha2000, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))
TSmn_H0_alloc_alpha3000 <- data.frame(matrix(unlist(lapply(TSmn_H0_mturk_alpha3000, function(x) x$all_prop)), nrow=Nsim*2, ncol=2, byrow=TRUE))


H0_x = c(mean(TSmn_H0_alloc_alpha1[,1]<=0.6 & TSmn_H0_alloc_alpha1[,1]>=0.4)*100,
         mean(TSmn_H0_alloc_alpha3[,1]<=0.6 & TSmn_H0_alloc_alpha3[,1]>=0.4)*100,
         mean(TSmn_H0_alloc_alpha5[,1]<=0.6 & TSmn_H0_alloc_alpha5[,1]>=0.4)*100,
         mean(TSmn_H0_alloc_alpha10[,1]<=0.6 & TSmn_H0_alloc_alpha10[,1]>=0.4)*100,
         
         # mean(TSmn_H0_alloc_alpha20[,1]<=0.6 & TSmn_H0_alloc_alpha20[,1]>=0.4)*100,
         # mean(TSmn_H0_alloc_alpha30[,1]<=0.6 & TSmn_H0_alloc_alpha30[,1]>=0.4)*100,
         # mean(TSmn_H0_alloc_alpha40[,1]<=0.6 & TSmn_H0_alloc_alpha40[,1]>=0.4)*100,
         # mean(TSmn_H0_alloc_alpha50[,1]<=0.6 & TSmn_H0_alloc_alpha50[,1]>=0.4)*100,
         # mean(TSmn_H0_alloc_alpha100[,1]<=0.6 & TSmn_H0_alloc_alpha100[,1]>=0.4)*100,
         
         mean(TSmn_H0_alloc_alpha150[,1]<=0.6 & TSmn_H0_alloc_alpha150[,1]>=0.4)*100,
         mean(TSmn_H0_alloc_alpha200[,1]<=0.6 & TSmn_H0_alloc_alpha200[,1]>=0.4)*100,
         mean(TSmn_H0_alloc_alpha250[,1]<=0.6 & TSmn_H0_alloc_alpha250[,1]>=0.4)*100,
         mean(TSmn_H0_alloc_alpha300[,1]<=0.6 & TSmn_H0_alloc_alpha300[,1]>=0.4)*100,
         mean(TSmn_H0_alloc_alpha400[,1]<=0.6 & TSmn_H0_alloc_alpha400[,1]>=0.4)*100,
         mean(TSmn_H0_alloc_alpha500[,1]<=0.6 & TSmn_H0_alloc_alpha500[,1]>=0.4)*100,
         mean(TSmn_H0_alloc_alpha700[,1]<=0.6 & TSmn_H0_alloc_alpha700[,1]>=0.4)*100,
         mean(TSmn_H0_alloc_alpha1000[,1]<=0.6 & TSmn_H0_alloc_alpha1000[,1]>=0.4)*100,
         mean(TSmn_H0_alloc_alpha2000[,1]<=0.6 & TSmn_H0_alloc_alpha2000[,1]>=0.4)*100,
         mean(TSmn_H0_alloc_alpha3000[,1]<=0.6 & TSmn_H0_alloc_alpha3000[,1]>=0.4)*100)


roc_plot = data.frame(alpha = c(1,3,5,10,
                                # 20,30,40,50,100,
                                150,200, 250, 300,400,500,700,1000, 2000, 3000),
              H0 = H0_x/100, 
              H1 = H1_y/100)

plot(roc_plot$H0, roc_plot$H1, lty=1, ylim=c(0,1), type = "l",
     ylab = expression(P(frac(T[1],T)~">"~0.9~"|"~H[1])),
     xlab = expression(P(0.4~"<"~frac(T[1],T)~"<"~0.6~"|"~H[0])))

ggplot(data = roc_plot, aes(x = roc_plot$H0, y = roc_plot$H1)) +
  geom_line(size=.7) + 
  ylab(expression(P(frac(T[1],T)>=0.9~"|"~H[1]))) + 
  xlab(expression(P(0.4~""<=""~frac(T[1],T)~""<=""~0.6~"|"~H[0]))) +
  geom_hline(yintercept=0.92, linetype="dashed", color = "gray", size=0.4) +
  geom_vline(xintercept=0.8, linetype="dashed", color = "gray", size=0.4) +
  annotate(geom="text", x=0.85, y=0.95, label=expression(alpha~"="~250),
             color="#595959", size = 4.5) +
  annotate(geom="text", x=0.2, y=0.01, label=expression(alpha~"="~1),
           color="#595959", size = 4.5) +
  annotate(geom="text", x=0.4, y=0.01, label=expression(alpha~"="~30),
           color="#595959", size = 4.5) +
  annotate(geom="text", x=0.6, y=0.01, label=expression(alpha~"="~100),
           color="#595959", size = 4.5) +
  annotate(geom="text", x=0.8, y=0.01, label=expression(alpha~"="~250),
           color="#595959", size = 4.5) +
  annotate(geom="text", x=0.99, y=0.01, label=expression(alpha~"="~3000),
           color="#595959", size = 4.5) +
  scale_y_continuous(breaks = c(0.00, 0.25, 0.50, 0.75, 0.92, 1.00), labels = c("0.00", "0.25", "0.50", "0.75", "0.92", "1.00")) +
  scale_x_continuous(breaks = c(0.20, 0.4, 0.6, 0.8, 1.00), labels = c("0.20", "0.40", "0.60", "0.80", "1.00")) +
  theme_bw()
