#####################################################################
#####  code for Simulation 1
library(rstan)
library(MASS)
library(msm) 
source("functions.R")

#####################################################################
n_student = 1000  
n_item = 20  
sigma_j=0.5  
rho=0.3 
######  generate the real parameters
a <- rtnorm(n_item, mean = 1.5, sd = 0.5, lower = 0, upper = Inf) 
b <- rnorm(n_item,0,1)
alpha1 <- rtnorm(n_item, mean = 0, sd = 0.5, lower = -1, upper = Inf) 
alpha2 <- rtnorm(n_item, mean = 0, sd = 0.5, lower = -1, upper = Inf) 
lambda <- rnorm(n_item,0.25,sqrt(0.2))
mu_theta_tao <- c(0,0)
rho_theta_tao <- 0.5*rho 
Sigma_theta_tao <- matrix(rho_theta_tao,2,2)
Sigma_theta_tao[1,1] <- 1
Sigma_theta_tao[2,2] <- 0.25
theta_tao <- mvrnorm(n_student,mu_theta_tao,Sigma_theta_tao)
theta <- theta_tao[,1]
tao <- theta_tao[,2]
#####################################################################
#####  main code
Y <- Generate_Y_2PGlogit(a,b,alpha1,alpha2,theta)###generate reponse
logT <- Generate_logT(tao,lambda,0.5)###generate response time
rstan_options(auto_write = TRUE)
data_irt = list(n_student = n_student,n_item = n_item, Y = Y, logT=logT)
fit = stan(file = "2PGlogit_RT.stan", data = data_irt, seed = 1, warmup = 1000, iter = 3000, chains = 4)###call the stan code
#####estimate the parameter from the stan result
post_mean_a = get_posterior_mean(fit, pars=c("a"))
post_mean_b = get_posterior_mean(fit, pars=c("b"))
post_mean_alpha1 = get_posterior_mean(fit, pars=c("alpha1"))
post_mean_alpha2 = get_posterior_mean(fit, pars=c("alpha2"))
post_mean_lambda = get_posterior_mean(fit, pars=c("lambda"))
post_mean_theta = get_posterior_mean(fit, pars=c("theta"))
post_mean_tao = get_posterior_mean(fit, pars=c("tao"))  
post_mean_sigma_tao = get_posterior_mean(fit, pars=c("sigma_tao"))
post_mean_sigma_tj = get_posterior_mean(fit, pars=c("sigma_tj"))
post_mean_Omega_person = get_posterior_mean(fit, pars=c("Omega_person"))
post_mean_Sigma_person = get_posterior_mean(fit, pars=c("Sigma_person"))