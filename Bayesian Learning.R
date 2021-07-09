# Bayesian Learning Exam 2021-01-15
# Author: Philip Dahlqvist-Sjöberg

library(latex2exp) # LaTeX in R
library(mvtnorm) # Multinormal
# library(dplyr)
# library(xtable)
# library(EnvStats) # To get truncated normal "rnormTrunc"
library(rstan)

# Clear all in r
rm(list = ls())
cat("\014")

wd = "C:/Users/Philip/Downloads/Bayesian Statistics/Exam/" # Working directory

###############################
########## Problem 1 ########## 
############################### 

#### a) ####

library(MASS);data(shrimp)

# param is the two parameters alpha_tilde and beta_tilde 
log_posterior = function(param, x){
  alpha = param[1]
  beta = param[2]
  prior_a = dnorm(alpha,0,1, log = T)
  prior_b = dnorm(beta,0,1, log = T)
  lik = sum(dbeta(x, exp(alpha), exp(beta), log = T))
  return(lik + prior_a + prior_b)
}

nPara = 2
initVal <- as.vector(rep(0.2,nPara))
X = shrimp/100 # To get percentage

OptimResults<-optim(initVal, log_posterior, gr = NULL, X, method=c("BFGS"), 
                    control = list(fnscale = -1), hessian = TRUE)

postMode = OptimResults$par
postCov = -solve(OptimResults$hessian)

# Print..

postMode
postCov

#### b) ####

nSim = 10000

# simulate alpha_tilde and beta_tilde
normal_approx_draw = rmvnorm(nSim, mean = postMode, sigma = postCov)

x_new = rep(0,nSim)

# simulate x_19 draw
for(i in 1:nSim){
  x_new[i] = rbeta(1, exp(normal_approx_draw[i,1]), exp(normal_approx_draw[i,2]))
}


hist(x_new, freq = F, breaks = 50, main = "Histogram of new x observation", xlab = TeX("$\\x_{19}$"))

#### c) ####

x_new = matrix(0, nSim, 10)

# Simulate 10000 new observations with 10 deliveries in each iteration
for(j in 1:10){
  
  for(i in 1:nSim){
    x_new[i, j] = rbeta(1, exp(normal_approx_draw[i,1]), exp(normal_approx_draw[i,2]))
  }
  
}
head(x_new)

prop_bad = apply(x_new, 1, function(i)sum(i<0.2))
hist(prop_bad, freq = F, breaks = 4, main = "Histogram of predictive distribution", xlab = "Number of bad cocktails")
table(prop_bad)
table(prop_bad)/nSim


###############################
########## Problem 2 ########## 
############################### 

# X adn y in the environment
load("C:/Users/Philip/Downloads/Bayesian Statistics/Exam/Sonar.RData")
# y: 1 = metal cylinder, 0 = rock
names(X[1]) = "Intercept"
#### a) ####

nPara <- ncol(X) # Number of covariates
tau = 10 # Prior standard deviations
mu <- as.vector(rep(0,nPara)) # Prior mean vector
Sigma <- diag(tau^2, nPara) # Prior var-cor matrix

# Setting up log-posterior

log_posterior = function(beta_vec, y, X, mu, Sigma, nPara){
  lin_pred = X %*% beta_vec
  log_lik = sum(lin_pred * y - log(1+exp(lin_pred)))
  log_prior = dmvnorm(beta_vec, rep(0,nPara), Sigma, log = T)
  return(log_lik + log_prior)
}

# Running optim

initVal <- as.vector(rep(0,nPara))

OptimResults<-optim(initVal, log_posterior, gr = NULL, y, X, mu, Sigma, nPara, method=c("BFGS"), 
                    control = list(fnscale = -1), hessian = TRUE)

# Picking the mode and obtaining the covariance matrix from the normal approximation.

postMode = OptimResults$par
postCov = -solve(OptimResults$hessian)

# Std for the beta vector and some naming..

postStd <- sqrt(diag(postCov)) 
names(postMode) <- c("Intercept", colnames(X[,2:nPara]))      
names(postStd) <- c("Intercept", colnames(X[,2:nPara]))

# Print..

postMode
postStd

# Plot marginals of beta

par(mfrow=c(2,3))
for(j in 1:6){
  beta_grid = seq(postMode[j]-4*postStd[j],postMode[j]+4*postStd[j], 0.1)
  plot(beta_grid, sapply(beta_grid, function(i)dnorm(i, postMode[j], postStd[j])), type = "l", ylab = "Density", xlab = TeX("$\\beta$"), 
       main = paste0("Posterior parameter ", j))  
}
par(mfrow=c(1,1))

#### b) ####

# I could not solve this....


# # Laplace approximation 
# p = 1
# beta_var = tau^2
# 
# lap_approx = sum(dnorm(X[,6], mean = postMode[6], sd = postStd[6], log = T)) + dnorm(postMode[6], mean = 0, sd = beta_var, log = T) +
#   (1/2)*log(det(postCov)) + (p/2)*log(2*pi)
# 
# M = exp(lap_approx)
# 
# 
# log_reg(X, postMode)

#### c) ####

# I could not solve this....

###############################
########## Problem 4 ########## 
############################### 

load("C:/Users/Philip/Downloads/Bayesian Statistics/Exam/Amazon.RData")

#### a) ####

stanModelNormal = '
// The input data is a vector y of length N.
data {
  // data
  int<lower=0> N;
  vector[N] y;
  // prior
  real mu0;
  real<lower=0> alpha0;
  real<lower=0> beta0;
  real<lower=0> nu0;
  real<lower=0> sigma20;
  real kappa0;
}

// The parameters in the model
parameters {
  real theta;
  real<lower=0> sigma2;
  real<lower=0> nu;
}

model {
  sigma2 ~ scaled_inv_chi_square(nu0, sqrt(sigma20));
  theta ~ normal(mu0,sqrt(sigma2/kappa0));
  nu ~ gamma(alpha0, beta0);
  for(i in 1:N)
  y[i] ~ student_t(nu, theta, sqrt(sigma2));
}
'

data = list(N=length(y), y = y)
prior = list(mu0 = 0, alpha0 = 4, beta0 = 1, nu0 = 5, sigma20 = 2^2, kappa0 = 1)
fit = stan(model_code = stanModelNormal, data = c(data, prior), iter = 10000)

# print(fit, digits_summary = 2) 
# traceplot(fit, pars = c("theta", "sigma2", "nu"), nrow = 2)  

library(RColorBrewer)
plotColors = brewer.pal(12, "Paired")

# Extract the posterior samples from stan's fit object
postDraws <- extract(fit, permuted = TRUE) # return a list of arrays 
thetaDraws <- postDraws$theta
sigmaDraws <- sqrt(postDraws$sigma2)
nuDraws <- postDraws$nu

# Plot marginals 
par(mfrow=c(2,2))
hist(thetaDraws, 30, main = expression(theta), xlab = "", ylab = "", yaxt='n', 
     col = plotColors[2], border = F)
hist(sigmaDraws, 30, main = expression(sigma), xlab = "", ylab = "", yaxt='n', 
     col = plotColors[2], border = F)
hist(nuDraws, 30, main = expression(nu), xlab = "", ylab = "", yaxt='n', 
     col = plotColors[2], border = F)
par(mfrow=c(1,1))

#### b) ####

nIter = length(thetaDraws)

pred_dist = rep(0,nIter)

for(i in 1:nIter){
  pred_dist[i] = rt(1, df = nuDraws[i])*sigmaDraws[i] + thetaDraws[i]
}

hist(pred_dist, freq = F, breaks = 50, xlab = TeX("$x_{new}$"), main = "Histogram of predictive distribution")

#### c) ####

bet = pred_dist
# Low return
bet[(bet/100)<0.03] = 0
# "High" return
bet[(bet/100)>=0.03] = (bet[(bet/100)>=0.03]/100 + 1)*1000*100

hist(bet)

# How much would I be willing to pay for the bet.
mean(bet)

########
