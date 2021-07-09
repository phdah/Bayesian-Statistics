###### Lab 2 ######
library(mvtnorm)
library(latex2exp)


###### 1 ######

rm(list = ls())
cat("\014")

## a)
tempLinkoping <- read.table(url("https://raw.githubusercontent.com/mattiasvillani/BayesLearnCourse/master/Labs/TempLinkoping.txt"), header = T)
Y = tempLinkoping$temp
X = tempLinkoping$time
X = matrix(c(rep(1, length(X)), X, X^2), ncol = 3, dimnames = list(1:length(X), c("X_0", "X_1", "X_2")))

n = nrow(tempLinkoping)
mu_0 = matrix(c(-10, 100, -100), ncol = 1)
omega_0 = 0.01*diag(3)
nu_0 = 4
sig2_0 = 1

n_draws = 10000
# Draw the joint prior

sigma_prior_draws = sig2_0*nu_0/rchisq(n_draws, nu_0)
beta_prior_draws = t(sapply(sigma_prior_draws , function(x)rmvnorm(1, mu_0, x*solve(omega_0))))

# Plot the prior distributions

par(mfrow=c(2,2))
plot(density(sigma_prior_draws)$x, density(sigma_prior_draws)$y, main = TeX("$\\sigma^2$ simulated prior distribution"),ylab = "Frequency", xlab = TeX("$\\sigma^2$"), xlim = c(0,10), type = "h")
legend("topright", legend = paste0("Mode: ", round(density(sigma_prior_draws)$x[which.max(density(sigma_prior_draws)$y)],2)))
hist(beta_prior_draws[,1], breaks = 50, main = TeX("$\\beta_0$ simulated prior distribution"), xlab = TeX("$\\beta_0$"))
legend("topright", legend = paste0("Mean: ", round(mean(beta_prior_draws[,1]),2)))
hist(beta_prior_draws[,2], breaks = 50, main = TeX("$\\beta_1$ simulated prior distribution"), xlab = TeX("$\\beta_1$"))
legend("topright", legend = paste0("Mean: ", round(mean(beta_prior_draws[,2]),2)))
hist(beta_prior_draws[,3], breaks = 50, main = TeX("$\\beta_2$ simulated prior distribution"), xlab = TeX("$\\beta_2$"))
legend("topright", legend = paste0("Mean: ", round(mean(beta_prior_draws[,3]),2)))
par(mfrow=c(1,1))

# The data

X = tempLinkoping$time

plot(X, Y, xlab = "X", main = "Temp over time")

# Polynomial regression curves


f = function(X, B){
  B[1] + B[2]*X + B[3]*X^2
}
f_it = function(X, B_it){
  B_it[,1] + B_it[,2]*X + B_it[,3]*X^2
}

plot_regression = function(X, Y, beta_prior_draws, sig2_0, prior_value){
  # s = paste0(round(mu_0[1]), ", ", round(mu_0[2]), ", ", round(mu_0[3]))
  # prio_s = ifelse(prior_value!=0, paste0("Temp over time Prior : ", s), paste0("Temp over time"))
  plot(X, Y, xlab = "X", main = paste0("Temp over time"), ylim = c(-50, 80))
  B = apply(beta_prior_draws , 2, function(z)density(z)$x[which.max(density(z)$y)])
  reg = f(X,B)
  B_it = beta_prior_draws
  q = apply(B_it, 1, function(B) f(X,B))
  apply(q, 2, function(C)lines(X, C, col = "blue"))
  lines(X, reg, col = "red", lwd = 3)
}
plot_regression(X, Y, beta_prior_draws[1:100,], sig2_0, 0)

# Not so good results, lets change the prior parameters:

priors = function(n_draws, sig2_0, nu_0, mu_0, omega_0){
  sigma_prior_draws = sig2_0*nu_0/rchisq(n_draws, nu_0)
  beta_prior_draws = t(sapply(sigma_prior_draws , function(x)rmvnorm(1, mu_0, x*solve(omega_0))))
  return(beta_prior_draws)
}

n_draws = 100 
mu_0 = matrix(c(-8, 105, -100), ncol = 1)
omega_0 = 0.01*diag(3)
nu_0 = 4
sig2_0 = 0.4

set.seed(1989);beta_prior_draws_opti = priors(n_draws, sig2_0, nu_0, mu_0, omega_0)
plot_regression(X, Y, beta_prior_draws_opti, sig2_0, 1)
cat("Prior: mu - ", mu_0, "\n    omega - ", omega_0, "\n       nu - ", nu_0, "\n    sigma - ", sig2_0)
# Tried to match the regression curves to my belives about the temp in linköping during the different months.

# Seems to be the sig^2 parameter that finds a good fit
# 
# n_draws = 10000
# sig2_0 = seq(0.01, 0.1, by = (0.1-0.01)/(6-1))
# 
# par(mfrow=c(3,3))
# sapply(sig2_0, function(sig2_0)plot_regression(X, Y, priors(n_draws, sig2_0, nu_0, mu_0, omega_0), sig2_0, 1))
# par(mfrow=c(1,1))

## b)
# Calculate the parameters 

X = tempLinkoping$time
X = matrix(c(rep(1, length(X)), X, X^2), ncol = 3, dimnames = list(1:length(X), c("X_0", "X_1", "X_2")))

n_draws = 1000
nu_n = nu_0 + n
B_hat = solve(t(X)%*%X)%*%t(X)%*%Y
mu_n = solve((t(X)%*%X + omega_0)) %*% (t(X)%*%X%*%B_hat + omega_0%*%mu_0)
omega_n = t(X)%*%X + omega_0
sig2_n = (nu_0*sig2_0 + (t(Y)%*%Y + t(mu_0)%*%omega_0%*%mu_0 - t(mu_n)%*%omega_n%*%mu_n)) / nu_n

# Draw the posteriors

sigma_posteriors_draws = sig2_n*nu_n/rchisq(n_draws, nu_n)
beta_posteriors_draws = t(sapply(sigma_posteriors_draws, function(x)rmvnorm(1, mu_n, x*solve(omega_n))))

# Plot the distributions 

par(mfrow=c(2,2))
plot(density(sigma_posteriors_draws)$x, density(sigma_posteriors_draws)$y, main = TeX("$\\sigma^2$ simulated posterior distribution"),ylab = "Frequency", xlab = TeX("$\\sigma^2$"), xlim = c(10,20), type = "h")
legend("topleft", legend = paste0("Mode: ", round(density(sigma_posteriors_draws)$x[which.max(density(sigma_posteriors_draws)$y)],2)))
hist(beta_posteriors_draws[,1], breaks = 50, main = TeX("$\\beta_0$ simulated posterior distribution"), xlab = TeX("$\\beta_0$"))
legend("topright", legend = paste0("Mean: ", round(mean(beta_posteriors_draws[,1]),2)))
hist(beta_posteriors_draws[,2], breaks = 50, main = TeX("$\\beta_1$ simulated posterior distribution"), xlab = TeX("$\\beta_1$"))
legend("topright", legend = paste0("Mean: ", round(mean(beta_posteriors_draws[,2]),2)))
hist(beta_posteriors_draws[,3], breaks = 50, main = TeX("$\\beta_2$ simulated posterior distribution"), xlab = TeX("$\\beta_2$"))
legend("topright", legend = paste0("Mean: ", round(mean(beta_posteriors_draws[,3]),2)))
par(mfrow=c(1,1))

# Plot the regression curves 

plot_regression_post = function(X, Y, beta_posteriors_draws, prior_value){
  # s = paste0(round(mu_n[1]), ", ", round(mu_n[2]), ", ", round(mu_n[3]))
  # paste0("Temp over time\nPrior sigma: ", s)
  plot(X, Y, xlab = "X", main = paste0("Temp over time"), ylim = c(-50, 80))
  B = apply(beta_posteriors_draws , 2, function(z)density(z)$x[which.max(density(z)$y)])
  reg = f(X,B)
  B_it = apply(beta_posteriors_draws, 2, function(x)as.vector(quantile(x, probs = c(0.025, 0.975)))) # Row1 2.5%, Row2 97.5%
  q = apply(B_it, 1, function(B) f(X,B))
  apply(q, 2, function(C)lines(X, C, col = "blue"))
  lines(X, reg, col = "red", lwd = 3)
  cat("Posterior: mu - ", mu_n, "\n    omega - ", omega_n, "\n       nu - ", nu_n, "\n    sigma - ", sig2_n)
}

n_draws = 366
sigma_posteriors_draws = sig2_n*nu_n/rchisq(n_draws, nu_n)
beta_posteriors_draws = t(sapply(sigma_posteriors_draws, function(x)rmvnorm(1, mu_n, x*solve(omega_n))))

X = tempLinkoping$time
plot_regression_post(X, Y, beta_posteriors_draws, 1)
# They don't contain all the points, and it shoulde'nt. 

## c)
# X = tempLinkoping$time

B = apply(beta_posteriors_draws , 2, function(z)density(z)$x[which.max(density(z)$y)])

f_d = function(B){
  X_max_temp = -B[,2]/(2*B[,3])
  return(X_max_temp)
}

x_tilde_all = f_d(beta_posteriors_draws)
hist(x_tilde_all)

x_tilde = X[which.max(f(X,B))]
abline(v = x_tilde, col = "red")

cat("The expected highest temperatur is at time:", round(x_tilde,2))

## d)

# A prior with low impact from the higher order terms, will have a smaller impact if the prior parameter (mu_0)
# for those values that we suspect have a small impact. Also, the orders that we are sure have a small impact,
# will have a small omega, since we are more sure they have a less impact. 

###### 2 ######

rm(list = ls())
cat("\014")

# a)

womenWork = read.table(url("https://raw.githubusercontent.com/mattiasvillani/BayesLearnCourse/master/Labs/WomenWork.dat"), header = T)

glmModel <- glm(Work ~ 0 + ., data = womenWork, family = binomial)

# b)

y = as.vector(womenWork[,1])
X = as.matrix(womenWork[,-1])
covNames = colnames(womenWork[,-1])                # ... and their names
nPara = dim(X)[2]

tau = 10

# Setting up the prior
mu = as.vector(rep(0,nPara)) # Prior mean vector
Sigma = tau^2*diag(nPara) # Prior variance 

# The log posterior function

LogPostLogistic <- function(betaVect,y,X,mu,Sigma){
  nPara <- length(betaVect)
  linPred <- X%*%betaVect
  logLik <- sum( linPred*y -log(1 + exp(linPred)))
  logPrior <- dmvnorm(betaVect, mu, Sigma, log=TRUE)
  return(logLik + logPrior)
}

initVal <- as.vector(rep(0,nPara)); 
OptimResults<-optim(initVal,LogPostLogistic,gr=NULL,y ,X,mu,Sigma,
                    method=c("BFGS"), control=list(fnscale=-1),hessian=TRUE)
postMode = OptimResults$par # Best set of parameter found
postCov = -solve(OptimResults$hessian) # inv(J) - Approx posterior normal covariance matrix
postStd <- sqrt(diag(postCov)) # Computing approximate stdev
names(postMode) <- covNames      # Naming the coefficient by covariates
names(postStd) <- covNames # Naming the coefficient by covariates

postMode
postStd

Small_child = c(postMode[7] - 1.96*postStd[7], postMode[7] + 1.96*postStd[7]) 
names(Small_child) = c("Lower 2.5%", "Upper 97.5 %")
Small_child
# Number of children age 6 or lower is the most important feature. The coefficient has the largest absolut value,
# with a relatively low standard deviation.

# c)

# a 40 year old woman, with two children (3 and 9 years old), 8 years of education, 
# 10 years of experience, and a husband with an income of 10

X_star = matrix(c(1, 10, 8, 10, (10/10)^2, 40, 1, 1), ncol = 8)

n_draws = 5000

pred_sim = function(n_draws, postMode, postCov, X_star){
  beta_draw = rmvnorm(n_draws, postMode, postCov)
  lin_pred = apply(beta_draw, 1, function(x) tcrossprod(X_star, x))
  prob_work = sapply(lin_pred,   function(x) exp(x)/(1+exp(x)))
  pred_work = sapply(prob_work,  function(x) rbinom(1, 1, x))
  return(list(prob_work, pred_work))
}

X_star_dist = pred_sim(n_draws, postMode, postCov, X_star)

hist(X_star_dist[[1]], xlab = "P(work|x)", main = "Posterior distribution", breaks = 40, col = "blue")
legend("topright", legend = paste0("Prop. work: ", round(sum(X_star_dist[[2]])/n_draws,2)))

