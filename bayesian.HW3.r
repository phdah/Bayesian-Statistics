#######################
library(mvtnorm)
library(latex2exp)
library(rstan)
library(xtable)

rm(list = ls())
cat("\014")
########## 1 ##########
# a)
rainfall = read.table(url("https://raw.githubusercontent.com/mattiasvillani/BayesLearnCourse/master/Labs/rainfall.dat"), header = F, col.names = "Rain")

n = nrow(rainfall)
mu_0 = 50
sig2_0 = 20
nu_0 = 1
tau_0 = 100
x_bar = mean(rainfall$Rain)
x = rainfall$Rain

n_draws = 10000

output = matrix(rep(NA, 2*n_draws), ncol = 2, dimnames = list(rep(1:n_draws), c("Mu", "Sigma")))

# Gibbs sampler

nu_n = nu_0 + n

for(i in 1:n_draws){

  s = (nu_0*sig2_0+sum((x-mu_0)^2))/(n+nu_0)
  
  sig2 = s*nu_n/rchisq(1, nu_n)

  w = (n/sig2)/(n/sig2 + 1/tau_0)
  mu_n = w*x_bar + (1-w)*mu_0
    
  tau_n = 1/(n/sig2 + 1/tau_0)
  
  mu = rnorm(1, mu_n, sqrt(tau_n))
  
  output[i, 1] = mu
  output[i, 2] = sig2
  
  mu_0 = mu
  sig2_0 = sig2
  
}
par(mfrow=c(2,1))
plot(1:n_draws, output[,1], type = "l", col = "blue", ylim = c(30,34), ylab = TeX("$\\mu$"), xlab = "Number of simulations", main = "Gibbs sample\nNormal Model")
abline(h=mean(output[,1]), col = "red")
plot(1:n_draws, output[,2], type = "l", col = "blue", ylim = c(1400, 1800), ylab = TeX("$\\sigma^2$"), xlab = "Number of simulations")
abline(h=mean(output[,2]), col = "red")
par(mfrow=c(1,1))
tail(output)
cat("Mean of mu: ", mean(output[,1]), "\nMean of variance: ", mean(output[,2]))

# b) 

rm(list = ls())
cat("\014")


# Defining a function that simulates from a Dirichlet distribution
rDirichlet <- function(param){
  nCat <- length(param)
  piDraws <- matrix(NA,nCat,1)
  for (j in 1:nCat){
    piDraws[j] <- rgamma(1,param[j],1)
  }
  piDraws = piDraws/sum(piDraws) # Diving every column of piDraws by the sum of the elements in that column.
  return(piDraws)
}

# Simple function that converts between two different representations of the mixture allocation
S2alloc <- function(S){
  n <- dim(S)[1]
  alloc <- rep(0,n)
  for (i in 1:n){
    alloc[i] <- which(S[i,] == 1)
  }
  return(alloc)
}


# Input 

rainfall = read.table(url("https://raw.githubusercontent.com/mattiasvillani/BayesLearnCourse/master/Labs/rainfall.dat"), header = F, col.names = "Rain")
x = as.matrix(rainfall$Rain)

# Model options
nComp <- 2    # Number of mixture components

# Prior options
alpha <- rep(1,nComp) # Dirichlet(alpha)
muPrior <- rep(50,nComp) # Prior mean of mu
tau2Prior <- rep(100,nComp) # Prior variance of mu
sigma2_0 <- rep(var(x),nComp) # s2_0 (best guess of sigma2)
nu0 <- rep(1,nComp) # degrees of freedom for prior on sigma2

# MCMC options
nIter <- 1000 # Number of Gibbs sampling draws

# Plotting options
plotFit <- T
lineColors <- c("blue", "green", "magenta", 'yellow')
sleepTime <- 0.1 # Adding sleep time between iterations for plotting

###### Defining a function that simulates from the 
rScaledInvChi2 <- function(n, df, scale){
  return((df*scale)/rchisq(n,df=df))
}

####### Defining a function that simulates from a Dirichlet distribution
rDirichlet <- function(param){
  nCat <- length(param)
  piDraws <- matrix(NA,nCat,1)
  for (j in 1:nCat){
    piDraws[j] <- rgamma(1,param[j],1)
  }
  piDraws = piDraws/sum(piDraws) # Diving every column of piDraws by the sum of the elements in that column.
  return(piDraws)
}

# Simple function that converts between two different representations of the mixture allocation
S2alloc <- function(S){
  n <- dim(S)[1]
  alloc <- rep(0,n)
  for (i in 1:n){
    alloc[i] <- which(S[i,] == 1)
  }
  return(alloc)
}

# Initial value for the MCMC
nObs <- length(x)
S <- t(rmultinom(nObs, size = 1 , prob = rep(1/nComp,nComp))) # nObs-by-nComp matrix with component allocations.
mu <- quantile(x, probs = seq(0,1,length = nComp))
sigma2 <- rep(var(x),nComp)
probObsInComp <- rep(NA, nComp)

# Setting up the plot
xGrid <- seq(min(x)-1*apply(x,2,sd),max(x)+1*apply(x,2,sd),length = 100)
xGridMin <- min(xGrid)
xGridMax <- max(xGrid)
mixDensMean <- rep(0,length(xGrid))
effIterCount <- 0
ylim <- c(0,2*max(hist(x, plot = F)$density))

output = matrix(rep(NA, 6*nIter), ncol = 6, dimnames = list(1:nIter, c("Mu_1", "Mu_2", "Sig2_1", "Sig2_2", "Pi", "1-Pi")))

for (k in 1:nIter){
  message(paste('Iteration number:',k))
  alloc <- S2alloc(S) # Just a function that converts between different representations of the group allocations
  nAlloc <- colSums(S)
  print(nAlloc)
  # Update components probabilities
  pi <- rDirichlet(alpha + nAlloc)
  output[k, 5] = pi[1]
  output[k, 6] = pi[2]
  
  # Update mu's
  for (j in 1:nComp){
    precPrior <- 1/tau2Prior[j]
    precData <- nAlloc[j]/sigma2[j]
    precPost <- precPrior + precData
    wPrior <- precPrior/precPost
    muPost <- wPrior*muPrior + (1-wPrior)*mean(x[alloc == j])
    tau2Post <- 1/precPost
    mu[j] <- rnorm(1, mean = muPost, sd = sqrt(tau2Post))
    output[k,j] = mu[j]
  }
  
  # Update sigma2's
  for (j in 1:nComp){
    sigma2[j] <- rScaledInvChi2(1, df = nu0[j] + nAlloc[j], scale = (nu0[j]*sigma2_0[j] + sum((x[alloc == j] - mu[j])^2))/(nu0[j] + nAlloc[j]))
    output[k,(j+2)] = sigma2[j]
  }
  
  # Update allocation
  for (i in 1:nObs){
    for (j in 1:nComp){
      probObsInComp[j] <- pi[j]*dnorm(x[i], mean = mu[j], sd = sqrt(sigma2[j]))
    }
    S[i,] <- t(rmultinom(1, size = 1 , prob = probObsInComp/sum(probObsInComp)))
  }
  
  # Printing the fitted density against data histogram
  if (plotFit && (k%%1 ==0)){
    effIterCount <- effIterCount + 1
    # hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = paste("Iteration number",k), ylim = ylim)
    mixDens <- rep(0,length(xGrid))
    components <- c()
    for (j in 1:nComp){
      compDens <- dnorm(xGrid,mu[j],sd = sqrt(sigma2[j]))
      mixDens <- mixDens + pi[j]*compDens
      # lines(xGrid, compDens, type = "l", lwd = 2, col = lineColors[j])
      components[j] <- paste("Component ",j)
    }
    mixDensMean <- ((effIterCount-1)*mixDensMean + mixDens)/effIterCount
    
    # lines(xGrid, mixDens, type = "l", lty = 2, lwd = 3, col = 'red')
    # legend("topright", box.lty = 1, legend = c("Data histogram",components, 'Mixture'),
    #        col = c("black",lineColors[1:nComp], 'red'), lwd = 2)
    # Sys.sleep(sleepTime)
  }
  
}

par(mfrow=c(1,2))
# Plot mu
plot(1:nIter, output[,1], type = "l", col = "blue", ylab = TeX("$\\mu_1$"), xlab = "Iteration", main = "Gibbs sample\nMix Normal Model")
abline(h=mean(output[,1]), col = "red")
legend("topright", legend = TeX("$\\mu_1$ "), col = "blue", lty = 1)

plot(1:nIter, output[,2], col = "orange", type = "l", ylab = TeX("$\\mu_2$"), xlab = "Iteration", main = "Gibbs sample\nMix Normal Model")
abline(h=mean(output[,2]), col = "red")
legend("topright", legend = TeX("$\\mu_2$ "), col = "orange", lty = 1)

# Plot sigma
plot(1:nIter, output[,3], type = "l", col = "blue", ylab = TeX("$\\sigma^2_1$"), xlab = "Iteration", main = "Gibbs sample\nMix Normal Model")
abline(h=mean(output[,3]), col = "red")
legend("topright", legend = TeX("$\\sigma^2_1$ "), col = "blue", lty = 1)

plot(1:nIter, output[,4], col = "orange", ylim = c(0,300), type = "l", ylab = TeX("$\\sigma^2_2$"), xlab = "Iteration", main = "Gibbs sample\nMix Normal Model")
abline(h=mean(output[,4]), col = "red")
legend("topright", legend = TeX("$\\sigma^2_1$ "), col = "orange", lty = 1)

par(mfrow=c(1,1))

# Plot pi
plot(1:nIter, output[,5], type = "l", ylim = c(0, 1), col = "blue", ylab = TeX("$\\pi$"), xlab = "Iteration", main = "Gibbs sample\nMix Normal Model")
abline(h=mean(output[,5]), col = "red")
legend("topright", legend = "Pi", col = "blue", lty = 1)

tail(output)

# c)

hist(x, breaks = 20, freq = FALSE, xlim = c(xGridMin,xGridMax), main = "Final fitted density")
lines(xGrid, mixDensMean, type = "l", lwd = 2, lty = 4, col = "red")
lines(xGrid, dnorm(xGrid, mean = 32.283, sd = sqrt(1547.349)), type = "l", lwd = 2, col = "blue")
legend("topright", box.lty = 1, legend = c("Data histogram","Mixture density","Normal density"), col=c("black","red","blue"), lwd = 2)


########## 2 ##########

rm(list = ls())
cat("\014")

# a)

eBayBidderData = read.table(url("https://raw.githubusercontent.com/mattiasvillani/BayesLearnCourse/master/Labs/eBayNumberOfBidderData.dat"), header = T)

head(eBayBidderData)

eBayBidderData_glm = eBayBidderData[,-2]

fit = glm(nBids ~.,data = eBayBidderData_glm, family = poisson(link=log))
summary(fit)

cat("Significant covariates are:\nIntercetp\nVerifyID\nSealed\nMajBlem\nLogBook\nMinBidShare")

# b)

y = as.matrix(eBayBidderData[,1])
X = as.matrix(eBayBidderData[,-1])
covNames = colnames(eBayBidderData[,-1])                # ... and their names
nPara = dim(X)[2]

# Setting up the prior
mu = as.vector(rep(0,nPara)) # Prior mean vector
Sigma = 100*solve(t(X)%*%X) # Prior variance 

# The log posterior function

log_posterior = function(beta_vec, y, X){
  n_para = ncol(X)
  beta_vec = matrix(beta_vec, nrow = n_para)
  log_prior = dmvnorm(as.vector(beta_vec), rep(0, n_para), 100*chol2inv(chol(crossprod(X))), log = T)
  lin_pred = X %*% beta_vec
  log_lik = sum(y * lin_pred - exp(lin_pred))
  return(log_prior + log_lik)
}

initVal <- as.vector(rep(0,nPara)); 
OptimResults<-optim(initVal,log_posterior,gr=NULL,y ,X,
                    method=c("BFGS"), control=list(fnscale=-1),hessian=TRUE)
postMode = OptimResults$par # Best set of parameter found
postCov = -solve(OptimResults$hessian) # inv(J) - Approx posterior normal covariance matrix
postStd <- sqrt(diag(postCov)) # Computing approximate stdev
names(postMode) <- covNames      # Naming the coefficient by covariates
names(postStd) <- covNames # Naming the coefficient by covariates

postMode
postStd

postTable = rbind(postMode, postStd)

xtable(postTable)

# Metropolis algorithm
mh_symmetric <- function(m = 10000, burn_in = 0, init_beta = postMode, init_Sigma = postCov, s = 0.65){
  # Create empty output list 
  output <- list(Chain = NA, Final = NA, Optim = NA, Accept_rate = NA, Iterations = NA)
  
  # Gibbs sampler
  burn <- burn_in # Possiblility to remove the first "r" observations in the chain
  e <- burn + m # Number of iterations, inncluding the burn in 
  
  # Initialize the chain. 
  beta = rep(0, length(init_beta))
  
  # Store the samples
  chain_beta <- matrix(rep(NA, (e-burn)*length(init_beta)), ncol = length(init_beta))
  
  accepted_number <- 0 # Initial value
  for(i in 1:e){
    
    # Draw a value from the proposal distribution.
    beta_candidate <- as.vector(rmvnorm(1, beta, s*init_Sigma))
   
    
    # Compute the acceptance probability, Equation 8 and Equation 6. 
    # We will do both equations in log domain here to avoid underflow. 
    accept <- log_posterior(beta_candidate, y, X) - log_posterior(beta, y, X)
    accept <- exp(accept)
    
    # LogPostPoisson <- function(betaVect,y,X,mu,Sigma){
    
    # Handle possible NaN 
    if(is.nan(accept)){
      accept <- 0
    }else{
      accept <- min(1, accept)
    }
    
    # Accept rho_candidate with probability accept.
    if (runif(1)<accept){
      beta <- beta_candidate
      accepted_number <- accepted_number+1
    }else{
      beta <- beta
    }
    # store without burn in 
    if(i >= burn){
      chain_beta[i-(burn),] <- beta
    }
  }
  output[[1]] = chain_beta # Store beta chain
  output[[2]] = colMeans(output[[1]]) # Final value
  names(output[[2]]) = names(init_beta) # Colnames final value
  output[[3]] = init_beta # Store optim_beta
  output[[4]] = accepted_number / (e + burn) # Store number of accepted candidates
  output[[5]] = e - burn # Store number of iterations
  return(output)
  
}
n_iterations = 10000
output = mh_symmetric(n_iterations, s = 0.76)

output[4] # Acceptance 0.234 optimal for random walk with more parameters
output[2]
output[3]

beta_chain = unlist(output[[1]])

par(mfrow=c(3,3))
for(i in 1:ncol(beta_chain)){
  b = paste0("$\\beta_", i, "$")
  plot(1:n_iterations, beta_chain[,i], type = "l", col = "blue", xlab = "", ylab = "", main = TeX(b))
}
par(mfrow=c(1,1))

########## 3 ##########

rm(list = ls())
cat("\014")

# a)
mu = 10; sig2 = 2; t = 200; phi = seq(-0.95,0.95,length.out = 5)

epsilon = rnorm(t-1, 0, sig2)

AR = function(mu, t, phi){
  
  x = matrix(rep(NA, t), ncol = 1)
  x[1] = mu
  for(i in 2:t){
    x[i] = mu + phi*(x[i-1]-mu) + epsilon[i-1]
  }
  return(x)
}

Iteration = 1:t

x_t = sapply(phi, function(phi)AR(mu, t, phi))
par(mfrow=c(3,2))
for(i in 1:ncol(x_t)){
  m = paste0("$\\phi$ = ", phi[i])
  plot(Iteration, x_t[,i], type = "l", ylim = c(-10,30), main = TeX(m), ylab = TeX("$\\x_t$"))
}
par(mfrow=c(1,1))
# phi represents the autocorrelation term. For large (+) phi, the previous observation will yield
# another observation in the same direction. Small (-) phi will give another observation in the 
# oposite direction, while phi close to 0 only be effcted by the epsilon term...

# b)

stanModel = "
data{
int<lower=0> N; // Number of datapoints
vector[N] y; // Data vector
}
parameters{
real mu;
real <lower=-0.999, upper=0.999> phi;
real <lower=0> sigma_sq;
}
model{
mu ~ normal(0, 100);
phi ~ uniform(-0.999, 0.999);
sigma_sq ~ scaled_inv_chi_square(1,1);
for(i in 2:N) 
y[i] ~ normal(mu + phi*(y[i-1]-mu), sqrt(sigma_sq));
}"

dataNormHiarch <- list(y = as.vector(AR(mu, t, 0.3)), N = t)

fit1<-stan(model_code=stanModel,
           data=dataNormHiarch,
           warmup=1000,
           iter=5000,
           chains=4)

print(fit1,digits_summary=3)

traceplot(fit1, pars = c("mu", "sigma_sq", "phi"), nrow = 3)

stan_hist(fit1, pars = c("mu", "phi"), nrow = 3, bins = 50)

xtable(summary(fit1)$summary)

# res = extract(fit1, permuted = F)
# res[, chain , parameter] <- parameters are ordered as stated in the model

# c)

campy = read.table(url("https://raw.githubusercontent.com/mattiasvillani/BayesLearnCourse/master/Labs/campy.dat"), header = T)

stanModel = "
data{
// data
int<lower=0> N;
int<lower=0> y[N];
// prior
real mu_0;
real<lower=0> sigma_sq_mu_0;
real a;
real<lower=a> b;
int<lower=0> nu_0;
real<lower=0> sigma_sq_0;
}
parameters{
real mu;
real <lower=-0.999, upper=0.999> phi;
real <lower=0> sigma_sq;
real x[N];
}
model{
mu ~ normal(mu_0, sigma_sq_mu_0);
sigma_sq ~ scaled_inv_chi_square(nu_0, sigma_sq_0);
phi ~ uniform(a, b);
for(i in 2:N) x[i] ~ normal(mu + phi*(x[i-1]-mu), sqrt(sigma_sq));
for(j in 1:N)
y[j] ~ poisson(exp(x[j]));
}"

data = list(N = nrow(campy), y = as.vector(campy$c))
prior = list(mu_0 = 0, sigma_sq_mu_0 = 50, a = -0.999, b = 0.999, nu_0 = 1, sigma_sq_0 = 5)

fit_poi = stan(model_code = stanModel, 
               data = c(data, prior),
               warmup=1000,
               iter=5000,
               chains=4)
print(fit_poi)

postDraws = extract(fit_poi)

xgrid = 1:nrow(campy)

plot(campy$c, ylab = TeX("$\\theta$"))
lines(xgrid, sapply(xgrid, function(i)exp(mean(postDraws$x[,i]))), type = "l", col = "darkorange")
lines(xgrid, sapply(xgrid, function(i)exp(quantile(postDraws$x[,i], probs = 0.025))), type = "l", col = "lightblue")
lines(xgrid, sapply(xgrid, function(i)exp(quantile(postDraws$x[,i], probs = 0.975))), type = "l", col = "lightblue")

# d)

data = list(N = nrow(campy), y = as.vector(campy$c))

# Updated prior belives about the vairance, a more informative prior model
prior = list(mu_0 = 0, sigma_sq_mu_0 = 50, a = -0.999, b = 0.999, nu_0 = 100, sigma_sq_0 = 0.1)

fit_poi = stan(model_code = stanModel, 
               data = c(data, prior),
               warmup=1000,
               iter=5000,
               chains=4)
print(fit_poi)

postDraws = extract(fit_poi)

xgrid = 1:nrow(campy)

plot(campy$c, ylab = TeX("$\\theta$"))
lines(xgrid, sapply(xgrid, function(i)exp(mean(postDraws$x[,i]))), type = "l", col = "darkorange")
lines(xgrid, sapply(xgrid, function(i)exp(quantile(postDraws$x[,i], probs = 0.025))), type = "l", col = "lightblue")
lines(xgrid, sapply(xgrid, function(i)exp(quantile(postDraws$x[,i], probs = 0.975))), type = "l", col = "lightblue")

# Smoother credible interval.

############






