##############################################
library(tidyverse)
library(gridExtra)
library(HDInterval)
library(latex2exp)
##############################################
set.seed(1995)
##############################################

# Question 1

##############################################

# a)

# Set values.

alpha_0 = beta_0 = 2; s = 14; f = 6; n_draws = 10000

# Draw from the posterior.

posterior_draws = rbeta(n_draws, alpha_0+s, beta_0 + f)

# Compute the theoretical mean and sd.

theoretical_mean = (alpha_0 + s)/(alpha_0 + s + beta_0 + f)
theoretical_sd = sqrt((alpha_0 + s)*(beta_0 + f)/((alpha_0 + s + beta_0 + f)^2*(alpha_0 + s + beta_0 + f+1)))

# Compute the mean and the sd as n increases to check the convergence. Lengthy expression for sd but avoids loops.

mean_convergence = cumsum(posterior_draws)/1:n_draws
sd_convergence = sqrt((cumsum(posterior_draws^2)[2:length(mean_convergence)]-mean_convergence[2:length(mean_convergence)]^2*
                         2:length(mean_convergence))/(1:(length(mean_convergence)-1)))

# Plot the convergence. Looks nice!

mean_convergence = as.data.frame(cbind(mean_convergence, draw = c(1:length(mean_convergence))))
sd_convergence = as.data.frame(cbind(sd_convergence, draw = c(2:(length(sd_convergence)+1))))
p1 = mean_convergence %>% ggplot(aes(1:length(mean_convergence), mean_convergence)) + geom_line(col = "red") +
  geom_line(aes(y = theoretical_mean)) + labs(x = "Number of posterior draws", y = "Posterior Mean") + theme_classic()
p2 = sd_convergence %>% ggplot(aes(draw, sd_convergence)) + geom_line(col = "blue") +
  geom_line(aes(y = theoretical_sd)) + labs(x = "Number of posterior draws", y = "Posterior SD") + theme_classic()
grid.arrange(p1, p2)


##############################################

# b)

# Compute MC-estimate of P(theta <= 0.4)

(mc_estimate_less0.4 = sum(posterior_draws <= 0.4)/n_draws)

# Compare with true prob.

(true_prob = pbeta(0.4, alpha_0 + s, beta_0 + f, lower.tail = T))

##############################################

# c)

# Create the log odds and make into df.

log_odds = log(posterior_draws) - log(1-posterior_draws)

# Plot posterior.

log_odds = as.data.frame(log_odds)
ggplot(data=log_odds, aes(x=log_odds)) + 
  geom_histogram(aes(y=..density..), color = "black", fill = "white" , bins = 20) +
  geom_density(alpha=.2, fill="#FF6666") + theme_classic() + labs(x = "Log-odds", y = "Density")


##############################################

# Question 2

##############################################

# a) 

# Set values.

y = c(14, 25, 45, 25, 30, 33, 19, 50, 34, 67)
n = length(y)
mu = 3.5

# Compute tau^2.

tau_sq = sum((log(y)-mu)^2)/n

# Draw from posterior. Compute the mean and sd from the draws.

sigma_posterior_draws = tau_sq*n/rchisq(n_draws, n)
(sigma_posterior_mean = mean(sigma_posterior_draws))
(sigma_posterior_sd = sd(sigma_posterior_draws))

# Compare with the theoretical values from scal-inv-chisquare. Looks similar.

(theoretical_mean = n*tau_sq/(n-2))
(theoretical_sd = sqrt(2*n^2*tau_sq^2/((n-2)^2*(n-4))))

xtable::xtable(cbind(c(sigma_posterior_mean, sigma_posterior_sd), c(theoretical_mean, theoretical_sd)))

# Plotting simulated posterior distribution.

sigma_posterior_draws_df <- data.frame(sigma_posterior_draws)
ggplot(data = sigma_posterior_draws_df) + 
  geom_density(alpha=.3, fill="#FF6666", aes(x = sigma_posterior_draws)) +
  theme_classic() + labs(x = TeX("$\\sigma^2$"), y = "Density")

##############################################

# b)

# Create the posterior draws of the gini-coefficient according to instructions.

gini_coeff = 2*pnorm(sqrt(as.numeric(sigma_posterior_draws)/2))-1

# Plot the posterior.

ggplot() +
  geom_density(aes(x = gini_coeff), alpha=.2, fill="#FF6666") +
  labs(x="Gini Coefficient", y = "Density") + theme_classic()

# Computed mean and standard deviation for the Gini Coefficient.

mean(gini_coeff)
sd(gini_coeff)

##############################################

# c)

# Equal-tail interval. Sort and cut away the 2.5% from every tail.

equal_tail_interval = c(sort(gini_coeff)[n_draws*0.025], sort(gini_coeff)[n_draws*0.975])

# HPDI using density function.

gini_coeff_density = density(gini_coeff) 

gini_coeff_density = as.data.frame(cbind(Density = gini_coeff_density$y, val = gini_coeff_density$x))

gini_coeff_density = gini_coeff_density[order(-gini_coeff_density$Density),]

gini_coeff_density$s <- cumsum(gini_coeff_density$Density)/sum(gini_coeff_density$Density)

hpdi = gini_coeff_density[gini_coeff_density$s < 0.95,]
hpdi = hpdi[order(hpdi$val),]
hpdi_val = hpdi$val[c(1, nrow(hpdi))]

# Check with reference..

hdi(gini_coeff, 0.95)

# Plot.

ggplot() +
  geom_density(aes(x = gini_coeff), alpha=.2, fill="#FF6666") +
  labs(x="Gini Coefficient", y = "Density") + theme_classic() + 
  geom_vline(xintercept = equal_tail_interval, color = "red", lty = 2) + 
  geom_vline(xintercept = hpdi_val, color = "blue", lty = 2)

##############################################

# Question 3

##############################################

# a)

# Create data.

convert_data = c(-2.44, 2.14, 2.54, 1.83, 2.02, 2.33, -2.79, 2.23, 2.07, 2.02)

# Create the grid.

kappa_grid = seq(0, 7.5, 0.01)

# Create the distribution(von Mises).

von_mises_distr = function(y, mu, kappa){
  exp(kappa*cos(y-mu))/(2*pi*besselI(kappa, 0))
}

# Create posterior distribution to feed the grid, i.e. product of likelihood and prior..

posterior_kappa = function(kappa, mu, lambda, data){
  likelihood = prod(sapply(data, von_mises_distr, mu = mu, kappa = kappa))
  priori = dexp(kappa, lambda)
  posterior = likelihood*priori
}

posterior_distribution = sapply(kappa_grid, posterior_kappa, mu = 2.39, lambda = 1, data = convert_data)

# Normalize posterior to sum to 1 and make into df.

result_posterior_k <- data.frame(k = kappa_grid, posterior = posterior_distribution/sum(posterior_distribution))
sum(result_posterior_k$posterior)
# Plot posterior over varying k.

ggplot(data=result_posterior_k, aes(x=k, y=posterior)) +
  geom_line(size=1) +
  geom_area(alpha=.2,fill="#FF6666") +
  theme_classic() + labs(x = TeX("$\\kappa$"), y = "Density") + 
  geom_vline(xintercept = kappa_grid[which.max(posterior_distribution)], linetype = "dashed")

##############################################

# b)

kappa_grid[which.max(posterior_distribution)]

##############################################

