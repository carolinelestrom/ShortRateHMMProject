####################################################################################################################
####################################################################################################################
#------------------------------------------ Fitting SSM ------------------------------------------------------------
####################################################################################################################
####################################################################################################################





#------------------------------------------ sigma & theta & kappa ------------------------------------------------------------

### Continuous time state space process for sigma, theta, and kappa

mllk <- function(theta.star, x, m, bm){
  ### S_t model formulation
  phi_state <- plogis(theta.star[1]) ### Mean reversion parameter of state process
  sigma_state <- exp(theta.star[2]) ### Vol/noise parameter of state process
  
  ### r_t model formulation
  kappa <- exp(theta.star[3]) ### r_t \sim N(theta, sigma^2/(2*kappa))
  theta <- theta.star[4] ### r_t \sim N(theta, sigma^2/(2*kappa))
  sigma <- exp(theta.star[5]) ### r_t \sim N(theta, sigma^2/(2*kappa))
  
  ### Intervals to approximate continuous state space
  b <- seq(-bm, bm, length = m + 1) ### Specify boundaries of m intervals
  h <- b[2] - b[1] ### h is the length of each interval
  bstar <- (b[-1] + b[-(m + 1)]) * 0.5 ### Midpoints of the m intervals
  
  ### Set up transition probability matrix
  Gamma <- matrix(0, m, m)
  for (i in 1:m){
    Gamma[i, ] <- h * dnorm(bstar, phi_state * bstar[i], sigma_state) ### m*m t.p.m. of the approx. HMM
  }
  
  ### Set up stationary distribution
  delta <- h * dnorm(bstar, 0, sigma_state / sqrt(1 - phi_state^2)) ### stat. initial distribution
  
  ### Forward algorithm
  foo <- delta * dnorm(x[1], exp(bstar) * theta, sqrt((exp(bstar) * sigma)^2/(2*kappa * exp(bstar)))) ### r_t \sim N(theta exp(theta_t), sigma^2/(2*kappa))
  l <- log(sum(foo))
  phi <- foo / sum(foo)
  for (t in 2:length(x)){
    foo <- phi %*% Gamma * dnorm(x[t], exp(bstar) * theta, sqrt((exp(bstar) * sigma)^2/(2*kappa * exp(bstar))))
    l <- l + log(sum(foo))
    phi <- foo / sum(foo)
  }
  return(-l)
}


theta.star <- c(qlogis(0.999), log(0.2), log(0.001), 0.05, log(0.0001))
# ssmmod_sigma_theta_kappa <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), 
#                                 m = 200, bm = 3, print.level = 2)
ssmmod_sigma_theta_kappa <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), 
                                m = 200, bm = 3, print.level = 2)

c(plogis(ssmmod_sigma_theta_kappa$estimate[1]), exp(ssmmod_sigma_theta_kappa$estimate[2]), exp(ssmmod_sigma_theta_kappa$estimate[3]), ssmmod_sigma_theta_kappa$estimate[4], exp(ssmmod_sigma_theta_kappa$estimate[5]))




