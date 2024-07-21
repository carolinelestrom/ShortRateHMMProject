####################################################################################################################
####################################################################################################################
#------------------------------------------ Viterbi SSM ------------------------------------------------------------
####################################################################################################################
####################################################################################################################

#------------------------------------------ sigma & theta & kappa ------------------------------------------------------------

### Continuous time state space process for sigma, theta, and kappa



c(plogis(ssmmod_sigma_theta_kappa$estimate[1]), exp(ssmmod_sigma_theta_kappa$estimate[2]), exp(ssmmod_sigma_theta_kappa$estimate[3]), ssmmod_sigma_theta_kappa$estimate[4], exp(ssmmod_sigma_theta_kappa$estimate[5]))


viterbi <- function(x, kappa, theta, sigma, phi_state, sigma_state, m, bm){ 
  n <- length(x)
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
  xi <- matrix(0, nrow = n, ncol = m) 
  foo <- delta * dnorm(x[1], exp(bstar) * theta, sqrt((exp(bstar) * sigma)^2/(2*kappa * exp(bstar))))
  xi[1, ] <- foo / sum(foo)
  for (t in 2:n){
    foo <- apply(xi[t - 1, ] * Gamma, 2, max) * dnorm(x[t], exp(bstar) * theta, sqrt((exp(bstar) * sigma)^2/(2*kappa * exp(bstar))))
    xi[t, ] <- foo / sum(foo) 
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ]) 
  for (t in (n - 1):1){
    iv[t] <- which.max(Gamma[, iv[t + 1]] * xi[t, ]) ### Index of states
  }
  
  s <-  rep(NA, length = n)
  for(t in 1:n){
    s[t] <- bstar[iv[t]]
  }
  return(s)
}



Viterbi_sigma_theta_kappa <- viterbi(x = as.numeric(yields_df$"3M"), 
                                     exp(ssmmod_sigma_theta_kappa$estimate[3]), 
                                     ssmmod_sigma_theta_kappa$estimate[4], 
                                     exp(ssmmod_sigma_theta_kappa$estimate[5]), 
                                     plogis(ssmmod_sigma_theta_kappa$estimate[1]), 
                                     exp(ssmmod_sigma_theta_kappa$estimate[2]), 
                                     m = 200, bm = 3)
plot(exp(Viterbi_sigma_theta_kappa)*ssmmod_sigma_theta_kappa$estimate[4])
points(as.numeric(yields_df$"3M"), col = "steelblue")

plot(Viterbi_sigma_theta_kappa)
