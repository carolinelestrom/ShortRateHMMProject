####################################################################################################################
####################################################################################################################
#---------------------------------------- Fitting Vasicek ----------------------------------------------------------
####################################################################################################################
####################################################################################################################


### Define the log-likelihood function for the Vasicek model
vasicek_loglik <- function(params, r) {
  a <- exp(params[1])
  b <- params[2]
  sigma <- exp(params[3])
  
  n <- length(r)
  dt <- 1  ### assuming unit time steps
  
  ### Stationary distribution
  exp_r <- b
  var_r <- sigma^2/(2*a)
  
  ### Calculate the log-likelihood
  ll <- 0
  for (i in 2:n) {
    #expected_r <- r[i-1] * exp(-a * dt) + b * (1 - exp(-a * dt))
    #var_r <- sigma^2 * (1 - exp(-2 * a * dt)) / (2 * a)
    ll <- ll + dnorm(r[i], mean = exp_r, sd = sqrt(var_r), log = TRUE)
  }
  
  return(-ll)  ### return negative log-likelihood for minimization
}



### Initial parameter estimates
init_params <- c(log(0.0001), 0.07, log(0.0005))

### Prediction data ~ First 9000 observations
PredDat <- as.numeric(yields_df$"3M")[1:6030]
PredDat <- rate_df$rate[1:6989] ### 7052 - 3 * 21
PredDat <- BondDat_df$Rate[1:5870]

### Use optimizer to fit the Vasicek model
mod <- nlm(f = vasicek_loglik, p = init_params, r = as.numeric(yields_df$"3M"), print.level = 2)


### Extract fitted parameters
a_hat <- exp(mod$estimate[1])
b_hat <- mod$estimate[2]
sigma_hat <- exp(mod$estimate[3])

### Print the fitted parameters
cat("Fitted parameters:\n")
cat("a:", a_hat, "\n") ### 0.0004137868 
cat("b:", b_hat, "\n") ### 0.0160007 
cat("sigma:", sigma_hat, "\n") ### 0.0005440518 

