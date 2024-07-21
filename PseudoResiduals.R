####################################################################################################################
####################################################################################################################
#----------------------------------------- Pseudo Residuals --------------------------------------------------------
####################################################################################################################
####################################################################################################################



#----------------------------------------- Continuous SSM --------------------------------------------------------

extract_parameters <- function(theta.star) {
  a <- exp(theta.star[3])
  b <- theta.star[4]
  sigma <- exp(theta.star[5])
  phi_state <- plogis(ssmmod_sigma_theta_kappa$estimate[1])
  sigma_state <- exp(ssmmod_sigma_theta_kappa$estimate[2])
  
  list(a = a, b = b, sigma = sigma, phi_state = phi_state, sigma_state = sigma_state)
}




## Compute the log-forward probabilities
lForward <- function(y, mod, m, bm){
  params <- extract_parameters(mod$estimate)  # Extract parameters from mod$estimate
  T <- length(y)
  
  ### Intervals to approximate continuous state space
  b <- seq(-bm, bm, length = m + 1) ### Specify boundaries of m intervals
  h <- b[2] - b[1] ### h is the length of each interval
  bstar <- (b[-1] + b[-(m + 1)]) * 0.5 ### Midpoints of the m intervals
  
  ### Set up transition probability matrix
  Gamma <- matrix(0, m, m)
  for (i in 1:m){
    Gamma[i, ] <- h * dnorm(bstar, params$phi_state * bstar[i], params$sigma_state) ### m*m t.p.m. of the approx. HMM
  }
  
  ### Set up stationary distribution
  delta <- h * dnorm(bstar, 0, params$sigma_state / sqrt(1 - params$phi_state^2)) ### stat. initial distribution
  
  mus <- matrix(exp(bstar) * params$b, nrow = T, ncol = m, byrow = TRUE)
  sigmas <- matrix(sqrt((exp(bstar) * params$sigma)^2 / (2 * params$a * exp(bstar))), nrow = T, ncol = m, byrow = TRUE)
  
  lalpha <- matrix(NA, m, T)
  P <- dnorm(y[1], mean = mus[1, ], sd = sigmas[1, ])
  foo <- delta * P
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  lalpha[, 1] <- lscale + log(foo)
  
  for (i in 2:T) {
    P <- dnorm(y[i], mean = mus[i, ], sd = sigmas[i, ])
    foo <- foo %*% Gamma * P
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
    lalpha[, i] <- log(foo) + lscale
  }
  
  return(lalpha)
}





## Compute the one-step-ahead forecast pseudo-residuals
PseudoResiduals <- function(y, mod, m, bm) {
  params <- extract_parameters(mod$estimate)  # Extract parameters from mod$estimate
  
  ### Intervals to approximate continuous state space
  b <- seq(-bm, bm, length = m + 1) ### Specify boundaries of m intervals
  h <- b[2] - b[1] ### h is the length of each interval
  bstar <- (b[-1] + b[-(m + 1)]) * 0.5 ### Midpoints of the m intervals
  
  ### Set up transition probability matrix
  Gamma <- matrix(0, m, m)
  for (i in 1:m){
    Gamma[i, ] <- h * dnorm(bstar, params$phi_state * bstar[i], params$sigma_state) ### m*m t.p.m. of the approx. HMM
  }
  
  ### Set up stationary distribution
  delta <- h * dnorm(bstar, 0, params$sigma_state / sqrt(1 - params$phi_state^2)) ### stat. initial distribution
  
  mus <- matrix(exp(bstar) * params$b, nrow = length(y), ncol = m, byrow = TRUE)
  sigmas <- matrix(sqrt((exp(bstar) * params$sigma)^2 / (2 * params$a * exp(bstar))), nrow = length(y), ncol = m, byrow = TRUE)
  
  la <- t(lForward(y = y, mod = mod, m = m, bm = bm))  # Compute log-forward probabilities using modified lForward function
  n <- length(y)
  Res <- rep(NA, n)
  pMat <- matrix(NA, nrow = n, ncol = m)
  
  pMat[1, ] <- pnorm(y[1], mean = mus[1, ], sd = (sigmas[1, ]))
  Res[1] <- qnorm(delta %*% pMat[1, ])
  
  for (i in 2:n) {
    pMat[i, ] <- pnorm(y[i], mean = mus[i, ], sd = (sigmas[i, ]))
    c <- max(la[i - 1, ])
    a <- exp(la[i - 1, ] - c) # - c
    weighted_Gamma <- Gamma / sum(a)
    Res[i] <- qnorm(a %*% weighted_Gamma %*% pMat[i, ])
  }
  
  return(list(Res = Res))
}









#----------------------------------------- N-state HMM --------------------------------------------------------

### Modify to take in correct number of parameters
### dependending on which parameters modeled as state dependent
extract_parameters <- function(theta.star, N) {
  a <- exp(theta.star[1]) #:N
  b <- theta.star[(N + 1)] #:(2 * N)
  sigma <- exp(theta.star[(2 * N + 1):(3 * N)]) #:(3 * N)
  Gamma <- diag(N)
  Gamma[!Gamma] <- exp(theta.star[(3 * N + 1):length(theta.star)])
  Gamma <- Gamma / rowSums(Gamma)
  delta <- solve(t(diag(N) - Gamma + 1), rep(1, N))
  
  list(a = a, b = b, sigma = sigma, Gamma = Gamma, delta = delta)
}





## Compute the log-forward probabilities
lForward <- function(y, mod, N) {
  params <- extract_parameters(mod$estimate, N)  # Extract parameters from mod$estimate
  T <- length(y)
  mus <- matrix(params$b, nrow = T, ncol = N, byrow = TRUE)
  sigmas <- matrix(sqrt(params$sigma^2 / (2 * params$a)), nrow = T, ncol = N, byrow = TRUE)
  
  lalpha <- matrix(NA, N, T)
  P <- dnorm(y[1], mean = mus[1, ], sd = sigmas[1, ])
  foo <- params$delta * P
  sumfoo <- sum(foo)
  lscale <- log(sumfoo)
  foo <- foo / sumfoo
  lalpha[, 1] <- lscale + log(foo)
  
  for (i in 2:T) {
    P <- dnorm(y[i], mean = mus[i, ], sd = sigmas[i, ])
    foo <- foo %*% params$Gamma * P
    sumfoo <- sum(foo)
    lscale <- lscale + log(sumfoo)
    foo <- foo / sumfoo
    lalpha[, i] <- log(foo) + lscale
  }
  
  return(lalpha)
}



## Compute the one-step-ahead forecast pseudo-residuals
PseudoResiduals <- function(y, mod, N) {
  params <- extract_parameters(mod$estimate, N)  # Extract parameters from mod$estimate
  mus <- matrix(params$b, nrow = length(y), ncol = N, byrow = TRUE)
  sigmas <- matrix(sqrt(params$sigma^2 / (2 * params$a)), nrow = length(y), ncol = N, byrow = TRUE)
  
  la <- t(lForward(y = y, mod = mod, N = N))  # Compute log-forward probabilities using modified lForward function
  n <- length(y)
  Res <- rep(NA, n)
  pMat <- matrix(NA, nrow = n, ncol = N)
  
  pMat[1, ] <- pnorm(y[1], mean = mus[1, ], sd = (sigmas[1, ]))
  Res[1] <- qnorm(params$delta %*% pMat[1, ])
  
  for (i in 2:n) {
    pMat[i, ] <- pnorm(y[i], mean = mus[i, ], sd = (sigmas[i, ]))
    c <- max(la[i - 1, ])
    a <- exp(la[i - 1, ] - c) # - c
    weighted_Gamma <- params$Gamma / sum(a)
    Res[i] <- qnorm(a %*% weighted_Gamma %*% pMat[i, ])
  }
  
  return(list(Res = Res))
}















#----------------------------------------- PseudoResiduals --------------------------------------------------------



### theta
pseudo_res_2_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2_theta, N = 2)
pseudo_res_3_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3_theta, N = 3)
pseudo_res_4_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4_theta, N = 4)
pseudo_res_5_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5_theta, N = 5)




### kappa
pseudo_res_2_kappa <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2_kappa, N = 2)
pseudo_res_3_kappa <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3_kappa, N = 3)
pseudo_res_4_kappa <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4_kappa, N = 4)
pseudo_res_5_kappa <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5_kappa, N = 5)





### sigma
pseudo_res_2_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2_sigma, N = 2)
pseudo_res_3_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3_sigma, N = 3)
pseudo_res_4_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4_sigma, N = 4)
pseudo_res_5_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5_sigma, N = 5)







### sigma theta
pseudo_res_2_sigma_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2_sigma_theta, N = 2)
pseudo_res_3_sigma_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3_sigma_theta, N = 3)
pseudo_res_4_sigma_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4_sigma_theta, N = 4)
pseudo_res_5_sigma_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5_sigma_theta, N = 5)





### sigma kappa
pseudo_res_2_kappa_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2_kappa_sigma, N = 2)
pseudo_res_3_kappa_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3_kappa_sigma, N = 3)
pseudo_res_4_kappa_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4_kappa_sigma, N = 4)
pseudo_res_5_kappa_sigma <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5_kappa_sigma, N = 5)





### kappa theta
pseudo_res_2_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod2_kappa_theta, N = 2)
pseudo_res_3_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3_kappa_theta, N = 3)
pseudo_res_4_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4_kappa_theta, N = 4)
pseudo_res_5_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5_kappa_theta, N = 5)







### sigma theta kappa
pseudo_res_2_sigma_kappa_theta <- PseudoResiduals(y = as.numeric(yields_df$"3M"), mod = mod2, params = params2, N = 2)
pseudo_res_3_sigma_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod3, N = 3)
pseudo_res_4_sigma_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod4, N = 4)
pseudo_res_5_sigma_kappa_theta <- PseudoResiduals(as.numeric(yields_df$"3M"), mod5, N = 5)






### Cont
theta.star <- ssmmod_sigma_theta_kappa$estimate
pseudo_res_SSM <- PseudoResiduals(as.numeric(yields_df$"3M"), ssmmod_sigma_theta_kappa, m = 200, bm = 3)











#----------------------------------------- SMALL PLOTS --------------------------------------------------------


par(oma=c(3,3,0,1),mar=c(3,3,2,3),mfrow=c(3,4))

### theta
qqnorm(pseudo_res_2_theta$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_theta$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_theta$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_theta$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_theta$Res, col = "steelblue", lwd = 3)



### kappa
qqnorm(pseudo_res_2_kappa$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_kappa$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_kappa$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_kappa$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_kappa$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_kappa$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_kappa$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_kappa$Res, col = "steelblue", lwd = 3)


### sigma
qqnorm(pseudo_res_2_sigma$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_sigma$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_sigma$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_sigma$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_sigma$Res, col = "steelblue", lwd = 3)




mtext(text="Theoretical Quantiles",side=1,line=1,outer=TRUE, cex = 2.7)
mtext(text="Sample Quantiles",side=2,line=0,outer=TRUE, cex = 2.7, adj = 0.5)
mtext(text=expression(theta),side=4,line=-1,outer=TRUE, cex = 2, adj = 0.85)
mtext(text=expression(kappa),side=4,line=-1,outer=TRUE, cex = 2, adj = 0.5)
mtext(text=expression(sigma),side=4,line=-1,outer=TRUE, cex = 2, adj = 0.15)













par(oma=c(3,3,0,1),mar=c(3,3,2,3),mfrow=c(3,4))

### theta sigma
qqnorm(pseudo_res_2_sigma_theta$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_sigma_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_sigma_theta$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_sigma_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_sigma_theta$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_sigma_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_sigma_theta$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_sigma_theta$Res, col = "steelblue", lwd = 3)



### kappa theta
qqnorm(pseudo_res_2_kappa_theta$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_kappa_theta$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_kappa_theta$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_kappa_theta$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_kappa_theta$Res, col = "steelblue", lwd = 3)


### sigma kappa
qqnorm(pseudo_res_2_kappa_sigma$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_kappa_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_kappa_sigma$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_kappa_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_kappa_sigma$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_kappa_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_kappa_sigma$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_kappa_sigma$Res, col = "steelblue", lwd = 3)




mtext(text="Theoretical Quantiles",side=1,line=1,outer=TRUE, cex = 2.7)
mtext(text="Sample Quantiles",side=2,line=0,outer=TRUE, cex = 2.7, adj = 0.5)
mtext(text=expression(theta * ", " * sigma),side=4,line=-1,outer=TRUE, cex = 2, adj = 0.85)
mtext(text=expression(theta * ", " * kappa),side=4,line=-1,outer=TRUE, cex = 2, adj = 0.5)
mtext(text=expression(kappa * ", " * sigma),side=4,line=-1,outer=TRUE, cex = 2, adj = 0.15)










par(oma=c(3,3,0,1),mar=c(3,3,2,3),mfrow=c(3,4))

### theta sigma kappa
qqnorm(pseudo_res_2_sigma_kappa_theta$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_sigma_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_sigma_kappa_theta$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_sigma_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_sigma_kappa_theta$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_sigma_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_sigma_kappa_theta$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_sigma_kappa_theta$Res, col = "steelblue", lwd = 3)



### SSM
qqnorm(pseudo_res_SSM$Res, main = "SSM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_SSM$Res, col = "steelblue", lwd = 3)



mtext(text="Theoretical Quantiles",side=1,line=1,outer=TRUE, cex = 2.7)
mtext(text="Sample Quantiles",side=2,line=0,outer=TRUE, cex = 2.7, adj = 0.5)
mtext(text=expression(sigma * ", " * kappa * ", " * theta),side=4,line=-1,outer=TRUE, cex = 2, adj = 0.85)
mtext(text="Continuous SSM",side=4,line=-1,outer=TRUE, cex = 2, adj = 0.5)





#----------------------------------------- BIG PLOT --------------------------------------------------------



par(mfrow = c(4, 4), mar = c(2, 2, 2, 1), oma = c(3, 3, 0, 2))


### theta
qqnorm(pseudo_res_2_theta$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_theta$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_theta$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_theta$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_theta$Res, col = "steelblue", lwd = 3)



### kappa
qqnorm(pseudo_res_2_kappa$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_kappa$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_kappa$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_kappa$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_kappa$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_kappa$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_kappa$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_kappa$Res, col = "steelblue", lwd = 3)


### sigma
qqnorm(pseudo_res_2_sigma$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_sigma$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_sigma$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_sigma$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_sigma$Res, col = "steelblue", lwd = 3)



### theta sigma
qqnorm(pseudo_res_2_sigma_theta$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_sigma_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_sigma_theta$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_sigma_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_sigma_theta$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_sigma_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_sigma_theta$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_sigma_theta$Res, col = "steelblue", lwd = 3)


mtext(text="",side=1,line=2,outer=TRUE, cex = 2.7)
mtext(text="Sample Quantiles",side=2,line=0.3,outer=TRUE, cex = 2.7, adj = -0.07)
mtext(text=expression(theta),side=4,line=0.3,outer=TRUE, cex = 2, adj = 0.87)
mtext(text=expression(kappa),side=4,line=-0.1,outer=TRUE, cex = 2, adj = 0.63)
mtext(text=expression(sigma),side=4,line=-0.1,outer=TRUE, cex = 2, adj = 0.37)
mtext(text=expression(theta * ", " * sigma),side=4,line=0.5,outer=TRUE, cex = 2, adj = 0.09)




par(mfrow = c(4, 4), mar = c(2, 2, 2, 1), oma = c(3, 3, 0, 2))



### kappa theta
qqnorm(pseudo_res_2_kappa_theta$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_kappa_theta$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_kappa_theta$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_kappa_theta$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_kappa_theta$Res, col = "steelblue", lwd = 3)


### sigma kappa
qqnorm(pseudo_res_2_kappa_sigma$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_kappa_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_kappa_sigma$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_kappa_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_kappa_sigma$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_kappa_sigma$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_kappa_sigma$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_kappa_sigma$Res, col = "steelblue", lwd = 3)





### theta sigma kappa
qqnorm(pseudo_res_2_sigma_kappa_theta$Res, main = "2-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_2_sigma_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_3_sigma_kappa_theta$Res, main = "3-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_3_sigma_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_4_sigma_kappa_theta$Res, main = "4-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_4_sigma_kappa_theta$Res, col = "steelblue", lwd = 3)

qqnorm(pseudo_res_5_sigma_kappa_theta$Res, main = "5-State HMM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_5_sigma_kappa_theta$Res, col = "steelblue", lwd = 3)



### SSM
qqnorm(pseudo_res_SSM$Res, main = "SSM", col = "#901a1E",
       xlab = "", ylab = "", cex.main = 2.7)
qqline(pseudo_res_SSM$Res, col = "steelblue", lwd = 3)


mtext(text="Theoretical Quantiles",side=1,line=2,outer=TRUE, cex = 2.7)
mtext(text="",side=2,line=0.3,outer=TRUE, cex = 2.7, adj = 0.5)
mtext(text=expression(theta * ", " * kappa),side=4,line=0.5,outer=TRUE, cex = 2, adj = 0.87)
mtext(text=expression(kappa * ", " * sigma),side=4,line=0.5,outer=TRUE, cex = 2, adj = 0.63)
mtext(text=expression(sigma * ", " * kappa * ", " * theta),side=4,line=0.5,outer=TRUE, cex = 2, adj = 0.37)
mtext(text=expression(sigma * ", " * kappa * ", " * theta),side=4,line=0.5,outer=TRUE, cex = 2, adj = 0.09)










