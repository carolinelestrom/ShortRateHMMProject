####################################################################################################################
####################################################################################################################
#------------------------------------------- State Prediction ------------------------------------------------------
####################################################################################################################
####################################################################################################################




predict <- function(theta.star, x, m, bm, steps = 1) {
  ### S_t model formulation
  phi_theta <- plogis(theta.star[1]) ### Mean reversion parameter of theta_t
  sigma_theta <- exp(theta.star[2]) ### Vol/noise parameter of theta_t
  
  ### r_t model formulation
  kappa <- exp(theta.star[3]) ### r_t \sim N(theta, sigma^2/(2*kappa))
  beta_theta <- theta.star[4] ### r_t \sim N(theta, sigma^2/(2*kappa))
  sigma <- exp(theta.star[5]) ### r_t \sim N(theta, sigma^2/(2*kappa))
  
  ### Intervals to approximate continuous state space
  b <- seq(-bm, bm, length = m + 1) ### Specify boundaries of m intervals
  h <- b[2] - b[1] ### h is the length of each interval
  bstar <- (b[-1] + b[-(m + 1)]) * 0.5 ### Midpoints of the m intervals
  
  ### Set up transition probability matrix
  Gamma <- matrix(0, m, m)
  for (i in 1:m){
    Gamma[i, ] <- h * dnorm(bstar, phi_theta * bstar[i], sigma_theta) ### m*m t.p.m. of the approx. HMM
  }
  
  ### Set up stationary distribution
  delta <- h * dnorm(bstar, 0, sigma_theta / sqrt(1 - phi_theta^2)) ### stat. initial distribution
  
  ### Forward algorithm
  foo <- delta * dnorm(x[1], exp(bstar) * beta_theta, sqrt(exp(bstar) * sigma)^2/(2* exp(bstar) * kappa)) ### r_t \sim N(theta exp(theta_t), sigma^2/(2*kappa))
  l <- log(sum(foo))
  phi <- foo / sum(foo)
  for (t in 2:length(x)){
    foo <- phi %*% Gamma * dnorm(x[t], exp(bstar) * beta_theta, sqrt(exp(bstar) * sigma)^2/(2* exp(bstar) *kappa))
    l <- l + log(sum(foo))
    phi <- foo / sum(foo)
  }
  prediction <- b[which.max(as.vector(phi %*% matrix.power(Gamma, steps)))]
  
  return(prediction)
}

predict(theta.star = theta.star, x = PredDat, m = 200, bm = 3, steps = 300)


####################################################################################################################
####################################################################################################################
#------------------------------------------ Predicting SSM ---------------------------------------------------------
####################################################################################################################
####################################################################################################################

PredDat <- as.numeric(yields_df$"3M")[1:9535] ### NERVOUS
PredDat <- as.numeric(yields_df$"3M")[1:8760] ###  MODERATE
PredDat <- as.numeric(yields_df$"3M")[1:9200] ###  CALM



### Model from simulated data
theta.star <- c(qlogis(0.999), log(0.2), log(0.001), 0.05, log(0.0001))
theta.star <- c(ssmmod_sigma_theta_kappa$estimate[1],
                ssmmod_sigma_theta_kappa$estimate[2],
                ssmmod_sigma_theta_kappa$estimate[3],
                ssmmod_sigma_theta_kappa$estimate[4],
                ssmmod_sigma_theta_kappa$estimate[5])
ssmmod_sigma_theta_kappa_pred2022 <- nlm(mllk, theta.star, x = PredDat,
                                         m = 200, bm = 3, print.level = 2)
ssmmod_sigma_theta_kappa_predmoderate <- nlm(mllk, theta.star, x = PredDat,
                                             m = 200, bm = 3, print.level = 2)
ssmmod_sigma_theta_kappa_predcalm <- nlm(mllk, theta.star, x = PredDat,
                                         m = 200, bm = 3, print.level = 2)






### Prediction of states
theta.star <- c(ssmmod_sigma_theta_kappa_predcalm$estimate[1],
                ssmmod_sigma_theta_kappa_predcalm$estimate[2],
                ssmmod_sigma_theta_kappa_predcalm$estimate[3],
                ssmmod_sigma_theta_kappa_predcalm$estimate[4],
                ssmmod_sigma_theta_kappa_predcalm$estimate[5])




future <- numeric(21*3)
for (i in 1:(21*3)){
  future[i] <- predict(theta.star = theta.star, x = PredDat, m = 200, bm = 3, steps = i)
}



### Prediction of the short rate
n <- 21*3 + 1
dt <- 1
kappa <- exp(ssmmod_sigma_theta_kappa_predcalm$estimate[3])
theta <- ssmmod_sigma_theta_kappa_predcalm$estimate[4]
sigma <- exp(ssmmod_sigma_theta_kappa_predcalm$estimate[5])
r0 <- PredDat[length(PredDat):length(PredDat)]
rSSMSim <- matrix(0, nrow = 1003, ncol = n)
rSSMSim[,1] <- r0

for (j in 1:1000){
  for (i in 2:n){
    rSSMSim[j,i] <- rSSMSim[j,i-1] + kappa * exp(future[i-1]) * (theta * exp(future[i-1]) - rSSMSim[j,i-1]) * dt + sigma * exp(future[i-1]) * sqrt(dt) * rnorm(1)
  }
}


for(i in 1:n){
  rSSMSim[1001,i] <- mean(rSSMSim[,i])
  rSSMSim[1002,i] <- quantile(rSSMSim[,i], 0.025)
  rSSMSim[1003,i] <- quantile(rSSMSim[,i], 0.975)
}



### Plot simulated short rate
plot(rSSMSim[1,], type = "l", xlab = "", ylab = "", 
     main = "", ylim = c(-0.003,0.005),
     col = "#00000025")
for (j in 2:1000){
  points(rSSMSim[j,], type = "l", col = "#00000025")
}
points(rSSMSim[1001,], type = "b", col = "#800000", lwd = 3)
points(rSSMSim[1002,], type = "b", col = "#CD5C5C", lwd = 3)
points(rSSMSim[1003,], type = "b", col = "#CD5C5C", lwd = 3)
points(as.numeric(yields_df$"3M")[length(PredDat):(length(PredDat) + 21*3)], type = "b", col = "#1034A6", lwd = 3)
mtext(text="Time",side=1,line=2.7,outer=FALSE, cex = 2)
mtext(text="Short Rate",side=2,line=2.3,outer=FALSE, cex = 2, adj = 0.5)
mtext(text="Prediction Continuous SSM - 2021",side=1,line=-39,outer=FALSE, cex = 3, adj = 0)







MSE_ContState <- rep(0, n)
for (i in 2:n){
  MSE_ContState[i] <- 1/999 * sum(as.numeric(yields_df$"3M")[length(PredDat) + (i - 1)] - r[,i])^2
}

round(MSE_ContState, 7)
mean(MSE_ContState)
