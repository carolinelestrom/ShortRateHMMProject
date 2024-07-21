####################################################################################################################
####################################################################################################################
#------------------------------------------- State Prediction ------------------------------------------------------
####################################################################################################################
####################################################################################################################


predict <- function(theta.star, x, N, steps = 1){
  a <- exp(theta.star[1:N])
  b <- theta.star[(N + 1):(2 * N)]
  sigma <- exp(theta.star[(2 * N + 1):(3 * N)])
  
  ### Initialize Gamma matrix
  Gamma <- diag(N)
  off_diag_elements <- exp(theta.star[(3 * N + 1):(length(theta.star))])
  Gamma[!Gamma] <- off_diag_elements
  Gamma <- Gamma / rowSums(Gamma)
  
  ### Compute the stationary distribution delta
  eig <- eigen(t(Gamma))
  delta <- Re(eig$vectors[,1])
  delta <- delta / sum(delta)
  #delta <- rep(1 / N, times = N)
  
  ### Initialize allprobs matrix
  allprobs <- matrix(1, length(x), N)
  ind <- which(!is.na(x))
  
  for (j in 1:N){
    allprobs[ind, j] <- dnorm(x[ind], mean = b[j], sd = sqrt(sigma[j]^2 / (2 * a[j])))
  }
  
  foo <- delta %*% diag(allprobs[1, ])
  l <- log(sum(foo))
  phi <- foo / sum(foo)
  
  for (t in 2:length(x)){
    if (t > nrow(allprobs)) {
      stop(paste("Subscript out of bounds: t =", t, "exceeds nrow(allprobs) =", nrow(allprobs)))
    }
    foo <- phi %*% Gamma %*% diag(allprobs[t, ])
    l <- l + log(sum(foo))
    phi <- foo / sum(foo)
  }
  prediction <- which.max(as.vector(phi %*% matrix.power(Gamma, steps)))
  return(prediction)
}




####################################################################################################################
####################################################################################################################
#------------------------------------------ Predicting HMM ---------------------------------------------------------
####################################################################################################################
####################################################################################################################

PredDat <- as.numeric(yields_df$"3M")[1:9535] ### NERVOUS
PredDat <- as.numeric(yields_df$"3M")[1:8760] ###  MODERATE
PredDat <- as.numeric(yields_df$"3M")[1:9200] ###  CALM


### Model from simulated data
N <- 5
theta.star <- c(mod5$estimate[1:N],
                mod5$estimate[(N + 1):(2 * N)],
                mod5$estimate[(2 * N + 1):(3 * N)],
                mod5$estimate[(3 * N + 1):(length(theta.star))])
theta.star <- mod5$estimate


mod5pred2022 <- nlm(mllk, theta.star, x = PredDat, N = N, print.level = 2, iterlim = 10000)
mod5predmoderate <- nlm(mllk, theta.star, x = PredDat, N = N, print.level = 2, iterlim = 10000)
mod5predcalm <- nlm(mllk, theta.star, x = PredDat, N = N, print.level = 2, iterlim = 10000)


Gamma5pred_hat <- diag(N)
Gamma5pred_hat[!Gamma5pred_hat] <- exp(mod5predcalm$estimate[(3 * N + 1):(length(theta.star))])
Gamma5pred_hat <- Gamma5pred_hat / rowSums(Gamma5pred_hat)

viterbi5 <- viterbi(x = PredDat, 
                    a = exp(mod5pred$estimate[1:N]),
                    b = mod5pred$estimate[(N + 1):(2 * N)],
                    sigma = exp(mod5pred$estimate[(2 * N + 1):(3 * N)]),
                    Gamma = Gamma5pred_hat,
                    delta = solve(t(diag(N) - Gamma5pred_hat + 1), rep(1, N)),
                    N = N)

plot(viterbi5)
table(viterbi5)



### Prediction of states
theta.star <- c(mod5predcalm$estimate[1:N],
                mod5predcalm$estimate[(N + 1):(2 * N)],
                mod5predcalm$estimate[(2 * N + 1):(3 * N)],
                mod5predcalm$estimate[(3 * N + 1):(length(theta.star))])



future <- numeric(21*3)
for (i in 1:(21*3)){
  future[i] <- predict(theta.star = theta.star, x = PredDat, steps = i, N = N)
}



### Prediction of the short rate
n <- 21*3 + 1
dt <- 1
N <- 5
kappa <- exp(mod5predcalm$estimate[1:N])
theta <- mod5predcalm$estimate[(N + 1):(2 * N)]
sigma <- exp(mod5predcalm$estimate[(2 * N + 1):(3 * N)])
Gamma <- Gamma5pred_hat
delta <- solve(t(diag(N) - Gamma5pred_hat + 1), rep(1, N)) ### Stationary distribution
r0 <- PredDat[length(PredDat):length(PredDat)]
r5Sim <- matrix(0, nrow = 1003, ncol = n)
r5Sim[,1] <- r0

for (j in 1:1000){
  for (i in 2:n){
    r5Sim[j,i] <- r5Sim[j,i-1] + kappa[future[i-1]] * (theta[future[i-1]] - r5Sim[j,i-1]) * dt + sigma[future[i-1]] * sqrt(dt) * rnorm(1)
  }
}

for(i in 1:n){
  r5Sim[1001,i] <- mean(r5Sim[,i])
  r5Sim[1002,i] <- quantile(r5Sim[,i], 0.025)
  r5Sim[1003,i] <- quantile(r5Sim[,i], 0.975)
}



### Plot simulated short rate
plot(r5Sim[1,], type = "l", xlab = "", ylab = "", 
     main = "", ylim = c(0.002,0.01),
     col = "#00000025")
for (j in 2:1000){
  points(r5Sim[j,], type = "l", col = "#00000025")
}
points(r5Sim[1001,], type = "b", col = "#800000", lwd = 3)
points(r5Sim[1002,], type = "b", col = "#CD5C5C", lwd = 3)
points(r5Sim[1003,], type = "b", col = "#CD5C5C", lwd = 3)
points(as.numeric(yields_df$"3M")[length(PredDat):(length(PredDat) + 21*3)], type = "b", col = "#1034A6", lwd = 3)
mtext(text="Time",side=1,line=2.7,outer=FALSE, cex = 2)
mtext(text="Short Rate",side=2,line=2.3,outer=FALSE, cex = 2, adj = 0.5)
mtext(text="Prediction w/ 5-state HMM - 2022",side=1,line=-39,outer=FALSE, cex = 3, adj = 0)






MSE_5State <- rep(0, n)
for (i in 2:n){
  MSE_5State[i] <- 1/999 * sum(as.numeric(yields_df$"3M")[length(PredDat) + (i - 1)] - r5Sim[,i])^2
}



round(MSE_5State, 7)
mean(MSE_5State)


