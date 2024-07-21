####################################################################################################################
####################################################################################################################
#---------------------------------------- Predicting Vasicek -------------------------------------------------------
####################################################################################################################
####################################################################################################################

PredDat <- as.numeric(yields_df$"3M")[1:9535] ### NERVOUS
PredDat <- as.numeric(yields_df$"3M")[1:8760] ###  MODERATE
PredDat <- as.numeric(yields_df$"3M")[1:9200] ###  CALM
mod1prednervous <- nlm(f = vasicek_loglik, p = init_params, r = PredDat, print.level = 2)
mod1predmoderate <- nlm(f = vasicek_loglik, p = init_params, r = PredDat, print.level = 2)
mod1predcalm <- nlm(f = vasicek_loglik, p = init_params, r = PredDat, print.level = 2)



### Prediction of the short rate
n <- 21*3 + 1
dt <- 1
kappa <- exp(mod1prednervous$estimate[1])
theta <- mod1prednervous$estimate[2]
sigma <- exp(mod1prednervous$estimate[3])
r0 <- PredDat[length(PredDat):length(PredDat)]
rVasiSim <- matrix(0, nrow = 1003, ncol = n)
rVasiSim[,1] <- r0

for (j in 1:1000){
  for (i in 2:n){
    rVasiSim[j,i] <- rVasiSim[j,i-1] + kappa * (theta - rVasiSim[j,i-1]) * dt + sigma * sqrt(dt) * rnorm(1)
  }
}

for(i in 1:n){
  rVasiSim[1001,i] <- mean(rVasiSim[,i])
  rVasiSim[1002,i] <- quantile(rVasiSim[,i], 0.025)
  rVasiSim[1003,i] <- quantile(rVasiSim[,i], 0.975)
}




par(mfrow = c(1,1), mar = c(2.1, 3.1, 4.1, 1.1))
### Plot simulated short rate
plot(rVasiSim[1,], type = "l", xlab = "", ylab = "", 
     main = "", ylim = c(-0.003,0.005),
     col = "#00000025")
for (j in 2:1000){
  points(rVasiSim[j,], type = "l", col = "#00000025")
}
points(rVasiSim[1001,], type = "b", col = "#800000", lwd = 3)
points(rVasiSim[1002,], type = "b", col = "#CD5C5C", lwd = 3)
points(rVasiSim[1003,], type = "b", col = "#CD5C5C", lwd = 3)
points(as.numeric(yields_df$"3M")[length(PredDat):(length(PredDat) + 21*3)], type = "b", col = "#1034A6", lwd = 3)
mtext(text="Time",side=1,line=2.7,outer=FALSE, cex = 2)
mtext(text="Short Rate",side=2,line=2.3,outer=FALSE, cex = 2, adj = 0.5)
mtext(text="Prediction Vasicek Model - 2021",side=1,line=-39,outer=FALSE, cex = 3, adj = 0)




MSE_Vasicek <- rep(0, n)
for (i in 2:n){
  MSE_Vasicek[i] <- 1/999 * sum(as.numeric(yields_df$"3M")[length(PredDat) + (i - 1)] - rVasiSim[,i])^2
}


round(MSE_Vasicek, 7)
mean(MSE_Vasicek)
