####################################################################################################################
####################################################################################################################
#---------------------------------------- Simulating Vasicek -------------------------------------------------------
####################################################################################################################
####################################################################################################################



### Prediction of the short rate
n <- 10000
T <- 10
dt <- T/n
kappa <- 0.6 
theta <- 0.03 
sigma <- 0.002 
r0 <- 0.027
VasiSim <- numeric(n)
VasiSim[1] <- r0

alpha <- 0.95
z <- qnorm((1 + alpha) / 2)
tau <- numeric(n)
for (i in 1:n){
  tau[i] <- i * T/n
}
CI_upper <- numeric(n)
CI_lower <- numeric(n)
CI_upper[1] <- (r0 * exp(-kappa * tau[1]) + theta * (1 - exp(-kappa * tau[1]))) + z * sqrt(sigma^2 * (1 - exp(-2 * kappa * tau[1])) / (2 * kappa))
CI_lower[1] <- (r0 * exp(-kappa * tau[1]) + theta * (1 - exp(-kappa * tau[1]))) - z * sqrt(sigma^2 * (1 - exp(-2 * kappa * tau[1])) / (2 * kappa))



for (i in 2:n){
  VasiSim[i] <- VasiSim[i-1] + kappa * (theta - VasiSim[i-1]) * dt + sigma * sqrt(dt) * rnorm(1)
  expected_r <- VasiSim[i-1] * exp(-kappa * dt) + theta * (1 - exp(-kappa * dt))
  var_r <- sigma^2 * (1 - exp(-2 * kappa * dt)) / (2 * kappa)
  CI_upper[i] <- (r0 * exp(-kappa * tau[i]) + theta * (1 - exp(-kappa * tau[i]))) + z * sqrt(sigma^2 * (1 - exp(-2 * kappa * tau[i])) / (2 * kappa))
  CI_lower[i] <- (r0 * exp(-kappa * tau[i]) + theta * (1 - exp(-kappa * tau[i]))) - z * sqrt(sigma^2 * (1 - exp(-2 * kappa * tau[i])) / (2 * kappa))
}

CI_upper_stat <- theta + z * sqrt(sigma^2 / (2 * kappa))
CI_lower_stat <- theta - z * sqrt(sigma^2 / (2 * kappa))



### Plot simulated short rate


data <- data.frame(
  Time = tau,
  SR = VasiSim,
  Lower = CI_lower,
  Upper = CI_upper,
  Lower_stat = CI_lower_stat,
  Upper_stat = CI_upper_stat
)

ggplot(data, aes(x = Time)) +
  geom_line(aes(y = SR), color = "#87A96B", size = 1) +
  geom_ribbon(aes(ymin = Lower_stat, ymax = Upper_stat), alpha = 0.2, fill = "#96C8A2", col = "#96C8A2") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.2, fill = "#801818", col = "#801818", linetype = "longdash") +
  theme_bw() +
  labs(title = "Simulated Short Rate with 95% Confidence Intervals",
       x = "Time in Years",
       y = "Short Rate") +
  scale_y_continuous(limits = c(0.023,0.039)) +
  theme(plot.title = element_text(size=37, hjust=0)) +
  theme(axis.title = element_text(size=27)) +
  theme(axis.text.x = element_text(size=23),
        axis.text.y = element_text(size=23)) +
  annotate(geom="text", x=8, y=0.039, label=paste("CI from cond. dist."),
           color="#801818", size = 17) +
  annotate(geom="text", x=7.85, y=0.037, label=paste("CI from stat. dist."),
           color="#96C8A2", size = 17)


