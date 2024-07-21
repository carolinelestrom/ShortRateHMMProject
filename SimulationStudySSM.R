####################################################################################################################
####################################################################################################################
#--------------------------------- Continuous State Space Simulation -----------------------------------------------
####################################################################################################################
####################################################################################################################



### Parameters for the simulation
n <- 10000
dt <- 1
kappa_true <- 0.003
theta_true <- 0.04
sigma_true <- 0.0002
phi_s_true <- 0.999  ### AR(1) parameter
sigma_s_true <- 0.002  ### Standard deviation/Volatility of the white noise
r0 <- theta_true

### Function to simulate the regime-switching Vasicek process
simulate_SSM <- function(n, kappa_true, theta_true, sigma_true, phi_s_true, sigma_s_true, r0) {
  r <- numeric(n)
  S <- numeric(n)
  S[1] <- 0.01  ### Initialize state sequence
  r[1] <- r0
  
  for (i in 2:n) {
    S[i] <- phi_s_true * S[i - 1] + sigma_s_true * rnorm(1, mean = 0, sd = 1)  ### Simulate state sequence
    r[i] <- r[i - 1] + kappa_true * exp(S[i]) * (theta_true * exp(S[i]) - r[i - 1]) * dt + sigma_true * exp(S[i]) * sqrt(dt) * rnorm(1)
  }
  
  return(list(r = r, S = S))
}

### Function to calculate the log-likelihood of the Vasicek process
mllk <- function(theta.star, x, m, bm){
  ### S_t model formulation
  phi_state <- plogis(theta.star[1]) ### Mean reversion parameter of state
  sigma_state <- exp(theta.star[2]) ### Vol/noise parameter of state
  
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

### Initial parameters for fitting
theta.star <- c(qlogis(phi_s_true), log(sigma_s_true), log(kappa_true), theta_true, log(sigma_true))

### Number of simulations
num_simulations <- 100

### Matrix to store the results
fitted_params_SSM <- matrix(NA, nrow = num_simulations, ncol = length(theta.star))


### Loop over simulations
for (i in 1:num_simulations) {
  set.seed(i)  ### For reproducibility
  
  ### Simulate the Vasicek process
  simulation <- simulate_SSM(n, kappa_true, theta_true, sigma_true, phi_s_true, sigma_s_true, r0)
  r <- simulation$r
  
  ### Fit the Vasicek model to the simulated series
  result <- nlm(mllk, theta.star, x = r, m = 200, bm = 3, print.level = 0)
  
  ### Store the fitted parameters
  fitted_params_SSM[i, ] <- result$estimate
}



### Loop over simulations
for (i in 35:num_simulations) {
  #set.seed(23+i)  ### For reproducibility
  
  ### Simulate the Vasicek process
  simulation <- simulate_SSM(n, kappa_true, theta_true, sigma_true, phi_s_true, sigma_s_true, r0)
  r <- simulation$r
  
  ### Fit the Vasicek model to the simulated series
  result <- nlm(mllk, theta.star, x = r, m = 200, bm = 3, print.level = 0)
  
  ### Store the fitted parameters
  fitted_params_SSM[i, ] <- result$estimate
}






### Convert the fitted parameters to their original scales
fitted_params_original_SSM <- cbind(
  plogis(fitted_params_SSM[, 1]), 
  exp(fitted_params_SSM[, 2]), 
  exp(fitted_params_SSM[, 3]), 
  fitted_params_SSM[, 4], 
  exp(fitted_params_SSM[, 5])
)

### Print the fitted parameters
print(fitted_params_original_SSM)

SSMFit <- numeric(5)
for (i in 1:5){
  SSMFit[i] <- mean(fitted_params_original_SSM[,i])
}


### Plot the fitted parameters
SimSSM <- as.data.frame(fitted_params_original_SSM)
colnames(SimSSM) <- c("phi_state", "sigma_state", "kappa", "theta", "sigma")
SimSSM$Index <- c(1:num_simulations)



SimSSM %>%
  ggplot(aes(x = Index, y = kappa)) +
  geom_point(col = "#87A96B", size = 3) +
  theme_bw() +
  xlab("Number of Fitted Models") +
  ylab("Kappa") +
  ggtitle("Simulated kappa models") +
  theme(plot.title = element_text(size=17, hjust=0)) +
  theme(axis.title = element_text(size=13)) +
  theme(axis.text.x = element_text(size=9, angle = 0, vjust = 0.7),
        axis.text.y = element_text(size=9)) +
  scale_color_manual(values=c("#CC5500", "#404080", "#87A96B", "#FFBCD9", "#996666")) + 
  geom_hline(yintercept = kappa_true[1], col = "#CB4154", linewidth = 1.3, linetype = "longdash") +
  annotate(geom="text", x=90.3, y=0.0033, label=paste("True Kappa ---"),
           color="#CB4154", size = 7) +
  annotate(geom="text", x=93, y=0.0032, label=paste("Estimated Kappa ---"),
           color="#87A96B", size = 7) +
  scale_y_continuous(breaks=c(0.001, 0.002, 0.003), limits = c(0.001, 0.0035))










####################################################################################################################
#------------------------------------------ Load RData -------------------------------------------------------------
####################################################################################################################




setwd("~/Documents/KU/Bielefeld/RCode")
load("fitted_params_SSM.RData")







kappa_true <- 0.02
theta_true <- 0.06
sigma_true <- 0.003
phi_s_true <- 0.999
sigma_s_true <- 0.009





### Convert the fitted parameters to their original scales
fitted_params_original_SSM <- cbind(
  plogis(fitted_params_SSM[, 1]), 
  exp(fitted_params_SSM[, 2]), 
  exp(fitted_params_SSM[, 3]), 
  fitted_params_SSM[, 4], 
  exp(fitted_params_SSM[, 5])
)



fitted_params_original_SSM



colMeans(fitted_params_original_SSM)
colSds(fitted_params_original_SSM)




### Plot the fitted parameters
SimSSM <- as.data.frame(fitted_params_original_SSM)
colnames(SimSSM) <- c("phi_state", "sigma_state", "kappa", "theta", "sigma")
SimSSM$Index <- c(1:num_simulations)






param_names <- c("phi_state", "sigma_state", "kappa", "theta", "sigma")
SimSSM_param <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(SimSSM$phi_state,
            SimSSM$sigma_state,
            SimSSM$kappa,
            SimSSM$theta,
            SimSSM$sigma)
)
true_vals <- data.frame(
  Group = c("phi_state", "sigma_state", "kappa", "theta", "sigma"),
  TrueValue = c(phi_s_true, sigma_s_true, kappa_true, theta_true, sigma_true)
)





#------------------------------------------ kappa -------------------------------------------------------------


param_names <- c("kappa")
SimSSM_param <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(SimSSM$kappa)
)
true_vals <- data.frame(
  Group = c("kappa"),
  TrueValue = c(kappa_true)
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(kappa))


ggplot(SimSSM_param, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 2) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, show.legend = FALSE, size = 7) +
  scale_color_manual(values=c("#404080", "#FFBCD9", "#87A96B")) +
  geom_hline(yintercept = kappa_true, col = "#404080", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_vals, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = expression("Estimated speed of mean reversion, " * kappa),
       x = "",
       y = "") +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = latex_labels) + 
  scale_y_continuous(limits = c(0.005,0.035)) + 
  theme(plot.margin=unit(c(.01,.3,-1.3,-.7),"cm")) ### c(top,right,bottom,left)






#------------------------------------------ theta -------------------------------------------------------------


param_names <- c("theta")
SimSSM_param <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(SimSSM$theta)
)
true_vals <- data.frame(
  Group = c("theta"),
  TrueValue = c(theta_true)
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(theta))


ggplot(SimSSM_param, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 2) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, show.legend = FALSE, size = 7) +
  scale_color_manual(values=c("#FFBCD9", "#FFBCD9", "#87A96B")) +
  geom_hline(yintercept = theta_true, col = "#FFBCD9", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_vals, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = expression("Estimated long-run mean, " * theta),
       x = "",
       y = "") +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = latex_labels) +
  scale_y_continuous(limits = c(0.03, 0.1)) + 
  theme(plot.margin=unit(c(.01,.3,-1.3,-.7),"cm")) ### c(top,right,bottom,left)





#------------------------------------------ sigma -------------------------------------------------------------


param_names <- c("sigma")
SimSSM_param <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(SimSSM$sigma)
)
true_vals <- data.frame(
  Group = c("sigma"),
  TrueValue = c(sigma_true)
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(sigma))


ggplot(SimSSM_param, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 2) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 7, show.legend = FALSE) +
  scale_color_manual(values=c("#87A96B", "#FFBCD9", "#87A96B")) +
  geom_hline(yintercept = sigma_true, col = "#87A96B", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_vals, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = expression("Estimated volatility, " * sigma),
       x = "",
       y = "") +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = latex_labels) +
  scale_y_continuous(limits = c(0.001, 0.01)) + 
  theme(plot.margin=unit(c(.01,.3,-1.3,-.7),"cm")) ### c(top,right,bottom,left)








#------------------------------------------ phi_state -------------------------------------------------------------


param_names <- c("phi_state")
SimSSM_param <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(SimSSM$phi_state)
)
true_vals <- data.frame(
  Group = c("phi_state"),
  TrueValue = c(phi_s_true)
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(phi["S"]))


ggplot(SimSSM_param, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 2) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 7, show.legend = FALSE) +
  scale_color_manual(values=c("#996666", "#FFBCD9", "#87A96B")) +
  geom_hline(yintercept = phi_s_true, col = "#996666", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_vals, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = expression("Estimated mean-reversion of state process, " * phi["S"]),
       x = "",
       y = "") +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = latex_labels) +
  scale_y_continuous(limits = c(0.99, 1.01)) + 
  theme(plot.margin=unit(c(.01,.3,-1.3,-.7),"cm")) ### c(top,right,bottom,left)











#------------------------------------------ sigma_state -------------------------------------------------------------


param_names <- c("sigma_state")
SimSSM_param <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(SimSSM$sigma_state)
)
true_vals <- data.frame(
  Group = c("sigma_state"),
  TrueValue = c(sigma_s_true)
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(sigma["S"]))


ggplot(SimSSM_param, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 2) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 7, show.legend = FALSE) +
  scale_color_manual(values=c("#CC5500", "#FFBCD9", "#87A96B")) +
  geom_hline(yintercept = sigma_s_true, col = "#CC5500", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_vals, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = expression("Estimated volatility of state process, " * sigma["S"]),
       x = "",
       y = "") +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = latex_labels) +
  scale_y_continuous(limits = c(0.00001, 0.0091)) + 
  theme(plot.margin=unit(c(.01,.3,-1.3,-.7),"cm")) ### c(top,right,bottom,left)








