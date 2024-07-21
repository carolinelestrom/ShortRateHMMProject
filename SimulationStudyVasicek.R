####################################################################################################################
####################################################################################################################
#---------------------------------------- Vasicek Simulation -------------------------------------------------------
####################################################################################################################
####################################################################################################################


### Set parameters for the Vasicek process
n <- 10000
dt <- 1
kappa_true <- 0.003
theta_true <- 0.04
sigma_true <- 0.0002
r0 <- theta_true

### Function to simulate Vasicek process
simulate_vasicek <- function(n, kappa, theta, sigma, r0) {
  r <- numeric(n)
  r[1] <- r0
  for (i in 2:n) {
    r[i] <- r[i-1] + kappa * (theta - r[i-1]) * dt + sigma * sqrt(dt) * rnorm(1)
  }
  return(r)
}

### Function to calculate the log-likelihood of the Vasicek process
vasicek_loglik <- function(params, r) {
  a <- exp(params[1])
  b <- params[2]
  sigma <- exp(params[3])
  
  n <- length(r)
  dt <- 1  ### Assuming unit time steps
  
  ll <- 0
  for (i in 2:n) {
    expected_r <- r[i-1] * exp(-a * dt) + b * (1 - exp(-a * dt))
    var_r <- sigma^2 * (1 - exp(-2 * a * dt)) / (2 * a)
    ll <- ll + dnorm(r[i], mean = expected_r, sd = sqrt(var_r), log = TRUE)
  }
  
  return(-ll)  ### Return negative log-likelihood for minimization
}

### Function to fit the Vasicek model
fit_vasicek <- function(r) {
  init_params <- c(log(kappa_true), theta_true, log(sigma_true))
  result <- nlm(f = vasicek_loglik, p = init_params, r = r, print.level = 0)
  return(c(exp(result$estimate[1]), result$estimate[2], exp(result$estimate[3])))
}

### Number of simulations
num_simulations <- 100

### List to store the results
simulated_series_Vasicek <- vector("list", num_simulations)
fitted_params_Vasicek <- matrix(NA, nrow = num_simulations, ncol = 3)

### Simulate and fit the Vasicek model


for (i in 1:num_simulations) {
  set.seed(i)  ### For reproducibility
  
  ### Simulate the Vasicek process
  simulated_series_Vasicek[[i]] <- simulate_vasicek(n, kappa_true, theta_true, sigma_true, r0)
  
  ### Fit the Vasicek model to the simulated series
  fitted_params_Vasicek[i, ] <- fit_vasicek(simulated_series_Vasicek[[i]])
}


VasicekFit <- numeric(3)
for (i in 1:3){
  VasicekFit[i] <- mean(fitted_params_Vasicek[,i])
}






####################################################################################################################
#------------------------------------------ Load RData -------------------------------------------------------------
####################################################################################################################




setwd("~/Documents/KU/Bielefeld/RCode")
load("fitted_params_Vasicek.RData")









fitted_params_Vasicek

colMeans(fitted_params_Vasicek)
colSds(fitted_params_Vasicek)


kappa_true <- 0.02
theta_true <- 0.06
sigma_true <- 0.003


### Plot the fitted parameters
SimVasicek <- as.data.frame(fitted_params_Vasicek)
colnames(SimVasicek) <- c("kappa", "theta", "sigma")
SimVasicek$Index <- c(1:num_simulations)






param_names <- c("kappa", "theta", "sigma")
SimVasicek_param <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(SimVasicek$kappa,
            SimVasicek$theta,
            SimVasicek$sigma)
)
true_vals <- data.frame(
  Group = c("kappa", "theta", "sigma"),
  TrueValue = c(kappa_true, theta_true, sigma_true)
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(kappa), expression(sigma), expression(theta))

ggplot(SimVasicek_param, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.5) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.5, show.legend = FALSE) +
  scale_color_manual(values=c("#404080", "#FFBCD9", "#87A96B")) +
  geom_hline(yintercept = kappa_true, col = "#404080", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = sigma_true, col = "#FFBCD9", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = theta_true, col = "#87A96B", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_vals, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = expression("Estimated " * kappa * ", " * theta * ", " * sigma),
       x = "",
       y = "") +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = latex_labels)








#------------------------------------------ kappa -------------------------------------------------------------


param_names <- c("kappa")
SimVasicek_param <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(SimVasicek$kappa)
)
true_vals <- data.frame(
  Group = c("kappa"),
  TrueValue = c(kappa_true)
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(kappa))


ggplot(SimVasicek_param, aes(x = Parameter, y = Value)) +
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
  scale_y_continuous(limits = c(0.015,0.025)) + 
  theme(plot.margin=unit(c(.01,.3,-1.3,-.7),"cm")) ### c(top,right,bottom,left)






#------------------------------------------ theta -------------------------------------------------------------


param_names <- c("theta")
SimVasicek_param <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(SimVasicek$theta)
)
true_vals <- data.frame(
  Group = c("theta"),
  TrueValue = c(theta_true)
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(theta))


ggplot(SimVasicek_param, aes(x = Parameter, y = Value)) +
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
  scale_y_continuous(limits = c(0.055, 0.065)) + 
  theme(plot.margin=unit(c(.01,.3,-1.3,-.7),"cm")) ### c(top,right,bottom,left)





#------------------------------------------ sigma -------------------------------------------------------------


param_names <- c("sigma")
SimVasicek_param <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(SimVasicek$sigma)
)
true_vals <- data.frame(
  Group = c("sigma"),
  TrueValue = c(sigma_true)
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(sigma))


ggplot(SimVasicek_param, aes(x = Parameter, y = Value)) +
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
  scale_y_continuous(limits = c(0.002, 0.004)) + 
  theme(plot.margin=unit(c(.01,.3,-1.3,-.7),"cm")) ### c(top,right,bottom,left)






