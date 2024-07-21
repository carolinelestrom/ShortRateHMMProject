####################################################################################################################
####################################################################################################################
#---------------------------------------- N-state Simulation -------------------------------------------------------
####################################################################################################################
####################################################################################################################



### Parameters
n <- 10000
dt <- 1
N <- 5
kappa_true <- c(0.002, 0.0025, 0.004, 0.0075, 0.005)
theta_true <- c(0.02, 0.03, 0.04, 0.035, 0.05)
sigma_true <- c(0.0002, 0.00017, 0.00018, 0.000099, 0.00035)
Gamma_true <- diag(N)
Gamma_true[!Gamma_true] <- seq(0.0001, 0.0005, length = N * (N - 1))
Gamma_true <- Gamma_true / rowSums(Gamma_true)
delta_true <- solve(t(diag(N) - Gamma_true + 1), rep(1, N))
r0 <- theta_true[1]

### Function to simulate the regime-switching Vasicek process
simulate_HMM <- function(n, N, kappa_true, theta_true, sigma_true, Gamma_true, delta_true, r0) {
  r <- numeric(n)
  S <- numeric(n)
  S[1] <- sample(x = 1:N, size = 1, prob = delta_true)
  r[1] <- r0
  
  for (i in 2:n) {
    S[i] <- sample(x = 1:N, size = 1, prob = Gamma_true[S[i - 1], ])
    r[i] <- r[i - 1] + kappa_true[S[i]] * (theta_true[S[i]] - r[i - 1]) * dt + sigma_true[S[i]] * sqrt(dt) * rnorm(1)
  }
  
  return(list(r = r, S = S))
}

### Function to calculate the log-likelihood of the Vasicek process
mllk <- function(theta.star, x, N){
  a <- exp(theta.star[1:N])
  b <- theta.star[(N + 1):(2 * N)]
  sigma <- exp(theta.star[(2 * N + 1):(3 * N)])
  
  ### Initialize Gamma matrix
  Gamma <- diag(N)
  off_diag_elements <- exp(theta.star[(3 * N + 1):(length(theta.star))])
  Gamma[!Gamma] <- off_diag_elements
  Gamma <- Gamma / rowSums(Gamma)
  
  ### Compute the stationary distribution delta
  #eig <- eigen(t(Gamma))
  #delta <- Re(eig$vectors[,1])
  #delta <- delta / sum(delta)
  delta <- rep(1 / N, times = N)
  
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
  
  return(-l)
}


### Function to fit the model
fit_HMM <- function(r, N, theta.star) {
  result <- nlm(mllk, theta.star, x = r, N = N, print.level = 0, iterlim = 10000)
  return(c(exp(result$estimate[1:N]), ### kappa
           result$estimate[(N + 1):(2 * N)], ### theta
           exp(result$estimate[(2 * N + 1):(3 * N)]), ### sigma
           exp(result$estimate[(3 * N + 1):(length(theta.star))]))) ### Gamma
}





### Number of simulations
num_simulations <- 100

### Initial parameters for fitting
theta.star <- c(log(kappa_true),
                theta_true,
                log(sigma_true),
                log(Gamma_true[row(Gamma_true) != col(Gamma_true)]))

### List to store the results
fitted_params_HMM <- matrix(NA, nrow = num_simulations, ncol = length(theta.star))



### Loop over simulations
for (i in 1:num_simulations) {
  set.seed(i)  ### For reproducibility
  print(i)
  
  ### Simulate the Vasicek process
  simulation <- simulate_HMM(n, N, kappa_true, theta_true, sigma_true, Gamma_true, delta_true, r0)
  r <- simulation$r
  
  ### Fit the Vasicek model to the simulated series
  fitted_params_HMM[i, ] <- fit_HMM(r, N, theta.star)
}

### Print the fitted parameters
print(fitted_params_HMM)

HMMFit <- numeric(35)
for (i in 1:35){
  HMMFit[i] <- mean(fitted_params_HMM[,i])
}


Gamma5Sim <- diag(N)
Gamma5Sim[!Gamma5Sim] <- HMMFit[(3 * N + 1):(length(theta.star))]
Gamma5Sim <- Gamma5Sim / rowSums(Gamma5Sim)

solve(t(diag(N) - Gamma5Sim + 1), rep(1, N))









####################################################################################################################
#------------------------------------------ Load RData -------------------------------------------------------------
####################################################################################################################




setwd("~/Documents/KU/Bielefeld/RCode")
load("fitted_params_HMM.RData")




### kappa
### 0.002, 0.004, 0.006, 0.008, 0.01
### 1.828782e-03 4.283450e-03 6.660169e-03 9.274132e-03 9.963341e-03


### theta
### 0.02, 0.04, 0.06, 0.08, 0.1
### 2.945431e-02 4.198366e-02 5.887665e-02 7.469502e-02 9.319463e-02

### sigma
### 0.0002, 0.0004, 0.0006, 0.0008, 0.001
### 2.412948e-04 3.473457e-04 4.812295e-04 5.840045e-04 9.787518e-04

### Gamma
### 0.998, 0.002,   0,         0,       0,
### 0.001, 0.998,   0.001,     0,       0,
### 0,     0.001,   0.998,     0.001,   0,
### 0,     0,       0.001,     0.998,   0.001,
### 0,     0,       0,         0.002,   0.998
### NA             1.955823e-03   1.993559e-04   4.794324e-05   4.238754e-05
### 1.131943e-03   NA             1.546756e-03   4.122288e-05   4.122288e-05
### 1.835459e-04   1.534977e-03   NA             2.371802e-03   4.105503e-05
### 4.896317e-05   4.088037e-05   2.265967e-03   NA             2.135142e-03
### 4.101800e-05   4.101800e-05   4.101804e-05   3.146751e-03   NA



fitted_params_HMM


colMeans(fitted_params_HMM)
colSds(fitted_params_HMM)


kappa_true <- c(0.002, 0.004, 0.006, 0.008, 0.01)
theta_true <- c(0.02, 0.04, 0.06, 0.08, 0.1)
sigma_true <- c(0.0002, 0.0004, 0.0006, 0.0008, 0.001)
Gamma_true <- matrix(c(0.998, 0.002, 0, 0, 0,
                       0.001, 0.998, 0.001, 0, 0,
                       0, 0.001, 0.998, 0.001, 0,
                       0, 0, 0.001, 0.998, 0.001,
                       0, 0, 0, 0.002, 0.998),
                     byrow = TRUE, ncol = 5)


### Plot the fitted parameters
SimHMM <- as.data.frame(fitted_params_HMM)
colnames(SimHMM) <- c("kappa1", "kappa2", "kappa3", "kappa4", "kappa5", 
                      "theta1", "theta2", "theta3", "theta4", "theta5", 
                      "sigma1","sigma2","sigma3","sigma4","sigma5",
                      "P(2->1)", "P(3->1)", "P(4->1)", "P(5->1)",
                      "P(1->2)", "P(3->2)", "P(4->2)", "P(5->2)",
                      "P(1->3)", "P(2->3)", "P(4->3)", "P(5->3)",
                      "P(1->4)", "P(2->4)", "P(3->4)", "P(5->4)",
                      "P(1->5)", "P(2->5)", "P(3->5)", "P(4->5)")
SimHMM$Index <- c(1:num_simulations)



#------------------------------------------ kappa -------------------------------------------------------------



param_names <- c("kappa1", "kappa2", "kappa3", "kappa4", "kappa5")
SimHMM_kappa <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(SimHMM$kappa1,
            SimHMM$kappa2,
            SimHMM$kappa3,
            SimHMM$kappa4,
            SimHMM$kappa5)
)
true_kappa <- data.frame(
  Group = c("kappa1", "kappa2", "kappa3", "kappa4", "kappa5"),
  TrueValue = kappa_true
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(kappa[1]), expression(kappa[2]), expression(kappa[3]), expression(kappa[4]), expression(kappa[5]))

ggplot(SimHMM_kappa, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 1) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 5, show.legend = FALSE) +
  scale_color_manual(values=c("#CC5500", "#404080", "#87A96B", "#FFBCD9", "#996666")) +
  geom_hline(yintercept = kappa_true[1], col = "#CC5500", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = kappa_true[2], col = "#404080", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = kappa_true[3], col = "#87A96B", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = kappa_true[4], col = "#FFBCD9", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = kappa_true[5], col = "#996666", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_kappa, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = expression("Estimated speed of mean reversion, " * kappa),
       x = "",
       y = "") +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = latex_labels) + 
  theme(plot.margin=unit(c(.01,.3,-1.3,-.7),"cm")) ### c(top,right,bottom,left)






#------------------------------------------ sigma -------------------------------------------------------------



param_names <- c("sigma1", "sigma2", "sigma3", "sigma4", "sigma5")
SimHMM_sigma <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(SimHMM$sigma1,
            SimHMM$sigma2,
            SimHMM$sigma3,
            SimHMM$sigma4,
            SimHMM$sigma5)
)
true_sigma <- data.frame(
  Group = c("sigma1", "sigma2", "sigma3", "sigma4", "sigma5"),
  TrueValue = sigma_true
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(sigma[1]), expression(sigma[2]), expression(sigma[3]), expression(sigma[4]), expression(sigma[5]))

ggplot(SimHMM_sigma, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 1) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 5, show.legend = FALSE) +
  scale_color_manual(values=c("#CC5500", "#404080", "#87A96B", "#FFBCD9", "#996666")) +
  geom_hline(yintercept = sigma_true[1], col = "#CC5500", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = sigma_true[2], col = "#404080", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = sigma_true[3], col = "#87A96B", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = sigma_true[4], col = "#FFBCD9", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = sigma_true[5], col = "#996666", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_sigma, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = expression("Estimated volatility, " * sigma),
       x = "",
       y = "") +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = latex_labels) + 
  theme(plot.margin=unit(c(.01,.3,-1.3,-.7),"cm")) ### c(top,right,bottom,left)





#------------------------------------------ theta -------------------------------------------------------------



param_names <- c("theta1", "theta2", "theta3", "theta4", "theta5")
SimHMM_theta <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(SimHMM$theta1,
            SimHMM$theta2,
            SimHMM$theta3,
            SimHMM$theta4,
            SimHMM$theta5)
)
true_theta <- data.frame(
  Group = c("theta1", "theta2", "theta3", "theta4", "theta5"),
  TrueValue = theta_true
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(theta[1]), expression(theta[2]), expression(theta[3]), expression(theta[4]), expression(theta[5]))

ggplot(SimHMM_theta, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 1) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 5, show.legend = FALSE) +
  scale_color_manual(values=c("#CC5500", "#404080", "#87A96B", "#FFBCD9", "#996666")) +
  geom_hline(yintercept = theta_true[1], col = "#CC5500", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = theta_true[2], col = "#404080", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = theta_true[3], col = "#87A96B", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = theta_true[4], col = "#FFBCD9", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = theta_true[5], col = "#996666", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_theta, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = expression("Estimated long-run mean, " * theta),
       x = "",
       y = "") +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = latex_labels) + 
  theme(plot.margin=unit(c(.01,.3,-1.3,-.7),"cm")) ### c(top,right,bottom,left)






#------------------------------------------ Gamma -------------------------------------------------------------

### Structure
### "P(1->2)", "P(1->3)", "P(1->4)", "P(1->5)"
### "P(2->1)", "P(2->3)", "P(2->4)", "P(2->5)"
### "P(3->1)", "P(3->2)", "P(3->4)", "P(3->5)"
### "P(4->1)", "P(4->2)", "P(4->3)", "P(4->5)"
### "P(5->1)", "P(5->2)", "P(5->3)", "P(5->4)"



### Gamma1

param_names <- c("P(1->1)", "P(1->2)", "P(1->3)", "P(1->4)", "P(1->5)")
SimHMM_Gamma1 <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(rep(NA, 100),
            as.numeric(SimHMM$"P(1->2)"),
            as.numeric(SimHMM$"P(1->3)"),
            as.numeric(SimHMM$"P(1->4)"),
            as.numeric(SimHMM$"P(1->5)"))
)
true_Gamma1 <- data.frame(
  Group = c("P(1->1)", "P(1->2)", "P(1->3)", "P(1->4)", "P(1->5)"),
  TrueValue = c(NA, Gamma_true[1,2], Gamma_true[1,3], Gamma_true[1,4], Gamma_true[1,5])
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(gamma[11]), expression(gamma[12]), expression(gamma[13]), expression(gamma[14]), expression(gamma[15]))

HMM_Gamma1 <- ggplot(SimHMM_Gamma1, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 1) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 5, show.legend = FALSE) +
  scale_color_manual(values=c("#CC5500", "#404080", "#87A96B", "#FFBCD9", "#996666")) +
  geom_hline(yintercept = Gamma_true[1,2], col = "#404080", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[1,3], col = "#87A96B", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[1,4], col = "#FFBCD9", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[1,5], col = "#996666", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_Gamma1, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = expression("Estimated transition probability matrix, " * Gamma),
       x = "",
       y = expression(gamma["1j"])) +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = NULL) +
  coord_cartesian(ylim = NULL) + 
  theme(plot.margin=unit(c(.01,.3,-1.3,.3),"cm")) ### c(bottom,left,top,right)


### Gamma2

param_names <- c("P(2->1)", "P(2->2)", "P(2->3)", "P(2->4)", "P(2->5)")
SimHMM_Gamma2 <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(as.numeric(SimHMM$"P(2->1)"),
            rep(NA, 100),
            as.numeric(SimHMM$"P(2->3)"),
            as.numeric(SimHMM$"P(2->4)"),
            as.numeric(SimHMM$"P(2->5)"))
)
true_Gamma2 <- data.frame(
  Group = c("P(2->1)", "P(2->2)", "P(2->3)", "P(2->4)", "P(2->5)"),
  TrueValue = c(Gamma_true[2,1], NA, Gamma_true[2,3], Gamma_true[2,4], Gamma_true[2,5])
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(gamma[21]), expression(gamma[22]), expression(gamma[23]), expression(gamma[24]), expression(gamma[25]))

HMM_Gamma2 <- ggplot(SimHMM_Gamma2, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 1) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 5, show.legend = FALSE) +
  scale_color_manual(values=c("#CC5500", "#404080", "#87A96B", "#FFBCD9", "#996666")) +
  geom_hline(yintercept = Gamma_true[2,1], col = "#CC5500", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[2,3], col = "#87A96B", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[2,4], col = "#FFBCD9", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[2,5], col = "#996666", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_Gamma2, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = "",
       x = "",
       y = expression(gamma["2j"])) +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = NULL) +
  coord_cartesian(ylim = NULL) + 
  theme(plot.margin=unit(c(-1.3,.3,-1.3,.3),"cm")) ### c(bottom,left,top,right)







### Gamma3

param_names <- c("P(3->1)", "P(3->2)", "P(3->3)", "P(3->4)", "P(3->5)")
SimHMM_Gamma3 <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(as.numeric(SimHMM$"P(3->1)"),
            as.numeric(SimHMM$"P(3->2)"),
            rep(NA, 100),
            as.numeric(SimHMM$"P(3->4)"),
            as.numeric(SimHMM$"P(3->5)"))
)
true_Gamma3 <- data.frame(
  Group = c("P(3->1)", "P(3->2)", "P(3->3)", "P(3->4)", "P(3->5)"),
  TrueValue = c(Gamma_true[3,1], Gamma_true[3,2], NA, Gamma_true[3,4], Gamma_true[3,5])
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(gamma[31]), expression(gamma[32]), expression(gamma[33]), expression(gamma[34]), expression(gamma[35]))

HMM_Gamma3 <- ggplot(SimHMM_Gamma3, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 1) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 5, show.legend = FALSE) +
  scale_color_manual(values=c("#CC5500", "#404080", "#87A96B", "#FFBCD9", "#996666")) +
  geom_hline(yintercept = Gamma_true[3,1], col = "#CC5500", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[3,2], col = "#404080", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[3,4], col = "#FFBCD9", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[3,5], col = "#996666", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_Gamma3, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = "",
       x = "",
       y = expression(gamma["3j"])) +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = NULL) +
  coord_cartesian(ylim = NULL) + 
  theme(plot.margin=unit(c(-1.3,.3,-1.3,.3),"cm")) ### c(bottom,left,top,right)








### Gamma4

param_names <- c("P(4->1)", "P(4->2)", "P(4->3)", "P(4->4)", "P(4->5)")
SimHMM_Gamma4 <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(as.numeric(SimHMM$"P(4->1)"),
            as.numeric(SimHMM$"P(4->2)"),
            as.numeric(SimHMM$"P(4->3)"),
            rep(NA, 100),
            as.numeric(SimHMM$"P(4->5)"))
)
true_Gamma4 <- data.frame(
  Group = c("P(4->1)", "P(4->2)", "P(4->3)", "P(4->4)", "P(4->5)"),
  TrueValue = c(Gamma_true[4,1], Gamma_true[4,2], Gamma_true[4,3], NA, Gamma_true[4,5])
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(gamma[41]), expression(gamma[42]), expression(gamma[43]), expression(gamma[44]), expression(gamma[45]))

HMM_Gamma4 <- ggplot(SimHMM_Gamma4, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 1) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 5, show.legend = FALSE) +
  scale_color_manual(values=c("#CC5500", "#404080", "#87A96B", "#FFBCD9", "#996666")) +
  geom_hline(yintercept = Gamma_true[4,1], col = "#CC5500", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[4,2], col = "#404080", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[4,3], col = "#87A96B", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[4,5], col = "#996666", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_Gamma4, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = "",
       x = "",
       y = expression(gamma["4j"])) +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = NULL) +
  coord_cartesian(ylim = NULL) + 
  theme(plot.margin=unit(c(-1.3,.3,-1.3,.3),"cm")) ### c(bottom,left,top,right)










### Gamma5

param_names <- c("P(5->1)", "P(5->2)", "P(5->3)", "P(5->4)", "P(5->5)")
SimHMM_Gamma5 <- data.frame(
  Parameter = rep(param_names, each = 100),
  Value = c(as.numeric(SimHMM$"P(5->1)"),
            as.numeric(SimHMM$"P(5->2)"),
            as.numeric(SimHMM$"P(5->3)"),
            as.numeric(SimHMM$"P(5->4)"),
            rep(NA, 100))
)
true_Gamma5 <- data.frame(
  Group = c("P(5->1)", "P(5->2)", "P(5->3)", "P(5->4)", "P(5->5)"),
  TrueValue = c(Gamma_true[5,1], Gamma_true[5,2], Gamma_true[5,3], Gamma_true[5,4], NA)
)

# Define LaTeX-style labels for the x-axis
latex_labels <- c(expression(gamma[51]), expression(gamma[52]), expression(gamma[53]), expression(gamma[54]), expression(gamma[55]))

HMM_Gamma5 <- ggplot(SimHMM_Gamma5, aes(x = Parameter, y = Value)) +
  geom_boxplot(aes(color = factor(Parameter)), show.legend = FALSE, alpha = 0.7, size = 1) +
  geom_jitter(aes(color = Parameter), width = 0.2, alpha = 0.7, size = 5, show.legend = FALSE) +
  scale_color_manual(values=c("#CC5500", "#404080", "#87A96B", "#FFBCD9", "#996666")) +
  geom_hline(yintercept = Gamma_true[5,1], col = "#CC5500", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[5,2], col = "#404080", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[5,3], col = "#87A96B", linewidth = 0.7, linetype = "longdash") +
  geom_hline(yintercept = Gamma_true[5,4], col = "#FFBCD9", linewidth = 0.7, linetype = "longdash") +
  geom_point(data = true_Gamma5, aes(x = Group, y = TrueValue), color = "black", size = 7, shape = 18) +
  labs(title = "",
       x = "",
       y = expression(gamma["5j"])) +
  theme_bw() +
  theme(plot.title = element_text(size=47, hjust=0)) +
  theme(axis.title = element_text(size=37)) +
  theme(axis.text.x = element_text(size=37),
        axis.text.y = element_text(size=23)) +
  scale_x_discrete(labels = c("j=1", "j=2", "j=3", "j=4", "j=5")) +
  coord_cartesian(ylim = NULL) + 
  theme(plot.margin=unit(c(-1.5,.3,-1.3,.3),"cm")) ### c(bottom,left,top,right)



### Collected plot

grid.arrange(HMM_Gamma1,
             HMM_Gamma2,
             ncol = 1)
grid.arrange(HMM_Gamma3,
             HMM_Gamma4,
             ncol = 1)
grid.arrange(HMM_Gamma5,
             HMM_Gamma1,
             ncol = 1)



