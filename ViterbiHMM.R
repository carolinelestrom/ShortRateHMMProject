####################################################################################################################
####################################################################################################################
#----------------------------------------- 2-state Viterbi ---------------------------------------------------------
####################################################################################################################
####################################################################################################################



viterbi <- function(x, a, b, sigma, Gamma, delta){ 
  n <- length(x)
  allprobs <- matrix(1, nrow = n, ncol = 2)
  ind <- which(!is.na(x))
  allprobs[ind, ] <- cbind(dnorm(x[ind], mean = b[1], sd = sqrt(sigma[1]^2/(2*a[1]))),
                           dnorm(x[ind], mean = b[2], sd = sqrt(sigma[2]^2/(2*a[2]))))
  xi <- matrix(0, nrow = n, ncol = 2) 
  foo <- delta * allprobs[1, ]
  xi[1, ] <- foo / sum(foo)
  for (t in 2:n){
    foo <- apply(xi[t - 1, ] * Gamma, 2, max) * allprobs[t, ]
    xi[t, ] <- foo / sum(foo) 
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ]) 
  for (t in (n - 1):1){
    iv[t] <- which.max(Gamma[, iv[t + 1]] * xi[t, ]) 
  }
  iv 
}

viterbi2 <- viterbi(x = as.numeric(yields_df$"3M"), a2_hat, b2_hat, sigma2_hat, Gamma2_hat, delta2_hat)
plot(viterbi2)



####################################################################################################################
####################################################################################################################
#----------------------------------------- N-state Viterbi ---------------------------------------------------------
####################################################################################################################
####################################################################################################################



viterbi <- function(x, a, b, sigma, Gamma, delta, N){ 
  n <- length(x)
  allprobs <- matrix(1, nrow = n, ncol = N)
  ind <- which(!is.na(x))
  for (j in 1:N){
    allprobs[ind, j] <- dnorm(x[ind], mean = b[j], sd = sqrt(sigma[j]^2 / (2 * a[j])))
  }
  xi <- matrix(0, nrow = n, ncol = N) 
  foo <- delta * allprobs[1, ]
  xi[1, ] <- foo / sum(foo)
  for (t in 2:n){
    foo <- apply(xi[t - 1, ] * Gamma, 2, max) * allprobs[t, ]
    xi[t, ] <- foo / sum(foo) 
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ]) 
  for (t in (n - 1):1){
    iv[t] <- which.max(Gamma[, iv[t + 1]] * xi[t, ]) 
  }
  iv 
}


### 3-state HMM
N <- 3
viterbi3 <- viterbi(x = as.numeric(yields_df$"3M"), a3_hat, b3_hat, sigma3_hat, Gamma3_hat, delta3_hat, N=N)
plot(viterbi3)

### 4-state HMM
N <- 4
viterbi4 <- viterbi(x = as.numeric(yields_df$"3M"), a4_hat, b4_hat, sigma4_hat, Gamma4_hat, delta4_hat, N=N)
plot(viterbi4)

### 4-state HMM
N <- 5
viterbi5 <- viterbi(x = as.numeric(yields_df$"3M"), a5_hat, b5_hat, sigma5_hat, Gamma5_hat, delta5_hat, N=N)
plot(viterbi5)

