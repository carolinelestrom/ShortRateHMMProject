
#------------------------------------------ theta ------------------------------------------------------------


### Define the log-likelihood function for the HMM with state-dependent normal distributions
mllk <- function(theta.star, x, N){
  a <- exp(theta.star[1])
  b <- theta.star[(N + 1):(2 * N)]
  sigma <- exp(theta.star[(2 * N + 1)])
  
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
    allprobs[ind, j] <- dnorm(x[ind], mean = b[j], sd = sqrt(sigma[1]^2 / (2 * a[1])))
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


### 2-state HMM
N <- 2
theta.star <- c(qlogis(c(0.9, 0.9)),
                log(c(0.0003, 0.0007)),
                0.01, 0.02,
                log(c(0.0005, 0.001))) 
mod2_theta <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)



### 3-state HMM
N <- 3
theta.star <- c(log(seq(0.0001, 0.0009, length = N)),
                seq(0.005, 0.08, length = N),
                log(seq(0.0001, 0.0015, length = N)),
                rep(-2, (N - 1) * N))
mod3_theta <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)

### 4-state HMM
N <- 4
theta.star <- c(log(seq(0.0001, 0.0009, length = N)),
                seq(0.005, 0.08, length = N),
                log(seq(0.0001, 0.0015, length = N)),
                rep(-2, (N - 1) * N))

mod4_theta <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)


### 5-state HMM
N <- 5
theta.star <- c(log(seq(0.0005, 0.001, length = N)),
                seq(0.03, 0.08, length = N),
                log(seq(0.0001, 0.001, length = N)),
                rep(-5, (N - 1) * N))



mod5_theta <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)






#------------------------------------------ kappa ------------------------------------------------------------


### Define the log-likelihood function for the HMM with state-dependent normal distributions
mllk <- function(theta.star, x, N){
  a <- exp(theta.star[1:N])
  b <- theta.star[(N + 1)]
  sigma <- exp(theta.star[(2 * N + 1)])
  
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
    allprobs[ind, j] <- dnorm(x[ind], mean = b[1], sd = sqrt(sigma[1]^2 / (2 * a[j])))
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


### 2-state HMM
N <- 2
theta.star <- c(qlogis(c(0.9, 0.9)),
                log(c(0.0003, 0.0007)),
                0.01, 0.02,
                log(c(0.0005, 0.001))) 
mod2_kappa <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)


### 3-state HMM
N <- 3
theta.star <- c(log(seq(0.0001, 0.0009, length = N)),
                seq(0.005, 0.08, length = N),
                log(seq(0.0001, 0.0015, length = N)),
                rep(-2, (N - 1) * N))
mod3_kappa <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)

### 4-state HMM
N <- 4
theta.star <- c(log(seq(0.0001, 0.0009, length = N)),
                seq(0.005, 0.08, length = N),
                log(seq(0.0001, 0.0015, length = N)),
                rep(-2, (N - 1) * N))

mod4_kappa <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)


### 5-state HMM
N <- 5
theta.star <- c(log(seq(0.0005, 0.001, length = N)),
                seq(0.03, 0.08, length = N),
                log(seq(0.0001, 0.001, length = N)),
                rep(-5, (N - 1) * N))



mod5_kappa <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)







#------------------------------------------ sigma ------------------------------------------------------------


### Define the log-likelihood function for the HMM with state-dependent normal distributions
mllk <- function(theta.star, x, N){
  a <- exp(theta.star[1])
  b <- theta.star[(N + 1)]
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
  delta <- rep(1 / N, times = N)
  
  ### Initialize allprobs matrix
  allprobs <- matrix(1, length(x), N)
  ind <- which(!is.na(x))
  
  for (j in 1:N){
    allprobs[ind, j] <- dnorm(x[ind], mean = b[1], sd = sqrt(sigma[j]^2 / (2 * a[1])))
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


### 2-state HMM
N <- 2
theta.star <- c(qlogis(c(0.9, 0.9)),
                log(c(0.0003, 0.0007)),
                0.01, 0.02,
                log(c(0.0005, 0.001))) 
mod2_sigma <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)


### 3-state HMM
N <- 3
theta.star <- c(log(seq(0.0001, 0.0009, length = N)),
                seq(0.005, 0.08, length = N),
                log(seq(0.0001, 0.0015, length = N)),
                rep(-2, (N - 1) * N))
mod3_sigma <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)

### 4-state HMM
N <- 4
theta.star <- c(log(seq(0.0001, 0.0009, length = N)),
                seq(0.005, 0.08, length = N),
                log(seq(0.0001, 0.0015, length = N)),
                rep(-2, (N - 1) * N))

mod4_sigma <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)


### 5-state HMM
N <- 5
theta.star <- c(log(seq(0.0005, 0.001, length = N)),
                seq(0.03, 0.08, length = N),
                log(seq(0.0001, 0.001, length = N)),
                rep(-5, (N - 1) * N))



mod5_sigma <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)





#------------------------------------------ sigma & theta ------------------------------------------------------------

### Define the log-likelihood function for the HMM with state-dependent normal distributions
mllk <- function(theta.star, x, N){
  a <- exp(theta.star[1])
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
    allprobs[ind, j] <- dnorm(x[ind], mean = b[j], sd = sqrt(sigma[j]^2 / (2 * a[1])))
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


### 2-state HMM
N <- 2
theta.star <- c(qlogis(c(0.9, 0.9)),
                log(c(0.0003, 0.0007)),
                0.01, 0.02,
                log(c(0.0005, 0.001))) 
mod2_sigma_theta <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)


### 3-state HMM
N <- 3
theta.star <- c(log(seq(0.0001, 0.0009, length = N)),
                seq(0.005, 0.08, length = N),
                log(seq(0.0001, 0.0015, length = N)),
                rep(-2, (N - 1) * N))
mod3_sigma_theta <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)

### 4-state HMM
N <- 4
theta.star <- c(log(seq(0.0001, 0.0009, length = N)),
                seq(0.005, 0.08, length = N),
                log(seq(0.0001, 0.0015, length = N)),
                rep(-2, (N - 1) * N))

mod4_sigma_theta <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)


### 5-state HMM
N <- 5
theta.star <- c(log(seq(0.0005, 0.001, length = N)),
                seq(0.03, 0.08, length = N),
                log(seq(0.0001, 0.001, length = N)),
                rep(-5, (N - 1) * N))



mod5_sigma_theta <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)


exp(mod5_sigma_theta$estimate[1:N])
mod5_sigma_theta$estimate[(N + 1):(2 * N)] 
exp(mod5_sigma_theta$estimate[(2 * N + 1):(3 * N)])

#------------------------------------------ kappa & theta ------------------------------------------------------------

### Define the log-likelihood function for the HMM with state-dependent normal distributions
mllk <- function(theta.star, x, N){
  a <- exp(theta.star[1:N])
  b <- theta.star[(N + 1):(2 * N)]
  sigma <- exp(theta.star[(2 * N + 1)])
  
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
    allprobs[ind, j] <- dnorm(x[ind], mean = b[j], sd = sqrt(sigma[1]^2 / (2 * a[j])))
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





### 2-state HMM
N <- 2
theta.star <- c(qlogis(c(0.9, 0.9)),
                log(c(0.0003, 0.0007)),
                0.01, 0.02,
                log(c(0.0005, 0.001))) 
mod2_kappa_theta <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)


### 3-state HMM
N <- 3
theta.star <- c(log(seq(0.0001, 0.0009, length = N)),
                seq(0.005, 0.08, length = N),
                log(seq(0.0001, 0.0015, length = N)),
                rep(-2, (N - 1) * N))
mod3_kappa_theta <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)

### 4-state HMM
N <- 4
theta.star <- c(log(seq(0.0001, 0.0009, length = N)),
                seq(0.005, 0.08, length = N),
                log(seq(0.0001, 0.0015, length = N)),
                rep(-2, (N - 1) * N))

mod4_kappa_theta <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)


### 5-state HMM
N <- 5
theta.star <- c(log(seq(0.0005, 0.001, length = N)),
                seq(0.03, 0.08, length = N),
                log(seq(0.0001, 0.001, length = N)),
                rep(-5, (N - 1) * N))



mod5_kappa_theta <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)





#------------------------------------------ sigma & kappa ------------------------------------------------------------

### Define the log-likelihood function for the HMM with state-dependent normal distributions
mllk <- function(theta.star, x, N){
  a <- exp(theta.star[1:N])
  b <- theta.star[(N + 1)]
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
    allprobs[ind, j] <- dnorm(x[ind], mean = b[1], sd = sqrt(sigma[j]^2 / (2 * a[j])))
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





### 2-state HMM
N <- 2
theta.star <- c(qlogis(c(0.9, 0.9)),
                log(c(0.0003, 0.0007)),
                0.01, 0.02,
                log(c(0.0005, 0.001))) 
mod2_kappa_sigma <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)


### 3-state HMM
N <- 3
theta.star <- c(log(seq(0.0001, 0.0009, length = N)),
                seq(0.005, 0.08, length = N),
                log(seq(0.0001, 0.0015, length = N)),
                rep(-2, (N - 1) * N))
mod3_kappa_sigma <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)

### 4-state HMM
N <- 4
theta.star <- c(log(seq(0.0001, 0.0009, length = N)),
                seq(0.005, 0.08, length = N),
                log(seq(0.0001, 0.0015, length = N)),
                rep(-2, (N - 1) * N))

mod4_kappa_sigma <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)


### 5-state HMM
N <- 5
theta.star <- c(log(seq(0.0005, 0.001, length = N)),
                seq(0.03, 0.08, length = N),
                log(seq(0.0001, 0.001, length = N)),
                rep(-5, (N - 1) * N))



mod5_kappa_sigma <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)


#------------------------------------------ sigma & kappa & theta ------------------------------------------------------------



### Define the log-likelihood function for the HMM with state-dependent normal distributions
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
  
  return(-l)
}

### 2-state HMM
N <- 2
theta.star <- c(qlogis(c(0.9, 0.9)),
                log(c(0.0003, 0.0007)),
                0.01, 0.02,
                log(c(0.0005, 0.001))) 
mod2 <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)

Gamma2_hat <- diag(plogis(mod2$estimate[1:2]))
Gamma2_hat[1,2] <- 1 - Gamma2_hat[1,1]
Gamma2_hat[2,1] <- 1 - Gamma2_hat[2,2]
###             [,1]         [,2]
### [1,] 0.999371654 0.0006283459
### [2,] 0.000468105 0.9995318950

delta2_hat <- solve(t(diag(2) - Gamma2_hat + 1), c(1, 1))
### 0.4269275 0.5730725


a2_hat <- exp(mod2$estimate[3:4]) ### 0.0004661632, 0.0007974269
b2_hat <- mod2$estimate[5:6] ### 0.00591042, 0.05428843
sigma2_hat <- exp(mod2$estimate[7:8]) ### 0.0002044272, 0.0007557500






### 3-state HMM
N <- 3
theta.star <- c(log(seq(0.0001, 0.0009, length = N)),
                seq(0.005, 0.08, length = N),
                log(seq(0.0001, 0.0015, length = N)),
                rep(-2, (N - 1) * N))
mod3 <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)



Gamma3_hat <- diag(N)
Gamma3_hat[!Gamma3_hat] <- exp(mod3$estimate[(3 * N + 1):(length(theta.star))])
Gamma3_hat <- Gamma3_hat / rowSums(Gamma3_hat)
###           [,1]      [,2]      [,3]
### [1,] 0.9981361 0.0018639 0.0000000
### [2,] 0.0017134 0.9976250 0.0006616
### [3,] 0.0000000 0.0002989 0.9997011

delta3_hat <- solve(t(diag(N) - Gamma3_hat + 1), rep(1, N))
### 0.2224343 0.2419723 0.5355934

a3_hat <- exp(mod3$estimate[1:N]) ### 0.0002516784 0.0008336284 0.0010366643
b3_hat <- mod3$estimate[(N + 1):(2 * N)] ### 0.0008120922 0.0139858111 0.0564406363
sigma3_hat <- exp(mod3$estimate[(2 * N + 1):(3 * N)]) ### 0.0000130 0.0002637 0.0008044






### 4-state HMM
N <- 4
theta.star <- c(log(seq(0.0005, 0.001, length = N)),
                seq(0.008, 0.07, length = N),
                log(seq(0.0005, 0.0015, length = N)),
                rep(-5, (N - 1) * N))
theta.star <- c(mod5$estimate[1:2], mod5$estimate[4:5],
                mod5$estimate[6:7], mod5$estimate[9:10],
                mod5$estimate[11:12], mod5$estimate[14:15],
                rep(-5, (N - 1) * N))

mod4 <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)

Gamma4_hat <- diag(N)
Gamma4_hat[!Gamma4_hat] <- exp(mod4$estimate[(3 * N + 1):(length(theta.star))])
Gamma4_hat <- Gamma4_hat / rowSums(Gamma4_hat)
###           [,1]      [,2]      [,3]      [,4]
### [1,] 0.9977881 0.0000487 0.0000520 0.0021112
### [2,] 0.0000116 0.9971555 0.0028202 0.0000127
### [3,] 0.0000269 0.0024303 0.9970804 0.0004624
### [4,] 0.0021344 0.0000503 0.0002800 0.9975353

delta4_hat <- solve(t(diag(N) - Gamma4_hat + 1), rep(1, N))
### 0.2742747 0.2105076 0.2350909 0.2801267

a4_hat <- exp(mod4$estimate[1:N]) ### 0.0002697959 0.0006844986 0.0006806138 0.0017699180
b4_hat <- mod4$estimate[(N + 1):(2 * N)] ### 0.0008119765 0.0526266806 0.0599571591 0.0139380071
sigma4_hat <- exp(mod4$estimate[(2 * N + 1):(3 * N)]) ### 0.0000134234 0.0001598952 0.0008763991 0.0003819540




### 5-state HMM
N <- 5
#theta.star <- c(log(seq(0.0001, 0.0009, length = N)),
#                seq(0.005, 0.08, length = N),
#                log(seq(0.0001, 0.0015, length = N)),
#                rep(-2, (N - 1) * N))
theta.star <- c(log(seq(0.0005, 0.001, length = N)),
                seq(0.03, 0.08, length = N),
                log(seq(0.0001, 0.001, length = N)),
                rep(-5, (N - 1) * N))
theta.star <- mod5$estimate

mod5 <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N = N, print.level = 2, iterlim = 10000)


Gamma5_hat <- diag(N)
Gamma5_hat[!Gamma5_hat] <- exp(mod5$estimate[(3 * N + 1):(length(theta.star))])
Gamma5_hat <- Gamma5_hat / rowSums(Gamma5_hat)
###           [,1]      [,2]      [,3]      [,4]
### [1,] 0.9977881 0.0000487 0.0000520 0.0021112
### [2,] 0.0000116 0.9971555 0.0028202 0.0000127
### [3,] 0.0000269 0.0024303 0.9970804 0.0004624
### [4,] 0.0021344 0.0000503 0.0002800 0.9975353

delta5_hat <- solve(t(diag(N) - Gamma5_hat + 1), rep(1, N))
### 0.2742747 0.2105076 0.2350909 0.2801267

a5_hat <- exp(mod5$estimate[1:N]) ### 0.0002697959 0.0006844986 0.0006806138 0.0017699180
b5_hat <- mod5$estimate[(N + 1):(2 * N)] ### 0.0008119765 0.0526266806 0.0599571591 0.0139380071
sigma5_hat <- exp(mod5$estimate[(2 * N + 1):(3 * N)]) ### 0.0000134234 0.0001598952 0.0008763991 0.0003819540













####################################################################################################################
####################################################################################################################
#-------------------------- Robustness test for various starting values --------------------------------------------
####################################################################################################################
####################################################################################################################


llks <- rep(NA, 100)
mods <- list()
for (k in 1:100){
  theta.star <- c(log(runif(N, 0, 1000)), runif(N, 0, 1000), log(runif(N, 0, 1000)), rep(-2, (N - 1) * N)) 
  mods[[k]] <- nlm(mllk, theta.star, x = as.numeric(yields_df$"3M"), N=N, stepmax = 5)
  llks[k] <- -mods[[k]]$minimum
}
plot(llks)





