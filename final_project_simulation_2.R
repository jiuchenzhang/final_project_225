# load required packages

library(MASS)
library(mvtnorm)
library(hexbin)
library(latex2exp)
library(ggplot2)
library(reshape2)

############ Simulation 2 ############
# Sample posterior distribution of   #
# 2d normal model with normal prior  #
# test running time, convergence     #
# of several algorithms              #
######################################

### Prior
mu <- rep(0, 2)
Tau <- 10 * diag(2)

### Data
normal_data <- function(n, seed){

  set.seed(seed)
  Sigma <- matrix(c(1, 0.7, 0.7, 1), ncol=2)
  Y <- mvrnorm(n=n, c(5, 3), Sigma)

  out <- list(Y, Sigma)
  names(out) <- c("Y", "Sigma")
  
  return(out)
}

### Common parameters used in all comparing algorithms
iter <- 5000
warmup <- 1000
epsilon <- 0.001
delta <- 0.1
leapfrog <- 50
batch_size <- 5

## Multivariate Metropolis-Hastings algorithm
MH <- function(Y, mu, Tau, Sigma, iter, warmup){

  start_time <- Sys.time()

  theta <- matrix(rep(0, iter * 2), ncol = 2)

  theta[1, ] <- mvrnorm(1, mu, Tau)

  for (i in 2:iter){
    old_state <- theta[i-1,]
    new_state <- mvrnorm(1, mu = old_state, Sigma = delta * diag(2)) # propose new state
    a <- min(1, dmvnorm(new_state, mu, Tau) / dmvnorm(old_state, mu, Tau) * 
              (exp(sum(log(dmvnorm(Y, new_state, Sigma))) - sum(log(dmvnorm(Y, old_state, Sigma))))))
    if (runif(n=1) < a){theta[i,] <- new_state} else {theta[i,] <- old_state}
  }

  end_time <- Sys.time()

  out <- list(theta[(warmup + 1): iter, ], end_time - start_time)
  names(out) <- c("theta", "time")
  return(out)
}

## Hamiltonian Monte Carlo

HMC <- function(Y, mu, Tau, Sigma, iter, warmup, leapfrog, MHC){

  n <- nrow(Y)

  start_time <- Sys.time()
  theta <- matrix(rep(0, iter * 2), ncol = 2) # parameter of interest
  r <- matrix(rep(0, iter * 2), ncol = 2) # momentum

  theta[1, ] <- mvrnorm(1, mu, Tau)

  for (i in 1:(iter - 1)){
    
    r[i, ] <- mvrnorm(n=1, mu = mu, Sigma = diag(2))
    
    theta_ <- theta[i, ]
    r_ <- r[i, ]
    
    dev <- - solve(Sigma) %*% (colSums(Y) - n * theta_) + solve(Tau) %*% (theta_ - mu)
    r_ <- r_ - epsilon / 2 * dev
    
    for (m in 1: leapfrog){
      theta_ <- theta_ + epsilon * r_
      dev <- - solve(Sigma) %*% (colSums(Y) - n * theta_) + solve(Tau) %*% (theta_ - mu)
      r_ <- r_ - epsilon * dev
    }
    
    theta[i+1, ] <- theta_
    
    if (MHC){

      dev <- - solve(Sigma) %*% (colSums(Y) - n * theta_) + solve(Tau) %*% (theta_ - mu)
      r_ <- r_ - epsilon / 2 * dev
      H1 <- - sum(log(dmvnorm(Y, mean = theta_, sigma = Sigma))) - 
            log(dmvnorm(t(theta_), mean = mu, sigma = Tau)) + 1/2 * t(r_) %*% r_
      H2 <- - sum(log(dmvnorm(Y, mean = theta[i, ], sigma = Sigma))) -
              log(dmvnorm(theta[i, ], mean = mu, sigma = Tau)) + 1/2 * t(r[i, ]) %*% r[i, ]
      rho <- exp(H1 - H2)
      u <- runif(n=1)
      if (u > min(1, rho)){
        theta[i+1, ] <- theta[i, ]
      }
    }
    
  }

  end_time <- Sys.time()

  out <- list(theta[(warmup + 1): iter, ], end_time - start_time)
  names(out) <- c("theta", "time")

  return(out)

}


## Naive Stochastic Gradient HMC
NSGHMC <- function(Y, mu, Tau, Sigma, iter, warmup, leapfrog, batch_size, MHC){
  
  n <- nrow(Y)
  start_time <- Sys.time()
  
  theta <- matrix(rep(0, iter * 2), ncol = 2) # parameter of interest
  r <- matrix(rep(0, iter * 2), ncol = 2) # momentum
  
  theta[1, ] <- mvrnorm(1, mu, Tau)
  
  for (i in 1:(iter - 1)){
    
    r[i, ] <- mvrnorm(n=1, mu = mu, Sigma = diag(2))
    
    theta_ <- theta[i, ]
    r_ <- r[i, ]
    
    batch <- Y[sample(c(1: n), batch_size), ]
    dev <- - n / batch_size * solve(Sigma) %*% (colSums(batch) - batch_size * theta_) + solve(Tau) %*% (theta_ - mu)
    r_ <- r_ - epsilon / 2 * dev
    
    for (m in 1: leapfrog){
      theta_ <- theta_ + epsilon * r_
      batch <- Y[sample(c(1: n), batch_size), ]
      dev <- - n / batch_size * solve(Sigma) %*% (colSums(batch) - batch_size * theta_) + solve(Tau) %*% (theta_ - mu)
      r_ <- r_ - epsilon * dev
    }
    
    theta[i+1, ] <- theta_
    
    if (MHC){
      
      batch <- Y[sample(c(1: n), batch_size), ]
      dev <- - n / batch_size * solve(Sigma) %*% (colSums(batch) - batch_size * theta_) + solve(Tau) %*% (theta_ - mu)
      r_ <- r_ - epsilon / 2 * dev
      
      H1 <- - sum(log(dmvnorm(Y, mean = theta_, sigma = Sigma))) -
        log(dmvnorm(t(theta_), mean = mu, sigma = Tau)) + 1/2 * t(r_) %*% r_
      H2 <- - sum(log(dmvnorm(Y, mean = theta[i, ], sigma = Sigma))) -
        log(dmvnorm(theta[i, ], mean = mu, sigma = Tau)) + 1/2 * t(r[i, ]) %*% r[i, ]
      rho <- exp(H1 - H2)
      u <- runif(n=1)
      if (u > min(1, rho)){
        theta[i+1, ] <- theta[i, ]
      }
    }
    
  }
  
  end_time <- Sys.time()
  
  out <- list(theta[(warmup + 1): iter, ], end_time - start_time)
  names(out) <- c("theta", "time")
  
  return(out)
}


## Stochastic Gradient HMC with Friction

SGHMC <- function(Y, mu, Tau, Sigma, iter, warmup, leapfrog){
  
  alpha <- 0.03
  V <- 1
  n <- nrow(Y)
  
  start_time <- Sys.time()
  
  
  theta <- matrix(rep(0, iter * 2), ncol = 2) # parameter of interest
  r <- matrix(rep(0, iter * 2), ncol = 2) # momentum
  beta <- V * epsilon^2 * 0.5
  
  theta[1, ] <- mvrnorm(n=1, mu, Tau)
  
  for (i in 1:(iter - 1)){
    
    r[i, ] <- mvrnorm(n=1, mu = mu, Sigma = diag(2))
    
    theta_ <- theta[i, ]
    r_ <- r[i, ]
    
    batch <- Y[sample(c(1: n), batch_size), ]
    dev <- - n / batch_size * solve(Sigma) %*% (colSums(batch) - batch_size * theta_) + solve(Tau) %*% (theta_ - mu)
    r_ <- r_ * epsilon
    sig <- 2 * epsilon^2 * (alpha - beta)
    
    for (m in 1: leapfrog){
      batch <- Y[sample(c(1: n), batch_size), ]
      dev <- - n / batch_size * solve(Sigma) %*% (colSums(batch) - batch_size * theta_) + solve(Tau) %*% (theta_ - mu)
      r_ <- r_ * (1 - alpha) - dev * epsilon^2 + mvrnorm(n=1, mu = mu, Sigma = sig * diag(2))
      theta_ <- theta_ + r_
    }
    theta[i+1, ] <- theta_
  }
  
  end_time <- Sys.time()
  
  out <- list(theta[(warmup + 1): iter, ], end_time - start_time)
  names(out) <- c("theta", "time")
  
  return(out)
}

## comparing running time

sample_size <- c(100, 500, 1000, 10000)

time_collection <- matrix(rep(0, 6 * 4), ncol = 6)

for (s in 1: 4){
  
  data <- normal_data(sample_size[s], 225)
  
  res1 <- MH(data$Y, mu, Tau, data$Sigma, iter, warmup)
  res2 <- HMC(data$Y, mu, Tau, data$Sigma, iter, warmup, leapfrog, TRUE)
  res3 <- HMC(data$Y, mu, Tau, data$Sigma, iter, warmup, leapfrog, FALSE)
  res4 <- NSGHMC(data$Y, mu, Tau, data$Sigma, iter, warmup, leapfrog, batch_size, TRUE)
  res5 <- NSGHMC(data$Y, mu, Tau, data$Sigma, iter, warmup, leapfrog, batch_size, FALSE)
  res6 <- SGHMC(data$Y, mu, Tau, data$Sigma, iter, warmup, leapfrog)
  
  time_collection[s, 1] <- res1$time
  time_collection[s, 2] <- res2$time
  time_collection[s, 3] <- res3$time
  time_collection[s, 4] <- res4$time
  time_collection[s, 5] <- res5$time
  time_collection[s, 6] <- res6$time
  
}

## visualize the time comparison

time_collection <- cbind(time_collection, sample_size)
time_collection <- data.frame(time_collection)
colnames(time_collection) <- c("MH", "HMC1", "HMC2", 
                               "NSGHMC1",
                               "NSGHMC2",
                               "SGMHC", "Sample_size")
time_collection <- melt(time_collection, id.vars = "Sample_size", measure.vars = c("MH", "HMC1", "HMC2", 
                                                                                   "NSGHMC1",
                                                                                   "NSGHMC2",
                                                                                   "SGMHC"))

ggplot(time_collection, aes(x=log(Sample_size), y=value, colour=variable)) + geom_line() + 
  labs(y="Running time(seconds)") + 
  scale_colour_manual(labels=c("Metropolis Hastings", "Standard HMC(with MH)", "Standard HMC(no MH)", "Naive stochastic gradient HMC(with MH)", 
                               "Naive stochastic gradient HMC(no MH)", "SGHMC"), values=c("red", "orange", "yellow", "green", "blue", "purple"))

