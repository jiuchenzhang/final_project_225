library(ggplot2)
library(latex2exp)
library(reshape2)

############# Simulation 1 ###############
# Reproduce the simulation 1 (figure 1)  #
# sample from the target distribution    #
# U(\theta) = -2\theta^2 + \theta^4      #
##########################################

## common parameters used in all comparing algorithms
epsilon <- 0.1
leapfrog <- 50
iter <- 80000 # total samples
warmup <- 10000 # warmup samples

theta_collection <- matrix(rep(0, 5 * (iter)), ncol = 5) # collect all samples generated from 5 comparing algorithms

## Hamiltonian Monte Carlo

theta <- rep(0, iter) # parameter of interest
r <- rep(0, iter) # momentum

theta[1] <- 0

for (i in 1:(iter - 1)){
  
  r[i] <- rnorm(n=1)
  
  theta_ <- theta[i]
  r_ <- r[i]
  
  dev <- - 4 * theta_ * (1 - theta_^2)
  r_ <- r_ - epsilon / 2 * dev
  
  for (m in 1: leapfrog){
    theta_ <- theta_ + epsilon * r_
    dev <- - 4 * theta_ * (1 - theta_^2)
    r_ <- r_ - epsilon * dev
  }
  dev <- - 4 * theta_ * (1 - theta_^2)
  r_ <- r_ - epsilon / 2 * dev
  
  H1 <- - 2 * theta_ ^ 2 + theta_ ^ 4 + 1 / 2 * r_ ^ 2
  H2 <- - 2 * theta[i] ^ 2 + theta[i] ^ 4 + 1 / 2 * r[i] ^ 2 
  rho <- exp(H1 - H2)
  u <- runif(n=1)
  if (u < min(1, rho)){
    theta[i+1] <- theta_
  } else {
    theta[i+1] <- theta[i]
  }
}

theta_collection[, 1] <- theta[1: iter]


## Hamiltonian Monte Carlo without MH correction

theta <- rep(0, iter) # parameter of interest
r <- rep(0, iter) # momentum

theta[1] <- 0

for (i in 1:(iter - 1)){
  
    
  r[i] <- rnorm(n=1)
  
  theta_ <- theta[i]
  r_ <- r[i]
  
  dev <- - 4 * theta_ * (1 - theta_^2)
  r_ <- r_ - epsilon / 2 * dev
  
  for (m in 1: leapfrog){
    theta_ <- theta_ + epsilon * r_
    dev <- - 4 * theta_ * (1 - theta_^2)
    r_ <- r_ - epsilon * dev
  }
  theta[i+1] <- theta_
}


theta_collection[, 2] <- theta[1: iter]


## Naive Stochastic Gradient HMC

theta <- rep(0, iter) # parameter of interest
r <- rep(0, iter) # momentum

theta[1] <- 0

for (i in 1:(iter - 1)){
  
  r[i] <- rnorm(n=1)
  
  theta_ <- theta[i]
  r_ <- r[i]
  
  dev <- - 4 * theta_ * (1 - theta_^2) + rnorm(n=1) * 2
  r_ <- r_ - epsilon / 2 * dev
  
  for (m in 1: leapfrog){
    theta_ <- theta_ + epsilon * r_
    dev <- - 4 * theta_ * (1 - theta_^2) + rnorm(n=1) * 2
    r_ <- r_ - epsilon * dev
  }
  dev <- - 4 * theta_ * (1 - theta_^2) + rnorm(n=1) * 2
  r_ <- r_ - epsilon / 2 * dev
  
  H1 <- - 2 * theta_ ^ 2 + theta_ ^ 4 + 1 / 2 * r_ ^ 2
  H2 <- - 2 * theta[i] ^ 2 + theta[i] ^ 4 + 1 / 2 * r[i] ^ 2 
  rho <- exp(H1 - H2)
  u <- runif(n=1)
  if (u < min(1, rho)){
    theta[i+1] <- theta_
  } else {
    theta[i+1] <- theta[i]
  }
}

theta_collection[, 3] <- theta[1: iter]



## Naive Stochastic Gradient HMC without MH correction

set.seed(225)

theta <- rep(0, iter) # parameter of interest
r <- rep(0, iter) # momentum

theta[1] <- 0

for (i in 1:(iter - 1)){
  
    
  r[i] <- rnorm(n=1)
  
  theta_ <- theta[i]
  r_ <- r[i]
  
  dev <- - 4 * theta_ * (1 - theta_^2) + rnorm(n=1) * 2
  r_ <- r_ - epsilon / 2 * dev
  
  for (m in 1: leapfrog){
    theta_ <- theta_ + epsilon * r_
    dev <- - 4 * theta_ * (1 - theta_^2) + rnorm(n=1) * 2
    r_ <- r_ - epsilon * dev
  }
  theta[i+1] <- theta_
}


theta_collection[, 4] <- theta[1: iter]

## Stochastic Gradient HMC with Friction
C <- 3
alpha <- epsilon * C
V <- 4

set.seed(225)

iter <- 80000 # total samples
warmup <- 10000 # warmup samples

theta <- rep(0, iter) # parameter of interest
r <- rep(0, iter) # momentum
beta <- V * epsilon^2 / 2

theta[1] <- 0

for (i in 1:(iter - 1)){
  
  r[i] <- rnorm(n=1)
  
  theta_ <- theta[i]
  r_ <- r[i]
  
  dev <- - 4 * theta_ * (1 - theta_^2) + rnorm(n=1) * 2
  r_ <- r_ * epsilon
  sig <- 2 * epsilon^2 * (alpha - beta)
  
  for (m in 1: leapfrog){
    theta_ <- theta_ + r_
    dev <- - 4 * theta_ * (1 - theta_^2) + rnorm(n=1) * 2
    r_ <- r_ * (1 - alpha) - dev * epsilon^2 + rnorm(n=1, mean=0, sd=sqrt(sig))
  }
  theta[i+1] <- theta_
}

theta_collection[, 5] <- theta[1: iter]

theta_collection <- data.frame(theta_collection)
colnames(theta_collection) <- c("HMC1", "HMC2", "NSGHMC1", "NSGHMC2", "SGHMC")
theta_collection <- melt(theta_collection, measure.vars=c("HMC1", "HMC2", "NSGHMC1", "NSGHMC2", "SGHMC")) # wide to long

ggplot(theta_collection, aes(x=value, colour=variable)) + geom_density() + 
  scale_colour_manual(labels=c("Standard HMC(with MH)", "Standard HMC(no MH)", "Naive stochastic gradient HMC(with MH)", 
                              "Naive stochastic gradient HMC(no MH)", "SGHMC"), values=c("red", "orange", "yellow", "green", "blue")) + 
  labs(x = TeX("$\\theta$"), y = "density", colour = "Algorithms", title=TeX(sprintf("$\\epsilon = %f$", epsilon))) 

ggsave(paste("simulation1_epsilon_", epsilon, ".jpg", sep=""), width = 20, height=12)
