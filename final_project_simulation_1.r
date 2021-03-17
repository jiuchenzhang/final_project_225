
############# Simulation 1 ###############
# Reproduce the simulation 1 (figure 1)  #
##########################################


## Hamiltonian Monte Carlo

### sample from the target distribution U(\theta) = -2\theta^2 + \theta^4

epsilon <- 1e-1
leapfrog <- 50

iter <- 80000 # total samples
warmup <- 10000 # warmup samples

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

### density plot
d_hmc1 <- density(theta[(warmup + 1): iter], bw=0.1)
plot(d_hmc1)



## Hamiltonian Monte Carlo without MH correction

### sample from the target distribution U(\theta) = -2\theta^2 + \theta^4

epsilon <- 1e-1
leapfrog <- 50

iter <- 80000 # total samples
warmup <- 10000 # warmup samples

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

### density plot
d_hmc2 <- density(theta[(warmup + 1): iter], bw=0.1)
plot(d_hmc2)


## Naive Stochastic Gradient HMC
### sample from the target distribution U(\theta) = -2\theta^2 + \theta^4

epsilon <- 1e-1
leapfrog <- 50

iter <- 80000 # total samples
warmup <- 10000 # warmup samples

set.seed(225)

theta <- rep(0, iter) # parameter of interest
r <- rep(0, iter) # momentum

theta[1] <- 0

for (i in 1:(iter - 1)){
  
  r[i] <- rnorm(n=1)
  
  theta_ <- theta[i]
  r_ <- r[i]
  
  dev <- - 4 * theta_ * (1 - theta_^2) + rnorm(n=1, sd=2)
  r_ <- r_ - epsilon / 2 * dev
  
  for (m in 1: leapfrog){
    theta_ <- theta_ + epsilon * r_
    dev <- - 4 * theta_ * (1 - theta_^2) + rnorm(n=1, sd=2)
    r_ <- r_ - epsilon * dev
  }
  dev <- - 4 * theta_ * (1 - theta_^2) + rnorm(n=1, sd=2)
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

### density plot
d_nsghmc1 <- density(theta[(warmup + 1): iter], bw=0.1)
plot(d_nsghmc1)



## Naive Stochastic Gradient HMC without MH correction

### sample from the target distribution U(\theta) = -2\theta^2 + \theta^4

epsilon <- 1e-1
leapfrog <- 50

iter <- 80000 # total samples
warmup <- 10000 # warmup samples

set.seed(225)

theta <- rep(0, iter) # parameter of interest
r <- rep(0, iter) # momentum

theta[1] <- 0

for (i in 1:(iter - 1)){
  
    
  r[i] <- rnorm(n=1)
  
  theta_ <- theta[i]
  r_ <- r[i]
  
  dev <- - 4 * theta_ * (1 - theta_^2) + rnorm(n=1, sd=2)
  r_ <- r_ - epsilon / 2 * dev
  
  for (m in 1: leapfrog){
    theta_ <- theta_ + epsilon * r_
    dev <- - 4 * theta_ * (1 - theta_^2) + rnorm(n=1, sd=2)
    r_ <- r_ - epsilon * dev
  }
  theta[i+1] <- theta_
}

### density plot
d_nsghmc2 <- density(theta[(warmup + 1): iter], bw=0.1)
plot(d_nsghmc2)


## Stochastic Gradient HMC with Friction

### sample from the target distribution U(\theta) = -2\theta^2 + \theta^4

epsilon <- 1e-1
leapfrog <- 50
alpha <- 0.03
V <- 1

iter <- 80000 # total samples
warmup <- 10000 # warmup samples

set.seed(225)

iter <- 80000 # total samples
warmup <- 10000 # warmup samples

theta <- rep(0, iter) # parameter of interest
r <- rep(0, iter) # momentum
beta <- V * epsilon^2 * 0.5

theta[1] <- 0

for (i in 1:(iter - 1)){
  
  r[i] <- rnorm(n=1)
  
  theta_ <- theta[i]
  r_ <- r[i]
  
  dev <- - 4 * theta_ * (1 - theta_^2) + rnorm(n=1, sd=2)
  r_ <- r_ * epsilon
  sig <- 2 * epsilon^2 * (alpha - beta)
  
  for (m in 1: leapfrog){
    dev <- - 4 * theta_ * (1 - theta_^2) + rnorm(n=1, sd=2)
    r_ <- r_ * (1 - alpha) - dev * epsilon^2 + rnorm(n=1, mean=0, sd=sqrt(2 * sig))
    theta_ <- theta_ + r_
  }
  theta[i+1] <- theta_
}

### density plot
d_sghmc <- density(theta[(warmup + 1): iter], bw=0.1)
plot(d_sghmc)

