## Metropolis-Hastings algorithm

### normal model with known variance
### normal conjugate prior for the mean
### proposal distribution: unif(old_state - radius, old_state + radius)

mu <- 0 # mean of prior normal
tau <- 10 # sd of prior normal
sigma <- 3 # sd of normal model
delta <- 5 # radius of proposal interval (uniform distribution around the current state)

set.seed(225)

Y <- rnorm(100, mean = 5, sd = sigma) # normal model: N(5, 9)

iter <- 2000 # total samples
warmup <- 1000 # warmup samples


theta <- rep(0, iter)

theta[1] <- rnorm(1, mean = mu, sd = tau)

for (i in 2:iter){
  old_state <- theta[i-1]
  new_state <- runif(1, min = old_state - delta, max = old_state + delta) # propose new state
  a <- min(1, exp(sum(log(dnorm(Y, mean = new_state, sd = sigma))) 
                  + log(dnorm(new_state, mu, tau)) 
                  - sum(log(dnorm(Y, mean = old_state, sd = sigma))) 
                  - log(dnorm(old_state, mu, tau))))
  if (runif(n=1) < a){theta[i] <- new_state} else {theta[i] <- old_state}
}

# trace plot
plot(c(1:2000), theta[1: 2000], type = "l", xlab = "index", ylab = "theta")

# histogram
hist(theta)

# sample average of posterior mean
# theoretical value: 5.566
print(mean(theta[1001: 2000])) 

## Multivariate Metropolis-Hastings algorithm

### normal model with known covariance matrix
### normal conjugate prior for the two-dimensional mean
library(MASS)
library(mvtnorm)
library(hexbin)
library(latex2exp)

mu <- rep(0, 2)
Tau <- 10 * diag(2)
Sigma <- 3 * diag(2)
delta <- 1
n <- 100

set.seed(225)

Y <- mvrnorm(n=n, c(5, 3), Sigma = Sigma) # normal model: N(c(5, 3), 9I_2)

iter <- 2000 # total samples
warmup <- 1000 # warmup samples

theta <- matrix(rep(0, iter * 2), ncol = 2)

theta[1, ] <- mvrnorm(1, mu, Tau)

for (i in 2:iter){
  old_state <- theta[i-1,]
  new_state <- mvrnorm(1, mu = old_state, Sigma = delta * diag(2)) # propose new state
  a <- min(1, dmvnorm(new_state, mu, Tau) / dmvnorm(old_state, mu, Tau) * 
             (exp(sum(log(dmvnorm(Y, new_state, Sigma))) - sum(log(dmvnorm(Y, old_state, Sigma))))))
  if (runif(n=1) < a){theta[i,] <- new_state} else {theta[i,] <- old_state}
}

### scatter plot
plot(theta[, 1], theta[, 2], type = "p", 
     xlab = TeX("$\\theta_1$"),
     ylab = TeX("$\\theta_2$"), 
     main = "Metropolis-Hastings") # very likely reject the new state
### posterior mean
print(colMeans(theta)) # theoretical value: 5.0376, 3.3197

## Hamiltonian Monte Carlo

### normal model with known covariance matrix
### normal conjugate prior for the two-dimensional mean

mu <- c(0, 0)
Tau <- 10 * diag(2)
Sigma <- 3 * diag(2)
epsilon <- 1e-2
n <- 100
leapfrog <- 100

set.seed(225)

Y <- mvrnorm(n = n, c(5, 3), Sigma = Sigma) # normal model: N(c(5, 3), 9I_2)

iter <- 2000 # total samples
warmup <- 1000 # warmup samples

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
  dev <- - solve(Sigma) %*% (colSums(Y) - n * theta_) + solve(Tau) %*% (theta_ - mu)
  r_ <- r_ - epsilon / 2 * dev
  
  H1 <- - sum(log(dmvnorm(Y, mean = theta_, sigma = Sigma))) - 
          log(dmvnorm(t(theta_), mean = mu, sigma = Tau)) + 1/2 * t(r_) %*% r_
  H2 <- - sum(log(dmvnorm(Y, mean = theta[i, ], sigma = Sigma))) -
          log(dmvnorm(theta[i, ], mean = mu, sigma = Tau)) + 1/2 * t(r[i, ]) %*% r[i, ]
  rho <- exp(H1 - H2)
  u <- runif(n=1)
  if (u < min(1, rho)){
    theta[i+1, ] <- theta_
  } else {
    theta[i+1, ] <- theta[i, ]
  }
}

### scatter plot
plot(theta[, 1], theta[, 2], type = "p", 
     xlab = TeX("$\\theta_1$"),
     ylab = TeX("$\\theta_2$"), 
     main = "Hamiltonian MC") 

### posterior mean
print(colMeans(theta)) # theoretical value: 5.0376, 3.3197


## Hamiltonian Monte Carlo without MH correction

### normal model with known covariance matrix
### normal conjugate prior for the two-dimensional mean

mu <- c(0, 0)
Tau <- 10 * diag(2)
Sigma <- 3 * diag(2)
epsilon <- 1e-2
n <- 100
leapfrog <- 100

set.seed(225)

Y <- mvrnorm(n = n, c(5, 3), Sigma = Sigma) # normal model: N(c(5, 3), 9I_2)

iter <- 2000 # total samples
warmup <- 1000 # warmup samples

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
}

### scatter plot
plot(theta[, 1], theta[, 2], type = "p", 
     xlab = TeX("$\\theta_1$"),
     ylab = TeX("$\\theta_2$"), 
     main = "Hamiltonian MC without MH correction") 

### posterior mean
print(colMeans(theta)) # theoretical value: 5.0376, 3.3197


## Naive Stochastic Gradient HMC

### normal model with known covariance matrix
### normal conjugate prior for the two-dimensional mean

mu <- c(0, 0)
Tau <- 10 * diag(2)
Sigma <- 3 * diag(2)
epsilon <- 5e-2
batch_size <- 5
n <- 100
leapfrog <- 100

set.seed(225)

Y <- mvrnorm(n = n, c(5, 3), Sigma = Sigma) # normal model: N(c(5, 3), 9I_2)

iter <- 2000 # total samples
warmup <- 1000 # warmup samples

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
  batch <- Y[sample(c(1: n), batch_size), ]
  dev <- - n / batch_size * solve(Sigma) %*% (colSums(batch) - batch_size * theta_) + solve(Tau) %*% (theta_ - mu)
  r_ <- r_ - epsilon / 2 * dev
  
  H1 <- - sum(log(dmvnorm(Y, mean = theta_, sigma = Sigma))) - 
          log(dmvnorm(t(theta_), mean = mu, sigma = Tau)) + 1/2 * t(r_) %*% r_
  H2 <- - sum(log(dmvnorm(Y, mean = theta[i, ], sigma = Sigma))) -
          log(dmvnorm(theta[i, ], mean = mu, sigma = Tau)) + 1/2 * t(r[i, ]) %*% r[i, ]
  rho <- exp(H1 - H2)
  u <- runif(n=1)
  if (u < min(1, rho)){
    theta[i+1, ] <- theta_
  } else {
    theta[i+1, ] <- theta[i, ]
  }
}

### scatter plot
plot(theta[, 1], theta[, 2], type = "p", 
     xlab = TeX("$\\theta_1$"),
     ylab = TeX("$\\theta_2$"), 
     main = "Naive Stochastic Gradient HMC") 

### posterior mean
print(colMeans(theta)) # theoretical value: 5.0376, 3.3197



## Naive Stochastic Gradient HMC without MH correction

### normal model with known covariance matrix
### normal conjugate prior for the two-dimensional mean

mu <- c(0, 0)
Tau <- 10 * diag(2)
Sigma <- 3 * diag(2)
epsilon <- 5e-2
batch_size <- 5
n <- 100
leapfrog <- 100

set.seed(225)

Y <- mvrnorm(n = n, c(5, 3), Sigma = Sigma) # normal model: N(c(5, 3), 9I_2)

iter <- 2000 # total samples
warmup <- 1000 # warmup samples

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
}

### scatter plot
plot(theta[, 1], theta[, 2], type = "p", 
     xlab = TeX("$\\theta_1$"),
     ylab = TeX("$\\theta_2$"), 
     main = "Naive Stochastic Gradient HMC without MH correction") 

### posterior mean
print(colMeans(theta)) # theoretical value: 5.0376, 3.3197


## Stochastic Gradient HMC with Friction

### normal model with known covariance matrix

mu <- c(0, 0)
Tau <- 10 * diag(2)
Sigma <- 3 * diag(2)
epsilon <- 5e-2
batch_size <- 5
n <- 100
leapfrog <- 100
alpha <- 0.03
V <- 1

set.seed(225)

Y <- mvrnorm(n = n, c(5, 3), Sigma = Sigma) # normal model: N(c(5, 3), 9I_2)

iter <- 2000 # total samples
warmup <- 1000 # warmup samples

theta <- matrix(rep(0, iter * 2), ncol = 2) # parameter of interest
r <- matrix(rep(0, iter * 2), ncol = 2) # momentum
beta <- V * epsilon^2 * 0.5

theta[1, ] <- c(0,0)

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

### scatter plot
plot(theta[, 1], theta[, 2], type = "p", 
     xlab = TeX("$\\theta_1$"),
     ylab = TeX("$\\theta_2$"), 
     main = "Naive Stochastic Gradient HMC without MH correction") 

### posterior mean
print(colMeans(theta)) # theoretical value: 5.0376, 3.3197

