library(faraway)
library(R2jags)

#focussing on one species, simulating total abundance across years and sites

sims <- matrix(0, nrow = 100, ncol = 10)
sims <- as.data.frame(sims)

x <- rep(0,10)

for (i in 1:10) {
  x[i] <- paste0("year", i)
}

colnames(sims) <- x

# real abundance

nsite <- 100
nyear <- 10

lambda <- matrix(0, nrow = nsite, ncol = nyear)

y <- matrix(0, nrow = nsite, ncol = nyear)

a <- rnorm(nsite, 0, 4)
b <- rnorm(nyear, 0, 4)

for (j in 1: nsite) {
  for (t in 1:nyear) {
    lambda[j,t] <- exp(a[j] + b[t])
    y[j,t] <- rpois(1, lambda[j,t])
  }
}

# simulate observed abundance

z <- ifelse(y==0, 0, 1)

p <- matrix(0, nrow = nsite, ncol = nyear)

for (j in 1: nsite) {
  for (t in 1:nyear) {
    # beta distribution from Isaac et al.
    p[j,t] <- runif(1,0,1)
    sims[j,t] <- rbinom(1, y[j,t], p[j,t]*z[j,t])
  }
}

data <- list(nyear=nyear, nsite=nsite, SI=sims)

out <- jags(data=data,
            parameters.to.save=c("y"),
            model.file="Site-abundance-model/SI_model.bug", 
            n.chains=4, n.iter=40000)

###################################################################

