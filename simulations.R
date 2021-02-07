require(faraway)
require(R2jags)
require(reshape2)

# focussing on one species
# simulating total abundance across years and sites

# real abundance

nsite <- 100
nyear <- 10

lambda <- matrix(0, nrow = nsite, ncol = nyear)

y <- matrix(0, nrow = nsite, ncol = nyear)

a <- rnorm(nsite, 0, 4)
b <- rnorm(nyear, 0, 4)
#b <- sort(b, decreasing = TRUE) # simulating a decline

for (j in 1: nsite) {
  for (t in 1:nyear) {
    lambda[j,t] <- exp(a[j] + b[t])
    y[j,t] <- rpois(1, lambda[j,t])
  }
}

# simulate observed abundance

sims <- matrix(0, nrow = nsite, ncol = nyear)
sims <- as.data.frame(sims)

x <- rep(0,nyear)

for (i in 1:nyear) {
  x[i] <- paste0("year", i)
}

colnames(sims) <- x

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
            parameters.to.save=c("y", "b"),
            model.file="Site-abundance-model/SI_model.bug", 
            n.chains=4, n.iter=40000)

###################################################################

# simulating a dataset with more species

# with one focal species - probability of detection=0.5

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

sims1 <- matrix(0, nrow = nsite, ncol = nyear)
sims1 <- as.data.frame(sims1)

x <- rep(0,nyear)

for (i in 1:nyear) {
  x[i] <- paste0("year", i)
}

colnames(sims1) <- x

x <- rep(0, nsite)

for (i in 1:nsite) {
  x[i] <- paste0("site", i)
}

sims1$site <- x

z <- ifelse(y==0, 0, 1)

p <- matrix(0, nrow = nsite, ncol = nyear)

for (j in 1: nsite) {
  for (t in 1:nyear) {
    # focal species
    p[j,t] <- 0.5
    sims1[j,t] <- rbinom(1, y[j,t], p[j,t]*z[j,t])
  }
}

sims1 <- melt(sims1, id.vars = "site")
sims1$species <- "species1"

# non-focal species

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

sims2 <- matrix(0, nrow = nsite, ncol = nyear)
sims2 <- as.data.frame(sims2)

x <- rep(0,nyear)

for (i in 1:nyear) {
  x[i] <- paste0("year", i)
}

colnames(sims2) <- x

x <- rep(0, nsite)

for (i in 1:nsite) {
  x[i] <- paste0("site", i)
}

sims2$site <- x

z <- ifelse(y==0, 0, 1)

p <- matrix(0, nrow = nsite, ncol = nyear)

for (j in 1: nsite) {
  for (t in 1:nyear) {
    # beta distribution from Isaac et al.
    p[j,t] <- rbeta(1, 2, 2)
    sims2[j,t] <- rbinom(1, y[j,t], p[j,t]*z[j,t])
  }
}

sims2 <- melt(sims2, id.vars = "site")
sims2$species <- "species2"

sim <- rbind(sims1, sims2)

######################################

# formatting the data to run JAGS

# for each species

source("functions.R")

simdata <- create_data(15,50,10)

simdata1 <- subset(simdata, species==species_name)
simdata1 <- dcast(data = simdata1,formula = site~year,fun.aggregate = sum,value.var = "observed")
rownames(simdata1) <- simdata1$site
simdata1 <- simdata1[,-1]

#data <- list(nsite = nrow(simdata1), nyear = ncol(simdata1), SI = simdata1)

# out <- jags(data=data,
#             parameters.to.save=c("y", "b"),
#             model.file="SI_model.bug", 
#             n.chains=4, n.iter=50000)

out <- run_model(simdata, species_list = c("species1"), n_iter=50000)
