require(faraway)
require(R2jags)
require(reshape2)
require(ggplot2)

# focussing on one species
# simulating total abundance across years and sites

# real abundance

nsite <- 50
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

sim <- create_data1(15,50,10, decline = TRUE)
simdata <- sim[[1]]

#data <- list(nsite = nrow(simdata1), nyear = ncol(simdata1), SI = simdata1)

# out <- jags(data=data,
#             parameters.to.save=c("y", "b"),
#             model.file="SI_model.bug", 
#             n.chains=4, n.iter=50000)

simdata1 <- subset(simdata, species==i)
simdata1 <- reshape2::dcast(data = simdata1,formula = site~year,fun.aggregate = sum,value.var = "observed")
rownames(simdata1) <- simdata1$site
simdata1 <- simdata1[,-1]

jags_data <- list(nsite = nrow(simdata1), nyear = ncol(simdata1), SI = simdata1)

output <- jags(data=jags_data,
               parameters.to.save=c("b", "r"),
               model.file="nb_model.bug", 
               n.chains=4, n.iter=10000)

out <- run_model(simdata, species_list = c("species1"), n_chains=3, n_iter=5000)


########################################

# plotting model output

# trend in actual abundance

trend <- sim[[2]]

actual <- subset(trend, species == "species1")$b_t

est <- out[["species1"]]$BUGSoutput$summary[paste0("b[", 1:nyear, "]"), c(1,3,7,8)]
est <- as.data.frame(est)
est$year <- 1:nyear
est$actual <- actual

p <- ggplot(est, aes(year, mean)) + geom_line(aes(color="Estimated")) +
  geom_ribbon(data=est, aes(ymin=`2.5%`, ymax=`97.5%`), alpha=0.3) + 
  geom_point(data=est, aes(year, actual, color="Actual")) + ylab("b_t") +
  scale_x_continuous(breaks=seq(1, 10, by=1)) +
  theme(legend.title = element_blank()) 


###################################################

# with visits

# for a given site

# given year

j=1
t=1
nspecies=20
mv=10

simdata <- data.frame()
a <- matrix(nrow = nsite, ncol=nspecies+1)
b <- matrix(nrow = nyear, ncol=nspecies+1)
lambda <- array(dim=c(nsite, nyear, nspecies+1))
y <- array(dim=c(nsite, nyear, nspecies+1))

d <- rep(0,nyear)
d[1] <- rnorm(1, 0, 3)
d[nyear] <- d[1] + log(0.7)
d[2:(nyear-1)] <- runif(nyear-2, d[nyear], d[1])
d[2:(nyear-1)] <- sort(d[2:(nyear-1)], decreasing = TRUE)

for (j in 1:nsite) {
  for (t in 1:nyear) {
    a[j,] <- rnorm(nspecies+1, 0, 3)
    b[t,] <- rnorm(nspecies+1, 0, 3)
    if (decline) {b[,1] <- d}
    lambda[j,t,] <- exp(a[j,]+b[t,])
    f1 <- function(x) rpois(1,x)
    y[j,t,] <- sapply(lambda[j,t,], f1)
    richness <- length(which(y[j,t,]!=0))/(nspecies+1)
    
    nvisits <- rbinom(1, mv, richness)
    
    f2 <- function(x, p_detect) rbinom(nvisits, x, p_detect)
    p_detect <- c(0.5, rbeta(20,2,2))
    sim <- sapply(y[j,t,], f2, p_detect=p_detect)
    if (nvisits==0) {
      sim <- data.frame(NULL)
    }
    else{
      sim <- reshape2::melt(sim)
      if (nvisits==1) {
        sim <- data.frame(rep(1,nspecies+1),1:(nspecies+1),sim[,1])
      }
      colnames(sim) <- c("visit", "species", "obs")
      f3 <- function(x) rep(x, nvisits)
      sim$actual <- sapply(y[j,t,], f3)[1:length(sim$visit)]
      sim$site <- j
      sim$year <- t
    }
    simdata <- rbind(simdata, sim)
  }
}

for (i in 1:10) {x[i] <- mean(subset(simdata, species==1 & year==i)$obs)}
for (i in 1:10) {x[i] <- mean(subset(simdata, species=="species1" & year==paste0("year", i))$actual)}

data <- create_data5(decline=FALSE)
simdata <- data[[1]]

out <- run_model(simdata, species_list = c(1), n_chains=3, n_iter=5000)

trend <- data[[2]]

actual <- trend[,1]

est <- out[[1]][[1]]$BUGSoutput$summary[paste0("b[", 1:nyear, "]"), c(1,3,7,8)]
est <- as.data.frame(est)
est$year <- 1:nyear
est$actual <- actual

p <- ggplot(est, aes(year, mean)) + geom_line(aes(color="Estimated")) +
  geom_ribbon(data=est, aes(ymin=`2.5%`, ymax=`97.5%`), alpha=0.3) + 
  geom_point(data=est, aes(year, actual, color="Actual")) + ylab("b_t") +
  scale_x_continuous(breaks=seq(1, 10, by=1)) +
  theme(legend.title = element_blank()) 
