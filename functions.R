require(reshape2)

create_data <- function(nspecies=20, nsite=100, nyear=10) {
  # one focal species - probability of detection=0.5
  
  lambda <- matrix(0, nrow = nsite, ncol = nyear)
  
  y <- matrix(0, nrow = nsite, ncol = nyear)
  
  a <- rnorm(nsite, 0, 3)
  b <- rnorm(nyear, 0, 3)
  #b <- sort(b, decreasing = TRUE) # simulating a decline
  
  for (j in 1: nsite) {
    for (t in 1:nyear) {
      lambda[j,t] <- exp(a[j] + b[t])
      y[j,t] <- rpois(1, lambda[j,t])
    }
  }
  y <- as.data.frame(y)
  
  # simulate observed abundance
  
  simfocal <- matrix(0, nrow = nsite, ncol = nyear)
  simfocal <- as.data.frame(simfocal)
  
  x <- rep(0,nyear)
  
  for (t in 1:nyear) {
    x[t] <- paste0("year", t)
  }
  
  colnames(simfocal) <- x
  colnames(y) <- x
  
  x <- rep(0, nsite)
  
  for (j in 1:nsite) {
    x[j] <- paste0("site", j)
  }
  
  simfocal$site <- x
  y$site <- x
  
  z <- ifelse(y==0, 0, 1)
  
  p <- matrix(0, nrow = nsite, ncol = nyear)
  
  for (j in 1: nsite) {
    for (t in 1:nyear) {
      # fixed probability of detection
      p[j,t] <- 0.5
      simfocal[j,t] <- rbinom(1, y[j,t], p[j,t]*z[j,t])
    }
  }
  
  simfocal <- melt(simfocal, id.vars = "site")
  simfocal$species <- "species1"
  colnames(simfocal) <- c("site", "year", "observed", "species")
  y <- melt(y, id.vars = "site")
  colnames(y) <- c("site", "year", "actual")
   
  simdata <- merge(simfocal, y)
  
  # non-focal species
  
  for (i in 1:nspecies) {
    
    lambda <- matrix(0, nrow = nsite, ncol = nyear)
    
    y <- matrix(0, nrow = nsite, ncol = nyear)
    
    a <- rnorm(nsite, 0, 3)
    b <- rnorm(nyear, 0, 3)
    
    for (j in 1: nsite) {
      for (t in 1:nyear) {
        lambda[j,t] <- exp(a[j] + b[t])
        y[j,t] <- rpois(1, lambda[j,t])
      }
    }
    y <- as.data.frame(y)
    
    # simulate observed abundance
    
    sims2 <- matrix(0, nrow = nsite, ncol = nyear)
    sims2 <- as.data.frame(sims2)
    
    x <- rep(0,nyear)
    
    for (t in 1:nyear) {
      x[t] <- paste0("year", t)
    }
    
    colnames(sims2) <- x
    colnames(y) <- x
    
    x <- rep(0, nsite)
    
    for (j in 1:nsite) {
      x[j] <- paste0("site", j)
    }
    
    sims2$site <- x
    y$site <- x
    
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
    sims2$species <- paste0("species", (i+1))
    
    colnames(sims2) <- c("site", "year", "observed", "species")
    y <- melt(y, id.vars = "site")
    colnames(y) <- c("site", "year", "actual")
    
    sims2 <- merge(sims2, y)

    simdata <- rbind(simdata, sims2)
  }

  simdata$site <- factor(simdata$site, levels = paste0("site", 1:nsite))
  simdata$species <- as.factor(simdata$species)
  
  return(simdata)
  
}

###############################

run_model <- function(data, species_list=unique(data$species), n_chains=4, n_iter=10000) {
  output <- list()
  
  for (i in species_list) {
    simdata1 <- subset(data, species==i)
    simdata1 <- dcast(data = simdata1,formula = site~year,fun.aggregate = sum,value.var = "observed")
    rownames(simdata1) <- simdata1$site
    simdata1 <- simdata1[,-1]
    
    jags_data <- list(nsite = nrow(simdata1), nyear = ncol(simdata1), SI = simdata1)
    
    output[[i]] <- jags(data=jags_data,
                        parameters.to.save=c("y", "b"),
                        model.file="SI_model.bug", 
                        n.chains=n_chains, n.iter=n_iter)
  }
  
  return(output)
}
