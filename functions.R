require(reshape2)

create_data <- function(nspecies=20, nsite=100, nyear=10) {
  # one focal species - probability of detection=0.5
  
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
  
  simfocal <- matrix(0, nrow = nsite, ncol = nyear)
  simfocal <- as.data.frame(simfocal)
  
  x <- rep(0,nyear)
  
  for (i in 1:nyear) {
    x[i] <- paste0("year", i)
  }
  
  colnames(simfocal) <- x
  
  x <- rep(0, nsite)
  
  for (i in 1:nsite) {
    x[i] <- paste0("site", i)
  }
  
  simfocal$site <- x
  
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
  
  simdata <- simfocal
  
  # non-focal species
  
  for (i in 1:nspecies) {
    
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
    
    simdata <- rbind(simdata, sims2)
  }
  colnames(simdata) <- c("site", "year", "abundance", "species")
  
  return(simdata)
}
