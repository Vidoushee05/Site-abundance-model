require(reshape2)
require(R2jags)

create_data <- function(nspecies=20, nsite=100, nyear=10, decline=FALSE) {
  # one focal species - probability of detection=0.5
  
  lambda <- matrix(0, nrow = nsite, ncol = nyear)
  
  y <- matrix(0, nrow = nsite, ncol = nyear)
  
  a <- rnorm(nsite, 0, 3)
  if (decline==FALSE) {b <- rnorm(nyear, 0, 3)}
  
  if (decline==TRUE) {
    # simulating a 30% decline in occupancy over 10 years
    b <- rep(0,nyear)
    b[1] <- rnorm(1, 0, 3)
    b[nyear] <- b[1] + log(0.7)
    b[2:(nyear-1)] <- runif(nyear-2, b[nyear], b[1])
    b[2:(nyear-1)] <- sort(b[2:(nyear-1)], decreasing = TRUE)
  }
  

  for (j in 1: nsite) {
    for (t in 1:nyear) {
      lambda[j,t] <- exp(a[j] + b[t])
      y[j,t] <- rpois(1, lambda[j,t])
    }
  }
  y <- as.data.frame(y)
  
  # store values of b
  trend <- data.frame(species="species1", year="year1", b_t=b[1], stringsAsFactors = FALSE)
  for (t in 2:nyear) {
    trend[t,] <- c("species1", paste0("year", t), b[t])
  }
  
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
  
  simfocal <- reshape2::melt(simfocal, id.vars = "site")
  simfocal$species <- "species1"
  colnames(simfocal) <- c("site", "year", "observed", "species")
  y <- reshape2::melt(y, id.vars = "site")
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
    
    trend2 <- data.frame(species=paste("species", i+1), year="year1", b_t=b[1], stringsAsFactors = FALSE)
    for (t in 2:nyear) {
      trend2[t,] <- c(paste0("species", i+1), paste0("year", t), b[t])
    }
    trend <- rbind(trend, trend2)
    
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
    
    sims2 <- reshape2::melt(sims2, id.vars = "site")
    sims2$species <- paste0("species", (i+1))
    
    colnames(sims2) <- c("site", "year", "observed", "species")
    y <- reshape2::melt(y, id.vars = "site")
    colnames(y) <- c("site", "year", "actual")
    
    sims2 <- merge(sims2, y)

    simdata <- rbind(simdata, sims2)
  }
  
  trend$b_t <- round(as.numeric(trend$b_t),2)
  simdata$site <- factor(simdata$site, levels = paste0("site", 1:nsite))
  simdata$species <- as.factor(simdata$species)
  
  simulations <- list(simdata, trend)
  
  return(simulations)
  
}

###############################

run_model <- function(data, species_list=unique(data$species), n_chains=4, n_iter=10000) {
  output <- list()
  
  for (i in species_list) {
    simdata1 <- subset(data, species==i)
    simdata1 <- reshape2::dcast(data = simdata1,formula = site~year,fun.aggregate = sum,value.var = "observed")
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

################################################

create_data.nb <- function(nspecies=20, nsite=100, nyear=10, decline=FALSE) {
  # one focal species - probability of detection=0.5
  
  lambda <- matrix(0, nrow = nsite, ncol = nyear)
  
  y <- matrix(0, nrow = nsite, ncol = nyear)
  
  a <- rexp(nsite, 1.5)
  if (decline==FALSE) {b <- rexp(nyear, 1.5)}
  
  if (decline==TRUE) {
    # simulating a 30% decline in occupancy over 10 years
    b <- rep(0,nyear)
    b[1] <- rexp(1, 1.5)
    b[nyear] <- b[1] + log(0.7)
    b[2:(nyear-1)] <- runif(nyear-2, b[nyear], b[1])
    b[2:(nyear-1)] <- sort(b[2:(nyear-1)], decreasing = TRUE)
  }
  
  
  for (j in 1: nsite) {
    for (t in 1:nyear) {
      lambda[j,t] <- exp(a[j] + b[t])
      y[j,t] <- rnbinom(1, lambda[j,t], 0.29)
    }
  }
  y <- as.data.frame(y)
  
  # store values of b
  trend <- data.frame(species="species1", year="year1", b_t=b[1], stringsAsFactors = FALSE)
  for (t in 2:nyear) {
    trend[t,] <- c("species1", paste0("year", t), b[t])
  }
  
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
  
  simfocal <- reshape2::melt(simfocal, id.vars = "site")
  simfocal$species <- "species1"
  colnames(simfocal) <- c("site", "year", "observed", "species")
  y <- reshape2::melt(y, id.vars = "site")
  colnames(y) <- c("site", "year", "actual")
  
  simdata <- merge(simfocal, y)
  
  # non-focal species
  
  for (i in 1:nspecies) {
    
    lambda <- matrix(0, nrow = nsite, ncol = nyear)
    
    y <- matrix(0, nrow = nsite, ncol = nyear)
    
    a <- rexp(nsite, 1.5)
    b <- rexp(nyear, 1.5)
    
    for (j in 1: nsite) {
      for (t in 1:nyear) {
        lambda[j,t] <- exp(a[j] + b[t])
        y[j,t] <- rnbinom(1, lambda[j,t], 0.29)
      }
    }
    y <- as.data.frame(y)
    
    trend2 <- data.frame(species=paste("species", i+1), year="year1", b_t=b[1], stringsAsFactors = FALSE)
    for (t in 2:nyear) {
      trend2[t,] <- c(paste0("species", i+1), paste0("year", t), b[t])
    }
    trend <- rbind(trend, trend2)
    
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
    
    sims2 <- reshape2::melt(sims2, id.vars = "site")
    sims2$species <- paste0("species", (i+1))
    
    colnames(sims2) <- c("site", "year", "observed", "species")
    y <- reshape2::melt(y, id.vars = "site")
    colnames(y) <- c("site", "year", "actual")
    
    sims2 <- merge(sims2, y)
    
    simdata <- rbind(simdata, sims2)
  }
  
  trend$b_t <- round(as.numeric(trend$b_t),2)
  simdata$site <- factor(simdata$site, levels = paste0("site", 1:nsite))
  simdata$species <- as.factor(simdata$species)
  
  simulations <- list(simdata, trend)
  
  return(simulations)
  
}

###############################################

