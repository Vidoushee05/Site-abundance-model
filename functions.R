require(reshape2)
require(R2jags)
require(zoo)
require(dplyr)

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

run_model <- function(simdata, species_list=unique(data$species), n_chains=4, n_iter=10000) {
  output <- list()
  missing <- list()
  
  for (i in species_list) {
    simdata1 <- subset(simdata, species==i)
    simdata1 <- reshape2::dcast(data = simdata1,formula = site~year,
                                fun.aggregate = function(x) round(mean(x)),
                                value.var = "obs")
    #rownames(simdata1) <- simdata1$site
    simdata1 <- simdata1[,-1]
    missing[[i]] <- sum(is.na(simdata1))/prod(dim(simdata1))
    simdata1 <- t(apply(simdata1, 1, na.locf, na.rm=FALSE))
    simdata1 <- t(apply(simdata1, 1, na.locf, fromLast=T))
    
    jags_data <- list(nsite = nrow(simdata1), nyear = ncol(simdata1), SI = simdata1)
    
    output[[i]] <- jags(data=jags_data,
                        parameters.to.save=c("a", "b", "c"),
                        model.file="SI_model.bug", 
                        n.chains=n_chains, n.iter=n_iter,
                        progress.bar = "none")
  }
  
  return(list(output, missing))
}

###############################################

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

# with visits
create_data2 <- function(nspecies=20, nsite=50, nyear=10, mv=10, decline=FALSE, nb=FALSE){
  simdata <- data.frame()
  a <- matrix(nrow = nsite, ncol=nspecies+1)
  b <- matrix(nrow = nyear, ncol=nspecies+1)
  p_detect <- rep(0, nspecies+1)
  lambda <- array(dim=c(nsite, nyear, nspecies+1))
  y <- array(dim=c(nsite, nyear, nspecies+1))

  if (decline) {
    # simulating a 30% decline in focal occupancy over 10 years
    d <- rep(0,nyear)
    d[1] <- rnorm(1, 0, 3)
    d[nyear] <- d[1] + log(0.7)
    d[2:(nyear-1)] <- runif(nyear-2, d[nyear], d[1])
    d[2:(nyear-1)] <- sort(d[2:(nyear-1)], decreasing = TRUE)
  }
  
  for (i in 1:(nspecies+1)) {
    a[,i] <- rnorm(nsite, 0, 3)
    b[,i] <- rnorm(nyear, 0, 3)
    p_detect[i] <- rbeta(1, 2, 2)
  }
  
  p_detect[1] <- 0.5  #fixed prob of detection for focal
  
  if (decline) {b[,1] <- d}
  
  for (j in 1:nsite) {
    for (t in 1:nyear) {
      lambda[j,t,] <- exp(a[j,]+b[t,])
      f1 <- function(x) rpois(1,x)
      y[j,t,] <- sapply(lambda[j,t,], f1)
      if (nb) {
        f11 <- function(x, r=0.425) rnbinom(1,x,r)
        y[j,t,] <- sapply(lambda[j,t,], f11)
      }
      richness <- length(which(y[j,t,]!=0))/(nspecies+1)
      
      nvisits <- rbinom(1, mv, richness)
      
      f2 <- function(x, p_detect) rbinom(nvisits, x, p_detect)
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
  out <- list()
  out[["data"]] <- simdata
  out[["params"]] <- b
  
  return(out)
}

##############################################

# simulating increase in prob of focal detection
create_data3 <- function(nspecies=20, nsite=50, nyear=10, mv=10, decline=FALSE, nb=FALSE){
  simdata <- data.frame()
  a <- matrix(nrow = nsite, ncol=nspecies+1)
  b <- matrix(nrow = nyear, ncol=nspecies+1)
  p_detect <- matrix(nrow = nyear, ncol=nspecies+1)
  lambda <- array(dim=c(nsite, nyear, nspecies+1))
  y <- array(dim=c(nsite, nyear, nspecies+1))
  
  if (decline) {
    # simulating a 30% decline in focal occupancy over 10 years
    d <- rep(0,nyear)
    d[1] <- rnorm(1, 0, 3)
    d[nyear] <- d[1] + log(0.7)
    d[2:(nyear-1)] <- runif(nyear-2, d[nyear], d[1])
    d[2:(nyear-1)] <- sort(d[2:(nyear-1)], decreasing = TRUE)
  }
  
  for (i in 1:(nspecies+1)) {
    a[,i] <- rnorm(nsite, 0, 3)
    b[,i] <- rnorm(nyear, 0, 3)
    p_detect[,i] <- rbeta(1, 2, 2)
  }
  
  for(t in 1:nyear) {
    p_detect[t,1] <- 0.4+(t-1)/((nyear-1)/0.1)
  }
  
  if (decline) {b[,1] <- d}
  
  for (j in 1:nsite) {
    for (t in 1:nyear) {
      lambda[j,t,] <- exp(a[j,]+b[t,])
      f1 <- function(x) rpois(1,x)
      y[j,t,] <- sapply(lambda[j,t,], f1)
      if (nb) {
        f11 <- function(x, r=0.425) rnbinom(1,x,r)
        y[j,t,] <- sapply(lambda[j,t,], f11)
      }
      richness <- length(which(y[j,t,]!=0))/(nspecies+1)
      
      nvisits <- rbinom(1, mv, richness)
      
      f2 <- function(x, p_detect) rbinom(nvisits, x, p_detect)
      sim <- sapply(y[j,t,], f2, p_detect=p_detect[t,])
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
  out <- list()
  out[["data"]] <- simdata
  out[["params"]] <- b
  
  return(out)
}

##############################################

# reduced recording effort
create_data4 <- function(nspecies=20, nsite=50, nyear=10, mv=10, decline=FALSE, nb=FALSE){
  simdata <- data.frame()
  a <- matrix(nrow = nsite, ncol=nspecies+1)
  b <- matrix(nrow = nyear, ncol=nspecies+1)
  p_detect <- rep(0, nspecies+1)
  lambda <- array(dim=c(nsite, nyear, nspecies+1))
  y <- array(dim=c(nsite, nyear, nspecies+1))
  
  if (decline) {
    # simulating a 30% decline in focal occupancy over 10 years
    d <- rep(0,nyear)
    d[1] <- rnorm(1, 0, 3)
    d[nyear] <- d[1] + log(0.7)
    d[2:(nyear-1)] <- runif(nyear-2, d[nyear], d[1])
    d[2:(nyear-1)] <- sort(d[2:(nyear-1)], decreasing = TRUE)
  }
  
  for (i in 1:(nspecies+1)) {
    a[,i] <- rnorm(nsite, 0, 3)
    b[,i] <- rnorm(nyear, 0, 3)
    p_detect[i] <- rbeta(1, 2, 2)
  }
  
  p_detect[1] <- 0.5  #fixed prob of detection for focal
  
  if (decline) {b[,1] <- d}
  
  for (j in 1:nsite) {
    for (t in 1:nyear) {
      lambda[j,t,] <- exp(a[j,]+b[t,])
      f1 <- function(x) rpois(1,x)
      y[j,t,] <- sapply(lambda[j,t,], f1)
      if (nb) {
        f11 <- function(x, r=0.425) rnbinom(1,x,r)
        y[j,t,] <- sapply(lambda[j,t,], f11)
      }
      richness <- length(which(y[j,t,]!=0))/(nspecies+1)
      
      nvisits <- rbinom(1, mv, richness)
      
      f2 <- function(x, p_detect) rbinom(nvisits, x, p_detect)
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
  
  #increase short lists from 10% to 30%
  for (t in 1:nyear) {
    single <- 0.1+(t-1)/((nyear-1)/0.2)
    sub <- subset(simdata, year==t)
    vis <- unique(sub[,c("visit", "site")])
    sel <- vis[sample(nrow(vis), single*nrow(vis)),]
    for (i in 1:nrow(sel)) {
      sub <- match_df(sub, sel[i,])
      sp <- sample(1:(nspecies+1), 1)
      sub <- subset(sub, species!=sp)
      simdata <- simdata[-as.integer(rownames(sub)),]
      sub <- subset(simdata, year==t)
    }
  }
  
  out <- list()
  out[["data"]] <- simdata
  out[["params"]] <- b
  
  return(out)
}

##############################################

# visits double over 10 years
create_data5 <- function(nspecies=20, nsite=50, nyear=10, mv=10, decline=FALSE, nb=FALSE){
  simdata <- data.frame()
  a <- matrix(nrow = nsite, ncol=nspecies+1)
  b <- matrix(nrow = nyear, ncol=nspecies+1)
  p_detect <- rep(0, nspecies+1)
  lambda <- array(dim=c(nsite, nyear, nspecies+1))
  y <- array(dim=c(nsite, nyear, nspecies+1))
  
  if (decline) {
    # simulating a 30% decline in focal occupancy over 10 years
    d <- rep(0,nyear)
    d[1] <- rnorm(1, 0, 3)
    d[nyear] <- d[1] + log(0.7)
    d[2:(nyear-1)] <- runif(nyear-2, d[nyear], d[1])
    d[2:(nyear-1)] <- sort(d[2:(nyear-1)], decreasing = TRUE)
  }
  
  for (i in 1:(nspecies+1)) {
    a[,i] <- rnorm(nsite, 0, 3)
    b[,i] <- rnorm(nyear, 0, 3)
    p_detect[i] <- rbeta(1, 2, 2)
  }
  
  p_detect[1] <- 0.5  #fixed prob of detection for focal
  
  if (decline) {b[,1] <- d}
  
  for (j in 1:nsite) {
    for (t in 1:nyear) {
      lambda[j,t,] <- exp(a[j,]+b[t,])
      f1 <- function(x) rpois(1,x)
      y[j,t,] <- sapply(lambda[j,t,], f1)
      if (nb) {
        f11 <- function(x, r=0.425) rnbinom(1,x,r)
        y[j,t,] <- sapply(lambda[j,t,], f11)
      }
      richness <- length(which(y[j,t,]!=0))/(nspecies+1)
      
      nvisits <- rbinom(1, mv, richness)
      
      f2 <- function(x, p_detect) rbinom(nvisits, x, p_detect)
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
  
  for (t in 1:(nyear-1)) {
    prop <- 0.5-(t-1)/((nyear-1)/0.5)
    sub <- subset(simdata, year==t)
    vis <- unique(sub[,c("visit", "site")])
    sel <- vis[sample(nrow(vis), prop*nrow(vis)),]
    for (i in 1:nrow(sel)) {
      sub <- match_df(sub, sel[i,])
      simdata <- simdata[-as.integer(rownames(sub)),]
      sub <- subset(simdata, year==t)
    }
  }
  
  out <- list()
  out[["data"]] <- simdata
  out[["params"]] <- b
  
  return(out)
}

##############################################

# increased visits biased towards focal
create_data6 <- function(nspecies=20, nsite=50, nyear=10, mv=10, decline=FALSE, nb=FALSE){
  simdata <- data.frame()
  a <- matrix(nrow = nsite, ncol=nspecies+1)
  b <- matrix(nrow = nyear, ncol=nspecies+1)
  p_detect <- rep(0, nspecies+1)
  lambda <- array(dim=c(nsite, nyear, nspecies+1))
  y <- array(dim=c(nsite, nyear, nspecies+1))
  
  if (decline) {
    # simulating a 30% decline in focal occupancy over 10 years
    d <- rep(0,nyear)
    d[1] <- rnorm(1, 0, 3)
    d[nyear] <- d[1] + log(0.7)
    d[2:(nyear-1)] <- runif(nyear-2, d[nyear], d[1])
    d[2:(nyear-1)] <- sort(d[2:(nyear-1)], decreasing = TRUE)
  }
  
  for (i in 1:(nspecies+1)) {
    a[,i] <- rnorm(nsite, 0, 3)
    b[,i] <- rnorm(nyear, 0, 3)
    p_detect[i] <- rbeta(1, 2, 2)
  }
  
  p_detect[1] <- 0.5  #fixed prob of detection for focal
  
  if (decline) {b[,1] <- d}
  
  for (j in 1:nsite) {
    for (t in 1:nyear) {
      lambda[j,t,] <- exp(a[j,]+b[t,])
      f1 <- function(x) rpois(1,x)
      y[j,t,] <- sapply(lambda[j,t,], f1)
      if (nb) {
        f11 <- function(x, r=0.425) rnbinom(1,x,r)
        y[j,t,] <- sapply(lambda[j,t,], f11)
      }
      richness <- length(which(y[j,t,]!=0))/(nspecies+1)
      
      nvisits <- rbinom(1, mv, richness)
      
      f2 <- function(x, p_detect) rbinom(nvisits, x, p_detect)
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
  
  for (t in 1:(nyear-1)) {
    prop <- 0.5-(t-1)/((nyear-1)/0.5)
    sub <- subset(simdata, year==t)
    vis <- unique(sub[,c("visit", "site")])
    vis$bias <- 2
    for (v in 1:nrow(vis)) {
      if (length(subset(match_df(sub, vis[v,1:2]), species==1)$actual>0)==1) {
        if (subset(match_df(sub, vis[v,1:2]), species==1)$actual>0) {
          vis$bias[v] <- 1
        }
      }
    }
    sel <- vis[sample(nrow(vis), prop*nrow(vis), prob = vis$bias),]
    for (i in 1:nrow(sel)) {
      sub <- match_df(sub, sel[i,1:2])
      simdata <- simdata[-as.integer(rownames(sub)),]
      sub <- subset(simdata, year==t)
    }
  }
  
  out <- list()
  out[["data"]] <- simdata
  out[["params"]] <- b
  
  return(out)
}

