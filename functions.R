require(reshape2)
require(R2jags)
require(zoo)
require(dplyr)
require(plyr)
require(MASS)

match_df <- function (x, y, on = NULL) {
  if (is.null(on)) {
    on <- intersect(names(x), names(y))
  }
  keys <- join.keys(x, y, on)
  x[keys$x %in% keys$y, , drop = FALSE]
}

###############################################

# with visits
create_data2 <- function(nspecies=20, nsite=50, nyear=10, mv=10, decline=FALSE, nb=FALSE, var_params=2){
  simdata <- data.frame()
  a <- matrix(nrow = nsite, ncol=nspecies+1)
  b <- matrix(nrow = nyear, ncol=nspecies+1)
  p_detect <- rep(0, nspecies+1)
  lambda <- array(dim=c(nsite, nyear, nspecies+1))
  y <- array(dim=c(nsite, nyear, nspecies+1))
  
  for (i in 1:(nspecies+1)) {
    a[,i] <- rnorm(nsite, 0, var_params)
    b[,i] <- rnorm(nyear, 0, var_params)
    p_detect[i] <- runif(1,0.16,0.88)
  }
  
  p_detect[1] <- 0.5  #fixed prob of detection for focal
  
  for (j in 1:nsite) {
    for (t in 1:nyear) {
      while (sum(a[j,]+b[t,]>12.2)>0) {
          a[j,] <- rnorm(nspecies+1, 0, var_params)
          b[t,] <- rnorm(nspecies+1, 0, var_params)
      }
    }
  }
  
  for (j in 1:nsite) {
    lambda[j,1,1] <- exp(a[j,1]+b[1,1])
    if (decline) {
      y[j,1,1] <- rpois(1, lambda[j,1,1])
      if (nb) {
        theta <- rnorm(1, 0.425, 8.82E-4)
        y[j,1,1] <- rnegbin(1, mu=lambda[j,1,1], theta)
      }
      y[j,nyear,1] <- round(0.7*y[j,1,1])
      y[j,2:(nyear-1),1] <- round(runif(nyear-2, y[j,nyear,1], y[j,1,1]))
      y[j,2:(nyear-1),1] <- sort(y[j,2:(nyear-1),1], decreasing = TRUE)
    }  else {for (t in 1:nyear) {
        lambda[j,t,1] <- exp(a[j,1]+b[t,1])
        f1 <- function(x) rpois(1,x)
        y[j,t,1] <- sapply(lambda[j,t,1], f1)
        if (nb) {
          theta <- rnorm(1, 0.425, 8.82E-4)
          f11 <- function(x) rnegbin(1,mu=x,theta)
          y[j,t,1] <- sapply(lambda[j,t,1], f11)
        }
    }}
    for (t in 1:nyear) {
        lambda[j,t,2:(nspecies+1)] <- exp(a[j,2:(nspecies+1)]+b[t,2:(nspecies+1)])
        f1 <- function(x) rpois(1,x)
        y[j,t,2:(nspecies+1)] <- sapply(lambda[j,t,2:(nspecies+1)], f1)
        if (nb) {
          theta <- rnorm(1, 0.425, 8.82E-4)
          f11 <- function(x) rnegbin(1,x,theta)
          y[j,t,2:(nspecies+1)] <- sapply(lambda[j,t,2:(nspecies+1)], f11)
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
create_data3 <- function(nspecies=20, nsite=50, nyear=10, mv=10, decline=FALSE, nb=FALSE, var_params=2){
  simdata <- data.frame()
  a <- matrix(nrow = nsite, ncol=nspecies+1)
  b <- matrix(nrow = nyear, ncol=nspecies+1)
  p_detect <- matrix(nrow = nyear, ncol=nspecies+1)
  lambda <- array(dim=c(nsite, nyear, nspecies+1))
  y <- array(dim=c(nsite, nyear, nspecies+1))
  
  if (decline) {
    # simulating a 30% decline in focal occupancy over 10 years
    d <- rep(0,nyear)
    d[1] <- rnorm(1, 0, var_params)
    d[nyear] <- d[1] + log(0.7)
    d[2:(nyear-1)] <- runif(nyear-2, d[nyear], d[1])
    d[2:(nyear-1)] <- sort(d[2:(nyear-1)], decreasing = TRUE)
  }
  
  for (i in 1:(nspecies+1)) {
    a[,i] <- rnorm(nsite, 0, var_params)
    b[,i] <- rnorm(nyear, 0, var_params)
    p_detect[,i] <- runif(1,0.16,0.88)
  }
  
  # detection increases
  for(t in 1:nyear) {
    p_detect[t,1] <- 0.4+(t-1)/((nyear-1)/0.1)
  }
  
  for (j in 1:nsite) {
    for (t in 1:nyear) {
      while (sum(a[j,]+b[t,]>12.2)>0) {
          a[j,] <- rnorm(nspecies+1, 0, var_params)
          b[t,] <- rnorm(nspecies+1, 0, var_params)
      }
    }
  }
  
  #if (decline) {b[,1] <- d}
  
  for (j in 1:nsite) {
    lambda[j,1,1] <- exp(a[j,1]+b[1,1])
    if (decline) {
      y[j,1,1] <- rpois(1, lambda[j,1,1])
      if (nb) {
        theta <- rnorm(1, 0.425, 8.82E-4)
        y[j,1,1] <- rnegbin(1, lambda[j,1,1], theta)
      }
      y[j,nyear,1] <- round(0.7*y[j,1,1])
      y[j,2:(nyear-1),1] <- round(runif(nyear-2, y[j,nyear,1], y[j,1,1]))
      y[j,2:(nyear-1),1] <- sort(y[j,2:(nyear-1),1], decreasing = TRUE)
    }  else {for (t in 1:nyear) {
      lambda[j,t,1] <- exp(a[j,1]+b[t,1])
      f1 <- function(x) rpois(1,x)
      y[j,t,1] <- sapply(lambda[j,t,1], f1)
      if (nb) {
        theta <- rnorm(1, 0.425, 8.82E-4)
        f11 <- function(x) rnegbin(1,x,theta)
        y[j,t,1] <- sapply(lambda[j,t,1], f11)
      }
    }}
    for (t in 1:nyear) {
      lambda[j,t,2:(nspecies+1)] <- exp(a[j,2:(nspecies+1)]+b[t,2:(nspecies+1)])
      f1 <- function(x) rpois(1,x)
      y[j,t,2:(nspecies+1)] <- sapply(lambda[j,t,2:(nspecies+1)], f1)
      if (nb) {
        theta <- rnorm(1, 0.425, 8.82E-4)
        f11 <- function(x) rnegbin(1,x,theta)
        y[j,t,2:(nspecies+1)] <- sapply(lambda[j,t,2:(nspecies+1)], f11)
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
create_data4 <- function(nspecies=20, nsite=50, nyear=10, mv=10, decline=FALSE, nb=FALSE, var_params=2){
  simdata <- data.frame()
  a <- matrix(nrow = nsite, ncol=nspecies+1)
  b <- matrix(nrow = nyear, ncol=nspecies+1)
  p_detect <- rep(0, nspecies+1)
  lambda <- array(dim=c(nsite, nyear, nspecies+1))
  y <- array(dim=c(nsite, nyear, nspecies+1))
  
  while (nrow(simdata)==0) {
    if (decline) {
      # simulating a 30% decline in focal occupancy over 10 years
      d <- rep(0,nyear)
      d[1] <- rnorm(1, 0, var_params)
      d[nyear] <- d[1] + log(0.7)
      d[2:(nyear-1)] <- runif(nyear-2, d[nyear], d[1])
      d[2:(nyear-1)] <- sort(d[2:(nyear-1)], decreasing = TRUE)
    }
    
    for (i in 1:(nspecies+1)) {
      a[,i] <- rnorm(nsite, 0, var_params)
      b[,i] <- rnorm(nyear, 0, var_params)
      p_detect[i] <- runif(1,0.16,0.88)
    }
    
    for (j in 1:nsite) {
      for (t in 1:nyear) {
        while (sum(a[j,]+b[t,]>12.2)>0) {
            a[j,] <- rnorm(nspecies+1, 0, var_params)
            b[t,] <- rnorm(nspecies+1, 0, var_params)
        }
      }
    }
    
    p_detect[1] <- 0.5  #fixed prob of detection for focal
    
    #if (decline) {b[,1] <- d}
    
    for (j in 1:nsite) {
      lambda[j,1,1] <- exp(a[j,1]+b[1,1])
      if (decline) {
        y[j,1,1] <- rpois(1, lambda[j,1,1])
        if (nb) {
          theta <- rnorm(1, 0.425, 8.82E-4)
          y[j,1,1] <- rnegbin(1, lambda[j,1,1], theta)}
        y[j,nyear,1] <- round(0.7*y[j,1,1])
        y[j,2:(nyear-1),1] <- round(runif(nyear-2, y[j,nyear,1], y[j,1,1]))
        y[j,2:(nyear-1),1] <- sort(y[j,2:(nyear-1),1], decreasing = TRUE)
      }  else {for (t in 1:nyear) {
        lambda[j,t,1] <- exp(a[j,1]+b[t,1])
        f1 <- function(x) rpois(1,x)
        y[j,t,1] <- sapply(lambda[j,t,1], f1)
        if (nb) {
          theta <- rnorm(1, 0.425, 8.82E-4)
          f11 <- function(x) rnegbin(1,x,theta)
          y[j,t,1] <- sapply(lambda[j,t,1], f11)
        }
      }}
      for (t in 1:nyear) {
        lambda[j,t,2:(nspecies+1)] <- exp(a[j,2:(nspecies+1)]+b[t,2:(nspecies+1)])
        f1 <- function(x) rpois(1,x)
        y[j,t,2:(nspecies+1)] <- sapply(lambda[j,t,2:(nspecies+1)], f1)
        if (nb) {
          theta <- rnorm(1, 0.425, 8.82E-4)
          f11 <- function(x) rnegbin(1,x,theta)
          y[j,t,2:(nspecies+1)] <- sapply(lambda[j,t,2:(nspecies+1)], f11)
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
  }
  
  out <- list()
  out[["data"]] <- simdata
  out[["params"]] <- b
  
  return(out)
}

##############################################

# visits double over 10 years
create_data5 <- function(nspecies=20, nsite=50, nyear=10, mv=10, decline=FALSE, nb=FALSE, var_params=2){
  simdata <- data.frame()
  a <- matrix(nrow = nsite, ncol=nspecies+1)
  b <- matrix(nrow = nyear, ncol=nspecies+1)
  p_detect <- rep(0, nspecies+1)
  lambda <- array(dim=c(nsite, nyear, nspecies+1))
  y <- array(dim=c(nsite, nyear, nspecies+1))
  
  if (decline) {
    # simulating a 30% decline in focal occupancy over 10 years
    d <- rep(0,nyear)
    d[1] <- rnorm(1, 0, var_params)
    d[nyear] <- d[1] + log(0.7)
    d[2:(nyear-1)] <- runif(nyear-2, d[nyear], d[1])
    d[2:(nyear-1)] <- sort(d[2:(nyear-1)], decreasing = TRUE)
  }
  
  for (i in 1:(nspecies+1)) {
    a[,i] <- rnorm(nsite, 0, var_params)
    b[,i] <- rnorm(nyear, 0, var_params)
    p_detect[i] <- runif(1,0.16,0.88)
  }
  
  p_detect[1] <- 0.5  #fixed prob of detection for focal
  
  for (j in 1:nsite) {
    for (t in 1:nyear) {
      while (sum(a[j,]+b[t,]>12.2)>0) {
          a[j,] <- rnorm(nspecies+1, 0, var_params)
          b[t,] <- rnorm(nspecies+1, 0, var_params)
      }
    }
  }
  
  #if (decline) {b[,1] <- d}
  
  for (j in 1:nsite) {
    lambda[j,1,1] <- exp(a[j,1]+b[1,1])
    if (decline) {
      y[j,1,1] <- rpois(1, lambda[j,1,1])
      if (nb) {
        theta <- rnorm(1, 0.425, 8.82E-4)
        y[j,1,1] <- rnegbin(1, lambda[j,1,1], theta)}
      y[j,nyear,1] <- round(0.7*y[j,1,1])
      y[j,2:(nyear-1),1] <- round(runif(nyear-2, y[j,nyear,1], y[j,1,1]))
      y[j,2:(nyear-1),1] <- sort(y[j,2:(nyear-1),1], decreasing = TRUE)
    }  else {for (t in 1:nyear) {
      lambda[j,t,1] <- exp(a[j,1]+b[t,1])
      f1 <- function(x) rpois(1,x)
      y[j,t,1] <- sapply(lambda[j,t,1], f1)
      if (nb) {
        theta <- rnorm(1, 0.425, 8.82E-4)
        f11 <- function(x) rnegbin(1,x,theta)
        y[j,t,1] <- sapply(lambda[j,t,1], f11)
      }
    }}
    for (t in 1:nyear) {
      lambda[j,t,2:(nspecies+1)] <- exp(a[j,2:(nspecies+1)]+b[t,2:(nspecies+1)])
      f1 <- function(x) rpois(1,x)
      y[j,t,2:(nspecies+1)] <- sapply(lambda[j,t,2:(nspecies+1)], f1)
      if (nb) {
        theta <- rnorm(1, 0.425, 8.82E-4)
        f11 <- function(x) rnegbin(1,x,theta)
        y[j,t,2:(nspecies+1)] <- sapply(lambda[j,t,2:(nspecies+1)], f11)
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
create_data6 <- function(nspecies=20, nsite=50, nyear=10, mv=10, decline=FALSE, nb=FALSE, var_params=2){
  simdata <- data.frame()
  a <- matrix(nrow = nsite, ncol=nspecies+1)
  b <- matrix(nrow = nyear, ncol=nspecies+1)
  p_detect <- rep(0, nspecies+1)
  lambda <- array(dim=c(nsite, nyear, nspecies+1))
  y <- array(dim=c(nsite, nyear, nspecies+1))
  
  if (decline) {
    # simulating a 30% decline in focal occupancy over 10 years
    d <- rep(0,nyear)
    d[1] <- rnorm(1, 0, var_params)
    d[nyear] <- d[1] + log(0.7)
    d[2:(nyear-1)] <- runif(nyear-2, d[nyear], d[1])
    d[2:(nyear-1)] <- sort(d[2:(nyear-1)], decreasing = TRUE)
  }
  
  for (i in 1:(nspecies+1)) {
    a[,i] <- rnorm(nsite, 0, var_params)
    b[,i] <- rnorm(nyear, 0, var_params)
    p_detect[i] <- runif(1,0.16,0.88)
  }
  
  p_detect[1] <- 0.5  #fixed prob of detection for focal
  
  for (j in 1:nsite) {
    for (t in 1:nyear) {
      while (sum(a[j,]+b[t,]>12.2)>0) {
          a[j,] <- rnorm(nspecies+1, 0, var_params)
          b[t,] <- rnorm(nspecies+1, 0, var_params)
      }
    }
  }
  
  #if (decline) {b[,1] <- d}
  
  for (j in 1:nsite) {
    lambda[j,1,1] <- exp(a[j,1]+b[1,1])
    if (decline) {
      y[j,1,1] <- rpois(1, lambda[j,1,1])
      if (nb) {
        theta <- rnorm(1, 0.425, 8.82E-4)
        y[j,1,1] <- rnegbin(1, lambda[j,1,1], theta)}
      y[j,nyear,1] <- round(0.7*y[j,1,1])
      y[j,2:(nyear-1),1] <- round(runif(nyear-2, y[j,nyear,1], y[j,1,1]))
      y[j,2:(nyear-1),1] <- sort(y[j,2:(nyear-1),1], decreasing = TRUE)
    }  else {for (t in 1:nyear) {
      lambda[j,t,1] <- exp(a[j,1]+b[t,1])
      f1 <- function(x) rpois(1,x)
      y[j,t,1] <- sapply(lambda[j,t,1], f1)
      if (nb) {
        theta <- rnorm(1, 0.425, 8.82E-4)
        f11 <- function(x) rnegbin(1,x,theta)
        y[j,t,1] <- sapply(lambda[j,t,1], f11)
      }
    }}
    for (t in 1:nyear) {
      lambda[j,t,2:(nspecies+1)] <- exp(a[j,2:(nspecies+1)]+b[t,2:(nspecies+1)])
      f1 <- function(x) rpois(1,x)
      y[j,t,2:(nspecies+1)] <- sapply(lambda[j,t,2:(nspecies+1)], f1)
      if (nb) {
        theta <- rnorm(1, 0.425, 8.82E-4)
        f11 <- function(x) rnegbin(1,x,theta)
        y[j,t,2:(nspecies+1)] <- sapply(lambda[j,t,2:(nspecies+1)], f11)
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

##############################################

# run abundance model
run_model <- function(simdata, species_list=1, model="SI_model.bug", n_chains=3, n_iter=5000) {
  output <- list()
  missing <- list()
  
  for (i in species_list) {
    simdata1 <- subset(simdata, species==i)
    if (nrow(simdata1)==0) {
      print("No records.") 
    }
    else {
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
                          model.file=model, 
                          n.chains=n_chains, n.iter=n_iter,
                          progress.bar = "none")
    }
  }
  
  return(list(output, missing))
}

##############################################

# function to run 100 simulations and record error rate/power
# using abundance model
assess_model <- function(nsims=100, scenarios="ABCDE", species_list=1, model="SI_model.bug", n_chains=3, n_iter=5000, nb=FALSE) {
  results <- list()
  create <- list(create_data2, create_data3, create_data4,
                 create_data5, create_data6)
  
  for (s in LETTERS[1:5]) {
    
    if (grepl(s, scenarios)) {
      test1 <- c()
      test2 <- c()
      miss <- c()
      
      for (i in 1:nsims) {
        data <- create[[which(LETTERS==s)]](decline=FALSE, nb=nb)
        simdata <- data[[1]]
        
        out <- run_model(simdata, species_list = species_list, model = model, n_chains=n_chains, n_iter=n_iter)
        
        if(sum(lengths(out))==0) {
          miss[i] <- NA
          test1[i] <- NA
          print(paste0("Scenario ", s, " + no trend: simulation ", i, " completed."))
        }
        
        else {
          miss[i] <- out[[2]][[1]]
        
          nyear <- formals(create[[which(LETTERS==s)]])$nyear
          est <- out[[1]][[1]]$BUGSoutput$summary[paste0("b[", 1:nyear, "]"), c(1,3,7,8)]
          est <- as.data.frame(est)
          est$year <- 1:nyear
          
          lntrend <- summary(lm(mean~year, data=est))
          test1[i] <- lntrend$coefficients[2,4]<0.05
          print(paste0("Scenario ", s, " + no trend: simulation ", i, " completed."))
        }
      }
        
     for (i in 1:nsims) {
        data <- create[[which(LETTERS==s)]](decline=TRUE, nb=nb)
        simdata <- data[[1]]
        
        out <- run_model(simdata, species_list = species_list, model = model, n_chains=n_chains, n_iter=n_iter)
        
        if(sum(lengths(out))==0) {
          miss[nsims+i] <- NA
          test2[i] <- NA
          print(paste0("Scenario ", s, " + decline: simulation ", i, " completed."))
        }
        
        else {
          miss[nsims+i] <- out[[2]][[1]]
          
          nyear <- formals(create[[which(LETTERS==s)]])$nyear
          est <- out[[1]][[1]]$BUGSoutput$summary[paste0("b[", 1:nyear, "]"), c(1,3,7,8)]
          est <- as.data.frame(est)
          est$year <- 1:nyear
          
          lntrend <- summary(lm(mean~year, data=est))
          test2[i] <- lntrend$coefficients[2,4]<0.05
          print(paste0("Scenario ", s, " + decline: simulation ", i, " completed."))
        }
      }
      
      test1 <- na.omit(test1)
      test2 <- na.omit(test2)
      errorI <- sum(test1)/length(test1)
      errorII <- 1 - sum(test2)/length(test2)
      power <- max(0, 1-(errorI+errorII))
      avgmiss <- mean(na.omit(miss))
      
      results[[s]] <- list(alpha=errorI, beta=errorII, power=power, avgmiss=avgmiss)
    }
  }
  
  return(results)
}

##############################################

# using occupancy detection model
assess_occmodel <- function(nsims=100, scenarios="ABCDE", species_list=1, 
                            model="SI_model.bug", n_chains=3, n_iter=5000, 
                            nb=FALSE) {
  results <- list()
  create <- list(create_data2, create_data3, create_data4,
                 create_data5, create_data6)
  
  for (s in LETTERS[1:5]) {
    
    if (grepl(s, scenarios)) {
      test1 <- c()
      test2 <- c()

      for (i in 1:nsims) {
        data <- create[[which(LETTERS==s)]](decline=FALSE, nb=nb)
        simdata <- data[[1]]
        
        vis <- unique(simdata[, c("visit", "year")])
        vis$dates <- 0
        nyear <- length(unique(simdata$year))
        sub1 <- data.frame()
        
        for (t in 1:nyear) {
          sub <- subset(vis, year==t)
          date_begin <- as.Date(paste0(2000+t, "-01-01"))
          date_end <- as.Date(paste0(2000+t,"-12-31"))
          rdates <- date_begin + (runif(nrow(sub))*(date_end-date_begin))
          sub$dates <- sort(rdates)
          sub1 <- rbind(sub1, sub)
        }
        
        simdata2 <- subset(join(simdata, sub1), obs>0)
        
        if(nrow(subset(simdata2, species==1))==0) {
          print(paste0("Scenario ", s, " + no trend: simulation ", i, " completed."))
        } else {
          out <- occDetModel(simdata2$species, simdata2$site, simdata2$dates, species_list=1,
                             write_results = FALSE)
          
          est <- out[[1]]$BUGSoutput$summary[paste0("alpha.p[", 1:nyear, "]"), c(1,3,7,8)]
          est <- as.data.frame(est)
          est$year <- 1:nyear
          
          lntrend <- summary(lm(mean~year, data=est))
          test1[i] <- lntrend$coefficients[2,4]<0.05
          print(paste0("Scenario ", s, " + no trend: simulation ", i, " completed."))
        }
      }
      
      for (i in 1:nsims) {
        data <- create[[which(LETTERS==s)]](decline=TRUE, nb=nb)
        simdata <- data[[1]]
        
        vis <- unique(simdata[, c("visit", "year")])
        vis$dates <- 0
        nyear <- length(unique(simdata$year))
        sub1 <- data.frame()
        
        for (t in 1:nyear) {
          sub <- subset(vis, year==t)
          date_begin <- as.Date(paste0(2000+t, "-01-01"))
          date_end <- as.Date(paste0(2000+t,"-12-31"))
          rdates <- date_begin + (runif(nrow(sub))*(date_end-date_begin))
          sub$dates <- sort(rdates)
          sub1 <- rbind(sub1, sub)
        }
        
        simdata2 <- subset(join(simdata, sub1), obs>0)
        
        if(nrow(subset(simdata2, species==1))==0) {
          print(paste0("Scenario ", s, " + no trend: simulation ", i, " completed."))
        } else {
          out <- occDetModel(simdata2$species, simdata2$site, simdata2$dates, species_list=1,
                             write_results = FALSE)
          
          est <- out[[1]]$BUGSoutput$summary[paste0("alpha.p[", 1:nyear, "]"), c(1,3,7,8)]
          est <- as.data.frame(est)
          est$year <- 1:nyear
          
          lntrend <- summary(lm(mean~year, data=est))
          test2[i] <- lntrend$coefficients[2,4]<0.05
          print(paste0("Scenario ", s, " + decline: simulation ", i, " completed."))
        }
      }
      
      test1 <- na.omit(test1)
      test2 <- na.omit(test2)
      errorI <- sum(test1)/length(test1)
      errorII <- 1 - sum(test2)/length(test2)
      power <- max(0, 1-(errorI+errorII))

      results[[s]] <- list(alpha=errorI, beta=errorII, power=power, avgmiss=avgmiss)
    }
  }
  
  return(results)
}
