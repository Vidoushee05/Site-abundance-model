source("functions.R")
library(ggplot2)
load("data1.Rdata")
library(knitr)
library(kableExtra)

par(mfrow=c(2,3))
for (i in seq(1,3.5,by=0.5)) {
  histdata <- list()
  for (f in 1:10) {
    data <- create_data2(decline=FALSE, var_params = i, nb=FALSE)
    simdata <- data[[1]]
    histdata[[f]] <- simdata$actual
  }
  min1 <- min(unlist(histdata))
  max1 <- max(unlist(histdata))
  col <- col2rgb(sample(colors(),10))
  transp <- rgb(col[1,], col[2,], col[3,], max = 255, alpha = 125)
  ax <- pretty(min1:max1, n=5)
  plot(hist(histdata[[1]], breaks=150, plot=FALSE), ylim=c(0,max(lengths(histdata))), 
       xlab = "Actual abundance",
       main=paste0("Variance=", i), col=transp[1])
  for (f in 2:10) {
    plot(hist(histdata[[f]], breaks=150, plot=FALSE), col=transp[f], add=TRUE)
  }
}

dev.off()

histdata <- list()

for (i in 1:6) {
  x <- seq(1,3.5,by=0.5)[i]
  histdata[[i]] <- list()
  for (f in 1:10) {
    data <- create_data2(decline=FALSE, var_params = x, nb=TRUE)
    simdata <- data[[1]]
    histdata[[i]][[f]] <- simdata$actual
  }
}

data2 <- subset(data1, common_name!="Large Blue")
sum <- summary(data1$site_index)[-7]
sum <- rbind(sum, summary(data2$site_index)[-7])
for (i in 1:6) {
  sum <- rbind(sum, summary(unlist(histdata[[i]])))
}

rownames(sum) <- c("UKBMS", "UKBMS (excluding Large Blue)", "Variance=1",
                   "Variance=1.5", "Variance=2", "Variance=2.5", "Variance=3",
                   "Variance=3.5")

sum[,4] <- round(sum[,4],2)

kable(sum, "latex")

data <- create_data2(decline=TRUE, var_params = 2, nb=FALSE)
simdata <- data[[1]]

x <- rep(0,10)
for (i in 1:10) {x[i] <- mean(subset(simdata2, species==1 & year==i)$actual)}
plot(x, type="l", xlab = "Year", ylab = "Mean abundance")


# probability of presence of focal
prob_occ <- rep(0,10)

data <- create_data2(decline=TRUE)
simdata <- data[[1]]

sub <- unique(subset(simdata[,c(2,4:6)], species==1))
for (t in 1:10) {
    sub1 <- subset(sub, year==t)
    prob_occ[t] <- proportion(sub1$actual>0)[2]
    x[t] <- mean(sub1$actual)
}

par(mfrow=c(1,2))
plot(x, type="l", xlab = "Year", ylab = "Mean abundance")
plot(prob_occ, type="l", xlab = "Year", ylab = "Probability of presence")

dev.off()

#################

# increased visits

data <- create_data2(decline=FALSE, var_params = 2, nb=FALSE)
simdata <- data[[1]]

records <- rep(0,10)

for (t in 1:10) {
  sub <- subset(simdata, year==t)
  records[t] <- nrow(sub)
}

plot(records, type="l", xlab="Year", ylab = "Number of records", ylim=c(2000, 7500))

col <- col2rgb(sample(colors(),10))


  data <- create_data5(decline=FALSE, var_params = 2, nb=FALSE)
  simdata <- data[[1]]
  
  records <- rep(0,10)
  
  for (t in 1:10) {
    sub <- subset(simdata, year==t)
    records[t] <- nrow(sub)
  }
  points(records, type="l", col = "red")
  
  legend("bottomright", legend = c("Control", "Increased Visits"), 
         col=c("black", "red"), lty=1)




