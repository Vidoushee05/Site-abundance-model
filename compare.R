source("functions.R")
library(sparta)

data <- create_data4(decline=TRUE)
simdata <- data[[1]]

vis <- unique(simdata[, c("visit", "year")])
vis$dates <- 0
nyear <- 10
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

#formatOccData(simdata2$species, simdata2$site, simdata2$dates)

out <- occDetModel(simdata2$species, simdata2$site, simdata2$dates, species_list=1)

est <- out[[1]]$BUGSoutput$summary[paste0("alpha.p[", 1:nyear, "]"), c(1,3,7,8)]
est <- as.data.frame(est)
est$year <- 1:nyear

summary(lm(mean~year, data=est))

out2 <- run_model(simdata, species_list = c(1), n_chains=3, n_iter=5000)
trend <- data[[2]]

actual <- trend[,1]

est <- out2[[1]][[1]]$BUGSoutput$summary[paste0("b[", 1:nyear, "]"), c(1,3,7,8)]
est <- as.data.frame(est)
est$year <- 1:nyear
est$actual <- actual

p <- ggplot(est, aes(year, mean)) + geom_line(aes(color="Estimated")) +
  geom_ribbon(data=est, aes(ymin=`2.5%`, ymax=`97.5%`), alpha=0.3) + 
  geom_point(data=est, aes(year, actual, color="Actual")) + ylab("b_t") +
  scale_x_continuous(breaks=seq(1, 10, by=1)) +
  theme(legend.title = element_blank()) 
