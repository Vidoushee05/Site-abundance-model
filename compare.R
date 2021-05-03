source("functions.R")
library(sparta)
library(gridExtra)
library(grid)
library(knitr)
library(kableExtra)

data <- create_data4(decline=TRUE, nb=TRUE)
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

out <- occDetModel(simdata2$species, simdata2$site, simdata2$dates, species_list=1,
                   write_results = FALSE)

est <- out[[1]]$BUGSoutput$summary[paste0("alpha.p[", 1:nyear, "]"), c(1,3,7,8)]
est <- as.data.frame(est)
est$year <- 1:nyear

summary(lm(mean~year, data=est))

out2 <- run_model(simdata, species_list = c(1),
                  n_chains=3, n_iter=10000)
trend <- data[[2]]

actual <- trend[,1]

est <- out2[[1]][[1]]$BUGSoutput$summary[paste0("b[", 1:nyear, "]"), c(1,3,7,8)]
est <- as.data.frame(est)
est$year <- 1:nyear
est$actual <- actual

p1 <- ggplot(est, aes(year, mean)) + geom_line(aes(color="Estimated")) +
  geom_point(aes(color="Estimated")) +
  geom_ribbon(data=est, aes(ymin=`2.5%`, ymax=`97.5%`), alpha=0.3) + 
  geom_point(data=est, aes(year, actual, color="Actual")) + ylab("b_t") +
  scale_x_continuous(breaks=seq(1, 10, by=1)) +
  theme(legend.title = element_blank()) + ggtitle("Negative Binomial simulation")

occA <- assess_occmodel(nsims=100, scenarios = "A")
save(occA, file="results/Control_occ.rdata")

occB <- assess_occmodel(nsims=100, scenarios = "B")
save(occB, file="results/Detect_occ.rdata")

occC <- assess_occmodel(nsims=100, scenarios = "C")
save(occC, file="results/Effort_occ.rdata")

occD <- assess_occmodel(nsims=100, scenarios = "D")
save(occD, file="results/Visits_occ.rdata")

occE <- assess_occmodel(nsims=100, scenarios = "E")
save(occE, file="results/Visits+Bias_occ.rdata")

####################################################

reslist1 <- list(resultsA, resultsB, resultsC, resultsD, resultsE)
reslist2 <- list(resultsAnb, resultsBnb, resultsCnb, resultsDnb, resultsEnb)
reslist3 <- list(occA, occB, occC, occD, occE)
 
p1 <- plot_alpha(reslist1) + theme(axis.title.y = element_blank())
p2 <- plot_alpha(reslist2) + theme(legend.position = "none") +
  theme(axis.title.y = element_blank())
p3 <- plot_alpha(reslist3) + theme(legend.position = "none") +
  theme(axis.title.y = element_blank())

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

legend <- get_legend(p1)

p1 <- p1 + theme(legend.position = "none")

grid.arrange(arrangeGrob(p1, top = 'Abundance model+Poisson'), 
             arrangeGrob(p2, top = 'Abundance model+NB'),
             arrangeGrob(p3, top = 'Occupancy detection model'),
             legend, ncol=4, widths = c(2,2,2,1),
             left = textGrob("Type I error rate", rot = 90),
             bottom = textGrob("Scenario", vjust=-1))

reslist <- list(resultsA, resultsB, resultsC, resultsD, resultsE, 
                resultsAnb, resultsBnb, resultsCnb, resultsDnb, resultsEnb)

nam <- names(unlist(reslist))
avgmiss <- unlist(reslist)[which(grepl("avgmiss", nam))]

x <- data.frame(names=names(avgmiss), avgmiss)
x$names <- c("Control", "Detect", "Effort", "Visits", "Visits+Bias",
             "Control+NB", "Detect+NB", "Effort+NB", "Visits+NB", "Visits+Bias+NB")
x$names <- factor(x$names, levels = c("Control", "Control+NB", "Detect", "Detect+NB", 
                                      "Effort", "Effort+NB", 
                                      "Visits", "Visits+NB", "Visits+Bias", "Visits+Bias+NB"))
                                      

p <- ggplot(data=x, aes(x=names, y=avgmiss, fill=names, order=names)) + 
  geom_bar(stat="identity") + 
  ylim(0,0.075) +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  xlab("Scenario") + ylab("Missingness") + labs(fill="Scenario")

#######################################################

# results tables

reslist4 <- list(resultsA, resultsB, resultsC, resultsD, resultsE, 
                resultsAnb, resultsBnb, resultsCnb, resultsDnb, resultsEnb,
                occA, occB, occC, occD, occE)

nam <- names(unlist(reslist1))
alpha1 <- unlist(reslist1)[which(grepl("alpha", nam))]

nam <- names(unlist(reslist2))
alpha2 <- unlist(reslist2)[which(grepl("alpha", nam))]

nam <- names(unlist(reslist3))
alpha3 <- unlist(reslist3)[which(grepl("alpha", nam))]

x <- data.frame(Scenario = c("Control", "Detect", "Effort", "Visits", "Visits+Bias"))
x[,2] <- alpha1 
x[,3] <- alpha2
x[,4] <- alpha3

x[,2:4] <- round(x[,2:4], 2)
colnames(x) <- c("Scenario", "AM+P", "AM+NB", "ODM")

kable(x, "latex")

nam <- names(unlist(reslist1))
power1 <- unlist(reslist1)[which(grepl("power", nam))]

nam <- names(unlist(reslist2))
power2 <- unlist(reslist2)[which(grepl("power", nam))]

nam <- names(unlist(reslist3))
power3 <- unlist(reslist3)[which(grepl("power", nam))]

x <- data.frame(Scenario = c("Control", "Detect", "Effort", "Visits", "Visits+Bias"))
x[,2] <- power1 
x[,3] <- power2
x[,4] <- power3

x[,2:4] <- round(x[,2:4], 2)
colnames(x) <- c("Scenario", "AM+P", "AM+NB", "ODM")

kable(x, "latex")
