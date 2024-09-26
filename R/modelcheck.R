#### ---- setup -------
#load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_simp.rdata")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_longrun.rdata")
#load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_longrun.rdata")
load("data/data.rdata")
library ('MCMCvis')
library ('coda')
library ('ggplot2')
library('reshape2')
library('bayestestR')
out <- lapply(post, as.mcmc)

# Identify chains with NAs that 
# failed to initialize
NAlist <- c()
for (i in 1:length(out)){
  NAlist[i] <- any (is.na(out[[i]][,1:286]) | out[[i]][,1:286]<0)
}
# Subset chains to those with good initial values
out <- out[!NAlist]
post2 <- post[!NAlist]
outp <- MCMCpstr(out, type="chains")

!NAlist
#### ---- pltfunction -------
# default settings for plots 
plt  <- function(object, params,...) {
  MCMCplot(object=out, 
           params=params, 
           guide_axis=TRUE, 
           HPD=TRUE, ci=c(80, 95), horiz=FALSE, 
           #ylim=c(-10,10),
           ...)
  }

#### ---- catplots1 -------
# Abundance of females at Los Haitises
par(mfrow=c(5,2))
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("NFY[",1:13, ", 1]"), 
    main="First-year (FY) minus hacked\n Los Haitises", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plot(2011:2023, datl$countsFY[,1], 
     ylab="Counts", xlab="Year", type="b")

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("NF[",1:13, ", 1]"), 
    main="Adult nonbreeder (NB)\n Los Haitises", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plot.new()

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("NB[",1:13, ", 1]"),
    main="Adult breeder (B)\n Los Haitises", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plot.new()

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("NAD[",1:13, ", 1]"), 
    main="Adult Breeders and Nonbreeders\n Los Haitises", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plot(2011:2023, datl$countsAdults[,1], 
     ylab="Counts", xlab="Year",  type="b")

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("Ntot[",1:13, ", 1]"), 
    main="All stages\n Los Haitises", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plot(2011:2023, datl$countsFY[,1]+datl$countsAdults[,1], 
     ylab="Counts", xlab="Year",  type="b")

# Abundance of females at Punta Cana
par(mfrow=c(5,2))
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("NFY[",1:13, ", 2]"), 
    main="First-year (FY) plus hacked\n Punta Cana", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plot(2011:2023, datl$countsFY[,2], 
     ylab="Counts", xlab="Year",  type="b")

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("NF[",1:13, ", 2]"), 
    main="Adult nonbreeder (NB)\n Punta Cana", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plot.new()

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("NB[",1:13, ", 2]"),
    main="Adult breeder (B)\n Punta Cana", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plot.new()

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("NAD[",1:13, ", 2]"), 
    main="Adult Breeders and Nonbreeders\n Punta Cana", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plot(2011:2023, datl$countsAdults[,2], 
     ylab="Counts", xlab="Year",  type="b")

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("Ntot[",1:13, ", 2]"), 
    main="All stages\n Punta Cana", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plot(2011:2023, datl$countsFY[,2]+datl$countsAdults[,2], 
     ylab="Counts", xlab="Year",  type="b")

#### ---- catplots2 -------
# Finer population segments
par(mfrow=c(4,2))
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("N[", 1, ", ", 1:13, ", 1]"), 
    main="Los Haitises\nFirst-years fledged", 
    labels = 2011:2023, 
    xlab = "Year", ylab= "Abundance")
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("N[", 2, ", ", 1:13, ", 1]"), 
    main="Los Haitises\nFY to NB", 
    labels = 2011:2023, 
    xlab = "Year", ylab= "Abundance")
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("N[", 3, ", ", 1:13, ", 1]"), 
    main="Los Haitises\nNB to NB", 
    labels = 2011:2023, 
    xlab = "Year", ylab= "Abundance")
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("N[", 4, ", ", 1:13, ", 1]"), 
    main="Los Haitises\nB to NB", 
    labels = 2011:2023, 
    xlab = "Year", ylab= "Abundance")
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("N[", 5, ", ", 1:13, ", 1]"), 
    main="Los Haitises\nFY to B", 
    labels = 2011:2023, 
    xlab = "Year", ylab= "Abundance")
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("N[", 6, ", ", 1:13, ", 1]"), 
    main="Los Haitises\nNB to B",
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("N[", 7, ", ", 1:13, ", 1]"), 
    main="Los Haitises\nB to B", 
    labels = 2011:2023, 
    xlab = "Year", ylab= "Abundance")


par(mfrow=c(4,2))
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("N[", 1, ", ", 1:13, ", 2]"), 
    main="Punta Cana\nFY fledged",
    labels = 2011:2023, 
    xlab = "Year", ylab= "Abundance")
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("N[", 2, ", ", 1:13, ", 2]"), 
    main="Punta Cana\nFY to NB", 
    labels = 2011:2023, 
    xlab = "Year", ylab= "Abundance")
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("N[", 3, ", ", 1:13, ", 2]"), 
    main="Punta Cana\nNB to NB", 
    labels = 2011:2023, 
    xlab = "Year", ylab= "Abundance")
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("N[", 4, ", ", 1:13, ", 2]"), 
    main="Punta Cana\nB to NB", 
    labels = 2011:2023, 
    xlab = "Year", ylab= "Abundance")
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("N[", 5, ", ", 1:13, ", 2]"), 
    main="Punta Cana\nFY to B", 
    labels = 2011:2023, 
    xlab = "Year", ylab= "Abundance")
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("N[", 6, ", ", 1:13, ", 2]"), 
    main="Punta Cana\nNB to B",
    labels = 2011:2023, 
    xlab = "Year", ylab= "Abundance")
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("N[", 7, ", ", 1:13, ", 2]"), 
    main="Punta Cana\nB to B", 
    labels = 2011:2023, 
    xlab = "Year", ylab= "Abundance")

#### ---- catplots3 -------
# I needed to abbreviate to save plot space
# FY=first-year, NB=Nonbreeder, B=Breeder
par(mfrow=c(1,2))
plt(object=out, 
    params=paste0("mus[",1:8, ", 1]"), 
    exact=TRUE, ISB=FALSE, 
    ylim=c(0,1),
    main="Overall means\n Los Haitises", 
    labels=c("FY survival", "NB survival", "B survival",
             "FY to B", "NB to B", "B to NB",
             "NB detection", "B detection")
    )

plt(object=out, 
    params=paste0("mus[",1:8, ", 2]"), 
    exact=TRUE, ISB=FALSE, 
    ylim=c(0,1),
    main="Overall means\n Punta Cana", 
    labels=c("FY survival", "NB survival", "B survival",
             "FY to B", "NB to B", "B to NB",
             "NB detection", "B detection"))

par(mfrow=c(1,1))
plt(object=out, 
    params="betas", 
    main= "Translocation effects",
    labels=c("FY survival", "NB survival", "B survival",
             "FY to B", "NB to B", "B to NB",
             "NB detection", "B detection"))


# Plot effort effects
plt(object=out, 
    params=paste0("deltas[",1:4,"]"),
    exact=TRUE, ISB=FALSE,
    main= "Effort effects", 
    labels=c("Number of Adults\nEffort", "Number of Adults\nEffort squared",
             "Number of FYs\nEffort", "Number of FYs\nEffort squared"),
    ylim=c(-0.5, 1))

plt(object=out, 
    params=paste0("deltas[",5:8,"]"),
    exact=TRUE, ISB=FALSE,
    main= "Effort effects continued", 
    labels=c("Nonbreeder detection\nEffort", "Nonbreeder detection\nEffort sq",
             "Breeder detection\nEffort", "Breeder detection\nEffort sq"),
    ylim=c(-5, 5))

# Fecundity
par(mfrow=c(1,1))
plt(object=out, 
    params=c("lmu.prod"), 
    labels= c("Productivity\n(log scale)\nLos Haitises",
              "Productivity\n(log scale)\nPunta Cana"))

f <- exp(outp$lmu.prod)

# Is fecundity at LHNP greater than PC
par(mfrow=c(1,1))
fdiff <- f[1,]-f[2,]
hist(fdiff, main="Productivity difference")
abline(v=0, lty=2)
# print probability of direction, similar to frequentist p-value
# so values <=0.025 and >=0.975
# Is the difference in fecundity >0 ?
mean(fdiff>0) 

# How many times greater is fecundity 
# at treated versus non-treated sites
# median(f.pred[1,]/f[1,])
# median(f.pred[2,]/f[2,])


# gamma = nest treatment effect on fecundity
par(mfrow=c(1,1))
plt(object=out, 
    params=c("gamma"), 
    main="Anti-Parasitic Fly\nTreatment Effects", ylim=c(0,3))

par(mfrow=c(1,1))
sds <- paste0("sds[", 1:9, "]")
plt(object=out, params=sds,
    exact=TRUE, ISB=FALSE,
    main="Temporal SDs (synchrony among sites)",
    labels=c("FY survival", "NB survival", "B survival",
             "FY to B", "NB to B", "B to NB",
             "NB detection", "B detection",
             "Productivity"))
sds2 <- paste0("sds2[", 1:9, "]")
plt(object=out, params=sds2,
    exact=TRUE, ISB=FALSE,
    main="Site-temporal SDs", 
    labels=c("FY survival", "NB survival", "B survival",
             "FY to B", "NB to B", "B to NB",
             "NB detection", "B detection",
             "Productivity"))
# Correlations among vital rates
# Plot is messy with only a few strong correlations
ind <- 1
Rs <- R2s <- c()
for (i in 1:(nrow(outp$R)-1)){
  for (j in (i+1):nrow(outp$R)){
  Rs[ind] <- paste0("R[",i,", ", j, "]")
  R2s[ind] <- paste0("R2[",i,", ", j, "]")
  ind <- ind+1
  }}
par(mfrow=c(2,1))
plt(object=out, params=Rs[1:18], exact=TRUE, ISB=FALSE,
    main="Correlations btw demographic rates\n over time (synchrony)",
    xlab = "Rhos", guide_lines=TRUE)
plt(object=out, params=Rs[19:36], exact=TRUE, ISB=FALSE,
    main="Correlations btw demographic rates\n over time (synchrony), continued...",
    xlab = "Rhos", guide_lines=TRUE)
par(mfrow=c(2,1))
plt(object=out, params=R2s[1:18], exact=TRUE, ISB=FALSE, 
    main="Correlations btw demographic rates\n over time and sites",
    xlab = "Rhos", guide_lines=TRUE)
plt(object=out, params=R2s[19:36], exact=TRUE, ISB=FALSE, 
    main="Correlations btw demographic rates\n over time and sites, continued ...",
    xlab = "Rhos", guide_lines=TRUE)

# Annual averages for integration into the population model
par(mfrow=c(1,1))
labs <- c(paste0("LH ",2011:2023), paste0("PC ",2011:2023))
plt(object=out, params="mn.phiFY", ylim=c(0,1),
    main="First-year survival", labels = labs,
    xlab = "Year", ylab= "Survival")
abline(v=13.5, lwd=2)
plt(object=out, params="mn.phiA", ylim=c(0,1),
    main="Adult nonbreeder", labels = labs,
    xlab = "Year", ylab= "Survival")
abline(v=13.5, lwd=2)
plt(object=out, params="mn.phiB", ylim=c(0,1),
    main="Breeder", labels = labs,
    xlab = "Year", ylab= "Survival")
abline(v=13.5, lwd=2)
plt(object=out, params="mn.psiFYB", ylim=c(0,1),
    main="First-year to breeder", labels = labs,
    xlab = "Year", ylab= "Recruitment")
abline(v=13.5, lwd=2)
plt(object=out, params="mn.psiAB", ylim=c(0,1),
    main="Adult nonbreeder to breeder", labels = labs,
    xlab = "Year", ylab= "Recruitment")
abline(v=13.5, lwd=2)
plt(object=out, params="mn.psiBA", ylim=c(0,1),
    main="Adult breeder to nonbreeder", labels = labs,
    xlab = "Year", ylab= "Recruitment")
abline(v=13.5, lwd=2)
plt(object=out, params="mn.pA", ylim=c(0,1),
    main="Nonbreeder", labels = labs,
    xlab = "Year", ylab= "Detection")
abline(v=13.5, lwd=2)
plt(object=out, params="mn.pB", ylim=c(0,1),
    main="Breeder", labels = labs,
    xlab = "Year", ylab= "Detection")
abline(v=13.5, lwd=2)
plt(object=out, params="mn.prod",
    main="", labels=labs,
    xlab = "Year", ylab= "Productivity")
abline(v=13.5, lwd=2)


#### ---- popgrowth --------
# plot population growth rates and
# correlations with demographics
# 2012 represents the 
# pop growth between breeding seasons of 2011-2012 etc.
# A population growth rate >1 means the population is growing
# while a pgr <1 means it's shrinking.
lam.m <- apply(outp$lambda, c(1,2), median)
lam.hdi <- apply(outp$lambda, c(1,2), HDInterval::hdi)
par(mfrow=c(1,2))
plot(2012:2023, lam.m[,1], type="b", pch=1, 
     ylab="Population growth rate", xlab="Year", 
     ylim=c(min(lam.hdi[,,1]), max(lam.hdi[,,1])),
     main="Los Haitises")
abline(h=1, lty=2)
segments(x0=2012:2023, x1=2012:2023, 
         y0 = lam.hdi[1,,1], y1= lam.hdi[2,,1])

plot(2012:2023, lam.m[,2], type="b", pch=1, lty=1,
     ylab="Population growth rate", xlab="Year",
     ylim=c(min(lam.hdi[,,2]), max(lam.hdi[,,2])),
     main="Punta Cana")
abline(h=1, lty=2)
segments(x0=2012:2023, x1=2012:2023, 
         y0 = lam.hdi[1,,2], y1= lam.hdi[2,,2])

# create a function to plot correlations
# between demographics and population growth rates
plot.cor <- function (lambda, x, x.lab, ind.x=1:12){
  # calculate the correlation coefficient 
  # over each iteration to propagate error
  cor.post <- array(NA, dim=dim(lambda)[-1])
  cor.est <- c(NA, NA)
  for (s in 1:2){
  for (i in 1:dim(lambda)[[3]]){
    cor.post[s,i] <- cor(lambda[1:12,s,i], x[ind.x,s,i])
  }}
  cor.df <- data.frame(median= apply(cor.post, 1, median) |> round(2),
             ldi=apply(cor.post, 1, HDInterval::hdi)[1,] |> round(2),
             hdi=apply(cor.post, 1, HDInterval::hdi)[2,] |> round(2),
             pd = apply(cor.post, 1, pd) |> round(2),
             row.names=c("LH", "PC"))
             
  
  lam.m <- apply(lambda, c(1,2), median)
  lam.hdi <- apply(lambda, c(1,2), HDInterval::hdi)
  x.m <- apply(x, c(1,2), median)
  x.hdi <- apply(x, c(1,2), HDInterval::hdi)

  par(mfrow=c(1,2))
  for (s in 1:2){
    x.lims <- c(min(x.hdi[,,s]), max(x.hdi[,,s]))
    y.lims <- c(min(lam.hdi[,,s]), max(lam.hdi[,,s]))
  plot(x.m[ind.x,s], lam.m[1:12,s], 
       xlim= x.lims,
       ylim= y.lims,
       type="n", ylab="Population growth rate", xlab=x.lab, 
       main=c("Los Haitises", "Punta Cana")[s])
  points(x.m[ind.x,s], lam.m[1:12,s], pch=1)
  segments(x0=x.hdi[1,ind.x,s], x1=x.hdi[2,ind.x,s], 
            y0 = lam.m[,s], y1= lam.m[,s])
  segments(x0=x.m[ind.x,s], x1=x.m[ind.x,s], 
            y0 = lam.hdi[1,,s], y1= lam.hdi[2,,s])
  text(x = x.lims[1], y = (y.lims[2]-y.lims[1])*0.9+y.lims[1], 
       paste("r = ", cor.df$median[s], " (",cor.df$ldi[s],", ", cor.df$hdi[s], ")", sep=""), 
       pos = 4, font = 3, cex = 1)
  text(x = x.lims[1], y = (y.lims[2]-y.lims[1])*0.8+y.lims[1], paste("P(r>0) = ", cor.df$pd[s], sep=""), 
       pos = 4, font = 3, cex = 1)
  }
}
# Plor correlations between population grrowth rates
# and demographics. 
# "r" is a correlation coefficient and represents
# the magnitude of the correlation
# P(r>0) is the probability of direction (similar to p-values)
# that is, the probability that an effect exists
plot.cor(outp$lambda, outp$mn.prod, x.lab="Productivity", ind.x=2:13)
plot.cor(outp$lambda, outp$mn.phiFY, x.lab="First-year Survival")
plot.cor(outp$lambda, outp$mn.phiA, x.lab="Nonbreeder Survival")
plot.cor(outp$lambda, outp$mn.phiB, x.lab="Breeder Survival")
plot.cor(outp$lambda, outp$mn.psiFYB, x.lab="First-year to Breeder")
plot.cor(outp$lambda, outp$mn.psiAB, x.lab="Nonbreeder to Breeder")
# Breeder to nonbreeder didn't vary over time
# Does the number of breeders correlate with pop growth?
plot.cor(outp$lambda, outp$NB, x.lab="Breeder Abundance", ind.x=2:13)
plot.cor(outp$lambda, outp$NF, x.lab="Nonbreeder Abundance", ind.x=2:13)
plot.cor(outp$lambda, outp$NFY, x.lab="First-year Abundance", ind.x=2:13)
plot.cor(outp$lambda, outp$Ntot, x.lab="All Stages Abundance", ind.x=2:13)

# Plot correlation between translocation 
# and population growth rates.
ind.x <- 1:12
lam.md <- apply(outp$lambda, c(1,2), median)
lam.hdis <- apply(outp$lambda, c(1,2), HDInterval::hdi)
plot(constl$hacked.counts[ind.x,1], lam.md[,1],
     xlab="Number translocated", 
     ylab="Population growth rate",
     type="n",
     xlim=c(min(constl$hacked.counts[,1]), max(constl$hacked.counts[,1])),
     ylim=c(min(lam.hdis[,,1]), max(lam.hdis[,,1])),
     main="Los Haitises")
points(constl$hacked.counts[ind.x,1], lam.md[,1])
segments(x0=constl$hacked.counts[ind.x,1], x1=constl$hacked.counts[ind.x,1], 
         y0 = lam.hdis[1,,1], y1= lam.hdis[2,,1])
plot(constl$hacked.counts[ind.x,2], lam.md[,2],
     xlab="Number translocated", 
     ylab="Population growth rate",
     type="n",
     xlim=c(min(constl$hacked.counts[,2]), max(constl$hacked.counts[,2])),
     ylim=c(min(lam.hdis[,,2]), max(lam.hdis[,,2])),
     main="Punta Cana")
points(constl$hacked.counts[ind.x,2], lam.md[,2])
segments(x0=constl$hacked.counts[ind.x,2], x1=constl$hacked.counts[ind.x,2], 
         y0 = lam.hdis[1,,2], y1= lam.hdis[2,,2])  


#### ---- paramests -------
pars1 <- c("sds", "sds2","mus", "betas",
           "NFY", "NF", "NB", "Ntot", 
           "mn.phiFY","mn.phiA", "mn.phiB", 
           "mn.psiFYB", "mn.psiAB", "mn.psiBA", 
           "mn.pA", "mn.pB")
# Estimates for the survival model 
# In this order: FY survival, NB survival, B survival, 
# FY to B recruitment, NB to B recruitent, B to NB recruitment,
# Detection NB, Detection B

MCMCsummary(post2, params = sds2, #c(sds, sds2), 
            exact=TRUE, ISB=FALSE,
            digits=2, HPD = T,
            hpd_prob = 0.95, pg0= TRUE,
            func=median, func_name="median")

# Mus are means for 
# mus[1, site] , where site=1 is LH and site=2 is PC
# Survival of first years = mus[1,]
# Survival of nonbreeders = mus[2,]
# Survival of breeders = mus[3,]
# Breeding propensity of first years = mus[4,]
# Breeding propensity of nonbreeders = mus[5,]
# Transition from breeder to nonbreeder = mus[6,]
# Detection probability of nonbreeders = mus[7,]
# Detection probability of breeders
# modeled as logit(probability) = lmus[x,site] + betas[x]*translocated + eps[x,t] + eta[x,s,t]
# except breeder to nonbreeder recruitment = lmus[6,site] 
MCMCsummary(post2, params = pars1[3], 
            digits=3, HPD = T, 
            hpd_prob = 0.95, pg0= TRUE,
            func=median, func_name="median")

MCMCsummary(post2, params = "lmus", 
            digits=3, HPD = T, 
            hpd_prob = 0.95, pg0= TRUE,
            func=median, func_name="median")

MCMCsummary(post2, params = pars1[4], 
            digits=3, HPD = T, 
            hpd_prob = 0.95, pg0= TRUE,
            func=median, func_name="median")

# Fecundity
# modeled as log(f) = lmu.f[site] + gamma*treatment + eps[x,t] + eta[x,s,t]
# log scale
MCMCsummary(post2, params = "lmu.prod", 
            digits=3, HPD = T, 
            hpd_prob = 0.95, pg0= TRUE,
            func=median, func_name="median")

# Estimates of population size
# NFY= first year, 
# NF = nonbreeders, 
# NB=breeders,
# Ntot=total
# All are presented as N[time, site], 
# where time 1=2011 ... 13=2023, 
# site 1 = LH and site 2=PC
MCMCsummary(post2, params = pars1[5:8], 
            digits=2, HPD = T, 
            hpd_prob = 0.95, pg0= TRUE,
            func=median, func_name="median")

MCMCsummary(post2, params = pars1[9:16], 
            digits=2, HPD = T, 
            hpd_prob = 0.95, pg0= TRUE,
            func=median, func_name="median")

MCMCsummary(post2, params = "deltas", 
            digits=2, HPD = T, 
            hpd_prob = 0.95, pg0= TRUE,
            func=median, func_name="median")
# Correlations among demographic rates time (synchrony)
MCMCsummary(post2, params = Rs,
            exact=TRUE, ISB=FALSE,
            digits=2, HPD = T,
            hpd_prob = 0.95, pg0= TRUE,
            func=median, func_name="median")
# Correlations among demographic rates site x time
MCMCsummary(post2, params = R2s, 
            exact=TRUE, ISB=FALSE,
            digits=2, HPD = T, 
            hpd_prob = 0.95, pg0= TRUE,
            func=median, func_name="median")

#### ---- traceplots ------
#MCMCtrace(post2, pdf=FALSE, params= "sds")
MCMCtrace(post2, pdf=FALSE, params= "sds2")
MCMCtrace(post2, pdf=FALSE, params= "mus")
MCMCtrace(post2, pdf=FALSE, params= "betas")
MCMCtrace(post2, pdf=FALSE, params= "NF")
MCMCtrace(post2, pdf=FALSE, params= "NFY")
MCMCtrace(post2, pdf=FALSE, params= "NB")
MCMCtrace(post2, pdf=FALSE, params= "R2")

#### ---- fit ------
# Goodness of fit check
fit.check <- function(out, ratio=FALSE, 
                      name.rep="f.dmape.rep", 
                      name.obs="f.dmape.obs",
                      jit=100,
                      ind=1,
                      lab=""){
  par(mfrow=c(1,1))
  # plot mean absolute percentage error
  samps <- MCMCpstr(out, "all", type="chains")
  rep <- samps[name.rep][[1]][ind,]
  obs <- samps[name.obs][[1]][ind,]
  mx <- max(c(rep, obs))
  mn <- min(c(rep, obs))
  plot(jitter(obs, amount=jit), 
       jitter(rep, amount=jit),
       main=paste0("Mean absolute percentage error\n",lab),
       ylab="Discrepancy replicate values",
       xlab="Discrepancy observed values", 
       xlim=c(mn,mx), ylim=c(mn,mx), 
       pch=16, cex=0.5, col="gray10")
  curve(1*x, from=mn, to=mx, add=T, lty=2, lwd=2, col="blue")
  bp1 <- round(mean(rep > obs),2)
  loc <- ifelse(bp1 < 0.5, "topleft", "bottomright")
  legend(loc, legend=bquote(p[B]==.(bp1)), bty="n", cex=2)
  
  if (ratio==TRUE){
    t.rep <- samps["tvm.rep"][[1]][ind,]
    t.obs <- samps["tvm.obs"][[1]][ind,]
    # plot variance/mean ratio
    hist(t.rep, nclass=50,
         xlab="variance/mean ", main=NA, axes=FALSE)
    abline(v=t.obs, col="red")
    axis(1); axis(2)
  }
  return(list('Bayesian p-value'=bp1))
}

# check goodness-of-fit for brood size
# breeder, ind=1
# fit.check(out, ratio=F,
#           name.rep="dmape.rep", 
#           name.obs="dmape.obs",
#           ind=1,
#           lab="Breeder counts- Poisson", jit=300)
# # nonbreeder, ind=2
# fit.check(out, ratio=F,
#           name.rep="dmape.rep", 
#           name.obs="dmape.obs",
#           ind=2,
#           lab="Nonbreeder counts- Poisson", jit=300)

fit.check(out, ratio=F,
          name.rep="dmape.rep", 
          name.obs="dmape.obs",
          ind=1,
          lab="Adults(Breeder+Nonbreeder)- Poisson", jit=300)
# first-year, ind=2
# poisson failed fit test bp=0
# Currently running models to try and fix
fit.check(out, ratio=F,
          name.rep="dmape.rep", 
          name.obs="dmape.obs",
          ind=2,
          lab="First-year counts\nNeg binomial-Poisson", jit=300)
# fecundity
fit.check(out, ratio=F,
          name.rep="f.dmape.rep", 
          name.obs="f.dmape.obs",
          ind=1,
          lab="Productivity-Neg binomial", jit=300)
