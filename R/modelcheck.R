#### ---- setup -------
#load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_short.rdata")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_da_longrun_2025_Apr_10.rdata")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_dd_longrun_2025_Apr_07.rdata")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_longrun_2025_Apr_01.rdata")
load("data/data-dd.rdata")
library ('MCMCvis')
library ('coda')
library ('ggplot2')
library('reshape2')
library('bayestestR')

# Identify chains with NAs that 
# failed to initialize
NAlist <- c()
for (i in 1:length(post)){
  NAlist[i] <- any (is.na(post[[i]][,1:311]) | post[[i]][,1:311]<0)
}
# Subset chains to those with good initial values
out <- post[!NAlist]
#out <- out[1:5] # omit chain 5 bc one parameter causes a convergence failure
!NAlist
outp <- MCMCpstr(out, type="chains")


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

# gamma = nest treatment effect on fecundity
par(mfrow=c(1,1))
plt(object=out, 
    params=c("gamma"), 
    main="Anti-Parasitic Fly\nTreatment Effects", ylim=c(0,3))

par(mfrow=c(1,1))
sds <- paste0("sds[", 1:9, "]")
plt(object=out, params=sds,
    exact=TRUE, ISB=FALSE,
    main="Site-temporal SDs", 
    labels=c("FY survival", "NB survival", "B survival",
             "FY to B", "NB to B", "B to NB",
             "NB detection", "B detection",
             "Productivity"))

#---- demo_avs ----
# Demographics- Annual Averages
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

#---- densitydependence ----
# Demographics- Annual Averages
# Annual averages for integration into the population model
plt(object=out, params="alphas",
    main="Coefficients for density dependence")

#### ---- paramests -------
# Estimates for the survival model 
# In this order: FY survival, NB survival, B survival, 
# FY to B recruitment, NB to B recruitent, B to NB recruitment,
# Detection NB, Detection B

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
MCMCsummary(out, params = "mus", 
            digits=3, HPD = T, 
            hpd_prob = 0.95, pg0= TRUE,
            func=median, func_name="median")

MCMCsummary(out, params = "betas", 
            digits=3, HPD = T, 
            hpd_prob = 0.95, pg0= TRUE,
            func=median, func_name="median")

# Fecundity
# modeled as log(f) = lmu.f[site] + gamma*treatment + eps[x,t] + eta[x,s,t]
# log scale
MCMCsummary(out, params = "lmu.prod", 
            digits=3, HPD = T, 
            hpd_prob = 0.95, pg0= TRUE,
            func=median, func_name="median")

# Effort effects on detections and resighting
MCMCsummary(out, params = "deltas", 
            digits=2, HPD = T, 
            hpd_prob = 0.95, pg0= TRUE,
            func=median, func_name="median")

# Random effects for temporal synchrony and site x year. 
MCMCsummary(out, params = sds, 
            exact=TRUE, ISB=FALSE,
            digits=2, HPD = T,
            hpd_prob = 0.95, pg0= TRUE,
            func=median, func_name="median")

# Correlations among demographic rates site x time
MCMCsummary(out, params = Rs, 
            exact=TRUE, ISB=FALSE,
            digits=2, HPD = T, 
            hpd_prob = 0.95, pg0= TRUE,
            func=median, func_name="median")

# Coeffieicents for density dependence
MCMCsummary(out, params = "alphas", 
            #exact=TRUE, ISB=FALSE,
            digits=2, HPD = T, 
            hpd_prob = 0.95, pg0= TRUE,
            func=median, func_name="median") |> round (3) 

#### ---- traceplots ------
out <- post
outp <- MCMCpstr(out, type="chains")
MCMCtrace(out, pdf=FALSE, params= "mus", Rhat=TRUE, priors=runif(20000, 0, 1), post_zm=FALSE)
MCMCtrace(out, pdf=FALSE, params= "lmu.prod", Rhat=TRUE, priors=rnorm(20000, 0, 5), post_zm=FALSE)
MCMCtrace(out, pdf=FALSE, params= "betas", Rhat=TRUE, priors=rnorm(20000, 0, 10), post_zm=FALSE)
MCMCtrace(out, pdf=FALSE, params= "deltas", Rhat=TRUE, priors=rnorm(20000, 0, 10), post_zm=FALSE)
MCMCtrace(out, pdf=FALSE, params= "gamma", Rhat=TRUE, priors=rnorm(20000, 0, 10), post_zm=FALSE)
MCMCtrace(out, pdf=FALSE, params= "alphas", Rhat=TRUE, priors=rnorm(20000, 0, 10), post_zm=FALSE)
MCMCtrace(out, pdf=FALSE, params= "sds", Rhat=TRUE, priors=rexp(20000, 1), post_zm=FALSE)
MCMCtrace(out, pdf=FALSE, params= "r", Rhat=TRUE, priors=rexp(20000, 0.05), post_zm=FALSE)
MCMCtrace(out, pdf=FALSE, params= "rr", Rhat=TRUE, priors=rexp(20000, 0.05), post_zm=FALSE)

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

# Counts of adults (nonbreeders+breeders)
# tiff(height=4, width=5, units="in", res=300,
#      filename= "C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//GOF_adults.tiff")
fit.check(out, ratio=F,
          name.rep="dmape.rep", 
          name.obs="dmape.obs",
          ind=1,
          lab="Adults(Breeder+Nonbreeder)- Poisson", jit=300)
# dev.off()
# Counts of first-year fledglings, ind=2
# Poisson failed fit test bp=0
# So we assigned negative binomial only for Punta Cana.
# Los Haitises excluded from neg bin because no zeroes in counts.
# tiff(height=4, width=5, units="in", res=300,
#      filename= "C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//GOF_FY.tiff")
fit.check(out, ratio=F,
          name.rep="dmape.rep", 
          name.obs="dmape.obs",
          ind=2,
          lab="First-year counts\nNeg binomial-Poisson", jit=300)
# dev.off()
# Productivity
# tiff(height=4, width=5, units="in", res=300,
#      filename= "C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//GOF_prod.tiff")
fit.check(out, ratio=F,
          name.rep="f.dmape.rep", 
          name.obs="f.dmape.obs",
          ind=1,
          lab="Productivity-Neg binomial", jit=300)
# dev.off()

