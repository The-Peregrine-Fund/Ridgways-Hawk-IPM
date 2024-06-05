#### ---- catplots -------
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_sites.rdata")
load("data/data.rdata")
library ('MCMCvis')
library ('coda')
out <- list(as.mcmc(post[[1]]), 
             as.mcmc(post[[2]]), 
             as.mcmc(post[[3]]),
             as.mcmc(post[[4]]))
outp <- MCMCpstr(out, type="chains")
# default settings for plots 
plt  <- function(object, params,...) {
  MCMCplot(object=out, 
           params=params, 
           guide_axis=TRUE, 
           HPD=TRUE, ci=c(80, 95), horiz=FALSE, ...)
  }

# I needed to abbreviate to save plot space
# FY=first-year, NB=Nonbreeder, B=Breeder
plt(object=out, params="sds", 
    main="Temporal SDs", 
    labels=c("FY survival", "NB survival", "B survival",
             "FY to B", "NB to B", "B to NB",
             "NB detection", "B detection",
             "Brood size", "Nest success"))
plt(object=out, params="sds2", 
    main="Site-temporal SDs", 
    labels=c("FY survival", "NB survival", "B survival",
             "FY to B", "NB to B", "B to NB",
             "NB detection", "B detection",
             "Brood size", "Nest success"))

plt(object=out, 
    params=paste0("mus[",1:13, ", 1]"), 
    exact=TRUE, ISB=FALSE, 
    ylim=c(0,1),
    main="Overall means\n Los Haitises", 
    labels=c("FY survival", "NB survival", "B survival",
             "FY to B", "NB to B", "B to NB",
             "NB detection", "B detection")
    )

plt(object=out, 
    params=paste0("mus[",1:13, ", 2]"), 
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

# Abundance of females at Los Haitises
par(mfrow=c(2,2))
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("NFY[",1:13, ", 1]"), 
    main="First-year (FY)\n Los Haitises", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("NF[",1:13, ", 1]"), 
    main="Adult nonbreeder (NB)\n Los Haitises", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("NB[",1:13, ", 1]"),
    main="Adult breeder (B)\n Los Haitises", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("Ntot[",1:13, ", 1]"), 
    main="All stages\n Los Haitises", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
# Abundance of females at Punta Cana
par(mfrow=c(2,2))
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("NFY[",1:13, ", 2]"), 
    main="First-year (FY)\n Punta Cana", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("NF[",1:13, ", 2]"), 
    main="Adult nonbreeder (NB)\n Punta Cana", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("NB[",1:13, ", 2]"),
    main="Adult breeder (B)\n Punta Cana", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("Ntot[",1:13, ", 2]"), 
    main="All stages\n Punta Cana", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")

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
plt(object=out, params=Rs, exact=TRUE, ISB=FALSE, 
    main="Temporal correlations\namong sites",
    xlab = "Rhos", guide_lines=TRUE)
plt(object=out, params=R2s, exact=TRUE, ISB=FALSE, 
    main="Site-temporal correlations",
    xlab = "Rhos", guide_lines=TRUE)
# lmu.brood = mean brood size (log scale), 
# sig.brood = SD among nests
# mu.nest = mean nest success
par(mfrow=c(1,1))
plt(object=out, 
    params=c("lmu.brood", "sig.brood", "mu.nest"), 
    labels= c("Brood size\n(log scale)\nLos Haitises",
              "Brood size\n(log scale)\nPunta Cana",
              "SD among females", 
              "Prob. of\nnest success\n Los Haitises", 
              "Prob. of\nnest success\n Punta Cana"))
# delta = nest treatment effect on brood size
# gamma = nest treatment effect on nest success
plt(object=out, 
    params=c("delta", "gamma"), 
    main="Treatment effects", 
    labels= c("Brood size", "Nest success"))

# Annual averages for integration into the population model
labs <- c(paste0("LH ",2011:2023), paste0("PC ",2011:2023))
plt(object=out, params="mn.phiFY", ylim=c(0,1),
    main="First-year survival", labels = labs,
    xlab = "Year", ylab= "Survival")
plt(object=out, params="mn.phiA", ylim=c(0,1),
    main="Adult nonbreeder", labels = labs,
    xlab = "Year", ylab= "Survival")
plt(object=out, params="mn.phiB", ylim=c(0,1),
    main="Breeder", labels = labs,
    xlab = "Year", ylab= "Survival")
plt(object=out, params="mn.psiFYB", ylim=c(0,1),
    main="First-year to breeder", labels = labs,
    xlab = "Year", ylab= "Recruitment")
plt(object=out, params="mn.psiAB", ylim=c(0,1),
    main="Adult nonbreeder to breeder", labels = labs,
    xlab = "Year", ylab= "Recruitment")
plt(object=out, params="mn.psiBA", ylim=c(0,1),
    main="Adult breeder to nonbreeder", labels = labs,
    xlab = "Year", ylab= "Recruitment")
plt(object=out, params="mn.pA", ylim=c(0,1),
    main="Nonbreeder", labels = labs,
    xlab = "Year", ylab= "Detection")
plt(object=out, params="mn.pB", ylim=c(0,1),
    main="Breeder", labels = labs,
    xlab = "Year", ylab= "Detection")
# plt(object=out, params="mn.f",
#     main="", labels = labs,
#     xlab = "Year", ylab= "Fecundity")
plt(object=out, params="mn.nest", ylim=c(0,1),
    main="", labels = labs,
    xlab = "Year", ylab= "Nest success")
# plt(object=out, params="mn.brood",
#     main="", labels = labs,
#     xlab = "Year", ylab= "Brood size")

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
MCMCsummary(post, params = pars1[1:2], 
            digits=2, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE)
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
MCMCsummary(post, params = pars1[3:4], 
            digits=3, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE)

# Estimates of population size
# NFY= first year, 
# NF = nonbreeders, 
# NB=breeders,
# Ntot=total
# All are presented as N[time, site], 
# where time 1=2011 ... 13=2023, 
# site 1 = LH and site 2=PC
MCMCsummary(post, params = pars1[5:8], 
            digits=2, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE)
MCMCsummary(post, params = pars1[9:16], 
            digits=2, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE)
MCMCsummary(post, params = "R", 
            digits=2, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE)
# Estimates of brood size, nest success, and fecundity
# indices are for sites, 1=LH and 2=PC
pars2 <- c("lmu.brood", "delta", "sig.brood",
           "mu.nest", "gamma") 
           #"mn.f", "mn.nest", "mn.brood" )
MCMCsummary(post, params = pars2, 
            digits=2, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE )

#### ---- traceplots ------
MCMCtrace(post, pdf=FALSE, params= "sds")
MCMCtrace(post, pdf=FALSE, params= "sds2")
MCMCtrace(post, pdf=FALSE, params= "mus")
MCMCtrace(post, pdf=FALSE, params= "betas")
MCMCtrace(post, pdf=FALSE, params= "NF")
MCMCtrace(post, pdf=FALSE, params= "NFY")
MCMCtrace(post, pdf=FALSE, params= "NB")
MCMCtrace(post, pdf=FALSE, params= "R")
# Brood size
## Mean brood size
MCMCtrace(post, pdf=FALSE, params= "lmu.brood")
## Effect from nest treatment on brood size
MCMCtrace(post, pdf=FALSE, params= "delta")
## SD among nests
MCMCtrace(post, pdf=FALSE, params= "sig.brood")
# Nest success
## Mean nest success
MCMCtrace(post, pdf=FALSE, params= "mu.nest")
## Effect from nest treatment on nest success
MCMCtrace(post, pdf=FALSE, params= "gamma")

#### ---- fit ------
# Goodness of fit check
plot.diag <- function(out, ratio=FALSE, 
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
plot.diag(out, ratio=F,
          name.rep="dmape.rep", 
          name.obs="dmape.obs",
          ind=1,
          lab="binomial", jit=300)
# nonbreeder, ind=2
plot.diag(out, ratio=F,
          name.rep="dmape.rep", 
          name.obs="dmape.obs",
          ind=2,
          lab="binomial", jit=300)
# first-year, ind=3
# poisson failed fit test bp=0
# Currently running models to try and fix
plot.diag(out, ratio=T,
          name.rep="dmape.rep", 
          name.obs="dmape.obs",
          ind=3,
          lab="Poisson", jit=300)

