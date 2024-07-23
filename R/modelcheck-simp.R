#### ---- setup -------
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_simp.Rdata")
load("data/data.rdata")
library ('MCMCvis')
library ('coda')

out <- list(as.mcmc(post[[1]]), 
            as.mcmc(post[[2]]), 
            as.mcmc(post[[3]]),
            as.mcmc(post[[4]]),
            as.mcmc(post[[5]]),
            as.mcmc(post[[6]]),
            as.mcmc(post[[7]]),
            as.mcmc(post[[8]]),
            as.mcmc(post[[9]]),
            as.mcmc(post[[10]]) )


outp <- MCMCpstr(out, type="chains")

#### ---- catplots -------
plt  <- function(object, params,...) {
  MCMCplot(object=out, 
           params=params, 
           guide_axis=TRUE, 
           HPD=TRUE, ci=c(80, 95), horiz=FALSE, 
           #ylim=c(-10,10),
           ...)
}

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


par(mfrow=c(1,1))
sds <- paste0("sds[", 1:9, "]")
plt(object=out, params=sds,
    exact=TRUE, ISB=FALSE,
    main="Temporal SDs (synchrony among sites)", 
    labels=c("FY survival", "NB survival", "B survival",
             "FY to B", "NB to B", "B to NB",
             "NB detection", "B detection",
             "Fecundity"))
sds2 <- paste0("sds2[", 1:9, "]")
plt(object=out, params=sds2,
    exact=TRUE, ISB=FALSE,
    main="Site-temporal SDs", 
    labels=c("FY survival", "NB survival", "B survival",
             "FY to B", "NB to B", "B to NB",
             "NB detection", "B detection",
             "Fecundity"))
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
# Both indices correspond to the parameters
# in the order listed above
# R[parameter1, parameter2], where R is the correlation between them
par(mfrow=c(1,1))
plt(object=out, params=Rs, exact=TRUE, ISB=FALSE, 
    main="Correlations btw demographic rates\n over time (synchrony)",
    xlab = "Rhos", guide_lines=TRUE)
# There's little temporal synchrony between sites
# R[1,3] appears slightly correlated and that suggests
# First year and breeder survival have some 
# synchrony. 

plt(object=out, params=R2s, exact=TRUE, ISB=FALSE, 
    main="Correlations btw demographic rates\n over time and sites",
    xlab = "Rhos", guide_lines=TRUE)
# Within each site 
# First year and breeder survival appear correlated (again R[1,3])
# as well as detection probability of breeders and fecundity (R[8,9])

# Annual averages for integration into the population model
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

#### ---- fecundity -------
# Averages
par(mfrow=c(1,1))
plt(object=out, 
    params=c("lmu.f"), 
    labels= c("Fecundity\n(log scale)\nLos Haitises",
              "Fecundity\n(log scale)\nPunta Cana"))

# gamma = nest treatment effect on fecundity
plt(object=out, 
    params=c("gamma"), 
    main="Anti-Parasitic Fly\nTreatment Effects", ylim=c(0,3))

plt(object=out, params="mn.f",
    main="", labels=labs,
    xlab = "Year", ylab= "Fecundity")
abline(v=13.5, lwd=2)

# So what's up with fecundity at LHNP? Let's examine the raw data.

# Are fewer nests being detected? No
# This would be consistent with lower detection prob of breeders.
nnests <- tapply(datl$f, list(constl$year.nest, constl$site.nest), length)
plot(2011:2023, nnests[,1], type="b",
     ylab="Number of nests detected", xlab="",
     main="Los Haitises")

# Are fewer fledglings being detected? Yes
ndet <- tapply(datl$f, list(constl$year.nest, constl$site.nest), sum)
plot(2011:2023, ndet[,1], type="b",
     ylab="Number of fledglings detected", xlab="",
     main="Los Haitises")

# Has the proportion of nests treated changed? Yes it has decreased.
ptreated <- tapply(constl$treat.nest, list(constl$year.nest, constl$site.nest), mean)
plot(2011:2023, ptreated[,1], type="b",
     ylab="Proportion of nests treated", xlab="",
     main="Los Haitises")

# Have the number of nests treated decreased? Maybe.
# The number of nests dropped slightly.
ntreated <- tapply(constl$treat.nest, list(constl$year.nest, constl$site.nest), sum)
plot(2011:2023, ntreated[,1], type="b",
     ylab="Number of nests treated", xlab="",
     main="Los Haitises")

# How do these line up with the number of breeders counted (Leah's dataset) at LH
# So the number of nests increased but the number of breeders decreased.
# That's confusing and doesn't make sense. Maybe that google sheet
# would be better to use given these discrepancies as you suggested. 
# No guarantees it will fix the average fecundity problem though.
# Perhaps it's density dependence kicking in? The number of breeders could
# be high enough that it's depressing fecundity.
plot(2011:2023, datl$counts[3,,1], type="b",
     ylab="Number of breeders detected", xlab="",
     main="Los Haitises")

# In summary the number of nest treatments remains relatively stable, 
# but it appears that the proportion of nests treated has
# declined because of an increase in the number of nests detected
# perhaps depressing the average affect on fecundity.
# Could it be that there are so many nests that the treatment effects
# are small at this large population level?
# Or density dependence may be depressing fecundity. 
# Count data of breeders do not align.



#### ---- tables -------
pars1 <- c("sds", "sds2","mus", "betas",
           "mn.phiFY","mn.phiA", "mn.phiB", 
           "mn.psiFYB", "mn.psiAB", "mn.psiBA", 
           "mn.pA", "mn.pB")
# Estimates of parameters in the survival model 
# In this order: 
# (1) FY survival 
# (2) NB survival 
# (3) B survival 
# (4) FY to B recruitment
# (5) NB to B recruitment
# (6) B to NB recruitment
# (7) Detection NB 
# (8) Detection B

# Temporal synchrony
# sds indices correspond to the order listed above. 
# This error term is for temporal variation among
# both sites during each year
# sds are on the logit scale
MCMCsummary(post, params = pars1[1], 
            digits=2, HPD = T,
            hpd_prob = 0.80, pg0= TRUE, 
            func = median,
            func_name = "median")
# Temporal by site interaction
# these sds allow sites to vary over time independently
# sds2 are on the logit scale
MCMCsummary(post, params = pars1[2], 
            digits=2, HPD = T,
            hpd_prob = 0.80, pg0= TRUE, 
            func = median,
            func_name = "median")
# Mus are overall means for survival, recruitment, and detection
# mus[parameter, site] , where site=1 is LH and site=2 is PC
# Again the order is listed above
# "mus" are on the probability scale, range[0,1]
# modeled as logit(probability) = lmus[x,site] + betas[x]*translocated + eps[x,t] + eta[x,s,t]
# except breeder to nonbreeder recruitment = lmus[6,site] 
MCMCsummary(post, params = pars1[3], 
            digits=3, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE, 
            func = median,
            func_name = "median")
# Same means but on the logit scale to match sds and sds2 
MCMCsummary(post, params = "lmus", 
            digits=3, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE, 
            func = median,
            func_name = "median")
# fecundity estimates correspond to sds[9] and sds2[9]
# these are on the log scale. Don't forget to exp(lmu.f).
# modeled as log(f) = lmu.f[site] + gamma*treatment + eps[x,t] + eta[x,s,t]
MCMCsummary(post, params = "lmu.f", 
            digits=2, HPD = T,
            hpd_prob = 0.80, pg0= TRUE, 
            func = median,
            func_name = "median")

# Effects for each parameter from
# translocation and hacking inidividuals
MCMCsummary(post, params = pars1[4], 
            digits=3, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE, 
            func = median,
            func_name = "median")

# Annual estimates of 
# (1) FY survival 
MCMCsummary(post, params = pars1[5], 
            digits=2, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE, 
            func = median,
            func_name = "median")
# (2) NB survival
MCMCsummary(post, params = pars1[6], 
            digits=2, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE, 
            func = median,
            func_name = "median")
# (3) B survival
MCMCsummary(post, params = pars1[7], 
            digits=2, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE, 
            func = median,
            func_name = "median")
# (4) FY to B recruitment
MCMCsummary(post, params = pars1[8], 
            digits=2, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE, 
            func = median,
            func_name = "median")
# (5) NB to B recruitment
MCMCsummary(post, params = pars1[9], 
            digits=2, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE, 
            func = median,
            func_name = "median")
# (6) B to NB recruitment
MCMCsummary(post, params = pars1[10], 
            digits=2, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE, 
            func = median,
            func_name = "median")
# (7) Detection NB 
MCMCsummary(post, params = pars1[11], 
            digits=2, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE, 
            func = median,
            func_name = "median")
# (8) Detection B
MCMCsummary(post, params = pars1[12], 
            digits=2, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE, 
            func = median,
            func_name = "median")

# Set up an index for temporal and site sds
ind <- 1
Rs <- R2s <- c()
for (i in 1:(nrow(outp$R)-1)){
  for (j in (i+1):nrow(outp$R)){
    Rs[ind] <- paste0("R[",i,", ", j, "]")
    R2s[ind] <- paste0("R2[",i,", ", j, "]")
    ind <- ind+1
  }}
# Correlations among demographic rates time (synchrony)
MCMCsummary(post, params = Rs, 
            digits=2, HPD = T,
            exact=TRUE, ISB=FALSE,
            hpd_prob = 0.80, pg0= TRUE)
# Correlations among demographic rates site x time
MCMCsummary(post, params = R2s,
            exact=TRUE, ISB=FALSE,
            digits=2, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE)

#### ---- traceplots ------
MCMCtrace(post, pdf=FALSE, params= "sds")
MCMCtrace(post, pdf=FALSE, params= "sds2")
MCMCtrace(post, pdf=FALSE, params= "mus")
MCMCtrace(post, pdf=FALSE, params= "betas")
MCMCtrace(post, pdf=FALSE, params= Rs, exact=TRUE, ISB=FALSE)
MCMCtrace(post, pdf=FALSE, params= R2s, exact=TRUE, ISB=FALSE)

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

# check goodness-of-fit for fecundity
fit.check(out, ratio=F,
          name.rep="f.dmape.rep", 
          name.obs="f.dmape.obs",
          ind=1,
          lab="Fecundity-Neg binomial", jit=300)

