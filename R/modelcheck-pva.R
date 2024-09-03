#### ---- setup -------
#load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_sites.rdata")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\pva.rdata")
load("data/data.rdata")
library ('MCMCvis')
library ('coda')
library ('ggplot2')
library('reshape2')
out <- list(as.mcmc(post[[1]]), 
             as.mcmc(post[[2]]), 
             as.mcmc(post[[3]]),
             as.mcmc(post[[4]]),
             as.mcmc(post[[5]]),
             as.mcmc(post[[6]]),
             as.mcmc(post[[7]]),
             as.mcmc(post[[8]]),
             as.mcmc(post[[9]]),
             as.mcmc(post[[10]]))

# Identify chains with NAs that 
# failed to initialize
NAlist <- c()
for (i in 1:length(out)){
  NAlist[i] <- any (is.na(out[[i]][,1:9493]) | out[[i]][,1:9493]<0)
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
library (viridis)
cols <- viridis(6)
pe <- apply(outp$extinct, c(1,2,3), mean)
par(mfrow=c(1,1))
plot(2024:2123, pe[1, ,1], 
     type="l", col=cols[1], lwd=2,
     ylim=c(0,1), main="Los Haitises",
     ylab="Quasi-extinction probability", xlab="Year")
lines(2024:2123, pe[2, ,1], col=cols[2], lwd=2)
lines(2024:2123, pe[3, ,1], col=cols[3], lwd=2)
lines(2024:2123, pe[4, ,1], col=cols[4], lwd=2, lty=2)
lines(2024:2123, pe[5, ,1], col=cols[5], lwd=2, lty=2)
lines(2024:2123, pe[6, ,1], col=cols[6], lwd=2, lty=2)
legend(x=2024, y=1, title="Scenarios",
        legend=c("hacked=0, nests=all", "hacked=0, nests=0", "hacked=0, nests=10",
                "hacked=10, nests=all", "hacked=10, nests=0", "hacked=10, nests=10"),
        xpd=NA, col=cols, lty=c(1,1,1,2,2,2), lwd=2)

par(mfrow=c(1,1))
plot(2024:2123, pe[1, ,2], 
     type="l", col=cols[1], lwd=2,
     ylim=c(0,1), main="Punta Cana",
     ylab="Quasi-extinction probability", xlab="Year")
lines(2024:2123, pe[2, ,2], col=cols[2], lwd=2)
lines(2024:2123, pe[3, ,2], col=cols[3], lwd=2)
lines(2024:2123, pe[4, ,2], col=cols[4], lwd=2, lty=2)
lines(2024:2123, pe[5, ,2], col=cols[5], lwd=2, lty=3)
lines(2024:2123, pe[6, ,2], col=cols[6], lwd=2, lty=4)
legend(x=2024, y=1, title="Scenarios",
       legend=c("hacked=0, nests=all", "hacked=0, nests=0", "hacked=0, nests=10",
                "hacked=10, nests=all", "hacked=10, nests=0", "hacked=10, nests=10"),
       xpd=NA, col=cols, lty=c(1,1,1,2,2,2), lwd=2)

plot(2024:2123, pe[3, ,1], type="l", 
     ylim=c(0,1),
     ylab="Prob. of extinction", xlab="Year", main="Scenario 3")
plot(2024:2123, pe[4, ,1], type="l", 
     ylim=c(0,1),
     ylab="Prob. of extinction", xlab="Year", main="Scenario 4")
plot(2024:2123, pe[5, ,1], type="l", 
     ylim=c(0,1),
     ylab="Prob. of extinction", xlab="Year", main="Scenario 5")
plot(2024:2123, pe[6, ,1], type="l", 
     ylim=c(0,1),
     ylab="Prob. of extinction", xlab="Year", main="Scenario 6")

plot(2024:2123, pe[1, ,2], type="l", 
     ylim=c(0,1),
     ylab="Prob. of extinction", xlab="Year", main="Scenario 1")
plot(2024:2123, pe[2, ,2], type="l", 
     ylim=c(0,1),
     ylab="Prob. of extinction", xlab="Year", main="Scenario 2")
plot(2024:2123, pe[3, ,2], type="l", 
     ylim=c(0,1),
     ylab="Prob. of extinction", xlab="Year", main="Scenario 3")
plot(2024:2123, pe[4, ,2], type="l", 
     ylim=c(0,1),
     ylab="Prob. of extinction", xlab="Year", main="Scenario 4")
plot(2024:2123, pe[5, ,2], type="l", 
     ylim=c(0,1),
     ylab="Prob. of extinction", xlab="Year", main="Scenario 5")
plot(2024:2123, pe[6, ,2], type="l", 
     ylim=c(0,1),
     ylab="Prob. of extinction", xlab="Year", main="Scenario 6")


par(mfrow=c(3,2))
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("extinct[1, ",1:100, ", 1]"), 
    main="All stages\n Los Haitises", 
    labels = 2024:2123,
    xlab = "Year", ylab= "Prob. of extinction",
    ylim=c(0,1))

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("extinct[2, ",1:100, ", 1]"), 
    main="All stages\n Los Haitises", 
    labels = 2024:2123,
    xlab = "Year", ylab= "Prob. of extinction",
    ylim=c(0,1))

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("extinct[3, ",1:100, ", 1]"), 
    main="All stages\n Los Haitises", 
    labels = 2024:2123,
    xlab = "Year", ylab= "Prob. of extinction",
    ylim=c(0,1))

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("extinct[4, ",1:100, ", 1]"), 
    main="All stages\n Los Haitises", 
    labels = 2024:2123,
    xlab = "Year", ylab= "Prob. of extinction",
    ylim=c(0,1))

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("extinct[5, ",1:100, ", 1]"), 
    main="All stages\n Los Haitises", 
    labels = 2024:2123,
    xlab = "Year", ylab= "Prob. of extinction",
    ylim=c(0,1))

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("extinct[6, ",1:100, ", 1]"), 
    main="All stages\n Los Haitises", 
    labels = 2024:2123,
    xlab = "Year", ylab= "Prob. of extinction",
    ylim=c(0,1))

##########################3
par(mfrow=c(3,2))
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("Ntot[1, ",1:113, ", 1]"), 
    main="All stages\n Los Haitises", 
    labels = 2011:2123,
    xlab = "Year", ylab= "Abundance")

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("Ntot[2, ",1:113, ", 1]"), 
    main="All stages\n Los Haitises", 
    labels = 2011:2123,
    xlab = "Year", ylab= "Abundance")

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("Ntot[3, ",1:113, ", 1]"), 
    main="All stages\n Los Haitises", 
    labels = 2011:2123,
    xlab = "Year", ylab= "Abundance")

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("Ntot[4, ",1:113, ", 1]"), 
    main="All stages\n Los Haitises", 
    labels = 2011:2123,
    xlab = "Year", ylab= "Abundance")

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("Ntot[5, ",1:113, ", 1]"), 
    main="All stages\n Los Haitises", 
    labels = 2011:2123,
    xlab = "Year", ylab= "Abundance")

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("Ntot[6, ",1:113, ", 1]"), 
    main="All stages\n Los Haitises", 
    labels = 2011:2123,
    xlab = "Year", ylab= "Abundance")

par(mfrow=c(3,2))
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("Ntot[1, ",1:113, ", 2]"), 
    main="All stages\n Punta Cana", 
    labels = 2011:2123,
    xlab = "Year", ylab= "Abundance")

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("Ntot[2, ",1:113, ", 2]"), 
    main="All stages\n Punta Cana", 
    labels = 2011:2123,
    xlab = "Year", ylab= "Abundance")

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("Ntot[3, ",1:113, ", 2]"), 
    main="All stages\n Punta Cana", 
    labels = 2011:2123,
    xlab = "Year", ylab= "Abundance")

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("Ntot[4, ",1:113, ", 2]"), 
    main="All stages\n Punta Cana", 
    labels = 2011:2123,
    xlab = "Year", ylab= "Abundance")

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("Ntot[5, ",1:113, ", 2]"), 
    main="All stages\n Punta Cana", 
    labels = 2011:2123,
    xlab = "Year", ylab= "Abundance")

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("Ntot[6, ",1:113, ", 2]"), 
    main="All stages\n Punta Cana", 
    labels = 2011:2123,
    xlab = "Year", ylab= "Abundance")
# plt(object=out, 
#     exact=TRUE, ISB=FALSE, 
#     params=paste0("NB[1, ",1:113, ", 1]"), 
#     main="Adult Breeder (NB)\n Los Haitises", 
#     labels = 2011:2123,
#     xlab = "Year", ylab= "Abundance")

# Abundance of females at Punta Cana
par(mfrow=c(5,2))
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("NFY[",1:13, ", 2]"), 
    main="First-year (FY)\n Punta Cana", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plot(2011:2023, datl$counts[1,,2]+constl$hacked.counts[,2], 
     ylab="Counts", type="b")

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
     ylab="Counts", type="b")

plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("Ntot[",1:13, ", 2]"), 
    main="All stages\n Punta Cana", 
    labels = 2011:2023,
    xlab = "Year", ylab= "Abundance")
plot(2011:2023, colSums(datl$counts[,,2]), 
     ylab="Counts", type="b")

#### ---- catplots2 -------
# Finer population segments
par(mfrow=c(4,2))
plt(object=out, 
    exact=TRUE, ISB=FALSE, 
    params=paste0("N[", 1, ", ", 1:13, ", 1]"), 
    main="Los Haitises\nFirst-years born", 
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
    main="Punta Cana\nFY born",
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
par(mfrow=c(1,1))
plt(object=out, params=Rs, exact=TRUE, ISB=FALSE, 
    main="Correlations btw demographic rates\n over time (synchrony)",
    xlab = "Rhos", guide_lines=TRUE)
plt(object=out, params=R2s, exact=TRUE, ISB=FALSE, 
    main="Correlations btw demographic rates\n over time and sites",
    xlab = "Rhos", guide_lines=TRUE)
# lmu.brood = mean brood size (log scale), 
# sig.brood = SD among nests
# mu.nest = mean nest success
par(mfrow=c(1,1))
plt(object=out, 
    params=c("lmu.f"), 
    labels= c("Fecundity\n(log scale)\nLos Haitises",
              "Fecundity\n(log scale)\nPunta Cana"))

# gamma = nest treatment effect on fecundity
plt(object=out, 
    params=c("gamma"), 
    main="Anti-Parasitic Fly\nTreatment Effects", ylim=c(0,3))

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
plt(object=out, params="mn.f",
    main="", labels=labs,
    xlab = "Year", ylab= "Fecundity")
abline(v=13.5, lwd=2)

#### ----popstructure -----

mdFY <-  apply(outp$NFY, c(1,2), median) 
mdB <-  apply(outp$NB, c(1,2), median) 
mdF <-  apply(outp$NF, c(1,2), median) 
lFY <- melt(mdFY)
lB <- melt(mdB)
lF <- melt(mdF)
lFY$Stage <- "First-year"
lB$Stage <- "Breeder"
lF$Stage <- "Nonbreeder"
ldat <- rbind(lFY, lB, lF)
colnames(ldat)[1:3] <- c("Year", "Sitenum", "Number") 
ldat$Site <- ifelse(ldat$Sitenum==1, "Los Haitises", "Punta Cana")

# Use median number of females in each stage
# to plot an approximate population structure
ggplot(ldat, aes(fill=Stage, y=as.numeric(Number), x=Year)) + 
  geom_bar(position="fill", stat="identity") +
  ylab("Proportion of population") + 
  facet_wrap("Site")

ggplot(ldat, aes(fill=Stage, y=as.numeric(Number), x=Year)) + 
  geom_bar(position="stack", stat="identity") +
  ylab("Numer of females") + 
  facet_wrap("Site", scales = "free_y")

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

MCMCsummary(post2, params = c(sds, sds2), 
            exact=TRUE, ISB=FALSE,
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
# modeled as logit(probability) = lmus[x,site] + betas[x]*translocated + eps[x,t] + eta[x,s,t]
# except breeder to nonbreeder recruitment = lmus[6,site] 
MCMCsummary(post2, params = pars1[3], 
            digits=3, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE)

MCMCsummary(post2, params = "lmus", 
            digits=3, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE)

MCMCsummary(post2, params = pars1[4], 
            digits=3, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE)

# Fecundity
# modeled as log(f) = lmu.f[site] + gamma*treatment + eps[x,t] + eta[x,s,t]
# log scale
MCMCsummary(post2, params = "lmu.f", 
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
MCMCsummary(post2, params = pars1[5:8], 
            digits=2, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE)

MCMCsummary(post2, params = pars1[9:16], 
            digits=2, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE)
# Correlations among demographic rates time (synchrony)
MCMCsummary(post2, params = "R", 
            digits=2, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE)
# Correlations among demographic rates site x time
MCMCsummary(post2, params = "R2", 
            digits=2, HPD = T, 
            hpd_prob = 0.80, pg0= TRUE)

#### ---- traceplots ------
MCMCtrace(post2, pdf=FALSE, params= "sds")
MCMCtrace(post2, pdf=FALSE, params= "sds2")
MCMCtrace(post2, pdf=FALSE, params= "mus")
MCMCtrace(post2, pdf=FALSE, params= "betas")
MCMCtrace(post2, pdf=FALSE, params= "NF")
MCMCtrace(post2, pdf=FALSE, params= "NFY")
MCMCtrace(post2, pdf=FALSE, params= "NB")
MCMCtrace(post2, pdf=FALSE, params= "R")

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
          lab="Fecundity-Neg binomial", jit=300)
