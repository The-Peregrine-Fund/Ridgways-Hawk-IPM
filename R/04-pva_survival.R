## ---- pva --------
library('nimble')
library('parallel')
library ('coda')

load("/bsuscratch/brianrolek/riha_ipm/data.rdata")
source("/bsuscratch/brianrolek/riha_ipm/MCMCvis.R")
load("/bsuscratch/brianrolek/riha_ipm/outputs/ipm_shortrun.rdata")

# load("C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//outputs//ipm_longrun.rdata")
# load("data//data.rdata")
# library ("MCMCvis")

#**********************
#* Set data to reflect PVA scenarios
#**********************
cpus <- 5 # number of processors

constl$K <- 50 # number of future years
constl$SC <- 45 # number of PVA scenarios
datl$constraint_data <- rbind(datl$constraint_data, array(1, dim=c(constl$K,2)) ) # to help constrain FYs born + hacked to be positive
constl$effort2 <- rbind(constl$effort2, array(0, dim=c(constl$K,2), dimnames=list(2024:(2024+constl$K-1), c("LH", "PC"))))

constl$num.treated <- rep( c(0, 15, 30, 45, 100, 0, 15, 30, 45, 100, 0, 15, 30, 45, 100), each=3 ) # 100s sub in for All but are over-ridden in model code
num.hacked <- rep( c(0, 0, 0, 0, 0, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10), each=3 )
surv.diff <- array(NA, dim=c(constl$SC, 3, 2), dimnames=list(scenario=1:constl$SC, stage=c('FY', 'NB', 'B'), site=c('LH', 'PC')))
surv.diff[,1,] <- matrix(c(rep( c(0, 0.130, 0.260), 15),
                           rep( c(0, 0.146, 0.292), 15)), ncol=2)
surv.diff[,2,] <- matrix(c(rep( c(0, 0.006, 0.012), 15),
                           rep( c(0, 0.022, 0.044), 15)), ncol=2)
surv.diff[,3,] <- matrix(c(rep( c(0, 0.016, 0.032), 15),
                           rep( c(0, 0.026, 0.052), 15)), ncol=2)
constl$surv.diff <- surv.diff 

hacked.counts <- array(0, dim=c(constl$SC, constl$nyr+constl$K, 2))
for (sc in 1:constl$SC){
  hacked.counts[sc,1:constl$nyr, 1:2] <- constl$hacked.counts[, 1:2] 
  hacked.counts[sc,14:(constl$nyr+constl$K), 1] <- -num.hacked[sc]
  hacked.counts[sc,14:(constl$nyr+constl$K), 2] <-  num.hacked[sc]
}
constl$hacked.counts <- hacked.counts

#**********************
#* Specify priors taken from posterior of IPM
#**********************
#* This helps ensure that all chains run
# Identify chains with NAs that failed to initialize
out <- lapply(post, as.mcmc)
NAlist <- c()
for (i in 1:length(out)){
  NAlist[i] <- any (is.na(out[[i]][,1:286]) | out[[i]][,1:286]<0)
}
out <- out[!NAlist]
outp <- MCMCpstr(out, type="chains")
!NAlist

Ni.func <- function (){
  # randomly select inits from ones that work
  # Ni <- outp$N[1:7,1:constl$nyr,1:2,
  #              sample(seq(1, 10000, by=100), 1, replace = F)]
  # take the means of the posterior for inits
  Ni <- apply(outp$N, c(1,2,3), mean) |> round()
  Ni.pva <- array(NA, dim=c(constl$SC, dim(Ni)+c(0,constl$K,0)))
  for (sc in 1:constl$SC){
    Ni.pva[sc, 1:7, 1:constl$nyr, 1:2] <- Ni
    for (t in constl$nyr:(constl$nyr+constl$K)){
      for (s in 1:2){
        Ni.pva[sc, 1:7, t, s] <- Ni[1:7, 13, s] 
      }}} # t sc
  return(Ni.pva)
} # function

u2 <- apply(outp$Ustar, c(1,2), mean)

inits.func.pva <- function (){
  list(  
    # fecundity 
    lmu.prod = apply(outp$lmu.prod, 1, mean),
    gamma = mean(outp$gamma), 
    rr = mean(outp$rr),
    # survival
    z = par_info[[1]]$inits$z, 
    mus = apply(outp$mus, c(1,2), mean), 
    betas = apply(outp$betas, 1, mean),
    deltas = apply(outp$deltas, 1, mean),
    sds =  apply(outp$sds, 1, mean),
    Ustar = u2,
    # counts
    countsAdults= matrix(c(374, 335, 305, 295, rep(NA, length(2015:2023)), rep(NA, length(2011:2023)) ), nrow=13), 
    r = mean(outp$r),
    N = Ni.func()
  )}

# Set seed then save for reproducibility
# then randomly draw for each chain
set.seed(1)
seeds <- sample(1:1000000, size=cpus, replace=FALSE)
par_info_pva <- list()
for (i in 1:cpus){
  par_info_pva[[i]] <- list(seed=seeds[i], inits = inits.func.pva())
}

#**********************
#* Parameter descriptions
#**********************
# Mark-resight-recovery data
#   Observations (po) = y  
#     1 seen first-year (age=0, just before 1st b-day)
#     2 seen nonbreeder
#     3 seen breeder
#     4 not seen
#   States (ps)
#     1 alive first-year
#     2 alive nonbreeder
#     3 alive breeder
#     4 dead
#   Groups
#     1 wild-born
#     2 translocated and hacked
###################################################
# PARAMETERS
#   phiFY: survival probability first year 
#   phiA: survival probability nonbreeders
#   phiB: survival probability breeders
#   psiFYB: recruitment probability from first-year to breeder
#   psiAB: recruitment probability from nonbreeders to breeder
#   psiBA: recruitment probability from breeder to nonbreeders 
#   pFY: resight probability first-years
#   pA: resight probability nonbreeders
#   pB: resight probability breeders
#   lmu.prod: mean productivity (males and females) per territory on the log scale
#   sds: standard deviations for multivariate normal random effects for time and site
#   R: correlation coefficients for multivariate normal random effects for time and site
#   lambda: population growth rate (derived)
#   extinct: binary indicator of extirpation at a site
#   gamma: coefficient of effect from nest treatments
#   betas: coefficient of effect from translocations
#   deltas: coefficient of effect from survey effort
#   mus: overall means for survival, recruitment, and detections
#   r and rr: "r" parameter for negative binomial distribution
#             also described as omega in manuscript

#**********************
#* Model code
#**********************
mycode <- nimbleCode(
  {
    ###################################################
    # Priors and constraints
    ###################################################
    # survival, recruitment, and detection can be correlated
    for (k in 1:8){ 
      betas[k] ~ dnorm(0, sd=20)  # prior for coefficients
      deltas[k] ~ dnorm(0, sd=10)
    } # k
    
    for (j in 1:8){
      for (s in 1:nsite){
        lmus[j,s] <- logit(mus[j,s])
        mus[j,s] ~ dbeta(1,1) # prior for means
      }}     # m population #s sex #h hacked
    
    # Temporal random effects and correlations between sites
    # Non-centered parameterization of the multivariate normal distribution to improve convergence
    for (jj in 1:p){ sds[jj] ~ dexp(1) }# prior for temporal variation
    # estimated using the multivariate normal distribution
    R[1:p,1:p] <- t(Ustar[1:p,1:p]) %*% Ustar[1:p,1:p] # calculate rhos, correlation coefficients
    Ustar[1:p,1:p] ~ dlkj_corr_cholesky(eta=1.1, p=p) # Ustar is the Cholesky decomposition of the correlation matrix
    # multivariate normal for temporal variance
    for (t in 1:(nyr+K)){ # survival params only have nyr-1, no problem to simulate from however
      for (s in 1:nsite){
        eta[1:p,s,t] <- diag(sds[1:p]) %*% t(Ustar[1:p,1:p]) %*% z.score[1:p,s,t]
        for(j in 1:p){
          z.score[j,s,t] ~ dnorm(0, sd=1)  # z-scores
        } # j
      } } # s t 
    
    #######################
    # Derived params
    #######################
    for(sc in 1:SC){
      for (s in 1:nsite){
        for (t in 1:(nyr-1+K)){
          lambda[sc, t, s] <-  Ntot[sc, t+1, s]/(Ntot[sc, t, s]) # population groewth rate
          loglambda[sc, t, s] <- log(lambda[sc, t, s]) # r
        }} #t
      
      for (s in 1:nsite){
        for (t in 1:K){
          # quasi-extinction probabilities given population segments, N total, adults, and breeders
          extinct[sc, t, s] <- equals(Ntot[sc, nyr+t, s], 0) 
          extinctAD[sc, t, s] <- equals(NAD[sc, nyr+t, s], 0)
          extinctB[sc, t, s] <- equals(NB[sc, nyr+t, s], 0)
        }}
    } # sc
    
    ###############################
    # Likelihood for productivity
    ###############################
    # Priors 
    for (s in 1:nsite){
      lmu.prod[s] ~ dnorm(0, sd=5)
    } # s
    gamma ~ dnorm(0, sd=10)
    rr ~ dexp(0.05)
    
    # Productivity likelihood       
    for (k in 1:npairsobs){
      prod[k] ~ dnegbin(ppp[k], rr)
      ppp[k] <- rr/(rr+mu.prod[k])
      log(mu.prod[k]) <- lmu.prod[site.pair[k]] +  
        gamma*treat.pair[k] + 
        eta[9, site.pair[k], year.pair[k] ] 
    } # k
    # Derive yearly productivity for population model
    # need to reorder because nimble doesn't 
    # handle nonconsecutive indices
    # yrind.pair is a matrix of indices for each site
    for (t in 1:nyr){
      for (s in 1:nsite){
        for (xxx in 1:pair.end[t,s]){
          prodmat[t,s,xxx] <- mu.prod[ yrind.pair[xxx,t,s] ]
        } # xxx
        mn.prod[1,t,s] <- mean( prodmat[t,s,1:pair.end[t,s]] )
        for(sc in 2:SC){
          mn.prod[sc, t, s] <- mn.prod[1,t,s]
        }}} # s t
    # Future projections- fecundity
    for(sc in 1:SC){
      for (t in (nyr+1):(nyr+K) ){
        for (s in 1:nsite){
          numer.perc.treat[sc, t, s] <-  # numerator
            step( NB[sc, t, s]-(num.treated[sc]+1) ) * num.treated[sc]  + # num.treated > NB, use num.treated
            ( 1-step( NB[sc, t, s]-(num.treated[sc]+1) )) * NB[sc, t, s]   #  num.treated < NB set to NB
          denom.perc.treat[sc, t, s] <- # denominator, total nests available for treatment
            step( NB[sc, t, s]-(num.treated[sc]+1) ) * NB[sc, t, s]  + # num.treated > NB, use num.treated
            ( 1-step( NB[sc, t, s]-(num.treated[sc]+1) )) * numer.perc.treat[sc, t, s]
          denom.perc.treat2[sc, t, s] <- equals(denom.perc.treat[sc, t, s], 0) +
            (1- equals(denom.perc.treat[sc, t, s], 0)) *denom.perc.treat[sc, t, s]
          perc.treat[sc, t, s] <- numer.perc.treat[sc, t, s] /
            denom.perc.treat2[sc, t, s] 
          
          log(mn.prod[sc, t, s]) <- lmu.prod[s] + 
            gamma*perc.treat[sc, t, s] + # treat.pair set to one here.
            eta[9, s, t] 
        }} } # s t sc
    
    ################################
    # Likelihood for counts
    ################################
    # Abundance for year=1
    for (v in 1:7){ 
      for (s in 1:nsite){
        # Subtract one to allow dcat to include zero
        N[1, v, 1, s] <- N2[1, v, 1, s] - 1 
        N2[1, v, 1, s] ~ dcat(pPrior[v, 1:s.end[v,s], s])
        # Assign estimates for years with data to all scenarios
        for (sc in 2:SC){
          N[sc, v, 1, s] <- N[1, v, 1, s]
        }
      }} # s t
    
    # Abundance for years > 1
    for (t in 1:(nyr-1)){
      for (s in 1:nsite){
        # Number of wild born juvs
        N[1, 1, t+1, s] ~ dpois( (NFY[1, t, s]*mn.phiFY[1, t, s]*mn.psiFYB[1, t, s] + # first year breeders
                                    NF[1, t, s]*mn.phiA[1, t, s]*mn.psiAB[1, t, s] + # nonbreeders to breeders
                                    NB[1, t, s]*mn.phiB[1, t, s]*(1-mn.psiBA[1, t, s])) # breeders remaining
                                 *(mn.prod[1, t+1, s]/2) ) # end Poisson
        # Abundance of nonbreeders
        N[1, 2, t+1, s] ~ dbin(mn.phiFY[1, t, s]*(1-mn.psiFYB[1, t, s]), NFY[1, t, s]) # Nestlings to second year nonbreeders
        N[1, 3, t+1, s] ~ dbin(mn.phiA[1, t, s]*(1-mn.psiAB[1, t, s]), NF[1, t, s]) # Nonbreeders to nonbreeders
        N[1, 4, t+1, s] ~ dbin(mn.phiB[1, t, s]*mn.psiBA[1, t, s], NB[1, t, s]) # Breeders to nonbreeders
        # Abundance of breeders
        N[1, 5, t+1, s] ~ dbin(mn.phiFY[1, t, s]*mn.psiFYB[1, t, s], NFY[1, t, s]) # Nestlings to second year breeders
        N[1, 6, t+1, s] ~ dbin(mn.phiA[1, t, s]*mn.psiAB[1, t, s], NF[1, t, s]) # Nonbreeder to breeder
        N[1, 7, t+1, s] ~ dbin(mn.phiB[1, t, s]*(1-mn.psiBA[1, t, s]), NB[1, t, s]) # Breeder to breeder
        # Assign estimates for years with data to 
        # all scenarios
        for (v in 1:7){ 
          for (sc in 2:SC){
            N[sc, v, t+1, s] <- N[1, v, t+1, s]
          } } # sc
      }} # s t
    
    # Future projections- population model 
    for (sc in 1:SC){
      for (t in nyr:(nyr+K-1)){
        for (s in 1:nsite){
          # Number of wild born juvs
          N[sc, 1, t+1, s] ~ dpois( (NFY[sc, t, s]*mn.phiFY[sc, t, s]*mn.psiFYB[sc, t, s] + # first year breeders
                                       NF[sc, t, s]*mn.phiA[sc, t, s]*mn.psiAB[sc, t, s] + # nonbreeders to breeders
                                       NB[sc, t, s]*mn.phiB[sc, t, s]*(1-mn.psiBA[sc, t, s])) # breeders remaining
                                    *(mn.prod[sc, t+1, s]/2) ) # end Poisson
          # Abundance of nonbreeders
          ## Second year nonbreeders
          N[sc, 2, t+1, s] ~ dbin(mn.phiFY[sc, t, s]*(1-mn.psiFYB[sc, t, s]), NFY[sc, t, s]) # Nestlings to second year nonbreeders
          ## Adult nonbreeders
          N[sc, 3, t+1, s] ~ dbin(mn.phiA[sc, t, s]*(1-mn.psiAB[sc, t, s]), NF[sc, t, s]) # Nonbreeders to nonbreeders
          N[sc, 4, t+1, s] ~ dbin(mn.phiB[sc, t, s]*mn.psiBA[sc, t, s], NB[sc, t, s]) # Breeders to nonbreeders
          # Abundance of breeders
          ## Second year breeders
          N[sc, 5, t+1, s] ~ dbin(mn.phiFY[sc, t, s]*mn.psiFYB[sc, t, s], NFY[sc, t, s]) # Nestlings to second year breeders
          ## Adult breeders
          N[sc, 6, t+1, s] ~ dbin(mn.phiA[sc, t, s]*mn.psiAB[sc, t, s], NF[sc, t, s]) # Nonbreeder to breeder
          N[sc, 7, t+1, s] ~ dbin(mn.phiB[sc, t, s]*(1-mn.psiBA[sc, t, s]), NB[sc, t, s]) # Breeder to breeder
        }} # s t
    } # sc
    
    # Count likelihoods, state-space model, and observation process 
    for (t in 1:nyr){
      for (s in 1:nsite){
        num.hacked[1, t, s] <- step( N[1, 1, t, 1] + hacked.counts[1, t, 1] ) * hacked.counts[1, t, s]  +
          (1- step( N[1, 1, t, 1] + hacked.counts[1, t, 1] )) * N[1, 1, t, 1] * 
          (-1*equals(s,1) + equals(s,2)) # change to 
        constraint_data[t, s] ~ dconstraint( (N[1, 1, t, s] + hacked.counts[1, t, s]) >= 0 ) # Transfers translocated first-year females
        NFY[1, t, s] <- N[1, 1, t, s] + hacked.counts[1, t, s] # Transfers translocated first-year females
        NF[1, t, s] <- sum(N[1, 2:4, t, s]) # number of adult nonbreeder females
        NB[1, t, s] <- sum(N[1, 5:7, t, s]) # number of adult breeder females
        NAD[1, t, s] <- NF[1, t, s] + NB[1, t, s] # number of adults
        Ntot[1, t, s] <- sum(N[1, 1:7, t, s]) # total number of females
        countsAdults[t, s] ~ dpois(lamAD[1, t, s]) # adult females
        # constrain N1+hacked.counts to be >=0
        log(lamAD[1, t, s]) <- log(NAD[1, t, s]) + deltas[1]*effort2[t, s] + deltas[2]*effort2[t, s]^2
        log(lamFY[1, t, s]) <- log(N[1, 1, t, s]) + deltas[3]*effort2[t, s] + deltas[4]*effort2[t, s]^2
        for (sc in 2:SC){
          lamAD[sc, t, s] <- lamAD[1, t, s]
          lamFY[sc, t, s] <- lamFY[1, t, s]
          NFY[sc, t, s] <- NFY[1, t, s]
          NF[sc, t, s] <- NF[1, t, s]
          NB[sc, t, s] <- NB[1, t, s]
          NAD[sc, t, s] <- NAD[1, t, s]
          Ntot[sc, t, s] <- Ntot[1, t, s]
          num.hacked[sc, t, s] <- num.hacked[1, t, s]
        } # sc
      }# s
      # First-years at different sites have different distributions
      # for better model fit
      countsFY[t, 1] ~ dpois(lamFY[1, t, 1]) # doesn't have any zeroes so poisson
      countsFY[t, 2] ~ dnegbin(pp[t], r) # first year females, includes translocated/hacked
      pp[t] <- r/(r+(lamFY[1, t, 2] ))
    } # t
    r ~ dexp(0.05)
    
    # Future projections of counts 
    for(sc in 1:SC){
      for (t in (nyr+1):(nyr+K)){
        for (s in 1:nsite){
          # constrain N1+hacked.counts to be >=0, and allow it to shrink as the population shrinks
          # If LH has fewer birds than translocations, none are translocated
          num.hacked[sc, t, s] <- step( N[sc, 1, t, 1] + hacked.counts[sc, t, 1] ) * hacked.counts[sc, t, s]  +
            (1- step( N[sc, 1, t, 1] + hacked.counts[sc, t, 1] )) * N[sc, 1, t, 1] * 
            (-1*equals(s,1) + equals(s,2)) # change to negative value for LH when using
          NFY[sc, t, s] <- N[sc, 1, t, s] + num.hacked[sc, t, s]
          NF[sc, t, s] <- sum(N[sc, 2:4, t, s]) # number of adult nonbreeder females
          NB[sc, t, s] <- sum(N[sc, 5:7, t, s]) # number of adult breeder females
          NAD[sc, t, s] <- NF[sc, t, s] + NB[sc, t, s] # number of adults
          Ntot[sc, t, s] <- sum(N[sc, 1:7, t, s]) # total number of females
        }# s
      }  }# t sc
    
    ################################
    # Likelihood for survival
    ################################ 
    # Calculate yearly averages for integration
    for (t in 1:(nyr-1)){
      for (s in 1:nsite){
        for (xxxx in 1:surv.end[t,s]){
          # need to reorder because nimble doesn't 
          # handle nonconsecutive indices
          # yrind.surv is a matrix of indices for each site
          phiFY2[ t, s, xxxx] <- phiFY[ yrind.surv[xxxx,t,s], t]
          phiA2[ t, s, xxxx] <- phiA[ yrind.surv[xxxx,t,s], t]
          phiB2[ t, s, xxxx] <- phiB[ yrind.surv[xxxx,t,s], t]
          psiFYB2[ t, s, xxxx] <- psiFYB[ yrind.surv[xxxx,t,s], t]
          psiAB2[ t, s, xxxx] <- psiAB[ yrind.surv[xxxx,t,s], t]
          psiBA2[ t, s, xxxx] <- psiBA[ yrind.surv[xxxx,t,s], t]
          pA2[ t, s, xxxx] <- pA[ yrind.surv[xxxx,t,s], t]
          pB2[ t, s, xxxx] <- pB[ yrind.surv[xxxx,t,s], t]
        } # xxxx
        mn.phiFY[1, t, s] <- mean( phiFY2[ t, s, 1:surv.end[t,s] ] ) 
        mn.phiA[1, t, s] <- mean( phiA2[ t, s, 1:surv.end[t,s] ] )
        mn.phiB[1, t, s] <- mean( phiB2[ t, s, 1:surv.end[t,s] ] )
        mn.psiFYB[1, t, s] <- mean( psiFYB2[ t, s, 1:surv.end[t,s] ] )
        mn.psiAB[1, t, s] <- mean( psiAB2[ t, s, 1:surv.end[t,s] ] )
        mn.psiBA[1, t, s] <- mean( psiBA2[ t, s, 1:surv.end[t,s] ] )
        mn.pA[1, t, s] <- mean( pA2[ t, s, 1:surv.end[t,s] ] )
        mn.pB[1, t, s] <- mean( pB2[ t, s, 1:surv.end[t,s] ] )
        for (sc in 2:SC){
          mn.phiFY[sc, t, s] <- mn.phiFY[1, t, s] 
          mn.phiA[sc, t, s] <- mn.phiA[1, t, s] 
          mn.phiB[sc, t, s] <- mn.phiB[1, t, s] 
          mn.psiFYB[sc, t, s] <- mn.psiFYB[1, t, s]
          mn.psiAB[sc, t, s] <- mn.psiAB[1, t, s]
          mn.psiBA[sc, t, s] <- mn.psiBA[1, t, s]
          mn.pA[sc, t, s] <- mn.pA[1, t, s]
          mn.pB[sc, t, s] <- mn.pB[1, t, s]
        }
      }}
    
    for (i in 1:nind){
      for (t in 1:(nyr-1)){
        #Survival
        logit(phiFY[i,t]) <- eta[1, site[i,t],t] + 
          lmus[1, site[i,t]] + betas[1]*hacked[i]  # first year
        logit(phiA[i,t]) <- eta[2, site[i,t],t] + 
          lmus[2, site[i,t]] #+  betas[2]*hacked[i] # nonbreeder
        logit(phiB[i,t]) <- eta[3, site[i,t],t] +  
          lmus[3, site[i,t]] + betas[3]*hacked[i] # breeder
        #Recruitment
        logit(psiFYB[i,t]) <- eta[4, site[i,t],t] +  
          lmus[4, site[i,t]] #+ betas[4]*hacked[i] # first year to breeder
        logit(psiAB[i,t]) <- eta[5, site[i,t],t] + 
          lmus[5, site[i,t]] + betas[5]*hacked[i] # nonbreeder to breeder
        logit(psiBA[i,t]) <-  
          lmus[6, site[i,t]] #+ betas[6]*hacked[i] # breeder to nonbreeder
        #Re-encounter
        logit(pA[i,t]) <- eta[7, site[i,t],t] + 
          lmus[7, site[i,t]] + betas[7]*hacked[i] +
          deltas[5]*effort2[t, site[i,t]] + deltas[6]*effort2[t, site[i,t]]^2 # resight of nonbreeders
        logit(pB[i,t]) <- eta[8, site[i,t],t] + 
          lmus[8, site[i,t]] + #betas[8]*hacked[i] + 
          deltas[7]*effort2[t, site[i,t]] + deltas[8]*effort2[t, site[i,t]]^2 # resight of breeders
      }#t
    }#i
    
    # Future projections of survival, recruitment, and detection
    for(sc in 1:SC){ 
      for (t in 7:(nyr+K)){
        for (s in 1:nsite){    
          # FYs we calculate the percent hacked for each year. 
          perc.hacked[sc, t, s] <- ( num.hacked[sc, t, s] / (NFY[sc, t, s]+0.0001) ) *
            equals(s,2) # set LH to zero
          # For adults we average over the previous 5 years as an approximation
          perc.hacked.5yr[sc, t, s] <- ( sum( num.hacked[sc, (t-6):(t-1), s] ) / 
                                           (sum( NFY[sc, (t-6):(t-1), s] )+0.0001) )*
            equals(s,2) # set LH to zero 
        }}}
    
    cap <- 0.99 # cap on maximum survival
    for(sc in 1:SC){ 
      for (t in nyr:(nyr+K)){
        for (s in 1:nsite){
          # enforce the survival cap
          mn.phiFY[sc, t, s] <- step(cap - mn.phiFY1[sc, t, s]) * mn.phiFY1[sc, t, s] + 
            (1 - step(cap - mn.phiFY1[sc, t, s])) * cap
          mn.phiA[sc, t, s] <- step(cap - mn.phiA1[sc, t, s]) * mn.phiA1[sc, t, s] + 
            (1 - step(cap - mn.phiA1[sc, t, s])) * cap
          mn.phiB[sc, t, s] <- step(cap - mn.phiB1[sc, t, s]) * mn.phiB1[sc, t, s] + 
            (1 - step(cap - mn.phiB1[sc, t, s])) * cap  
          # Survival
          mn.phiFY1[sc, t, s] <- ilogit( eta[1, s, t] + 
                                           lmus[1, s] + betas[1]*perc.hacked[sc, t, s] ) + # change perc.hacked to zero at LH bc no birds translocated there
            surv.diff[sc, 1, s]  
          mn.phiA1[sc, t, s] <- ilogit( eta[2, s, t] + 
                                          lmus[2, s] ) + surv.diff[sc, 2, s]
          mn.phiB1[sc, t, s] <- ilogit( eta[3, s, t] +  
                                          lmus[3, s] + betas[3]*perc.hacked.5yr[sc, t, s] ) + # change perc.hacked to zero at LH bc no birds translocated there
            surv.diff[sc, 3, s]
          
          #Recruitment
          logit( mn.psiFYB[sc, t, s] ) <- eta[4, s, t] +  
            lmus[4, s]    # first year to breeder
          logit( mn.psiAB[sc, t, s] ) <- eta[5, s, t] +  # nonbreeder to breeder
            lmus[5, s] + betas[5]*perc.hacked.5yr[sc, t, s] # change perc.hacked to zero at LH bc no birds translocated there
          logit( mn.psiBA[sc, t, s] ) <- 
            lmus[6, s]  # breeder to nonbreeder
          #Re-encounter
          logit( mn.pA[sc, t, s] ) <- eta[7, s, t] + 
            lmus[7, s] + betas[7]*perc.hacked.5yr[sc, t, s] # resight of nonbreeders
          logit( mn.pB[sc, t, s] ) <- eta[8, s, t] +  
            lmus[8, s]  # resight of breeders
        } # s
      } # t
    } # sc
    
    # Define state-transition and observation matrices
    for (i in 1:nind){  
      # Define probabilities of state S(t+1) given S(t)
      for (t in first[i]:(nyr-1)){
        ps[1,i,t,1] <- 0
        ps[1,i,t,2] <- phiFY[i,t] * (1-psiFYB[i,t]) 
        ps[1,i,t,3] <- phiFY[i,t] * psiFYB[i,t]
        ps[1,i,t,4] <- (1-phiFY[i,t])
        
        ps[2,i,t,1] <- 0
        ps[2,i,t,2] <- phiA[i,t] * (1-psiAB[i,t])
        ps[2,i,t,3] <- phiA[i,t] * psiAB[i,t]
        ps[2,i,t,4] <- (1-phiA[i,t])
        
        ps[3,i,t,1] <- 0
        ps[3,i,t,2] <- phiB[i,t] * psiBA[i,t]
        ps[3,i,t,3] <- phiB[i,t] * (1-psiBA[i,t])
        ps[3,i,t,4] <- (1-phiB[i,t])
        
        ps[4,i,t,1] <- 0
        ps[4,i,t,2] <- 0
        ps[4,i,t,3] <- 0
        ps[4,i,t,4] <- 1
        
        # Define probabilities of O(t) given S(t)
        po[1,i,t,1] <- 1 
        po[1,i,t,2] <- 0
        po[1,i,t,3] <- 0
        po[1,i,t,4] <- 0
        
        po[2,i,t,1] <- 0
        po[2,i,t,2] <- pA[i,t]
        po[2,i,t,3] <- 0
        po[2,i,t,4] <- (1-pA[i,t])
        
        po[3,i,t,1] <- 0
        po[3,i,t,2] <- 0
        po[3,i,t,3] <- pB[i,t]
        po[3,i,t,4] <- (1-pB[i,t])
        
        po[4,i,t,1] <- 0
        po[4,i,t,2] <- 0
        po[4,i,t,3] <- 0
        po[4,i,t,4] <- 1
      } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
      # Define latent state at first capture
      z[i,first[i]] <- y[i,first[i]]
      for (t in (first[i]+1):nyr){
        # State process: draw S(t) given S(t-1)
        z[i,t] ~ dcat(ps[z[i,t-1], i, t-1, 1:4])
        # Observation process: draw O(t) given S(t)
        y[i,t] ~ dcat(po[z[i,t], i, t-1, 1:4])
      } #t
    } #i
  } )

#**********************
#* Function to run model in NIMBLE
#**********************
run_pva <- function(info, datl, constl, code){
  library('nimble')
  library('coda')
  # function for multivariate normal
  uppertri_mult_diag <- nimbleFunction(
    run = function(mat = double(2), vec = double(1)) {
      returnType(double(2))
      p <- length(vec)
      out <- matrix(nrow = p, ncol = p, init = FALSE)
      for(i in 1:p){
        out[ ,i] <- mat[ ,i] * vec[i]
      }
      return(out)
    })
  assign('uppertri_mult_diag', uppertri_mult_diag, envir = .GlobalEnv)
  
  params <- c(
    # pop growth 
    "lambda",
    # prod
    "lmu.prod", "gamma", "rr", "mn.prod", 
    # survival 
    "mus", "lmus", "betas", "deltas",
    # abundance
    "NB", "NF", "NFY", "N", "NAD", "Ntot",
    "r",
    # error terms
    "eta", "sds", "Ustar", "R",
    # yearly summaries
    'mn.phiFY', 'mn.phiA', 'mn.phiB',
    'mn.psiFYB', 'mn.psiAB', 'mn.psiBA',
    'mn.pA', 'mn.pB',
    # pva
    "extinct", "extinctAD", "extinctB"#,
    # "perc.hacked.5yr", "perc.hacked", "num.hacked", "perc.treat",
    # "numer.perc.treat", "denom.perc.treat"
  )
  #n.chains=1; n.thin=1; n.iter=50; n.burnin=25
  #n.chains=1; n.thin=50; n.iter=100000; n.burnin=50000
  n.chains=1; n.thin=100; n.iter=500000; n.burnin=400000
  
  mod <- nimbleModel(code, 
                     constants = constl, 
                     data = datl, 
                     inits = info$inits, 
                     buildDerivs = FALSE,
                     calculate=T 
  ) 
  
  cmod <- compileNimble(mod, showCompilerOutput = TRUE)
  confhmc <- configureMCMC(mod)
  confhmc$setMonitors(params)
  hmc <- buildMCMC(confhmc)
  chmc <- compileNimble(hmc, project = mod, 
                        resetFunctions = TRUE,
                        showCompilerOutput = TRUE)
  
  post <- runMCMC(chmc,
                  niter = n.iter, 
                  nburnin = n.burnin,
                  nchains = n.chains,
                  thin = n.thin,
                  samplesAsCodaMCMC = T, 
                  setSeed = info$seed)
  
  return(post)
} # run_pva function end

#*****************
#* Run chains in parallel
#*****************

this_cluster <- makeCluster(cpus)
post <- parLapply(cl = this_cluster, 
                  X = par_info_pva[1:cpus], 
                  fun = run_pva, 
                  datl = datl, 
                  constl = constl, 
                  code = mycode)
stopCluster(this_cluster)
save(post, mycode, seeds, cpus,
     file="/bsuscratch/brianrolek/riha_ipm/outputs/pva_longrun.rdata")
# save(post, mycode, seeds, cpus,
#      file="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//outputs//pva_survival_longrun.rdata")
