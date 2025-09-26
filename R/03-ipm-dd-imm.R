## ---- ipm --------
library('nimble')
library('parallel')
library ('coda')
library ('nimbleEcology')
source("/bsuscratch/brianrolek/riha_ipm/MCMCvis.R")
load("/bsuscratch/brianrolek/riha_ipm/data-dd.rdata")
#load("N-inits.rdata")
#load("/bsuscratch/brianrolek/riha_ipm/outputs/ipm-find-inits.rdata")
load("/bsuscratch/brianrolek/riha_ipm/outputs/ipm-2025-09-21.rdata")
outp <- MCMCpstr(post, type="chains")
#load("/bsuscratch/brianrolek/riha_ipm/outputs/ipm_longrun_2025_Apr_03.rdata")
#load("/bsuscratch/brianrolek/riha_ipm/outputs/ipm_shortrun.rdata")
#load("/bsuscratch/brianrolek/riha_ipm/outputs/ipm_dd_long_2025_Apr_03.rdata")
#outp2 <- MCMCpstr(post[-5], type="chains")
#constl$Ntot <- apply(outp2$Ntot, c(1,2), mean)
constl$p <- 8
cpus <- 4

# library (MCMCvis)
# load("data/data-dd.rdata")
# load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm-find-inits.rdata")
# load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_shortrun.rdata")
# outp <- screen(post)
# rm(list=c("mycode", "post"))
# cpus <- 5
# info <- par_info[[1]]

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
    for (k in 1:7){
      #betas[k] ~ dnorm(0, sd=10)  # prior for translocations coefficients
      betas[k] <- 0
    }
    for (ww in 1:8){
      deltas[ww] ~ dnorm(0, sd=10) # prior for survey effort coefficients
      #deltas[ww] <- 0
    } # k
    
    for (j in 1:8){
      for (s in 1:nsite){
        lmus[j,s] <- logit(mus[j,s])
        mus[j,s] ~ dbeta(1,1) # prior for overall means
      }}     # 
    
    for (qq in 1:6){
      for (s in 1:nsite){
        alphas[qq,s] ~ dnorm(0, sd=5)
      }} # qq s
    
    # Temporal random effects and correlations between sites
    # Non-centered parameterization of the multivariate normal distribution to improve convergence
    for (jj in 1:p){ sds[jj] ~ T(dnorm(0, sd=2), 0, ) }# prior for temporal variation
    R[1:p,1:p] <- t(Ustar[1:p,1:p]) %*% Ustar[1:p,1:p] # calculate rhos, correlation coefficients
    Ustar[1:p,1:p] ~ dlkj_corr_cholesky(eta=1.1, p=p) # Ustar is the Cholesky decomposition of the correlation matrix
    for (t in 1:nyr){ # survival params only have nyr-1, no problem to simulate from however
      for (s in 1:nsite){
        eta[1:p,s,t] <- diag(sds[1:p]) %*% t(Ustar[1:p,1:p]) %*% z.score[1:p,s,t]
        for(j in 1:p){
          z.score[j,s,t] ~ dnorm(0, sd=1)  # z-scores
        } # j
      } } # s t 
    
    #######################
    # Derived params
    #######################
    for (s in 1:nsite){
      for (t in 1:(nyr-1)){
        lambda[t, s] <-  Ntot[t+1, s]/(Ntot[t, s])
        loglambda[t, s] <- log(lambda[t, s])
      }} #t
    
    ###############################
    # Likelihood for productivity
    ###############################
    # Priors
    for (s in 1:nsite){
      lmu.prod[s] ~ dnorm(0, sd=5)
    } # s
    gamma ~ dnorm(0, sd=10)
    rr ~ dgamma(0.01,0.01)
    for (s in 1:nsite){
    lmu.treat[s] <- logit(mu.treat[s])
    mu.treat[s] ~ dbeta(1,1)
    }
    
    # Productivity likelihood      
    for (k in 1:npairsobs){
      prod[k] ~ dnegbin(ppp[k], rr)
      ppp[k] <- rr/(rr+mu.prod[k])
      log(mu.prod[k]) <- mn.prod[treat.pair[k] + 1, year.pair[k], site.pair[k], ] +
        # covariate as sum to zero contrasts
        # intercept remains the overall mean
        gamma*c(-1,1)[ treat.pair[k] + 1]
    } # k
    for (t in 1:(nyr-1)){ # this index for dens dep ipm bc mn.prod[t-1] would cause a zero index
      for (s in 1:nsite){
        for (tr in 1:2){
        mn.prod[tr, t, s] <-  lmu.prod[s] + 
            #alphas[6, s]*(Ntot[t, s]-mnC[s]) + 
            eta[7, s, t]
        }
        num.treat[t,s] ~ dpois( ptreat[t,s]*NB[t,s] ) # binomial doesn't sample
        logit(ptreat[t,s]) <- lmu.treat[s]
      }} # s t
    
    # GOF for productivity
    for (k in 1:npairsobs){
      f.obs[k] <- prod[k] # observed counts
      f.exp[k] <- mu.prod[k] # expected counts adult breeder
      f.rep[k] ~ dnegbin(ppp[k], rr) # expected counts
      f.dssm.obs[k] <- abs( ( f.obs[k] - f.exp[k] ) / (f.obs[k]+0.001) )
      f.dssm.rep[k] <- abs( ( f.rep[k] - f.exp[k] ) / (f.rep[k]+0.001) )
    } # k
    f.dmape.obs <- sum(f.dssm.obs[1:npairsobs])
    f.dmape.rep <- sum(f.dssm.rep[1:npairsobs])
    f.tvm.obs <- sd(f.obs[1:npairsobs])^2/mean(f.obs[1:npairsobs])
    f.tvm.rep <- sd(f.rep[1:npairsobs])^2/mean(f.rep[1:npairsobs])

    ################################
    # Immigration
    ################################
    omega ~ dunif(0,1) # assumes a maximum of 1:1 immigrant per adult
    ################################
    # Likelihood for counts
    ################################
    # Abundance for year=1
    for (s in 1:nsite){
    for (v in 1:10){
        # subtract one to allow dcat to include zero
        N[v, 1, s] <- N2[v, 1, s] - 1
        N2[v, 1, s] ~ dcat(pPrior[v, 1:s.end[v,s], s]) # Priors differ for FYs and adults
      } # t
      } # s
    # Abundance for years > 1
    for (t in 1:(nyr-1)){
      for (s in 1:nsite){
        # Number of wild born juvs
        # treated nests
        N[1, t+1, s] ~ dpois(
          (NFY[t, s]*mn.phiFY[t, s]*mn.psiFYB[t, s] + # first year breeders
             NF[t, s]*mn.phiA[t, s]*mn.psiAB[t, s] + # nonbreeders to breeders
             NB[t, s]*mn.phiB[t, s]*(1-mn.psiBA[t, s])) * # breeders remaining
            ptreat[t, s]*
            mn.prod[2, t, s]/2 )
        # untreated nests
        #N[2, t+1, s] <- 0
        N[2, t+1, s] ~ dpois(
          (NFY[t, s]*mn.phiFY[t, s]*mn.psiFYB[t, s] + # first year breeders
          NF[t, s]*mn.phiA[t, s]*mn.psiAB[t, s] + # nonbreeders to breeders
          NB[t, s]*mn.phiB[t, s]*(1-mn.psiBA[t, s])) * # breeders remaining
                 (1-ptreat[t, s]) *
                mn.prod[1, t, s]/2 )
         # end Poisson
        # mn.prod should be t+1 in a postbreeding ipm
        # but here we specify as t to prevent an acyclic graph
        # and we adjust time indices in the productivity model
        # to account for this
        # Abundance of nonbreeders
        N[3, t+1, s] ~ dbin(mn.phiFY[t, s]*(1-mn.psiFYB[t, s]), NFY[t, s]) # Nestlings to nonbreeders
        N[4, t+1, s] ~ dbin(mn.phiA[t, s]*(1-mn.psiAB[t, s]), NF[t, s]) # Nonbreeders to nonbreeders
        N[5, t+1, s] ~ dbin(mn.phiB[t, s]*mn.psiBA[t, s], NB[t, s]) # Breeders to nonbreeders
        # skips N[5,t] see immigrants below
        # Abundance of breeders
        N[7, t+1, s] ~ dbin(mn.phiFY[t, s]*mn.psiFYB[t, s], NFY[t, s]) # Nestlings to second year breeders
        N[8, t+1, s] ~ dbin(mn.phiA[t, s]*mn.psiAB[t, s], NF[t, s]) # Nonbreeder to breeder
        N[9, t+1, s] ~ dbin(mn.phiB[t, s]*(1-mn.psiBA[t, s]), NB[t, s]) # Breeder to breeder
      } # s
      N[6, t+1, 1] ~ dpois( NF[t, 1]*omega ) # number of nonbreeder immigrants, assumes similar recruitment as natal birds
      N[10, t+1, 1] ~ dpois( NB[t, 1]*omega ) # number of breeder immigrants, assumes similar recruitment as natal birds
      #N[6, t+1, 1] <- 0 # immigrants to nonbreeder
      N[6, t+1, 2] <- 0 # set to zero for PC, isolated
      #N[10, t+1, 1] <- 0 # immigrants to breeder
      N[10, t+1, 2] <- 0 # set to zero for PC, isolated
    } # t

    # Count likelihoods, state-space model, and observation process
    for (t in 1:nyr){
      for (s in 1:nsite){
        NFY2[t,s] <- N[1, t, s] + N[2, t, s]
        constraint_data[t, s] ~ dconstraint( (NFY2[t,s] + hacked.counts[t, s]) >= 0 ) # constrain N1+hacked.counts to be >=0
        NFY[t, s] <- NFY2[t,s] + hacked.counts[t, s] # Transfers translocated first-year females
        
        # constraint_data[t, s] ~ dconstraint( (N[1, t, s] + hacked.counts[t, s]) >= 0 ) # constrain N1+hacked.counts to be >=0
        # NFY[t, s] <- N[1, t, s] + hacked.counts[t, s] # Transfers translocated first-year females
        NF[t, s] <- sum(N[3:6, t, s]) # number of adult nonbreeder females
        NB[t, s] <- sum(N[7:10, t, s]) # number of adult breeder females
        NAD[t, s] <- NF[t, s] + NB[t, s]  # number of adults
        Ntot[t, s] <- sum(N[1:10, t, s]) # total number of females
        countsAdults[t, s] ~ dpois(lamAD[t, s]) # adult females, includes nonbreeders and breeders but its a temp approx. to help with inits
        log(lamAD[t, s]) <- log(NAD[t, s]) #+ deltas[1]*effort2[t, s] + deltas[2]*effort2[t, s]^2
        log(lamFY[t, s]) <- log(NFY2[ t, s]) #+ deltas[3]*effort2[t, s] + deltas[4]*effort2[t, s]^2
      }# s
      # First-years at different sites have different distributions
      # for better model fit
      countsFY[t, 1] ~ dpois(lamFY[t, 1]) # doesn't have any zeroes so poisson fits
      countsFY[t, 2] ~ dnegbin(pp[t], r) # first year females, includes translocated/hacked
      pp[t] <- r/(r+(lamFY[t, 2] ))
    } # t
    r ~ dgamma(0.01,0.01)
    ###################
    # Assess GOF of the state-space models for counts
    # Step 1: Compute statistic for observed data
    # Step 2: Use discrepancy measure: mean absolute error
    # Step 3: Use test statistic: number of turns
    ###################
    for (t in 1:nyr){
      c.repFY[t, 1] ~ dpois( lamFY[t, 1] )
      c.repFY[t, 2] ~ dnegbin( pp[t], r )
      for (s in 1:nsite){
        c.expAD[t, s] <- lamAD[t, s]  # expected counts adult breeder
        c.expFY[t, s] <- lamFY[t, s]
        c.obsAD[t, s] <- countsAdults[t, s]
        c.obsFY[t, s] <- countsFY[t, s]  # first year
        c.repAD[t, s] ~ dpois( lamAD[t, s] ) # simulated counts
        dssm.obsAD[t, s] <- abs( ( (c.obsAD[t, s]) - (c.expAD[t, s]) ) / (c.obsAD[t, s]+0.001)  )
        dssm.obsFY[t, s] <- abs( ( (c.obsFY[t, s]) - (c.expFY[t, s]) ) / (c.obsFY[t, s]+0.001)  )
        dssm.repAD[t, s] <- abs( ( (c.repAD[t, s]) - (c.expAD[t, s]) ) / (c.repAD[t, s]+0.001) )
        dssm.repFY[t, s] <- abs( ( (c.repFY[t, s]) - (c.expFY[t, s]) ) / (c.repFY[t, s]+0.001) )
      }} # t
    dmape.obs[1] <- sum(dssm.obsAD[1:nyr, 1:nsite])
    dmape.obs[2] <- sum(dssm.obsFY[1:nyr, 1:nsite])
    dmape.rep[1] <- sum(dssm.repAD[1:nyr, 1:nsite])
    dmape.rep[2] <- sum(dssm.repFY[1:nyr, 1:nsite])
    tvm.obs[1] <- sd(c.obsAD[1:nyr, 1:nsite])^2/mean(c.obsAD[1:nyr, 1:nsite])
    tvm.obs[2] <- sd(c.obsFY[1:nyr, 1:nsite])^2/mean(c.obsFY[1:nyr, 1:nsite])
    tvm.rep[1] <- sd(c.repAD[1:nyr, 1:nsite])^2/mean(c.repAD[1:nyr, 1:nsite])
    tvm.rep[2] <- sd(c.repFY[1:nyr, 1:nsite])^2/mean(c.repFY[1:nyr, 1:nsite])
    # # Test statistic for number of turns
    for (s in 1:nsite){
      for (t in 1:(nyr-2)){
        tt1.obsAD[t, s] <- step(c.obsAD[t+2, s] - c.obsAD[t+1, s])
        tt2.obsAD[t, s] <- step(c.obsAD[t+1, s] - c.obsAD[t, s])
        tt3.obsAD[t, s] <- equals(tt1.obsAD[t, s] + tt2.obsAD[t, s], 1)
        tt1.obsFY[t, s] <- step(c.obsFY[t+2, s] - c.obsFY[t+1, s])
        tt2.obsFY[t, s] <- step(c.obsFY[t+1, s] - c.obsFY[t, s])
        tt3.obsFY[t, s] <- equals(tt1.obsFY[t, s] + tt2.obsFY[t, s], 1)
      }} # t
    tturn.obs[1] <- sum(tt3.obsAD[1:(nyr-2), 1:nsite])
    tturn.obs[2] <- sum(tt3.obsFY[1:(nyr-2), 1:nsite])
    for (s in 1:nsite){
      for (t in 1:(nyr-2)){
        tt1.repAD[t, s] <- step(c.repAD[t+2, s] - c.repAD[t+1, s])
        tt2.repAD[t, s] <- step(c.repAD[t+1, s] - c.repAD[t, s])
        tt3.repAD[t, s] <- equals(tt1.repAD[t, s] + tt2.repAD[t, s], 1)
        tt1.repFY[t, s] <- step(c.repFY[t+2, s] - c.repFY[t+1, s])
        tt2.repFY[t, s] <- step(c.repFY[t+1, s] - c.repFY[t, s])
        tt3.repFY[t, s] <- equals(tt1.repFY[t, s] + tt2.repFY[t, s], 1)
      }} # t
    tturn.rep[1] <- sum(tt3.repAD[1:(nyr-2), 1:nsite])
    tturn.rep[2] <- sum(tt3.repFY[1:(nyr-2), 1:nsite])

    # ###############################
    # Likelihood for survival
    # ###############################
    # Calculate averages for sites each year for integration
    # for (t in 1:nyr){
    #   for (s in 1:nsite){
    #     for (xxxx in 1:surv.end[t,s]){
    #       # Reorder because nimble doesn't
    #       # handle nonconsecutive indices
    #       # yrind.surv is a matrix of indices for each site
    #       phiFY2[ t, s, xxxx] <- phiFY[ yrind.surv[xxxx,t,s], t]
    #       phiA2[ t, s, xxxx] <- phiA[ yrind.surv[xxxx,t,s], t]
    #       phiB2[ t, s, xxxx] <- phiB[ yrind.surv[xxxx,t,s], t]
    #       psiFYB2[ t, s, xxxx] <- psiFYB[ yrind.surv[xxxx,t,s], t]
    #       psiAB2[ t, s, xxxx] <- psiAB[ yrind.surv[xxxx,t,s], t]
    #       psiBA2[ t, s, xxxx] <- psiBA[ yrind.surv[xxxx,t,s], t]
    #       pA2[ t, s, xxxx] <- pA[ yrind.surv[xxxx,t,s], t]
    #       pB2[ t, s, xxxx] <- pB[ yrind.surv[xxxx,t,s], t]
    #     } # xxxx
    #     mn.phiFY[ t, s] <- mean( phiFY2[ t, s, 1:surv.end[t,s] ] )
    #     mn.phiA[ t, s] <- mean( phiA2[ t, s, 1:surv.end[t,s] ] )
    #     mn.phiB[ t, s] <- mean( phiB2[ t, s, 1:surv.end[t,s] ] )
    #     mn.psiFYB[ t, s] <- mean( psiFYB2[ t, s, 1:surv.end[t,s] ] )
    #     mn.psiAB[ t, s] <- mean( psiAB2[ t, s, 1:surv.end[t,s] ] )
    #     mn.psiBA[ t, s] <- mean( psiBA2[ t, s, 1:surv.end[t,s] ] )
    #     mn.pA[ t, s] <- mean( pA2[ t, s, 1:surv.end[t,s] ] )
    #     mn.pB[ t, s] <- mean( pB2[ t, s, 1:surv.end[t,s] ] )
    #   }}
    
    for (t in 1:(nyr-1)){
      for (s in 1:nsite){
        logit(mn.phiFY[t,s]) <- lmn.phiFY[t,s]
        logit(mn.phiA[t,s]) <- lmn.phiA[t,s]
        logit(mn.phiB[t,s]) <- lmn.phiB[t,s]
        logit(mn.psiFYB[t,s]) <- lmn.psiFYB[t,s]
        logit(mn.psiAB[t,s]) <- lmn.psiAB[t,s]
        logit(mn.psiBA[t,s]) <- lmn.psiBA[t,s]

        lmn.phiFY[t,s] <- lmus[1,s] + eta[1,s,t] 
          #alphas[1,s]*(Ntot[t, s]-mnC[s])
        lmn.phiA[t,s] <- lmus[2,s]  #eta[2,s,t] + KEEP OFF
          #alphas[2,s]*(Ntot[t, s]-mnC[s])
        lmn.phiB[t,s] <- lmus[3,s] + eta[3,s,t] 
          #alphas[3,s]*(Ntot[t, s]-mnC[s])
        lmn.psiFYB[t,s] <- lmus[4,s] + eta[4,s,t] 
          #alphas[4,s]*(Ntot[t, s]-mnC[s])
        lmn.psiAB[t,s] <- lmus[5,s] + eta[5,s,t] 
          #alphas[5,s]*(Ntot[t, s]-mnC[s])
        lmn.psiBA[t,s] <- lmus[6,s]
      } # s
    for (i in 1:nind){
      # include individual effects here
        #Survival
        logit(phiFY[i,t]) <- 
          # lmus[1,site[i,t]] + 
          # eta[1,site[i,t],t] +
          # alphas[1, site[i,t]]*(Ntot[t, site[i,t]]-mnC[site[i,t]]) +
          # betas[1]*hacked[i] # first year
          lmn.phiFY[t, site[i,t]] +
          betas[1]*c(-1, 1)[ hacked[i]+1 ]
        logit(phiA[i,t]) <- 
          # lmus[2,site[i,t]] + 
          # eta[2,site[i,t],t] +
          # alphas[2, site[i,t]]*(Ntot[t, site[i,t]]-mnC[site[i,t]])  +
          # betas[2]*hacked[i] # nonbreeder
          lmn.phiA[t, site[i,t]] 
          #betas[2]*c(-1, 1)[ hacked[i]+1 ]
        logit(phiB[i,t]) <- 
          # lmus[3,site[i,t]] + 
          # eta[3,site[i,t],t] +
          # alphas[3, site[i,t]]*(Ntot[t, site[i,t]]-mnC[site[i,t]])  +
          # betas[3]*hacked[i]  # breeder
          lmn.phiB[t, site[i,t]] +
          betas[3]*c(-1, 1)[ hacked[i]+1 ]
          
        #Recruitment
        logit(psiFYB[i,t]) <- 
          # lmus[4,site[i,t]] + 
          # eta[4,site[i,t],t] +
          # alphas[4, site[i,t]]*(Ntot[t, site[i,t]]-mnC[site[i,t]])  +
          # betas[4]*hacked[i]  # first year to breeder
          lmn.psiFYB[t, site[i,t]] 
          #betas[4]*c(-1, 1)[ hacked[i]+1 ]
        
        logit(psiAB[i,t]) <- 
          # lmus[5,site[i,t]] + 
          # eta[5,site[i,t],t] +
          # alphas[5, site[i,t]]*(Ntot[t, site[i,t]]-mnC[site[i,t]]) +
          # betas[5]*hacked[i] # nonbreeder to breeder
          lmn.psiAB[t, site[i,t]] +
          betas[5]*c(-1, 1)[ hacked[i]+1 ]
        logit(psiBA[i,t]) <- 
          #lmus[6, site[i,t]]  # breeder to nonbreeder
          lmn.psiBA[t, site[i,t]]
        #Re-encounter
        logit(pA[i,t]) <- 
          lmus[7, site[i,t]] +  
          #betas[6]*c(-1, 1)[ hacked[i]+1 ] + KEEP OFF
          deltas[5]*effort2[t, site[i,t]] + 
          deltas[6]*effort2[t, site[i,t]]^2# resight of nonbreeders
        logit(pB[i,t]) <- 
          lmus[8, site[i,t]] +
          #eta[6, site[i,t],t] +
          betas[7]*c(-1, 1)[ hacked[i]+1 ] + 
          deltas[7]*effort2[t, site[i,t]] + 
          deltas[8]*effort2[t, site[i,t]]^2# resight of breeders
      }#t
    }#i
    
    # Define state-transition and observation matrices
    for (i in 1:nind){  
      # Define probabilities of state S(t+1) given S(t)
      for (t in first[i]:(nyr-1)){
        ps[i,1,1,t] <- 0
        ps[i,1,2,t] <- phiFY[i,t] * (1-psiFYB[i,t]) 
        ps[i,1,3,t] <- phiFY[i,t] * psiFYB[i,t]
        ps[i,1,4,t] <- (1-phiFY[i,t])
        
        ps[i,2,1,t] <- 0
        ps[i,2,2,t] <- phiA[i,t] * (1-psiAB[i,t])
        ps[i,2,3,t] <- phiA[i,t] * psiAB[i,t]
        ps[i,2,4,t] <- (1-phiA[i,t])
        
        ps[i,3,1,t] <- 0
        ps[i,3,2,t] <- phiB[i,t] * psiBA[i,t]
        ps[i,3,3,t] <- phiB[i,t] * (1-psiBA[i,t])
        ps[i,3,4,t] <- (1-phiB[i,t])
        
        ps[i,4,1,t] <- 0
        ps[i,4,2,t] <- 0
        ps[i,4,3,t] <- 0
        ps[i,4,4,t] <- 1
        
        # Define probabilities of O(t) given S(t)
        po[i,1,1,t] <- 1 
        po[i,1,2,t] <- 0
        po[i,1,3,t] <- 0
        po[i,1,4,t] <- 0
        
        po[i,2,1,t] <- 0
        po[i,2,2,t] <- pA[i,t]
        po[i,2,3,t] <- 0
        po[i,2,4,t] <- (1-pA[i,t])
        
        po[i,3,1,t] <- 0
        po[i,3,2,t] <- 0
        po[i,3,3,t] <- pB[i,t]
        po[i,3,4,t] <- (1-pB[i,t])
        
        po[i,4,1,t] <- 0
        po[i,4,2,t] <- 0
        po[i,4,3,t] <- 0
        po[i,4,4,t] <- 1
      } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
      ps[i,1:4,1:4,nyr] <- matrix(0,4,4)
      y[i, (first[i] + 1):nyr] ~ 
        dDHMMo(init = ps[i, y.first[i], 1:4, first[i]], 
               probObs = po[i, 1:4, 1:4, first[i]:(nyr - 1)], 
               probTrans = ps[i, 1:4, 1:4, (first[i] + 1):(nyr - 1)], 
               len = nyr - first[i], 
               checkRowSums = 0)
    } #i
  } )

#**********************
#* Function to run model in NIMBLE
#**********************
run_ipm <- function(info, datl, constl, code, outp){
  library('nimble')
  library('coda')
  library ('nimbleEcology')
  #library ('MCMCvis')
  source("/bsuscratch/brianrolek/riha_ipm/MCMCvis.R")
  source("/bsuscratch/brianrolek/riha_ipm/check_N_inits.R")
  
  params <- c(
    # fecundity
    "lmu.prod", "gamma", "rr", "mn.prod", 
    #"mu.treat", "ptreat",
    # survival 
    "mus", "lmus", "betas", "deltas", "alphas",
    # abundance
    #"lNAD", "lN", "lamAD",
    "NB", "NF", "NFY", "N", "NAD",
    "r",
    "N", "Ntot",
    "lambda",
    "omega",
    # error terms
    "eta", "sds", "Ustar", "R", "z.score",
    # yearly summaries
    'mn.phiFY', 'mn.phiA', 'mn.phiB',
    'mn.psiFYB', 'mn.psiAB', 'mn.psiBA',
    #'mn.pA', 'mn.pB',
    # goodness of fit
    "f.dmape.obs", "f.dmape.rep",
    "f.tvm.obs", "f.tvm.rep",
    "dmape.obs", "dmape.rep",
    "tvm.obs", "tvm.rep",
    "tturn.obs", "tturn.rep"
  )
  
  # initial values  
  inits.func <- function (){
    # sample from working inits
    # to find good inits for N
    # see function provided in file check.N.inits
    # which specifies the get_inits() function
    list(
      # fecundity
      lmu.prod = apply(outp$lmu.prod, 1, mean),
      gamma = mean(outp$gamma),
      rr = mean(outp$rr),
      mu.treat = c(0.80, 0.55),
      # survival
      mus = apply(outp$mus, c(1,2), mean),
      #betas = apply(outp$betas, 1, mean),
      deltas = c(0,0,0,0, apply(outp$deltas, 1, mean)[5:8]),
      alphas = apply(outp$alphas, c(1,2), mean),
      sds = apply(outp$sds, 1, mean),
      Ustar = apply(outp$Ustar, c(1,2), mean),
      z.score = apply(outp$z.score , c(1,2,3), mean),
      # counts
      r = mean(outp$r),
      N = get_inits2() , # sample from inits of chains that worked
      # immigration
      omega = runif(1, 0.01, 0.1)
    )}
  inits <- inits.func()
  
  mod <- nimbleModel(code, 
                     constants = constl, 
                     data = datl, 
                     inits = inits, 
                     buildDerivs = FALSE, # doesn't work when TRUE, no hope for HMC
                     calculate=T ) 
  cmod <- compileNimble(mod, showCompilerOutput = TRUE)
  confhmc <- configureMCMC(cmod)
  confhmc$setMonitors(params)
  hmc <- buildMCMC(confhmc)
  chmc <- compileNimble(hmc, project = cmod, 
                        resetFunctions = TRUE,
                        showCompilerOutput = TRUE)
  #nc=1; nt=100; ni=50000; nb=25000
  #nc=1; nt=100; ni=100000; nb=50000
  #nc=1; nt=100; ni2=250000; nb=150000
  #nc=1; nt=1; ni=3; nb=1
  nc=1; nt=1; ni=200; nb=0
  
  post <- runMCMC(chmc,
                   niter = ni,
                   nburnin = nb,
                   nchains = nc,
                   thin = nt,
                   samplesAsCodaMCMC = T,
                   setSeed = info$seed,
                   inits = inits)
  
  return(post)
} # run_ipm function end

# #*****************
# #* Run chains in parallel
# #*****************
this_cluster <- makeCluster(cpus)
post <- parLapply(cl = this_cluster,
                  X = par_info[1:cpus],
                  fun = run_ipm,
                  datl = datl,
                  constl = constl,
                  code = mycode,
                  outp = outp)
stopCluster(this_cluster)

flnm2 <- paste0("/bsuscratch/brianrolek/riha_ipm/outputs/ipm-imm",
                substr(Sys.time(), 1, 10), ".rdata")
save(post, mycode, 
     file = flnm2)
# 

# save(post, mycode,
#   file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm-simp-2025-05-02.rdata")
