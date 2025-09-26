## ---- ipm --------
library('nimble')
library('parallel')
library ('coda')
load("/bsuscratch/brianrolek/riha_ipm/data.rdata")
# load("data/data.rdata")
cpus <- 5

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
      betas[k] ~ dnorm(0, sd=20)  # prior for translocations coefficients
      deltas[k] ~ dnorm(0, sd=10) # prior for survey effort coefficients
    } # k
    
    for (j in 1:8){
      for (s in 1:nsite){
        lmus[j,s] <- logit(mus[j,s])
        mus[j,s] ~ dbeta(1,1) # prior for overall means
      }}     # 
    
    # Temporal random effects and correlations between sites
    # Non-centered parameterization of the multivariate normal distribution to improve convergence
    for (jj in 1:p){ sds[jj] ~ dexp(1) }# prior for temporal variation
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
        mn.prod[t,s] <- mean( prodmat[t,s,1:pair.end[t,s]] )
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
    f.tvm.obs <- sd(brood[1:npairsobs])^2/mean(brood[1:npairsobs])
    f.tvm.rep <- sd(f.rep[1:npairsobs])^2/mean(f.rep[1:npairsobs])
    
    ################################
    # Likelihood for counts
    ################################
    # Abundance for year=1
    for (v in 1:7){ 
      for (s in 1:nsite){
        # subtract one to allow dcat to include zero
        N[v, 1, s] <- N2[v, 1, s] - 1 
        N2[v, 1, s] ~ dcat(pPrior[v, 1:s.end[v,s], s]) # Priors differ for FYs and adults
      }} # s t
    # Abundance for years > 1
    for (t in 1:(nyr-1)){
      for (s in 1:nsite){
        # Number of wild born juvs
        N[1, t+1, s] ~ dpois( (NFY[t, s]*mn.phiFY[t, s]*mn.psiFYB[t, s] + # first year breeders
                                 NF[t, s]*mn.phiA[t, s]*mn.psiAB[t, s] + # nonbreeders to breeders
                                 NB[t, s]*mn.phiB[t, s]*(1-mn.psiBA[t, s])) # breeders remaining
                              *mn.prod[t+1, s]/2 ) # end Poisson
        # Abundance of nonbreeders
        N[2, t+1, s] ~ dbin(mn.phiFY[t, s]*(1-mn.psiFYB[t, s]), NFY[t, s]) # Nestlings to nonbreeders
        N[3, t+1, s] ~ dbin(mn.phiA[t, s]*(1-mn.psiAB[t, s]), NF[t, s]) # Nonbreeders to nonbreeders
        N[4, t+1, s] ~ dbin(mn.phiB[t, s]*mn.psiBA[t, s], NB[t, s]) # Breeders to nonbreeders
        # Abundance of breeders
        N[5, t+1, s] ~ dbin(mn.phiFY[t, s]*mn.psiFYB[t, s], NFY[t, s]) # Nestlings to second year breeders
        N[6, t+1, s] ~ dbin(mn.phiA[t, s]*mn.psiAB[t, s], NF[t, s]) # Nonbreeder to breeder
        N[7, t+1, s] ~ dbin(mn.phiB[t, s]*(1-mn.psiBA[t, s]), NB[t, s]) # Breeder to breeder
      }} # s t
    
    # Count likelihoods, state-space model, and observation process    
    for (t in 1:nyr){
      for (s in 1:nsite){
        countsAdults[t, s] ~ dpois(lamAD[t, s]) # adult females 
        constraint_data[t, s] ~ dconstraint( (N[1, t, s] + hacked.counts[t, s]) >= 0 ) # constrain N1+hacked.counts to be >=0
        NFY[t, s] <- N[1, t, s] + hacked.counts[t, s] # Transfers translocated first-year females
        NF[t, s] <- sum(N[2:4, t, s]) # number of adult nonbreeder females
        NB[t, s] <- sum(N[5:7, t, s]) # number of adult breeder females
        NAD[t, s] <- NF[t, s] + NB[t, s] # number of adults
        Ntot[t, s] <- sum(N[1:7, t, s]) # total number of females
        log(lamAD[t, s]) <- log(NAD[t, s]) + deltas[1]*effort2[t, s] + deltas[2]*effort2[t, s]^2  
        log(lamFY[t, s]) <- log(N[1, t, s]) + deltas[3]*effort2[t, s] + deltas[4]*effort2[t, s]^2
      }# s
      # First-years at different sites have different distributions
      # for better model fit
      countsFY[t, 1] ~ dpois(lamFY[t, 1]) # doesn't have any zeroes so poisson fits
      countsFY[t, 2] ~ dnegbin(pp[t], r) # first year females, includes translocated/hacked
      pp[t] <- r/(r+(lamFY[t, 2] ))
    } # t
    r ~ dexp(0.05)
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
    tvm.obs[1] <- sd(c.obsB[1:nyr, 1:nsite])^2/mean(c.obsB[1:nyr, 1:nsite])
    tvm.obs[2] <- sd(c.obsFY[1:nyr, 1:nsite])^2/mean(c.obsFY[1:nyr, 1:nsite])
    tvm.rep[1] <- sd(c.repB[1:nyr, 1:nsite])^2/mean(c.repB[1:nyr, 1:nsite])
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
    
    ################################
    # Likelihood for survival
    ################################ 
    # Calculate averages for sites each year for integration
    for (t in 1:nyr){
      for (s in 1:nsite){
        for (xxxx in 1:surv.end[t,s]){
          # Reorder because nimble doesn't 
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
        mn.phiFY[ t, s] <- mean( phiFY2[ t, s, 1:surv.end[t,s] ] ) 
        mn.phiA[ t, s] <- mean( phiA2[ t, s, 1:surv.end[t,s] ] )
        mn.phiB[ t, s] <- mean( phiB2[ t, s, 1:surv.end[t,s] ] )
        mn.psiFYB[ t, s] <- mean( psiFYB2[ t, s, 1:surv.end[t,s] ] )
        mn.psiAB[ t, s] <- mean( psiAB2[ t, s, 1:surv.end[t,s] ] )
        mn.psiBA[ t, s] <- mean( psiBA2[ t, s, 1:surv.end[t,s] ] )
        mn.pA[ t, s] <- mean( pA2[ t, s, 1:surv.end[t,s] ] )
        mn.pB[ t, s] <- mean( pB2[ t, s, 1:surv.end[t,s] ] )
      }}
    
    for (i in 1:nind){
      for (t in 1:nyr){
        #Survival
        logit(phiFY[i,t]) <- eta[1, site[i,t],t] + 
          lmus[1, site[i,t]] + betas[1]*hacked[i]  # first year
        logit(phiA[i,t]) <- eta[2, site[i,t],t] +  
          lmus[2, site[i,t]] +  betas[2]*hacked[i] # nonbreeder
        logit(phiB[i,t]) <- eta[3, site[i,t],t] +  
          lmus[3, site[i,t]] + betas[3]*hacked[i] # breeder
        #Recruitment
        logit(psiFYB[i,t]) <- eta[4, site[i,t],t] +  
          lmus[4, site[i,t]] + betas[4]*hacked[i] # first year to breeder
        logit(psiAB[i,t]) <- eta[5, site[i,t],t] +  
          lmus[5, site[i,t]] + betas[5]*hacked[i] # nonbreeder to breeder
        logit(psiBA[i,t]) <-  
          lmus[6, site[i,t]] # breeder to nonbreeder
        #Re-encounter
        logit(pA[i,t]) <- eta[7, site[i,t],t] + 
          lmus[7, site[i,t]] + betas[7]*hacked[i] +
          deltas[5]*effort2[t, site[i,t]] + deltas[6]*effort2[t, site[i,t]]^2# resight of nonbreeders
        logit(pB[i,t]) <- eta[8, site[i,t],t] + 
          lmus[8, site[i,t]] + betas[8]*hacked[i] + 
          deltas[7]*effort2[t, site[i,t]] + deltas[8]*effort2[t, site[i,t]]^2# resight of breeders
      }#t
    }#i
    
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
run_ipm <- function(info, datl, constl, code){
  library('nimble')
  library('coda')

  params <- c(
    # pop growth 
    "lambda",
    # fecundity
    "lmu.prod", "gamma", "rr", "mn.prod", 
    # survival 
    "mus", "lmus", "betas", "deltas",
    # abundance
    "NB", "NF", "NFY", "N", "NAD",
    "r",
    "N", "Ntot",
    # error terms
    "eta", "sds", "Ustar", "R", "z.score",
    # yearly summaries
    'mn.phiFY', 'mn.phiA', 'mn.phiB',
    'mn.psiFYB', 'mn.psiAB', 'mn.psiBA',
    'mn.pA', 'mn.pB',
    # goodness of fit
    "f.dmape.obs", "f.dmape.rep",
    "f.tvm.obs", "f.tvm.rep",
    "dmape.obs", "dmape.rep",
    "tvm.obs", "tvm.rep",
    "tturn.obs", "tturn.rep"
  )

  #n.chains=1; n.thin=1; n.iter=100; n.burnin=50
  n.chains=1; n.thin=100; n.iter=100000; n.burnin=50000
  
  mod <- nimbleModel(code, 
                     constants = constl, 
                     data = datl, 
                     inits = info$inits, 
                     buildDerivs = FALSE, # doesn't work when TRUE, no hope for HMC
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
} # run_ipm function end

#*****************
#* Run chains in parallel
#*****************
this_cluster <- makeCluster(cpus)
post <- parLapply(cl = this_cluster, 
                  X = par_info[1:cpus], 
                  fun = run_ipm, 
                  datl = datl, 
                  constl = constl, 
                  code = mycode)
stopCluster(this_cluster)

save(post, mycode,
     file="/bsuscratch/brianrolek/riha_ipm/outputs/ipm_longrun_2025_Apr_03.rdata")

# save(post, mycode,
#      file="/bsuscratch/brianrolek/riha_ipm/outputs/ipm.rdata")