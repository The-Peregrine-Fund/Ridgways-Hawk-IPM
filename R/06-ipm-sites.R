## ---- ipm1 --------
library('nimble')
library('parallel')
library ('coda')
load("/bsuscratch/brianrolek/riha_ipm/data.rdata")
#load("data/data.rdata")
#################################
# The model
####################################################
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

mycode <- nimbleCode(
  {
    ###################################################
    # Priors and constraints
    ###################################################
    # survival, recruitment, and detection can be correlated
    for (k in 1:8){ 
      betas[k] ~ dunif(-10, 10)  # prior for coefficients
    } # k
    for (j in 1:8){ 
      for (s in 1:nsite){
        lmus[j,s] <- logit(mus[j,s])
        mus[j,s] ~ dbeta(1,1) # prior for means
      }}     # m population #s sex #h hacked 
    
    # # Temporal random effects and correlations among all sites
    # for (j in 1:p){ sds[j] ~ dexp(1) }# prior for temporal variation
    # # estimated using the multivariate normal distribution
    # R[1:p,1:p] <- t(Ustar[1:p,1:p]) %*% Ustar[1:p,1:p] # calculate rhos, correlation coefficients
    # Ustar[1:p,1:p] ~ dlkj_corr_cholesky(eta=1.3, p=p) # Ustar is the Cholesky decomposition of the correlation matrix
    # U[1:p,1:p] <- uppertri_mult_diag(Ustar[1:p, 1:p], sds[1:p])
    # # multivariate normal for temporal variance
    # for (t in 1:nyr){ # survival params only have nyr-1, no problem to simulate from however
    #   eps[1:p,t] ~ dmnorm(mu.zeroes[1:p],
    #                       cholesky = U[1:p, 1:p], prec_param = 0)
    # }

    # Temporal random effects and correlations between sites
    for (jj in 1:p2){ sds2[jj] ~ dexp(1) }# prior for temporal variation
    # estimated using the multivariate normal distribution
    R2[1:p2,1:p2] <- t(Ustar2[1:p2,1:p2]) %*% Ustar2[1:p2,1:p2] # calculate rhos, correlation coefficients
    Ustar2[1:p2,1:p2] ~ dlkj_corr_cholesky(eta=1.3, p=p2) # Ustar is the Cholesky decomposition of the correlation matrix
    U2[1:p2,1:p2] <- uppertri_mult_diag(Ustar2[1:p2, 1:p2], sds2[1:p2])
    # multivariate normal for temporal variance
    for (t in 1:nyr){ # survival params only have nyr-1, no problem to simulate from however
      for (s in 1:nsite){
      eta[1:p2,s,t] ~ dmnorm(mu.zeroes2[1:p2],
                          cholesky = U2[1:p2, 1:p2], prec_param = 0)
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
    # Likelihood for fecundity
    ###############################
    # Priors for number of fledglings
    # and nest success
    for (s in 1:nsite){
      lmu.brood[s] ~ T(dnorm(0, sd=5), 0, ) # restrict to >=1
      lmu.nest[s] <- logit(mu.nest[s])
      mu.nest[s] ~ dbeta(1, 1)
    } # s
    delta ~ dunif(-10, 10)
    sig.brood ~ dexp(1)
    gamma ~ dunif(-10, 10)
    
    # Two models for fecundity, (1) brood size and (2) nest success  
    # Brood size
    for (k in 1:nbrood){
      brood[k] ~ dlnorm(lam[k], sdlog=sig.brood) # truncating seems to break something
      lam[k] <- lmu.brood[site.brood[k]] +
                delta*treat.brood[k] +
                eps[9, year.brood[k] ] + eta[9, site.brood[k], year.brood[k] ]
    } # k
    # Nest success       
    for (n in 1:nnest){
      nest.success[n] ~ dbern( nu[n] )
      logit(nu[n]) <- lmu.nest[site.nest[n]] + 
                      gamma*treat.nest[n] + 
                      eps[10, year.nest[n] ] + eta[10, site.nest[n], year.nest[n] ]
    } # n 
    # derive yearly brood size for population model
    for (t in 1:nyr){
      for (s in 1:nsite){
      for (xx in 1:brood.end[t,s]){
        broodmat[t,s,xx] <- lam[yrind.brood[xx,t,s]]
      } # xx
      for (xxx in 1:nest.end[t,s]){
        nestmat[t,s,xxx] <- nu[yrind.nest[xxx,t,s]]
      } # xxx
      mn.brood[t,s] <- exp( mean(broodmat[t,s,1:brood.end[t,s]]) )
      mn.nest[t,s] <- mean( nestmat[t,s,1:nest.end[t,s]] )
      mn.f[t,s] <- mn.brood[t,s]*mn.nest[t,s] # calc fecundity
    }} # s t
    
    # GOF for number of fledglings
    for (k in 1:nbrood){
      f.obs[k] <- brood[k] # observed counts
      f.exp[k] <- lam[k] # expected counts adult breeder
      f.rep[k] ~ dlnorm(lam[k], sdlog=sig.brood) # expected counts
      f.dssm.obs[k] <- abs( ( f.obs[k] - f.exp[k] ) / (f.obs[k]+0.001) )
      f.dssm.rep[k] <- abs( ( f.rep[k] - f.exp[k] ) / (f.rep[k]+0.001) )
    } # k
    f.dmape.obs <- sum(f.dssm.obs[1:nbrood])
    f.dmape.rep <- sum(f.dssm.rep[1:nbrood])
    f.tvm.obs <- sd(brood[1:nbrood])^2/mean(brood[1:nbrood])
    f.tvm.rep <- sd(f.rep[1:nbrood])^2/mean(f.rep[1:nbrood])
    
    ################################
    # Likelihood for counts
    ################################
    # Abundance for year=1
    for (v in 1:7){ 
      for (s in 1:nsite){
      # subtract one because dcat can't be zero
      N[v, 1, s] <- N2[v, 1, s] - 1 
      N2[v, 1, s] ~ dcat(pPrior[1:s.end[s], s]) 
        }} # s t
    # Abundance for years > 1
    for (t in 1:(nyr-1)){
      for (s in 1:nsite){
        # Number of wild born juvs
        N[1, t+1, s] ~ dpois( (NFY[t, s]*mn.phiFY[t, s]*mn.psiFYB[t, s] + # first year breeders
                            NF[t, s]*mn.phiA[t, s]*mn.psiAB[t, s] + # nonbreeders to breeders
                            NB[t, s]*mn.phiB[t, s]*(1-mn.psiBA[t, s])) # breeders remaining
                             *mn.f[t+1, s] ) # end Poisson
        # Abundance of nonbreeders
        ## Second year nonbreeders
        N[2, t+1, s] ~ dbin(mn.phiFY[t, s]*(1-mn.psiFYB[t, s]), NFY[t, s]) # Nestlings to second year nonbreeders
        ## Adult nonbreeders
        N[3, t+1, s] ~ dbin(mn.phiA[t, s]*(1-mn.psiAB[t, s]), NF[t, s]) # Nonbreeders to nonbreeders
        N[4, t+1, s] ~ dbin(mn.phiB[t, s]*mn.psiBA[t, s], NB[t, s]) # Breeders to nonbreeders
        # Abundance of breeders
        ## Second year breeders
        N[5, t+1, s] ~ dbin(mn.phiFY[t, s]*mn.psiFYB[t, s], NFY[t, s]) # Nestlings to second year breeders
        ## Adult breeders
        N[6, t+1, s] ~ dbin(mn.phiA[t, s]*mn.psiAB[t, s], NF[t, s]) # Nonbreeder to breeder
        N[7, t+1, s] ~ dbin(mn.phiB[t, s]*(1-mn.psiBA[t, s]), NB[t, s]) # Breeder to breeder
    }} # s t
    
    # Observation process    
    for (t in 1:nyr){
      for (s in 1:nsite){
        #counts[1, t, s] ~ dpois(NFY[t, s]) # first year males, includes translocated/hacked
        counts[2, t, s] ~ dbin(mn.pA[t,s], NF[t, s]) # nonbreeding adult females 
        counts[3, t, s] ~ dbin(mn.pB[t,s], NB[t, s]) # breeding females 
        NFY[t, s] <- N[1, t, s] + hacked.counts[t, s] # Includes translocated first-year females
        NF[t, s] <- sum(N[2:4, t, s])  # number of adult nonbreeder females
        NB[t, s] <- sum(N[5:7, t, s]) # number of adult breeder females
        Ntot[t, s] <- sum(N[1:7, t, s]) # total number of females
    }# s
      counts[1, t, 1] ~ dpois(NFY[t, 1])
      counts[1, t, 2] ~ dnegbin(pp[t], r) # first year females, includes translocated/hacked
      pp[t] <- r/(r+NFY[t, 2])
      } # t
    r ~ dexp(0.1)
    ###################
    # Assess GOF of the state-space models for counts
    # Step 1: Compute statistic for observed data
    # Step 2: Use discrepancy measure: mean absolute error
    # Step 3: Use test statistic: number of turns
    ###################
    for (t in 1:nyr){
      c.repFY[t, 1] ~ dpois( NFY[t, 1] )
      c.repFY[t, 2] ~ dnegbin( pp[t], r )
    for (s in 1:nsite){
      c.expB[t, s] <- NB[t, s]  # expected counts adult breeder
      c.expA[t, s] <- NF[t, s] # nonbreeder
      c.expFY[t, s] <- NFY[t, s]
      c.obsB[t, s] <- counts[3, t, s]
      c.obsA[t, s] <- counts[2, t, s]
      c.obsFY[t, s] <- counts[1, t, s]  # first year
      c.repB[t, s] ~ dbin( mn.pB[t,s], NB[t,s] ) # simulated counts
      c.repA[t, s] ~ dbin( mn.pA[t,s], NF[t,s] )
      dssm.obsB[t, s] <- abs( ( (c.obsB[t, s]) - (c.expB[t, s]) ) / (c.obsB[t, s]+0.001)  )
      dssm.obsA[t, s] <- abs( ( (c.obsA[t, s]) - (c.expA[t, s]) ) / (c.obsA[t, s]+0.001)  )
      dssm.obsFY[t, s] <- abs( ( (c.obsFY[t, s]) - (c.expFY[t, s]) ) / (c.obsFY[t, s]+0.001)  )
      dssm.repB[t, s] <- abs( ( (c.repB[t, s]) - (c.expB[t, s]) ) / (c.repB[t, s]+0.001) )
      dssm.repA[t, s] <- abs( ( (c.repA[t, s]) - (c.expA[t, s]) ) / (c.repA[t, s]+0.001) )
      dssm.repFY[t, s] <- abs( ( (c.repFY[t, s]) - (c.expFY[t, s]) ) / (c.repFY[t, s]+0.001) )
    }} # t
    dmape.obs[1] <- sum(dssm.obsB[1:nyr, 1:nsite])
    dmape.obs[2] <- sum(dssm.obsA[1:nyr, 1:nsite])
    dmape.obs[3] <- sum(dssm.obsFY[1:nyr, 1:nsite])
    dmape.rep[1] <- sum(dssm.repB[1:nyr, 1:nsite])
    dmape.rep[2] <- sum(dssm.repA[1:nyr, 1:nsite])
    dmape.rep[3] <- sum(dssm.repFY[1:nyr, 1:nsite])
    tvm.obs[1] <- sd(c.obsB[1:nyr, 1:nsite])^2/mean(c.obsB[1:nyr, 1:nsite])
    tvm.obs[2] <- sd(c.obsA[1:nyr, 1:nsite])^2/mean(c.obsA[1:nyr, 1:nsite])
    tvm.obs[3] <- sd(c.obsFY[1:nyr, 1:nsite])^2/mean(c.obsFY[1:nyr, 1:nsite])
    tvm.rep[1] <- sd(c.repB[1:nyr, 1:nsite])^2/mean(c.repB[1:nyr, 1:nsite])
    tvm.rep[2] <- sd(c.repA[1:nyr, 1:nsite])^2/mean(c.repA[1:nyr, 1:nsite])
    tvm.rep[3] <- sd(c.repFY[1:nyr, 1:nsite])^2/mean(c.repFY[1:nyr, 1:nsite])
    # # Test statistic for number of turns
    for (s in 1:nsite){
    for (t in 1:(nyr-2)){
      tt1.obsB[t, s] <- step(c.obsB[t+2, s] - c.obsB[t+1, s])
      tt2.obsB[t, s] <- step(c.obsB[t+1, s] - c.obsB[t, s])
      tt3.obsB[t, s] <- equals(tt1.obsB[t, s] + tt2.obsB[t, s], 1)
      tt1.obsA[t, s] <- step(c.obsA[t+2, s] - c.obsA[t+1, s])
      tt2.obsA[t, s] <- step(c.obsA[t+1, s] - c.obsA[t, s])
      tt3.obsA[t, s] <- equals(tt1.obsA[t, s] + tt2.obsA[t, s], 1)
      tt1.obsFY[t, s] <- step(c.obsFY[t+2, s] - c.obsFY[t+1, s])
      tt2.obsFY[t, s] <- step(c.obsFY[t+1, s] - c.obsFY[t, s])
      tt3.obsFY[t, s] <- equals(tt1.obsFY[t, s] + tt2.obsFY[t, s], 1)
    }} # t
    tturn.obs[1] <- sum(tt3.obsB[1:(nyr-2), 1:nsite])
    tturn.obs[2] <- sum(tt3.obsA[1:(nyr-2), 1:nsite])
    tturn.obs[3] <- sum(tt3.obsFY[1:(nyr-2), 1:nsite])
    for (s in 1:nsite){
    for (t in 1:(nyr-2)){
      tt1.repB[t, s] <- step(c.repB[t+2, s] - c.repB[t+1, s])
      tt2.repB[t, s] <- step(c.repB[t+1, s] - c.repB[t, s])
      tt3.repB[t, s] <- equals(tt1.repB[t, s] + tt2.repB[t, s], 1)
      tt1.repA[t, s] <- step(c.repA[t+2, s] - c.repA[t+1, s])
      tt2.repA[t, s] <- step(c.repA[t+1, s] - c.repA[t, s])
      tt3.repA[t, s] <- equals(tt1.repA[t, s] + tt2.repA[t, s], 1)
      tt1.repFY[t, s] <- step(c.repFY[t+2, s] - c.repFY[t+1, s])
      tt2.repFY[t, s] <- step(c.repFY[t+1, s] - c.repFY[t, s])
      tt3.repFY[t, s] <- equals(tt1.repFY[t, s] + tt2.repFY[t, s], 1)
    }} # t
    tturn.rep[1] <- sum(tt3.repB[1:(nyr-2), 1:nsite])
    tturn.rep[2] <- sum(tt3.repA[1:(nyr-2), 1:nsite])
    tturn.rep[3] <- sum(tt3.repFY[1:(nyr-2), 1:nsite])
    
    ################################
    # Likelihood for survival
    ################################ 
    # Calculate an averages for sites
    # each year 
    for (t in 1:nyr){
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
        logit(phiFY[i,t]) <- eta[1, site[i,t],t] + eps[1,t] + 
                              lmus[1, site[i,t]] + betas[1]*hacked[i]  # first year
        logit(phiA[i,t]) <- eta[2, site[i,t],t] + eps[2,t] + 
                              lmus[2, site[i,t]] +  betas[2]*hacked[i] # nonbreeder
        logit(phiB[i,t]) <- eta[3, site[i,t],t] + eps[3,t] + 
                              lmus[3, site[i,t]] + betas[3]*hacked[i] # breeder
        #Recruitment
        logit(psiFYB[i,t]) <- eta[4, site[i,t],t] + eps[4,t] + 
                              lmus[4, site[i,t]] + betas[4]*hacked[i] # first year to breeder
        logit(psiAB[i,t]) <- eta[5, site[i,t],t] + eps[5,t] + 
                              lmus[5, site[i,t]] + betas[5]*hacked[i] # nonbreeder to breeder
        logit(psiBA[i,t]) <- eta[6, site[i,t],t] + eps[6,t] + 
                              lmus[6, site[i,t]] + betas[6]*hacked[i] # breeder to nonbreeder
        #Re-encounter
        logit(pA[i,t]) <- eta[7, site[i,t],t] + eps[7,t] + 
                              lmus[7, site[i,t]] + betas[7]*hacked[i] # resight of nonbreeders
        logit(pB[i,t]) <- eta[8, site[i,t],t] + eps[8,t] + 
                              lmus[8, site[i,t]] + betas[8]*hacked[i] # resight of breeders
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

run_ipm <- function(info, datl, constl, code){
  library('nimble')
  library('coda')
  # helper function for multivariate normal
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
  
params <- c(# pop growth 
            "lambda",
            # fecundity
            "lmu.brood", "delta", "sig.brood", 
            "lmu.nest", "mu.nest", "gamma", 
            "mn.brood", "mn.nest", "mn.f", 
             # survival 
             "mus", "lmus", "betas",
             # abundance
             "NB", "NF", "NFY", "N",
             "r",
             "N", "Ntot",
             # error terms
             "eps", "sds", "Ustar", "U", "R",
             "eta", "sds2", "Ustar2", "U2", "R2",
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

# n.chains=1; n.thin=200; n.iter=500000; n.burnin=300000
#n.chains=1; n.thin=50; n.iter=100000; n.burnin=50000
n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
# n.chains=1; n.thin=1; n.iter=200; n.burnin=100

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

this_cluster <- makeCluster(4)
post <- parLapply(cl = this_cluster, 
                  X = par_info, 
                  fun = run_ipm, 
                  datl = datl, 
                  constl = constl, 
                  code = mycode)
stopCluster(this_cluster)

save(post, mycode,
     file="/bsuscratch/brianrolek/riha_ipm/outputs/ipm_sites_hmc.rdata")

# save(post,  
#      file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm-sites.Rdata")
