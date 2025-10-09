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
constl$p <- 9
cpus <- 4

library (MCMCvis)
load("data/data-dd.rdata")
# load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm-find-inits.rdata")
# load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_shortrun.rdata")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm-2025-09-21.rdata")
outp <- MCMCpstr(post, type="chains")
rm(list=c("mycode", "post"))
# cpus <- 5
# info <- par_info[[1]]
constl$p <- 9
cpus <- 4
datl$hacked <- constl$hacked
constl <- constl[-18]
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
    for (j in 1:8){
      for (s in 1:nsite){
        lmus[j,s] <- logit(mus[j,s])
        mus[j,s] ~ dbeta(1,1) # prior for overall means
      }}     
    
    for (k in 1:7){
      betas[k] ~ dnorm(0, sd=2)  # prior for translocations coefficients
      #betas[k] ~ dunif(-12, 5)
      #betas[k] <- 0
    }
    for (ww in 1:8){
      deltas[ww] ~ dnorm(0, sd=2) # prior for survey effort coefficients
      #deltas[ww] ~ dunif(-3, 3)
      #deltas[ww] <- 0
    } # k
    
    for (qq in 1:6){
      for (s in 1:nsite){
        alphas[qq,s] ~ dnorm(0, sd=2)
        #alphas[qq,s] ~ dunif(-2, 2)
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
      lmu.prod[s] ~ dunif(-4.6, 0.5)
    } # s
    #gamma ~ dnorm(0, sd=2)
    gamma ~ dunif(-2,2)
    rr ~ dexp(0.05)
    for (s in 1:nsite){
    lmu.treat[s] <- logit(mu.treat[s])
    mu.treat[s] ~ dbeta(1,1)
    lmu.hacked[s] <- logit(mu.hacked[s])
    mu.hacked[s] ~ dbeta(1,1)
    }
    
    # Productivity likelihood      
    for (k in 1:npairsobs){
      prod[k] ~ dnegbin(ppp[k], rr)
      ppp[k] <- rr/(rr+mu.prod[k])
      log(mu.prod[k]) <- mn.prod[treat.pair[k] + 1,  year.pair[k], site.pair[k]]
        # covariate as sum to zero contrasts
        # intercept remains the overall mean
        #gamma*c(-1,1)[ treat.pair[k] + 1]
    } # k
    for (t in 1:(nyr-1)){ # this index for dens dep ipm bc mn.prod[t-1] would cause a zero index
      for (s in 1:nsite){
        for (tr in 1:2){
        mn.prod[tr, t, s] <-  lmu.prod[s] + 
            gamma*c(-1,1)[ tr ] #+
            #alphas[6, s]*(Ntot[t, s]-mnC[s]) + 
            #eta[7, s, t]
        }
        # percent of breeding territories treated for parasitic bot flies
        num.treat[t,s] ~ dpois( ptreat[t,s]*NB[t,s] ) # binomial doesn't sample
        logit(ptreat[t,s]) <- lmu.treat[s] #+ 
          #eta[8, s, t]
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
    omega ~ dbeta(1,1) # assumes a maximum of 1:1 immigrant per adult
    ################################
    # Likelihood for counts
    ################################
    # Abundance for year=1
    for (s in 1:nsite){
    for (v in 1:8){
        # subtract one to allow dcat to include zero
        N[v, 1, s] <- N2[v, 1, s] - 1
        N2[v, 1, s] ~ dcat(pPrior[v, 1:s.end[v,s], s]) # Priors differ for FYs and adults
      } # t
      } # s
    # Abundance for years > 1
    for (t in 1:(nyr-1)){
      for (s in 1:nsite){
        # Abundance of nonbreeders
        N[2, t+1, s] ~ dbin(mn.phiA[t, s]*(1-mn.psiAB[t, s]), NF[t, s]) # Nonbreeders to nonbreeders
        N[3, t+1, s] ~ dbin(mn.phiB[t, s]*mn.psiBA[t, s], NB[t, s]) # Breeders to nonbreeders
        N[4, t+1, s] <- 0 # Immigrants to nonbreeders
        # Abundance of breeders
        N[6, t+1, s] ~ dbin(mn.phiA[t, s]*mn.psiAB[t, s], NF[t, s]) # Nonbreeder to breeder
        N[7, t+1, s] ~ dbin(mn.phiB[t, s]*(1-mn.psiBA[t, s]), NB[t, s]) # Breeder to breeder
        N[8, t+1, s] <- 0 # immigrants to breeders
      } # s
      # FY to nonbreeder
      N[1, t+1, 1] ~ dpois( # born wild to nonbreeder
          NB[t,1]*(mn.phiFY[t, 1]*(1-mn.psiFYB[t, 1])*
                 (mn.prod[1, t, 1]/2)*(1-phacked.lh[t])* # not hacked
                   (1-ptreat[t,1]) ) + # not treated
          NB[t,1]*(mn.phiFY[t, 1]*(1-mn.psiFYB[t, 1])*
                   (mn.prod[2, t, 1]/2)*(1-phacked.lh[t])* # not hacked
                   ptreat[t,1]) ) # treated
      N[1, t+1, 2] ~ dpois( # born hacked to nonbreeder
        # Prod of LH to PC
          NB[t,1]*(mn.phiFY[t, 2]*(1-mn.psiFYB[t, 2])*
                 (mn.prod[1, t, 1]/2)*phacked.pc[t]* # hacked 
                   (1-ptreat[t,1])) + # not treated
          NB[t,1]*(mn.phiFY[t, 2]*(1-mn.psiFYB[t, 2])*
                   (mn.prod[2, t, 1]/2)*phacked.pc[t]* # hacked
                    ptreat[t,1])  + # treated
        # Prod of PC birds
          NB[t,2]*(mn.phiFY[t, 2]*(1-mn.psiFYB[t, 2])*
                       (mn.prod[1, t, 2]/2)*
                       (1-ptreat[t,2])) + # not treated
          NB[t,2]*(mn.phiFY[t, 2]*(1-mn.psiFYB[t, 2])*
                       (mn.prod[2, t, 2]/2)*
                      ptreat[t,2])  # treated
      )
      
      
      
      N[1, t+1, 1] ~ dpois( # born wild to nonbreeder
        NFledged[t,1]*(mn.phiFY[t, 1]*(1-mn.psiFYB[t, 1])*
                   (1-phacked.lh[t]))
        )
      N[1, t+1, 2] ~ dpois( # born hacked to nonbreeder
        # Prod of LH to PC
        NFledged[t,1]*(mn.phiFY[t, 2]*(1-mn.psiFYB[t, 2])*
                   phacked.pc[t]) + 
        NFledged[t,2]*(mn.phiFY[t, 2]*(1-mn.psiFYB[t, 2])
      )
      
      # FY to breeder
      N[5, t+1, 1] ~ dpois( # born wild to breeder
          NB[t,1]* (mn.phiFY[t, 1]*mn.psiFYB[t, 1]*
                  (mn.prod[1, t, 1]/2)*(1-phacked.lh[t])*
                    (1-ptreat[t,1])) +
          NB[t,1]* (mn.phiFY[t, 1]*mn.psiFYB[t, 1]*
                    (mn.prod[2, t, 1]/2)*(1-phacked.lh[t]))*
                    ptreat[t,1] )
      N[5, t+1, 2] ~ dpois( # born hacked to breeder
          # Prod of LH to PC
          NB[t,1]* (mn.phiFY[t, 2]*mn.psiFYB[t, 2]*
                    (mn.prod[1, t, 1]/2)*phacked.pc[t]*
                    (1-ptreat[t,1])) +
          NB[t,2]* (mn.phiFY[t, 2]*mn.psiFYB[t, 2]*
                      (mn.prod[2, t, 2]/2)*phacked.pc[t])*
                  ptreat[t,2] +  
          # Prod of PC
          NB[t,2]* (mn.phiFY[t, 2]*mn.psiFYB[t, 2]*
                  (mn.prod[1, t, 2]/2)*phacked.pc[t]*
                    (1-ptreat[t,2])) +
          NB[t,2]* (mn.phiFY[t, 2]*mn.psiFYB[t, 2]*
                    (mn.prod[2, t, 2]/2)*phacked.pc[t])*
                    ptreat[t,2] )
    } # t

    # Count likelihoods, state-space model, and observation process
    for (t in 1:nyr){
      for (s in 1:nsite){
        Nfledged[t, s] ~ dpois( NB[t,s]*(mn.prod[1, t, s]/2)*(1-ptreat[t,s]) +
                                NB[t,s]*(mn.prod[2, t, s]/2)*ptreat[t,s] )
        #NFY[t, s] <- sum(N[find[1:2], t, s]) # number of adult nonbreeder females
        NF[t, s] <- sum(N[1:4, t, s]) # number of adult nonbreeder females
        NB[t, s] <- sum(N[5:8, t, s]) # number of adult breeder females
        NAD[t, s] <- NF[t, s] + NB[t, s]  # number of adults
        Ntot[t, s] <- sum(N[1:8, t, s]) # total number of females
        countsAdults[t, s] ~ dpois(lamAD[t, s]) # adult females, includes nonbreeders and breeders but its a temp approx. to help with inits
        log(lamAD[t, s]) <- log(NAD[t, s]) #+ deltas[1]*effort2[t, s] + deltas[2]*effort2[t, s]^2
      }# s
    } # t
    
    ###################
    # Assess GOF of the state-space models for counts
    # Step 1: Compute statistic for observed data
    # Step 2: Use discrepancy measure: mean absolute error
    # Step 3: Use test statistic: number of turns
    ###################
    for (t in 1:nyr){
      for (s in 1:nsite){
        c.expAD[t, s] <- lamAD[t, s]  # expected counts adult breeder
        c.obsAD[t, s] <- countsAdults[t, s]
        c.repAD[t, s] ~ dpois( lamAD[t, s] ) # simulated counts
        dssm.obsAD[t, s] <- abs( ( (c.obsAD[t, s]) - (c.expAD[t, s]) ) / (c.obsAD[t, s]+0.001)  )
        dssm.repAD[t, s] <- abs( ( (c.repAD[t, s]) - (c.expAD[t, s]) ) / (c.repAD[t, s]+0.001) )
      }} # t
    dmape.obs[1] <- sum(dssm.obsAD[1:nyr, 1:nsite])
    dmape.rep[1] <- sum(dssm.repAD[1:nyr, 1:nsite])
    tvm.obs[1] <- sd(c.obsAD[1:nyr, 1:nsite])^2/mean(c.obsAD[1:nyr, 1:nsite])
    tvm.rep[1] <- sd(c.repAD[1:nyr, 1:nsite])^2/mean(c.repAD[1:nyr, 1:nsite])
    # # Test statistic for number of turns
    for (s in 1:nsite){
      for (t in 1:(nyr-2)){
        tt1.obsAD[t, s] <- step(c.obsAD[t+2, s] - c.obsAD[t+1, s])
        tt2.obsAD[t, s] <- step(c.obsAD[t+1, s] - c.obsAD[t, s])
        tt3.obsAD[t, s] <- equals(tt1.obsAD[t, s] + tt2.obsAD[t, s], 1)
      }} # t
    tturn.obs[1] <- sum(tt3.obsAD[1:(nyr-2), 1:nsite])
    tturn.obs[2] <- sum(tt3.obsFY[1:(nyr-2), 1:nsite])
    for (s in 1:nsite){
      for (t in 1:(nyr-2)){
        tt1.repAD[t, s] <- step(c.repAD[t+2, s] - c.repAD[t+1, s])
        tt2.repAD[t, s] <- step(c.repAD[t+1, s] - c.repAD[t, s])
        tt3.repAD[t, s] <- equals(tt1.repAD[t, s] + tt2.repAD[t, s], 1)
      }} # t
    tturn.rep[1] <- sum(tt3.repAD[1:(nyr-2), 1:nsite])


    # Probability a FY was hacked
    # Probability of an individual at a site
    # has a history of being hacked
    lmu.hacked.lh <- logit(mu.hacked.lh)
    mu.hacked.lh ~ dbeta(1,1)
    lmu.hacked.pc <- logit(mu.hacked.pc)
    mu.hacked.pc ~ dbeta(1,1)
    
    for (h in 1:nhacked.lh){
      hacked.lh[h] ~ dbern( phacked.lh[ lh.year[h] ] )
    } # h
    for (hh in 1:nhacked.pc){
      hacked.pc[hh] ~ dbern( phacked.pc[ pc.year[hh] ] )

    } # h
    for (t in 1:nyr){
      hacked.lh[t] ~ dbin( phacked.lh[ t ], 
                           Nfledged[t, s] )
      hacked.pc[t] ~ dbin( phacked.pc[ t ], 
                            Nfledged[t, s] )
      # probability that a fledgling at LH
      # is translocated from the site to
      # Aniana Vargas or Punta Cana
      logit(phacked.lh[t]) <- lmu.hacked.lh #+ eta[9, s, t]
      # probability that a fledgling at LH
      # is translocated to PC
      logit(phacked.pc[t]) <- lmu.hacked.pc #+ eta[9, s, t]
    } # t

    # ###############################
    # Likelihood for survival
    # ###############################
    
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
          lmn.phiFY[t, site[i,t]] #+
          #betas[1]*c(-1, 1)[ hacked[i]+1 ]
        logit(phiA[i,t]) <- 
          lmn.phiA[t, site[i,t]] #+
          #betas[2]*c(-1, 1)[ hacked[i]+1 ]
        logit(phiB[i,t]) <- 
          lmn.phiB[t, site[i,t]] #+
          #betas[3]*c(-1, 1)[ hacked[i]+1 ]
        #Recruitment
        logit(psiFYB[i,t]) <- 
          lmn.psiFYB[t, site[i,t]] #+
          #betas[4]*c(-1, 1)[ hacked[i]+1 ]
        
        logit(psiAB[i,t]) <- 
          lmn.psiAB[t, site[i,t]] #+
          #betas[5]*c(-1, 1)[ hacked[i]+1 ]
        logit(psiBA[i,t]) <- 
          lmn.psiBA[t, site[i,t]]
        #Re-encounter
        logit(pA[i,t]) <- 
          lmus[7, site[i,t]] #+  
          #deltas[5]*effort2[t, site[i,t]] + 
          #deltas[6]*effort2[t, site[i,t]]^2# resight of nonbreeders
        logit(pB[i,t]) <- 
          lmus[8, site[i,t]] #+
          #eta[6, site[i,t],t] #+
          #betas[7]*c(-1, 1)[ hacked[i]+1 ] + 
          #deltas[7]*effort2[t, site[i,t]] + 
          #deltas[8]*effort2[t, site[i,t]]^2# resight of breeders
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
  library ('MCMCvis')
  #source("/bsuscratch/brianrolek/riha_ipm/MCMCvis.R")
  #source("/bsuscratch/brianrolek/riha_ipm/check_N_inits.R")
  
  params <- c(
    # fecundity
    "lmu.prod", "gamma", "rr", "mn.prod", 
    "mu.treat", "ptreat",
    # survival 
    "mus", "lmus", "betas", "deltas", "alphas",
    "mu.hacked", "lmu.hacked",
    # abundance
    #"lNAD", "lN", "lamAD",
    "NB", "NF", #"NFY", 
    "N", "NAD",
    #"r",
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
      betas = apply(outp$betas, 1, mean),
      deltas = c(0,0,0,0, apply(outp$deltas, 1, mean)[5:8]),
      alphas = apply(outp$alphas, c(1,2), mean),
      sds = runif(constl$p),
      Ustar = diag(constl$p),
      z.score = array(runif(constl$p*2*13), dim=c(constl$p,2,13)),
      mu.hacked = runif(2),
      # counts
      #r = mean(outp$r),
      #N = get_inits2(hpc=TRUE) , # sample from inits of chains that worked
      # immigration
      omega = runif(1, 0.01, 0.1)
    )}
  inits2 <- inits.func()
  
  mod <- nimbleModel(code, 
                     constants = constl, 
                     data = datl, 
                     inits = inits2, 
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
  #nc=1; nt=1; ni=10000; nb=5000
  nc=1; nt=1; ni=200; nb=0
  
  post <- runMCMC(chmc,
                   niter = ni,
                   nburnin = nb,
                   nchains = nc,
                   thin = nt,
                   samplesAsCodaMCMC = T,
                   setSeed = info$seed,
                   inits = inits2)
  
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

flnm2 <- paste0("/bsuscratch/brianrolek/riha_ipm/outputs/ipm-poistreat",
                substr(Sys.time(), 1, 10), ".rdata")
save(post, mycode, 
     file = flnm2)
# 

# save(post, mycode,
#   file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm-simp-2025-05-02.rdata")
