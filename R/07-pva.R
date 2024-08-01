## ---- ipm1 --------
library('nimble')
library('parallel')
library ('coda')
load("/bsuscratch/brianrolek/riha_ipm/data.rdata")
#load("data/data.rdata")
constl$K <- 100

hacked.counts <- array(0, dim=c(6, 113, 3))
for (sc in 1:6){
hacked.counts[sc,1:13,1:3] <- constl$hacked.counts 
hacked.counts[sc,14:113,1:3] <- c(0, 0, 0, 1, 1, 1)[sc]*cbind(rep(-10, 100), rep(10, 100), rep(0, 100) )
}

constl$hacked.counts <- hacked.counts
datl$constraint_data <- rbind(datl$constraint_data, array(1, dim=c(constl$K,2)) )
constl$treat.nest2 <- c(1, 0, 1, 1, 0, 1) # treated or not
constl$hacked2 <- c(0, 0, 0, 1, 1, 1) # scenario- hacked or not

# index to help add informative priors for transition NB to B
#constl$j2 <- c(1,2,3,4,5,7,8)
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
      betas[k] ~ dunif(-20, 20)  # prior for coefficients
    } # k
    # for (j in 1:7){
    #   for (s in 1:nsite){
    #     lmus[j2[j],s] <- logit(mus[j2[j],s])
    #     mus[j2[j],s] ~ dbeta(1,1) # prior for means
    #   }}     # m population #s sex #h hacked
    # mus[6,1] ~ dbeta(1,15) # informative prior for transition NB to B
    # mus[6,2] ~ dbeta(1,1)
    
    for (j in 1:8){
      for (s in 1:nsite){
        lmus[j,s] <- logit(mus[j,s])
        mus[j,s] ~ dbeta(1,1) # prior for means
      }}     # m population #s sex #h hacked
    
    # Temporal random effects and correlations among all sites
    for (j in 1:p){ sds[j] ~ dexp(1) }# prior for temporal variation
    # estimated using the multivariate normal distribution
    R[1:p,1:p] <- t(Ustar[1:p,1:p]) %*% Ustar[1:p,1:p] # calculate rhos, correlation coefficients
    Ustar[1:p,1:p] ~ dlkj_corr_cholesky(eta=1.1, p=p) # Ustar is the Cholesky decomposition of the correlation matrix
    U[1:p,1:p] <- uppertri_mult_diag(Ustar[1:p, 1:p], sds[1:p])
    # multivariate normal for temporal variance
    for (t in 1:(nyr+K) ){ # survival params only have nyr-1, no problem to simulate from however
      eps[1:p,t] ~ dmnorm(mu.zeroes[1:p],
                          cholesky = U[1:p, 1:p], prec_param = 0)
    }
    
    # Temporal random effects and correlations between sites
    for (jj in 1:p2){ sds2[jj] ~ dexp(1) }# prior for temporal variation
    # estimated using the multivariate normal distribution
    R2[1:p2,1:p2] <- t(Ustar2[1:p2,1:p2]) %*% Ustar2[1:p2,1:p2] # calculate rhos, correlation coefficients
    Ustar2[1:p2,1:p2] ~ dlkj_corr_cholesky(eta=1.1, p=p2) # Ustar is the Cholesky decomposition of the correlation matrix
    U2[1:p2,1:p2] <- uppertri_mult_diag(Ustar2[1:p2, 1:p2], sds2[1:p2])
    # multivariate normal for temporal variance
    for (t in 1:(nyr+K) ){ # survival params only have nyr-1, no problem to simulate from however
      for (s in 1:nsite){
        eta[1:p2,s,t] ~ dmnorm(mu.zeroes2[1:p2],
                               cholesky = U2[1:p2, 1:p2], prec_param = 0)
      } } # s t 
    #######################
    # Derived params
    #######################
    for(sc in 1:6){
      
    for (s in 1:nsite){
      for (t in 1:(nyr-1+K)){
        lambda[sc, t, s] <-  Ntot[sc, t+1, s]/(Ntot[sc, t, s])
        loglambda[sc, t, s] <- log(lambda[sc, t, s])
      }} #t
    
    for (s in 1:nsite){
      for (t in 1:K){
        extinct[sc, t, s] <- equals(Ntot[sc, nyr+t, s], 0)
        extinctAD[sc, t, s] <- equals(NAD[sc, nyr+t, s], 0)
        extinctB[sc, t, s] <- equals(NB[sc, nyr+t, s], 0)
      }}
    } # sc
    
    ###############################
    # Likelihood for fecundity
    ###############################
    # Priors for number of fledglings
    # and nest success
    for (s in 1:nsite){
      lmu.f[s] ~ dnorm(0, sd=5)
    } # s
    gamma ~ dunif(-20, 20)
    rr ~ dexp(0.05)
    
    # Fecundity       
    for (k in 1:nnest){
      f[k] ~ dnegbin(ppp[k], rr)
      ppp[k] <- rr/(rr+mu.f[k])
      log(mu.f[k]) <- lmu.f[site.nest[k]] +  
        gamma*treat.nest[k] + 
        eps[9, year.nest[k] ] + 
        eta[9, site.nest[k], year.nest[k] ] 
    } # k
    # derive yearly brood size for population model
    for (t in 1:nyr){
      for (s in 1:nsite){
        for (xxx in 1:nest.end[t,s]){
          fecmat[t,s,xxx] <- mu.f[ yrind.nest[xxx,t,s] ]
        } # xxx
        mn.f[1,t,s] <- mean( fecmat[t,s,1:nest.end[t,s]] )
        for(sc in 2:6){
          mn.f[sc,t,s] <- mn.f[1,t,s]
      }}} # s t
    # Future fecundity

      for(sc in 1:6){
        for (t in (nyr+1):(nyr+K) ){
          for (s in 1:nsite){
        # calc percent of nests treated so it cannot exceed 100% 
        perc.treat[sc, t, s] <- step( 1- 10/NB[sc, t, s] ) * 10/NB[sc, t, s] +
                                (1-step( 1- 10/NB[sc, t, s] ))
        log(mn.f[sc,t,s]) <- lmu.f[s] +  
                              gamma*treat.nest2[sc]*perc.treat[sc, t, s] + 
                              eps[9, t] + 
                              eta[9, s, t] 
  }} } # s t sc
    
    ################################
    # Likelihood for counts
    ################################
    # Abundance for year=1
    for(sc in 1:1){
      for (v in 1:7){ 
      for (s in 1:nsite){
        # subtract one because to allow dcat to include zero
        N[sc, v, 1, s] <- N2[sc, v, 1, s] - 1 
        N2[sc, v, 1, s] ~ dcat(pPrior[1:s.end[s], s])
        # Assign estimates for years with data to 
        # all scenarios
        for (sc2 in 2:6){
          N[sc2, v, 1, s] <- N[1, v, 1, s]
        }
      }} # s t

    # Abundance for years > 1
    for (t in 1:(nyr-1)){
      for (s in 1:nsite){
        # Number of wild born juvs
        N[sc, 1, t+1, s] ~ dpois( (NFY[sc, t, s]*mn.phiFY[sc, t, s]*mn.psiFYB[sc, t, s] + # first year breeders
                                 NF[sc, t, s]*mn.phiA[sc, t, s]*mn.psiAB[sc, t, s] + # nonbreeders to breeders
                                 NB[sc, t, s]*mn.phiB[sc, t, s]*(1-mn.psiBA[sc, t, s])) # breeders remaining
                              *mn.f[sc, t+1, s] ) # end Poisson
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
        # Assign estimates for years with data to 
        # all scenarios
        for (sc2 in 2:6){
          N[sc2, 1, t+1, s] <- N[sc, 1, t+1, s]
          N[sc2, 2, t+1, s] <- N[sc, 2, t+1, s]
          N[sc2, 3, t+1, s] <- N[sc, 3, t+1, s]
          N[sc2, 4, t+1, s] <- N[sc, 4, t+1, s]
          N[sc2, 5, t+1, s] <- N[sc, 5, t+1, s]
          N[sc2, 6, t+1, s] <- N[sc, 6, t+1, s]
          N[sc2, 7, t+1, s] <- N[sc, 7, t+1, s]
        } # sc2
        }} # s t
    } # sc
    
  # Abundance for years after data collection
  # PVA portion   
    for (sc in 1:6){
      for (t in nyr:(nyr+K-1)){
        for (s in 1:nsite){
          # Number of wild born juvs
          N[sc, 1, t+1, s] ~ dpois( (NFY[sc, t, s]*mn.phiFY[sc, t, s]*mn.psiFYB[sc, t, s] + # first year breeders
                                       NF[sc, t, s]*mn.phiA[sc, t, s]*mn.psiAB[sc, t, s] + # nonbreeders to breeders
                                       NB[sc, t, s]*mn.phiB[sc, t, s]*(1-mn.psiBA[sc, t, s])) # breeders remaining
                                    *mn.f[sc, t+1, s] ) # end Poisson
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
    
    # Observation process 
  for(sc in 1:6){
    for (t in 1:(nyr+K)){
      for (s in 1:nsite){
        # constrain N1+hacked.counts to be >=0
        NFY[sc, t, s] <- step( (N[sc, 1, t, s] + hacked.counts[sc, t, s]) ) * (N[sc, 1, t, s] + hacked.counts[sc, t, s]) 
        #constraint_data[t, s] ~ dconstraint( (N[sc, 1, t, s] + hacked.counts[sc, t, s]) >= 0 ) # Transfers translocated first-year females
        #NFY[sc, t, s] <- N[sc, 1, t, s] + hacked.counts[sc, t, s] # Transfers translocated first-year females
        NF[sc, t, s] <- sum(N[sc, 2:4, t, s]) # number of adult nonbreeder females
        NB[sc, t, s] <- sum(N[sc, 5:7, t, s]) # number of adult breeder females
        NAD[sc, t, s] <- NF[sc, t, s] + NB[sc, t, s] # number of adults
        Ntot[sc, t, s] <- sum(N[sc, 1:7, t, s]) # total number of females
      }# s
    }  }# t sc
    
  for(sc in 1:6){
    for (t in 1:nyr){
      for (s in 1:nsite){
        countsAdults[t, s] ~ dpois(NAD[sc, t, s]) # adult females 
      } # s
      # First-years at different sites have different distributions
      # for better model fit
      counts[1, t, 1] ~ dpois(N[sc, 1, t, 1]) # doesn't have any zeroes so poisson
      counts[1, t, 2] ~ dnegbin(pp[sc, t], r) # first year females, includes translocated/hacked
      pp[sc, t] <- r/(r+N[sc, 1, t, 2])
      }} # t
    r ~ dexp(0.05)
    
    ################################
    # Likelihood for survival
    ################################ 
    # Calculate an averages for sites
    # each year for integration
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
        for (sc in 1:6){
        mn.phiFY[sc, t, s] <- mean( phiFY2[ t, s, 1:surv.end[t,s] ] ) 
        mn.phiA[sc, t, s] <- mean( phiA2[ t, s, 1:surv.end[t,s] ] )
        mn.phiB[sc, t, s] <- mean( phiB2[ t, s, 1:surv.end[t,s] ] )
        mn.psiFYB[sc, t, s] <- mean( psiFYB2[ t, s, 1:surv.end[t,s] ] )
        mn.psiAB[sc, t, s] <- mean( psiAB2[ t, s, 1:surv.end[t,s] ] )
        mn.psiBA[sc, t, s] <- mean( psiBA2[ t, s, 1:surv.end[t,s] ] )
        mn.pA[sc, t, s] <- mean( pA2[ t, s, 1:surv.end[t,s] ] )
        mn.pB[sc, t, s] <- mean( pB2[ t, s, 1:surv.end[t,s] ] )
      }}}
    
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
        logit(psiBA[i,t]) <- #eta[6, site[i,t],t] + eps[6,t] + 
          lmus[6, site[i,t]] #+ betas[6]*hacked[i] # breeder to nonbreeder
        #Re-encounter
        logit(pA[i,t]) <- eta[7, site[i,t],t] + eps[7,t] + 
          lmus[7, site[i,t]] + betas[7]*hacked[i] # resight of nonbreeders
        logit(pB[i,t]) <- eta[8, site[i,t],t] + eps[8,t] + 
          lmus[8, site[i,t]] + betas[8]*hacked[i] # resight of breeders
      }#t
    }#i
    
    # Future survival, recruitment, and detection
    for(sc in 1:6){ 
    for (t in (nyr+1):(nyr+K)){
        for (s in 1:nsite){
          # calculate percent hacked
          # assigns a one (1.0 or 100%) if abundance is less than 10
          p.hackedFY[sc, t, s] <- step( 1- 10/NFY[sc, t, s] ) * 10/NFY[sc, t, s] +
                                      (1-step( 1- 10/NFY[sc, t, s] )) 
          p.hackedA[sc, t, s] <- step( 1- 10/NF[sc, t, s] ) * 10/NF[sc, t, s] +
                                      (1-step( 1- 10/NF[sc, t, s] ))
          p.hackedB[sc, t, s] <- step( 1- 10/NB[sc, t, s] ) * 10/NB[sc, t, s] +
                                      (1-step( 1- 10/NB[sc, t, s] ))
        # Hacking/translocation only affects birds
        # moved to Punta Cana (site 2) hence- equals(s, 2)
        #Survival
        logit(mn.phiFY[sc, t, s]) <- eta[1, s, t] + eps[1, t] + 
          lmus[1, s] + betas[1]*hacked2[sc]*p.hackedFY[sc, t, s]*equals(s, 2)  # first year
        logit(mn.phiA[sc, t, s]) <- eta[2, s, t] + eps[2, t] + 
          lmus[2, s] +  betas[2]*hacked2[sc]*p.hackedA[sc, t, s]*equals(s, 2) # nonbreeder
        logit(mn.phiB[sc, t, s]) <- eta[3, s, t] + eps[3, t] + 
          lmus[3, s] + betas[3]*hacked2[sc]*p.hackedB[sc, t, s]*equals(s, 2) # breeder
        #Recruitment
        logit(mn.psiFYB[sc, t, s]) <- eta[4, s, t] + eps[4, t] + 
          lmus[4, s] + betas[4]*hacked2[sc]*p.hackedFY[sc, t, s]*equals(s, 2) # first year to breeder
        logit(mn.psiAB[sc, t, s]) <- eta[5, s, t] + eps[5, t] + 
          lmus[5, s] + betas[5]*hacked2[sc]*p.hackedA[sc, t, s]*equals(s, 2) # nonbreeder to breeder
        logit(mn.psiBA[sc, t, s]) <- #eta[6, site[i,t],t] + eps[6,t] + 
          lmus[6, s] #+ betas[6]*hacked[i] # breeder to nonbreeder
        #Re-encounter
        logit(mn.pA[sc, t, s]) <- eta[7, s, t] + eps[7, t] + 
          lmus[7, s] + betas[7]*hacked2[sc]*p.hackedA[sc, t, s]*equals(s, 2) # resight of nonbreeders
        logit(mn.pB[sc, t, s]) <- eta[8, s, t] + eps[8, t] + 
          lmus[8, s] + betas[8]*hacked2[sc]*p.hackedB[sc, t, s]*equals(s, 2) # resight of breeders
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
    "lmu.f", "gamma", "rr", "mn.f", 
    # survival 
    "mus", "lmus", "betas",
    # abundance
    "NB", "NF", "NFY", "N", "NAD", "Ntot",
    "r",
    # error terms
    "eps", "sds", "Ustar", "U", "R",
    "eta", "sds2", "Ustar2", "U2", "R2",
    # yearly summaries
    'mn.phiFY', 'mn.phiA', 'mn.phiB',
    'mn.psiFYB', 'mn.psiAB', 'mn.psiBA',
    'mn.pA', 'mn.pB',
    # pva
    "extinct", "extinctAD", "extinctB"
    #"p.hackedFY", "p.hackedA", "p.hackedB"
  )
  
  #n.chains=1; n.thin=200; n.iter=500000; n.burnin=300000
  #n.chains=1; n.thin=50; n.iter=100000; n.burnin=50000
  n.chains=1; n.thin=1; n.iter=2000; n.burnin=0
  #n.chains=1; n.thin=200; n.iter=600000; n.burnin=400000
  
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


this_cluster <- makeCluster(5)
post <- parLapply(cl = this_cluster, 
                  X = par_info_pva[1:5], 
                  fun = run_ipm, 
                  datl = datl, 
                  constl = constl, 
                  code = mycode)
stopCluster(this_cluster)

save(post, mycode,
     file="/bsuscratch/brianrolek/riha_ipm/outputs/pva.rdata")

# save(post, mycode,
#      file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_statespace.rdata")
