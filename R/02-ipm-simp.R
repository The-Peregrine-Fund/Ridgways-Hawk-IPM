## ---- ipm-simp --------
library('nimble')
library('parallel')
library ('coda')
load("/bsuscratch/brianrolek/riha_ipm/data.rdata")

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
    # Survival, recruitment, and detection can be correlated
    for (k in 1:8){ 
      betas[k] ~ dunif(-20, 20)  # prior for coefficients
    } # k
    for (j in 1:8){ 
      for (s in 1:nsite){
        lmus[j,s] <- logit(mus[j,s])
        mus[j,s] ~ dbeta(1,1) # prior for means
      }}     # m population #s sex #h hacked 
    
    # # Temporal random effects and correlations among all sites, synchrony
    # for (j in 1:p){ sds[j] ~ dexp(1) }# prior for temporal variation estimated using the multivariate normal distribution
    # R[1:p,1:p] <- t(Ustar[1:p,1:p]) %*% Ustar[1:p,1:p] # calculate rhos, correlation coefficients
    # Ustar[1:p,1:p] ~ dlkj_corr_cholesky(eta=1.1, p=p) # Ustar is the Cholesky decomposition of the correlation matrix
    # U[1:p,1:p] <- uppertri_mult_diag(Ustar[1:p, 1:p], sds[1:p])
    # # Multivariate normal for temporal variance
    # for (t in 1:nyr){ 
    #   eps[1:p,t] ~ dmnorm(mu.zeroes[1:p],
    #                       cholesky = U[1:p, 1:p], prec_param = 0)
    # }
    
    # Temporal random effects and correlations between sites
    for (jj in 1:p2){ sds2[jj] ~ dexp(1) }# prior for temporal variation estimated using the multivariate normal distribution
    R2[1:p2,1:p2] <- t(Ustar2[1:p2,1:p2]) %*% Ustar2[1:p2,1:p2] # calculate rhos, correlation coefficients
    Ustar2[1:p2,1:p2] ~ dlkj_corr_cholesky(eta=1.1, p=p2) # Ustar is the Cholesky decomposition of the correlation matrix
    U2[1:p2,1:p2] <- uppertri_mult_diag(Ustar2[1:p2, 1:p2], sds2[1:p2])
    # multivariate normal for temporal and site variance
    for (t in 1:nyr){ 
      for (s in 1:nsite){
        eta[1:p2,s,t] ~ dmnorm(mu.zeroes2[1:p2],
                               cholesky = U2[1:p2, 1:p2], prec_param = 0)
      } } # s t 

    ###############################
    # Likelihood for productivity
    ###############################
    # Priors for productivity
    for (s in 1:nsite){
      lmu.f[s] ~ dnorm(0, sd=5)
    } # s
    gamma ~ dunif(-20, 20)
    rr ~ dexp(0.05)
    
    # Productivity likelihood      
    for (k in 1:npairsobs){
      f[k] ~ dnegbin(ppp[k], rr)
      ppp[k] <- rr/(rr+mu.f[k])
      log(mu.f[k]) <- lmu.f[site.pair[k]] +  
        gamma*treat.pair[k] + 
        #eps[9, year.pair[k] ] + 
        eta[9, site.pair[k], year.pair[k] ] 
    } # k
    # Derive yearly brood size for population model
    # Need to reorder because nimble doesn't 
    # handle nonconsecutive indices
    # yrind.pair is a matrix of indices for each site
    for (t in 1:nyr){
      for (s in 1:nsite){
        for (xxx in 1:pair.end[t,s]){
          fecmat[t,s,xxx] <- mu.f[ yrind.pair[xxx,t,s] ]
        } # xxx
        mn.f[t,s] <- mean( fecmat[t,s,1:pair.end[t,s]] )
      }} # s t
    
    # GOF for number of fledglings
    for (k in 1:npairsobs){
      f.obs[k] <- f[k] # observed counts
      f.exp[k] <- mu.f[k] # expected counts adult breeder
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
# Omitted. See IPM or PVA.
    
    ################################
    # Likelihood for survival
    ################################ 
    # Calculate yearly averages for sites for integration
    for (t in 1:nyr){
      for (s in 1:nsite){
        for (xxxx in 1:surv.end[t,s]){
          # Need to reorder because nimble doesn't 
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
        logit(phiFY[i,t]) <- eta[1, site[i,t],t] + # eps[1,t] + 
          lmus[1, site[i,t]] + betas[1]*hacked[i]  # first year
        logit(phiA[i,t]) <- eta[2, site[i,t],t] +# eps[2,t] + 
          lmus[2, site[i,t]] +  betas[2]*hacked[i] # nonbreeder
        logit(phiB[i,t]) <- eta[3, site[i,t],t] +# eps[3,t] + 
          lmus[3, site[i,t]] + betas[3]*hacked[i] # breeder
        #Recruitment
        logit(psiFYB[i,t]) <- eta[4, site[i,t],t] +# eps[4,t] + 
          lmus[4, site[i,t]] + betas[4]*hacked[i] # first year to breeder
        logit(psiAB[i,t]) <- eta[5, site[i,t],t] +# eps[5,t] + 
          lmus[5, site[i,t]] + betas[5]*hacked[i] # nonbreeder to breeder
        logit(psiBA[i,t]) <- #eta[6, site[i,t],t] + eps[6,t] + 
          lmus[6, site[i,t]] #+ betas[6]*hacked[i] # breeder to nonbreeder
        #Re-encounter
        logit(pA[i,t]) <- eta[7, site[i,t],t] +# eps[7,t] + 
          lmus[7, site[i,t]] + betas[7]*hacked[i] # resight of nonbreeders
        logit(pB[i,t]) <- eta[8, site[i,t],t] +# eps[8,t] + 
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

#**********************
#* Function to run model in NIMBLE
#**********************
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
    #"lambda",
    # fecundity
    "lmu.f", "gamma", "rr", "mn.f", 
    # survival 
    "mus", "lmus", "betas",
    # abundance
    # "NB", "NF", "NFY", "N",
    # "r",
    # "N", "Ntot",
    # error terms
    #"eps", "sds", "Ustar", "U", "R",
    "eta", "sds2", "Ustar2", "U2", "R2",
    # yearly summaries
    'mn.phiFY', 'mn.phiA', 'mn.phiB',
    'mn.psiFYB', 'mn.psiAB', 'mn.psiBA',
    'mn.pA', 'mn.pB',
    # goodness of fit
    "f.dmape.obs", "f.dmape.rep",
    "f.tvm.obs", "f.tvm.rep"
    # "dmape.obs", "dmape.rep",
    # "tvm.obs", "tvm.rep",
    # "tturn.obs", "tturn.rep"
  )
  
  #n.chains=1; n.thin=200; n.iter=500000; n.burnin=300000
  #n.chains=1; n.thin=50; n.iter=100000; n.burnin=50000
  #n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
  n.chains=1; n.thin=200; n.iter=500000; n.burnin=300000
  
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
this_cluster <- makeCluster(10)
post <- parLapply(cl = this_cluster, 
                  X = par_info, 
                  fun = run_ipm, 
                  datl = datl, 
                  constl = constl, 
                  code = mycode)
stopCluster(this_cluster)

save(post, mycode,
     file="/bsuscratch/brianrolek/riha_ipm/outputs/ipm_simp.rdata")

# save(post, mycode,
#      file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_statespace.rdata")