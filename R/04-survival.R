library ('nimble')
library ('MCMCvis')
library ('coda')
library('parallel')

setwd("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\GitHub\\Ridgways-Hawk-IPM\\")
load("data\\data.rdata")
set.seed(131517)

constl$p2 <- constl$p <- 8
par_info[[1]]$inits$sds2 <- par_info[[1]]$inits$sds <- par_info[[1]]$inits$sds[1:constl$p]
par_info[[2]]$inits$sds2 <- par_info[[2]]$inits$sds <- par_info[[2]]$inits$sds[1:constl$p]
par_info[[3]]$inits$sds2 <- par_info[[3]]$inits$sds <- par_info[[3]]$inits$sds[1:constl$p]
par_info[[4]]$inits$sds2 <- par_info[[4]]$inits$sds <- par_info[[4]]$inits$sds[1:constl$p]

par_info[[1]]$inits$Ustar2 <- par_info[[1]]$inits$Ustar <- par_info[[1]]$inits$Ustar[1:constl$p, 1:constl$p]
par_info[[2]]$inits$Ustar2 <- par_info[[2]]$inits$Ustar <- par_info[[2]]$inits$Ustar[1:constl$p, 1:constl$p]
par_info[[3]]$inits$Ustar2 <- par_info[[3]]$inits$Ustar <- par_info[[3]]$inits$Ustar[1:constl$p, 1:constl$p]
par_info[[4]]$inits$Ustar2 <- par_info[[4]]$inits$Ustar <- par_info[[4]]$inits$Ustar[1:constl$p, 1:constl$p]

# runif(constl$p) ->
# par_info[[1]]$inits$mus -> par_info[[2]]$inits$mus -> par_info[[3]]$inits$mus ->par_info[[4]]$inits$mus
par_info[[1]]$inits$betas <- par_info[[2]]$inits$betas <- par_info[[3]]$inits$betas <- par_info[[4]]$inits$betas <- rep(0, constl$p)

datl$mu.zeroes <- rep(0, constl$p)
datl$mu.zeroes2 <- rep(0, constl$p2)
####################################################
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
    ################################
    # Likelihood for survival
    ################################ 
    # survival, recruitment, and detection can be correlated
    
    # priors for means
    for (j in 1:p){ # coefficient
      sds[j] ~ dexp(1)
      betas[j] ~ dunif(-10,10)
      for (s in 1:nsite){
      lmus[j,s] <- logit(mus[j,s])
      mus[j,s] ~ dbeta(1,1)
    }} # s, p
    
    # estimated using the multivariate normal distribution
    R[1:p,1:p] <- t(Ustar[1:p,1:p]) %*% Ustar[1:p,1:p] # calculate rhos, correlation coefficients
    Ustar[1:p,1:p] ~ dlkj_corr_cholesky(eta=1.1, p=p) # Ustar is the Cholesky decomposition of the correlation matrix
    U[1:p,1:p] <- uppertri_mult_diag(Ustar[1:p, 1:p], sds[1:p])
    # multivariate normal for temporal variance
    for (t in 1:(nyr-1)){
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
    for (t in 1:nyr){ # survival params only have nyr-1, no problem to simulate from however
      for (s in 1:nsite){
        eta[1:p2,s,t] ~ dmnorm(mu.zeroes2[1:p2],
                               cholesky = U2[1:p2, 1:p2], prec_param = 0)
      } } # s t 
    
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
      for (t in 1:(nyr-1)){
        #Survival
        logit(phiFY[i,t]) <- eps[1,t] + eta[1, site[i,t], t] + 
          lmus[1, site[i,t]] + betas[1]*hacked[i]  # first year
        logit(phiA[i,t]) <- eps[2,t] + eta[2, site[i,t], t] + 
          lmus[2, site[i,t]] +  betas[2]*hacked[i] # nonbreeder
        logit(phiB[i,t]) <- eps[3,t]  + eta[3, site[i,t], t] + 
          lmus[3, site[i,t]] + betas[3]*hacked[i] # breeder
        #Recruitment
        logit(psiFYB[i,t]) <- eps[4,t] + eta[4, site[i,t], t] + 
          lmus[4, site[i,t]] + betas[4]*hacked[i] # first year to breeder
        logit(psiAB[i,t]) <- eps[5,t] + eta[5, site[i,t], t] + 
          lmus[5, site[i,t]] + betas[5]*hacked[i] # nonbreeder to breeder
        logit(psiBA[i,t]) <- eps[6,t] + eta[6, site[i,t], t] + 
          lmus[6, site[i,t]] + betas[6]*hacked[i] # breeder to nonbreeder
        #Re-encounter
        logit(pA[i,t]) <- eps[7,t] + eta[7, site[i,t], t] + 
          lmus[7, site[i,t]] + betas[7]*hacked[i] # resight of nonbreeders
        logit(pB[i,t]) <- eps[8,t] + eta[8, site[i,t], t] + 
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
  } ) #model

run_surv <- function(info, datl, constl, code){
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
  
  
  params <- c( "mus", "lmus", "betas", "sds", "sds2",
               "Ustar", "U", "R", 
               "Ustar2", "U2", "R2",
               'mn.phiFY', 'mn.phiA', 'mn.phiB',
               'mn.psiFYB', 'mn.psiAB', 'mn.psiBA',
               'mn.pA', 'mn.pB', 
               "eps", "eta")
  
  #nc=1; nt=1; ni=10000; nb=0
  nc=1; nt=100; ni=200000; nb=100000
  
  # can't use enableWAIC=TRUE 
  mod <- nimbleModel(code, calculate=T, 
                     constants = constl, 
                     data = datl, 
                     inits = info$inits)
  cmod <- compileNimble(mod, showCompilerOutput = TRUE )
  confhmc <- configureMCMC(mod)
  confhmc$setMonitors(params)
  hmc <- buildMCMC(confhmc)
  chmc <- compileNimble(hmc, project=mod, 
                        resetFunctions = TRUE,
                        showCompilerOutput = TRUE )
  
  post <- runMCMC(chmc,
                  niter = ni, 
                  nburnin = nb,
                  nchains = nc,
                  thin = nt,
                  samplesAsCodaMCMC = T,
                  setSeed = info$seed)
  
  return(post)
}

this_cluster <- makeCluster(4)
post <- parLapply(cl = this_cluster, 
                  X = par_info, 
                  fun = run_surv, 
                  datl = datl, 
                  constl = constl, 
                  code = mycode)
stopCluster(this_cluster)

save(post, run_surv,  
     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\surv_nimble.Rdata")
