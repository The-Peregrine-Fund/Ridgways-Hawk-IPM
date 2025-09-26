## ---- ipm --------
library('nimble')
library('parallel')
library ('coda')

# Screen posteriors to retain good inital values
screen <- function(x) {
  NAlist <- c()
  for (i in 1:length(x)){ # are there NAs or negative values in N
    NAlist[i] <- any (is.na(x[[i]][,1:311]) | x[[i]][,1:311]<0)
  } # end loop 
  # Subset chains to those with good initial values
  #print(!NAlist)
  outp <- MCMCpstr(x[!NAlist], type="chains")
  return(outp)
} # end function

# source("/bsuscratch/brianrolek/riha_ipm/MCMCvis.R")
# load("/bsuscratch/brianrolek/riha_ipm/data-dd.rdata")
# #load("/bsuscratch/brianrolek/riha_ipm/outputs/ipm-find-inits.rdata")
# load("/bsuscratch/brianrolek/riha_ipm/outputs/ipm_shortrun.rdata")
# outp <- screen(post)
# rm(list=c("mycode", "post"))
# cpus <- 5

library (MCMCvis)
load("data/data-dd.rdata")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm-find-inits.rdata")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_shortrun.rdata")
outp <- screen(post)
rm(list=c("mycode", "post"))
cpus <- 1
info <- par_info[[1]]

constl$p <- 10

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
      deltas[k] ~ dnorm(0, sd=10) # prior for survey effort coefficients
    } # k
    
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
    
    ###############################
    # Likelihood for productivity
    ###############################
    # Priors
    for (s in 1:nsite){
      lmu.prod[s] ~ dnorm(0, sd=5)
      lmu.treat[s] <- logit(mu.treat[s])
      mu.treat[s] ~ dbeta(1,1)
    } # s
    gamma ~ dnorm(0, sd=10)
    rr ~ dexp(0.05)
    sd.treat ~ dexp(1)
    
    # Productivity likelihood      
    for (k in 1:npairsobs){
      prod[k] ~ dnegbin(ppp[k], rr)
      ppp[k] <- rr/(rr+mu.prod[k])
      mu.prod[k] <- mn.prod[treat.pair[k]+1, site.pair[k], year.pair[k]] 
    } # k
    for (t in 1:(nyr-1)){ # this index for dens dep ipm bc mn.prod[t-1] would cause a zero index
      for (s in 1:nsite){
        for (tr in 1:2){
          log(mn.prod[tr, s, t]) <- lmu.prod[s] + 
            gamma*c(0,1)[tr] + 
            eta[9, s, t]
        } # tr
        num.treat[t, s] ~ dbinom( ptreat[t, s], lamAD[t, s] )
        logit(ptreat[t, s]) <-  lmu.treat[s] + eta[10,s,t]
      }} # s t
    
    # Count likelihoods, state-space model, and observation process    
    for (s in 1:nsite){    
      lNAD[s] ~ dunif(-10,7)
      lN[s] ~ dunif(-10,7)
    }
    
    for (t in 1:nyr){
      for (s in 1:nsite){
        countsAdults[t, s] ~ dpois(lamAD[t, s]) # adult females
        log(lamAD[t, s]) <- lNAD[s] + deltas[1]*effort2[t, s] + deltas[2]*effort2[t, s]^2  
        log(lamFY[t, s]) <- lN[s] + deltas[3]*effort2[t, s] + deltas[4]*effort2[t, s]^2
      }# s
      # First-years at different sites have different distributions
      # for better model fit
      countsFY[t, 1] ~ dpois(lamFY[t, 1]) # doesn't have any zeroes so poisson fits
      countsFY[t, 2] ~ dnegbin(pp[t], r) # first year females, includes translocated/hacked
      pp[t] <- r/(r+(lamFY[t, 2] ))
    } # t
    r ~ dexp(0.05)

  } )

#**********************
#* Function to run model in NIMBLE
#**********************
run_ipm <- function(info, datl, constl, code, outp){
  library('nimble')
  library('coda')
  #source("/bsuscratch/brianrolek/riha_ipm/MCMCvis.R")
  
  params <- c(
    # pop growth 
    #"lambda",
    # fecundity
    "lmu.prod", "gamma", "rr", "mn.prod", 
    "mu.treat", "ptreat",
    # survival 
    #"mus", "lmus", "betas", 
    "deltas", #"alphas",
    # abundance
    "lNAD", "lN",
    #"NB", "NF", "NFY", "N", "NAD",
    "r",
    #"N", "Ntot",
    # immigration 
    #"omega",
    # error terms
    "eta", "sds", "Ustar", "R", "z.score"
    # yearly summaries
  )
  
  # initial values  
  inits.func <- function (){
    # sample from working inits
    ind <- sample(1:ncol(outp$gamma), size=1, replace=TRUE)
    
    U <- outp$Ustar[,,ind]
    U <- rbind(U, rep(0, constl$p-1) )
    U <- cbind(U, c(runif(constl$p-1, -0.1, 0.1), runif(1, 0, 0.3) ))
    list(
      # fecundity
      lmu.prod = outp$lmu.prod[,ind],
      gamma = outp$gamma[,ind],
      rr = outp$rr[,ind],
      mu.treat = runif(2),
      #sd.treat = runif(1),
      # survival
      z = info$inits$z,
      mus = outp$mus[,,ind],
      betas = outp$betas[,ind],
      deltas = c(0,0, outp$deltas[3:8,ind]),
      sds = c(outp$sds[,ind], runif(1)),
      Ustar = U,
      # counts
      r = outp$r[,ind],
      #lNAD= runif(constl$nyr*constl$nsite, 4, 6.5), dim=c(constl$nyr,constl$nsite)),
      #lN= array(runif(constl$nyr*constl$nsite, 3, 6), dim=c(constl$nyr,constl$nsite)),
      lNAD= c(log(80), log(43)),
      lN= runif(2, 3, 6),
      #N = N.ar, # sample from inits of chains that worked
      z.score = array(runif(constl$p*constl$nsite*constl$nyr, -0.1, 0.1), 
                      dim=c(constl$p, constl$nsite, constl$nyr))#,

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
  #nc=1; nt=100; ni=100000; nb=50000
  #nc=1; nt=100; ni2=250000; nb=150000
  #nc=1; nt=1; ni=3; nb=1
  nc=1; nt=1; ni=50000; nb=5000
  
  post3 <- runMCMC(chmc,
                   niter = ni,
                   nburnin = nb,
                   nchains = nc,
                   thin = nt,
                   samplesAsCodaMCMC = T,
                   setSeed = info$seed,
                   inits = inits)
  
  return(post3)
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

flnm2 <- paste0("/bsuscratch/brianrolek/riha_ipm/outputs/ipm-simp-",
                substr(Sys.time(), 1, 10), ".rdata")
save(post, file = flnm2)
# 
# save(post, mycode,
#      file="/bsuscratch/brianrolek/riha_ipm/outputs/ipm_2025-11-04-11i.rdata")