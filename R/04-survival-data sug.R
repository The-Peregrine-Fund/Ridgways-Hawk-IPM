library ('nimble')
library ('MCMCvis')
library ('coda')
library('parallel')

setwd("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\GitHub\\Ridgways-Hawk-IPM\\")
load("data\\data.rdata")
set.seed(5757575)

constl$p <- 8
par_info[[1]]$inits$sds <- par_info[[1]]$inits$sds[1:constl$p]
par_info[[2]]$inits$sds <- par_info[[2]]$inits$sds[1:constl$p]
par_info[[3]]$inits$sds <- par_info[[3]]$inits$sds[1:constl$p]
par_info[[4]]$inits$sds <- par_info[[4]]$inits$sds[1:constl$p]

par_info[[1]]$inits$Ustar <- par_info[[1]]$inits$Ustar[1:constl$p, 1:constl$p]
par_info[[2]]$inits$Ustar <- par_info[[2]]$inits$Ustar[1:constl$p, 1:constl$p]
par_info[[3]]$inits$Ustar <- par_info[[3]]$inits$Ustar[1:constl$p, 1:constl$p]
par_info[[4]]$inits$Ustar <- par_info[[4]]$inits$Ustar[1:constl$p, 1:constl$p]

runif(constl$p) ->
par_info[[1]]$inits$mus -> par_info[[2]]$inits$mus -> par_info[[3]]$inits$mus ->par_info[[4]]$inits$mus

datl$mu.zeroes <- rep(0, constl$p)

da.matrix <- array(4, dim=c(100, constl$nyr))
datl$yaug <- rbind(datl$y, da.matrix)
par_info[[1]]$inits$z <- rbind(par_info[[1]]$inits$z, da.matrix)
par_info[[2]]$inits$z <- rbind(par_info[[2]]$inits$z, da.matrix)
par_info[[3]]$inits$z <- rbind(par_info[[3]]$inits$z, da.matrix)
par_info[[4]]$inits$z <- rbind(par_info[[4]]$inits$z, da.matrix)

constl$M <- nrow(datl$yaug)
constl$first <- c(constl$first, rep(1, 100))
# To Do 
# x(1) Get basic model working and converging
# x(2) Add sites as index
# (3) Add inter-site correlations as an another level of correlations
#       in addition to intra-site demographic correlations. 
#       See Kery and Schaub IPM book ch16

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
      betas[j] ~ dnorm(0, sd=2)
      lmus[j] <- logit(mus[j])
      mus[j] ~ dbeta(1,1)
    } # p
      # for (m in 1:4){ # population
      #   lmus[j,m] <- logit(mus[j,m])
      #   mus[j,m] ~ dbeta(1,1)
      # }  }   # m population #s sex #h hacked 
    
    # estimated using the multivariate normal distribution
    R[1:p,1:p] <- t(Ustar[1:p,1:p]) %*% Ustar[1:p,1:p] # calculate rhos, correlation coefficients
    Ustar[1:p,1:p] ~ dlkj_corr_cholesky(eta=1.0, p=p) # Ustar is the Cholesky decomposition of the correlation matrix
    U[1:p,1:p] <- uppertri_mult_diag(Ustar[1:p, 1:p], sds[1:p])
    # multivariate normal for temporal variance
    for (t in 1:(nyr-1)){
      eps[1:p,t] ~ dmnorm(mu.zeroes[1:p],
                          cholesky = U[1:p, 1:p], prec_param = 0)
    }
    # Calculate counts from HMM 
    for (i in 1:M){
    for (t in first[i]:nyr){
      zFY[i,t] <- equals(z[i,t], 1)
      zF[i,t] <- equals(z[i,t], 2)
      zB[i,t] <- equals(z[i,t], 3)
    } }
    for (t in 1:nyr){
      NFY[t] <- sum(zFY[1:M,t])
      NF[t] <- sum(zF[1:M,t])
      NB[t] <- sum(zB[1:M,t])
      } 

    # Calculate an average for each year 
    # add site here later
    for (t in 1:(nyr-1)){
      eta.phiFY[t] <- mean(phiFY[1:nind,t])
      eta.phiA[t] <- mean(phiA[1:nind,t])
      eta.phiB[t] <- mean(phiB[1:nind,t])
      eta.psiFYB[t] <- mean(psiFYB[1:nind,t])
      eta.psiAB[t] <- mean(psiAB[1:nind,t])
      eta.psiBA[t] <- mean(psiBA[1:nind,t])
      eta.pA[t] <- mean(pA[1:nind,t])
      eta.pB[t] <- mean(pB[1:nind,t])
    } # t
    
    for (i in 1:M){
      for (t in 1:(nyr-1)){
        #Survival
        logit(phiFY[i,t]) <- eps[1,t] + lmus[1] #+ betas[1]*hacked[i]  # first year
        logit(phiA[i,t]) <- eps[2,t] + lmus[2] #+  betas[2]*hacked[i] # nonbreeder
        logit(phiB[i,t]) <- eps[3,t]  + lmus[3] #+ betas[3]*hacked[i] # breeder
        #Recruitment
        logit(psiFYB[i,t]) <- eps[4,t] + lmus[4] #+ betas[4]*hacked[i] # first year to breeder
        logit(psiAB[i,t]) <- eps[5,t] + lmus[5] #+ betas[5]*hacked[i] # nonbreeder to breeder
        logit(psiBA[i,t]) <- eps[6,t] + lmus[6] #+ betas[6]*hacked[i] # breeder to nonbreeder
        #Re-encounter
        logit(pA[i,t]) <- eps[7,t] + lmus[7] #+ betas[7]*hacked[i] # resight of nonbreeders
        logit(pB[i,t]) <- eps[8,t] + lmus[8] #+ betas[8]*hacked[i] # resight of breeders
      }#t
    }#i
    
    # Define state-transition and observation matrices
    for (i in 1:M){  
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
    for (i in 1:M){
      # Define latent state at first capture
      z[i,first[i]] <- yaug[i,first[i]]
      for (t in (first[i]+1):nyr){
        # State process: draw S(t) given S(t-1)
        z[i,t] ~ dcat(ps[z[i,t-1], i, t-1, 1:4])
        # Observation process: draw O(t) given S(t)
        yaug[i,t] ~ dcat(po[z[i,t], i, t-1, 1:4])
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
  
  
  params <- c( "mus", "lmus", "eps", "betas", "sds",
               "Ustar", "U", "R",              
               'eta.phiFY', 'eta.phiA', 'eta.phiB',
               'eta.psiFYB', 'eta.psiAB', 'eta.psiBA',
               'eta.pA', 'eta.pB',
               "NFY", "NF", "NB", "z")
  
  nc=1; nt=1; ni=20000; nb=10000
  #nc=1; nt=50; ni=50000; nb=25000
  
  # build derivatives and enable WAIC breaks this
  mod <- nimbleModel(code, calculate=T, 
                     constants = constl, 
                     data = datl, 
                     inits = info$inits)
  mod$simulate()      
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

# constl <- constl[-5] # remove duplicate nyr

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
