## ---- ipm1 --------
#################################
# The model
################################
library('nimble')
library('parallel')
library ('coda')
load("/bsuscratch/brianrolek/riha_ipm/data.rdata")
#load("data/data.rdata")

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
    ###################################################
    # Priors and constraints
    ###################################################
    # survival, recruitment, and detection can be correlated
    for (j in 1:8){ 
        betas[j] ~ dunif(-10, 10)  # prior for coefficients
        lmus[j] <- logit(mus[j])
        mus[j] ~ dbeta(1,1) # prior for means
      }     # m population #s sex #h hacked 
    for (j in 1:p){ sds[j] ~ dexp(1) }# prior for temporal variation
    
    # estimated using the multivariate normal distribution
    R[1:p,1:p] <- t(Ustar[1:p,1:p]) %*% Ustar[1:p,1:p] # calculate rhos, correlation coefficients
    Ustar[1:p,1:p] ~ dlkj_corr_cholesky(eta=1.3, p=p) # Ustar is the Cholesky decomposition of the correlation matrix
    U[1:p,1:p] <- uppertri_mult_diag(Ustar[1:p, 1:p], sds[1:p])
    # multivariate normal for temporal variance
    for (t in 1:nyr){ # survival params only have nyr-1, no problem to simulate from however
      eps[1:p,t] ~ dmnorm(mu.zeroes[1:p],
                          cholesky = U[1:p, 1:p], prec_param = 0)
    }
    #######################
    # Derived params
    #######################
    for (t in 1:(nyr-1)){
      lambda[t] <-  Ntot[t+1]/(Ntot[t])
      loglambda[t] <- log(lambda[t])
    } #t
    
    ###############################
    # Likelihood for fecundity
    ###############################
    # Priors for number of fledglings
    lmu.brood ~ dnorm(0, sd=5)
    delta ~ dunif(-10, 10)
    sig.brood ~ dexp(1)
    
    # priors for nest success
    lmu.nest <- logit(mu.nest)
    mu.nest ~ dbeta(1, 1) 
    gamma ~ dunif(-10, 10)
    
    # Two models for fecundity, (1) brood size and (2) nest success  
    # Brood size
    for (k in 1:nbrood){
      brood[k] ~ dlnorm(lam[k], sdlog=sig.brood) # truncating seems to break something
      lam[k] <- lmu.brood +
                delta*treat.brood[k] +
                eps[9,year.brood[k] ]
    } # k
    # Nest success       
    for (n in 1:nnest){
      nest.success[n] ~ dbern( nu[n] )
      logit(nu[n]) <- lmu.nest + 
                      gamma*treat.nest[n] + 
                      eps[10,year.nest[n] ]
    } # n 
    # derive yearly brood size for population model
    for (t in 1:nyr){
      for (xx in 1:brood.end[t]){
        broodmat[t,xx] <- lam[yrind.brood[xx,t]]
      } # xx
      for (xxx in 1:nest.end[t]){
        nestmat[t,xxx] <- nu[yrind.nest[xxx,t]]
      } # xxx
      mn.brood[t] <- exp( mean(broodmat[t,1:brood.end[t]]) )
      mn.nest[t] <- mean( nestmat[t,1:nest.end[t]] )
      mn.f[t] <- mn.brood[t]*mn.nest[t] # calc fecundity
    }
    
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
    for (v in 1:7){ N[v,1] ~ dcat(pPrior[1:102]) }
    # Abundance for years > 1
    for (t in 1:(nyr-1)){
      # Number of wild born juvs
      N[1, t+1] ~ dpois( (NFY[t]*eta.phiFY[t]*eta.psiFYB[t] + # first year breeders
                          NF[t]*eta.phiA[t]*eta.psiAB[t] + # nonbreeders to breeders
                          NB[t]*eta.phiB[t]*(1-eta.psiBA[t])) # breeders remaining
                           *mn.f[t+1] ) # end Poisson
      # Abundance of nonbreeders
      ## Second year nonbreeders
      N[2, t+1] ~ dbin(eta.phiFY[t]*(1-eta.psiFYB[t]), NFY[t]) # Nestlings to second year nonbreeders
      ## Adult nonbreeders
      N[3, t+1] ~ dbin(eta.phiA[t]*(1-eta.psiAB[t]), NF[t]) # Nonbreeders to nonbreeders
      N[4, t+1] ~ dbin(eta.phiB[t]*eta.psiBA[t], NB[t]) # Breeders to nonbreeders
      # Abundance of breeders
      ## Second year breeders
      N[5, t+1] ~ dbin(eta.phiFY[t]*eta.psiFYB[t], NFY[t]) # Nestlings to second year breeders
      ## Adult breeders
      N[6, t+1] ~ dbin(eta.phiA[t]*eta.psiAB[t], NF[t]) # Nonbreeder to breeder
      N[7, t+1] ~ dbin(eta.phiB[t]*(1-eta.psiBA[t]), NB[t]) # Breeder to breeder
    } # t

    for (t in 1:nyr){
      NFY[t] <- N[1, t] #+ aug[t] # Includes translocated first years
      NF[t] <- sum(N[2:4, t])  # number of adult nonbreeders
      NB[t] <- sum(N[5:7, t]) # number of adult breeders
      Ntot[t] <- sum(N[1:7, t]) # total number
    } # t
    
    # Observation process    
    for (t in 1:nyr){
      # counts.marked[1, t] ~ dpois(NFY[t]) # first year males, includes translocated/hacked
      # counts.marked[2, t] ~ dpois(NF[t]) # nonbreeding adult males 
      # counts.marked[3, t] ~ dpois(NB[t]) # breeding males 
      # log(NFY.mu[t]) ~ dnorm(0, sd = 5)
      # log(NF.mu[t]) ~ dnorm(0, sd = 5)
      # log(NB.mu[t]) ~ dnorm(0, sd = 5)
      # NF[t] ~ dpois(NF.mu[t])
      # NB[t] ~ dpois(NB.mu[t])
      counts[1, t] ~ dpois(NFY[t]) # first year males, includes translocated/hacked
      counts[2, t] ~ dbin(eta.pA[t], NF[t]) # nonbreeding adult males 
      #NF.mu[t] ~ dpois(NF[t])
      counts[3, t] ~ dbin(eta.pB[t], NB[t]) # breeding males 
      #NB.mu[t] ~ dpois(NB[t])
    } # t
    
    ###################
    # Assess GOF of the state-space models for counts
    # Step 1: Compute statistic for observed data
    # Step 2: Use discrepancy measure: mean absolute error
    # Step 3: Use test statistic: number of turns
    ###################
    for (t in 1:nyr){
      c.expB[t] <- NB[t]  # expected counts adult breeder
      c.expA[t] <- NF[t] # nonbreeder
      c.expFY[t] <- NFY[t]
      c.obsB[t] <- counts[3, t]
      c.obsA[t] <- counts[2, t]
      c.obsFY[t] <- counts[1, t]  # first year
      c.repB[t] ~ dbin( eta.pB[t], NB[t] ) # simulated counts
      c.repA[t] ~ dbin( eta.pA[t], NF[t] )
      c.repFY[t] ~ dpois( NFY[t] )
      dssm.obsB[t] <- abs( ( (c.obsB[t]) - (c.expB[t]) ) / (c.obsB[t]+0.001)  )
      dssm.obsA[t] <- abs( ( (c.obsA[t]) - (c.expA[t]) ) / (c.obsA[t]+0.001)  )
      dssm.obsFY[t] <- abs( ( (c.obsFY[t]) - (c.expFY[t]) ) / (c.obsFY[t]+0.001)  )
      dssm.repB[t] <- abs( ( (c.repB[t]) - (c.expB[t]) ) / (c.repB[t]+0.001) )
      dssm.repA[t] <- abs( ( (c.repA[t]) - (c.expA[t]) ) / (c.repA[t]+0.001) )
      dssm.repFY[t] <- abs( ( (c.repFY[t]) - (c.expFY[t]) ) / (c.repFY[t]+0.001) )
    } # t
    dmape.obs[1] <- sum(dssm.obsB[1:nyr])
    dmape.obs[2] <- sum(dssm.obsA[1:nyr])
    dmape.obs[3] <- sum(dssm.obsFY[1:nyr])
    dmape.rep[1] <- sum(dssm.repB[1:nyr])
    dmape.rep[2] <- sum(dssm.repA[1:nyr])
    dmape.rep[3] <- sum(dssm.repFY[1:nyr])
    tvm.obs[1] <- sd(c.obsB[1:nyr])^2/mean(c.obsB[1:nyr])
    tvm.obs[2] <- sd(c.obsA[1:nyr])^2/mean(c.obsA[1:nyr])
    tvm.obs[3] <- sd(c.obsFY[1:nyr])^2/mean(c.obsFY[1:nyr])
    tvm.rep[1] <- sd(c.repB[1:nyr])^2/mean(c.repB[1:nyr])
    tvm.rep[2] <- sd(c.repA[1:nyr])^2/mean(c.repA[1:nyr])
    tvm.rep[3] <- sd(c.repFY[1:nyr])^2/mean(c.repFY[1:nyr])
    # # Test statistic for number of turns
    for (t in 1:(nyr-2)){
      tt1.obsB[t] <- step(c.obsB[t+2] - c.obsB[t+1])
      tt2.obsB[t] <- step(c.obsB[t+1] - c.obsB[t])
      tt3.obsB[t] <- equals(tt1.obsB[t] + tt2.obsB[t], 1)
      tt1.obsA[t] <- step(c.obsA[t+2] - c.obsA[t+1])
      tt2.obsA[t] <- step(c.obsA[t+1] - c.obsA[t])
      tt3.obsA[t] <- equals(tt1.obsA[t] + tt2.obsA[t], 1)
      tt1.obsFY[t] <- step(c.obsFY[t+2] - c.obsFY[t+1])
      tt2.obsFY[t] <- step(c.obsFY[t+1] - c.obsFY[t])
      tt3.obsFY[t] <- equals(tt1.obsFY[t] + tt2.obsFY[t], 1)
    } # t
    tturn.obs[1] <- sum(tt3.obsB[1:(nyr-2)])
    tturn.obs[2] <- sum(tt3.obsA[1:(nyr-2)])
    tturn.obs[3] <- sum(tt3.obsFY[1:(nyr-2)])

    for (t in 1:(nyr-2)){
      tt1.repB[t] <- step(c.repB[t+2] - c.repB[t+1])
      tt2.repB[t] <- step(c.repB[t+1] - c.repB[t])
      tt3.repB[t] <- equals(tt1.repB[t] + tt2.repB[t], 1)
      tt1.repA[t] <- step(c.repA[t+2] - c.repA[t+1])
      tt2.repA[t] <- step(c.repA[t+1] - c.repA[t])
      tt3.repA[t] <- equals(tt1.repA[t] + tt2.repA[t], 1)
      tt1.repFY[t] <- step(c.repFY[t+2] - c.repFY[t+1])
      tt2.repFY[t] <- step(c.repFY[t+1] - c.repFY[t])
      tt3.repFY[t] <- equals(tt1.repFY[t] + tt2.repFY[t], 1)
    } # t
    tturn.rep[1] <- sum(tt3.repB[1:(nyr-2)])
    tturn.rep[2] <- sum(tt3.repA[1:(nyr-2)])
    tturn.rep[3] <- sum(tt3.repFY[1:(nyr-2)])
    
    ################################
    # Likelihood for survival
    ################################ 
    # Calculate an average for each year 
    # add site here later
    for (t in 1:nyr){
      eta.phiFY[t] <- mean(phiFY[1:nind,t])
      eta.phiA[t] <- mean(phiA[1:nind,t])
      eta.phiB[t] <- mean(phiB[1:nind,t])
      eta.psiFYB[t] <- mean(psiFYB[1:nind,t])
      eta.psiAB[t] <- mean(psiAB[1:nind,t])
      eta.psiBA[t] <- mean(psiBA[1:nind,t])
      eta.pA[t] <- mean(pA[1:nind,t])
      eta.pB[t] <- mean(pB[1:nind,t])
    } # t
    
    for (i in 1:nind){
      for (t in 1:nyr){
        #Survival
        logit(phiFY[i,t]) <- eps[1,t] + lmus[1] + betas[1]*hacked[i]  # first year
        logit(phiA[i,t]) <- eps[2,t] + lmus[2] +  betas[2]*hacked[i] # nonbreeder
        logit(phiB[i,t]) <- eps[3,t] + lmus[3] + betas[3]*hacked[i] # breeder
        #Recruitment
        logit(psiFYB[i,t]) <- eps[4,t] + lmus[4] + betas[4]*hacked[i] # first year to breeder
        logit(psiAB[i,t]) <- eps[5,t] + lmus[5] + betas[5]*hacked[i] # nonbreeder to breeder
        logit(psiBA[i,t]) <- eps[6,t] + lmus[6] + betas[6]*hacked[i] # breeder to nonbreeder
        #Re-encounter
        logit(pA[i,t]) <- eps[7,t] + lmus[7] + betas[7]*hacked[i] # resight of nonbreeders
        logit(pB[i,t]) <- eps[8,t] + lmus[8] + betas[8]*hacked[i] # resight of breeders
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
            "lmu.brood", "delta", "sig.brood", #"sig.brood.t",
            "lmu.nest", "mu.nest", "gamma", #"sig.nest",
            "mn.brood", "mn.nest", "mn.f", 
             # survival 
             "mus", "lmus", "betas",
             # abundance
             "NB", "NF", "NFY", 
             #"NB.mu", "NF.mu", "NFY.mu",
             "N", "Ntot",
             # error terms
             "eps", "sds", "Ustar", "U", "R",
             # yearly summaries
             'eta.phiFY', 'eta.phiA', 'eta.phiB',
             'eta.psiFYB', 'eta.psiAB', 'eta.psiBA',
             'eta.pA', 'eta.pB',
             # goodness of fit
             "f.dmape.obs", "f.dmape.rep",
             "f.tvm.obs", "f.tvm.rep",
             "dmape.obs", "dmape.rep",
             "tvm.obs", "tvm.rep",
             "tturn.obs", "tturn.rep"
)

 #n.chains=1; n.thin=200; n.iter=500000; n.burnin=300000
#n.chains=1; n.thin=50; n.iter=100000; n.burnin=50000
n.chains=1; n.thin=10; n.iter=10000; n.burnin=0
# n.chains=1; n.thin=1; n.iter=200; n.burnin=100

mod <- nimbleModel(code, 
                   constants = constl, 
                   data = datl, 
                   inits = info$inits, 
                   buildDerivs = FALSE,
                   calculate=T 
                   # nodes=c('eta.phiFY', 'eta.phiA', 'eta.phiB',
                   #         'eta.psiFYB', 'eta.psiAB', 'eta.psiBA',
                   #         'eta.pA', 'eta.pB',"NB", 
                   #         "NF", "NFY", "Ntot",
                   #         "lmus", "lambda", 
                   #         "mn.brood", "mn.nest", "mn.f",
                   #         "f.dmape.obs", "f.dmape.rep",
                   #         "f.tvm.obs", "f.tvm.rep",
                   #         "dmape.obs", "dmape.rep",
                   #         "tvm.obs", "tvm.rep",
                   #         "tturn.obs", "tturn.rep"
                   #         )
                   ) # doesn't work when TRUE, no hope for HMC


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
     file="/bsuscratch/brianrolek/riha_ipm/outputs/ipm_noburnin.rdata")

#save(post,  
#     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm.Rdata")