## ---- ipm1 --------
#################################
# The model
################################
library('nimble')
library('nimbleHMC')
library('MCMCvis')
library('parallel')
library ('coda')
load("/bsuscratch/brianrolek/riha_ipm/data.RData")
load("data/data.RData")

# from IPMbook package
dUnif <- function (lower, upper) 
{
  A <- round(lower)
  B <- round(upper)
  nrow <- length(lower)
  out <- matrix(0, nrow = nrow, ncol = max(B))
  n <- B - A + 1
  for (i in 1:nrow) out[i, A[i]:B[i]] <- rep(1/n[i], n[i])
  return(drop(out))
}

run_ipm <- function(seed, datl, constl){
  library('nimble')
  library('coda')
  library ('nimbleHMC')
  
  code <- nimbleCode(
      {
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
  ###################################################
  # Priors and constraints
  ###################################################
  # survival, recruitment, and detection can be correlated
  for (j in 1:p){ # coefficient
    sds[j] ~ dexp(1) # prior for temporal variation
    betas[j] ~ dnorm(0, sd=2) # prior for coefficients
    for (m in 1:4){ # population
      lmus[j,m] <- logit(mus[j,m])
      mus[j,m] ~ dbeta(1,1) # prior for means
    }  }   # m population #s sex #h hacked 
  
  # estimated using the multivariate normal distribution
  R[1:p,1:p] <- t(Ustar[1:p,1:p]) %*% Ustar[1:p,1:p] # calculate rhos, correlation coefficients
  Ustar[1:p,1:p] ~ dlkj_corr_cholesky(eta=1.0, p=p) # Ustar is the Cholesky decomposition of the correlation matrix
  U[1:p,1:p] <- uppertri_mult_diag(Ustar[1:p, 1:p], sds[1:p])
  # multivariate normal for temporal variance
  for (t in 1:(nyr-1)){
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
  lmu.brood ~ T(dnorm(0, sd=10), log(mintrunc), log(maxtrunc))
  delta ~ dunif(-5, 5)
  gamma ~ dunif(-5, 5)
  # priors for nest success
  lmu.nest <- logit(mu.nest)
  mu.nest ~ dbeta(1,1) 
  sig.brood ~ dexp(1)
  sig.nest ~ dexp(1)
  sig ~ dexp(1)
  
  # Two models for fecundity, (1) brood size and (2) nest success  
  # Brood size
  for (k in 1:nbrood){
    brood[k] ~ dlnorm(lam[k], sdlog=sig) # truncating seems to break something
    lam[k] <- lmu.brood +
              delta*treat.brood[k] +
              eps[9, year.brood[k] ]
  } # k
  # Nest success       
  for (n in 1:nnest){
    nest.success[n] ~ dbern( nu[n] )
    logit(nu[n]) <- lmu.nest + 
                    gamma*treat.nest[n] + 
                    eps[10, year.nest[n] ]
  } # n 
  
  # summarize yearly brood size for population model
  # accounts for nest treatments
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
  
  # # # GOF for number of fledglings
  for (k in 1:nbrood){
    f.obs[k] <- brood[k] # observed counts
    f.exp[k] <- lam[k] # expected counts adult breeder
    f.rep[k] ~ dlnorm(lam[k], sdlog=sig) # expected counts
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
  N[v,1] ~ dcat(pInit) 
  }
  # Abundance for years > 1
  for (t in 2:(nyr-1)){
  # Number of wild born juvs
  N[1,t+1] ~ dpois( 
                    (NFY[t]*eta.phiFY[t]*eta.psiFYB[t])*mn.f[t]/2 + # first year breeders
                    (NF[t]*eta.phiA[t]*eta.psiAB[t] + # nonbreeders to breeders
                      NB[t]*eta.phiB[t]*(1-eta.psiBA[t]) # breeders remaining
                      ) *mn.f[t]/2 # Adult breeders
                    ) # end Poisson
  # Abundance of nonbreeders
  ## Second year nonbreeders
  N[2,t+1] ~ dbin(eta.phiFY[t]*(1-eta.psiNB[t]), NFY[t]) # Nestlings to second year nonbreeders
  ## Adult nonbreeders
  N[3,t+1] ~ dbin(eta.phiA[t]*(1-eta.psiAB[t]), NF[t]) # Nonbreeders to nonbreeders
  N[4,t+1] ~ dbin(eta.phiB[t]*eta.psiBA[t], NB[t]) # Breeders to nonbreeders
  # Abundance of breeders
  ## Second year breeders
  N[5,t+1] ~ dbin(eta.phiFY[t]*eta.psiNB[t], NFY[t]) # Nestlings to second year breeders
  ## Adult breeders
  N[6,t+1] ~ dbin(eta.phiA[t]*eta.psiAB[t], NF[t]) # Nonbreeder to breeder
  N[7,t+1] ~ dbin(eta.phiB[t]*(1-eta.psiBA[t]), NB[t]) # Breeder to breeder
  } # t
  
  for (t in 1:nyr){
  NFY[t] <- N[1,t] + aug[t] # Includes translocated first years
  NF[t] <- sum(N[c(2,3,4),t])  # number of adult nonbreeders
  NB[t] <- sum(N[c(5,6,7),t]) # number of adult breeders
  Ntot[t] <- sum(N[c(1,2,3,4,5,6,7),t]) # total number
  } # t
  
  # Observation process    
  for (t in 2:nyr){
  counts.marked[1, t] ~ dpois(NFY[t]) # first year males, includes translocated/hacked
  counts.marked[2, t] ~ dpois(NF[t]) # nonbreeding adult males    
  counts.marked[3, t] ~ dpois(NB[t]) # breeding males 
  } # t
  
  ###################
  # Assess GOF of the state-space models for counts
  # Step 1: Compute statistic for observed data
  # Step 2: Use discrepancy measure: mean absolute error
  # Step 3: Use test statistic: number of turns
  ###################
  # for (t in 2:nyr){ 
  # c.expB[t] <- NB[t] + 0.001 # expected counts adult breeder
  # c.expA[t] <- NF[t] + 0.001 # nonbreeder
  # c.expFYW[t] <- countFYW[t] + 0.001 # first year
  # c.obsB[t] <- countB[t] + 0.001
  # c.obsA[t] <- countA[t] + 0.001
  # c.obsFYW[t] <- countNFYW[t] + 0.001
  # dssm.obsB[t] <- abs( ( (c.obsB[t]) - (c.expB[t]) ) / (c.obsB[t]+0.001)  )
  # dssm.obsA[t] <- abs( ( (c.obsA[t]) - (c.expA[t]) ) / (c.obsA[t]+0.001)  )
  # dssm.obsFYW[t] <- abs( ( (c.obsFYW[t]) - (c.expFYW[t]) ) / (c.obsFYW[t]+0.001)  )
  # } # t
  # dmape.obs[1] <- sum(dssm.obsB[2:n.yr])
  # dmape.obs[2] <- sum(dssm.obsA[2:n.yr])
  # dmape.obs[3] <- sum(dssm.obsFYW[2:n.yr])
  # # Compute fit statistic for replicate data
  # # Mean absolute error
  # for (t in 2:nyr){ 
  # c.repB[t] ~ dpois( NB[t] ) # expected counts
  # c.repA[t] ~ dpois( NF[t] ) 
  # c.repFYW[t] ~ dpois( NFYW[1,t] )
  # } # t
  # for (t in 2:nyr){ 
  # dssm.repB[t] <- abs( ( (c.repB[t]) - (c.expB[t]) ) / (c.repB[t]+0.001) )
  # dssm.repA[t] <- abs( ( (c.repA[t]) - (c.expA[t]) ) / (c.repA[t]+0.001) )
  # dssm.repFYW[t] <- abs( ( (c.repFYW[t]) - (c.expFYW[t]) ) / (c.repFYW[t]+0.001) )
  # } # t
  # dmape.rep[1] <- sum(dssm.repB[2:n.yr])
  # dmape.rep[2] <- sum(dssm.repA[2:n.yr])
  # dmape.rep[3] <- sum(dssm.repFYW[2:n.yr])
  # 
  # # Test statistic for number of turns
  # for (t in 2:(nyr-2)){
  # tt1.obsB[t] <- step(countB[t+2] - countB[t+1])
  # tt2.obsB[t] <- step(countB[t+1] - countB[t])
  # tt3.obsB[t] <- equals(tt1.obsB[t] + tt2.obsB[t], 1)
  # tt1.obsA[t] <- step(countA[t+2] - countA[t+1])
  # tt2.obsA[t] <- step(countA[t+1] - countA[t])
  # tt3.obsA[t] <- equals(tt1.obsA[t] + tt2.obsA[t], 1)
  # tt1.obsFYW[t] <- step(countFYW[t+2] - countFYW[t+1])
  # tt2.obsFYW[t] <- step(countFYW[t+1] - countFYW[t])
  # tt3.obsFYW[t] <- equals(tt1.obsFYW[t] + tt2.obsFYW[t], 1)
  # } # t
  # tturn.obs[1] <- sum(tt3.obsB[2:(n.yr-2)])
  # tturn.obs[2] <- sum(tt3.obsA[2:(n.yr-2)])
  # tturn.obs[3] <- sum(tt3.obsFYW[2:(n.yr-2)])
  # 
  # for (t in 2:(nyr-2)){
  # tt1.repB[t] <- step(c.repB[t+2] - c.repB[t+1])
  # tt2.repB[t] <- step(c.repB[t+1] - c.repB[t])
  # tt3.repB[t] <- equals(tt1.repB[t] + tt2.repB[t], 1)
  # tt1.repA[t] <- step(c.repA[t+2] - c.repA[t+1])
  # tt2.repA[t] <- step(c.repA[t+1] - c.repA[t])
  # tt3.repA[t] <- equals(tt1.repA[t] + tt2.repA[t], 1)
  # tt1.repFYW[t] <- step(c.repFYW[t+2] - c.repFYW[t+1])
  # tt2.repFYW[t] <- step(c.repFYW[t+1] - c.repFYW[t])
  # tt3.repFYW[t] <- equals(tt1.repFYW[t] + tt2.repFYW[t], 1)
  # } # t
  # tturn.rep[1] <- sum(tt3.repB[2:(n.yr-2)])
  # tturn.rep[2] <- sum(tt3.repA[2:(n.yr-2)])
  # tturn.rep[3] <- sum(tt3.repFYW[2:(n.yr-2)])
  
  
  ################################
  # Likelihood for survival
  ################################ 
  # Calculate an average for each year 
  # add site here later
  for (t in 1:(nyr-1)){
    eta.phiFY[t] <- mean(phiFY[1:nind,t])
    eta.phiA[t] <- mean(phiFY[1:nind,t])
    eta.phiB[t] <- mean(phiFY[1:nind,t])
    eta.psiNB[t] <- mean(psiFY[1:nind,t])
    eta.psiAB[t] <- mean(psiAB[1:nind,t])
    eta.psiBA[t] <- mean(psiBA[1:nind,t])
    eta.pA[t] <- mean(pA[1:nind,t])
    eta.pB[t] <- mean(pB[1:nind,t])
  } # t
  
  for (i in 1:nind){
    for (t in 1:(nyr-1)){
      #Survival
      logit(phiFY[i,t]) <- eps[1,t] + lmus[1, site[i,t]] + betas[1]*hacked[i]  # first year
      logit(phiA[i,t]) <- eps[2,t] + lmus[2, site[i,t]] +  betas[2]*hacked[i] # nonbreeder
      logit(phiB[i,t]) <- eps[3,t]  + lmus[3, site[i,t]] + betas[3]*hacked[i] # breeder
      #Recruitment
      logit(psiFYB[i,t]) <- eps[4,t] + lmus[4, site[i,t]] + betas[4]*hacked[i] # first year to breeder
      logit(psiAB[i,t]) <- eps[5,t] + lmus[5, site[i,t]] + betas[5]*hacked[i] # nonbreeder to breeder
      logit(psiBA[i,t]) <- eps[6,t] + lmus[6, site[i,t]] + betas[6]*hacked[i] # breeder to nonbreeder
      #Re-encounter
      logit(pA[i,t]) <- eps[7,t] + lmus[7, site[i,t]] + betas[7]*hacked[i] # resight of nonbreeders
      logit(pB[i,t]) <- eps[8,t] + lmus[8, site[i,t]] + betas[8]*hacked[i] # resight of breeders
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

get.first <- function(x) min(which(x!=5))
f <- apply(datl$y, 1, get.first)
get.last <- function(x) max(which(x!=5))
l <- apply(datl$y, 1, get.last)
TFmat <- is.na(z.inits) & is.na(datl$z)
for (i in 1:dim(TFmat)[1]){  TFmat[i,1:f[i]] <- FALSE }
z.inits[TFmat] <- sample(size=445, c(2,3), replace=T, prob=c(0.5, 0.5) ) 

# create inits for Rhos
p <- 11
Ustar <- array(runif(p*p, -0.5, 0.5), dim=c(p,p))
diag(Ustar) <- 1 # set diagonal to 1
Ustar[lower.tri(Ustar)] <- 0 # set lower diag to zero
t(Ustar)%*%Ustar

# Set initial values of N 
# to small numbers for immigrant and emigrants
# because preliminary analyses suggested low rates
# This is necessary to run model and avoid initial values error.
inits <- function(){list(
  # fecundity inits
  lmu.f.mean = runif(1, -2, 2),
  beta = runif(1, -2, 2),
  mu.nu = runif(1), 
  # survival
  lmu = runif(9), 
  Ustar = Ustar
  )} 

params <- c( # fecundity
            "lmu.brood", "delta", "sig.brood", "sig",
            "lmu.nest", "mu.nest", "gamma", "sig.nest",
            "eps.brood", "eps.nest",
            "mn.brood", "mn.nest", "mn.f", 
            "lam", 
            "f.dmape.obs", "f.dmape.rep",
            "f.tvm.obs", "f.tvm.rep",
             # survival 
             "mu", "lmu", 
             # abundance
             "N", "NB", "NF", "NFY", "NFYW", "NH", "Ntot",
             # error terms
             "eps", "sds", "Ustar", "U", "R",
             # yearly summaries
             'eta.phiFY', 'eta.phiA', 'eta.phiB',
             'eta.psiNB', 'eta.psiAB', 'eta.psiBA',
             'eta.pA', 'eta.pB',
             # goodness of fit
             "dmape.obs", "dmape.rep", 
             "tvm.obs", "tvm.rep",
             "tturn.obs", "tturn.rep",
             "f.dmape.obs", "f.dmape.rep",
             "f.tvm.obs", "f.tvm.rep"
)

n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000

mod <- nimbleModel(code, calculate=T, 
                   constants = constl, 
                   data = datl, 
                   inits = inits(), 
                   buildDerivs = TRUE)

mod$calculate()

cmod <- compileNimble(mod)
confhmc <- configureMCMC(mod)
#confhmc <- configureHMC(mod)
#confhmc$addSampler(target = 'mu.nu', type = 'RW')
#confhmc$addSampler(target = 'sigma.nu', type = 'RW')

confhmc$setMonitors(params)
hmc <- buildMCMC(confhmc)
chmc <- compileNimble(hmc, project = mod, resetFunctions = TRUE)

post <- runMCMC(chmc,
                niter = n.iter, 
                nburnin = n.burnin,
                nchains = n.chains,
                thin = n.thin,
                samplesAsCodaMCMC = T)

return(post)
} # run_ipm function end

this_cluster <- makeCluster(4)
post <- parLapply(cl = this_cluster, 
                       X = 1:4, 
                       fun = run_ipm, 
                       dat = datl, const = constl)
stopCluster(this_cluster)



#save(file=paste("./", m, ".Rdata", sep=""), list="out")