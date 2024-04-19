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
  # and this is estimated using the multivariate normal distribution
  p <- nimDim(eps)[1] # extract the number of intercepts/coefficients
  R <- t(Ustar)%*%Ustar # calculate rhos, correlation coefficients
  Ustar[1:p,1:p] ~ dlkj_corr_cholesky(eta=1, p=p) # Ustar is the Cholesky decomposition of the correlation matrix
  U[1:p,1:p] <- uppertri_mult_diag(Ustar[1:p, 1:p], sds[1:p])
  # multivariate normal for temporal variance
  for (t in 1:(nyr-1)){
    eps[1:p,t] ~ dmnorm(0, cholesky = U[1:p, 1:p], prec_param = 0)
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
  # Priors for fecundity
  lmu.f.mean ~ dnorm(0, sd=2)
  beta[1] ~ dnorm(0, sd=2)
  beta[2] ~ dnorm(0, sd=2)
  logit(lmu.nu) <- mu.nu
  mu.nu ~ dbeta(1,1) # conjugate to Bernoulli
  
  for (k in 1:K){
    f[k] ~ T(dpois(z[k]*mu.f[k]), ,3)
    log(mu.f[k]) <- lmu.f[year[k]] + beta[1]*treat[k] + beta[2]*num_transl[k] 
    z[k] ~ dbern(nu[year[k]])
  }
  
  for (t in 1:nyr){
    logit(nu[t]) <- lmu.nu + eps[10,t]
    lmu.f[t] <- lmu.f.mean + eps[11,t]
  }
    
  # GOF for fecundity
  for (k in 1:K){
    f.obs[k] <- f[k] # observed counts
    f.exp[k] <- z[k]*mu.f[k] # expected counts adult breeder
    f.rep[k] ~ dpois( mu.f[k] * nu[year[k]] ) # expected counts
    f.dssm.obs[k] <- abs( ( (f.obs[k]) - (f.exp[k]) ) / (f.obs[k]+0.001) )
    f.dssm.rep[k] <- abs( ( (f.rep[k]) - (f.exp[k]) ) / (f.rep[k]+0.001) )
  } # k
  f.dmape.obs <- sum(dssm.obs[1:K])
  f.dmape.rep <- sum(dssm.rep[1:K])
  f.tvm.obs <- sd(f[1:K])^2/mean(f[1:K])
  f.tvm.rep <- sd(f.rep[1:K])^2/mean(f.rep[1:K])
      
  ################################
  # Likelihood for counts
  ################################
  # Abundance for year=1
  for (v in 1:5){
  N[v,1] ~ dcat(pInit) 
  }
  # Abundance for years >2
  for (t in 2:(nyr-1)){
  # Number of wild born juvs
  N[1,t+1] ~ dpois(NFY[t]*phiFY[t]*fFY[1,t]/2 + 
                  NB[t]*phiFY[t]*fAD[2,t]/2)
  # Abundance of nonbreeders
  ## Nonbreeders to nonbreeders
  N[2,t+1] ~ dbin(phiNB[t]*(1-psiNB_B[t]), NF[t])
  ## Breeders to nonbreeders
  N[3,t+1] ~ dbin(phiB[t]*psiB_NB[t], NB[t])
  # Abundance of breeders
  ## Nonbreeders to breeders
  N[4,t+1] ~ dbin(phi_NB[t]*psi_NB_B[t], NF[t])
  ## Breeders to breeders
  N[5,t+1] ~ dbin(phi_B[t]*(1-psi_B_NB[t]), NB[t])
  } # t
  
  for (t in 1:nyr){
  Ntot[t] <- sum(N[c(1,2,3,4,5),t]) # total number
  NB[t] <- sum(N[c(4,5),t]) # number of breeders
  NF[t] <- sum(N[c(2,3),t])  # number of nonbreeders
  NFY[t] <- N[1,t] # Includes translocated first years
  NFYW[t] <- N[1,t] # Number of wild born first years
  } # t
  
  # Observation process    
  for (t in 2:nyr){
  countB[t] ~ dpois(NB[t]) # breeding males 
  countA[t] ~ dpois(NF[t]) # nonbreeding adult males
  countFYW[t] ~ dpois(NFYW[t]) # first year wild males, includes translocated/hacked
  } # t
  
  ###################
  # Assess GOF of the state-space models for counts
  # Step 1: Compute statistic for observed data
  # Step 2: Use discrepancy measure: mean absolute error
  # Step 3: Use test statistic: number of turns
  ###################
  for (t in 2:nyr){ 
  c.expB[t] <- NB[t] + 0.001 # expected counts adult breeder
  c.expA[t] <- NF[t] + 0.001 # nonbreeder
  c.expFYW[t] <- countFYW[t] + 0.001 # first year
  c.obsB[t] <- countB[t] + 0.001
  c.obsA[t] <- countA[t] + 0.001
  c.obsFYW[t] <- countNFYW[t] + 0.001
  dssm.obsB[t] <- abs( ( (c.obsB[t]) - (c.expB[t]) ) / (c.obsB[t]+0.001)  )
  dssm.obsA[t] <- abs( ( (c.obsA[t]) - (c.expA[t]) ) / (c.obsA[t]+0.001)  )
  dssm.obsFYW[t] <- abs( ( (c.obsFYW[t]) - (c.expFYW[t]) ) / (c.obsFYW[t]+0.001)  )
  } # t
  dmape.obs[1] <- sum(dssm.obsB[2:n.yr])
  dmape.obs[2] <- sum(dssm.obsA[2:n.yr])
  dmape.obs[3] <- sum(dssm.obsFYW[2:n.yr])
  # Compute fit statistic for replicate data
  # Mean absolute error
  for (t in 2:nyr){ 
  c.repB[t] ~ dpois( NB[t] ) # expected counts
  c.repA[t] ~ dpois( NF[t] ) 
  c.repFYW[t] ~ dpois( NFYW[1,t] )
  } # t
  for (t in 2:nyr){ 
  dssm.repB[t] <- abs( ( (c.repB[t]) - (c.expB[t]) ) / (c.repB[t]+0.001) )
  dssm.repA[t] <- abs( ( (c.repA[t]) - (c.expA[t]) ) / (c.repA[t]+0.001) )
  dssm.repFYW[t] <- abs( ( (c.repFYW[t]) - (c.expFYW[t]) ) / (c.repFYW[t]+0.001) )
  } # t
  dmape.rep[1] <- sum(dssm.repB[2:n.yr])
  dmape.rep[2] <- sum(dssm.repA[2:n.yr])
  dmape.rep[3] <- sum(dssm.repFYW[2:n.yr])
  
  # Test statistic for number of turns
  for (t in 2:(nyr-2)){
  tt1.obsB[t] <- step(countB[t+2] - countB[t+1])
  tt2.obsB[t] <- step(countB[t+1] - countB[t])
  tt3.obsB[t] <- equals(tt1.obsB[t] + tt2.obsB[t], 1)
  tt1.obsA[t] <- step(countA[t+2] - countA[t+1])
  tt2.obsA[t] <- step(countA[t+1] - countA[t])
  tt3.obsA[t] <- equals(tt1.obsA[t] + tt2.obsA[t], 1)
  tt1.obsFYW[t] <- step(countFYW[t+2] - countFYW[t+1])
  tt2.obsFYW[t] <- step(countFYW[t+1] - countFYW[t])
  tt3.obsFYW[t] <- equals(tt1.obsFYW[t] + tt2.obsFYW[t], 1)
  } # t
  tturn.obs[1] <- sum(tt3.obsB[2:(n.yr-2)])
  tturn.obs[2] <- sum(tt3.obsA[2:(n.yr-2)])
  tturn.obs[3] <- sum(tt3.obsFYW[2:(n.yr-2)])
  
  for (t in 2:(nyr-2)){
  tt1.repB[t] <- step(c.repB[t+2] - c.repB[t+1])
  tt2.repB[t] <- step(c.repB[t+1] - c.repB[t])
  tt3.repB[t] <- equals(tt1.repB[t] + tt2.repB[t], 1)
  tt1.repA[t] <- step(c.repA[t+2] - c.repA[t+1])
  tt2.repA[t] <- step(c.repA[t+1] - c.repA[t])
  tt3.repA[t] <- equals(tt1.repA[t] + tt2.repA[t], 1)
  tt1.repFYW[t] <- step(c.repFYW[t+2] - c.repFYW[t+1])
  tt2.repFYW[t] <- step(c.repFYW[t+1] - c.repFYW[t])
  tt3.repFYW[t] <- equals(tt1.repFYW[t] + tt2.repFYW[t], 1)
  } # t
  tturn.rep[1] <- sum(tt3.repB[2:(n.yr-2)])
  tturn.rep[2] <- sum(tt3.repA[2:(n.yr-2)])
  tturn.rep[3] <- sum(tt3.repFYW[2:(n.yr-2)])
  
  
  ################################
  # Likelihood for survival
  ################################
  # priors for means
  for (pp in 1:p){ # coefficient
    for (h in 1:2){ # translocated
        for (m in 1:4){ # population
          lmu[pp,h,m] <- logit(mu[pp,h])
          mu[pp,h,m] <- dbeta(1,1)
        }  } }   # m population #s sex #h hacked #t time
  
  for (i in 1:nind){
    for (t in 1:(nyr-1)){
      #Survival
      logit(phiFY[i,t]) <- lmu[1, transl[i], pop[i]] + eps[1,t]  # first year
      logit(phiA[i,t]) <- lmu[2, transl[i], pop[i]] + eps[2,t] # nonbreeder
      logit(phiB[i,t]) <- lmu[3, transl[i], pop[i]] + eps[3,t]  # breeder
      #Recruitment
      logit(psiFYB[i,t]) <- lmu[4, transl[i], pop[i]] + eps[4,t]  # first year to breeder
      logit(psiAB[i,t]) <- lmu[5, transl[i], pop[i]] + eps[5,t]  # nonbreeder to breeder
      logit(psiBA[i,t]) <- lmu[6, transl[i], pop[i]] + eps[6,t]  # breeder to nonbreeder
      #Re-encounter
      logit(pFY[i,t]) <- lmu[7, transl[i], pop[i]] + eps[7,t]  # resight of nonbreeders
      logit(pA[i,t]) <- lmu[8, transl[i], pop[i]] + eps[8,t]  # resight of nonbreeders
      logit(pB[i,t]) <- lmu[9, transl[i], pop[i]] + eps[9,t]  # resight of breeders
    }#t
  }#i
  
  # Define state-transition and observation matrices
  for (i in 1:nind){  
    # Define probabilities of state S(t+1) given S(t)
    for (t in first[i]:(nyr-1)){
      ps[1,i,t,1] <- 0
      ps[1,i,t,2] <- sFY[i,t] * (1-psiFYB[i,t]) 
      ps[1,i,t,3] <- sFY[i,t] * psiFYB[i,t]
      ps[1,i,t,4] <- (1-sFY[i,t])
      
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
      po[1,i,t,1] <- pFY[i,t] 
      po[1,i,t,2] <- 0
      po[1,i,t,3] <- 0
      po[1,i,t,4] <- 0
      
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- pA[i,t]
      po[2,i,t,3] <- 0
      po[2,i,t,4] <- 0
      
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 0
      po[3,i,t,3] <- pB[i,t]
      po[3,i,t,4] <- 0
      
      po[4,i,t,1] <- (1-pFY[i,t])
      po[4,i,t,2] <- (1-pA[i,t])
      po[4,i,t,3] <- (1-pB[i,t])
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
             "lmu.f.mean", "mu.f","lmu.f", "beta",
             "lmu.nu", "mu.nu", "nu", 
             # survival 
             "mu", "lmu", 
             # abundance
             "N", "NB", "NF", "NFY", "NFYW", "NH", "Ntot",
             # error terms
             "eps", "sds", "Ustar", "U", "R",
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