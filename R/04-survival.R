library ('nimble')
library ('nimbleHMC')
library ('MCMCvis')
library ('coda')
library('parallel')
load("data/data.Rdata")
set.seed(5757575)

run_surv <- function(seed, datl, constl){
  library('nimble')
  library('coda')
  library ('nimbleHMC')

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
  assign('uppertri_mult_diag', ifgreaterFun, envir = .GlobalEnv)
  
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
    #   Sex
    #     1 female
    #     2 male
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
################################
# Likelihood for survival
################################ 
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
# priors for means
      for (pp in 1:p){ # coefficient
        for (h in 1:2){ # translocated
          for (s in 1:2){ # sex
            for (m in 1:4){ # population
            lmu[pp,s,h,m] <- logit(mu[pp,s,h])
            mu[pp,s,h,m] <- dbeta(1,1)
          }  } } }  # m population #s sex #h hacked #t time
      
for (i in 1:nind){
  for (t in 1:(nyr-1)){
    #Survival
    logit(phiFY[i,t]) <- lmu[1, sex[i], transl[i], pop[i]] + eps[1,t]  # first year
    logit(phiA[i,t]) <- lmu[2, sex[i], transl[i], pop[i]] + eps[2,t] # nonbreeder
    logit(phiB[i,t]) <- lmu[3, sex[i], transl[i], pop[i]] + eps[3,t]  # breeder
    #Recruitment
    logit(psiFYB[i,t]) <- lmu[4, sex[i], transl[i], pop[i]] + eps[4,t]  # first year to breeder
    logit(psiAB[i,t]) <- lmu[5, sex[i], transl[i], pop[i]] + eps[5,t]  # nonbreeder to breeder
    logit(psiBA[i,t]) <- lmu[6, sex[i], transl[i], pop[i]] + eps[6,t]  # breeder to nonbreeder
    #Re-encounter
    logit(pFY[i,t]) <- lmu[7, sex[i], transl[i], pop[i]] + eps[7,t]  # resight of nonbreeders
    logit(pA[i,t]) <- lmu[8, sex[i], transl[i], pop[i]] + eps[8,t]  # resight of nonbreeders
    logit(pB[i,t]) <- lmu[9, sex[i], transl[i], pop[i]] + eps[9,t]  # resight of breeders
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
} ) #model
  
      params <- c( "mu", "lmu", "eps", 
                   "Ustar", "U", "R", "sds"
      )
      
      # create initial values for missing y data
      get.first <- function(x) min(which(x!=4))
      f <- apply(datl$y, 1, get.first)
      get.last <- function(x) max(which(x!=4))
      l <- apply(datl$y, 1, get.last)
      
      TFmat <- is.na(z.inits) & is.na(datl$z)
      for (i in 1:nrow(TFmat) ){  TFmat[i,1:f[i]] <- FALSE }
      z.inits[TFmat] <- sample(size=445, c(2,3), replace=T, prob=c(0.5, 0.5) ) 
      # create inits for rhos
      p <- 9
      Ustar <- array(runif(p*p, -0.5, 0.5), dim=c(p,p))
      diag(Ustar) <- 1 # set diagonal to 1
      Ustar[lower.tri(Ustar)] <- 0 # set lower diag to zero
      t(Ustar)%*%Ustar
      
      inits <- function(){ list(z = z.inits,
                                Ustar = Ustar) }
      
      n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
      
      mod <- nimbleModel(code, calculate=T, 
                         constants = constl, 
                         data = datl, 
                         inits = inits(), 
                         buildDerivs = TRUE)
      
      mod$calculate()
      
      cmod <- compileNimble(mod )
      
      confhmc <- configureHMC(mod)
      confhmc$setMonitors(params)
      hmc <- buildMCMC(confhmc)
      chmc <- compileNimble(hmc, project=mod, resetFunctions = TRUE)
      
      post <- runMCMC(chmc,
                      niter = n.iter, 
                      nburnin = n.burnin,
                      nchains = n.chains,
                      thin = n.thin,
                      samplesAsCodaMCMC = T)
      
      return(post)
    }

this_cluster <- makeCluster(4)
post <- parLapply(cl = this_cluster, 
                        X = 1:4, 
                        fun = run_surv, 
                        dat = datl, const = constl)
stopCluster(this_cluster)

extr_fun <- function(x) {list(as.mcmc(x[[1]]), 
                              as.mcmc(x[[2]]), 
                              as.mcmc(x[[3]]),
                              as.mcmc(x[[4]]))}

out <- extr_fun(post)
save(post, out, run_surv,  
     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\surv.Rdata")
