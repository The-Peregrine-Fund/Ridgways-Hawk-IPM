#******************
#* Check goodness-of-fit 
#* for population count data
#******************
## ---- start --------
library ('nimble')
library ('nimbleHMC')
library ('MCMCvis')
library ('coda')
library('parallel')
load("data/data.Rdata")
set.seed(5757575)

## ---- models --------
#***********************
#* Create models with 3 distributions
#* poisson, negative binomial, and zero-inflated Poisson
#***********************

### ---- pois --------
#***********************
#* Poisson model
#***********************
run_p <- function(seed, datl, constl){
  library('nimble')
  library('coda')
  library ('nimbleHMC')
  code <- nimbleCode(
    {
      # priors
      mu ~ dnorm(0, sd=5) # constrain to reasonable values <exp(5) or <148
      # likelihood
      for (t in 1:ntime){
        count[t] ~ dpois(lambda[t])
        # abundance
        log(lambda[t]) <- mu 
      } # t
      # Assess GOF of the models for counts
      # Step 1: Compute statistic for observed data
      # Step 2: Use discrepancy measure: mean absolute error
      # Step 3: Use test statistic: number of turns
      for (t in 1:ntime){
        c.obs[t] <- count[t] # observed counts
        c.exp[t] <- lambda[t] # expected counts adult breeder
        c.rep[t] ~ dpois(lambda[t]) # expected counts
        # Compute fit statistics, Mean absolute error
        dssm.obs[t] <- abs( ( (c.obs[t]) - (c.exp[t]) ) / (c.obs[t]+0.001) )
        dssm.rep[t] <- abs( ( (c.rep[t]) - (c.exp[t]) ) / (c.rep[t]+0.001) )
      } # t
      dmape.obs <- sum(dssm.obs[1:ntime])
      dmape.rep <- sum(dssm.rep[1:ntime])
      # variance-mean ratio
      tvm.rep <- sd(c.rep[1:ntime])^2/mean(c.rep[1:ntime])
      tvm.obs <- sd(y[1:ntime])^2/mean(y[1:ntime])
    }
  ) # end model
  
  params <- c(    "mu",  
                  "dssm.obs", "dssm.rep",
                  "dmape.obs", "dmape.rep",
                  "tvm.rep", "tvm.obs"
  )
  
  # create initial values for missing y data
  c.inits <- array(NA, dim(datl$count))
  c.inits[is.na(datl$count)] <- rpois(sum(is.na(datl$count)), 2)
  inits <- function(){ list (count = c.inits,
                             mu = runif(constl$nsite, -2, 2)
  )}
  n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
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


### ---- odpois --------
#***********************
#* Overdispersed poisson model
#***********************
run_odp <- function(seed, datl, constl){
  library('nimble')
  library('coda')
  library ('nimbleHMC')
  code <- nimbleCode(
    {
      # priors
      mu ~ dnorm(0, sd=5) # constrain to reasonable values <exp(5) or <148
      sigma.time ~ dexp(2)
      # likelihood
      for (t in 1:ntime){
        count[t] ~ dpois(lambda[t])
        # abundance
        log(lambda[t]) <- mu + eps.time[t] 
        eps.time[t] ~ dnorm(0, sd=sigma.time)
      } # t
      # Assess GOF of the models for counts
      # Step 1: Compute statistic for observed data
      # Step 2: Use discrepancy measure: mean absolute error
      # Step 3: Use test statistic: number of turns
      for (t in 1:ntime){
        c.obs[t] <- count[t] # observed counts
        c.exp[t] <- lambda[t] # expected counts adult breeder
        c.rep[t] ~ dpois(lambda[t]) # expected counts
        # Compute fit statistics, Mean absolute error
        dssm.obs[t] <- abs( ( (c.obs[t]) - (c.exp[t]) ) / (c.obs[t]+0.001) )
        dssm.rep[t] <- abs( ( (c.rep[t]) - (c.exp[t]) ) / (c.rep[t]+0.001) )
      } # t
      dmape.obs <- sum(dssm.obs[1:ntime])
      dmape.rep <- sum(dssm.rep[1:ntime])
      # variance-mean ratio
      tvm.rep <- sd(c.rep[1:ntime])^2/mean(c.rep[1:ntime])
      tvm.obs <- sd(y[1:ntime])^2/mean(y[1:ntime])
    }
  ) # end model
  
  params <- c(    "sigma.time",
                  "mu",  
                  "dssm.obs", "dssm.rep",
                  "dmape.obs", "dmape.rep",
                  "tvm.rep", "tvm.obs"
  )
  
  # create initial values for missing y data
  c.inits <- array(NA, dim(datl$count))
  c.inits[is.na(datl$count)] <- rpois(sum(is.na(datl$count)), 2)
  inits <- function(){ list (count = c.inits,
                             sigma.time = rexp(1),
                             mu = runif(constl$nsite, -2, 2)
  )}
  n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
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

### ---- nb --------
#***********************
#* Negative binomial model
#* Best fitting model
#***********************
run_nb <- function(seed, datl, constl){
  library('nimble')
  library('coda')
  library ('nimbleHMC')
  code <- nimbleCode(
    {
      # priors
      mu ~ dnorm(0, sd=5) # constrain to reasonable values <exp(5) or <148
      r ~ dexp(0.2) 
      # likelihood
      for (t in 1:ntime){
        p[t] <- r/(r+lambda[t])
        count[t] ~ dnegbin(p[t], r)
        # abundance
        log(lambda[t]) <- mu 
      } # t
      # Assess GOF of the models for counts
      # Step 1: Compute statistic for observed data
      # Step 2: Use discrepancy measure: mean absolute error
      # Step 3: Use test statistic: number of turns
      for (t in 1:ntime){
        c.obs[t] <- count[t] # observed counts
        c.exp[t] <- lambda[t] # expected counts adult breeder
        c.rep[t] ~ dnegbin(p[t],r) # expected counts
        # Compute fit statistics, Mean absolute error
        dssm.obs[t,j] <- abs( ( (c.obs[t]) - (c.exp[t]) ) / (c.obs[t]+0.001) )
        dssm.rep[t,j] <- abs( ( (c.rep[t]) - (c.exp[t]) ) / (c.rep[t]+0.001) )
      } # t
      dmape.obs <- sum(dssm.obs[1:ntime])
      dmape.rep <- sum(dssm.rep[1:ntime])
      # variance-mean ratio
      tvm.rep <- sd(c.rep[1:ntime])^2/mean(c.rep[1:ntime])
      tvm.obs <- sd(y[1:ntime])^2/mean(y[1:ntime])
    }
  ) # end model
  
  params <- c(    "mu", "r", "p",
                  "dssm.obs", "dssm.rep",
                  "dmape.obs", "dmape.rep",
                  "tvm.rep", "tvm.obs"
  )
  
  # create initial values for missing y data
  c.inits <- array(NA, dim(datl$count))
  c.inits[is.na(datl$count)] <- rpois(sum(is.na(datl$count)), 2)
  inits <- function(){ list (count = c.inits,
                             sigma.time = rexp(1),
                             mu = runif(1, -2, 2),
                             r = runif(1)
  )}
  n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE)
  
  mod$calculate()
  
  cmod <- compileNimble(mod)
  
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

### ---- zip --------
#***********************
#* Zero-inflated Poisson model
#***********************
run_zip <- function(seed, datl, constl){
  library('nimble')
  library('coda')
  library ('nimbleHMC')
  code <- nimbleCode(
    {
      # priors
        mu ~ dnorm(0, sd=5) # constrain to reasonable values <exp(5) or <148
        psi ~ dunif(0, 1) 
      # likelihood
      for (t in 1:ntime){
          count[t] ~ dpois(lambda.star[t])
          lambda.star[t] <- lambda[t]*z[t]
          z[t] ~ dbern(psi)
          # abundance
          log(lambda[t]) <- mu 
      } # t
      
      # Assess GOF of the models for counts
      # Step 1: Compute statistic for observed data
      # Step 2: Use discrepancy measure: mean absolute error
      # Step 3: Use test statistic: number of turns
      for (t in 1:ntime){
          c.obs[t] <- count[t] # observed counts
          c.exp[t] <- z[t]*lambda[t] # expected counts adult breeder
          c.rep[t] ~ dpois(lambda[t]*psi) # expected counts
          # Compute fit statistics, Mean absolute error
          dssm.obs[t] <- abs( ( (c.obs[t]) - (c.exp[t]) ) / (c.obs[t]+0.001) )
          dssm.rep[t] <- abs( ( (c.rep[t]) - (c.exp[t]) ) / (c.rep[t]+0.001) )
        } # t
      dmape.obs <- sum(dssm.obs[1:ntime])
      dmape.rep <- sum(dssm.rep[1:ntime])
      # variance-mean ratio
      tvm.rep <- sd(c.rep[1:ntime])^2/mean(c.rep[1:ntime])
      tvm.obs <- sd(y[1:ntime])^2/mean(y[1:ntime])
    }
  ) # end model
  
  params <- c(    "sigma.time",
                  "mu", "psi",
                  "dssm.obs", "dssm.rep",
                  "dmape.obs", "dmape.rep",
                  "tvm.rep", "tvm.obs"
  )
  
  # create initial values for missing y data
  c.inits <- array(NA, dim(datl$count))
  c.inits[is.na(datl$count)] <- rpois(sum(is.na(datl$count)), 2)
  inits <- function(){ list (count = c.inits,
                             mu = runif(constl$nsite, -2, 2),
                             psi = runif(constl$nsite)
  )}
  
  n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
                     buildDerivs = TRUE)
  
  mod$calculate()
  
  cmod <- compileNimble(mod )
  # Use MCMC rather than HMC 
  # because HMC samplers are stuck
  confhmc <- configureMCMC(mod)
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

### ---- norm --------
#***********************
#* Normal model
#***********************
run_norm <- function(seed, datl, constl){
  library('nimble')
  library('coda')
  library ('nimbleHMC')
  code <- nimbleCode(
    {
      # priors
      mu ~ dnorm(0, sd=5) # constrain to reasonable values <exp(5) or <148
      sigma ~ dexp(2)
      # likelihood
      for (t in 1:ntime){
        count[t] ~ dnorm(mu[t], sd=sigma)
      } # t
      # Assess GOF of the models for counts
      # Step 1: Compute statistic for observed data
      # Step 2: Use discrepancy measure: mean absolute error
      # Step 3: Use test statistic: number of turns
      for (t in 1:ntime){
        c.obs[t] <- count[t] # observed counts
        c.exp[t] <- mu[t] # expected counts adult breeder
        c.rep[t] ~ dnorm(mu[t], sigma) # expected counts
        # Compute fit statistics, Mean absolute error
        dssm.obs[t] <- abs( ( (c.obs[t]) - (c.exp[t]) ) / (c.obs[t]+0.001) )
        dssm.rep[t] <- abs( ( (c.rep[t]) - (c.exp[t]) ) / (c.rep[t]+0.001) )
      } # t
      dmape.obs <- sum(dssm.obs[1:ntime])
      dmape.rep <- sum(dssm.rep[1:ntime])
      # variance-mean ratio
      tvm.rep <- sd(c.rep[1:ntime])^2/mean(c.rep[1:ntime])
      tvm.obs <- sd(y[1:ntime])^2/mean(y[1:ntime])
    }
  ) # end model
  
  params <- c(    "sigma",
                  "mu",  
                  "dssm.obs", "dssm.rep",
                  "dmape.obs", "dmape.rep",
                  "tvm.rep", "tvm.obs"
  )
  
  # create initial values for missing y data
  c.inits <- array(NA, dim(datl$count))
  c.inits[is.na(datl$count)] <- rpois(sum(is.na(datl$count)), 2)
  inits <- function(){ list (count = c.inits,
                             sigma = rexp(1),
                             mu = runif(constl$nsite, -2, 2)
  )}
  n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
  
  mod <- nimbleModel(code, calculate=T, constants = constl, 
                     data = datl, inits = inits(), 
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

# ---- runmodels --------
#**********************
#* Run models
#**********************
# approximately XX minute runtime per model
this_cluster <- makeCluster(4)
postB <- postNB <- postFY <- list()
postB[[1]] <- parLapply(cl = this_cluster, 
                  X = 1:4, 
                  fun = run_p, 
                  dat = datl, 
                  const = constl)
postB[[2]] <- parLapply(cl = this_cluster, 
                        X = 1:4, 
                        fun = run_odp, 
                        dat = datl, 
                        const = constl)
postB[[3]] <- parLapply(cl = this_cluster, 
                       X = 1:4, 
                       fun = run_nb, 
                       dat = datl, 
                       const = constl)
postB[[4]] <- parLapply(cl = this_cluster, 
                       X = 1:4, 
                       fun = run_zip, 
                       dat = datl, 
                       const = constl)
postB[[5]] <- parLapply(cl = this_cluster, 
                        X = 1:4, 
                        fun = run_norm, 
                        dat = datl, 
                        const = constl)

postNB[[1]] <- parLapply(cl = this_cluster, 
                       X = 1:4, 
                       fun = run_p, 
                       dat = datl, 
                       const = constl)
postNB[[2]] <- parLapply(cl = this_cluster, 
                        X = 1:4, 
                        fun = run_odp, 
                        dat = datl, 
                        const = constl)
postNB[[3]] <- parLapply(cl = this_cluster, 
                       X = 1:4, 
                       fun = run_nb, 
                       dat = datl, 
                       const = constl)
postNB[[4]] <- parLapply(cl = this_cluster, 
                       X = 1:4, 
                       fun = run_zip, 
                       dat = datl, 
                       const = constl)
postNB[[5]] <- parLapply(cl = this_cluster, 
                        X = 1:4, 
                        fun = run_norm, 
                        dat = datl, 
                        const = constl)

postFY[[1]] <- parLapply(cl = this_cluster, 
                       X = 1:4, 
                       fun = run_p, 
                       dat = datl, 
                       const = constl)
postFY[[2]] <- parLapply(cl = this_cluster, 
                         X = 1:4, 
                         fun = run_odp, 
                         dat = datl, 
                         const = constl)
postFY[[3]] <- parLapply(cl = this_cluster, 
                       X = 1:4, 
                       fun = run_nb, 
                       dat = datl, 
                       const = constl)
postFY[[4]] <- parLapply(cl = this_cluster, 
                       X = 1:4, 
                       fun = run_zip, 
                       dat = datl, 
                       const = constl)
postFY[[5]] <- parLapply(cl = this_cluster, 
                         X = 1:4, 
                         fun = run_norm, 
                         dat = datl, 
                         const = constl)
stopCluster(this_cluster)

extr_fun <- function(x) {list(as.mcmc(x[[1]]), 
                             as.mcmc(x[[2]]), 
                             as.mcmc(x[[3]]),
                             as.mcmc(x[[4]]))}

outB <- lapply(postB, extr_fun)
outNB <- lapply(postNB, extr_fun)
outFY <- lapply(postFY, extr_fun)
save(postB, postNB, postFY, 
     outB, outNB, outFY, 
     run_nb, run_nb, run_zip, 
     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\GOF.Rdata")

# ---- PPCfunction --------
# Function for posterior predictive checks
# to assess goodness-of-fit
plot.diag <- function(out, ratio=FALSE, lab=""){
  par(mfrow=c(1,1))
  # plot mean absolute percentage error
  samps <- MCMCpstr(out, "all", type="chains")
  mx <- max(c(samps$dmape.rep, samps$dmape.obs))
  mn <- min(c(samps$dmape.rep, samps$dmape.obs))
  plot(jitter(samps$dmape.obs, amount=300), 
       jitter(samps$dmape.rep, amount=300),
       main=paste0("Mean absolute percentage error\nmodel\n",lab),
       ylab="Discrepancy replicate values",
       xlab="Discrepancy observed values", 
       xlim=c(mn,mx), ylim=c(mn,mx), 
       pch=16, cex=0.5, col="gray10")
  curve(1*x, from=mn, to=mx, add=T, lty=2, lwd=2, col="blue")
  bp1 <- round(mean(samps$dmape.rep > samps$dmape.obs),2)
  loc <- ifelse(bp1 < 0.5, "topleft", "bottomright")
  legend(loc, legend=bquote(p[B]==.(bp1)), bty="n", cex=2)
  
  if (ratio==TRUE){
    # plot variance/mean ratio
    hist(samps$tvm.rep, nclass=50,
         xlab="variance/mean ", main=NA, axes=FALSE)
    abline(v=samps$tvm.obs, col="red")
    axis(1); axis(2)
  }
  return(list('Bayesian p-value'=bp1))
}

## ---- plotppcs -----


## ---- ppc --------
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\GOF.Rdata")
pars_pois <-c( "sigma.time", "mu")
pars_nb <-c( "sigma.time", "mu", "r")
pars_zip <- c( "sigma.time", "psi", "mu")
# check convergence
MCMCtrace(nb, params_nb, pdf=F,
          ind = TRUE, Rhat = TRUE, n.eff = TRUE)
par(mfrow=c(1,1))
# plot point estimates
MCMCplot(object = nb, params = pars_nb)
# plot posterior predicitive checks
plot.diag(nb) # posterior predictive check