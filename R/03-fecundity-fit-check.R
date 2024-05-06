library('nimble')
library('nimbleHMC')
library('MCMCvis')
library('parallel')
library ('coda')

datl <- list( # productivity data
              f = dfp$'# Fledged',
              z = ifelse(dfp$'# Fledged'>0, 1, NA)
)

constl <- list( nyr = nyr,
                K = nrow(dfp),
                year = as.numeric(factor(dfp$Year)), 
                tr = dfp$tr
)



run_p <- function(seed, datl, constl){
  library('nimble')
  library('coda')
  library ('nimbleHMC')
  
  code <- nimbleCode(
    {
      ###############################
      # Likelihood for productivity
      ###############################
      lmu.f.mean ~ dnorm(0, sd=2)
      sigma.f ~ dexp(1)
      beta ~ dnorm(0, sd=2) 
      for (k in 1:K){
        f[k] ~ T(dpois(mu.f[k]),,3)
        log(mu.f[k]) <- lmu.f[year[k]] + beta*tr[k] 
      } # k
      
      for (t in 1:nyr){
        lmu.f[t] <- lmu.f.mean + eps.f[t]
        eps.f[t] ~ dnorm(0, sd=sigma.f)
      } # t
      
      # GOF fecundity- Mean absolute percentage
      for (k in 1:K){
        f.obs[k] <- f[k] # observed counts
        f.exp[k] <- mu.f[k] # expected counts adult breeder
        f.rep[k] ~ dpois( mu.f[k] ) # expected counts
        f.dssm.obs[k] <- abs( ( (f.obs[k]) - (f.exp[k]) ) / (f.obs[k]+0.001) )
        f.dssm.rep[k] <- abs( ( (f.rep[k]) - (f.exp[k]) ) / (f.rep[k]+0.001) )
        
      } # k
      dmape.obs <- sum(f.dssm.obs[1:K])
      dmape.rep <- sum(f.dssm.rep[1:K])
      tvm.obs <- sd(f[1:K])^2/mean(f[1:K])
      tvm.rep <- sd(f.rep[1:K])^2/mean(f.rep[1:K])
    }) # model end
  
  params <- c( "lmu.f.mean", "sigma.f", "mu.f", "beta",
               "eps.f", "lmu.f",
               "dmape.obs", "dmape.rep",
               "tvm.obs", "tvm.rep"
  )
  
  inits <- function(){ list(lmu.f.mean = runif(1, -2, 2),
                            beta = runif(1, -2, 2), 
                            sigma.f = runif(1) )}
  
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
} # model function end

run_zip_trunc <- function(seed, datl, constl){
  library('nimble')
  library('coda')
  library ('nimbleHMC')
  
  code <- nimbleCode(
    {
      ###############################
      # Likelihood for productivity
      ###############################
      lmu.f.mean ~ dnorm(0, sd=2)
      beta ~ dnorm(0, sd=2)
      sigma.f ~ dexp(1)
      logit(lmu.nu) <- mu.nu
      mu.nu ~ dunif(0,1)
      sigma.nu ~ dexp(1)
      
      for (k in 1:K){
        f[k] ~ T(dpois(z[k]*mu.f[k]),,3)
        log(mu.f[k]) <- lmu.f[year[k]] + beta*tr[k] 
        z[k] ~ dbern(nu[year[k]])
      }
      
      for (t in 1:nyr){
        logit(nu[t]) <- lmu.nu + eps.nu[t]
        eps.nu[t] ~ dnorm(0, sd=sigma.nu)
        lmu.f[t] <- lmu.f.mean + eps.f[t]
        eps.f[t] ~ dnorm(0, sd=sigma.f)
      }
      
      # GOF fecundity
      for (k in 1:K){
        f.obs[k] <- f[k] # observed counts
        f.exp[k] <- z[k]*mu.f[k] # expected counts adult breeder
        f.rep[k] ~ dpois( mu.f[k] * nu[year[k]] ) # expected counts
        dssm.obs[k] <- abs( ( (f.obs[k]) - (f.exp[k]) ) / (f.obs[k]+0.001) )
        dssm.rep[k] <- abs( ( (f.rep[k]) - (f.exp[k]) ) / (f.rep[k]+0.001) )
        
      } # k
      dmape.obs <- sum(dssm.obs[1:K])
      dmape.rep <- sum(dssm.rep[1:K])
      tvm.obs <- sd(f[1:K])^2/mean(f[1:K])
      tvm.rep <- sd(f.rep[1:K])^2/mean(f.rep[1:K])
    }) # model end
  
  params <- c( "lmu.f.mean", "sigma.f", "mu.f","lmu.f", "beta",
               "lmu.nu", "sigma.nu", "mu.nu", "nu", 
               "eps.f", "eps.nu",
               "dmape.obs", "dmape.rep",
               "tvm.obs", "tvm.rep"
  )
  
  inits <- function(){ list(lmu.f.mean = runif(1, -2, 2),
                            beta = runif(1, -2, 2),
                            mu.nu = runif(1), 
                            sigma.f = runif(1),
                            sigma.nu = runif(1)) }
  
  n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
  
  mod <- nimbleModel(code, calculate=T, 
                     constants = constl, 
                     data = datl, 
                     inits = inits(), 
                     buildDerivs = TRUE)
  
  mod$calculate()
  
  cmod <- compileNimble(mod)
  #confhmc <- configureMCMC(mod)
  confhmc <- configureHMC(mod)
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
} # model function end

run_zip_hmc <- function(seed, datl, constl){
  library('nimble')
  library('coda')
  library ('nimbleHMC')
  
code <- nimbleCode(
    {
###############################
# Likelihood for productivity
###############################
lmu.f.mean ~ dnorm(0, sd=2)
beta ~ dnorm(0, sd=2)
sigma.f ~ dexp(1)
logit(lmu.nu) <- mu.nu
mu.nu ~ dunif(0,1)
sigma.nu ~ dexp(1)

for (k in 1:K){
  f[k] ~ dpois(z[k]*mu.f[k])
  log(mu.f[k]) <- lmu.f[year[k]] + beta*tr[k] 
  z[k] ~ dbern(nu[year[k]])
}

for (t in 1:nyr){
  logit(nu[t]) <- lmu.nu + eps.nu[t]
  eps.nu[t] ~ dnorm(0, sd=sigma.nu)
  lmu.f[t] <- lmu.f.mean + eps.f[t]
  eps.f[t] ~ dnorm(0, sd=sigma.f)
}

# GOF fecundity
for (k in 1:K){
  f.obs[k] <- f[k] # observed counts
  f.exp[k] <- z[k]*mu.f[k] # expected counts adult breeder
  f.rep[k] ~ dpois( mu.f[k] * nu[year[k]] ) # expected counts
  dssm.obs[k] <- abs( ( (f.obs[k]) - (f.exp[k]) ) / (f.obs[k]+0.001) )
  dssm.rep[k] <- abs( ( (f.rep[k]) - (f.exp[k]) ) / (f.rep[k]+0.001) )

} # k
dmape.obs <- sum(dssm.obs[1:K])
dmape.rep <- sum(dssm.rep[1:K])
tvm.obs <- sd(f[1:K])^2/mean(f[1:K])
tvm.rep <- sd(f.rep[1:K])^2/mean(f.rep[1:K])
}) # model end
  
  params <- c( "lmu.f.mean", "sigma.f", "mu.f","lmu.f", "beta",
               "lmu.nu", "sigma.nu", "mu.nu", "nu", 
               "eps.f", "eps.nu",
               "dmape.obs", "dmape.rep",
               "tvm.obs", "tvm.rep"
  )
  
  inits <- function(){ list(lmu.f.mean = runif(1, -2, 2),
                            beta = runif(1, -2, 2),
                            mu.nu = runif(1), 
                            sigma.f = runif(1),
                            sigma.nu = runif(1)) }
  
  n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
  
  mod <- nimbleModel(code, calculate=T, 
                     constants = constl, 
                     data = datl, 
                     inits = inits(), 
                     buildDerivs = TRUE)
  
  mod$calculate()
  
  cmod <- compileNimble(mod)
  #confhmc <- configureMCMC(mod)
  confhmc <- configureHMC(mod)
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
} # model function end


run_zip_hmc_marginalized <- function(seed, datl, constl){
  library('nimble')
  library('coda')
  library ('nimbleHMC')
  
  dZIP <- nimbleFunction(
    run = function(x = double(), lambda = double(),
                   zeroProb = double(), log = logical(0, default = 0)) {
      returnType(double())
      ## First handle non-zero data
      prob <- zeroProb * dbinom(x, size = 1, prob = 0) + (1 - zeroProb) * dpois(x, lambda)
      if (log) return(log(prob))
      return(prob)
    },
    buildDerivs = 'run'
  )
  assign('dZIP', dZIP, envir = .GlobalEnv)
  
  rZIP <- nimbleFunction(
    run = function(n = integer(), lambda = double(), zeroProb = double()) {
      returnType(double())
      isStructuralZero <- rbinom(1, prob = zeroProb, size = 1)
      if (isStructuralZero) return(0)
      return(rpois(1, lambda))
    })
  assign('rZIP', rZIP, envir = .GlobalEnv)
  
  registerDistributions(list(
    dZIP = list(
      BUGSdist = "dZIP(lambda, zeroProb)",
      discrete = TRUE,
      range = c(0, Inf),
      types = c('value = double()', 'lambda = double()', 'zeroProb = double()')
    )))
  
  code <- nimbleCode(
    {
      ###############################
      # Likelihood for productivity
      ###############################
      lmu.f.mean ~ dnorm(0, sd=2)
      beta ~ dnorm(0, sd=2)
      sigma.f ~ dexp(1)
      logit(lmu.nu) <- mu.nu
      mu.nu ~ dunif(0,1)
      sigma.nu ~ dexp(1)
      
      for (k in 1:K){
        f[k] ~ dZIP(mu.f[k], zeroProb = 1-nu[year[k]])
        log(mu.f[k]) <- lmu.f[year[k]] + beta*tr[k] 
      }
      
      for (t in 1:nyr){
        logit(nu[t]) <- lmu.nu + eps.nu[t]
        eps.nu[t] ~ dnorm(0, sd=sigma.nu)
        lmu.f[t] <- lmu.f.mean + eps.f[t]
        eps.f[t] ~ dnorm(0, sd=sigma.f)
      }
      
      # GOF fecundity
      # for (k in 1:K){
      #   f.obs[k] <- f[k] # observed counts
      #   f.exp[k] <- z[k]*mu.f[k] # expected counts adult breeder
      #   f.rep[k] ~ dpois( mu.f[k] * nu[year[k]] ) # expected counts
      #   dssm.obs[k] <- abs( ( (f.obs[k]) - (f.exp[k]) ) / (f.obs[k]+0.001) )
      #   dssm.rep[k] <- abs( ( (f.rep[k]) - (f.exp[k]) ) / (f.rep[k]+0.001) )
      #   
      # } # k
      # dmape.obs <- sum(dssm.obs[1:K])
      # dmape.rep <- sum(dssm.rep[1:K])
      # tvm.obs <- sd(f[1:K])^2/mean(f[1:K])
      # tvm.rep <- sd(f.rep[1:K])^2/mean(f.rep[1:K])
    }) # model end
  
  params <- c( "lmu.f.mean", "sigma.f", "mu.f","lmu.f", "beta",
               "lmu.nu", "sigma.nu", "mu.nu", "nu", 
               "eps.f", "eps.nu"
               # "dmape.obs", "dmape.rep",
               # "tvm.obs", "tvm.rep"
               #"dd.obs", "dd.rep"
  )
  
  inits <- function(){ list(lmu.f.mean = runif(1, -2, 2),
                            beta = runif(1, -2, 2),
                            mu.nu = runif(1), 
                            sigma.f = runif(1),
                            sigma.nu = runif(1) )}
                            #z= rbinom(length(datl$f), size=1, prob=0.9)) }
  
  n.chains=1; n.thin=25; n.iter=50000; n.burnin=25000
  
  mod <- nimbleModel(code, calculate=T, 
                     constants = constl, 
                     data = datl, 
                     inits = inits(), 
                     buildDerivs = TRUE)
  
  mod$calculate()
  
  cmod <- compileNimble(mod)
  #confhmc <- configureMCMC(mod)
  confhmc <- configureHMC(mod)
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
} # model function end


this_cluster <- makeCluster(4)
post <- list()
post[[1]] <- parLapply(cl = this_cluster, 
                  X = 1:4, 
                  fun = run_p, 
                  dat = datl, const = constl)

post[[2]] <- parLapply(cl = this_cluster, 
                       X = 1:4, 
                       fun = run_zip_trunc, 
                       dat = datl, const = constl)

post[[3]] <- parLapply(cl = this_cluster, 
                  X = 1:4, 
                  fun = run_zip_hmc, 
                  dat = datl, const = constl)

post[[4]] <- parLapply(cl = this_cluster, 
                       X = 1:4, 
                       fun = run_zip_hmc_marginalized, 
                       dat = datl, const = constl)

stopCluster(this_cluster)

pars <- c( "lmu.f.mean", "sigma.f", "beta",
           "lmu.nu", "mu.nu","sigma.nu", "nu")

MCMCsummary(post[[1]], digits=2, pars)
MCMCsummary(post[[2]], digits=2, pars)

MCMCtrace(post[[1]], pars, pdf=F)
MCMCtrace(post[[2]], pars, pdf=F)
MCMCtrace(post[[3]], pars, pdf=F)


extr_fun <- function(x) {list(as.mcmc(x[[1]]), 
                              as.mcmc(x[[2]]), 
                              as.mcmc(x[[3]]),
                              as.mcmc(x[[4]]))}
out <- list()
out[[1]] <- extr_fun(post[[1]])
out[[2]] <- extr_fun(post[[2]])

out[[1]] <- do.call(rbind, post[[1]])
MCMCplot( out[[1]] , pars)
MCMCplot( MCMCpstr(out[[2]][[1]], type = 'chains') , pars)

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
       main=paste0("Mean absolute percentage error\n",lab),
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

plot.diag(post[[1]], ratio=T, lab="Poisson")
plot.diag(post[[2]], lab="ZIP")

save(post, run_f,  
     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\fec-fit.Rdata")
