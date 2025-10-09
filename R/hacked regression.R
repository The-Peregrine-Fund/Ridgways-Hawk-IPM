## ---- ipm --------
library('nimble')
library('parallel')
library ('coda')
library ('nimbleEcology')
library ('MCMCvis')
load("data/data-dd.rdata")
cpus <- 4

datlist <- list(hacked= constl$hacked)
constlist <- list(
  nyr=constl$nyr,
  nsite=constl$nsite,
  nind=constl$nind,
  first=constl$first,
  site=constl$site
)


mycode <- nimbleCode(
  {

# Probability of an individual at a site
# has a history of being hacked
lmu.hacked.lh <- logit(mu.hacked.lh)
mu.hacked.lh ~ dbeta(1,1)
for (h in 1:nhacked.lh){
  hacked.lh[h] ~ dbern( phacked.lh[ lh.year[h] ] )
} # h
  for (t in 1:nyr){
    logit(phacked.lh[t]) <- lmu.hacked.lh #+ eta[9, s, t]
  } # t

lmu.hacked.pc <- logit(mu.hacked.pc)
mu.hacked.pc ~ dbeta(1,1)
for (hh in 1:nhacked.pc){
  hacked.pc[hh] ~ dbern( phacked.pc[ pc.year[hh] ] )
} # h
for (t in 1:nyr){
  logit(phacked.pc[t]) <- lmu.hacked.pc #+ eta[9, s, t]
} # t

  }
)

run_ipm <- function(info, datl, constl, code){
library('nimble')
library('coda')
library ('nimbleEcology')
library ('MCMCvis')

params <- c(
  "mu.hacked.lh", "lmu.hacked.lh", "phacked.lh",
  "mu.hacked.pc", "lmu.hacked.pc", "phacked.pc"
)

# initial values  
inits.func <- function (){

  list(
mu.hacked = runif(2)
  )}
inits2 <- inits.func()

mod <- nimbleModel(code, 
                   constants = constl, 
                   data = datl, 
                   inits = inits2, 
                   buildDerivs = FALSE, # doesn't work when TRUE, no hope for HMC
                   calculate=T ) 
cmod <- compileNimble(mod, showCompilerOutput = TRUE)
confhmc <- configureMCMC(cmod)
confhmc$setMonitors(params)
hmc <- buildMCMC(confhmc)
chmc <- compileNimble(hmc, project = cmod, 
                      resetFunctions = TRUE,
                      showCompilerOutput = TRUE)
#nc=1; nt=100; ni=50000; nb=25000
#nc=1; nt=100; ni=100000; nb=50000
#nc=1; nt=100; ni2=250000; nb=150000
#nc=1; nt=1; ni=3; nb=1
nc=1; nt=1; ni=10000; nb=5000
#nc=1; nt=1; ni=200; nb=0

post <- runMCMC(chmc,
                niter = ni,
                nburnin = nb,
                nchains = nc,
                thin = nt,
                samplesAsCodaMCMC = T,
                setSeed = info$seed)

return(post)
} # run_ipm function end

# #*****************
# #* Run chains in parallel
# #*****************
this_cluster <- makeCluster(cpus)
post <- parLapply(cl = this_cluster,
                  X = par_info[1:cpus],
                  fun = run_ipm,
                  datl = datlist,
                  constl = constlist,
                  code = mycode)
stopCluster(this_cluster)

save(post, mycode, 
     file = flnm2)