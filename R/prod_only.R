library('nimble')
library('parallel')
library ('coda')
#load("/bsuscratch/brianrolek/riha_ipm/data-dd.rdata")
load("data/data.rdata")
cpus <- 5

# Add data augmentation to productivity
# Add 200 unobserved breeders for each site/year
len <- 200*constl$nsite*constl$nyr
prod.da <- rep(NA, len)
site.da <- rep_len(1:constl$nsite, length.out=len)
yr.da <- rep_len(1:constl$nyr, length.out=len)
trt.da <- rep_len(0, length.out=len)

datl$prod <- c(datl$prod, prod.da)
constl$site.pair <- c(constl$site.pair, site.da)
constl$year.pair <- c(constl$year.pair, yr.da)
constl$treat.pair <- c(constl$treat.pair, trt.da)

for (ii in 1:length(par_info)){
  par_info[[ii]]$inits$prod <- c(rep(NA, constl$npairsobs), 
                                 rep(0, len))
}
constl$npairstot <- length(datl$prod) 

# Create a mtrix to deal with nonconsecutive indices
constl$pair.end <- table(constl$year.pair, constl$site.pair)
yrind.pair <- array(NA, dim = c(max(constl$pair.end), constl$nyr, constl$nsite) )
for (t in 1:(constl$nyr-1)){
  for (s in 1:(constl$nsite)){
    if(constl$pair.end[t, s]==0){next} else{
      yrind.pair[1:constl$pair.end[t, s],t,s] <- which(constl$year.pair==t & constl$site.pair==s)
    } # else
  }} # s t
constl$yrind.pair <- yrind.pair

datl$Ntot <- matrix( rpois(constl$nyr*constl$nsite, lambda=50), nrow=constl$nyr )
datl$NB <- matrix( rpois(constl$nyr*constl$nsite, lambda= 25), nrow=constl$nyr )

str(par_info[[1]]$inits)
par_info[[1]]$inits$alphas <- runif(2, -1, 1)

#**********************
#* Model code
#**********************
mycode <- nimbleCode(
  {
    ###################################################
    # Priors and constraints
    ##################################################
    
    ###############################
    # Likelihood for productivity
    ###############################
    # Priors
    for (s in 1:nsite){
      lmu.prod[s] ~ dnorm(0, sd=5)
    } # s
    gamma ~ dnorm(0, sd=10)
    alphas[1] ~ dnorm(0, sd=10)
    alphas[2] ~ dnorm(0, sd=10)
    rr ~ dexp(0.05)
    
    # Productivity likelihood      
    for (k in 1:npairstot){
      prod[k] ~ dnegbin(ppp[k], rr)
      ppp[k] <- rr/(rr+mu.prod[k])
      log(mu.prod[k]) <- lmu.prod[site.pair[k]] +  
        gamma*treat.pair[k] +
        alphas[site.pair[k]]*(Ntot[year.pair[k], site.pair[k]]-mnC[site.pair[k]]) 
    } # k
    # Derive yearly productivity for population model
    # need to reorder because nimble doesn't 
    # handle nonconsecutive indices
    # yrind.pair is a matrix of indices for each site
    for (t in 1:(nyr-1)){ # this index for dens dep ipm bc mn.prod[t-1] would cause a zero index
      for (s in 1:nsite){
        for (xxx in 1:pair.end[t,s]){ # reorder to deal with nonconsecutive indices
          prodmat[xxx,t,s] <- mu.prod[ yrind.pair[xxx,t,s] ]
        } # xxx
        mn.prod[t, s] <- mean.dyn( mat= prodmat[1:358 , t, s], 
                                   end=NB[t,s] )
      }} # s t
    
    # # GOF for productivity
    # for (k in 1:npairsobs){
    #   f.obs[k] <- prod[k] # observed counts
    #   f.exp[k] <- mu.prod[k] # expected counts adult breeder
    #   f.rep[k] ~ dnegbin(ppp[k], rr) # expected counts
    #   f.dssm.obs[k] <- abs( ( f.obs[k] - f.exp[k] ) / (f.obs[k]+0.001) )
    #   f.dssm.rep[k] <- abs( ( f.rep[k] - f.exp[k] ) / (f.rep[k]+0.001) )
    # } # k
    # f.dmape.obs <- sum(f.dssm.obs[1:npairsobs])
    # f.dmape.rep <- sum(f.dssm.rep[1:npairsobs])
    # f.tvm.obs <- sd(brood[1:npairsobs])^2/mean(brood[1:npairsobs])
    # f.tvm.rep <- sd(f.rep[1:npairsobs])^2/mean(f.rep[1:npairsobs])
    
})

#**********************
#* Function to run model in NIMBLE
#**********************
  
  params <- c(
    # pop growth 
    "lambda",
    # fecundity
    "lmu.prod", "gamma", "rr", "mn.prod", 
    # survival 
    "mus", "lmus", "betas", "deltas", "alphas",
    # abundance
    "NB", "NF", "NFY", "N", "NAD",
    "r",
    "N", "Ntot",
    # error terms
    "eta", "sds", "Ustar", "R", "z.score",
    # yearly summaries
    'mn.phiFY', 'mn.phiA', 'mn.phiB',
    'mn.psiFYB', 'mn.psiAB', 'mn.psiBA',
    'mn.pA', 'mn.pB',
    "prod.meas", "prod.unmeas",
    # goodness of fit
    "f.dmape.obs", "f.dmape.rep",
    "f.tvm.obs", "f.tvm.rep",
    "dmape.obs", "dmape.rep",
    "tvm.obs", "tvm.rep",
    "tturn.obs", "tturn.rep"
  )
  
  n.chains=1; n.thin=1; n.iter=100; n.burnin=50
  #n.chains=1; n.thin=100; n.iter=100000; n.burnin=50000
  #n.chains=1; n.thin=100; n.iter=250000; n.burnin=150000
  # create a nimble function to handle dynamic indices
  mean.dyn <- nimbleFunction(
    run = function(mat = double(), end = integer(1)) {
      return( mean(mat[1:end]) )
      returnType(double(1))
    })
  assign('mean.dyn', mean.dyn, envir = .GlobalEnv)
  
  mod <- nimbleModel(mycode, 
                     constants = constl, 
                     data = datl, 
                     inits = par_info[[1]]$inits, 
                     buildDerivs = FALSE, # doesn't work when TRUE, no hope for HMC
                     calculate=T 
  ) 
  
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
