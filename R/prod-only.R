load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_dd_longrun_2025_Apr_03.rdata")
library('nimble')
library('parallel')
library ('coda')
library (MCMCvis)
load("data/data-dd.rdata")
outp <- MCMCpstr(post, type="chains")
# Simulate the Number of breeders because IPM estimates this
# eta <- array( rnorm(constl$nyr*2, mean=0, sd=2), dim=c(2, constl$ny) )
# mns <- apply(datl$num.treat, 2, mean)
# NB <- array(NA, c(constl$nyr, 2))
# 
# for (s in 1:2){
#   for (t in 1:constl$nyr){
# NB[t, s] <- rpois(1, mns[s] + eta[t])
#   }}
datl$NB <- apply(outp$NB, c(1,2), mean)[-1,] |> ceiling()
datl$NB[c(6,10),1] <- c(110,99)

mycode <- nimbleCode(
  {
for (s in 1:nsite){
  lmu.prod[s] ~ dnorm(0, sd=5)
  lmu.treat[s] <- logit(mu.treat[s])
  mu.treat[s] ~ dbeta(1,1)
} # s
gamma ~ dnorm(0, sd=10)
rr ~ dexp(0.05)
sd.treat ~ dexp(1)
sd.eta ~ dexp(1) 

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
        gamma*c(0,1)[tr] + #alphas[6, s]*(Ntot[t, s]-mnC[s]) + 
        eta[s, t]
    } # tr
    eta[s, t] ~ dnorm( 0, sd = sd.eta)
    num.treat[t, s] ~ dbinom( ptreat[t, s], NB[t,s] )
    logit(ptreat[t, s]) ~ dnorm( lmu.treat[s], sd = sd.treat)
  }} # s t
}
) # end model code

params <- c(
  # fecundity
  "lmu.prod", "gamma", "rr",  
  "mu.treat", "sd.treat", "sd.eta", 
  "mn.prod", "ptreat"
)

inits.func <- function (){
  # sample from working inits
  list(
    # fecundity
    lmu.prod = runif(2, -2, 2),
    gamma = runif(1, -1, 1),
    rr = runif(1, 40, 80),
    mu.treat = runif(2, 0.4, 0.6),
    sd.treat = runif(1),
    sd.eta = runif(1)
  )}
inits <- inits.func()

mod <- nimbleModel(mycode, 
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
nc=3; nt=10; ni=20000; nb=15000
#nc=1; nt=100; ni2=250000; nb=150000
#nc=1; nt=1; ni=3; nb=1

post2 <- runMCMC(chmc,
                 niter = ni,
                 nburnin = nb,
                 nchains = nc,
                 thin = nt,
                 samplesAsCodaMCMC = T,
                 setSeed = par_info[[1]]$seed,
                 inits = inits)
MCMCtrace(post2, params[-c(7,8)], pdf=F, Rhat=TRUE)

outp2 <- MCMCpstr(post2, type="chains", exact=FALSE, ISB=TRUE)
out2 <- do.call(rbind, post2)
MCMCplot(out2, params="ptreat")
