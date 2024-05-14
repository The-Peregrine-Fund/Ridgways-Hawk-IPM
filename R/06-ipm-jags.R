################
# JAGS version
###############
setwd("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\GitHub\\Ridgways-Hawk-IPM\\")
load("data\\data.rdata")
library (jagsUI)
#load("/bsuscratch/brianrolek/riha_ipm/data.rdata")

m<- c("ipm-jags")
modfl <- paste("./", m, ".txt", sep="")
sink(modfl)
cat("
    model{
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
   
    #######################
    # Derived params
    #######################
    # for (t in 1:(nyr-1)){    
    #   lambda[t] <-  Ntot[t+1]/(Ntot[t])
    #   loglambda[t] <- log(lambda[t])
    # } #t
    
    ###############################
    # Likelihood for fecundity
    ###############################
    # Priors for number of fledglings
    lmu.brood ~ dnorm(0, 1/(5*5))
    delta ~ dunif(-5, 5)
    sig.brood ~ dexp(1)
    sig.brood.t ~ dexp(1)
    
    # priors for nest success
    lmu.nest <- logit(mu.nest)
    mu.nest ~ dbeta(1, 1) 
    gamma ~ dunif(-5, 5)
    sig.nest ~ dexp(1)
    
    # Two models for fecundity, (1) brood size and (2) nest success  
    # Brood size
    for (k in 1:nbrood){
      brood[k] ~ dlnorm(lam[k], 1/(sig.brood*sig.brood) ) # truncating seems to break something
      lam[k] <- lmu.brood +
                delta*treat.brood[k] +
                eps.brood[year.brood[k] ]
    } # k
    # Nest success       
    for (n in 1:nnest){
      nest.success[n] ~ dbern( nu[n] )
      logit(nu[n]) <- lmu.nest + 
                      gamma*treat.nest[n] + 
                      eps.nest[year.nest[n] ]
    } # n 
    for ( t in 1:nyr){
      eps.brood[t] ~ dnorm(0, 1/(sig.brood.t*sig.brood.t) )
      eps.nest[t] ~ dnorm(0, 1/(sig.nest*sig.nest) )
    } # t
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
      f.rep[k] ~ dlnorm(lam[k], 1/(sig.brood*sig.brood) ) # expected counts
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
    # # Abundance for year=1
    for (v in 1:7){
      N[v,1] ~ dcat(pInit[1:88])
    }
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
    # sig.NFY ~ dexp(1)
    # sig.NF ~ dexp(1)
    # sig.NB ~ dexp(1)
    
    for (t in 1:nyr){
      # counts.marked[1, t] ~ dnorm(NFY[t], 1/(sig.NFY*sig.NFY) ) # first year males, includes translocated/hacked
      # counts.marked[2, t] ~ dnorm(NF[t], 1/(sig.NF*sig.NF) ) # nonbreeding adult males 
      # counts.marked[3, t] ~ dnorm(NB[t], 1/(sig.NB*sig.NB) ) # breeding males 
      
      counts.marked[1, t] ~ dpois(NFY[t]) # first year males, includes translocated/hacked
      counts.marked[2, t] ~ dpois(NF[t]) # nonbreeding adult males 
      counts.marked[3, t] ~ dpois(NB[t]) # breeding males 
      
      # log(NFY[t]) <- lNFY[t]
      # lNFY[t] ~ dnorm( 0, 1/(5*5) )
      # log(NF[t]) <- lNF[t]
      # lNF[t] ~ dnorm( 0, 1/(5*5) )
      # log(NB[t]) <- lNB[t]
      # lNB[t] ~ dnorm( 0, 1/(5*5) )
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
      eta.phiA[t] <- mean(phiA[1:nind,t])
      eta.phiB[t] <- mean(phiB[1:nind,t])
      eta.psiFYB[t] <- mean(psiFYB[1:nind,t])
      eta.psiAB[t] <- mean(psiAB[1:nind,t])
      eta.psiBA[t] <- mean(psiBA[1:nind,t])
      eta.pA[t] <- mean(pA[1:nind,t])
      eta.pB[t] <- mean(pB[1:nind,t])
    } # t
    
    # multivariate normal for temporal variance
  for (t in 1:nyr){
    eps[1:p,t] ~ dmnorm.vcov(mu.zeroes[1:p], sigma2[1:p, 1:p])
  }

# priors and diagonals
for (j in 1:p){ # coefficients
    sigma[j] ~ dexp(1)
    sigma2[j,j] <- sigma[j]*sigma[j]*rho[j,j]
    rho[j,j] <- 1
    #betas[j] ~ dnorm( 0, 1/(2*2) )
    lmus[j] <- logit(mus[j])
    mus[j] ~ dbeta(1,1)
  } # j

# upper diagonal  
  for (j in 1:(p-1)){
    for (k in (j+1):p){
      sigma2[j,k] <- sigma[j] * sigma[k] * rho[j,k] 
      sigma2[k,j] <- sigma2[j,k]
      rho[j,k] ~ dunif(-1,1)
      rho[k,j] <- rho[j,k]
    }}
  
    for (i in 1:nind){
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
} # model
",fill = TRUE)
sink() 

n.chains=1; n.thin=50; n.iter=100000; n.burnin=50000
n.chains=3; n.thin=5; n.iter=10000; n.burnin=5000; n.adapt=1000
n.chains=1; n.thin=1; n.iter=200; n.burnin=100

pars <- c(# pop growth 
  #"lambda",
  # fecundity
  "lmu.brood", "delta", "sig.brood", 
  "lmu.nest", "mu.nest", "gamma", 
  "sig.brood.t", "sig.nest",
  "mn.brood", "mn.nest", "mn.f", 
  #"lam", 
  # survival 
  "mus", "lmus", 
  # abundance
  "N", "NB", "NF", "NFY", "Ntot",
  # error terms
  "sigma2", "sigma", "rho",
  # yearly summaries
  'eta.phiFY', 'eta.phiA', 'eta.phiB',
  'eta.psiNB', 'eta.psiAB', 'eta.psiBA',
  'eta.pA', 'eta.pB',
  # goodness of fit
  # "dmape.obs", "dmape.rep", 
  # "tvm.obs", "tvm.rep",
  # "tturn.obs", "tturn.rep",
  "f.dmape.obs", "f.dmape.rep",
  "f.tvm.obs", "f.tvm.rep"
)

datl <- c(datl, constl)
# create initial values for y data
get.first <- function(x) min(which(x!=4), na.rm=T)
f <- apply(datl$y, 1, get.first)
get.last <- function(x) max(which(x!=4), na.rm=T)
l <- apply(datl$y, 1, get.last)

z.inits <- array(NA, dim=dim(datl$y), dimnames=dimnames(datl$y))
for (i in 1:constl$nind){
  if(l[i]==constl$nyr) { next } else{ 
    z.inits[i, (l[i]+1):constl$nyr] <- 4 }}

TFmat <- is.na(z.inits) & is.na(z)
for (i in 1:nrow(TFmat) ){  TFmat[i,1:f[i]] <- FALSE }

# we just want the 2s and 3s
suby <- array(NA, dim(datl$y))
suby[datl$y %in% c(2,3)] <- datl$y[datl$y %in% c(2,3)]
for (i in 1:nrow(z.inits) ){
  for (t in f[i]:ncol(z.inits) ){
    mx <- ifelse( max(suby[i,1:t], na.rm=T)=="-Inf", 2, max(suby[i,1:t], na.rm=T))
    if (TFmat[i,t]==TRUE & mx==2){
      z.inits[i,t] <- 2 }
    if (TFmat[i,t]==TRUE & mx==3){
      z.inits[i,t] <- mx }
  }}

# create inits for rhos
z.inits[ z %in% c(2,3) ] <- z[ z %in% c(2,3) ]
z.inits [datl$y %in% c(2,3)] <- NA

# Abundance
N <- array(NA, dim=c(7, constl$nyr) )
#N[] <- 0
N[1,] <- datl$counts.marked[1,]
# N[2,] <- round(datl$counts.marked[2,]/2)
# N[3,] <- round(datl$counts.marked[2,]/2)
N[4,] <- 0
# N[5,] <- round(datl$counts.marked[3,]/3)
# N[6,] <- round(datl$counts.marked[3,]/3)
# N[7,] <- round(datl$counts.marked[3,]/3)

inits.func <- function(){
  list(z = z.inits, 
       sigma = runif(p)#,
       #N = N)
)}

datl$z <- z
for (i in 1:datl$nind){
  datl$z[i,1:f[i]] <- NA
}

out <- jags( datl, 
            inits= inits.func,
            pars, modfl,  
            n.chains = n.chains, n.thin = n.thin, 
            n.burnin = n.burnin, n.adapt=n.adapt, n.iter=n.iter, 
            parallel=T, module=c("glm", "bugs"))

save(out,  
     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm-jags.Rdata")

#save(post,  
#     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm.Rdata")

#save(file=paste("./", m, ".Rdata", sep=""), list="out")