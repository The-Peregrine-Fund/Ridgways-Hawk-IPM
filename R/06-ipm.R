## ---- ipm1 --------
#################################
# The model
################################
library (jagsUI)
load("/scratch/brolek/aplo_ipm/data/final-data.Rdata")
datl$pp <- c(0, 0, tapply(datl$prod, datl$year.p, sum, na.rm=T), 42)
datl$pp[c(15,17)] <- NA
datl$countOM <- datl$countJM-datl$aug

m<- c("ipm")
modfl <- paste(".\\", m, ".txt", sep="")
sink(modfl)
cat("
    model{
    ####################################################
    ####################################################
    # Mark-resight-recovery data
    #   Observations (po) = y  
    #     1 seen first year (age 0)
    #     2 seen nonbreeder
    #     3 seen breeder
    #     4 seen dead
    #     5 not seen
    #   States (ps)
    #     1 alive first year
    #     2 alive nonbreeder
    #     3 alive breeder
    #     4 Recovered recently dead
    #     5 Dead not recovered/long dead
    #     6 Emigrated and alive
    #     7 Emigrated and dead
    #   Groups
    #     1 wild born
    #     2 hacked
    #   Sex
    #     1 female
    #     2 male
    #   Effort
    #     1 low
    #     2 high
    ###################################################
    # PARAMETERS
    #   sO: survival probability first year 
    #       (this is the letter 'O' rather than zero so jags can parse)
    #   sA: survival probability nonbreeders
    #   sB: survival probability breeders
    #   psiOB: recruitment probability from first-year to breeder
    #   psiAB: recruitment probability from nonbreeders to breeder
    #   psiBA: recruitment probability from breeder to nonbreeders 
    #   pA: resight probability nonbreeders
    #   pB: resight probability breeder
    #   F: fecundity
    #   omega: immigration
    ###################################################
    # Priors and constraints
    ###################################################
    rNB ~ dunif(0,50)   
    for (m in 1:2){ mu.F[m] ~ dunif(0,5) } # m # limits to help run model
    sigma.F ~ dnorm(0, 1/(2*2) )T(0,)
    r ~ dunif(0,1)
    
    # Survival loops for demographic categories by sex, hacked, effort 
    
    sigma.AS.s ~ dnorm(0, 1/(2*2) )T(0,)
    ASalpha <- logit(ASalpha1)
    ASalpha1 ~ dunif(0, 1)
    
    sigma.BS.s ~ dnorm(0, 1/(2*2) )T(0,)
    BSalpha <- logit(BSalpha1) 
    BSalpha1 ~ dunif(0, 1)
    
    sigma.BAR.psi ~ dnorm(0, 1/(2*2) )T(0,)
    BARalpha <- logit(BARalpha1)
    BARalpha1 ~ dunif(0, 1)
    
    sigma.pB ~ dnorm(0, 1/(2*2) )T(0,)
    
    for (k in 1:2){
    mu.pB[k]<- logit(mu.pB1[k])
    mu.pB1[k] ~ dunif(0, 1)
    } # k
    
    for (t in 1:(n.yr-1)){
    logit(eta.ASalpha[t]) <- ASalpha + eps.AS.phi[t]
    eps.AS.phi[t] ~ dnorm(0, 1/(sigma.AS.s*sigma.AS.s) )
    
    logit(eta.BSalpha[t]) <- BSalpha + eps.BS.phi[t]
    eps.BS.phi[t] ~ dnorm(0, 1/(sigma.BS.s*sigma.BS.s) )    
    
    logit(eta.BARalpha[t]) <- BARalpha + eps.BAR.psi[t]
    eps.BAR.psi[t] ~ dnorm(0, 1/(sigma.BAR.psi*sigma.BAR.psi))
    
    logit(eta.pB[t]) <- mu.pB[effort[t]] + eps.pB[t]
    eps.pB[t] ~ dnorm(0, 1/(sigma.pB*sigma.pB) )
    } #t
    
    for(s in 1:2){
    
    sigma.OBR.psi[s] ~ dnorm(0, 1/(2*2) )T(0,)
    OBRalpha[s] <- logit(OBRalpha1[s]) 
    OBRalpha1[s] ~ dunif(0, 1)
    
    for (t in 1:(n.yr-1)){
    logit(eta.OBRalpha[s,t]) <- OBRalpha[s] + eps.OBR.psi[s,t]
    eps.OBR.psi[s,t] ~ dnorm(0, 1/(sigma.OBR.psi[s]*sigma.OBR.psi[s]) )
    } #t
    
    for (h in 1:2){
    sigma.ABR.psi[s,h] ~ dnorm(0, 1/(2*2) )T(0,)
    ABRalpha[s,h] <- logit(ABRalpha1[s,h])
    ABRalpha1[s,h] ~ dunif(0, 1)
    
    sigma.OS.s[s,h] ~ dnorm(0, 1/(2*2) )T(0,)
    OSalpha[s,h] <- logit(OSalpha1[s,h])
    OSalpha1[s,h] ~ dunif(0, 1)
    
    for (t in 1:(n.yr-1)){
    logit(eta.ABRalpha[s,h,t]) <- ABRalpha[s,h] + eps.ABR.psi[s,h,t]
    eps.ABR.psi[s,h,t] ~ dnorm(0, 1/(sigma.ABR.psi[s,h]*sigma.ABR.psi[s,h]) )
    
    logit(eta.OSalpha[s,h,t]) <- OSalpha[s,h] + eps.OS.s[s,h,t]
    eps.OS.s[s,h,t] ~ dnorm(0, 1/(sigma.OS.s[s,h]*sigma.OS.s[s,h]) )
    } #t
    } } #h #s
    
    for (h in 1:2){
    sigma.pA[h] ~ dnorm(0, 1/(2*2) )T(0,)
    for (k in 1:2){
    mu.pA[k,h] <- logit(mu.pA1[k,h])
    mu.pA1[k,h] ~ dunif(0,1)
    } # k
    
    for (t in 1:(n.yr-1)){
    logit(eta.pA[h,t]) <- mu.pA[effort[t], h] + eps.pA[h,t]
    eps.pA[h,t] ~ dnorm(0, 1/(sigma.pA[h]*sigma.pA[h]) )
    }} #t #h
    
    #######################
    # Derived params
    #######################
    for (t in 3:(nyr-1)){    
    lambda[t] <-  Ntot[t+1]/(Ntot[t])
    loglambda[t] <- log(lambda[t])
    } #t
    
    ###############################
    # Likelihood for productivity
    ###############################
    for (t in 1:nyr){ 
    J[t] ~ dpois( f[t]*NB[t] )
    log(f[t]) <- log(mu.f[manage[t]]) + eps.f[t]
    eps.f[t] ~ dnorm (0, sd=sigma.f )
    }
    # GOF fecundity
    for (t in 1:nyr){ 
    J.rep[t] ~ dpois( f[t]*NB[t] )
    J.exp[t] <- f[t]*NB[t]
    d.obs[t] <- J[t]* log((J[t]+0.001)/(J.exp[t]+0.001)) - (J[t]-J.exp[t])
    d.rep[t] <- J.rep[t]*log((J.rep[t]+0.001)/(J.exp[t]+0.001)) - (J.rep[t]-J.exp[t])
    } # t
    dd.obs <- sum(d.obs)
    tvm.obs <- sd(J)^2/mean(J)
    dd.rep <- sum(d.rep)
    tvm.rep <- sd(J.rep)^2/mean(J.rep)
    
    ################################
    # Likelihood for counts
    ################################
    for (t in 1:(nyr-1)){
    # Number of wild born juvs
    N[1,t+1] ~ dpois(NFY[t]*phiFY[t]*fFY[1,t]/2 + 
                    NB[t]*phiFY[t]*fAD[2,t]/2)
    
    # Number of nonbreeders
    # nonbreeders to nonbreeders
    N[2,t+1] ~ dbin(phiNB[t]*(1-psiNB_B[t]), NA[t])
    # breeders to nonbreeders
    N[3,t+1] ~ dbin(phiB[t]*psiB_NB[t], NB[t])
    
    # Number of breeders
    # nonbreeders to breeders
    N[4,t+1] ~ dbin(phi_NB[t]*psi_NB_B[t], NA[t])
    # breeders to breeders
    N[5,t+1] ~ dbin(phi_B[t]*(1-psi_B_NB[t]), NB[t])
    } # t
    
    for (t in 1:nyr){
    Ntot[t] <- sum(N[c(1,2,3,4,5),t]) + aug[t] # total number
    NB[t] <- sum(N[c(4,5),t]) # number of breeders
    NA[t] <- sum(N[c(2,3),t])  # number of nonbreeders
    NFY[t] <- N[1,t] + aug[t] # number of total first years. Includes translocated first years
    NFYW[t] <- N[1,t] # Number of wild born first years
    NH[t] <- aug[t] # Number of hacked first years
    } # t
    
    # Observation process    
    for (t in 2:nyr){
    countB[t] ~ dnegbin(pNB[t],rNB) # breeding males negative binomial
    pNB[t] <- rNB/(rNB+NB[t])
    countA[t] ~ dpois(NA[t]) # nonbreeding adult males
    countFYW[t] ~ dpois(NW[t]) # first year wild males 
    } # t
    
    ###################
    # Assess GOF of the state-space models for counts
    # Step 1: Compute statistic for observed data
    # Step 2: Use discrepancy measure: mean absolute error
    # Step 3: Use test statistic: number of turns
    ###################
    for (t in 2:nyr){ 
    c.expB[t] <- NB[t] + 0.001 # expected counts adult breeder
    c.expA[t] <- NA[t] + 0.001 # nonbreeder
    c.expFYW[t] <- NFYW[t] + 0.001 # first year
    c.obsB[t] <- countB[t] + 0.001
    c.obsA[t] <- countA[t] + 0.001
    c.obsFYW[t] <- countFYW[t] + 0.001
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
    c.repB[t] ~ dnegbin(pNB[t],rNB) # expected counts
    c.repA[t] ~ dpois(NA[t] ) 
    c.repFYW[t] ~ dpois(NFYW[1,t] )
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
    for (i in 1:nind){
    for (t in 1:(nyr-1)){
    #Survival
    sO[i,t] <- eta.OSalpha[sex[i],hacked[i],t] # first year
    sA[i,t] <- eta.ASalpha[t] # nonbreeder
    sB[i,t] <- eta.BSalpha[t] # breeder
    #Recruitment
    psiOB[i,t] <- eta.OBRalpha[sex[i],t] # first year to breeder
    psiAB[i,t] <- eta.ABRalpha[sex[i], hacked[i],t] # nonbreederto breeder
    psiBA[i,t] <- eta.BARalpha[t] # breeder to nonbreeder
    #Re-encounter
    pA[i,t] <- eta.pA[hacked[i],t] # resight of nonbreeders
    pB[i,t] <- eta.pB[t]  # resight of breeders
    }#t
    }#i
    
    # Define state-transition and observation matrices
    for (i in 1:nind){  
    # Define probabilities of state S(t+1) given S(t)
    for (t in first[i]:(nyr-1)){
    ps[1,i,t,1] <- 0
    ps[1,i,t,2] <- sO[i,t] * (1-psiOB[i,t]) * (1-em[t])
    ps[1,i,t,3] <- sO[i,t] * psiOB[i,t] * (1-em[t])
    ps[1,i,t,4] <- (1-sO[i,t]) * r * (1-em[t])
    ps[1,i,t,5] <- (1-sO[i,t]) * (1-r) * (1-em[t])
    ps[1,i,t,6] <- sO[i,t] * em[t]
    ps[1,i,t,7] <- (1-sO[i,t]) * (1-r) * em[t]
    
    ps[2,i,t,1] <- 0
    ps[2,i,t,2] <- sA[i,t] * (1-psiAB[i,t]) * (1-em[t])
    ps[2,i,t,3] <- sA[i,t] * psiAB[i,t] * (1-em[t])
    ps[2,i,t,4] <- (1-sA[i,t]) * r * (1-em[t])
    ps[2,i,t,5] <- (1-sA[i,t]) * (1-r) * (1-em[t])
    ps[2,i,t,6] <- sA[i,t] * em[t]
    ps[2,i,t,7] <- (1-sA[i,t]) * (1-r) * em[t]
    
    ps[3,i,t,1] <- 0
    ps[3,i,t,2] <- sB[i,t] * psiBA[i,t] * (1-em[t])
    ps[3,i,t,3] <- sB[i,t] * (1-psiBA[i,t]) * (1-em[t])
    ps[3,i,t,4] <- (1-sB[i,t]) * r * (1-em[t])
    ps[3,i,t,5] <- (1-sB[i,t]) * (1-r) * (1-em[t])
    ps[3,i,t,6] <- sB[i,t] * em[t]
    ps[3,i,t,7] <- (1-sB[i,t]) * (1-r) * em[t]
    
    ps[4,i,t,1] <- 0
    ps[4,i,t,2] <- 0
    ps[4,i,t,3] <- 0
    ps[4,i,t,4] <- 0
    ps[4,i,t,5] <- 1
    ps[4,i,t,6] <- 0
    ps[4,i,t,7] <- 0
    
    ps[5,i,t,1] <- 0
    ps[5,i,t,2] <- 0
    ps[5,i,t,3] <- 0
    ps[5,i,t,4] <- 0
    ps[5,i,t,5] <- 1
    ps[5,i,t,6] <- 0
    ps[5,i,t,7] <- 0
    
    ps[6,i,t,1] <- 0
    ps[6,i,t,2] <- 0
    ps[6,i,t,3] <- 0
    ps[6,i,t,4] <- 0
    ps[6,i,t,5] <- 0
    ps[6,i,t,6] <- 1
    ps[6,i,t,7] <- 0
    
    ps[7,i,t,1] <- 0
    ps[7,i,t,2] <- 0
    ps[7,i,t,3] <- 0
    ps[7,i,t,4] <- 0
    ps[7,i,t,5] <- 0
    ps[7,i,t,6] <- 0
    ps[7,i,t,7] <- 1
    
    # Define probabilities of O(t) given S(t)
    po[1,i,t,1] <- 1 
    po[1,i,t,2] <- 0
    po[1,i,t,3] <- 0
    po[1,i,t,4] <- 0
    po[1,i,t,5] <- 0
    
    po[2,i,t,1] <- 0
    po[2,i,t,2] <- pA[i,t]
    po[2,i,t,3] <- 0
    po[2,i,t,4] <- 0
    po[2,i,t,5] <- 1-pA[i,t]
    
    po[3,i,t,1] <- 0
    po[3,i,t,2] <- 0
    po[3,i,t,3] <- pB[i,t]
    po[3,i,t,4] <- 0
    po[3,i,t,5] <- 1-pB[i,t]
    
    po[4,i,t,1] <- 0
    po[4,i,t,2] <- 0
    po[4,i,t,3] <- 0
    po[4,i,t,4] <- 1
    po[4,i,t,5] <- 0
    
    po[5,i,t,1] <- 0
    po[5,i,t,2] <- 0
    po[5,i,t,3] <- 0
    po[5,i,t,4] <- 0
    po[5,i,t,5] <- 1
    
    po[6,i,t,1] <- 0
    po[6,i,t,2] <- 0
    po[6,i,t,3] <- 0
    po[6,i,t,4] <- 0
    po[6,i,t,5] <- 1
    
    po[7,i,t,1] <- 0
    po[7,i,t,2] <- 0
    po[7,i,t,3] <- 0
    po[7,i,t,4] <- 0
    po[7,i,t,5] <- 1
    } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
    # Define latent state at first capture
    z[i,first[i]] <- y[i,first[i]]
    for (t in (first[i]+1):nyr){
    # State process: draw S(t) given S(t-1)
    z[i,t] ~ dcat(ps[z[i,t-1], i, t-1, 1:7])
    # Observation process: draw O(t) given S(t)
    y[i,t] ~ dcat(po[z[i,t], i, t-1, 1:5])
    } #t
    } #i
    } #model
    ",fill = TRUE)
sink() 

get.first <- function(x) min(which(x!=5))
f <- apply(datl$y, 1, get.first)
get.last <- function(x) max(which(x!=5))
l <- apply(datl$y, 1, get.last)
TFmat <- is.na(z.inits) & is.na(datl$z)
for (i in 1:dim(TFmat)[1]){  TFmat[i,1:f[i]] <- FALSE }
z.inits[TFmat] <- sample(size=445, c(2,3), replace=T, prob=c(0.5, 0.5) ) 

# Set initial values of N 
# to small numbers for immigrant and emigrants
# because preliminary analyses suggested low rates
# This is necessary to run model and avoid initial values error.
Ni <- array(NA, dim=c(17,26))
Ni[12:17, 2:26] <- 0

inits <- function(){list(z = z.inits, N=Ni)} 

params <- c( "rNB", 
             # fecundity
             "sigma.f", "eps.f", "mu.f",
             # survival 
             "OSalpha", "ASalpha", "BSalpha", "OBRalpha", "ABRalpha", "BARalpha",  "mu.pA", "mu.pB",
             "OSalpha1", "ASalpha1", "BSalpha1","OBRalpha1", "ABRalpha1", "BARalpha1","mu.pA1", "mu.pB1",
             "sigma.OS.s", "sigma.AS.s", "sigma.BS.s", "sigma.OBR.psi", "sigma.ABR.psi", "sigma.BAR.psi", "sigma.pA", "sigma.pB",
             # abundance
             "N", "NB", "NF", "NFY", "NFYW", "NH", "Ntot",
             # goodness of fit
             "dmape.obs", "dmape.rep", "tvm.obs", "tvm.rep", "dd.obs", "dd.rep",
             "tturn.obs", "tturn.rep",
             # error terms
             "eps.OS.s", "eps.AS.s", "eps.BS.s", "eps.OBR.psi", "eps.ABR.psi", "eps.BAR.psi", "eps.pA", "eps.pB", 
             "eta.OSalpha", "eta.ASalpha", "eta.BSalpha", "eta.OBRalpha", "eta.ABRalpha", "eta.BARalpha", "eta.pA",  "eta.pB", 
             "F"
)

# MCMC settings
ni <- 200000; nt <- 50; nb <- 150000; nc <- 3; na <- 1000 
out <- jags(datl, inits, params, modfl,  
            n.chains = nc, n.thin = nt, n.burnin = nb, n.adapt=na, n.iter=ni, 
            parallel=T, module=c("glm", "bugs"))
#save(file=paste("./", m, ".Rdata", sep=""), list="out")