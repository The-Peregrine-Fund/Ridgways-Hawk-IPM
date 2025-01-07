## ---- pva --------
start <- Sys.time()
library('nimble')
library ('coda')

load("/bsuscratch/brianrolek/riha_ipm/data.rdata")
source("/bsuscratch/brianrolek/riha_ipm/MCMCvis.R")
load("/bsuscratch/brianrolek/riha_ipm/outputs/ipm_longrun.rdata")

# load("C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//outputs//ipm_longrun.rdata")
# load("data//data.rdata")
# library ("MCMCvis")

#**********************
#* Set data to reflect PVA scenarios
#**********************
constl$K <- 50 # number of future years
constl$SC <- 45 # number of PVA scenarios
datl$constraint_data <- rbind(datl$constraint_data, array(1, dim=c(constl$K,2)) ) # to help constrain FYs born + hacked to be positive
constl$effort2 <- rbind(constl$effort2, array(0, dim=c(constl$K,2), dimnames=list(2024:(2024+constl$K-1), c("LH", "PC"))))

num.treated <- rep( c(0, 15, 30, 45, 100, 0, 15, 30, 45, 100, 0, 15, 30, 45, 100), each=3 ) # 100s sub in for All but are over-ridden in model code
n.hacked <- rep( c(0, 0, 0, 0, 0, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10), each=3 )
surv.diff <- array(NA, dim=c(constl$SC, 3, 2), dimnames=list(scenario=1:constl$SC, stage=c('FY', 'NB', 'B'), site=c('LH', 'PC')))
surv.diff[,1,] <- matrix(c(rep( c(0, 0.130, 0.260), 15),
                           rep( c(0, 0.146, 0.292), 15)), ncol=2)
surv.diff[,2,] <- matrix(c(rep( c(0, 0.006, 0.012), 15),
                           rep( c(0, 0.022, 0.044), 15)), ncol=2)
surv.diff[,3,] <- matrix(c(rep( c(0, 0.016, 0.032), 15),
                           rep( c(0, 0.026, 0.052), 15)), ncol=2)
constl$surv.diff <- surv.diff 

hacked.counts <- array(0, dim=c(constl$SC, constl$nyr+constl$K, 2))
for (sc in 1:constl$SC){
  hacked.counts[sc,1:constl$nyr, 1:2] <- constl$hacked.counts[, 1:2] 
  hacked.counts[sc,14:(constl$nyr+constl$K), 1] <- -n.hacked[sc]
  hacked.counts[sc,14:(constl$nyr+constl$K), 2] <-  n.hacked[sc]
}
constl$hacked.counts <- hacked.counts

#**********************
#* Extract posterior of IPM
#**********************
out <- lapply(post, as.mcmc)
outp <- MCMCpstr(out[-5], type="chains") # omit chain 5, convergence problem

#**********************
#* Simulation code
#**********************
nsims <- 2000 # simulate for 1000 iter
npost <- ncol(outp$gamma) # size of IPM posterior iterations
p <- constl$p
nyr <- constl$nyr
K <- 50
nsite <- constl$nsite
nsc <- constl$SC
betas <- outp$betas
deltas <- outp$deltas
lmus <- outp$lmus
Ustar <- outp$Ustar[1:p,1:p, ]
sds  <- outp$sds # prior for temporal variation

# Create arrays needed to store nodes.
eta <- array(NA, dim=c(p, nsite, nyr+K, nsims))
z.score <- array( rnorm(p*nsite*(nyr+K)*nsims, 0, sd=1), # generate random z.score values
                       dim=c(p, nsite, nyr+K, nsims))  
N <- array(NA, dim=c(nsc , 7, nyr+K, nsite, nsims))
NFY <- NF <- NB <- NAD <- Ntot <- 
  array(NA, dim=c(nsc , nyr+K, nsite, nsims))
lmu.prod  <- outp$lmu.prod
gamma  <- outp$gamma[1, ]
mn.prod <- mn.phiFY <- mn.phiA <- mn.phiB <- 
  mn.phiFY1 <- mn.phiA1 <- mn.phiB1 <-  
  mn.psiFYB <- mn.psiAB <- mn.psiBA <- 
  array(NA, dim=c(nsc, nyr+K, nsite, nsims))
num.hacked <- perc.hacked <- perc.hacked.5yr <- array(NA, dim=c(nsc, nyr+K, nsite, nsims) )
numer.perc.treat <- denom.perc.treat <- denom.perc.treat2 <- 
  perc.treat <- mn.prod <- array(NA, dim=c(nsc, nyr+K, nsite, nsims))
lambda <- array(NA, dim=c(nsc, nyr+K-1, nsite, nsims))
extinct <- extinctAD <- extinctB <- array(NA, dim=c(nsc, K, nsite, nsims)) 

# Backfill arrays with IPM posteriors in years with data
z.score[1:p, 1:nsite, 1:nyr, 1:nsims] <- outp$z.score  
N[1, 1:7, 1:nyr, 1:nsite, 1:nsims] <- outp$N
NFY[1, 1:nyr, 1:nsite, 1:nsims] <- outp$NFY
NF[1, 1:nyr, 1:nsite, 1:nsims] <- outp$NF
NB[1, 1:nyr, 1:nsite, 1:nsims] <- outp$NB
NAD[1, 1:nyr, 1:nsite, 1:nsims] <- outp$NAD
Ntot[1, 1:nyr, 1:nsite, 1:nsims] <- outp$Ntot
mn.prod[1, 1:nyr, 1:nsite, ] <- outp$mn.prod
mn.phiFY[1, 1:nyr, 1:nsite, ] <- outp$mn.phiFY
mn.phiA[1, 1:nyr, 1:nsite, ] <- outp$mn.phiA
mn.phiB[1, 1:nyr, 1:nsite, ] <- outp$mn.phiB
mn.psiFYB[1, 1:nyr, 1:nsite, ] <- outp$mn.psiFYB
mn.psiAB[1, 1:nyr, 1:nsite, ] <- outp$mn.psiAB
mn.psiBA[1, 1:nyr, 1:nsite, ] <- outp$mn.psiBA

#*********************
#* Simulations
#*********************
# Site and temporal random effects
    for (t in 1:(nyr+K) ){ 
      for (s in 1:nsite){
        for (ii in 1:nsims){
        eta[1:p,s,t, ii] <- diag(sds[1:p, ii]) %*% t(Ustar[1:p,1:p, ii]) %*% z.score[1:p,s,t, ii]
      } } }  # i s t

# Abundance
# Assign IPM posteriors to all scenarios
for (sc in 2:nsc){
  N[sc, 1:7, 1:nyr, 1:nsite, ] <- N[1, 1:7, 1:nyr, 1:nsite, ]
  NFY[sc, 1:nyr, 1:nsite, ] <- NFY[1, 1:nyr, 1:nsite, ]
  NF[sc, 1:nyr, 1:nsite, ] <- NF[1, 1:nyr, 1:nsite, ]
  NB[sc, 1:nyr, 1:nsite, ] <- NB[1, 1:nyr, 1:nsite, ]
  NAD[sc, 1:nyr, 1:nsite, ] <- NAD[1, 1:nyr, 1:nsite, ]
  Ntot[sc, 1:nyr, 1:nsite, ] <- Ntot[1, 1:nyr, 1:nsite, ]
  num.hacked[sc, 1:nyr, 1:nsite, ] <- num.hacked[1, 1:nyr, 1:nsite, ]
  
  mn.prod[sc, 1:nyr, 1:nsite, ] <- mn.prod[1, 1:nyr, 1:nsite, ]
  mn.phiFY[sc, 1:nyr, 1:nsite, ] <- mn.phiFY[1, 1:nyr, 1:nsite, ]
  mn.phiA[sc, 1:nyr, 1:nsite, ] <- mn.phiA[1, 1:nyr, 1:nsite, ]
  mn.phiB[sc, 1:nyr, 1:nsite, ] <- mn.phiB[1, 1:nyr, 1:nsite, ]
  mn.psiFYB[sc, 1:nyr, 1:nsite, ] <- mn.psiFYB[1, 1:nyr, 1:nsite, ]
  mn.psiAB[sc, 1:nyr, 1:nsite, ] <- mn.psiAB[1, 1:nyr, 1:nsite, ]
  mn.psiBA[sc, 1:nyr, 1:nsite, ] <- mn.psiBA[1, 1:nyr, 1:nsite, ]
} # sc
 
# Calculate the number hacked, so it can decrease
# as populations at Los Haitises decrease
for (t in 1:nyr){
  for (s in 1:nsite){
    num.hacked[1, t, s, ] <- nimStep( N[1, 1, t, 1, ] + 
                                        hacked.counts[1, t, 1] ) * hacked.counts[1, t, s]  +
      (1- nimStep( N[1, 1, t, 1, ] + hacked.counts[1, t, 1] )) * N[1, 1, t, 1, ] * 
      (-1*nimEquals(s,1) + nimEquals(s,2)) # change to 
  }# s
} # t
for (sc in 2:nsc){
  num.hacked[sc, 1:nyr, 1:2, ] <- num.hacked[1, 1:nyr, 1:2, ]
}

# Future projections of survival, recruitment, and detection
  for (t in 7:nyr){
    for (s in 1:nsite){    
      # FYs we calculate the percent hacked for each year. 
      perc.hacked[1, t, s, ] <- ifelse(NFY[1, t, s, ]>0, 
                                        num.hacked[1, t, s, ] / NFY[1, t, s, ] *
                                          nimEquals(s,2), 0) # set LH to zero
      # For adults we average over the previous 5 years as an approximation
      perc.hacked.5yr[1, t, s, ] <- ifelse( sum( NFY[1, (t-6):(t-1), s, ])>0,
                                             ( sum( num.hacked[1, (t-6):(t-1), s, ] ) / 
                                                 sum( NFY[1, (t-6):(t-1), s, ] ) )*
                                               nimEquals(s,2), 0) # set LH to zero 
    }}
for (sc in 2:nsc){
  perc.hacked[sc, 1:nyr, 1:2, ] <- perc.hacked[1, 1:nyr, 1:2, ]
  perc.hacked.5yr[sc, 1:nyr, 1:2, ] <- perc.hacked.5yr[1, 1:nyr, 1:2, ]
}

cap <- 0.99 # cap on maximum survival
    # Future projections- population model 
    for (sc in 1:nsc){
      for (t in nyr:(nyr+K-1)){
        for (s in 1:nsite){

          # Calculate the percent hacked for each year for survival and recruitment 
          perc.hacked[sc, t, s, ] <- ifelse(NFY[sc, t, s, ]>0, 
                                            num.hacked[sc, t, s, ] / NFY[sc, t, s, ] *
                                            nimEquals(s,2), 0) # set LH to zero
          # For adults we average over the previous 5 years as an approximation
          perc.hacked.5yr[sc, t, s, ] <- ifelse( sum( NFY[sc, (t-6):(t-1), s, ])>0,
                                                 ( sum( num.hacked[sc, (t-6):(t-1), s, ] ) / 
                                                  sum( NFY[sc, (t-6):(t-1), s, ] ) )*
                                                  nimEquals(s,2), 0) # set LH to zero 
          
          # Survival
          mn.phiFY1[sc, t, s, ] <- ilogit( eta[1, s, t, ] + 
                                             lmus[1, s, ] + betas[1, ]*perc.hacked[sc, t, s, ] ) + # change perc.hacked to zero at LH bc no birds translocated there
            surv.diff[sc, 1, s]  
          mn.phiA1[sc, t, s, ] <- ilogit( eta[2, s, t, ] + 
                                            lmus[2, s, ] ) + 
            surv.diff[sc, 2, s]
          mn.phiB1[sc, t, s, ] <- ilogit( eta[3, s, t, ] +  
                                            lmus[3, s, ] + betas[3, ]*perc.hacked.5yr[sc, t, s, ] ) + # change perc.hacked to zero at LH bc no birds translocated there
            surv.diff[sc, 3, s]
          
          # enforce the survival cap
          mn.phiFY[sc, t, s, ] <- nimStep(cap - mn.phiFY1[sc, t, s, ]) * mn.phiFY1[sc, t, s, ] + 
            (1 - nimStep(cap - mn.phiFY1[sc, t, s, ])) * cap
          mn.phiA[sc, t, s, ] <- nimStep(cap - mn.phiA1[sc, t, s, ]) * mn.phiA1[sc, t, s, ] + 
            (1 - nimStep(cap - mn.phiA1[sc, t, s, ])) * cap
          mn.phiB[sc, t, s, ] <- nimStep(cap - mn.phiB1[sc, t, s, ]) * mn.phiB1[sc, t, s, ] + 
            (1 - nimStep(cap - mn.phiB1[sc, t, s, ])) * cap  
          
          #Recruitment
          mn.psiFYB[sc, t, s, ] <- ilogit( eta[4, s, t, ] +  
                                             lmus[4, s, ] )    # first year to breeder
          mn.psiAB[sc, t, s, ] <- ilogit( eta[5, s, t, ] +  # nonbreeder to breeder
                                            lmus[5, s, ] + betas[5, ]*perc.hacked.5yr[sc, t, s, ] ) # change perc.hacked to zero at LH bc no birds translocated there
          mn.psiBA[sc, t, s, ] <- ilogit( lmus[6, s, ]  ) # breeder to nonbreeder
          
          # Abundance of nonbreeders
          ## Second year nonbreeders
          N[sc, 2, t+1, s, ] <-  rbinom(n=nsims,
                                      p=mn.phiFY[sc, t, s, ]*(1-mn.psiFYB[sc, t, s, ]), 
                                      size=NFY[sc, t, s, ]) # Nestlings to second year nonbreeders
          ## Adult nonbreeders
          N[sc, 3, t+1, s, ] <-  rbinom(n=nsims,
                                      p=mn.phiA[sc, t, s, ]*(1-mn.psiAB[sc, t, s, ]), 
                                      size=NF[sc, t, s, ]) # Nonbreeders to nonbreeders
          N[sc, 4, t+1, s, ] <-  rbinom(n=nsims,
                                      p=mn.phiB[sc, t, s, ]*mn.psiBA[sc, t, s, ], 
                                      size=NB[sc, t, s, ]) # Breeders to nonbreeders
          # Abundance of breeders
          ## Second year breeders
          N[sc, 5, t+1, s, ] <-  rbinom(n=nsims,
                                      p=mn.phiFY[sc, t, s, ]*mn.psiFYB[sc, t, s, ], 
                                      size=NFY[sc, t, s, ]) # Nestlings to second year breeders
          ## Adult breeders
          N[sc, 6, t+1, s, ] <-  rbinom(n=nsims,
                                      p=mn.phiA[sc, t, s, ]*mn.psiAB[sc, t, s, ], 
                                      size=NF[sc, t, s, ]) # Nonbreeder to breeder
          N[sc, 7, t+1, s, ] <-  rbinom(n=nsims,
                                      p=mn.phiB[sc, t, s, ]*(1-mn.psiBA[sc, t, s, ]), 
                                      size=NB[sc, t, s, ]) # Breeder to breeder

          # Add the population segments into life stages
          NF[sc, t+1, s, ] <- sum(N[sc, 2:4, t+1, s, ]) # number of adult nonbreeder females
          NB[sc, t+1, s, ] <- sum(N[sc, 5:7, t+1, s, ]) # number of adult breeder females
          NAD[sc, t+1, s, ] <- NF[sc, t+1, s, ] + NB[sc, t+1, s, ] # number of adults
          
          # calc precent of territories treated for parasitic nest flies 
          # which affects productivity
          numer.perc.treat[sc, t+1, s, ] <-  # numerator
            nimStep( NB[sc, t+1, s, ]-(num.treated[sc]+1) ) * num.treated[sc]  + # num.treated > NB, use num.treated
            ( 1-nimStep( NB[sc, t+1, s, ]-(num.treated[sc]+1) )) * NB[sc, t+1, s, ]   #  num.treated < NB set to NB
          denom.perc.treat[sc, t+1, s, ] <- # denominator, total nests available for treatment
            nimStep( NB[sc, t+1, s, ]-(num.treated[sc]+1) ) * NB[sc, t+1, s, ]  + # num.treated > NB, use num.treated
            ( 1-nimStep( NB[sc, t+1, s, ]-(num.treated[sc]+1) )) * numer.perc.treat[sc, t+1, s, ]
          denom.perc.treat2[sc, t+1, s, ] <- ifelse(denom.perc.treat[sc, t+1, s, ]==0, 1, 0) +
            (1- ifelse(denom.perc.treat[sc, t+1, s, ]==0, 1, 0)) *denom.perc.treat[sc, t+1, s, ]
          perc.treat[sc, t+1, s, ] <- numer.perc.treat[sc, t+1, s, ] /
            denom.perc.treat2[sc, t+1, s, ] 
          # productivity
          mn.prod[sc, t+1, s, ] <- exp(lmu.prod[s, ] + 
                                       gamma*perc.treat[sc, t+1, s, ] + # treat.pair set to one here.
                                       eta[9, s, t+1, ] )
          
          # Number of wild born juvs
          N[sc, 1, t+1, s, ] <-  rpois(nsims, 
                                     lambda= (NFY[sc, t, s, ]*mn.phiFY[sc, t, s, ]*mn.psiFYB[sc, t, s, ] + # first year breeders
                                                NF[sc, t, s, ]*mn.phiA[sc, t, s, ]*mn.psiAB[sc, t, s, ] + # nonbreeders to breeders
                                                NB[sc, t, s, ]*mn.phiB[sc, t, s, ]*(1-mn.psiBA[sc, t, s, ])) # breeders remaining
                                     *(mn.prod[sc, t+1, s, ]/2) ) # end Poisson
          # Add the population segments into life stages
          # constrain N1+hacked.counts to be >=0, and allow it to shrink as the population shrinks
          # If LH has fewer birds than translocations, none are translocated
          num.hacked[sc, t+1, s, ] <- nimStep( N[sc, 1, t+1, 1, ] + hacked.counts[sc, t+1, 1] ) * hacked.counts[sc, t+1, s]  +
            (1- nimStep( N[sc, 1, t+1, 1, ] + hacked.counts[sc, t+1, 1] )) * N[sc, 1, t+1, 1, ] * 
            (-1*nimEquals(s,1) + nimEquals(s,2)) # change to negative value for LH when using
          NFY[sc, t+1, s, ] <- N[sc, 1, t+1, s, ] + num.hacked[sc, t+1, s, ]
          Ntot[sc, t+1, s, ] <- sum(N[sc, 1:7, t+1, s, ]) # total number of females 
                  }} # s t
    } # sc

    #######################
    # Derived params
    #######################
    for(sc in 1:nsc){
      for (s in 1:nsite){
        for (t in 1:(nyr+K-1)){
          lambda[sc, t, s, ] <-  Ntot[sc, t+1, s, ]/(Ntot[sc, t, s, ]) # population growth rate
        }} #t
      
      for (s in 1:nsite){
        for (t in 1:K){
          # quasi-extinction probabilities given population segments, N total, adults, and breeders
          extinct[sc, t, s, ] <- ifelse(Ntot[sc, nyr+t, s, ]==0, 1, 0) 
          extinctAD[sc, t, s, ] <- ifelse(NAD[sc, nyr+t, s, ]==0, 1, 0)
          extinctB[sc, t, s, ] <- ifelse(NB[sc, nyr+t, s, ]==0, 1, 0)
        }}
    } # sc

end <- Sys.time()

post <- list(runtime=end-start,
             NFY=NFY,
             NF=NF,
             NB=NB,
             NAD=NAD,
             Ntot=Ntot,
             N=N,
             lambda=lambda,
             extinct=extinct,
             extinctAD=extinctAD,
             extinctB=extinctB,
             mn.prod=mn.prod,
             mn.phiFY=mn.phiFY,
             mn.phiA=mn.phiA,
             mn.phiB=mn.phiB,
             mn.psiFYB=mn.psiFYB,
             mn.psiAB=mn.psiAB,
             mn.psiBA=mn.psiBA,
             num.hacked=num.hacked,
             perc.hacked=perc.hacked,
             perc.hacked.5yr=perc.hacked.5yr
)    

save(post,
     file="/bsuscratch/brianrolek/riha_ipm/outputs/pva-r.rdata")
# save(post, mycode, seeds, cpus,
#      file="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//outputs//pva_survival_longrun.rdata")

# postprocess
load('C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\pva-r.rdata')