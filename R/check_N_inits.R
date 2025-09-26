# library(MCMCvis)
# load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_dd_longrun_2025_Apr_03.rdata")
# load("data/data-dd.rdata")

get_inits <- function(){
source("/bsuscratch/brianrolek/riha_ipm/MCMCvis.R")
load("/bsuscratch/brianrolek/riha_ipm/data-dd.rdata")
load("/bsuscratch/brianrolek/riha_ipm/outputs/ipm_dd_long_2025_Apr_03.rdata")

outp <- MCMCpstr(post, type="chains")
ind <- sample(1:ncol(outp$gamma), size=1, replace=TRUE)
N.ar <- array(NA, dim=c(9, 13, 2))
# Rearrange for updated model that
# allows immigration
N.ar[1,,1] <- outp$N[1, 1:13, 1, ind] |> rpois(13)
N.ar[2,,1] <- outp$N[2, 1:13, 1, ind] |> rpois(13) 
N.ar[3,,1] <- outp$N[3, 1:13, 1, ind] |> rpois(13)
N.ar[4,,1] <- outp$N[4, 1:13, 1, ind] |> rpois(13)
N.ar[6,,1] <- outp$N[5, 1:13, 1, ind] |> rpois(13) 
N.ar[7,,1] <- outp$N[6, 1:13, 1, ind] |> rpois(13)
N.ar[8,,1] <- outp$N[7, 1:13, 1, ind] |> rpois(13)

N.ar[1,,2] <- outp$N[1, 1:13, 2, ind] |> rpois(13)
N.ar[2,,2] <- outp$N[2, 1:13, 2, ind] |> rpois(13)
N.ar[3,,2] <- outp$N[3, 1:13, 2, ind] |> rpois(13)
N.ar[4,,2] <- outp$N[4, 1:13, 2, ind] |> rpois(13)
N.ar[6,,2] <- outp$N[5, 1:13, 2, ind] |> rpois(13)
N.ar[7,,2] <- outp$N[6, 1:13, 2, ind] |> rpois(13)
N.ar[8,,2] <- outp$N[7, 1:13, 2, ind] |> rpois(13)

N.ar[5,,1:2] <- 0
N.ar[9,,1:2] <- 0
# add them all up
nyr <- constl$nyr
nsite <- constl$nsite
hc <- constl$hacked.counts[,-3]
# SETUP COMPLETE

# A function to sum totals 
get_sums <- function(N.ar){
  NFY <- NF <- NB <- NAD <- Ntot <- array(NA, dim=c(nyr, nsite))
  for (t in 1:nyr){
  for (s in 1:nsite){
    NFY[t, s] <- N.ar[1, t, s] + hc[t, s] # Transfers translocated first-year females
    NF[t, s] <- sum(N.ar[2:5, t, s]) # number of adult nonbreeder females
    NB[t, s] <- sum(N.ar[6:9, t, s]) # number of adult breeder females
    NAD[t, s] <- NF[t, s] + NB[t, s]  # number of adults
    Ntot[t, s] <- sum(N.ar[1:9, t, s]) # total number of females
  }} # t,s
  return(list(NFY=NFY, NF=NF, NB=NB, NAD=NAD, Ntot=Ntot))
} # function

# START THE TESTS AND REASSIGNMENTS
# check N[1] is > num hacked
pos.N1 <- (N.ar[1, , ] + hc) >= 0  # constrain N1+hacked.counts to be >=0
while(any(pos.N1==FALSE)){
  # (1) First substitute FALSE values 
  # so N1 is greater than num.hacked.
  for (t in 1:(nyr-1)){
    for (s in 1:nsite){
  N.ar[1,t,s] <-  ifelse(pos.N1[t,s]==FALSE, abs(constl$hacked.counts[t,s])+10, N.ar[1,t,s]+10)
    }}
  gsums <- get_sums(N.ar)
  list2env(gsums, envir = .GlobalEnv)
  # retest
  pos.N1 <- (N.ar[1, , ] + hc) >= 0 
} # while

# (2) Next check that the dynamics add up
# N2 and N6 at t+1 cannot exceed NFY at t
# Test N.ar <= NFY
N2 <-  NFY[1:12,1:2]>=N.ar[2,2:13,1:2]
N6 <- NFY[1:12,1:2]>=N.ar[6,2:13,1:2]
while(any(N2==FALSE) | any(N6==FALSE)){
for (t in 1:(nyr-1)){
  for (s in 1:nsite){
    N.ar[2,t+1,s] <- ifelse(N2[t,s]==FALSE, rpois(1, NFY[t,s]/2),  N.ar[2,t+1,s])
    N.ar[6,t+1,s] <- ifelse(N6[t,s]==FALSE, rpois(1, NFY[t,s]/2),  N.ar[6,t+1,s])
}} # t s 
  gsums <- get_sums(N.ar)
  list2env(gsums, envir = .GlobalEnv)
  N2 <-  NFY[1:12,1:2]>=N.ar[2,2:13,1:2]
  N6 <- NFY[1:12,1:2]>=N.ar[6,2:13,1:2]
}

# Check nonbreeders
N3 <-  NF[1:12,1:2]>=N.ar[3,2:13,1:2]
N7 <- NF[1:12,1:2]>=N.ar[7,2:13,1:2]
while(any(N3==FALSE) | any(N7==FALSE)){
  for (t in 1:(nyr-1)){
    for (s in 1:nsite){
      N.ar[3,t+1,s] <- ifelse(N3[t,s]==FALSE, rpois(1, NF[t,s]/2),  N.ar[3,t+1,s])
      N.ar[7,t+1,s] <- ifelse(N7[t,s]==FALSE, rpois(1, NF[t,s]/2),  N.ar[7,t+1,s])
    }} # t s 
  gsums <- get_sums(N.ar)
  list2env(gsums, envir = .GlobalEnv)
  N3 <-  NF[1:12,1:2]>=N.ar[3,2:13,1:2]
  N7 <- NF[1:12,1:2]>=N.ar[7,2:13,1:2]
}

# check breeders
N4 <-  NB[1:12,1:2]>=N.ar[4,2:13,1:2]
N8 <- NB[1:12,1:2]>=N.ar[8,2:13,1:2]
while(any(N4==FALSE) | any(N8==FALSE)){
  for (t in 1:(nyr-1)){
    for (s in 1:nsite){
      N.ar[4,t+1,s] <- ifelse(N4[t,s]==FALSE, rpois(1, NB[t,s]/2),  N.ar[4,t+1,s])
      N.ar[8,t+1,s] <- ifelse(N8[t,s]==FALSE, rpois(1, NB[t,s]/2),  N.ar[8,t+1,s])
    }} # t s 
  gsums <- get_sums(N.ar)
  list2env(gsums, envir = .GlobalEnv)
  N4 <-  NB[1:12,1:2]>=N.ar[4,2:13,1:2]
  N8 <- NB[1:12,1:2]>=N.ar[8,2:13,1:2]
}

pos.N1 <- (N.ar[1, , ] + hc) >= 0
N2 <-  NFY[1:12,1:2]>=N.ar[2,2:13,1:2]
N6 <- NFY[1:12,1:2]>=N.ar[6,2:13,1:2]
N3 <-  NF[1:12,1:2]>=N.ar[3,2:13,1:2]
N7 <- NF[1:12,1:2]>=N.ar[7,2:13,1:2]
N4 <-  NB[1:12,1:2]>=N.ar[4,2:13,1:2]
N8 <- NB[1:12,1:2]>=N.ar[8,2:13,1:2]

return(N.ar)
} # end function
#save(file = "data/N-inits.rdata", N.ar)



# Create a separate function
# For model including probability 
# of nest treatment
get_inits2 <- function(){
  source("/bsuscratch/brianrolek/riha_ipm/MCMCvis.R")
  load("/bsuscratch/brianrolek/riha_ipm/data-dd.rdata")
  load("/bsuscratch/brianrolek/riha_ipm/outputs/ipm_dd_long_2025_Apr_03.rdata")
  
  outp <- MCMCpstr(post, type="chains")
  ind <- sample(1:ncol(outp$gamma), size=1, replace=TRUE)
  N.ar <- array(NA, dim=c(10, 13, 2))
  # Rearrange for updated model that
  # allows immigration
  N.ar[1,,1] <- (outp$N[1, 1:13, 1, ind]*1) |> rpois(13)
  N.ar[2,,1] <- 0
    #(outp$N[1, 1:13, 1, ind]*0.2) |> rpois(13)
  
  N.ar[3,,1] <- outp$N[2, 1:13, 1, ind] |> rpois(13) 
  N.ar[4,,1] <- outp$N[3, 1:13, 1, ind] |> rpois(13)
  N.ar[5,,1] <- outp$N[4, 1:13, 1, ind] |> rpois(13)
  N.ar[7,,1] <- outp$N[5, 1:13, 1, ind] |> rpois(13) 
  N.ar[8,,1] <- outp$N[6, 1:13, 1, ind] |> rpois(13)
  N.ar[9,,1] <- outp$N[7, 1:13, 1, ind] |> rpois(13)
  
  N.ar[1,,2] <- (outp$N[1, 1:13, 2, ind]*0.8) |> rpois(13)
  N.ar[2,,2] <- (outp$N[1, 1:13, 2, ind]*0.2) |> rpois(13)
  
  N.ar[3,,2] <- outp$N[2, 1:13, 2, ind] |> rpois(13)
  N.ar[4,,2] <- outp$N[3, 1:13, 2, ind] |> rpois(13)
  N.ar[5,,2] <- outp$N[4, 1:13, 2, ind] |> rpois(13)
  N.ar[7,,2] <- outp$N[5, 1:13, 2, ind] |> rpois(13)
  N.ar[8,,2] <- outp$N[6, 1:13, 2, ind] |> rpois(13)
  N.ar[9,,2] <- outp$N[7, 1:13, 2, ind] |> rpois(13)
  
  N.ar[6,,1:2] <- 0
  N.ar[10,,1:2] <- 0
  # add them all up
  nyr <- constl$nyr
  nsite <- constl$nsite
  hc <- constl$hacked.counts[,-3]
  # SETUP COMPLETE
  
  # A function to sum totals 
  get_sums <- function(N.ar){
    NFY <- NF <- NB <- NAD <- Ntot <- array(NA, dim=c(nyr, nsite))
    for (t in 1:nyr){
      for (s in 1:nsite){
        NFY[t, s] <- N.ar[1, t, s] + N.ar[2, t, s] + hc[t, s] # Transfers translocated first-year females
        NF[t, s] <- sum(N.ar[3:6, t, s]) # number of adult nonbreeder females
        NB[t, s] <- sum(N.ar[7:10, t, s]) # number of adult breeder females
        NAD[t, s] <- NF[t, s] + NB[t, s]  # number of adults
        Ntot[t, s] <- sum(N.ar[1:10, t, s]) # total number of females
      }} # t,s
    return(list(NFY=NFY, NF=NF, NB=NB, NAD=NAD, Ntot=Ntot))
  } # function
  
  # START THE TESTS AND REASSIGNMENTS
  # check N[1] is > num hacked
  pos.N1 <- (N.ar[1, , ] + N.ar[2, , ] + hc) >= 0  # constrain N1+hacked.counts to be >=0
  while(any(pos.N1==FALSE)){
    # (1) First substitute FALSE values 
    # so N1 is greater than num.hacked.
    for (t in 1:(nyr-1)){
      for (s in 1:nsite){
        N.ar[1,t,s] <-  ifelse(pos.N1[t,s]==FALSE, abs(constl$hacked.counts[t,s])+10, N.ar[1,t,s]+10)
        #N.ar[2,t,s] <-  ifelse(pos.N1[t,s]==FALSE, abs(constl$hacked.counts[t,s])+2, N.ar[1,t,s]+2)
      }}
    gsums <- get_sums(N.ar)
    list2env(gsums, envir = .GlobalEnv)
    # retest
    pos.N1 <- (N.ar[1, , ] + N.ar[2, , ] + hc) >= 0 
  } # while
  gsums <- get_sums(N.ar)
  list2env(gsums, envir = .GlobalEnv)
  
  # (2) Next check that the dynamics add up
  # N2 and N6 at t+1 cannot exceed NFY at t
  # Test N.ar <= NFY
  N3 <-  NFY[1:12,1:2]>=N.ar[3,2:13,1:2]
  N7 <- NFY[1:12,1:2]>=N.ar[7,2:13,1:2]
  while(any(N3==FALSE) | any(N7==FALSE)){
    for (t in 1:(nyr-1)){
      for (s in 1:nsite){
        N.ar[3,t+1,s] <- ifelse(N3[t,s]==FALSE, rpois(1, NFY[t,s]/2),  N.ar[3,t+1,s])
        N.ar[7,t+1,s] <- ifelse(N7[t,s]==FALSE, rpois(1, NFY[t,s]/2),  N.ar[7,t+1,s])
      }} # t s 
    gsums <- get_sums(N.ar)
    list2env(gsums, envir = .GlobalEnv)
    N3 <-  NFY[1:12,1:2]>=N.ar[3,2:13,1:2]
    N7 <- NFY[1:12,1:2]>=N.ar[7,2:13,1:2]
  }
  gsums <- get_sums(N.ar)
  list2env(gsums, envir = .GlobalEnv)
  
  # Check nonbreeders
  N4 <-  NF[1:12,1:2]>=N.ar[4,2:13,1:2]
  N8 <- NF[1:12,1:2]>=N.ar[8,2:13,1:2]
  while(any(N4==FALSE) | any(N8==FALSE)){
    for (t in 1:(nyr-1)){
      for (s in 1:nsite){
        N.ar[4,t+1,s] <- ifelse(N4[t,s]==FALSE, rpois(1, NF[t,s]/2),  N.ar[4,t+1,s])
        N.ar[8,t+1,s] <- ifelse(N8[t,s]==FALSE, rpois(1, NF[t,s]/2),  N.ar[8,t+1,s])
      }} # t s 
    gsums <- get_sums(N.ar)
    list2env(gsums, envir = .GlobalEnv)
    N4 <-  NF[1:12,1:2]>=N.ar[4,2:13,1:2]
    N8 <- NF[1:12,1:2]>=N.ar[8,2:13,1:2]
  }
  gsums <- get_sums(N.ar)
  list2env(gsums, envir = .GlobalEnv)
  
  # check breeders
  N5 <-  NB[1:12,1:2]>=N.ar[5,2:13,1:2]
  N9 <- NB[1:12,1:2]>=N.ar[9,2:13,1:2]
  while(any(N5==FALSE) | any(N9==FALSE)){
    for (t in 1:(nyr-1)){
      for (s in 1:nsite){
        N.ar[5,t+1,s] <- ifelse(N5[t,s]==FALSE, rpois(1, NB[t,s]/2),  N.ar[5,t+1,s])
        N.ar[9,t+1,s] <- ifelse(N9[t,s]==FALSE, rpois(1, NB[t,s]/2),  N.ar[9,t+1,s])
      }} # t s 
    gsums <- get_sums(N.ar)
    list2env(gsums, envir = .GlobalEnv)
    N5 <-  NB[1:12,1:2]>=N.ar[5,2:13,1:2]
    N9 <- NB[1:12,1:2]>=N.ar[9,2:13,1:2]
  }
  gsums <- get_sums(N.ar)
  list2env(gsums, envir = .GlobalEnv)
  
  pos.N1 <- (N.ar[1, , ] + N.ar[2, , ] + hc) >= 0
  N3 <-  NFY[1:12,1:2]>=N.ar[3,2:13,1:2]
  N7 <- NFY[1:12,1:2]>=N.ar[7,2:13,1:2]
  N4 <-  NF[1:12,1:2]>=N.ar[4,2:13,1:2]
  N8 <- NF[1:12,1:2]>=N.ar[8,2:13,1:2]
  N5 <-  NB[1:12,1:2]>=N.ar[5,2:13,1:2]
  N9 <- NB[1:12,1:2]>=N.ar[9,2:13,1:2]
  N.ar[2,,] <- 0
  return(N.ar)
} # end function
#save(file = "data/N-inits.rdata", N.ar)
