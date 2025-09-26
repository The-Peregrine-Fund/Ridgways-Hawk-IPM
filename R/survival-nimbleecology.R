## ---- ipm --------
library('nimble')
library('parallel')
library ('coda')
library ('nimbleEcology')
load("data/data-dd.rdata")
y.first <- c()
for (i in 1:constl$nind){
  y.first[i] <- datl$y[i,constl$first[i]]  
}
constl$y.first <- y.first

cpus <- 3
constl$p <- 8
# library (MCMCvis)
# load("data/data-dd.rdata")
# load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm-find-inits.rdata")
# load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_shortrun.rdata")
# outp <- screen(post)
# rm(list=c("mycode", "post"))
# cpus <- 5
# info <- par_info[[1]]

#**********************
#* Model code
#**********************
mycode <- nimbleCode(
  {
    ###################################################
    # Priors and constraints
    ###################################################
    # survival, recruitment, and detection can be correlated
    for (k in 1:8){
      betas[k] ~ dnorm(0, sd=20)  # prior for translocations coefficients
      deltas[k] ~ dnorm(0, sd=10) # prior for survey effort coefficients
    } # k
    
    for (j in 1:8){
      for (s in 1:nsite){
        lmus[j,s] <- logit(mus[j,s])
        mus[j,s] ~ dbeta(1,1) # prior for overall means
      }}     # 
    
    # Temporal random effects and correlations between sites
    # Non-centered parameterization of the multivariate normal distribution to improve convergence
    for (jj in 1:p){ sds[jj] ~ dexp(1) }# prior for temporal variation
    R[1:p,1:p] <- t(Ustar[1:p,1:p]) %*% Ustar[1:p,1:p] # calculate rhos, correlation coefficients
    Ustar[1:p,1:p] ~ dlkj_corr_cholesky(eta=1.1, p=p) # Ustar is the Cholesky decomposition of the correlation matrix
    for (t in 1:nyr){ # survival params only have nyr-1, no problem to simulate from however
      for (s in 1:nsite){
        eta[1:p,s,t] <- diag(sds[1:p]) %*% t(Ustar[1:p,1:p]) %*% z.score[1:p,s,t]
        for(j in 1:p){
          z.score[j,s,t] ~ dnorm(0, sd=1)  # z-scores
        } # j
      } } # s t 
    

    for (i in 1:nind){
      for (t in 1:nyr){
        #Survival
        logit(phiFY[i,t]) <- eta[1, site[i,t],t] + 
          lmus[1, site[i,t]] + betas[1]*c(-1,1)[hacked[i]+1] #+  # first year
          #alphas[1, site[i,t]]*(Ntot[t, site[i,t]]-mnC[site[i,t]])
        logit(phiA[i,t]) <- eta[2, site[i,t],t] +  
          lmus[2, site[i,t]]  +  betas[2]*c(-1,1)[hacked[i]+1] #+ # nonbreeder
          #alphas[2, site[i,t]]*(Ntot[t, site[i,t]]-mnC[site[i,t]])
        logit(phiB[i,t]) <- eta[3, site[i,t],t] +  
          lmus[3, site[i,t]] + betas[3]*c(-1,1)[hacked[i]+1] #+ # breeder
          #alphas[3, site[i,t]]*(Ntot[t, site[i,t]]-mnC[site[i,t]])
        #Recruitment
        logit(psiFYB[i,t]) <- eta[4, site[i,t],t] +  
          lmus[4, site[i,t]] + betas[4]*c(-1,1)[hacked[i]+1] #+ # first year to breeder
          #alphas[4, site[i,t]]*(Ntot[t, site[i,t]]-mnC[site[i,t]])
        logit(psiAB[i,t]) <- eta[5, site[i,t],t] +  
          lmus[5, site[i,t]] + betas[5]*c(-1,1)[hacked[i]+1] #+ # nonbreeder to breeder
          #alphas[5, site[i,t]]*(Ntot[t, site[i,t]]-mnC[site[i,t]])
        logit(psiBA[i,t]) <-  
          lmus[6, site[i,t]]  # breeder to nonbreeder
        #Re-encounter
        logit(pA[i,t]) <- eta[7, site[i,t],t] + 
          lmus[7, site[i,t]] + betas[7]*c(-1,1)[hacked[i]+1] +
          deltas[5]*effort2[t, site[i,t]] + deltas[6]*effort2[t, site[i,t]]^2# resight of nonbreeders
        logit(pB[i,t]) <- eta[8, site[i,t],t] + 
          lmus[8, site[i,t]] + betas[8]*c(-1,1)[hacked[i]+1] + 
          deltas[7]*effort2[t, site[i,t]] + deltas[8]*effort2[t, site[i,t]]^2# resight of breeders
      }#t
    }#i
    
    # Define state-transition and observation matrices
    for (i in 1:nind){  
      # Define probabilities of state S(t+1) given S(t)
      for (t in first[i]:(nyr-1)){
        ps[i,1,1,t] <- 0
        ps[i,1,2,t] <- phiFY[i,t] * (1-psiFYB[i,t]) 
        ps[i,1,3,t] <- phiFY[i,t] * psiFYB[i,t]
        ps[i,1,4,t] <- (1-phiFY[i,t])
        
        ps[i,2,1,t] <- 0
        ps[i,2,2,t] <- phiA[i,t] * (1-psiAB[i,t])
        ps[i,2,3,t] <- phiA[i,t] * psiAB[i,t]
        ps[i,2,4,t] <- (1-phiA[i,t])
        
        ps[i,3,1,t] <- 0
        ps[i,3,2,t] <- phiB[i,t] * psiBA[i,t]
        ps[i,3,3,t] <- phiB[i,t] * (1-psiBA[i,t])
        ps[i,3,4,t] <- (1-phiB[i,t])
        
        ps[i,4,1,t] <- 0
        ps[i,4,2,t] <- 0
        ps[i,4,3,t] <- 0
        ps[i,4,4,t] <- 1
        
        # Define probabilities of O(t) given S(t)
        po[i,1,1,t] <- 1 
        po[i,1,2,t] <- 0
        po[i,1,3,t] <- 0
        po[i,1,4,t] <- 0
        
        po[i,2,1,t] <- 0
        po[i,2,2,t] <- pA[i,t]
        po[i,2,3,t] <- 0
        po[i,2,4,t] <- (1-pA[i,t])
        
        po[i,3,1,t] <- 0
        po[i,3,2,t] <- 0
        po[i,3,3,t] <- pB[i,t]
        po[i,3,4,t] <- (1-pB[i,t])
        
        po[i,4,1,t] <- 0
        po[i,4,2,t] <- 0
        po[i,4,3,t] <- 0
        po[i,4,4,t] <- 1
      } #t
    } #i
    
    # Likelihood 
    for (i in 1:nind){
      ps[i,1:4,1:4,nyr] <- matrix(0,4,4)
      y[i, (first[i] + 1):nyr] ~ 
        dDHMMo(init = ps[i, y.first[i], 1:4, first[i]], 
               probObs = po[i, 1:4, 1:4, first[i]:(nyr - 1)], 
               probTrans = ps[i, 1:4, 1:4, (first[i] + 1):(nyr - 1)], 
               len = nyr - first[i], 
               checkRowSums = 0)
    } #i
  } )

#**********************
#* Function to run model in NIMBLE
#**********************
run_ipm <- function(info, datl, constl, code, outp){
  library('nimble')
  library('nimbleEcology')
  library('coda')
  #library ('MCMCvis')
  #source("/bsuscratch/brianrolek/riha_ipm/MCMCvis.R")
  
  params <- c(
    # survival 
    "mus", "lmus", "betas", "deltas",
    # error terms
    "eta", "sds", "Ustar", "R", "z.score"
  )
  
  # initial values  
  inits.func <- function (){
    list(
      mus = array(runif(8*2), dim=c(8,2)),
      betas = rep(0, 8),
      deltas = rep(0, 8),
      sds = rep(0.01, 8),
      Ustar = diag(constl$p),
      z.score = array(runif(constl$p*constl$nsite*constl$nyr, -0.1, 0.1), 
                      dim=c(constl$p, constl$nsite, constl$nyr))
    )}
  inits <- inits.func()
  
  mod <- nimbleModel(code, 
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
  nc=1; nt=100; ni=100000; nb=50000
  #nc=1; nt=100; ni2=250000; nb=150000
  #nc=1; nt=1; ni=3; nb=1
  #nc=1; nt=1; ni=50000; nb=25000
  
  post <- runMCMC(chmc,
                  niter = ni,
                  nburnin = nb,
                  nchains = nc,
                  thin = nt,
                  samplesAsCodaMCMC = T,
                  setSeed = info$seed,
                  inits = inits)
  
  return(post)
} # run_ipm function end

# #*****************
# #* Run chains in parallel
# #*****************
this_cluster <- makeCluster(cpus)
post <- parLapply(cl = this_cluster,
                  X = par_info[1:cpus],
                  fun = run_ipm,
                  datl = datl,
                  constl = constl,
                  code = mycode,
                  outp = outp)
stopCluster(this_cluster)

# flnm2 <- paste0("/bsuscratch/brianrolek/riha_ipm/outputs/ipm-simp-",
#                 substr(Sys.time(), 1, 10), ".rdata")
# save(post, file = flnm2)
save(post, mycode,
  file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\survival-2025-09-18.rdata")

MCMCtrace(post, pdf=F, params="mus", Rhat=T )
MCMCtrace(post, pdf=F, params="betas", Rhat=T )
MCMCtrace(post, pdf=F, params="deltas", Rhat=T )
MCMCtrace(post, pdf=F, params="sds", Rhat=T )
MCMCtrace(post, pdf=F, params="R", Rhat=T )
