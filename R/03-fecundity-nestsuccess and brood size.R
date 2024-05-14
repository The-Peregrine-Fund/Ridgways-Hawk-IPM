library('nimble')
library('nimbleHMC')
library('MCMCvis')
library('parallel')
library ('coda')

datl <- list( # productivity data
  brood = lb$brood,
  nest.success = lp$nestsuccess
)

constl <- list( 
  nyr = nyr,
  treat.brood = lb$treat,
  nbrood = nrow(lb),
  year.brood = lb$year2, 
  yrind.brood = yrind.brood,
  brood.end = brood.end,
  mintrunc = min(lb$brood),
  maxtrunc = max(lb$brood),
  
  treat.nest = lp$treat,
  nnest = nrow(lp),
  year.nest = lp$year2,
  yrind.nest = yrind.nest,
  nest.end = nest.end
)

run_f <- function(seed, datl, constl){
  library('nimble')
  library('coda')
  library ('nimbleHMC')
  
  code <- nimbleCode(
    {
      # Priors for fecundity
      #lmu.brood ~ T(dnorm(0, sd=10), log(mintrunc), log(maxtrunc))
      #lmu.brood ~ T(dnorm(0, sd=10), mintrunc, maxtrunc)
      lmu.brood ~ dnorm(0, sd=10)
      
      delta ~ dunif(-5, 5)
      gamma ~ dunif(-5, 5)
      
      lmu.nest <- logit(mu.nest)
      mu.nest ~ dbeta(1,1) 
      
      sig.brood ~ dexp(1)
      sig.nest ~ dexp(1)
      sig ~ dexp(1)
      
      # Two models for fecundity, (1) brood size and (2) nest success  
      # Brood size
      for (k in 1:nbrood){
        brood[k] ~ dlnorm(lam[k], sdlog=sig) # truncating seems to break something
        #brood[k] ~ T(dnorm(lam[k], sd=sig), mintrunc, ) # truncated at the max observed fledglings per territory (3) 
        lam[k] <- lmu.brood +
                 delta*treat.brood[k] +
                 eps.brood[ year.brood[k] ]
        
        #brood[k] ~ dpois(lam[k]) # truncated at the max observed fledglings per territory (3) 
        # log(lam[k]) <- lmu.brood +
        #                 delta*treat.brood[k] +
        #                 eps.brood[ year.brood[k] ]
      } # k
      # Nest success       
      for (n in 1:nnest){
        nest.success[n] ~ dbern( nu[n] )
        logit(nu[n]) <- lmu.nest + gamma*treat.nest[n] + eps.nest[ year.nest[n] ]
      } # n 
      
      for ( t in 1:nyr){
        eps.brood[t] ~ dnorm(0, sd=sig.brood)
        eps.nest[t] ~ dnorm(0, sd=sig.nest)
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
        f.rep[k] ~ dlnorm(lam[k], sdlog=sig) # expected counts
        f.dssm.obs[k] <- abs( ( f.obs[k] - f.exp[k] ) / (f.obs[k]+0.001) )
        f.dssm.rep[k] <- abs( ( f.rep[k] - f.exp[k] ) / (f.rep[k]+0.001) )
      } # k
      f.dmape.obs <- sum(f.dssm.obs[1:nbrood])
      f.dmape.rep <- sum(f.dssm.rep[1:nbrood])
      f.tvm.obs <- sd(brood[1:nbrood])^2/mean(brood[1:nbrood])
      f.tvm.rep <- sd(f.rep[1:nbrood])^2/mean(f.rep[1:nbrood])
    }) # end code
  
  params <- c( "lmu.brood", "delta", "sig.brood", "sig",
               "lmu.nest", "mu.nest", "gamma", "sig.nest",
               "eps.brood", "eps.nest",
               "mn.brood", "mn.nest", "mn.f", 
               "lam", 
               "f.dmape.obs", "f.dmape.rep",
               "f.tvm.obs", "f.tvm.rep"
  )
  
  inits <- function(){ list(lmu.brood = runif(1, 1, 1.5),
                            delta = runif(1, -1, 1), 
                            sig.brood = rexp(1),
                            sig = rexp(1),
                            mu.nest = runif(1),
                            gamma = runif(1, -2, 2), 
                            sig.nest = rexp(1)
                            )}
  
  #n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
  n.chains=1; n.thin=10; n.iter=50000; n.burnin=25000
  #n.chains=1; n.thin=10; n.iter=200; n.burnin=100
  
  mod <- nimbleModel(code, calculate=T, 
                     constants = constl, 
                     data = datl, 
                     inits = inits(), 
                     buildDerivs = TRUE)
  
  mod$calculate()
  
  cmod <- compileNimble(mod)
  #confhmc <- configureMCMC(mod)
  confhmc <- configureHMC(mod)
  
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
post <- parLapply(cl = this_cluster, 
                  X = 1:4, 
                  fun = run_f, 
                  dat = datl, const = constl)
stopCluster(this_cluster)

save(post,  
     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\fecundity-lnormal.Rdata")


pars <- c( "lmu.brood", "delta", "sig.brood", "sig",
           "lmu.nest", "mu.nest", "gamma", "sig.nest"#,
           #"mn.brood", "mn.nest"
)

out <- MCMCpstr(post, type="chains")

MCMCsummary(post, digits=2, pars)
MCMCtrace(post, pars, pdf=F)

hist( exp(out$lam) )
range( exp(out$lam) )
hist( exp(out$lmu.brood) )

# plot number of fledglings from data versus 
# model estimates
modelest <- apply(out$mn.brood, 1, mean)
datest <- tapply(datl$brood, constl$year.brood, mean)
par(mfrow=c(1,1))
plot(2011:2023, datest, type="l", ylim=c(0,2))
lines(2011:2023, modelest, lty=2)

apply(out$mn.nest, 1, mean)
apply(out$mn.nest, 1, quantile, probs=c(0.025, 0.975))
apply(out$mn.f, 1, mean)

# compare data and model estimates directly 
# also compare with Woolaver et al. 2015
# 2005-2009
# Fledging rate of 0.64 fledglings per 
# active nest (fledgling nest -1 )
dataest <- tapply(lp$fledged, lp$year2, mean)
#modelest1 <- apply(out$mn.brood, 1, mean)
#modelest2 <- apply(out$mn.nest, 1, mean)
modelest3 <- apply(out$mn.f, 1, mean)
par(mfrow=c(1,1))
plot(2011:2023, dataest, type="l", 
     ylim=c(0.1,0.6), 
     xlab="Year", ylab="Fledglings per nest")
#lines(2011:2023, modelest1*modelest2, lty=2)
lines(2011:2023, modelest3, lty=3)


proptr <- tapply(constl$treat.nest, constl$year.nest, mean)
plot(2011:2023, proptr, type="l",
     ylab="Prop of nests treated", xlab="Year")
numnests <- tapply(datl$nest.success, constl$year.nest, length)
plot(2011:2023, numnests, type="l",
     ylab="Num of nests monitored", xlab="Year")

# check fit
plot.diag <- function(out, ratio=FALSE, lab=""){
  par(mfrow=c(1,1))
  # plot mean absolute percentage error
  samps <- MCMCpstr(out, "all", type="chains")
  mx <- max(c(samps$f.dmape.rep, samps$f.dmape.obs))
  mn <- min(c(samps$f.dmape.rep, samps$f.dmape.obs))
  plot(jitter(samps$f.dmape.obs), 
       jitter(samps$f.dmape.rep),
       main=paste0("Mean absolute percentage error\n",lab),
       ylab="Discrepancy replicate values",
       xlab="Discrepancy observed values", 
       xlim=c(mn,mx), ylim=c(mn,mx), 
       pch=16, cex=0.5, col="gray10")
  curve(1*x, from=mn, to=mx, add=T, lty=2, lwd=2, col="blue")
  bp1 <- round(mean(samps$f.dmape.rep > samps$f.dmape.obs),2)
  loc <- ifelse(bp1 < 0.5, "topleft", "bottomright")
  legend(loc, legend=bquote(p[B]==.(bp1)), bty="n", cex=2)
  
  if (ratio==TRUE){
    # plot variance/mean ratio
    xs <- c(samps$f.tvm.rep[1,], samps$f.tvm.obs[1,])
    hist(samps$f.tvm.rep, nclass=50,
         xlab="variance/mean ", main=NA, axes=FALSE,
         xlim=c(min(xs), max(xs)))
    abline(v=samps$f.tvm.obs[1,], col="red")
    axis(1); axis(2)
  }
  return(list('Bayesian p-value'=bp1))
}


plot.diag(post, ratio=T, lab="Poisson")
plot.diag(post, ratio=F, lab="log-Normal")








#############################
# JAGS
############################
library (jagsUI)
m<- c("fec-jags")
modfl <- paste("./", m, ".txt", sep="")
datl <- list( # productivity data
  brood = lb$brood,
  nest.success = lp$nestsuccess,
  nyr = nyr,
  treat.brood = lb$treat,
  nbrood = nrow(lb),
  year.brood = lb$year2, 
  yrind.brood = yrind.brood,
  brood.end = brood.end,
  mintrunc = min(lb$brood),
  maxtrunc = max(lb$brood),
  
  treat.nest = lp$treat,
  nnest = nrow(lp),
  year.nest = lp$year2,
  yrind.nest = yrind.nest,
  nest.end = nest.end
)

sink(modfl)
cat("
    model{
      lmu.brood ~ dnorm(0, 1/(10*10) )
      delta ~ dnorm(0, 1/(10*10))
      gamma ~ dnorm(0, 1/(10*10))
      logit(lmu.nest) <- mu.nest
      mu.nest ~ dbeta(1,1) 
      sig.brood ~ dnorm(0, 1/(2*2)) T(0,)
      sig.nest ~ dnorm(0, 1/(2*2)) T(0,)
      sig ~ dnorm(0, 1/(2*2)) T(0,)
      
      # Two models for fecundity, (1) brood size and (2) nest success  
      # Brood size
      for (k in 1:nbrood){
        brood[k] ~ dnorm(lam[k], 1/(sig*sig)) T( mintrunc, maxtrunc) # truncated at the max observed fledglings per territory (3) 
        #brood[k] ~ dpois(lam[k]) T( mintrunc, maxtrunc) # truncated at the max observed fledglings per territory (3) 
        lam[k] <- lmu.brood + 
                        delta*treat.brood[k] + 
                        eps.brood[ year.brood[k] ]
      } # k
      # Nest success       
      for (n in 1:nnest){
        nest.success[n] ~ dbern( nu[n] )
        logit(nu[n]) <- lmu.nest + gamma*treat.nest[n] + eps.nest[ year.nest[n] ]
      } # n 
      
      for ( t in 1:nyr){
        eps.brood[t] ~ dnorm(0, 1/(sig.brood*sig.brood) )
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
        mn.brood[t] <- mean(broodmat[t,1:brood.end[t]])
        mn.nest[t] <- mean(nestmat[t,1:nest.end[t]])
        mn.f[t] <- mn.brood[t]*mn.nest[t] # calc fecundity
      } # t
} # model
", fill = TRUE)
sink() 

params <- c( "lmu.brood", "delta", "sig.brood", "sig",
             "lmu.nest", "mu.nest", "gamma", "sig.nest",
             "eps.brood", "eps.nest",
             "mn.brood", "mn.nest", "mn.f", 
             "lam"
)

inits <- function(){ list(lmu.brood = runif(1, 0, 1.09),
                          delta = runif(1, -2, 2), 
                          sig.brood = runif(1),
                          eps.brood = runif(constl$nyr, -0.1, 0.1),
                          sig = runif(1),
                          mu.nest = runif(1),
                          gamma = runif(1, -2, 2), 
                          sig.nest = runif(1),
                          eps.nest = runif(constl$nyr, -0.1, 0.1)
)}

n.chains=1; n.thin=10; n.iter=20000; n.burnin=10000
out <- jags(datl, inits= inits, params, modfl,  
            n.chains = nc, n.thin = nt, n.burnin = nb, n.adapt=na, n.iter=ni, 
            parallel=T, module=c("glm", "bugs"))

save(out,  
     file="C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\surv.Rdata")

library (MCMCvis)
MCMCtrace(out, "mu", pdf=FALSE)
MCMCpstr