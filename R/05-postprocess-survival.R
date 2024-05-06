library ('MCMCvis')
library('bayestestR')
library ('HDInterval')
library ('reshape2')
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\surv_nimble.Rdata")
post2 <- MCMCpstr(post, type="chains")

out <- MCMCpstr(post, type="chains")
MCMCvis::MCMCtrace(post, "R", pdf=FALSE)
MCMCvis::MCMCtrace(post, "mus", pdf=FALSE)
MCMCvis::MCMCtrace(post, "sds", pdf=FALSE)
MCMCvis::MCMCtrace(post, "betas", pdf=FALSE)
MCMCsummary(post, "R", HPD=TRUE, pg0=TRUE, round=3)
MCMCsummary(post, "mus", HPD=TRUE, pg0=TRUE, round=3)
MCMCsummary(post, "sds", HPD=TRUE, pg0=TRUE, round=3)
MCMCsummary(post, "betas", HPD=TRUE, pg0=TRUE, round=3)

eps <- apply(out$eps, c(1, 2), median) 
plot(eps[3,], eps[4,])


# Create a function to compare means
compare.means <- function(x, par) {
  # get the dimensions
  x2 <- x[par][[1]]
  ld <- melt(x2)
  ld$par.num <- paste0(ld$Var1, "_", ld$Var2)
  ld$value2 <- plogis(ld$value)
  npars <- length(unique(ld$par.num))
  niter <- length(unique(ld$Var3))
  
  nms <- sort(unique(ld$par.num))
  diff.ar <- diff.ar2 <- array(NA, dim=c(npars, npars, niter), 
                               dimnames=list(nms, nms))
  prd <- prd2 <- array(NA, dim=c( npars, npars ), dimnames=list(nms, nms))
  for (i in 1:(npars-1)){
    for (j in (i+1):npars){ # upper diagonal
      par1 <- nms[i]
      par2 <- nms[j]  
      p1 <- ld[ld$par.num == par1, "value"]
      p2 <- ld[ld$par.num == par2, "value"]
      p11 <- ld[ld$par.num == par1, "value2"]
      p22 <- ld[ld$par.num == par2, "value2"]
      diff.ar[i,j,] <- p1 - p2
      prd[i,j] <- round(pd(p1 - p2),3)
      diff.ar2[i,j,] <- p11 - p22
      prd2[i,j] <- round(pd(p11 - p22),3)
    }}
  dmn <- round(apply(diff.ar, c(1,2), mean),3)
  dlhd <- round(apply(diff.ar, c(1,2), HDInterval::hdi),3)[1,,]
  duhd<- round(apply(diff.ar, c(1,2), HDInterval::hdi),3)[2,,]
  
  df <- data.frame(Var1=melt(dmn)$Var1, 
                   Var2=melt(dmn)$Var2,
                   Mean=melt(dmn)$value, 
                   LHDI=melt(dlhd)$value, UHDI=melt(duhd)$value, 
                   PD=melt(prd)$value)
  df <- df[!is.na(df$Mean),]
  df <- df[order(df$Var1, df$Var2),]

return(df)
}

df <- compare.means(post2, "lmus")

#*****************************
#* Compare prior and posterior densities
#*****************************
prior <- rbeta(100000, 1, 1)
priord <- density(prior)
postd <- density(post2$mus[1,4,])
plot(postd, xlim=c(0,1))
lines(priord, lty=2)
