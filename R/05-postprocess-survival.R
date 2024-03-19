library ('MCMCvis')
library('bayestestR')
# Create a function to compare means
compare.means <- function(x) {
  # get the dimensions
  npars <- nrow(x)
  niter <- ncols(x)
  
  diff.ar <- array(NA, dim=c(npars, npars, niter))
  prd <- array(NA, dim=c( npars, npars ))
  for (i in 1:(npars-1)){
    for (j in (i+1):npars){ # upper diagonal
      diff.ar[i,j,] <- x$sims.list$a1[,i] - x$sims.list$a1[,j]
      prd[i,j] <- round(pd(x$sims.list$a1[,i] - x$sims.list$a1[,j])[[1]],3)
    }}
  dmn <- round(apply(diff.ar, c(1,2), mean),3)
  dlhd <- round(apply(diff.ar, c(1,2), HDInterval::hdi),3)[1,,]
  duhd<- round(apply(diff.ar, c(1,2), HDInterval::hdi),3)[2,,]
  
  df <- cbind(melt(dmn), melt(dlhd), melt(duhd),melt(prd))
  df <- df[!is.na(df[,3]), c(1,2,3,6,9,12)]
  colnames(df) <- c("v1", "v2", 
                    "dmn", "dlhd", "duhd",
                    "pd")
  keep <- ifelse( (df$v1 %in% 1:7) & (df$v2 %in% 1:7) |
                    (df$v1 %in% 8:10) & (df$v2 %in% 8:10) |
                    (df$v1 %in% 11:12) & (df$v2 %in% 11:12), 
                  TRUE, FALSE)
  df <- df[keep,]
  df <- df[order(df$v1),]
}

df <- compare.means(x)
