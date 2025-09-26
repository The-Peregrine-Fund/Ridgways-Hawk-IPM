## ---- postprocess --------
library ('MCMCvis')
library ('coda')
library ('ggplot2')
library('reshape2')
library ('tidybayes')
library ('bayestestR')
library ('ggpubr')
library('viridis')
library ('HDInterval')
library ('abind')
library ('ggimage')
load("data/data.rdata")
#load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_longrun.rdata")
#load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_longrun_2025_Apr_01.rdata")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_dd_longrun_2025_Apr_03.rdata")
out <- lapply(post[1:5], as.mcmc) # omit chain 5, causing lack of convergence in 1 param
NAlist <- c()
for (i in 1:length(out)){
  NAlist[i] <- any (is.na(out[[i]][,1:286]) | out[[i]][,1:286]<0)
}
# Subset chains to those with good initial values
out <- out[!NAlist]
outp <- MCMCpstr(out, type="chains")
niter <- dim(outp$N)[4]


