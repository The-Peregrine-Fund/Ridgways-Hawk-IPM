# find initial values that work
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_da_longrun_2025_Apr_10.rdata")
library ('MCMCvis')
colnames(post[[1]][, 1:182])

wf <- list()
for(i in 1:length(post)){
wf[[i]]<- which((is.na(post[[i]][1,1:182]) | post[[i]][1,1:182]<0))
}

outp <- lapply(post, MCMCpstr, type="chains")
mnl <- which.min(lapply(wf, length ))
Ni <- outp[[mnl]]$N[,,,1]
# find which nodes have broken chains for each posterior
broken <- list()
for (j in 1:length(wf)){
  broken[[j]] <- names(wf[[mnl]]) %in% names(wf[[j]]) 
}


# replace N[5,2,1 and N[6,2,1]

# 5,2,1 s that worked
outp[[5]]$N[5,2,1,1]
outp[[7]]$N[5,2,1,1]
outp[[8]]$N[5,2,1,1]
outp[[10]]$N[5,2,1,1]
outp[[11]]$N[5,2,1,1]
outp[[13]]$N[5,2,1,1]
outp[[16]]$N[5,2,1,1]
outp[[18]]$N[5,2,1,1]
outp[[19]]$N[5,2,1,1]
outp[[20]]$N[5,2,1,1]

outp[[5]]$N[5,2,1,1]
outp[[7]]$N[5,2,1,1]
outp[[8]]$N[5,2,1,1]
outp[[10]]$N[5,2,1,1]
outp[[11]]$N[5,2,1,1]
outp[[13]]$N[5,2,1,1]
outp[[16]]$N[5,2,1,1]
outp[[18]]$N[5,2,1,1]
outp[[19]]$N[5,2,1,1]
outp[[20]]$N[5,2,1,1]

Ni[5,2,1] <- outp[[5]]$N[5,2,1,1]
Ni[6,2,1] <- outp[[5]]$N[6,2,1,1]

save(file="data//N_inits.rdata",
     Ni=Ni )
