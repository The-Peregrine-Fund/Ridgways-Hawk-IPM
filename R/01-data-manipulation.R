#************************************
#* Data manipulation script for 
#* Ridgway's Hawk IPM in DR and Haiti
#* 
#* Author: Brian Rolek
#* Email: rolek.brian@peregrinefund.org
#* Organization: The Peregrine Fund
#************************************
# hacking = # hacked
# success= successful hack

# Link history from banding data 
# to territory name, number and year combined in nest success
# los anianas vargas (2018), los brazos (hack site 2019), 
# los cameita (2023), and punta cana = all hacked to start

#**********************
#* 1. Load packages and data ----
#**********************
library('lubridate')
library('readxl')
library ('dplyr')
library('ggplot2')
library ('reshape2')

dfl <- "data\\RIHA data 2011-2023.xlsx"
df <- read_excel(dfl, sheet='Individuals', na=c('#N/A', 'n/a','N/A', '', 'NA', '-') ) 
df$year.resighted <- as.numeric(df$year.resighted)
df$origin[ df$origin=="lHNP" ] <- "LHNP"

# Reassign bands to abbreviate
df$right.color[df$right.color=="Not recorded"] <- "NR"
df$right.color[df$right.color=="No band"] <- "NB"
df$left.color[df$left.color=="Not recorded"] <- "NR"
df$left.color[df$left.color=="No band"] <- "NB"
df$r.code[df$r.code=="Not recorded"] <- "NR"
df$r.code[df$r.code=="None"] <- "NB"
df$l.code[df$l.code=="Not recorded"] <- "NR"
df$l.code[df$l.code=="None"] <- "NB"
df$ID <- paste(df$l.code, df$r.code, 
               df$`left.color`, df$`right.color`, sep="-")
# subset
df <- df[df$sex=="Female",]
df <- df[df$current.population!="AV",] # remove because small sample size

# separate marked and unmarked birds
not.banded <- c("NB-NB-NB-NB", "NB-NR-NB-NR", "NR-NR-NR-NR")
df.unmarked <- df[df$ID %in% not.banded,] # unmarked df
df.unmarked.nodups <- df.unmarked[!duplicated(df.unmarked$territory.name, df.unmarked$year.resighted), ]

df.marked <- df[!df$ID %in% not.banded,] # marked df
df.marked <- df.marked[!is.na(df.marked$ID),]
# subset marked birds by the max stage
# for each bird and for each year
df.marked.max <- df.marked %>% 
  group_by(data.frame(ID, year.resighted)) %>%
  slice(which.max(Stage))

#**********************
#* 2. Extract counts 
#**********************
# 'B NB N' - 0= Not seen, 1= Nestling, 2=Nonbreeder, 3= Breeder
counts.marked <- table(df.marked.max$Stage, df.marked.max$year.resighted, df.marked.max$current.population)

counts.unmarked <- table(df.unmarked.nodups$Stage, df.unmarked.nodups$year.resighted, df.unmarked.nodups$current.population)
dimnames(counts.marked)[[1]] <- c('nestling', 'nonbreeder', 'breeder')
dimnames(counts.unmarked)[[1]] <- c('nonbreeder', 'breeder')

#**********************
#* 3. Territory occupancy 
#**********************
# maybe needs work
# data missing zeroes for surveyed but unoccupied sites
occ <- tapply(df$territory.name, list(df$territory.name, df$year.resighted), 
       function(x) {ifelse(length(x)>0,1,NA)}, default=NA)
occ <- occ[order( as.numeric(rownames(occ)) ), ]

#**********************
#* 2. Clean up dataset ----
#**********************
# Ridgway's Hawk phenology notes
# nest construction Feb-March, monitoring Jan-Jul
# Eggs- mid-Feb first clutch mean lay date
# Nestlings- early Mar-late Dec.
# Fledglings- mid-April to May first nests, or July for 3rd nest attempts
# Cutoff dates for annual cycle
fyr <- min(df.marked.max$year.resighted, na.rm=T)
lyr <- max(df.marked.max$year.resighted, na.rm=T)
df.cut <- data.frame(
  start = paste0(fyr:(lyr-1), "-07-14"), # cutoffs around first fledging
  end = paste0((fyr+1):lyr, "-07-14"),
  nms = paste0("Y", fyr:(lyr-1),"_", (fyr+1):lyr) )
df.cut$yr_int <- df.cut$start %--% df.cut$end 

nyr <- lyr-fyr+1
nsite <- length(unique(df$current.population))
#**********************
#* 3. Mark-recapture/resight ----
#**********************
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
#   Sex
#     1 female
#     2 male

# remove individuals only detected during last yr
min.yr.ind <- tapply(df.marked.max$year.resighted, df.marked.max$ID, min, na.rm=T)
last.year.inds <- rownames(min.yr.ind)[min.yr.ind==lyr]
dfm <- df.marked.max[!df.marked.max$ID %in% last.year.inds, ]

nind <- length(unique(dfm$ID))
ind.nu <- sort(unique(dfm$ID))

# detected/not detected
det <- tapply(dfm$ID, list(dfm$ID, dfm$year.resighted), 
              function(x) {ifelse(length(x)>0, 1, 0)}, default=0)
# by stage
y <- tapply(dfm$Stage, list(dfm$ID, dfm$year.resighted), 
              max,  na.rm=T, default=0)
y[y==0] <- 4 # assign 'not seen' as a four

# create z data when latent states are known
z <- array(NA, dim = dim(y), dimnames=dimnames(y))
z[] <- ifelse(y %in% c(2,3), y, NA)

# site covariate 
pop <- tapply(as.numeric(factor(dfm$current.population)), 
              list(dfm$ID, dfm$year.resighted), 
            max,  na.rm=T, default=NA)
pop[pop=='-Inf'] <- NA
# fill in NAs during unobserved years
for (i in 1:nind){
  vals <- which( !is.na(pop[i,]) )
  novals <- which( is.na(pop[i,]) )
  if (length(novals)==0) { next } else{
  for (j in 1:length(novals) ){
    min.val <-  vals[ which.min(vals-novals[j]) ]
    max.val <-  vals[ which.max(vals-novals[j]) ]
    if (novals[j] < max.val){
    # which value is nearest in time
    pop[ i, novals[j] ] <- pop[i,min.val] 
    } else {
      pop[ i, novals[j] ] <- pop[i, max.val]
    } # else
  } # j
  } # else
} # i

# did any birds switch sites?
pop[apply(pop, 1, function(x){ y <- min(x); any(x!=y) } ),]
# none switched

pop2 <- pop
pop2[y==4]
lpop <- melt(pop)
colnames(lpop) <- c("ID", "year", "site")
lpop$year2 <- lpop$year-2010
surv.end <- table(lpop$year2, lpop$site)
yrind.surv <- array(NA, dim = c(max(surv.end), nyr, nsite) )
for (t in 1:nyr){
  for (s in 1:nsite){
    spop <- pop[,t]
    if(surv.end[t, s]==0){next} else{
      yrind.surv[1:surv.end[t, s],t,s] <- which(spop==s)
    } # else
  }} # s t

#**********************
#* 4. Calculate minimum ages ----
#**********************
firstfun <- function(x){ min(which(x %in% c(1,2,3))) }
frst <- apply(y, 1, firstfun)
frst <- ifelse(frst=="Inf" | frst=="-Inf", nyr, frst)

lastfun <- function(x){ max(which(x %in% c(1,2,3))) }
lst <- apply(y, 1, lastfun)
lst <- ifelse(lst=="Inf" | lst=="-Inf", nyr, lst)

age2 <- age <- array(NA, dim=c(nind, nyr),
             dimnames=list(ind.nu, fyr:lyr ))
# set ages for first year of detection
for (i in 1:nind){
  if(y[i, frst[i]]==1){
    age[i, frst[i]] <- 0
    } else{
      age[i, frst[i]] <- 2 # assign 2 if not newborn
    } # else
  } # i
# project ages into future
  for (i in 1:nind){
    for (t in (frst[i]+1):nyr){ 
      if (frst[i] != nyr) {
      age[i, t] <- age[i, t-1] + 1
      }
    }}
# project ages into the past
for (i in 1:nind){
  for (t in (frst[i]):1){ 
    if (frst[i] != 1) {
      age[i, t-1] <- age[i, t] - 1
    }
  }}

for (i in 1:nind){
  for (t in (frst[i]:lst[i]) ){
  age2[i,t] <- ifelse( age[i,t]>=0, age[i,t], NA)
    }}
    
mage <- melt(age2)
colnames(mage) <- c("ID", "Year", "Age")
mage <- mage[!is.na(mage$Age),]
mage$Age <- factor(mage$Age, levels=0:8)
ggplot(mage, aes(fill=Age, y=as.numeric(Age), x=Year)) + 
  geom_bar(position="fill", stat="identity") +
  ylab("Proportion of population")

#**********************
#* reassign stages to 
#* include SY stage
#**********************
y2 <- y
y2[y2==4] <- 6 # not seen becomes state 6
y2[y2==3] <- 5 # breeder becomes 5
y2[y2==2] <- 4 # nonbreeder becomes 4

# age==1 then reassign as FY
y2[age==1 & y==2] <- 2 # FY nonbrreder 
y2[age==1 & y==3] <- 3 # FY breeder

# not many FYs documented
table(y2)

# Need to do the same to zs
# NEEDS WORK
z2 <-  z   

#**********************
#* 7. Productivity data ----
#**********************
# Productivity from mark resight data
# Stage=1 is a nestling
df$prod <- ifelse(df$Stage==1, 1, 
                  ifelse(df$Stage==3, 0, NA))
prod <- tapply(df$prod, list(df$territory.name, df$year.resighted), sum, na.rm=TRUE)
table(prod)

df$group <- ifelse(df$group=="NR", NA, as.numeric(df$group) )
levels(factor(df$group))
treat <- tapply(df$group, list(df$territory.name, df$year.resighted), max, na.rm=TRUE)
treat <- ifelse(treat==3, 1, ifelse(treat==1,0,NA))
# plot treatments over time
plot(2011:2023, colSums(treat, na.rm=T), 
     type="b", xlab="Year", ylab="Number treated",
     main="All Sites")
# plot just Los Haitises
keepersLH <- df$territory.name[df$current.population=="LHNP"]
treatLH <- treat[rownames(treat) %in% keepersLH,]
plot(2011:2023, colSums(treatLH, na.rm=T), 
     type="b", xlab="Year", ylab="Number treated",
     main="Los Haitises")
# plot punta cana
keepersPC <- df$territory.name[df$current.population=="PC"]
treatPC <- treat[rownames(treat) %in% keepersPC,]
plot(2011:2023, colSums(treatPC, na.rm=T), 
     type="b", xlab="Year", ylab="Number treated",
     main="Punta Cana")
# make long form to speed up bc many NAs
lprod <- melt(prod)
ltreat <- melt(treat)
lp <- merge(lprod, ltreat, by=c("Var1", "Var2"))
colnames(lp) <- c("ter","year", "fledged", "treat")
lp <- lp[!is.na(lp$fledged) & !is.na(lp$treat),]

lp$nestsuccess <- ifelse(lp$fledged>0, 1, 0) 
lp$brood <- ifelse(lp$fledged>0, lp$fledged, NA)
lp$year2 <- lp$year-(min(lp$year)-1)

lb <- lp[!is.na(lp$brood),]

# add sites
sites <- df[!duplicated(df$territory.name), c("current.population", "territory.name")]
lp <- merge(lp, sites, by.x="ter", by.y="territory.name")
lp$site2 <- as.numeric(factor(lp$current.population))

lb <- merge(lb, sites, by.x="ter", by.y="territory.name")
lb$site2 <- as.numeric(factor(lb$current.population))
# create matrixes to summarize by year and site
brood.end <- table(lb$year2, lb$site2)
yrind.brood <- array(NA, dim = c(max(brood.end), nyr, nsite) )
for (t in 1:nyr){
  for (s in 1:nsite){
    if(brood.end[t, s]==0){next} else{
  yrind.brood[1:brood.end[t, s],t,s] <- which(lb$year2==t & lb$site2==s)
} # else
  }} # s t

nest.end <- table(lp$year2, lp$site2)
yrind.nest <- array(NA, dim = c(max(nest.end), nyr, nsite) )
for (t in 1:nyr){
  for (s in 1:nsite){
    if(nest.end[t, s]==0){next} else{
  yrind.nest[1:nest.end[t, s],t,s] <- which(lp$year2==t & lp$site2==s)
} # else
  }} # s t

#**********************
#* 8. Hacked birds ----
#**********************
# History
# LHNP=wild and Adult birds in Los Haitises, 
# WPC=wild fledged in Punta Cana, 
# HPC=hacked in Punta Cana (brought from LHNP), 
# FPC=fostered in Punta Cana (brought from LHNP), 
# PC = Adults birds in Puntacana 
# HPS= hacked in Pedro Sanchez (brought from LHNP),

# Hacked birds might not matter if we assume one big population
# ie no immigration or emigration
# Fostered birds will matter because they will directly influence 
# productivity overall
# or try individual data for comparison
df.hacked <- df[df$history == "Hack" & df$Stage==1,]
hacked.to <- tapply(df.hacked$ID, list(df.hacked$year.resighted, df.hacked$current.population), length, default=0)
hacked.from <- tapply(df.hacked$ID, list(df.hacked$year.resighted, df.hacked$origin), length, default=0)
hacked <- data.frame(yr=2011:2023)
hacked <- merge(hacked, hacked.from, by.x="yr", by.y=0, all.x=T)
hacked <- merge(hacked, hacked.to, by.x="yr", by.y=0, all.x=T)
hacked[is.na(hacked$LHNP),"LHNP"] <- 0
hacked[is.na(hacked$PC),"PC"] <- 0
hacked$LHNP <- hacked$LHNP*-1

#create a hacked or not hacked covariate for each individual
# this places fostered eggs into wild birds
hacked.cov.survival <- ifelse(rownames(y) %in% df.hacked$ID, "hacked", "wild") |>
  factor(levels=c("wild","hacked"))
names(hacked.cov.survival) <-  rownames(y)

# check for changes in hacked status over time to
# check data
tapply(dfm$history, list(dfm$ID, dfm$year.resighted), min)
#**********************
#* 9. Fostered ----
#**********************
#* Brought from los haistes to puntacana
#* only one bird for now
df.foster <- df[df$history =="Foster",]

#**********************
#* 10. Covariates ----
#**********************
p <- p2 <- 10
counts <- array(NA, dim=dim(counts.marked), dimnames=dimnames(counts.marked))
counts[1,,] <- counts.marked[1, , ]
counts[c(2,3),,] <- counts.marked[c(2,3), ,] + counts.unmarked[c(1,2), ,]
#########################################
# FY hacked birds need to be moved from PC back 
# to LHNP where they were born
# and subtracted from PC

# add to LHNP
counts[1,,1] <- counts[1,,1] + hacked$LHNP*-1
# subtract from PC
counts[1,,2] <- counts[1,,2] - hacked$PC
#####################################
# from IPMbook package
dUnif <- function (lower, upper) 
{
  A <- round(lower)
  B <- round(upper)
  nrow <- length(lower)
  out <- matrix(0, nrow = nrow, ncol = max(B))
  n <- B - A + 1
  for (i in 1:nrow) out[i, A[i]:B[i]] <- rep(1/n[i], n[i])
  return(drop(out))
}
pPrior <- array(NA, dim=c(max(counts)*2,2))
s.end <- c(max(counts[,,1])*2, max(counts[,,2])*2)
pp1 <- dUnif(1, s.end[1])
pp2 <- dUnif(1, s.end[2])
pPrior[1:s.end[1],1] <- pp1
pPrior[1:s.end[2],2] <- pp2

#**********************
#* 11. Data list for analysis ----
#**********************
datl <- list( # productivity data
              brood = lb$brood,
              nest.success = lp$nestsuccess,
              # survival data
              y = y,
              #z = z, 
              mu.zeroes = rep(0,p),
              mu.zeroes2 = rep(0,p2),
              # count data
              counts = counts
              )

constl <- list( # survival 
                first = frst, 
                nind = nind,
                nyr = nyr,
                site = pop, 
                nsite = nsite,
                yrind.surv = yrind.surv,
                surv.end = surv.end,
                # productivity 
                treat.brood = lb$treat,
                nbrood = nrow(lb),
                year.brood = lb$year2, 
                yrind.brood = yrind.brood,
                brood.end = brood.end,
                site.brood = as.numeric(factor(lb$current.population)), 
                mintrunc = min(lb$brood),
                maxtrunc = max(lb$brood),
                
                treat.nest = lp$treat,
                nnest = nrow(lp),
                year.nest = lp$year2,
                yrind.nest = yrind.nest,
                nest.end = nest.end,
                site.nest = as.numeric(factor(lp$current.population)),
                pPrior = pPrior,
                p = p, # number of random yr effects
                p2 = p2, # number of random yr effects
                
                s.end=s.end,
                hacked = as.numeric(hacked.cov.survival)-1,
                hacked.counts = hacked[, -1]
                # #pop_f = as.numeric(factor(dfp$site)),
)

#*******************
#* Initial values
#*******************
# create initial values for missing y data
get.first <- function(x) min(which(x!=4), na.rm=T)
f <- apply(datl$y, 1, get.first)
get.last <- function(x) max(which(x!=4), na.rm=T)
l <- apply(datl$y, 1, get.last)

# These inits worked but excludes direct input of
# z as data
z.inits <- array(NA, dim=dim(datl$y), dimnames=dimnames(datl$y))
for (i in 1:constl$nind){
  if(l[i]==constl$nyr) { next } else{
    z.inits[i, (l[i]+1):constl$nyr] <- 4 }}
TFmat <- is.na(z.inits) & is.na(z)
for (i in 1:nrow(TFmat) ){  TFmat[i,1:f[i]] <- FALSE }
# For "not seen" this part subs in 
# 2s or 3s depending on which was last observed
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
z.inits[datl$y %in% c(1:3)] <- datl$y[datl$y %in% c(1:3)]
# combine z.inits and z so we can drop z from data 
z.inits <- ifelse(is.na(z.inits), z, z.inits)

# # create inits for rhos
Ustar <- array( runif(p*p, 0.01, 0.1), dim=c(p,p))
diag(Ustar) <- 1 # set diagonal to 1
Ustar[lower.tri(Ustar)] <- 0 # set lower diag to zero
t(Ustar)%*%Ustar

# Abundance
NFYi <- datl$counts[1,,]
NFi <- datl$counts[2,,] + 10
NBi <- datl$counts[3,,] + 16

# In n~binom(p,N) ensure that is always n<N
# tedious but necessary
# FYs
N22 <-  rbind( c(0, 2), round(NFYi[1:(nyr-1),]/2) )
N55 <- rbind( c(0, 2), NFYi[1:(nyr-1),]  )-N22
# Nonbreeders
N33 <- rbind( c(4, 4), round(NFi[1:(nyr-1),]/2) )
N66 <- rbind( c(10, 10), NFi[1:(nyr-1),])-N33
# breeders
N77 <- rbind( c(36, 12) , NBi[1:(nyr-1),] )

#N11_1 <- ifelse(datl$counts[1,,1]<hacked.from, hacked.from, datl$counts[1,,1])
#N11_2 <- cbind(N11_1, datl$counts[1,,2])
                                
N <- array(NA, dim=c(7, constl$nyr, constl$nsite) )
N[1,,] <- datl$counts[1,,1]
N[2,,] <- N22
N[3,,] <- N33
N[4,,] <- 0
N[5,,] <- N55
N[6,,] <- N66
N[7,,] <- N77

N[3,2:4,1] <- 3
N[6,2:4,1] <- 3

NFYi2 <-  N[1,,]
NFi2 <-  apply(N[2:4,,], c(2,3), sum)
NBi2 <-  apply(N[5:7,,], c(2,3), sum)

# Lots of trial and error because of cascading numbers
# Make a function to check for valid inits
Ncheck <- function(x, NFYi2, NFi2, NBi2){
Nc <- array(NA, dim=dim(x))
for (s in 1:constl$nsite){
  for(t in 1: (nyr-1) ){
    Nc[2,t,s] <-  N[2,t+1,s] <= NFYi2[t,s]
    Nc[3,t,s] <-  N[3,t+1,s] <= NFi2[t,s]
    Nc[4,t,s] <-  N[4,t+1,s] <= NBi2[t,s]
    
    Nc[5,t,s] <-  N[5,t+1,s] <= NFYi2[t,s]
    Nc[6,t,s] <-  N[6,t+1,s] <= NFi2[t,s]
    Nc[7,t,s] <-  N[7,t+1,s] <= NBi2[t,s]
  } }
return (Nc)
}

Ncheck(N, NFYi2=NFYi2, NFi2=NFi2, NBi2=NBi2)

N[7,2,1] <- 40
N[7,6,1] <- 70

N[3,6,2] <- 5
N[6,6,2] <- 5
N[7,2,2] <- 16
N[7,9,2] <- 32

NFYi2 <-  N[1,,]
NFi2 <-  apply(N[2:4,,], c(2,3), sum)
NBi2 <-  apply(N[5:7,,], c(2,3), sum)
Ncheck(N, NFYi2=NFYi2, NFi2=NFi2, NBi2=NBi2)

N[3,7,2] <- 4

NFYi2 <-  N[1,,]
NFi2 <-  apply(N[2:4,,], c(2,3), sum)
NBi2 <-  apply(N[5:7,,], c(2,3), sum)
Ncheck(N, NFYi2=NFYi2, NFi2=NFi2, NBi2=NBi2)

N[3,8,2] <- 4
N[6,8,2] <- 4
NFYi2 <-  N[1,,]
NFi2 <-  apply(N[2:4,,], c(2,3), sum)
NBi2 <-  apply(N[5:7,,], c(2,3), sum)
Ncheck(N, NFYi2=NFYi2, NFi2=NFi2, NBi2=NBi2)

N[7,9,2] <- 29
NFYi2 <-  N[1,,]
NFi2 <-  apply(N[2:4,,], c(2,3), sum)
NBi2 <-  apply(N[5:7,,], c(2,3), sum)
Ncheck(N, NFYi2=NFYi2, NFi2=NFi2, NBi2=NBi2)  
  
# check that counts are < totals
counts[2,,1] < NFi2[,1]
counts[2,,2] < NFi2[,2]
counts[3,,1] < NBi2[,1]
counts[3,,2] < NBi2[,2]

# We need good initial values to prevent chains from getting stuck. 
# Use results from a non-integrated analysis get better starting values 
library ('MCMCvis')
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm-simp.rdata")
mus <- MCMCsummary(post, params = "mus", 
            digits=2, HPD = T, 
            hpd_prob = 0.85, pg0= TRUE)
sds <- MCMCsummary(post, params = "sds", 
            digits=2, HPD = T, 
            hpd_prob = 0.85, pg0= TRUE)
betas <- MCMCsummary(post, params = "betas", 
                   digits=2, HPD = T, 
                   hpd_prob = 0.85, pg0= TRUE)

repro <- MCMCsummary(post, 
                          params = c("lmu.brood", "delta",
                          "sig.brood",
                          "mu.nest", "gamma"), 
                          digits=2, HPD = T, 
                          hpd_prob = 0.85, pg0= TRUE)

inits.func1 <- function (){
  list(  
  # fecundity inits
  lmu.brood = c(repro$mean[1], repro$mean[1]),
  delta = repro$mean[2], 
  sig.brood = repro$mean[3],
  mu.nest = c(repro$mean[4], repro$mean[4]),
  gamma = repro$mean[5], 
  # survival
  z = z.inits, 
  mus = cbind(mus$mean, mus$mean), # values from non-integrated run
  betas = betas$mean,
  sds = sds$mean,
  Ustar = Ustar,
  sds2 = rexp(p2, 3),
  Ustar2 = Ustar,
  # counts
  N = N,
  NFY = NFYi2,
  NF = NFi2,
  NB = NBi2, 
  Ntot = NFYi2+NFi2+NBi2, 
  r = rexp(1)
  )}

par_info <- # allows for different seed for each chain
  list(
    list(seed=1, inits = inits.func1()),
    list(seed=2, inits = inits.func1()),
    list(seed=3, inits = inits.func1()),
    list(seed=4, inits = inits.func1())
  )

inits.func <- function() {
  list(  
  # fecundity inits
  # lmu.brood = runif(1, 0, 0.5),
  # delta = runif(1, -0.5, 0.5), 
  # sig.brood = rexp(1),
  # sig.brood.t = rexp(1),
  # sig.nest = rexp(1),
  # mu.nest = runif(1),
  # gamma = runif(1, -0.5, 0.5), 
  # survival
  z = z.inits
  # mus = runif(8),
  # betas = runif(8,-0.5, 0.5), 
  # N=N,
  # NFY = datl$counts.marked[1,],
  # NF = datl$counts.marked[2,],
  # NB = datl$counts.marked[3,]
) }

save(datl, constl, par_info, inits.func, z,
     file="data\\data.rdata")

#*********************
#* Data checks
#*********************

# Check count data
colSums(counts.marked[2:3,,1]) + colSums(counts.marked[2:3,,1])
colSums(counts.marked[2:3,,2]) + colSums(counts.marked[2:3,,2])
#colSums(counts.marked[2:3,,3]) + colSums(counts.marked[2:3,,3])
colSums(table(df$ID, df$year.resighted))

# Check survival data
# How many transitions were observed?
# First year to breeder
# includes unobserved gaps eg 1,4,3
# 1 and 3 not necessarily consecutive
FYB <- FYBT <- c()
ones <- y==1
twos <- y==2
threes <- y==3
fours <- y==4

for (i in 1:nind){
  if ( any(ones[i,]) & any(threes[i,]) ){
    o <- which(ones[i,])
    tw <- which(twos[i,])
    th <- which(threes[i,])
    f <- which(fours[i,])
    if( length(tw)>0) {
       if( any(tw < th) ){ # If there are twos we need to check if its position is btw 1 and 3
         FYB[i] <- 0 
         FYBT[i] <- NA 
         o <- tw <- th <- f <- c()
         next
       }} else{ FYB[i] <- 1
                FYBT[i] <- min(th)
                o <- tw <- th <- f <- c()
               next }
  } else { FYB[i] <- 0
            FYBT[i] <- NA
            o <- tw <- th <- f <- c()
            next }
  }
sum(FYB)
y[which(FYB==1),]
dfFY <- data.frame(trans=FYB, yr=FYBT)

# First year to breeder, observed with no gaps
FYB2 <- c()
for (i in 1:nind){
  o <- ifelse( length(which(y[i,]==1))==0,
              nyr, which(y[i,]==1))
  if (o!=nyr){
    FYB2[i] <- ifelse( y[i,o+1]==3, 1, 0)
  } else{ FYB2[i] <- 0 }
} # i
sum(FYB2)
y[which(FYB2==1),]

# First year to nonbreeder, observed with no gaps
FYA <- c()
for (i in 1:nind){
  o <- ifelse( length(which(y[i,]==1))==0,
               nyr, which(y[i,]==1))
  if (o!=nyr){
    FYA[i] <- ifelse( y[i,o+1]==2, 1, 0)
  } else{ FYA[i] <- 0 }
} # i
sum(FYA)
y[which(FYA==1),]

# Nonbreeder to Breeder, observed with no gaps
AB <- c()
for (i in 1:nind){
  o <- ifelse( length(which(y[i,]==2))==0,
               nyr, which(y[i,]==2))
  if (o!=nyr){
    AB[i] <- ifelse( y[i,o+1]==3, 1, 0)
  } else{ AB[i] <- 0 }
} # i
sum(AB)
y[which(AB==1),]

# Nonbreeder to Breeder, observed with no gaps
BA <- c()
for (i in 1:nind){
  o <- ifelse( length(which(y[i,]==3))==0,
               nyr, which(y[i,]==3))
  if (o!=nyr){
    BA[i] <- ifelse( y[i,o+1]==2, 1, 0)
  } else{ BA[i] <- 0 }
} # i
sum(BA)
y[which(BA==1),]

# compare 
sum(FYA)
sum(FYB2)
sum(AB)
sum(BA)
