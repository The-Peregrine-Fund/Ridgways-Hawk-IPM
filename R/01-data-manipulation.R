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
df <- df[df$sex=="Female",]
df <- df[df$current.population!="AV",] # remove because small sample size

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
# separate marked and unmarked birds
not.banded <- c("NB-NB-NB-NB", "NB-NR-NB-NR", "NR-NR-NR-NR")
df.unmarked <- df[df$ID %in% not.banded,] # unmarked df
df.unmarked.nodups <- df.unmarked[!duplicated(df.unmarked$territory.number, df.unmarked$year.resighted), ]

df.marked <- df[!df$ID %in% not.banded,] # marked df
df.marked <- df.marked[!is.na(df.marked$ID),]
# subset marked birds by the max stage
# for each bird and for each year
df.marked.max <- df.marked %>% 
  group_by(data.frame(ID,year.resighted)) %>%
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
occ <- tapply(df$territory.number, list(df$territory.number, df$year.resighted), 
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
keepers <- df$territory.name[df$current.population=="LHNP"]
treatLH <- treat[rownames(treat) %in% keepers,]
plot(2011:2023, colSums(treatLH, na.rm=T), 
     type="b", xlab="Year", ylab="Number treated",
     main="Los Haitises only")
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

brood.end <- table(lb$year2)
yrind.brood <- array(NA, dim = c(max(brood.end), nyr) )
for (t in 1:nyr){
  yrind.brood[1:brood.end[t],t] <- which(lb$year2==t)
}

nest.end <- table(lp$year2)
yrind.nest <- array(NA, dim = c(max(nest.end), nyr) )
for (t in 1:nyr){
  yrind.nest[1:nest.end[t],t] <- which(lp$year2==t)
}

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
df.hacked <- df[df$history == "Hack",]
hacked.to <- tapply(df.hacked$ID, list(df.hacked$year.resighted, df.hacked$current.population), length, default=0)
hacked.from <- tapply(df.hacked$ID, list(df.hacked$year.resighted, df.hacked$origin), length, default=0)

# compare individual data 
data.frame(hacked.to, hacked.from)

#create a hacked or not hacked covariate for each individual
# this places fostered eggs into wild birds
hacked.cov.survival <- ifelse(rownames(y) %in% df.hacked$ID, "hacked", "wild") |>
  factor(levels=c("wild","hacked"))
names(hacked.cov.survival) <-  rownames(y)

# check for changes in hacked status over time to
# check data
tapply(dfm$history, list(dfm$ID, dfm$year.resighted), min)
#**********************
#* 8. Fostered ----
#**********************
#* Brought from los haistes to puntacana
#* only one bird for now
df.foster <- df[df$history =="Foster",]

#**********************
#* 8. Covariates ----
#**********************

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

p <- 8
#**********************
#* 9. Data list for analysis ----
#**********************
datl <- list( # productivity data
              brood = lb$brood,
              nest.success = lp$nestsuccess,
              # survival data
              y = y,
              #z = z, 
              mu.zeroes = rep(0,p),
              # count data
              counts.unmarked = apply(counts.unmarked, c(1,2), sum),
              counts.marked = apply(counts.marked, c(1,2), sum)
              )

constl <- list( # survival 
                first = frst, 
                #hacked =  as.numeric(hacked.cov.survival)-1, # 1 = wild, 2 = hacked
                nind = nind,
                nyr = nyr,
                site = pop, 
                # productivity 
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
                nest.end = nest.end,
                pInit = dUnif(1, max(datl$counts.marked)),
                p = p, # number of random yr effects
                hacked = as.numeric(hacked.cov.survival)-1
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
Ustar <- array( runif(p*p, 0.1, 0.5), dim=c(p,p))
diag(Ustar) <- 1 # set diagonal to 1
Ustar[lower.tri(Ustar)] <- 0 # set lower diag to zero
t(Ustar)%*%Ustar

# Abundance
N <- array(NA, dim=c(7, constl$nyr) )
#N[] <- 0
# N[1,] <- datl$counts.marked[1,]
# N[2,] <- round(datl$counts.marked[2,]/2)
# N[3,] <- round(datl$counts.marked[2,]/2)
N[4,] <- 0
# N[5,] <- round(datl$counts.marked[3,]/3)
# N[6,] <- round(datl$counts.marked[3,]/3)
# N[7,] <- round(datl$counts.marked[3,]/3)

#ifelse(N[2,]>datl$counts.marked[2,], datl$counts.marked[2,], N[2,])

inits.func1 <- function (){
  list(  
  # fecundity inits
  lmu.brood = runif(1, 0, 0.5),
  delta = runif(1, -0.1, 0.1), 
  sig.brood = rexp(1),
  sig.brood.t = rexp(1),
  sig.nest = rexp(1),
  mu.nest = runif(1),
  gamma = runif(1, -0.1, 0.1), 
  # survival
  z = z.inits, 
  mus = runif(8),
  betas = runif(8,-0.1, 0.1),
  sds = rexp(p, 6),
  Ustar = Ustar,
  # eta.phiFY = runif(constl$nyr),
  # eta.phiA = runif(constl$nyr),
  # eta.phiB = runif(constl$nyr),
  # eta.psiFYB = runif(constl$nyr),
  # eta.psiAB = runif(constl$nyr),
  # eta.psiBA = runif(constl$nyr),
  # eta.pA = runif(constl$nyr),
  # eta.pB = runif(constl$nyr),
  # counts
  # N = N,
  NFY = datl$counts.marked[1,],
  NF = datl$counts.marked[2,],
  NB = datl$counts.marked[3,]
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
colSums(counts.marked[2:3,,3]) + colSums(counts.marked[2:3,,3])
colSums(table(df$ID, df$year.resighted))

# Check survival data
# How many transitions were observed?
# First year to breeder
# includes unobserved gaps eg 1,4,3
# 1 and 3 not necessarily consecutive
FYB <- c()
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
         o <- tw <- th <- f <- c()
         next
       }} else{ FYB[i] <- 1
                o <- tw <- th <- f <- c()
               next }
  } else { FYB[i] <- 0 
            o <- tw <- th <- f <- c()
            next }
  }
sum(FYB)
y[which(FYB==1),]

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
