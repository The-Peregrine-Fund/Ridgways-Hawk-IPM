#************************************
#* Data manipulation script for 
#* Ridgway's Hawk IPM in DR and Haiti
#* 
#* Author: Brian Rolek
#* Email: rolek.brian@peregrinefund.org
#* Organization: The Peregrine Fund
#************************************
# hacking = # hacked
# suvccess= successful hack

# Link history from banding data 
# to territory name, number and year combined in nest success
# los anianas vargas (2018), los brazos (hack site 2019), 
# los cameita (2023), and punta cana = all hacked to start

#**********************
#* 1. Load packages and data ----
#**********************
library('lubridate')
library('readxl')
dfl <- "data\\RIHA data 2015-2017_v2.xlsx"
df <- read_excel(dfl, sheet='2015_2017_Banding data', na=c('#N/A', 'n/a','N/A', '', 'NA', '-') ) # live/dead tibble
# separate unmarked birds
df.um <- df[df$OID==999,] # unmarked df
df.m <- df[df$OID!=999 & df$Code...10!='none',] # marked df
df.m <- df.m[!is.na(df.m$OID),]

df.m.nodups <- df.m[!duplicated(df.m$OID), ]

#**********************
#* 2. Extract counts 
#**********************
# 'B NB N' - 0= Not seen, 1= Nestling, 2=Nonbreeder, 3= Breeder
counts.marked <- table(df.m.nodups$'B NB N', df.m.nodups$'Year banded/Resighted', df.m.nodups$Sex)
counts.unmarked <- table(df.um$'B NB N', df.um$'Year banded/Resighted', df.um$Sex)

#**********************
#* 3. Territory occupancy 
#**********************
# maybe needs work
# data missing zeroes for surveyed but unoccupied sites
occ <- tapply(df$Territoy_number, list(df$Territoy_number, df$`Year banded/Resighted`), 
       function(x) {ifelse(length(x)>0,1,0)}, default=0)

#**********************
#* 2. Clean up dataset ----
#**********************
# Ridgway's Hawk phenology notes
# nest construction Feb-March, monitoring Jan-Jul
# Eggs- mid-Feb first clutch mean lay date
# Nestlings- early Mar-late Dec.
# Fledglings- mid-April to May first nests, or July for 3rd nest attempts
# Cutoff dates for annual cycle
fyr <- min(df.m$`Year banded/Resighted`, na.rm=T)
lyr <- max(df.m$`Year banded/Resighted`, na.rm=T)
df.cut <- data.frame(
  start = paste0(fyr:(lyr-1), "-07-14"), # cutoffs around first fledging
  end = paste0((fyr+1):lyr, "-07-14"),
  nms = paste0("Y", fyr:(lyr-1),"_", (fyr+1):lyr) )
df.cut$yr_int <- df.cut$start %--% df.cut$end 

nyr <- lyr-fyr+1
nind <- length(unique(df.m$OID))
ind.nu <- sort(unique(df.m$OID))

#**********************
#* 3. Exposure lengths in days ----
#**********************
# From APLO
# Mark-resight-recovery data
#   Observations (po) = y  
#     1 seen first year (age 0)
#     2 seen nonbreeder
#     3 seen breeder
#     4 recovered dead
#     5 not seen
#   States (ps)
#     1 alive first year
#     2 alive nonbreeder
#     3 alive breeder
#     4 dead
#     5 dead not recovered
#     6 emigrated alive
#     7 emigrated dead

# detected/not detected
det <- tapply(df.m$OID, list(df.m$OID, df.m$`Year banded/Resighted`), 
              function(x) {ifelse(length(x)>0, 1, 0)}, default=0)
# by stage
y <- tapply(df.m$'B NB N', list(df.m$OID, df.m$`Year banded/Resighted`), 
              max,  na.rm=T, default=0)

#**********************
#* 4. Calculate ages ----
#**********************
firstfun <- function(x){ min(which(x==1)) }
frst <- apply(det, 1, firstfun)
age <- array(NA, dim=c(nind, nyr),
             dimnames=list(ind.nu, 2015:2017 ))
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
      if (frst[i]!= nyr)
      age[i, t] <- age[i, t-1] + 1
    }}

#**********************
#* 7. Productivity data ----
#**********************
dflp <- "data\\RIHA data 2015-2017.xlsx"
dfp <- read_excel(dfl, sheet='2015_2017_nest summary', na=c('#N/A', 'n/a','N/A', '', 'NA', '-') ) # live/dead tibble
dfp$fledged <- ifelse(dfp$`# Fledged`!="NR", as.numeric(dfp$`# Fledged`), NA)
fl <- tapply(dfp$fledged, list(dfp$`Territory Number`, dfp$Year), sum, na.rm=T)
colSums(fl, na.rm=T)
colMeans(fl, na.rm=T)

# was a nest treated with insecticide?
# 1=No, 3=yes
dfp$tr <- ifelse(dfp$Group=="NR", NA,  
                 ifelse(dfp$Group=="3", 1, 0))
tr <- tapply(dfp$tr , list(dfp$`Territory Number`, dfp$Year), max, na.rm=T)

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
# 
dfp$num.hacked <- ifelse(dfp$Hacking=="NR", NA, as.numeric(dfp$Hacking))
tapply(dfp$num.hacked, list(dfp$Year), sum, na.rm=T)

hacked <- dfp[dfp$Hacking!=0,]
dfp$`Territory Number`

#**********************
#* 8. Covariates ----
#**********************


#**********************
#* 9. Data list for analysis ----
#**********************
datl <- list( # productivity data
              prod_nest=fl,
              prod_tot=colSums(fl, na.rm=T),
              # survival data
              y=y,
              # count data
              # indices

              # deterministic data
              )

constl <- list( # survival 
                sex=ld.nodups$sex, # manually check "unknowns" and "NA"s
                transl=, # translocated
                nind=nind,
                nyr=nyr,
                # productivity 
                fledged=dfp$fledged,
                nnests=nrow(dfp),
                nter=,
                fosters=, # number fostered
                site_f=as.numeric(factor(dfp$Population)), # 
                tr= dfp$tr, # nest treatment

                J=

)