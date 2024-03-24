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
dfl <- "data\\RIHA data 2015-2017_v2.xlsx"
df <- read_excel(dfl, sheet='2015_2017_Banding data', na=c('#N/A', 'n/a','N/A', '', 'NA', '-', 
                          "None", "Not Recorded", "Not recorded") ) # live/dead tibble
df$ID <- paste(df$Code...8, df$Code...10, 
               df$`Left Color`, df$`Right color`, sep="-")

# separate unmarked birds
not.banded <- c("NA-NA-NA-NA", "NA-NA-No band-No band")
df.unmarked <- df[df$ID %in% not.banded,] # unmarked df

df.marked <- df[df$Sex=="Female",]
df.marked <- df.marked[!df.marked$ID %in% not.banded,] # marked df
df.marked <- df.marked[!is.na(df.marked$ID),]
df.marked.nodups <- df.marked[!duplicated(df.marked$ID, df.marked$`Year banded/Resighted`), ]

#**********************
#* 2. Extract counts 
#**********************
# 'B NB N' - 0= Not seen, 1= Nestling, 2=Nonbreeder, 3= Breeder
counts.marked <- table(df.marked.nodups$'B NB N', df.marked.nodups$'Year banded/Resighted', df.marked.nodups$Sex)
counts.unmarked <- table(df.unmarked$'B NB N', df.unmarked$'Year banded/Resighted', df.unmarked$Sex)

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
fyr <- min(df.marked$`Year banded/Resighted`, na.rm=T)
lyr <- max(df.marked$`Year banded/Resighted`, na.rm=T)
df.cut <- data.frame(
  start = paste0(fyr:(lyr-1), "-07-14"), # cutoffs around first fledging
  end = paste0((fyr+1):lyr, "-07-14"),
  nms = paste0("Y", fyr:(lyr-1),"_", (fyr+1):lyr) )
df.cut$yr_int <- df.cut$start %--% df.cut$end 

nyr <- lyr-fyr+1
nind <- length(unique(df.marked$ID))
ind.nu <- sort(unique(df.marked$ID))

#**********************
#* 3. Exposure lengths in days ----
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

# detected/not detected
det <- tapply(df.marked$ID, list(df.marked$ID, df.marked$`Year banded/Resighted`), 
              function(x) {ifelse(length(x)>0, 1, 0)}, default=0)
# by stage
y <- tapply(df.marked$'B NB N', list(df.marked$ID, df.marked$`Year banded/Resighted`), 
              max,  na.rm=T, default=0)
y[y==0] <- 4 # assign 'not seen' as a four

# create z data when latent states are known
z <- array(NA, dim = dim(y), dimnames=dimnames(y))
z[] <- ifelse(y %in% c(2,3), y, NA)

# population covariate 
pop <- tapply(as.numeric(factor(df.marked$Population)), list(df.marked$ID, df.marked$`Year banded/Resighted`), 
            max,  na.rm=T, default=NA)
pop[pop=='-Inf'] <- NA

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

#**********************
#* 7. Productivity data ----
#**********************
dflp <- "data\\RIHA data 2015-2017.xlsx"
dfp <- read_excel(dfl, sheet='2015_2017_nest summary', na=c('#N/A', 'n/a','N/A', '', 'NA', '-', "NR") ) # live/dead tibble
# omit "NO ID"s
dfp <- dfp[substr(dfp$`Territory Number`, 1, 5)!="No ID", ]
# omit NAs
dfp <- dfp[!is.na(dfp$`# Fledged` ), ]
fl <- tapply(dfp$`# Fledged`, list(dfp$`Territory Number`, dfp$Year), sum, na.rm=T)
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
hacked <- tapply(dfp$Hacking, list(dfp$Year), sum, na.rm=T)

#**********************
#* 8. Covariates ----
#**********************


#**********************
#* 9. Data list for analysis ----
#**********************
datl <- list( # productivity data
              f=dfp$`# Fledged`,
              # survival data
              y=y,
              # count data
              countA = counts.marked[3,,1], # counts.marked[stage, year, sex]
              countB = counts.marked[4,,1], 
              countFYW = counts.marked[2,,1]
              )

constl <- list( # survival 
                transl=, # translocated
                nind=nind,
                nyr=nyr,
                # productivity 
                nyr=nyr,
                K=nrow(dfp),
                year=as.numeric(factor(dfp$Year)), 
                tr=dfp$tr,
                fosters=, # number fostered
                site_f=as.numeric(factor(dfp$Population)),
                NH = hacked
)
