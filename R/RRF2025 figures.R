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

#*******************
#### ---- pva_ext1 -----
#* Plot PVA results
#*******************
library ('viridis')
library ("tidybayes")
library ('coda')
library ('MCMCvis')
library ('reshape2')
library('ggplot2')
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\pva_shortrun100k.rdata")
out <- lapply(post[-5], as.mcmc) 
outp <- MCMCpstr(out, type="chains")
scen <- data.frame('Scenarios' = 1:45,
                   'Number translocated' = rep(c(0,5,10), each=15), 
                   'Number of nests treated' = rep( rep(c(0, 15, 30, 45, 100), each=3), 3),
                   'Mortality reduction' = rep(c('100% Mortality', '80% Mortality', '60% Mortality'), 15) )

#*********************
#* Line plot of extinction probability over time
#*********************
pe <- apply(outp$extinct, c(1,2,3), mean)
lpe <- melt(pe)
colnames(lpe) <- c('Scenarios', 'Time', 'Site2', 'Extinct')
pe.df <- merge(lpe, scen, by='Scenarios', all.x=TRUE)
pe.df$Year <- pe.df$Time+2023
pe.df$Site <- ifelse(pe.df$Site2==1, "Los Haitises", "Punta Cana")
pe.df$Mortality.reduction <- factor(pe.df$Mortality.reduction, 
                                    levels= c('100% Mortality', '80% Mortality', '60% Mortality'))

p6 <- ggplot(pe.df, aes(x=Year, y=Extinct, group=Scenarios, 
                        color = factor(Number.translocated), 
                        linetype = factor(Number.of.nests.treated))) +
  geom_line( linewidth = 1) +
  facet_grid(rows= vars(Site, factor(Number.translocated)), 
             cols= vars(Mortality.reduction)) +
  scale_color_viridis_d(option="viridis", direction=1,
                        name="Number translocated") +
  scale_linetype(name="Number of nests treated") +
  scale_x_continuous( name="Future year",
                      breaks=c(2020, 2030, 2040, 2050, 2060, 2070),
                      labels=c('', '2030', '', '2050', '', '2070'), 
                      minor_breaks=NULL) +
  scale_y_continuous( name="Extinction probability",
                      breaks=seq(0, 1, by=0.1),
                      labels=c('0', '', '', '', '', '0.5', '', '', '', '', '1'), 
                      limits=c(0,1),
                      minor_breaks=NULL) +
  theme_minimal()  
p6
ggsave(filename="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Extinction_lines.tiff",
       plot=p6,
       device="tiff",
       width=6, height=6, dpi=300)

#### ---- pva_ext2 -----
#********************
#* PVA heatmap of extinction probability projected 50 years
#********************
pe.df2 <- pe.df[pe.df$Time==50, ]
pe.df2 <- pe.df2[pe.df2$Site=="Los Haitises",]


p7 <- ggplot(pe.df2, aes(x = as.factor(Number.translocated), 
                         y = as.factor(Number.of.nests.treated))) +
  geom_tile( aes(fill = Extinct ) ) + 
  geom_text(aes(label=round(Extinct,2)), color="gray70") +
  scale_fill_viridis_c(option="plasma", direction=1,
                       name = "Extinction\nprobability after\n50 years",
                       begin=0, end=1, limits=c(0,1)) +
  facet_wrap(~  Site + factor(Mortality.reduction, 
                              levels=c("100% Mortality", "80% Mortality", "60% Mortality"))) +
  scale_y_discrete(breaks=c(0,15,30,45,100), 
                   name="Number of territories treated") +
  scale_x_discrete(name="Number of females translocated") +
  theme_minimal( ) + theme(strip.background = element_rect(size = 0.5)) 
  
p7  
ggsave(filename="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Extinction_50yrHeatmap.tiff",
       plot=p7,
       device="tiff",
       width=6, height=4, dpi=300)
pe.df3 <- pe.df2[!is.na(pe.df2$Extinct),]
pe.df3 <- pe.df3[rev(order(pe.df3$Site2, pe.df3$Extinct)),]
pe.df3

#### ---- pva_abund -----
#***************************
#* Plot total Abundance
#***************************
lnt1 <-  outp$Ntot[,1:13, , ] |> melt()
colnames(lnt1) <- c("Scenarios", "Time", "Site2", "Iter", "Abund")
lnt1$Year <- lnt1$Time + 2010
lnt1$Site <- factor(ifelse(lnt1$Site2==1, 'Los Haitises', 'Punta Cana'))
nt.df <- merge(lnt1, scen, by='Scenarios', all.x=TRUE)
nt.df$Mortality.reduction <- factor(nt.df$Mortality.reduction,
                                    levels=c('100% Mortality','80% Mortality', '60% Mortality' ))

# calculate future means
lnt2 <- apply(outp$Ntot[,13:63, , ], c(1, 2, 3), mean) |> melt()
colnames(lnt2) <- c("Scenarios", "Time", "Site2", "Abund")
lnt2$Year <- lnt2$Time + 2022
lnt2$Site <- factor(ifelse(lnt2$Site2==1, 'Los Haitises', 'Punta Cana'))
nt.df2 <- merge(lnt2, scen, by='Scenarios', all.x=TRUE)
nt.df2$Mortality.reduction <- factor(nt.df2$Mortality.reduction,
                                     levels=c('100% Mortality','80% Mortality', '60% Mortality' ))

p8 <- ggplot() +
  stat_lineribbon(data= nt.df,  
                  aes(x=Year, y=Abund, group=Scenarios),
                  fill='gray60',
                  linewidth=0.5,
                  point_interval = 'median_hdci', .width = 0.95,
                  na.rm=TRUE) +
  geom_line(data= nt.df2,
            aes(x=Year, y=Abund, #group=Scenarios,
                color = factor(Number.translocated),
                linetype = factor(Number.of.nests.treated) ),
            linewidth=0.5) +
  facet_grid(rows = vars(Site, Number.translocated), cols = vars(Mortality.reduction),
             scales ='free_y') +
  ylab('Number of females') +
  scale_color_viridis_d(option="viridis", direction=1,
                        name="Number translocated") +
  scale_linetype(name="Number of nests treated") +
  theme_minimal()

p8
ggsave(filename="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Abundance_projections.tiff",
       plot=p8,
       device="tiff",
       width=8, height=6, dpi=300)

#### ---- pva_pc -----
#******************************
#* Future percent change of population
#* between 2023 and 2073
#******************************
# calculate percent change and make a database for plotting
pcl <- ((( outp$Ntot[,63,,]/outp$Ntot[,13,,] )-1)*100) |> melt() 
colnames(pcl) <- c("Scenario", "Site.num", "Iter", "PC")
pc.mns <- tapply(pcl$PC, list(pcl$Scenario, pcl$Site.num), mean) |> melt()
colnames(pc.mns) <- c("Scenario", "Site.num", "PC")
pc.df <- merge(pc.mns, scen, by.x="Scenario", by.y="Scenarios")
pc.df$Site <- ifelse(pc.df$Site.num==1, "Los Haitises", "Punta Cana")
# scale color so that percent change=0 is the middle color (gray),
# -100 is the extreme lower bound and max is the extreme upper bound
pc.df$pc.sc <- NA
for (i in 1:length(pc.df$PC)){
  if(pc.df$PC[i]<=0){
    pc.df$pc.sc[i] <- (pc.df$PC[i]-0)/sd(pc.df$PC[pc.df$PC<=0])}
  else{ pc.df$pc.sc[i] <- (pc.df$PC[i]-0)/sd(pc.df$PC[pc.df$PC>0]) }
}

labs <- c(-90, 0, 700, 1400)
labs.sc <- c( (-90-0)/sd(pc.df$PC[pc.df$PC<=0]),
              (0-0)/sd(pc.df$PC[pc.df$PC<=0]),
              (700-0)/sd(pc.df$PC[pc.df$PC>0]),
              (1400-0)/sd(pc.df$PC[pc.df$PC>0]))

pc.df <- pc.df[pc.df$Site=="Los Haitises", ]
# plot
p11 <- ggplot(pc.df, aes(x = as.factor(Number.translocated), 
                         y = as.factor(Number.of.nests.treated))) +
  geom_tile( aes(fill = pc.sc ) ) + 
  geom_text(aes(label=round(PC,0)), color="gray50") +
  scale_fill_viridis_c(option="turbo", direction=-1,
                       name = "Percent change\nafter 50 years",
                       breaks = labs.sc, labels = labs,
                       begin=0, end=1, limits=c(-2.898255, 4.863298)) +
  facet_wrap(~  Site + factor(Mortality.reduction, 
                              levels=c("100% Mortality", "80% Mortality", "60% Mortality"))) +
  scale_y_discrete(breaks=c(0, 15, 30, 45, 100), 
                   name="Number of territories treated") +
  scale_x_discrete(name="Number of females translocated") +
  theme_minimal( ) + theme(strip.background = element_rect(size = 0.5))
p11 
ggsave(filename="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Abundance_50yrHeatmap.tiff",
       plot=p11,
       device="tiff",
       width=6, height=4, dpi=300)

#### ---- pva_pd -----
#**********************************
# Probability of Future Population Decline
#**********************************
# p.decline <- apply( (outp$Ntot[,63,,]-outp$Ntot[,13,,])<0, c(1,2), mean)
# pdl <- melt(p.decline)
# colnames(pdl) <- c("Scenario", "Site.num", "pd")
# pdl$Site <- ifelse(pdl$Site.num==1, "Los Haitises", "Punta Cana")
# pdl <- merge(pdl, scen, by.x="Scenario", by.y="Scenarios")

p.decline <- tapply(pcl$PC<0, list(pcl$Scenario, pcl$Site.num), mean)
pdl <- melt(p.decline)
colnames(pdl) <- c("Scenario", "Site.num", "pd")
pdl$Site <- ifelse(pdl$Site.num==1, "Los Haitises", "Punta Cana")
pd <- merge(pdl, scen, by.x="Scenario", by.y="Scenarios")

pd <- pd[pd$Site=="Los Haitises", ]
p12 <- ggplot(pd, aes(x = as.factor(Number.translocated), 
                      y = as.factor(Number.of.nests.treated))) +
  geom_tile( aes(fill = pd ) ) + 
  geom_text(aes(label=round(pd,2)), color="gray50") +
  scale_fill_viridis_c(option="magma", direction=1,
                       name = "Prob. of decline\nafter 50 years",
                       breaks = c(0,0.25,0.50,0.75,1.0), 
                       labels = c(0,0.25,0.50,0.75,1.0),
                       begin=0, end=1, limits=c(0, 1)) +
  facet_wrap(~  Site + factor(Mortality.reduction, 
                              levels=c("100% Mortality", "80% Mortality", "60% Mortality"))) +
  scale_y_discrete(breaks=c(0, 15, 30, 45, 100), 
                   name="Number of territories treated") +
  scale_x_discrete(name="Number of females translocated") +
  theme_minimal( ) + theme(strip.background = element_rect(size = 0.5))
p12 
ggsave(filename="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//ProbDecline_50yrHeatmap.tiff",
       plot=p12,
       device="tiff",
       width=6, height=4, dpi=300)

#### ---- pva_tt_ext -----
#*****************************
# Time to extinction
#*****************************
# conditional on extinction threshold (D=X)

D <- 1
abund <- outp$Ntot[,14:(13+50),,]
dimnames(abund) <- list('Scenarios'=1:45, 
                        'Time'= 2024:(2024+49), 
                        'Site'=c("Los Haitises", "Punta Cana"), 
                        'Iter'=1: dim(outp$extinct)[[4]])
qextinct <- abund == D
ext <- qextinct[,50,,] <= D
lte <- melt( qextinct )
colnames(lte)[5] <- "Extinct"
le <- melt( ext )
colnames(le)[4] <- "EverExtinct"
lte2 <- merge(lte, le, by=c("Scenarios", "Site", "Iter")  )
lte3<- lte2[lte2$EverExtinct==TRUE, ]
lte4 <- lte3[lte3$Extinct==TRUE, ]
lte5 <- tapply(lte4$Time-2023, list(lte4$Scenarios, lte4$Site, lte4$Iter), min)
lte6 <- melt(lte5)
lte7 <- lte6[!is.na(lte6$value), ]
colnames(lte7) <- c('Scenarios', 'Site', 'Iter', 'TTE')
te <- table(lte7$TTE, lte7$Scenarios, lte7$Site) |> melt()
colnames(te) <- c('TE', 'Scen', 'Site', 'Count')
te$Scenarios <- factor(paste0("Scenario ", te$Scen), 
                       levels=paste0('Scenario ', 1:45)) 
teLH <- te[te$Site=="Los Haitises",]
tePC <- te[te$Site=="Punta Cana",]

p9 <- ggplot() +
  geom_col(data = teLH, aes(x= TE, y=Count)) + 
  facet_wrap(vars(Scenarios),
             nrow=9, ncol=5, scales="free_y", shrink=FALSE)

p10 <- ggplot() +
  geom_col(data = tePC, aes(x= TE, y=Count)) + 
  facet_wrap(vars(Scenarios),
             nrow=9, ncol=5, scales="free_y", shrink=FALSE)
p9
p10

# plot Abundance siumulations
library (gganimate)

lN <- melt(outp$Ntot[1,,1,1:100])
colnames(lN) <- c("time", "Iteration", "Abundance")
lN$Year <- lN$time+2010

lN1 <- lN[lN$Iteration==1,]
ggplot(lN1, aes(x = Year, y = Abundance, 
               group=Iteration, color=Iteration)) +
  geom_point() +
  geom_line() +
  geom_vline(aes(xintercept=2023), linetype="dashed")+
  scale_color_viridis_c(option="D", alpha=0.5) +
  transition_reveal(Year) + 
  theme_bw(base_size = 22) +
  ylab("Number of females") +
  ylim(c(0, 1200)) 

anim_save("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\figs\\PVA_Abundance-1Iter.gif")


ggplot(lN, aes(x = Year, y = Abundance, 
               group=Iteration, color=Iteration)) +
  geom_point() +
  geom_line() +
  geom_vline(aes(xintercept=2023), linetype="dashed")+
  scale_color_viridis_c(option="D", alpha=0.5) +
  transition_reveal(Year) + 
  theme_bw(base_size = 22) +
  ylab("Number of females") +
  ylim(c(0, 1200))

anim_save("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\figs\\PVA_Abundance.gif")


  transition_states(Year,
                    transition_length = 2,
                    state_length = 1)
