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
load("data/data.rdata")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_longrun.rdata")
out <- lapply(post[1:4], as.mcmc) # omit chain 5, causing lack of convergence in 1 param
outp <- MCMCpstr(out, type="chains")

#********************
#* Calculate summary stats
#********************
#* Abundance at both sites
Nall <- apply(outp$Ntot, c(1,3) , sum)
l.all.post <- array(NA, dim=c(constl$nyr-1, ncol(Nall)) ) 
for (t in 2:13){
l.all.post[t-1, ]  <- Nall[t,] / Nall[t-1,]
}
l.mn <- apply(l.all.post, 2, mean)
# Average annual pop growth
median(l.mn)
HDInterval::hdi(l.mn)

# 2011 to 2023 pop growth, growth between years
l.all.post2  <- Nall[13,] / Nall[1,]
median(l.all.post2)
HDInterval::hdi(l.all.post2)

# Population growth rates
mn.lambda <- apply(outp$lambda, c(2,3), mean)
(md.lambda <- apply(mn.lambda, 1, median))
(hdi95.lambda <- apply(mn.lambda, 1, HDInterval::hdi))

# Overall averages
mn.mus <- apply(outp$mus, c(1,3), mean)
(md.mus <- apply(mn.mus, 1, median)) |> round(2)
(hdi.mus <- apply(mn.mus, 1, HDInterval::hdi)) |> round(2)
# Does detection differ by life stage?
pd(mn.mus[8,]-mn.mus[7,])

# Compare between sites
df.sites <- data.frame(mus=1:8, pd=rep(NA,8), greater=NA)
for (i in 1:8){
  first <- median(outp$mus[i, 1, ]) > median(outp$mus[i, 2, ])
  df.sites$greater[i] <- ifelse(first, "LH", "PC")
  if (first){
    diff <- outp$mus[i, 1, ]-outp$mus[i, 2, ]
    df.sites$pd[i] <- mean(diff>0) |> round(2) 
  } else{
    diff <- outp$mus[i, 2, ]-outp$mus[i, 1, ]
    df.sites$pd[i] <- mean(diff>0) |> round(2)
  }
}
df.sites

(md.mus <- apply(outp$mus, c(1,2), median)) |> round(2)
(hdi.mus <- apply(outp$mus, c(1,2), HDInterval::hdi)) |> round(2)
# productivity
# see section below

# Compare survival, recruitment, and detection rates
# among life stages within each site
df.comp <- data.frame(param1 = c('phiA', 'phiA', 'phiB', 'psiA', 'pB'), 
                 param2 = c('phiB', 'phiFY', 'phiFY', 'psiFY', 'pA'), 
                 pd_LH = NA, pd_PC = NA)
## Los Haitises, survival
df.comp$pd_LH[1] <- (outp$lmus[2,1,]-outp$lmus[3,1,]) |> pd() |> round(2)
df.comp$pd_LH[2] <- (outp$lmus[2,1,]-outp$lmus[1,1,]) |> pd() |> round(2)
df.comp$pd_LH[3] <- (outp$lmus[3,1,]-outp$lmus[1,1,]) |> pd() |> round(2)
## Los Haitises, recruitment rates
df.comp$pd_LH[4] <- (outp$lmus[5,1,]-outp$lmus[4,1,]) |> pd() |> round(2)
## Los Haitises, resight probability
df.comp$pd_LH[5] <- (outp$lmus[8,1,]-outp$lmus[7,1,]) |> pd() |> round(2)
## Punta Canas, survival
df.comp$pd_PC[1] <- (outp$lmus[2,2,]-outp$lmus[3,2,]) |> pd() |> round(2)
df.comp$pd_PC[2] <- (outp$lmus[2,2,]-outp$lmus[1,2,]) |> pd() |> round(2)
df.comp$pd_PC[3] <- (outp$lmus[3,2,]-outp$lmus[1,2,]) |> pd() |> round(2)
## Punta Cana, recruitment rates
df.comp$pd_PC[4] <- (outp$lmus[5,2,]-outp$lmus[4,2,]) |> pd() |> round(2)
## Punta Cana, resight probability
df.comp$pd_PC[5] <- (outp$lmus[8,2,]-outp$lmus[7,2,]) |> pd() |> round(2)
df.comp 

#*******************
#* Plot IPM results
#*******************
#Plot population size
lNall <- melt(Nall)
colnames(lNall) <- c("Time", "Iter", "Abundance")
lNall$Site.ind <- 3; lNall$Site <- "Both"
lNall <- lNall[,c(1,4,2,3,5)]

Ntot <- outp$Ntot[1:13,,] # thin by 10 to speed up plotting
lN <- melt(Ntot)
colnames(lN) <- c("Time", "Site.ind", "Iter", "Abundance")
lN$Site <- ifelse(lN$Site.ind==1, "Los Haitises", "Punta Cana")
lN <- rbind(lN, lNall)

lN$Year <- lN$Time + 2010

med_hdis <- function(x, cm){
            df <- data.frame(y=median(x),
                             ymin=HDInterval::hdi(x, credMass=cm)[1],
                             ymax=HDInterval::hdi(x, credMass=cm)[2] )
}

counts.tot <- datl$countsAdults+datl$countsFY+constl$hacked.counts[,-3]
df.counts <- melt(counts.tot)
colnames(df.counts) <- c("Year", "Site_ab", "Count")
df.counts$Site <- ifelse(df.counts$Site_ab=="LHNP", "Los Haitises", "Punta Cana")

df.counts.all <- data.frame(Year=2011:2023, Site_ab="Both", Count=counts.tot[,1] + counts.tot[,2], Site="Both" )
df.counts <- rbind(df.counts, df.counts.all)

p1 <- ggplot() + theme_minimal(base_size=14) + 
  geom_line(data=lN, aes(x=Year, y=Abundance, group=Iter), 
            color="gray30", linewidth=0.1, alpha=0.01 ) +
  stat_summary(data=lN, aes(x=Year, y=Abundance), 
               geom="errorbar", width=0, 
               fun.data =med_hdis, fun.args=list(cm=0.95) ) +
  stat_summary(data=lN, aes(x=Year, y=Abundance), 
               geom="errorbar", width=0, linewidth=1,
               fun.data =med_hdis, fun.args=list(cm=0.80) ) +
  stat_summary(data=lN, aes(x=Year, y=Abundance), 
               geom="line", linewidth=1,
               fun.data =function(x){data.frame(y=median(x))} ) +
  geom_line(data=df.counts, aes(x=Year, y=Count), 
            col="black", linewidth=1, linetype=2 ) +
  facet_wrap("Site", scales="free_y") +
  xlab("Year") + ylab("Number")
p1

# ggsave(filename="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Abundance and counts.tiff",
#        plot=p1,
#        device="tiff",
#        width=9, height=4, dpi=300)

#### ----popstructure -----
# Plot life stage structure
mdFY <-  apply(outp$NFY[1:13,,], c(1,2), median) 
mdB <-  apply(outp$NB[1:13,,], c(1,2), median) 
mdF <-  apply(outp$NF[1:13,,], c(1,2), median) 
lFY <- melt(mdFY)
lB <- melt(mdB)
lF <- melt(mdF)
lFY$Stage <- "First-year"
lB$Stage <- "Breeder"
lF$Stage <- "Nonbreeder"
ldat <- rbind(lFY, lB, lF)
colnames(ldat)[1:3] <- c("Year", "Sitenum", "Number") 
ldat$Site <- ifelse(ldat$Sitenum==1, "Los Haitises", "Punta Cana")
ldat$year <- ldat$Year+2010

# Use median number of females in each stage
# to plot an approximate population structure
p3 <- ggplot(ldat, aes(fill=Stage, y=as.numeric(Number), x=year)) + 
  geom_bar(position="fill", stat="identity") +
  ylab("Proportion of population") + xlab("Year") +
  facet_wrap("Site") + 
  theme_minimal(base_size=14)+ theme(legend.position="top") +
  scale_fill_viridis_d(option = "mako")

p4 <- ggplot(ldat, aes(fill=Stage, y=as.numeric(Number), x=year)) + 
  geom_bar(position="stack", stat="identity") +
  ylab("Median number of females") + xlab("Year") +
  facet_wrap("Site", scales = "free_y") + 
  theme_minimal(base_size=14 ) + theme(legend.position="top") + 
  scale_fill_viridis_d(option = "mako") 

# plot with uncertainty
lFY <- melt(outp$NFY)
lB <- melt(outp$NB)
lF <- melt(outp$NF)
colnames(lFY) <- colnames(lB) <- colnames(lF) <- c("Time", "Site", "Iteration", "Abundance")
lFY$Year <- lFY$Time + 2010
lB$Year <- lB$Time + 2010
lF$Year <- lF$Time + 2010
lFY$Stage <- "First-year"
lB$Stage <- "Breeder"
lF$Stage <- "Nonbreeder"
lALL <- rbind(lFY, lB, lF)
lALL$Site <- ifelse(lALL$Site==1, "Los Haitises", "Punta Cana")

mdFY <-  apply(outp$NFY[1:13,,], c(1,2), median) 
mdB <-  apply(outp$NB[1:13,,], c(1,2), median) 
mdF <-  apply(outp$NF[1:13,,], c(1,2), median) 
hdiFY <-  apply(outp$NFY[1:13,,], c(1,2), HDInterval::hdi) 
hdiB <-  apply(outp$NB[1:13,,], c(1,2), HDInterval::hdi) 
hdiF <-  apply(outp$NF[1:13,,], c(1,2), HDInterval::hdi) 
medFY <- melt(mdFY)
medB <- melt(mdB)
medF <- melt(mdF)
hdisFY <- melt(hdiFY)
hdisB <- melt(hdiB)
hdisF <- melt(hdiF)
dfstages <- rbind(medFY, medB, medF)
dfhdis <- rbind(hdisFY, hdisB, hdisF)
dfstages$LHDI <- dfhdis[dfhdis$Var1=="lower",]
dfstages$UHDI <- dfhdis[dfhdis$Var1=="upper",]
cols <- mako(3)

med_hdis <- function(x, cm){
  df <- data.frame(y=median(x),
                   ymin=HDInterval::hdi(x, credMass=cm)[1],
                   ymax=HDInterval::hdi(x, credMass=cm)[2] )
}

p5 <- ggplot() +  theme_minimal(base_size=14) +
  theme(legend.position="none") +
  geom_line(data=lALL, aes(x=Year, y=Abundance, 
                           group=Iteration, color=Stage,
                           alpha=Stage), 
            linewidth=0.1) +
  stat_summary(data=lALL, aes(x=Year, y=Abundance),
               geom="errorbar", width=0,
               fun.data =med_hdis, fun.args=list(cm=0.95) ) +
  stat_summary(data=lALL, aes(x=Year, y=Abundance),
               geom="errorbar", width=0, linewidth=1,
               fun.data =med_hdis, fun.args=list(cm=0.80) ) +
  stat_summary(data=lALL, aes(x=Year, y=Abundance),
               geom="line", linewidth=1,
               fun.data =function(x){data.frame(y=median(x))} ) +
  scale_color_manual( values = c("First-year"=cols[2], "Breeder"=cols[1], "Nonbreeder"=cols[3])  ) +
  scale_alpha_manual( values = c("First-year"=0.015, "Breeder"=0.0075, "Nonbreeder"=0.15) ) +
  facet_wrap(Stage~Site, scales="free_y", 
             shrink=TRUE, ncol=2) +
  xlab("Year") + ylab("Number of females") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )

p35 <- ggarrange(p3, p5, ncol=1, nrow=2)
p35
# ggsave(filename="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Stage structure_abundance.tiff",
#        plot=p35,
#        device="tiff",
#        width=6, height=8, dpi=300)

#### ---- survival -----
# plot survival, recruitment, and detection
mus.mat <- outp$mus
dimnames(mus.mat)[[1]] <- c("FY", 
                            "A\nSurvival", "B", 
                            "FY:B",
                            "A:B\nRecruitment", 
                            "B:A",
                            "A", "B\nDetection")
dimnames(mus.mat)[[2]] <- c("Los Haitises", "Punta Cana")
lmus <- melt(mus.mat)
colnames(lmus) <- c("Parameter", "Site", "Iter", "value")

p6 <- ggplot(data= lmus, aes(x=value, y=Parameter )) +
  stat_pointinterval(point_interval=median_hdci, .width = c(.95, .80)) +
  coord_flip() +
  facet_wrap("Site") +
  xlab("Probability") + ylab("Parameter") +
  theme_bw(base_size=14) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) 

p6
# ggsave(filename="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Survival Recruitment Detection.tiff",
#        plot=p6,
#        device="tiff",
#        width=8, height=4, dpi=300)

# print estimates
mds <- tapply(lmus$value, list(lmus$Parameter, lmus$Site), median) |> round(2)
hdi_LH <- tapply(lmus$value, list(lmus$Parameter, lmus$Site), HDInterval::hdi)[,1]
hdi_LH2 <- do.call(rbind, hdi_LH) |> round(2)                                                                  
hdi_PC <- tapply(lmus$value, list(lmus$Parameter, lmus$Site), HDInterval::hdi)[,2]
hdi_PC2 <- do.call(rbind, hdi_PC) |> round(2) 
df <- data.frame(param= c("phiFY", "phiA", "phiB", "psiFY:B", "psiA:B", "psiB:A", "pA", "pB"),
  md_LH=mds[,1], lhdi_LH=hdi_LH2[,1], uhdi_LH=hdi_LH2[,2], 
  md_PC=mds[,2], lhdi_PC=hdi_PC2[,1], uhdi_PC=hdi_PC2[,2])
df
# Plot predicted survival, recruitment, and detection
# with effects from translocation and hacking
# Birds were only hacked from LHNP to PC here
# so we only predict values for PC

betas.pd <- apply(outp$betas, 1, pd)
ind <- which (betas.pd>0.975) # which are significant
pred.mus <- array(NA, dim=c(dim(outp$mus)))
for (m in ind){
  for (tr in 1:2){
    pred.mus[m,tr,] <- plogis( outp$lmus[m,2,] + outp$betas[m,]*c(0,1)[tr] ) 
  }}
# treatment, mu, site, iter
mus.md <- apply(pred.mus, c(1,2), median)
mus.HDI80 <- apply(pred.mus, c(1,2), HDInterval::hdi, credMass=0.8)
mus.HDI95 <- apply(pred.mus, c(1,2), HDInterval::hdi)
dimnames(pred.mus) <- list(Params=c("First-year\nSurvival", 
                                    "Nonbreeder\nSurvival",
                                    "Breeder\nSurvival",
                                    "First-year\nRecruitment",
                                    "Nonbreeder\nRecruitment",
                                    "Breeder to Nonbreeder\nTransition",
                                   "Nonbreeder\nDetection",
                                   "Breeder\nDetection"),
                           Translocated=c("Not Translocated",
                                          "Translocated"),
                           Iter=1:dim(pred.mus)[[3]])
lpms <- melt(pred.mus)
colnames(lpms)[[4]] <- "Probability"
lpms <- lpms[!is.na(lpms$Probability),]

p.betas <- ggplot(lpms, aes(x=Translocated, y=Probability) ) +
  stat_pointinterval(point_interval=median_hdci, .width = c(.95, .80)) +
  facet_wrap("Params", nrow=1) +
  theme_bw(base_size=14) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  xlab("") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), 
                     minor_breaks = seq(0, 1, by = 0.1)) +
  scale_x_discrete(labels=c("Not Translocated"="Not\nTranslocated",
                            "Translocated"="Trans-\nlocated"))

p.betas
# ggsave(filename="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Translocation effects_ggplot.tiff",
#        plot=p.betas,
#        device="tiff",
#        width=8, height=4, dpi=300)

#### ---- productivity -----
# Plot fecundity in response to nest treatments
f <- exp(outp$lmu.prod)
f.md <- apply(f, 1, median)
f.HDI80 <- apply(f, 1, HDInterval::hdi, credMass=0.8)
f.HDI95 <- apply(f, 1, HDInterval::hdi)
rownames(f) <- c("Los Haitises", "Punta Cana")
lf <- melt(f)
colnames(lf) <- c("Site", "Iter", "Productivity")
lf$Treatment <- "Untreated"

# Calculate treatment effects on fecundity
f.pred <- array(NA, dim=dim(outp$lmu.prod))
for (s in 1:2){
  f.pred[s,] <- exp(outp$lmu.prod[s,] + outp$gamma[,1])
} # s
f2.md <- apply(f.pred, 1, median)
f2.HDI80 <- apply(f.pred, 1, HDInterval::hdi, credMass=0.8)
f2.HDI95 <- apply(f.pred, 1, HDInterval::hdi)
rownames(f.pred) <- c("Los Haitises", "Punta Cana")
lf.pred <- melt(f.pred)
colnames(lf.pred) <- c("Site", "Iter", "Productivity")
lf.pred$Treatment <- "Treated"
lprod <- rbind(lf,lf.pred)  
lprod$Treatment <- factor(lprod$Treatment, levels=c('Untreated', 'Treated'))

p.prod <- ggplot(lprod, aes(x=Treatment, y=Productivity) ) +
            stat_pointinterval(point_interval=median_hdci, .width = c(.95, .80), expand=T) +
            facet_wrap("Site") +
            theme_bw(base_size=14) +
            theme(panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank()) +
            xlab("")+ 
            scale_y_continuous(breaks = seq(0, 1.5, by = 0.4), 
                               minor_breaks = seq(0, 1.5, by = 0.1))

p.prod
# ggsave(filename="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Fecundity_treatment_ggplot.tiff",
#        plot=p.prod,
#        device="tiff",
#        width=4, height=3, dpi=300)

# Calculate probability of direction
f.diff <- f[1,]-f[2,]
pd(f.diff) |> round(2)
f.md |> round(2)
f.HDI95 |> round(2)

f.diff2 <- f.pred[1,]-f.pred[2,]
pd(f.diff2) |> round(2)
f2.md  |> round(2)
f2.HDI95 |> round(2)

median(f.pred[1,]/f[1,])
median(f.pred[2,]/f[2,])

#### ---- cors_demo -----
#*******************
# Correlations among demographic rates and detection
#*******************
pdR <- apply(outp$R, c(1,2), pd) |> round(2)
max(pdR[pdR<1])
apply(outp$R, c(1,2), median) |> round(2)
apply(outp$R, c(1,2), HDInterval::hdi)[1,,] |> round(2)
apply(outp$R, c(1,2), HDInterval::hdi)[2,,] |> round(2)

plt  <- function(object, params,...) {
  MCMCplot(object=out, 
           params=params, 
           guide_axis=TRUE, 
           HPD=TRUE, ci=c(80, 95), horiz=FALSE, 
           #ylim=c(-10,10),
           ...)
}

ind <- 1
Rs <- c()
for (i in 1:(nrow(outp$R)-1)){
  for (j in (i+1):nrow(outp$R)){
    Rs[ind] <- paste0("R[",i,", ", j, "]")
    ind <- ind+1
  }}
par(mfrow=c(2,1))
plt(object=out, params=Rs[1:18], exact=TRUE, ISB=FALSE, 
    main="Correlations btw demographic rates\n over time and sites",
    xlab = "Rhos", guide_lines=TRUE)
plt(object=out, params=Rs[19:36], exact=TRUE, ISB=FALSE, 
    main="Correlations btw demographic rates\n over time and sites, continued ...",
    xlab = "Rhos", guide_lines=TRUE)

#### ---- popgr -----
#*******************
# Population growth rates over time
#*******************
lam.m <- apply(outp$lambda[1:12,,], c(1,2), median)
lam.hdi <- apply(outp$lambda[1:12,,], c(1,2), HDInterval::hdi)
par(mfrow=c(1,2))
plot(2012:2023, lam.m[,1], type="b", pch=1, 
     ylab="Population growth rate", xlab="Year", 
     ylim=c(min(lam.hdi[,,1]), max(lam.hdi[,,1])),
     main="Los Haitises")
abline(h=1, lty=2)
segments(x0=2012:2023, x1=2012:2023, 
         y0 = lam.hdi[1,,1], y1= lam.hdi[2,,1])

plot(2012:2023, lam.m[,2], type="b", pch=1, lty=1,
     ylab="Population growth rate", xlab="Year",
     ylim=c(min(lam.hdi[,,2]), max(lam.hdi[,,2])),
     main="Punta Cana")
abline(h=1, lty=2)
segments(x0=2012:2023, x1=2012:2023, 
         y0 = lam.hdi[1,,2], y1= lam.hdi[2,,2])

#### ---- cors1 -----
#*******************
# Correlations between demographics and population growth rates
#*******************
# create a function to plot correlations
# between demographics and population growth rates
plot.cor <- function (lam.post, x.post, x.lab, ind.x=1:12, ind.y=1:12, site=1, yaxt="b"){
  # calculate the correlation coefficient 
  # over each iteration to propagate error
  lambda <- lam.post[ind.y, site,]
  x <- x.post[ind.x, site, ]
  df <- data.frame(ats1=c(0.8, 0.9, 1.0, 1.1, 1.2, 1.3),
                   ats2=c(1.0, 1.5, 2.0, 2.5, 3.0, 3.5),
                   labs1=c("0.8", NA, "1.0", NA, "1.2", NA),
                   labs2=c("1.0", NA, "2.0", NA, "3.0", NA)
  )
  if(site==1){ df <- df[, c(1,3)]} else{
    df <- df[, c(2,4)]
  }
  lwd <- 1
  
  cor.post <- array(NA, dim=dim(lambda)[-1])
    for (i in 1:dim(lambda)[[2]]){
      cor.post[i] <- cor(lambda[,i], x[,i])
    }
  cor.df <- data.frame(median= median(cor.post) |> round(digits=2),
                       ldi=HDInterval::hdi(cor.post)[1] |> round(digits=2),
                       hdi=HDInterval::hdi(cor.post)[2] |> round(digits=2),
                       pd = pd(cor.post) |> round(digits=2))
  
  lam.m <- apply(lambda, 1, median)
  lam.hdi <- apply(lambda, 1, HDInterval::hdi)
  x.m <- apply(x, 1, median)
  x.hdi <- apply(x, 1, HDInterval::hdi)
    x.lims <- c(min(x.hdi[,]), max(x.hdi[,]))
    y.lims <- c(min(lam.hdi[,]), max(lam.hdi[,]))
    if(yaxt=="n"){
    plot(x.m, lam.m[ind.y], 
         xlim= x.lims,
         ylim= y.lims,
         type="n", 
         ylab="", xlab=x.lab, cex.lab=1.5, 
         main=c("", "")[site],
         yaxt=yaxt, xaxt="n")
      axis(2, at=df[,1], labels=c(rep(NA, 5)))
    } else {
      plot(x.m, lam.m[ind.y], 
           xlim= x.lims,
           ylim= y.lims,
           type="n", 
           ylab="", xlab=x.lab, cex.lab=1.5,  
           main=c("", "")[site], 
           yaxt="n", xaxt="n")
      axis(2, at=df[,1], labels=df[,2], las=1, cex.axis=1.5)
    }
    axis(1, cex.axis=1.5)
    points(x.m, lam.m, pch=1)
    segments(x0=x.hdi[1,], x1=x.hdi[2,], 
             y0 = lam.m, y1= lam.m, lwd=lwd)
    segments(x0=x.m, x1=x.m, 
             y0 = lam.hdi[1,], y1= lam.hdi[2,], lwd=lwd)
    xs <- c(x.lims[1], x.lims[1]+(x.lims[2]-x.lims[1]*1.5))
    ys <- c( (y.lims[2]-y.lims[1])+y.lims[1], (y.lims[2]-y.lims[1])+y.lims[1], 
             (y.lims[2]-y.lims[1])*0.6+y.lims[1], (y.lims[2]-y.lims[1])*0.6+y.lims[1])
    polygon(x=c(xs, rev(xs)), y=ys, col=alpha("white", 0.7), border=NA )
    text(x = x.lims[1], y = (y.lims[2]-y.lims[1])*0.9+y.lims[1], 
         paste("r = ", cor.df$median, " (",cor.df$ldi,", ", cor.df$hdi, ")", sep=""), 
         pos = 4, font = 3, cex = 1.5)
    text(x = x.lims[1], y = (y.lims[2]-y.lims[1])*0.7+y.lims[1], paste("P(r>0) = ", cor.df$pd, sep=""), 
         pos = 4, font = 3, cex = 1.5)
    
  }

# Plot correlations between population growth rates
# and demographics. 
# "r" is a correlation coefficient and represents
# the magnitude of the correlation
# P(r>0) is the probability of direction (similar to p-values)
# that is, the probability that an effect exists
# tiff(height=8, width=8, units="in", res=300,
#      filename= "C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//PopGrowthand Demographics1.tiff")
par(mfrow=c(4,3), mar=c(4,1,3,1), oma=c(1,5,1,1))
plot.cor(outp$lambda, outp$mn.prod, x.lab="Productivity", ind.x=2:13)
plot.cor(outp$lambda, outp$mn.phiFY, x.lab="First-year Survival", yaxt="n")
mtext("Los Haitises", side=3, line=1, cex=2)
plot.cor(outp$lambda, outp$mn.phiA, x.lab="Nonbreeder Survival", yaxt="n")
plot.cor(outp$lambda, outp$mn.phiB, x.lab="Breeder Survival")
plot.cor(outp$lambda, outp$mn.psiFYB, x.lab="First-year to Breeder", yaxt="n")
plot.cor(outp$lambda, outp$mn.psiAB, x.lab="Nonbreeder to Breeder", yaxt="n")

# tiff(height=4, width=8, units="in", res=300,
#      filename= "C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//PopGrowthand Demographics2.tiff")
# par(mfrow=c(2,3), mar=c(4,1,1,1), oma=c(0,5,3,0))
plot.cor(outp$lambda, outp$mn.prod, x.lab="Productivity", ind.x=2:13, site=2)
plot.cor(outp$lambda, outp$mn.phiFY, x.lab="First-year Survival", yaxt="n", site=2)
mtext("Punta Cana", side=3, line=1, cex=2)
plot.cor(outp$lambda, outp$mn.phiA, x.lab="Nonbreeder Survival", yaxt="n", site=2)
plot.cor(outp$lambda, outp$mn.phiB, x.lab="Breeder Survival", site=2)
plot.cor(outp$lambda, outp$mn.psiFYB, x.lab="First-year to Breeder", yaxt="n", site=2)
plot.cor(outp$lambda, outp$mn.psiAB, x.lab="Nonbreeder to Breeder", yaxt="n", site=2)

mtext("Population growth rate", side=2, outer=TRUE, line=2.5, cex=2)
# dev.off()

#### ---- cors2 -----
# Is there evidence of density dependence
par(mfrow=c(2,3), mar=c(4,1,1,1), oma=c(0,5,3,0))
plot.cor(outp$Ntot, outp$mn.prod, x.lab="Productivity", ind.x=2:13)
plot.cor(outp$Ntot, outp$mn.phiFY, x.lab="First-year Survival", yaxt="n")
plot.cor(outp$Ntot, outp$mn.phiA, x.lab="Nonbreeder Survival", yaxt="n")
plot.cor(outp$Ntot, outp$mn.phiB, x.lab="Breeder Survival")
plot.cor(outp$Ntot, outp$mn.psiFYB, x.lab="First-year to Breeder", yaxt="n")
plot.cor(outp$Ntot, outp$mn.psiAB, x.lab="Nonbreeder to Breeder", yaxt="n")
mtext("Abundance", side=2, outer=TRUE, line=2.5, cex=2)
mtext("Los Haitises", side=3, outer=TRUE, cex=2)

par(mfrow=c(2,3), mar=c(4,1,1,1), oma=c(0,5,3,0))
plot.cor(outp$Ntot, outp$mn.prod, x.lab="Productivity", ind.x=2:13, site=2)
plot.cor(outp$Ntot, outp$mn.phiFY, x.lab="First-year Survival", yaxt="n", site=2)
plot.cor(outp$Ntot, outp$mn.phiA, x.lab="Nonbreeder Survival", yaxt="n", site=2)
plot.cor(outp$Ntot, outp$mn.phiB, x.lab="Breeder Survival", site=2)
plot.cor(outp$Ntot, outp$mn.psiFYB, x.lab="First-year to Breeder", yaxt="n", site=2)
plot.cor(outp$Ntot, outp$mn.psiAB, x.lab="Nonbreeder to Breeder", yaxt="n", site=2)
mtext("Abundance", side=2, outer=TRUE, line=2.5, cex=2)
mtext("Punta Cana", side=3, outer=TRUE, cex=2)

#### ---- contributions -----
#********************
# Transient Life Table Response Experiment (tLTRE)
# Overall contributions
#********************
lambda <- expression( ( # numerator start
                      (nfy*phiFY*psiFYB*f*0.5 + # repro of first year breeders
                       nf*phiA*psiAB*f*0.5 + # repro of nonbreeders to breeders
                       nb*phiB*(1-psiBA)*f*0.5) + # repro of breeders remaining breeders
                       nfy*phiFY*(1-psiFYB)  +
                       nf*phiA*(1-psiAB) +
                       nb*phiB*psiBA +
                       nfy*phiFY*psiFYB +
                       nb*phiB*(1-psiBA) ) / # numerator end
                       (nfy+nf+nb) # denominator
                      ) # expression end

# Calculate proportional population sizes
nyr <- nrow(outp$NFY)-1
niter <- dim(outp$NFY)[[3]]
n1 <- outp$NFY[1:nyr, ,] / outp$Ntot[1:nyr, ,]
n2 <- outp$NF[1:nyr, ,] / outp$Ntot[1:nyr, ,]
n3 <- outp$NB[1:nyr, ,] / outp$Ntot[1:nyr, ,]

# Extract the mean demographic rates and population sizes and store them in a list
mu <- list(phiFY=outp$mus[1,,], phiA=outp$mus[2,,], phiB=outp$mus[3,,], 
           psiFYB=outp$mus[4,,], psiAB=outp$mus[5,,], psiBA=outp$mus[6,,], 
           f=exp(outp$lmu.prod),
           nfy=apply(n1, c(2,3), mean), nf=apply(n2, c(2,3), mean), nb=apply(n3, c(2,3), mean))

# Calculate sensitivities
nms <- c('phiFY', 'phiA', 'phiB',
          'psiFYB', 'psiAB', 'psiBA', 
          'f', 
          'nfy', 'nf', 'nb' )
sens.func <- function( name, expr=lambda, envir=mu) { eval(D(expr=expr, name=name), envir=envir) }
sensl <- lapply(nms, sens.func)
names(sensl) <- nms
sens <- abind(sensl, along=3)

# Define matrix to store results
cont <- array(NA, dim=c(2, niter, 10) )
dimnames(cont)[[3]] <- nms

# Calculate contributions for each demographic rate and stage-structured population size at each MCMC draw
for (i in 1:niter){
  for (s in 1:2){
  dp_stoch <- cbind(outp$mn.phiFY[1:12,s,i], outp$mn.phiF[1:12,s,i], outp$mn.phiB[1:12,s,i], 
                    outp$mn.psiFYB[1:12,s,i], outp$mn.psiAB[1:12,s,i], outp$mn.psiBA[1:12,s,i],
                    outp$mn.prod[2:13,s,i],
                    n1[1:12,s,i], n2[1:12,s,i], n3[1:12,s,i])
  # Derive process variance and covariance among demographic parameters using shrunk estimates of
  #     demographic rates and proportional pop. sizes
  dp_varcov <- var(dp_stoch)
  sensvec <- sens[s, i, ]
  # Calculate demographic contributions
  contmatrix <- dp_varcov * outer(sensvec, sensvec)
  cont[s, i,  ] <- rowSums(contmatrix)
}}

CRI <- apply(cont[,,], c(1,3), quantile, probs=c(0.025, 0.975))
names.arg <- c(expression(italic(phi)[italic(FY)]),
               expression(italic(phi)[italic(A)]),
               expression(italic(phi)[italic(B)]),
               expression(italic(psi)[italic(FY:B)]),
               expression(italic(psi)[italic(A:B)]),
               expression(italic(psi)[italic(B:A)]),
               expression(italic(f)),
               expression(italic(n)[italic(FY)]),
               expression(italic(n)[italic(A)]), 
               expression(italic(n)[italic(B)]))
# tiff(height=4, width=6, units="in", res=300,
#      filename= "C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//LTRE_overall.tiff")
par(mfrow=c(2,1), mar=c(3,5,1,1))
a <- barplot(colMeans(cont[1,,]), names.arg=names.arg, 
             ylab=expression( paste("Contribution to ", lambda) ), las=1,
             ylim=range(CRI[,1,]), col=rep("grey65", 10),
             border=rep("grey65", 10),
             main="Los Haitises", font.main=1, 
             yaxt="n")
axis(2, at=c(0, 0.01, 0.02, 0.03), las=1)
segments(a, CRI[1,1,], a, CRI[2,1,], lwd=1.5)

b <- barplot(colMeans(cont[2,,]), names.arg=names.arg, 
             ylab=expression( paste("Contribution to ", lambda) ), las=1,
             ylim=range(CRI[,2,]), col=rep("grey65", 10),
             border=rep("grey65", 10), 
             main="Punta Cana", font.main=1,
             yaxt="n")
axis(2, at=c(0, 0.01, 0.02, 0.03), las=1)
segments(a, CRI[1,2,], a, CRI[2,2,], lwd=1.5)
# dev.off()

#*******************
# Transient Life Table Response Experiment (tLTRE) 
# Inter-annual Contributions
#*******************
sens <- diff <- array(NA, dim=c(constl$nyr-1, 2, 10, dim(outp$mn.phiA)[[3]]),
              dimnames=list(Year=2012:2023, 
                            site=c('Los Haitises', 'Punta Cana'), 
                            Param=c('phiFY', 'phiA', 'phiB',
                                      'psiFYB', 'psiAB', 'psiBA', 
                                      'f', 
                                      'nfy', 'nf', 'nb' ),
                            Iter=1:dim(outp$mn.phiA)[[3]]) )

getDiff <- function(x){
            x[2:constl$nyr, , ] - x[1:(constl$nyr-1), , ]
}
diff[,,'phiFY',] <- getDiff(outp$mn.phiFY)
diff[,,'phiA',] <- getDiff(outp$mn.phiA)
diff[,,'phiB',] <- getDiff(outp$mn.phiB)
diff[,,'psiFYB',] <- getDiff(outp$mn.psiFYB)
diff[,,'psiAB',] <- getDiff(outp$mn.psiAB)
diff[,,'psiBA',] <- getDiff(outp$mn.psiBA)
diff[,,'f',] <- getDiff(outp$mn.prod)
diff[,,'nfy',] <- getDiff( outp$NFY/outp$Ntot )
diff[,,'nf',] <- getDiff( outp$NF/outp$Ntot )
diff[,,'nb',] <- getDiff( outp$NB/outp$Ntot )

getMn <- function(x) {
  ( x[2:constl$nyr, , ] + x[1:(constl$nyr-1), , ] ) / 2
} 
means <- list('phiFY' = getMn(outp$mn.phiFY), 
              'phiA' = getMn(outp$mn.phiA), 
              'phiB' = getMn(outp$mn.phiB),
              'psiFYB' = getMn(outp$mn.psiFYB), 
              'psiAB' = getMn(outp$mn.psiAB),
              'psiBA' = getMn(outp$mn.psiBA),
              'f' = getMn(outp$mn.prod), 
              'nfy' = getMn(outp$NFY), 
              'nf' = getMn(outp$NF), 
              'nb' =  getMn(outp$NB) )

sens[,,'phiFY',] <- eval( D(lambda, "phiFY" ), envir = means )
sens[,,'phiA',] <- eval( D(lambda, "phiA" ), envir = means )
sens[,,'phiB',] <- eval( D(lambda, "phiB" ), envir = means )
sens[,,'psiFYB',] <- eval( D(lambda, "phiFYB" ), envir = means )
sens[,,'psiAB',] <- eval( D(lambda, "psiAB" ), envir = means )
sens[,,'psiBA',] <- eval( D(lambda, "psiBA" ), envir = means )
sens[,,'f',] <- eval( D(lambda, "f" ), envir = means )
sens[,,'nfy',] <- eval( D(lambda, "nfy" ), envir = means )
sens[,,'nf',] <- eval( D(lambda, "nf" ), envir = means )
sens[,,'nb',] <- eval( D(lambda, "nb" ), envir = means )

conts <- diff*sens

# ~~~~ figure 9.7 ~~~~
# Create a plot for the sequential differences in finite growth rates
niter <- dim(outp$mn.phiA)[[3]]
nyrs <- constl$nyr-1
diff.lam <- outp$lambda[2:nyrs, , ] - outp$lambda[1:(nyrs-1), ,]
differences <- cbind(2012:2022, apply(diff.lam, c(1,2), mean))
differences <- rbind(differences, c(2023, NA, NA))

# Mean contributions
V <- apply(conts, 3:1, mean)
V1 <- V2 <- V
V1[V1<0] <- 0
V2[V2>0] <- 0

# Make figure
colors <- c(colorRampPalette(c("blue4", "lightblue"))(3),
             colorRampPalette(c("darkred", "yellow"))(3),
             colorRampPalette(c("darkviolet"))(1), 
             colorRampPalette(c("lightgreen", "darkgreen"))(3))

legend.text <- c(expression(italic(phi)[italic(FY)]),
                 expression(italic(phi)[italic(A)]),
                 expression(italic(phi)[italic(B)]),
                 expression(italic(psi)[italic(FY:B)]),
                 expression(italic(psi)[italic(A:B)]),
                 expression(italic(psi)[italic(B:A)]),
                 expression(italic(f)),
                 expression(italic(n)[italic(FY)]),
                 expression(italic(n)[italic(A)]), 
                 expression(italic(n)[italic(B)]))

# tiff(height=6.5, width=8, units="in", res=300,
#      filename= "C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//LTRE.tiff")
par(mfrow=c(2,2), mar=c(1, 2, 0, 1), oma=c(7,7,1,1))
txt.size <- 1.5
yl <- list(bquote( paste("Difference in population")),
         bquote( paste("growth rate (", Delta, lambda, ")" )) ) 
barplot(differences[,2], ylim=c(-0.16, 0.16),
        #ylab=expression( paste("Difference in population", "\ngrowth rate (", Delta, lambda, ")" )), 
        ylab= "",
        xpd=NA,
        axes=FALSE, border=NA,
        cex.lab=txt.size, main="Los Haitises", cex.main=txt.size, 
        font.main=1 )
mtext(do.call(expression,yl), side=2, cex=txt.size,
      line=c(5,3.5) )
axis(2, las=1)
abline(h=0)

barplot(differences[,3], ylim=c(-1.5, 1.5),
        ylab="", 
        xpd=NA,
        axes=FALSE, border=NA, main="Punta Cana",
        cex.main=txt.size, 
        font.main=1)
axis(2, las=1)
abline(h=0)

barplot(V1[,1,], ylim=c(-0.19, 0.19), 
        col=colors, ylab="",
        xlab="", 
        xpd=NA,
        axes=FALSE,
        border=NA,
        cex.lab=txt.size)
mtext(expression( paste("Contribution to ", Delta, lambda) ), side=2, cex=txt.size,
      line=4.25)
barplot(V2[,1,], add=TRUE, col=colors, axes=FALSE, border=NA)
axis(2, las=1)
abline(h=0)

barplot(V1[,2,], ylim=c(-0.12, 0.12), 
        col=colors, ylab="",
        xlab="", axes=FALSE,  
        names.arg=2012:2023,
        border=NA)
barplot(V2[,2,], add=TRUE, col=colors, axes=FALSE, border=NA)
axis(2, las=1)
abline(h=0)
mtext(side=1 , outer=TRUE, "Beginning year", line=2.5, cex=txt.size)
reset <- function() {
  par(mfrow=c(1, 1), oma=rep(0, 4), mar=rep(0, 4), new=TRUE)
  plot(0:1, 0:1, type="n", xlab="", ylab="", axes=FALSE)
}
reset()
legend(legend = legend.text, fill=colors,
       bty="n", x="bottom", border=NA, 
       xpd=NA, x.intersp=0.005, text.width=NA,
       cex=1.5, horiz=TRUE)
# dev.off()

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
  # facet_grid(rows=  vars(Site) , cols= vars(Mortality.reduction), 
  #            scales='free_y') +
  #facet_wrap(~  Site + Mortality.reduction, scales='free_y') +
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

p7 <- ggplot(pe.df2, aes(x = as.factor(Number.translocated), 
                  y = as.factor(Number.of.nests.treated))) +
  geom_tile( aes(fill = Extinct ) ) + 
  geom_text(aes(label=round(Extinct,2)), color="gray70") +
  scale_fill_viridis_c(option="plasma", direction=1,
                       name = "Extinction\nprobability after\n50 years") +
  facet_wrap(~  Site + factor(Mortality.reduction, 
                              levels=c("100% Mortality", "80% Mortality", "60% Mortality"))) +
  scale_y_discrete(breaks=c(0,15,30,45, 100), 
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

# calculate future medians
lnt2 <- apply(outp$Ntot[,13:63, , ], c(1, 2, 3), median) |> melt()
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

#******************************
#* Population change
#* between 2023 and 2073
#******************************
ab2023 <- nt.df2[nt.df2$Year=='2023', ]
ab2023 <- ab2023[order(ab2023$Scenarios, ab2023$Site2),]
ab2073 <- nt.df2[nt.df2$Year=='2073', ]
ab2073 <- ab2073[order(ab2073$Scenarios, ab2073$Time, ab2073$Site2),]

ab2073$lambda <- ab2073$Abund / ab2023$Abund
ab2073$pc <- (ab2073$lambda-1)*100

# scale color so that percent change=0 is the middle color (gray),
# -100 is the extreme lower bound and 811 is the extreme upper bound
for (i in 1:length(ab2073$pc)){
  if(ab2073$pc[i]<=0){  
  ab2073$pc.sc[i] <- (ab2073$pc[i]-0)/sd(ab2073$pc[ab2073$pc<=0])}
  else{ ab2073$pc.sc[i] <- (ab2073$pc[i]-0)/sd(ab2073$pc[ab2073$pc>0]) }
}

labs <- c(-100, 0, 400, 800)
labs.sc <- c( (-100-0)/sd(ab2073$pc[ab2073$pc<=0]),
              (0-0)/sd(ab2073$pc[ab2073$pc<=0]),
              (400-0)/sd(ab2073$pc[ab2073$pc>0]),
              (800-0)/sd(ab2073$pc[ab2073$pc>0]) )
# plot
p11 <- ggplot(ab2073, aes(x = as.factor(Number.translocated), 
                         y = as.factor(Number.of.nests.treated))) +
  geom_tile( aes(fill = pc.sc ) ) + 
  geom_text(aes(label=round(pc,0)), color="white") +
  scale_fill_viridis_c(option="cividis", direction=-1,
                       name = "Percent change\nafter 50 years",
                       breaks = labs.sc, labels = labs,
                       begin=0, end=1) +
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
