## ---- postprocess --------
library ('MCMCvis')
library ('coda')
library ('ggplot2')
library('reshape2')
library ("tidybayes")
library ('bayestestR')
library ("ggpubr")
library("viridis")
library ("HDInterval")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_longrun.rdata")
load("data/data.rdata")
out <- lapply(post, as.mcmc) 
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
(outp$lmu.prod[1,]-outp$lmu.prod[2,]) |> pd()

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
  stat_halfeye(point_interval=median_hdi, .width = c(.95, .80)) +
  coord_flip() +
  facet_wrap("Site") +
  xlab("Probability") + ylab("Parameter") +
  theme_bw(base_size=14) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) 

p6
# ggsave(filename="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Survival Recruitment Detection.tiff",
#        plot=p5,
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
# tiff(height=4, width=6, units="in", res=300,
#      filename= "C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Translocation effects.tiff")
par(mar=c(3,4,4,1), mfrow=c(1,1))
for (tr in 1:2){
  if (tr==1){
    plot(1:4, c(0,0,1,1), 
         type="n",
         pch=c(16, 17)[tr], cex=2,
         ylim=c(0,1),
         xlim=c(0.5, 4.5), 
         cex.axis=1, cex.lab=1.25,
         ylab="Probability", xlab="",
         main="",
         xaxt="n", yaxt="n")
  abline(h= seq(0,1,by=0.1), col="gray90")
  abline(v= seq(0.5,6,by=1), col="gray90")
  points(1:4+c(-0.1,0.1)[tr], mus.md[ind,tr], 
         pch=c(16, 17)[tr], cex=1.5)
  }else{
    points(1:4+c(-0.1,0.1)[tr], mus.md[ind,tr], 
           pch=c(16, 15)[tr], cex=1.5)
  }}

axis(1, at=1:4, cex.axis=0.85,
     labels=c("First-year\nSurvival", 
              "Breeder\nSurvival", 
              "Nonbreeder\n to Breeder\nRecruitment", 
              "Nonbreeder\nDetection"), las=1, padj=0.5)
axis(2, at=seq(0,1, by=0.1), cex.axis=1,
     labels=c(0, NA,NA,NA,NA,0.5, NA,NA,NA,NA, 1))
for (m in 1:length(ind)){
  for (tr in 1:2){
    lines(rep(c(1:4)[m] + c(-0.1,0.1)[tr],2), mus.HDI95[1:2,ind[m],tr], lwd=2)
    lines(rep(c(1:4)[m] + c(-0.1,0.1)[tr],2), mus.HDI80[1:2,ind[m],tr], lwd=4)
  }}
legend(x=2.5,y=1.25, pch=c(16,15), pt.cex=1.5, cex=1, xpd=NA, horiz=T, 
       legend=c("Not translocated", "Translocated" ),
       xjust=0.5)
#dev.off()

#### ---- productivity -----
# Plot fecundity in response to nest treatments
f <- exp(outp$lmu.prod)
f.md <- apply(f, 1, median)
f.HDI80 <- apply(f, 1, HDInterval::hdi, credMass=0.8)
f.HDI95 <- apply(f, 1, HDInterval::hdi)

# Calculate treatment effects on fecundity
f.pred <- array(NA, dim=dim(outp$lmu.prod))
for (s in 1:2){
  f.pred[s,] <- exp(outp$lmu.prod[s,] + outp$gamma[,1])
} # s
f2.md <- apply(f.pred, 1, median)
f2.HDI80 <- apply(f.pred, 1, HDInterval::hdi, credMass=0.8)
f2.HDI95 <- apply(f.pred, 1, HDInterval::hdi)

# tiff(height=4, width=6, units="in", res=300,
#      filename= "C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Fecundity_treatment.tiff")
par(mar=c(3,5,3,1), mfrow=c(1,1))
plot(c(1, 2)-0.1, f.md, 
     ylim=c(0, 1.55),
     xlim=c(0.5, 2.5), 
     ylab="Productivity", xlab="",
     main="",
     cex.axis=1.5, cex.lab=1.5,
     xaxt="n", yaxt="n", type="n")
points(c(1, 2)-0.1, f.md, pch=16, cex=1.5)
abline(h=seq(0, 1.6, by=0.1), col="gray90")
abline(v=1.5, col="black")
axis(1, at=c(1,2), cex.axis=1,
     labels=c("Los Haitises\nNational Park","Punta Cana"),
     padj=0.5)
axis(2, at=seq(0, 1.6, by=0.1), 
     cex.axis=1,
     labels=c("0", NA, NA, NA, NA, "0.5", NA, NA, NA, NA, "1.0", NA, NA, NA, NA, "1.5", NA))
lines(c(1,1)-0.1, f.HDI95[,1], lwd=2)
lines(c(2,2)-0.1, f.HDI95[,2], lwd=2)
lines(c(1,1)-0.1, f.HDI80[,1], lwd=4)
lines(c(2,2)-0.1, f.HDI80[,2], lwd=4)

points(c(1.1, 2.1), f2.md, 
       pch=18, cex=2.5)

lines(c(1,1) +0.1, f2.HDI95[,1], lwd=3)
lines(c(2,2) +0.1, f2.HDI95[,2], lwd=3)
lines(c(1,1) +0.1, f2.HDI80[,1], lwd=6)
lines(c(2,2) +0.1, f2.HDI80[,2], lwd=6)
legend(x=1.5,y=1.9, pch=c(16,18), pt.cex=c(1.5,2.5), cex=1,
       legend=c("Untreated", "Treated" ), 
       horiz=TRUE, xjust=0.5, xpd=NA)
#dev.off()

# Calculate probability of direction
f.diff <- f[1,]-f[2,]
pd(f.diff)
f.md 
f.HDI95 

f.diff2 <- f.pred[1,]-f.pred[2,]
pd(f.diff2)
f2.md 
f2.HDI95

median(f.pred[1,]/f[1,])
median(f.pred[2,]/f[2,])

#*******************
# Correlations between demographic rates
#*******************
pdR2 <- apply(outp$R2, c(1,2), pd) |> round(2)
max(pdR2[pdR2<1])
apply(outp$R2, c(1,2), median) |> round(2)
apply(outp$R2, c(1,2), HDInterval::hdi)[1,,] |> round(2)
apply(outp$R2, c(1,2), HDInterval::hdi)[2,,] |> round(2)

#### ---- cors -----
#*******************
# Correlations between demographics and population growth rates
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

# Plor correlations between population grrowth rates
# and demographics. 
# "r" is a correlation coefficient and represents
# the magnitude of the correlation
# P(r>0) is the probability of direction (similar to p-values)
# that is, the probability that an effect exists
# tiff(height=8, width=8, units="in", res=300,
#      filename= "C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//PopGrowthand Demographics1.tiff")
# par(mfrow=c(4,3), mar=c(4,1,1,1), oma=c(0,5,3,0))
par(mfrow=c(2,3), mar=c(4,1,1,1), oma=c(0,5,3,0))
plot.cor(outp$lambda, outp$mn.prod, x.lab="Fecundity", ind.x=2:13)
plot.cor(outp$lambda, outp$mn.phiFY, x.lab="First-year Survival", yaxt="n")
plot.cor(outp$lambda, outp$mn.phiA, x.lab="Nonbreeder Survival", yaxt="n")
plot.cor(outp$lambda, outp$mn.phiB, x.lab="Breeder Survival")
plot.cor(outp$lambda, outp$mn.psiFYB, x.lab="First-year to Breeder", yaxt="n")
plot.cor(outp$lambda, outp$mn.psiAB, x.lab="Nonbreeder to Breeder", yaxt="n")
mtext("Population growth rate", side=2, outer=TRUE, line=2.5, cex=2)
mtext("Los Haitises", side=3, outer=TRUE, cex=2)
#dev.off()

# tiff(height=4, width=8, units="in", res=300,
#      filename= "C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//PopGrowthand Demographics2.tiff")
par(mfrow=c(2,3), mar=c(4,1,1,1), oma=c(0,5,3,0))
plot.cor(outp$lambda, outp$mn.prod, x.lab="Fecundity", ind.x=2:13, site=2)
plot.cor(outp$lambda, outp$mn.phiFY, x.lab="First-year Survival", yaxt="n", site=2)
plot.cor(outp$lambda, outp$mn.phiA, x.lab="Nonbreeder Survival", yaxt="n", site=2)
plot.cor(outp$lambda, outp$mn.phiB, x.lab="Breeder Survival", site=2)
plot.cor(outp$lambda, outp$mn.psiFYB, x.lab="First-year to Breeder", yaxt="n", site=2)
plot.cor(outp$lambda, outp$mn.psiAB, x.lab="Nonbreeder to Breeder", yaxt="n", site=2)
mtext("Abundance", side=2, outer=TRUE, line=2.5, cex=2)
mtext("Punta Cana", side=3, outer=TRUE, cex=2)
#dev.off()

# Is there evidence of density dependence
par(mfrow=c(2,3), mar=c(4,1,1,1), oma=c(0,5,3,0))
plot.cor(outp$Ntot, outp$mn.prod, x.lab="Fecundity", ind.x=2:13)
plot.cor(outp$Ntot, outp$mn.phiFY, x.lab="First-year Survival", yaxt="n")
plot.cor(outp$Ntot, outp$mn.phiA, x.lab="Nonbreeder Survival", yaxt="n")
plot.cor(outp$Ntot, outp$mn.phiB, x.lab="Breeder Survival")
plot.cor(outp$Ntot, outp$mn.psiFYB, x.lab="First-year to Breeder", yaxt="n")
plot.cor(outp$Ntot, outp$mn.psiAB, x.lab="Nonbreeder to Breeder", yaxt="n")
mtext("Abundance", side=2, outer=TRUE, line=2.5, cex=2)
mtext("Los Haitises", side=3, outer=TRUE, cex=2)
#dev.off()

# tiff(height=4, width=8, units="in", res=300,
#      filename= "C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//PopGrowthand Demographics2.tiff")
par(mfrow=c(2,3), mar=c(4,1,1,1), oma=c(0,5,3,0))
plot.cor(outp$Ntot, outp$mn.prod, x.lab="Fecundity", ind.x=2:13, site=2)
plot.cor(outp$Ntot, outp$mn.phiFY, x.lab="First-year Survival", yaxt="n", site=2)
plot.cor(outp$Ntot, outp$mn.phiA, x.lab="Nonbreeder Survival", yaxt="n", site=2)
plot.cor(outp$Ntot, outp$mn.phiB, x.lab="Breeder Survival", site=2)
plot.cor(outp$Ntot, outp$mn.psiFYB, x.lab="First-year to Breeder", yaxt="n", site=2)
plot.cor(outp$Ntot, outp$mn.psiAB, x.lab="Nonbreeder to Breeder", yaxt="n", site=2)
mtext("Population growth rate", side=2, outer=TRUE, line=2.5, cex=2)
mtext("Punta Cana", side=3, outer=TRUE, cex=2)
#dev

#### ---- pva -----
#*******************
#* Plot PVA results
#*******************
library ('viridis')
library ("tidybayes")
library ('coda')
library ('MCMCvis')
library ('reshape2')
library('ggplot2')
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\pva_longrun.rdata")
out <- lapply(post, as.mcmc) 
outp <- MCMCpstr(out, type="chains")
scen <- data.frame('Scenarios' = 1:36,
                   'Number translocated' = rep(c(0,5,10), each=12), 
                   'Number of nests treated' = rep( rep(c(0, 15, 30, 45), each=3), 3),
                   'Mortality reduction' = rep(c('100% Mortality', '80% Mortality', '60% Mortality'), 12) )

#*********************
#* Line plot of extinction probability over time
#*********************
pe <- apply(outp$extinct, c(1,2,3), mean)
lpe <- melt(pe)
colnames(lpe) <- c('Scenarios', 'Time', 'Site2', 'Extinct')
pe.df <- merge(lpe, scen, by='Scenarios', all.x=TRUE)
pe.df$Year <- pe.df$Time+2023
pe.df$Site <- ifelse(pe.df$Site2==1, "Los Haitises", "Punta Cana")
pe.df[(pe.df$Site2==2 & pe.df$Number.translocated>0),]$Extinct <- NA
pe.df$Mortality.reduction <- factor(pe.df$Mortality.reduction, 
                                    levels= c('100% Mortality', '80% Mortality', '60% Mortality'))

p6 <- ggplot(pe.df, aes(x=Year, y=Extinct, group=Scenarios, 
                  color = factor(Number.translocated), 
                  linetype = factor(Number.of.nests.treated))) +
  geom_line( linewidth = 1) +
  facet_grid(rows=  vars(Site), cols= vars(Mortality.reduction), 
             scales='free_y') +
  #facet_wrap(~  Site + Mortality.reduction) +
  scale_color_viridis_d(option="viridis", direction=1,
                        name="Number translocated") +
  scale_linetype(name="Number of nests treated") +
  scale_x_continuous( name="Future year",
                      breaks=c(2020, 2030, 2040, 2050, 2060, 2070),
                     labels=c('', '2030', '', '2050', '', '2070'), 
                     minor_breaks=NULL) +
  theme_minimal() + 
  ylab("Extinction probability")
p6
# ggsave(filename="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Extinction_lines.tiff",
#        plot=p6,
#        device="tiff",
#        width=6, height=4, dpi=300)

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
  scale_y_discrete(breaks=c(0,15,30,45), 
                   name="Number of nests treated") +
  scale_x_discrete(name="Number translocated") +
  theme_minimal( ) + theme(strip.background = element_rect(size = 0.5))
p7  
# ggsave(filename="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Extinction_50yrHeatmap.tiff",
#        plot=p7,
#        device="tiff",
#        width=6, height=4, dpi=300)
pe.df3 <- pe.df2[!is.na(pe.df2$Extinct),]
pe.df3 <- pe.df3[rev(order(pe.df3$Site2, pe.df3$Extinct)),]
pe.df3

#***************************
# Plot total Abundance
lnt1 <-  outp$Ntot[,1:13, , ] |> melt()
colnames(lnt1) <- c("Scenarios", "Time", "Site2", "Iter", "Abund")
lnt1$Year <- lnt1$Time + 2010
lnt1$Site <- factor(ifelse(lnt1$Site2==1, 'Los Haitises', 'Punta Cana'))
nt.df <- merge(lnt1, scen, by='Scenarios', all.x=TRUE)
nt.df$Mortality.reduction <- factor(nt.df$Mortality.reduction,
                                    levels=c('100% Mortality','80% Mortality', '60% Mortality' ))
nt.df[(nt.df$Site2==2 & nt.df$Number.translocated>0),]$Abund <- NA
nt.df <- nt.df[nt.df$Scenarios %in% c(1,2,3), ]

# calculate future means
lnt2 <- apply(outp$Ntot[,13:63, , ], c(1, 2, 3), median) |> melt()
colnames(lnt2) <- c("Scenarios", "Time", "Site2", "Abund")
lnt2$Year <- lnt2$Time + 2022
lnt2$Site <- factor(ifelse(lnt2$Site2==1, 'Los Haitises', 'Punta Cana'))
nt.df2 <- merge(lnt2, scen, by='Scenarios', all.x=TRUE)
nt.df2$Mortality.reduction <- factor(nt.df2$Mortality.reduction,
                                    levels=c('100% Mortality','80% Mortality', '60% Mortality' ))
nt.df2[(nt.df2$Site2==2 & nt.df2$Number.translocated>0),]$Abund <- NA

p8 <- ggplot() +
      stat_lineribbon(data= nt.df,  
                      aes(x=Year, y=Abund, group=Scenarios),
                      fill='gray40',
                      linewidth=0.5,
                      point_interval = 'median_hdci', .width = 0.95,
                       na.rm=TRUE) +
      geom_line(data= nt.df2,
                aes(x=Year, y=Abund, group=Scenarios,
                    color = factor(Number.translocated),
                    linetype = factor(Number.of.nests.treated) ),
                linewidth=0.5) +
      facet_grid(rows = vars(Site), cols = vars(Mortality.reduction),
                 scales ='free_y') +
      ylab('Number of females') +
      scale_color_viridis_d(option="viridis", direction=1,
                            name="Number translocated") +
      scale_linetype(name="Number of nests treated") +
      theme_minimal()
  
p8
# ggsave(filename="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Abundance_projections.tiff",
#        plot=p9,
#        device="tiff",
#        width=8, height=4, dpi=300)

# Time to extinction
# conditional on extinction threshold (D=X)
#### ---- pva2 -----
D <- 1
abund <- outp$Ntot[,14:(13+50),,]
dimnames(abund) <- list('Scenarios'=1:36, 
                        'Time'= 2024:(2024+49), 
                        'Site'=c("Los Haitises", "Punta Cana"), 
                        'Iter'=1:10000)
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
                       levels=paste0('Scenario ', 1:36)) 
teLH <- te[te$Site=="Los Haitises",]
tePC <- te[te$Site=="Punta Cana",]

p9 <- ggplot() +
      geom_col(data = teLH, aes(x= TE, y=Count)) + 
      facet_wrap(vars(Scenarios),
                 nrow=12, ncol=3, scales="free_y", shrink=FALSE)

p10 <- ggplot() +
  geom_col(data = tePC, aes(x= TE, y=Count)) + 
  facet_wrap(vars(Scenarios),
             nrow=12, ncol=3, scales="free_y", shrink=FALSE)
p9
p10