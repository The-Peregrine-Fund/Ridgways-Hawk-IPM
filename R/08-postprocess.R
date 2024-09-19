library ('MCMCvis')
library ('coda')
library ('ggplot2')
library('reshape2')
library ("tidybayes")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\ipm_longrun.rdata")
load("C:\\Users\\rolek.brian\\OneDrive - The Peregrine Fund\\Documents\\Projects\\Ridgways IPM\\outputs\\pva.rdata")
load("data/data.rdata")
out <- lapply(post, as.mcmc) 
outp <- MCMCpstr(out, type="chains")

#********************
#*Calculate summary stats
#********************
# Population growth rates
mn.lambda <- apply(outp$lambda[1,1:12,,], c(2,3), mean)
(md.lambda <- apply(mn.lambda, 1, median))
(hdi95.lambda <- apply(mn.lambda, 1, HDInterval::hdi))
# Overall averages
mn.mus <- apply(outp$mus, c(1,3), mean)
(md.mus <- apply(mn.mus, 1, median))
(hdi.mus <- apply(mn.mus, 1, HDInterval::hdi))
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
#*******************
#* Plot IPM results
#*******************
#Plot population size
Ntot <- outp$Ntot[1,1:13,,] # thin by 10 to speed up plotting
lN <- melt(Ntot)
colnames(lN) <- c("Time", "Site.ind", "Iter", "Abundance")
lN$Year <- lN$Time +2010
lN$Site <- ifelse(lN$Site.ind==1, "Los Haitises", "Punta Cana")

med_hdis <- function(x, cm){
            df <- data.frame(y=median(x),
                             ymin=HDInterval::hdi(x, credMass=cm)[1],
                             ymax=HDInterval::hdi(x, credMass=cm)[2] )
}

counts.tot <- datl$countsAdults+datl$countsFY+constl$hacked.counts[,-3]
df.counts <- melt(counts.tot)
colnames(df.counts) <- c("Year", "Site_ab", "Count")
df.counts$Site <- ifelse(df.counts$Site_ab=="LHNP", "Los Haitises", "Punta Cana")

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

ggsave(filename="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Abundance and counts.tiff",
       plot=p1,
       device="tiff",
       width=6, height=4, dpi=300)

#### ----popstructure -----
# Plot life stage structure
mdFY <-  apply(outp$NFY[1,1:13,,], c(1,2), median) 
mdB <-  apply(outp$NB[1,1:13,,], c(1,2), median) 
mdF <-  apply(outp$NF[1,1:13,,], c(1,2), median) 
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
  ylab("Proportion of population") + 
  facet_wrap("Site") + 
  theme_minimal(base_size=14)

p4 <- ggplot(ldat, aes(fill=Stage, y=as.numeric(Number), x=year)) + 
  geom_bar(position="stack", stat="identity") +
  ylab("Median number of females") + xlab("Year") +
  facet_wrap("Site", scales = "free_y") + 
  theme_minimal(base_size=14 ) + theme(legend.position="top") + 
  scale_fill_viridis_d(option = "mako") 

ggsave(filename="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Stage structure.tiff",
       plot=p4,
       device="tiff",
       width=6, height=4, dpi=300)

# PLot predicted survival, recruitment, and detection
# with effects from translocation and hacking
# Birds were only hacked from LHNP to PC here
# so we only predict values for PC
# pred.mus <- array(NA, dim=dim(outp$mus))
# for (m in 1:8){
#   for (tr in 1:2){
#     pred.mus[m,tr,] <- plogis( outp$lmus[m,2,] + outp$betas[m,]*c(0,1)[tr] ) 
#   }}
# 
# dimnames(pred.mus) <- list(param=c("FY", "A", "B", "FY:B", "A:B", "B:A", "A", "B"),
#                            trans=c("Not translocated", "Translocated"), 
#                            iter=1:10000)
# lm <- melt(pred.mus)
# 
# p5 <- ggplot(data= lm, aes(x=value, y=param )) +
#   stat_halfeye(point_interval=median_hdi, .width = c(.95, .80)) +
#   coord_flip() +
#   facet_wrap("trans") +
#   xlab("Probability") + ylab("Parameter") +
#   theme_bw(base_size=14) +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank()) 

pred.mus <- array(NA, dim=c(dim(outp$mus)))
for (m in c(1,2,3,5,7)){
  for (tr in 1:2){
    pred.mus[m,tr,] <- plogis( outp$lmus[m,2,] + outp$betas[m,]*c(0,1)[tr] ) 
  }}
# treatment, mu, site, iter
mus.md <- apply(pred.mus, c(1,2), median)
mus.HDI80 <- apply(pred.mus, c(1,2), HDInterval::hdi, credMass=0.8)
mus.HDI95 <- apply(pred.mus, c(1,2), HDInterval::hdi)
tiff(height=4, width=6, units="in", res=300,
     filename= "C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Translocation effects.tiff")
par(mar=c(3,4,4,1), mfrow=c(1,1))
for (tr in 1:2){
  if (tr==1){
    plot(1:5, c(0,0,1,1,1), 
         type="n",
         pch=c(16, 17)[tr], cex=2,
         ylim=c(0,1),
         xlim=c(0.5, 5.5), 
         cex.axis=1, cex.lab=1.25,
         ylab="Probability", xlab="",
         main="",
         xaxt="n", yaxt="n")
  abline(h= seq(0,1,by=0.1), col="gray90")
  abline(v= seq(0.5,6,by=1), col="gray90")
  points(1:5+c(-0.1,0.1)[tr], mus.md[c(1,2,3,5,7),tr], 
         pch=c(16, 17)[tr], cex=1.5)
  }else{
    points(1:5+c(-0.1,0.1)[tr], mus.md[c(1,2,3,5,7),tr], 
           pch=c(16, 15)[tr], cex=1.5)
  }}

axis(1, at=1:5, cex.axis=0.85,
     labels=c("First-year\nSurvival", 
              "Nonbreeder\nSurvival", "Breeder\nSurvival", 
              "Nonbreeder\n to Breeder\nRecruitment", 
              "Nonbreeder\nDetection"), las=1, padj=0.5)
axis(2, at=seq(0,1, by=0.1), cex.axis=1,
     labels=c(0, NA,NA,NA,NA,0.5, NA,NA,NA,NA, 1))
for (m in 1:5){
  for (tr in 1:2){
    lines(rep(c(1:5)[m] + c(-0.1,0.1)[tr],2), mus.HDI95[1:2,c(1,2,3,5,7)[m],tr], lwd=2)
    lines(rep(c(1:5)[m] + c(-0.1,0.1)[tr],2), mus.HDI80[1:2,c(1,2,3,5,7)[m],tr], lwd=4)
  }}
legend(x=3,y=1.25, pch=c(16,15), pt.cex=1.5, cex=1, xpd=NA, horiz=T, 
       legend=c("Not translocated", "Translocated" ),
       xjust=0.5)
dev.off()

# Plot fecundity in response to nest treatments
f <- exp(outp$lmu.prod)
f.md <- apply(f, 1, median)
f.HDI80 <- apply(f, 1, HDInterval::hdi, credMass=0.8)
f.HDI95 <- apply(f, 1, HDInterval::hdi)

# Calculate treatment effects on fecundity
f.pred <- array(NA, dim=dim(outp$lmu.f))
for (s in 1:2){
  f.pred[s,] <- exp(outp$lmu.f[s,] + outp$gamma[,1])
} # s
f2.md <- apply(f.pred, 1, median)
f2.HDI80 <- apply(f.pred, 1, HDInterval::hdi, credMass=0.8)
f2.HDI95 <- apply(f.pred, 1, HDInterval::hdi)

tiff(height=4, width=6, units="in", res=300,
     filename= "C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Fecundity_treatment.tiff")
par(mar=c(3,5,3,1), mfrow=c(1,1))
plot(c(1, 2)-0.1, f.md, 
     ylim=c(min(f.HDI95), max(f2.HDI95)),
     xlim=c(0.5, 2.5), 
     ylab="Fecundity", xlab="",
     main="",
     cex.axis=1.5, cex.lab=1.5,
     xaxt="n", yaxt="n", type="n")
points(c(1, 2)-0.1, f.md, pch=16, cex=1.5)
abline(h=seq(0, 0.8, by=0.1), col="gray90")
abline(v=1.5, col="black")
axis(1, at=c(1,2), cex.axis=1,
     labels=c("Los Haitises\nNational Park","Punta Cana"),
     padj=0.5)
axis(2, at=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8), 
     cex.axis=1,
     labels=c(0, NA, NA, NA, 0.4, NA, NA, NA, 0.8))
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
legend(x=1.5,y=1.05, pch=c(16,18), pt.cex=c(1.5,2.5), cex=1,
       legend=c("Untreated", "Treated" ), 
       horiz=TRUE, xjust=0.5, xpd=NA)
dev.off()

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

p5 <- ggplot(data= lmus, aes(x=value, y=Parameter )) +
      stat_halfeye(point_interval=median_hdi, .width = c(.95, .80)) +
      coord_flip() +
      facet_wrap("Site") +
      xlab("Probability") + ylab("Parameter") +
      theme_bw(base_size=14) +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()) 
  
ggsave(filename="C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Survival Recruitment Detection.tiff",
       plot=p5,
       device="tiff",
       width=8, height=4, dpi=300)

#*******************
# Correlations between demographics and population growth rates
#*******************
lam.m <- apply(outp$lambda[1,1:12,,], c(1,2), median)
lam.hdi <- apply(outp$lambda[1,1:12,,], c(1,2), HDInterval::hdi)
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
plot.cor <- function (lam.post, x.post, x.lab, ind.x=1:12, site=1, yaxt="b"){
  # calculate the correlation coefficient 
  # over each iteration to propagate error
  lambda <- lam.post[1, 1:12, site,]
  x <- x.post[1,ind.x, site, ]
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
  cor.df <- data.frame(median= median(cor.post) |> round(2),
                       ldi=HDInterval::hdi(cor.post)[1] |> round(2),
                       hdi=HDInterval::hdi(cor.post)[2] |> round(2),
                       pd = pd(cor.post) |> round(2))
  
  lam.m <- apply(lambda, 1, median)
  lam.hdi <- apply(lambda, 1, HDInterval::hdi)
  x.m <- apply(x, 1, median)
  x.hdi <- apply(x, 1, HDInterval::hdi)
    x.lims <- c(min(x.hdi[,]), max(x.hdi[,]))
    y.lims <- c(min(lam.hdi[,]), max(lam.hdi[,]))
    if(yaxt=="n"){
    plot(x.m, lam.m[1:12], 
         xlim= x.lims,
         ylim= y.lims,
         type="n", 
         ylab="", xlab=x.lab, cex.lab=1.5, 
         main=c("", "")[site],
         yaxt=yaxt, xaxt="n")
      axis(2, at=df[,1], labels=c(rep(NA, 5)))
    } else {
      plot(x.m, lam.m[1:12], 
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
tiff(height=8, width=8, units="in", res=300,
     filename= "C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//PopGrowthand Demographics1.tiff")
par(mfrow=c(4,3), mar=c(4,1,1,1), oma=c(0,5,3,0))
plot.cor(outp$lambda, outp$mn.f, x.lab="Fecundity", ind.x=2:13)
plot.cor(outp$lambda, outp$mn.phiFY, x.lab="First-year Survival", yaxt="n")
plot.cor(outp$lambda, outp$mn.phiA, x.lab="Nonbreeder Survival", yaxt="n")
plot.cor(outp$lambda, outp$mn.phiB, x.lab="Breeder Survival")
plot.cor(outp$lambda, outp$mn.psiFYB, x.lab="First-year to Breeder", yaxt="n")
plot.cor(outp$lambda, outp$mn.psiAB, x.lab="Nonbreeder to Breeder", yaxt="n")
mtext("Population growth rate", side=2, outer=TRUE, line=2.5, cex=2)
mtext("Los Haitises", side=3, outer=TRUE, cex=2)
dev.off()

tiff(height=4, width=8, units="in", res=300,
     filename= "C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//PopGrowthand Demographics2.tiff")
par(mfrow=c(2,3), mar=c(4,1,1,1), oma=c(0,5,3,0))
plot.cor(outp$lambda, outp$mn.f, x.lab="Fecundity", ind.x=2:13, site=2)
plot.cor(outp$lambda, outp$mn.phiFY, x.lab="First-year Survival", yaxt="n", site=2)
plot.cor(outp$lambda, outp$mn.phiA, x.lab="Nonbreeder Survival", yaxt="n", site=2)
plot.cor(outp$lambda, outp$mn.phiB, x.lab="Breeder Survival", site=2)
plot.cor(outp$lambda, outp$mn.psiFYB, x.lab="First-year to Breeder", yaxt="n", site=2)
plot.cor(outp$lambda, outp$mn.psiAB, x.lab="Nonbreeder to Breeder", yaxt="n", site=2)
mtext("", side=2, outer=TRUE, line=2.5, cex=2)
mtext("Punta Cana", side=3, outer=TRUE, cex=2)
dev.off()
# Breeder to nonbreeder didn't vary over time
# Does the number of breeders correlate with pop growth?
plot.cor(outp$lambda, outp$NB, x.lab="Breeder Abundance", ind.x=2:13)
plot.cor(outp$lambda, outp$NF, x.lab="Nonbreeder Abundance", ind.x=2:13)
plot.cor(outp$lambda, outp$NFY, x.lab="First-year Abundance", ind.x=2:13)
plot.cor(outp$lambda, outp$Ntot, x.lab="All Stages Abundance", ind.x=2:13)


#*******************
#* Plot PVA results
#*******************
# Abundance of females at Los Haitises
library (viridis)
cols <- viridis(6)
pe <- apply(outp$extinct, c(1,2,3), mean)

tiff(height=3, width=8, units="in", res=300,
     filename= "C://Users//rolek.brian//OneDrive - The Peregrine Fund//Documents//Projects//Ridgways IPM//figs//Quasiextinction.tiff")
par(mar=c(2, 2, 1, 1), oma=c(3,3,0,0))
layout(matrix(c(1,1,2,2,3), 1, 5, byrow = TRUE))

plot(2024:2123, pe[1, ,1], 
     type="n", col=cols[1], lwd=2,
     ylim=c(0,1), main="Los Haitises",
     ylab="Quasi-extinction probability", xlab="Year",
     xaxt="n")
axis(1, at=seq(2020, 2120, by=10), 
     labels=c(2020, NA, 2040, NA, 2060, NA, 2080, NA, 2100, NA, 2120))
abline(h=seq(0, 1, by=0.1), col="gray90")
abline(v=seq(2030, 2130, by=10), col="gray90")
lines(2024:2123, pe[2, ,1], col=cols[1], lwd=2)
lines(2024:2123, pe[3, ,1], col=cols[2], lwd=2)
lines(2024:2123, pe[1, ,1], col=cols[3], lwd=2)
lines(2024:2123, pe[5, ,1], col=cols[4], lwd=2, lty=2)
lines(2024:2123, pe[6, ,1], col=cols[5], lwd=2, lty=2)
lines(2024:2123, pe[4, ,1], col=cols[6], lwd=2, lty=2)

plot(2024:2123, pe[1, ,2], 
     type="n", col=cols[1], lwd=2,
     ylim=c(0,1), main="Punta Cana",
     ylab="", xlab="Year",
     xaxt="n")
axis(1, at=seq(2020, 2120, by=10), 
     labels=c(2020, NA, 2040, NA, 2060, NA, 2080, NA, 2100, NA, 2120))
abline(h=seq(0, 1, by=0.1), col="gray90")
abline(v=seq(2030, 2130, by=10), col="gray90")
lines(2024:2123, pe[2, ,2], col=cols[1], lwd=2)
lines(2024:2123, pe[3, ,2], col=cols[2], lwd=2)
lines(2024:2123, pe[1, ,2], col=cols[3], lwd=2)
lines(2024:2123, pe[5, ,2], col=cols[4], lwd=2, lty=2)
lines(2024:2123, pe[6, ,2], col=cols[5], lwd=2, lty=3)
lines(2024:2123, pe[4, ,2], col=cols[6], lwd=2, lty=4)
plot.new()
legend(x=-0.3, y=0.5, title="Scenarios",
       legend=c("1 trans=0, nests=0", "2 trans=0, nests=10", "3 trans=0, nests=all",
                "4 trans=10, nests=0", "5 trans=10, nests=10", "6 trans=10, nests=all"),
       xpd=NA, col=cols, lty=c(1,1,1,2,2,2), lwd=2,
       xjust=0, yjust=0.5)
mtext("Future year", 1,  outer=TRUE, line=1, adj=0.39)
mtext("Quasi-extinction probability", 2, outer=TRUE, line=1)
dev.off()

# Time to extinction
# conditional on extinction threshold (D=X)
D <- 3
qextinct <- outp$Ntot <=3
ext <- outp$extinct[,100,,]==1
t.extinction <- apply(qextinct[])