library (HDInterval)
# 2023 values
## Breeder survival
### LH and PC
apply(outp$mn.phiB[12,,], c(1), median)
apply(outp$mn.phiB[12,,], c(1), hdi)
### both sites averaged
phiB.mn <- apply(outp$mn.phiB[13,,], 2, mean)
median(phiB.mn)
hdi(phiB.mn)

## Non-Breeder survival
### LH and PC
apply(outp$mn.phiA[12,,], c(1), median)
apply(outp$mn.phiA[12,,], c(1), hdi)
### both sites averaged
phiA.mn <- apply(outp$mn.phiA[13,,], 2, mean)
median(phiA.mn)
hdi(phiA.mn)

## FY survival
### LH and PC
apply(outp$mn.phiFY[12,,], c(1), median)
apply(outp$mn.phiFY[12,,], c(1), hdi)
### both sites averaged
phiFY.mn <- apply(outp$mn.phiFY[12,,], 2, mean)
median(phiFY.mn)
hdi(phiFY.mn)

## Fecundity
### LH and PC
apply(outp$mn.f[13,,]*2, c(1), median)
apply(outp$mn.f[13,,]*2, c(1), hdi)
### both sites averaged
f23.mn <- apply(outp$mn.f[13,,]*2, 2, mean)
median(f23.mn)
hdi(f23.mn)

## Total abundance
### LH and PC
apply(outp$Ntot[13,,], c(1), median)
apply(outp$Ntot[13,,], c(1), hdi)
### both sites averaged
median(outp$Ntot[13,1,]+outp$Ntot[13,2,])
hdi(outp$Ntot[13,1,]+outp$Ntot[13,2,])

## Mature abundance
### LH and PC
apply(outp$NAD[13,,], c(1), median)
apply(outp$NAD[13,,], c(1), hdi)
### both sites averaged
median(outp$NAD[13,1,]+outp$NAD[13,2,])
hdi(outp$NAD[13,1,]+outp$NAD[13,2,])
