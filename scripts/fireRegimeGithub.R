##################################################################################
## hypothesis test of shifting fire intensity from 'natural' to 'anthropogenic' ## 
## fire regime in the Late Pleistocene of Sahul                                 ##
##                                                                              ##
## Corey J. A. Bradshaw, Flinders University/CABAH                              ##
## January 2023                                                                 ##
##                                                                              ##
## code for paper:                                                              ##
## Bird, M, CJA Bradshaw, M Brand, R Comley, X Hadeen, Z Jacobs, CR Rowe,       ##
##    CM Wurster, C Zwart. Pleistocene establishment of an Anthropogenic        ##
##    fire regime in Australia’s tropical savannas.                             ##
##################################################################################


## import data
d2c <- read.table("d2c.csv", header=T, sep=",") # distance to coast
toc <- read.table("toc.csv", header=T, sep=",") # total organic carbon
spacmar <- read.table("SPACMAR.csv", header=T, sep=",") # mass accumulation rate of stable polycyclic aromatic carbon 
charpar <- read.table("charpar.csv", header=T, sep=",") # particle accumulation rate of micro-charcoal 
pcC4 <- read.table("pcC4.csv", header=T, sep=",") # % grass pollen (of total dryland pollen)
spacC13 <- read.table("SPACC13.csv", header=T, sep=",") # δ13C of the stable polycyclic aromatic carbon

# interval for standardised ages
agest <- seq(0,150000,250)

# d2c vs. toc
plot(d2c$age,d2c$d2c, pch=19, cex=0.3, type="l", col="red")
par(new=T,xaxt="n")
plot(toc$age,toc$TOC, pch=19, cex=0.3, type="l")

d2c.intp <- approx(x=d2c$age, y=d2c$d2c, xout=agest)
toc.intp <- approx(x=toc$age, y=toc$TOC, xout=agest)
d2ctoc.dat <- na.omit(data.frame(d2c.intp,toc.intp$y))
colnames(d2ctoc.dat) <- c("age","d2c","toc")
head(d2ctoc.dat)
plot(d2ctoc.dat$d2c, d2ctoc.dat$toc, pch=19, cex=0.2, xlab="distance to coast", ylab="TOC")

# spacmar vs. charpar
spacmar.intp <- approx(x=spacmar$age, y=spacmar$SPACMAR, xout=agest, ties="mean")
charpar.intp <- approx(x=charpar$age, y=charpar$charpar, xout=agest, ties="mean")
spacmarcharpar.dat <- data.frame(spacmar.intp,charpar.intp$y)
colnames(spacmarcharpar.dat) <- c("age","spacmar","charpar")
head(spacmarcharpar.dat)
plot(spacmarcharpar.dat$spacmar, spacmarcharpar.dat$charpar, pch=19, cex=0.2, xlab="SPAC MAR", ylab="char PAR")

# pcC4 vs. spacC13
pcC4.intp <- approx(x=pcC4$age, y=pcC4$pcC4, xout=agest, ties="mean")
spacC13.intp <- approx(x=spacC13$age, y=spacC13$SPACC13, xout=agest, ties="mean")
pcC4spacC13.dat <- data.frame(pcC4.intp,spacC13.intp$y)
colnames(pcC4spacC13.dat) <- c("age","pcC4","spac13")
head(pcC4spacC13.dat)
plot(pcC4spacC13.dat$pcC4, pcC4spacC13.dat$spac13, pch=19, cex=0.2, xlab="SPAC MAR", ylab="char PAR")
write.table(pcC4spacC13.dat, "pcC4spacC13.csv", sep=",")

# split by pre and post 30 ka
# pcC4 vs. spacC13
pcC4spacC13old <- subset(pcC4spacC13.dat, age > 30000)
pcC4spacC13young <- subset(pcC4spacC13.dat, age <= 30000)
par(mfrow=c(1,2))
plot((pcC4spacC13old$pcC4),pcC4spacC13old$spac13, pch=19, cex=0.3)
plot((pcC4spacC13young$pcC4),pcC4spacC13young$spac13, pch=19, cex=0.3)
par(mfrow=c(1,1))

cor.young <- cor(pcC4spacC13young$pcC4,pcC4spacC13young$spac13, method="spearman")
cor.old <- cor(pcC4spacC13old$pcC4,pcC4spacC13old$spac13, method="spearman")

lyoung <- dim(pcC4spacC13young)[1]
seqold <- 1:dim(pcC4spacC13old)[1]
iter <- 10000
corRan <- rep(0,iter)
for (i in 1:iter) {
  ran.sub <- sample(seqold,lyoung,replace=T)  
  corRan.old <- cor((pcC4spacC13old[ran.sub,]$pcC4), (pcC4spacC13old[ran.sub,]$spac13), method="spearman", use="pairwise.complete.obs")
  corRan[i] <- ifelse(corRan.old < cor.young, 1, 0)
}
sum(corRan)/iter

# spacmar vs. charpar
spacmarcharparold <- subset(spacmarcharpar.dat, age > 30000 & spacmar > 0 & charpar > 0)
spacmarcharparyoung <- subset(spacmarcharpar.dat, age <= 30000 & spacmar > 0 & charpar > 0)
par(mfrow=c(1,2))
plot(log10(spacmarcharparold$spacmar),log10(spacmarcharparold$charpar), pch=19, cex=0.3)
plot(log10(spacmarcharparyoung$spacmar),log10(spacmarcharparyoung$charpar), pch=19, cex=0.3)
par(mfrow=c(1,1))

cor.young <- cor(log10(spacmarcharparyoung$spacmar),log10(spacmarcharparyoung$charpar), method="spearman")
cor.young
cor.old <- cor(log10(spacmarcharparold$spacmar),log10(spacmarcharparold$charpar), method="spearman")
cor.old

lyoung <- dim(spacmarcharparyoung)[1]
seqold <- 1:dim(spacmarcharparold)[1]
iter <- 10000
corRan <- rep(0,iter)
for (i in 1:iter) {
  ran.sub <- sample(seqold,lyoung,replace=T)  
  corRan.old <- cor(log10(spacmarcharparold[ran.sub,]$spacmar), log10(spacmarcharparold[ran.sub,]$charpar), method="spearman")
  corRan[i] <- ifelse(corRan.old < cor.young, 1, 0)
}
sum(corRan)/iter


## only 'wet' periods (> 10% total organic carbon)
comb.dat <- na.omit(data.frame(d2c.intp,toc.intp$y, spacmar.intp$y, charpar.intp$y, pcC4.intp$y, spacC13.intp$y))
colnames(comb.dat) <- c("age","d2c","toc","spacmar","charpar","pcC4","spacC13")
head(comb.dat)
combwet <- subset(comb.dat, toc > 10)  
combwetyoung <- subset(combwet, age < 30000)  
combwetold <- subset(combwet, age >= 30000)  

# spacmar vs. charpar
par(mfrow=c(1,2))
plot(log10(combwetold$spacmar),log10(combwetold$charpar), pch=19, cex=0.3)
plot(log10(combwetyoung$spacmar),log10(combwetyoung$charpar), pch=19, cex=0.3)
par(mfrow=c(1,1))

cor.young <- cor(log10(combwetyoung$spacmar),log10(combwetyoung$charpar), method="spearman")
cor.young
cor.old <- cor(log10(combwetold$spacmar),log10(combwetold$charpar), method="spearman")
cor.old

lyoung <- dim(combwetyoung)[1]
seqold <- 1:dim(combwetold)[1]
iter <- 10000
corRan <- corRan.old <- rep(0,iter)
for (i in 1:iter) {
  ran.sub <- sample(seqold,lyoung,replace=T)  
  corRan.old[i] <- cor(log10(combwetold[ran.sub,]$spacmar), log10(combwetold[ran.sub,]$charpar), method="spearman")
  corRan[i] <- ifelse(corRan.old[i] < cor.young, 1, 0)
}
sum(corRan)/iter
hist(corRan.old,main="")
abline(v=cor.young, col="red", lty=2)


# pcC4 vs spacC13
par(mfrow=c(1,2))
plot((combwetold$pcC4),(combwetold$spacC13), pch=19, cex=0.3)
plot((combwetyoung$pcC4),(combwetyoung$spacC13), pch=19, cex=0.3)
par(mfrow=c(1,1))

cor.young <- cor(combwetyoung$pcC4,(combwetyoung$spacC13), method="spearman")
cor.young
cor.old <- cor(combwetold$pcC4,(combwetold$spacC13), method="spearman")
cor.old

lyoung <- dim(combwetyoung)[1]
seqold <- 1:dim(combwetold)[1]
iter <- 10000
corRan <- corRan.old <- rep(0,iter)
for (i in 1:iter) {
  ran.sub <- sample(seqold,lyoung,replace=T)  
  corRan.old[i] <- cor(combwetold[ran.sub,]$pcC4, (combwetold[ran.sub,]$spacC13), method="spearman")
  corRan[i] <- ifelse(corRan.old[i] > cor.young, 1, 0)
}
sum(corRan)/iter
hist(corRan.old,main="")
abline(v=cor.young, col="red", lty=2)


## wetness threshold sensitivity analysis
iter <- 10000
wet.thresh.vec <- seq(1,10,1)

pr.ran1 <- pr.ran2 <- rep(NA,length(wet.thresh.vec))

for (w in 1:length(wet.thresh.vec)) {
  combwet <- subset(comb.dat, toc > wet.thresh.vec[w])  
  combwetyoung <- subset(combwet, age < 30000)  
  combwetold <- subset(combwet, age >= 30000)  

  # MARspac vs. PARchar
  cor.young1 <- cor(log10(combwetyoung$spacmar),log10(combwetyoung$charpar), method="spearman")
  cor.old1 <- cor(log10(combwetold$spacmar),log10(combwetold$charpar), method="spearman")

  lyoung1 <- dim(combwetyoung)[1]
  seqold1 <- 1:dim(combwetold)[1]
  corRan1 <- corRan.old1 <- rep(0,iter)
  for (i in 1:iter) {
    ran.sub <- sample(seqold1,lyoung1,replace=T)  
    corRan.old1[i] <- cor(log10(combwetold[ran.sub,]$spacmar), log10(combwetold[ran.sub,]$charpar), method="spearman")
    corRan1[i] <- ifelse(corRan.old1[i] < cor.young1, 1, 0)
  }
  pr.ran1[w] <- sum(corRan1)/iter

  
  # pcC4 vs δC13spac
  cor.young2 <- cor(combwetyoung$pcC4,combwetyoung$spacC13, method="spearman")
  cor.old2 <- cor(combwetold$pcC4,combwetold$spacC13, method="spearman")

  lyoung2 <- dim(combwetyoung)[1]
  seqold2 <- 1:dim(combwetold)[1]
  corRan2 <- corRan.old2 <- rep(0,iter)
  for (i in 1:iter) {
    ran.sub <- sample(seqold2,lyoung2,replace=T)  
    corRan.old2[i] <- cor(combwetold[ran.sub,]$pcC4, combwetold[ran.sub,]$spacC13, method="spearman")
    corRan2[i] <- ifelse(corRan.old2[i] > cor.young2, 1, 0)
  }
  pr.ran2[w] <- sum(corRan2)/iter
  
  print(wet.thresh.vec[w])
  
} # end w

par(mfrow=c(1,2))
plot(wet.thresh.vec, pr.ran1, type="l", xlab="wetness threshold", ylab="Pr(random) MARspac vs. PARchar")
plot(wet.thresh.vec, pr.ran2, type="l", xlab="wetness threshold", ylab="Pr(random) pcC4 vs δC13spac")
par(mfrow=c(1,1))



## temporal split sensitivity analysis
iter <- 10000
tsplit.thresh.vec <- seq(50000,10000,-5000)

pr.ran1 <- pr.ran2 <- rep(NA,length(tsplit.thresh.vec))

for (t in 1:length(tsplit.thresh.vec)) {
  combwet <- subset(comb.dat, toc > 10)  
  combwetyoung <- subset(combwet, age < tsplit.thresh.vec[t])  
  combwetold <- subset(combwet, age >= tsplit.thresh.vec[t])  
  
  # MARspac vs. PARchar
  cor.young1 <- cor(log10(combwetyoung$spacmar),log10(combwetyoung$charpar), method="spearman")
  cor.old1 <- cor(log10(combwetold$spacmar),log10(combwetold$charpar), method="spearman")
  
  lyoung1 <- dim(combwetyoung)[1]
  seqold1 <- 1:dim(combwetold)[1]
  corRan1 <- corRan.old1 <- rep(0,iter)
  for (i in 1:iter) {
    ran.sub <- sample(seqold1,lyoung1,replace=T)  
    corRan.old1[i] <- cor(log10(combwetold[ran.sub,]$spacmar), log10(combwetold[ran.sub,]$charpar), method="spearman")
    corRan1[i] <- ifelse(corRan.old1[i] < cor.young1, 1, 0)
  }
  pr.ran1[t] <- sum(corRan1)/iter
  
  
  # pcC4 vs δC13spac
  cor.young2 <- cor(combwetyoung$pcC4,combwetyoung$spacC13, method="spearman")
  cor.old2 <- cor(combwetold$pcC4,combwetold$spacC13, method="spearman")
  
  lyoung2 <- dim(combwetyoung)[1]
  seqold2 <- 1:dim(combwetold)[1]
  corRan2 <- corRan.old2 <- rep(0,iter)
  for (i in 1:iter) {
    ran.sub <- sample(seqold2,lyoung2,replace=T)  
    corRan.old2[i] <- cor(combwetold[ran.sub,]$pcC4, combwetold[ran.sub,]$spacC13, method="spearman")
    corRan2[i] <- ifelse(corRan.old2[i] > cor.young2, 1, 0)
  }
  pr.ran2[t] <- sum(corRan2)/iter
  
  print(tsplit.thresh.vec[t])
  
} # end t

par(mfrow=c(1,2))
plot(tsplit.thresh.vec, pr.ran1, type="l", xlab="temporal split", ylab="Pr(random) MARspac vs. PARchar")
plot(tsplit.thresh.vec, pr.ran2, type="l", xlab="temporal split", ylab="Pr(random) pcC4 vs δC13spac")
par(mfrow=c(1,1))


