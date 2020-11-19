### R code from vignette source 'd:/courses/scr/lectures/open-JS/Open-JS.Rnw'

###################################################
### code chunk number 1: Open-JS.Rnw:207-215
###################################################
T <- 10      # years/primary periods
K <- 3       # 3 secondary sampling occasion
N0 <- 25     # Abundance in year 1
M <- 500     # Easiest way to simulate data is using data augmentation
phi <- 0.7   # Apparent survival
gamma <- 0.3 # Per-capital recruitment rate
p0 <- 0.4
sigma <- 0.1


###################################################
### code chunk number 2: Open-JS.Rnw:219-230
###################################################
set.seed(340)
co <- seq(0.25, 0.75, length=5)
X <- cbind(rep(co, each=5), rep(co, times=5))
J <- nrow(X)
xlim <- ylim <- c(0,1)
s <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
d <- p <- matrix(NA, M, J)
for(i in 1:M) {
    d[i,] <- sqrt((s[i,1]-X[,1])^2 + (s[i,2]-X[,2])^2)
    p[i,] <- p0*exp(-d[i,]^2/(2*sigma^2))
}


###################################################
### code chunk number 3: Open-JS.Rnw:242-260
###################################################
set.seed(034)
z <- recruitable <- died <- recruited <- matrix(0, M, T)
z[1:N0,1] <- 1 # First N0 are alive
recruitable[(N0+1):M,1] <- 1
for(t in 2:T) {
    prevN <- sum(z[,t-1]) # number alive at t-1
    ER <- prevN*gamma # expected number of recruits
    prevA <- sum(recruitable[,t-1]) # Number available to be recruited
    gammaPrime <- ER/prevA
    if(gammaPrime > 1)
        stop("M isn't big enough")
    for(i in 1:M) {
        z[i,t] <- rbinom(1, 1, z[i,t-1]*phi + recruitable[i,t-1]*gammaPrime)
        recruitable[i,t] <- 1 - max(z[i,1:(t)]) # to be recruited
        died[i,t] <- z[i,t-1]==1 & z[i,t]==0
        recruited[i,t] <- z[i,t]==1 & z[i,t-1]==0
    }
}


###################################################
### code chunk number 4: Open-JS.Rnw:265-269
###################################################
N <- colSums(z) # Population size
Deaths <- colSums(died)
Recruits <- colSums(recruited)
everAlive <- sum(rowSums(z)>0)


###################################################
### code chunk number 5: Open-JS.Rnw:285-293
###################################################
yall <- array(NA, c(M, J, K, T))
for(i in 1:M) {
    for(t in 1:T) {
        for(j in 1:J) {
            yall[i,j,1:K,t] <- rbinom(K, 1, z[i,t]*p[i,j])
        }
    }
}


###################################################
### code chunk number 6: Open-JS.Rnw:298-301
###################################################
detected <- rowSums(yall) > 0
y <- yall[detected,,,]
str(y)


###################################################
### code chunk number 7: NDR
###################################################
plot(1:T, N, ylim=c(0, 50), type="o", pch=16,
     xlab="Year", ylab="")
lines(2:T, Deaths[-1], col="red", type="o", pch=16)
lines(2:T, Recruits[-1], col="blue", type="o", pch=16)
legend(1, 50, c("Population size", "Deaths", "Recruits"),
       col=c("black", "red", "blue"), pch=16, lty=1)


###################################################
### code chunk number 8: Open-JS.Rnw:344-347
###################################################
M <- nrow(y) + 50
yz <- array(0, c(M, J, K, T))
yz[1:nrow(y),,,] <- y


###################################################
### code chunk number 9: Open-JS.Rnw:352-355
###################################################
zi <- matrix(0, M, T)
zi[1:nrow(y),] <- 1
ji1 <- function() list(phi=0.01, gamma=0.01, z=zi)


###################################################
### code chunk number 10: Open-JS.Rnw:360-366 (eval = FALSE)
###################################################
## jd1 <- list(y=yz, M=M, X=X,
##             J=J, K=K, T=T, xlim=xlim, ylim=ylim)
## jp1 <- c("phi", "gamma", "p0", "sigma", "N", "Deaths", "Recruits", "Ntot")
## library(rjags)
## jm1 <- jags.model("JS-spatial.jag", jd1, ji1, n.chains=1, n.adapt=500)
## jc1 <- coda.samples(jm1, jp1, 1000)


###################################################
### code chunk number 11: Ntot
###################################################
hist(as.matrix(jc1[,"Ntot"]), xlab="Total population size", ylab="", main="", freq=FALSE, xlim=c(nrow(y), M))
abline(v=M, lwd=3, col="blue")


###################################################
### code chunk number 12: jc1
###################################################
plot(jc1[,c("phi", "gamma", "p0", "sigma")])


###################################################
### code chunk number 13: jcN1-4
###################################################
plot(jc1[,c("N[1]", "N[2]", "N[3]", "N[4]")])


###################################################
### code chunk number 14: jcN5-8
###################################################
plot(jc1[,c("N[5]", "N[6]", "N[7]", "N[8]")])


###################################################
### code chunk number 15: dd1
###################################################
N <- 0:50
nu0 <- 2
nu1 <- 0.05
plot(N, nu0*exp(-nu1*N), type="l", ylim=c(0,2), ylab="Per-capita Recruitment")


###################################################
### code chunk number 16: Open-JS.Rnw:509-521
###################################################
T <- 10      # years/primary periods
K <- 3       # 3 secondary sampling occasion
## Is it necessary to be far from equilibrium to detect density-dependence?
## Equilibrium here is where (1-phi) == gamma, where gamma is function of N
N0 <- 10     # Abundance in year 1
M <- 500     # Easiest way to simulate data is using data augmentation
phi <- 0.7   # Apparent survival
##gamma <- 0.3 # Per-capital recruitment rate
nu0 <- 2
nu1 <- 0.05
p0 <- 0.4
sigma <- 0.1


###################################################
### code chunk number 17: Open-JS.Rnw:525-536
###################################################
set.seed(3479)
co <- seq(0.25, 0.75, length=5)
X <- cbind(rep(co, each=5), rep(co, times=5))
J <- nrow(X)
xlim <- ylim <- c(0,1)
s <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
d <- p <- matrix(NA, M, J)
for(i in 1:M) {
    d[i,] <- sqrt((s[i,1]-X[,1])^2 + (s[i,2]-X[,2])^2)
    p[i,] <- p0*exp(-d[i,]^2/(2*sigma^2))
}


###################################################
### code chunk number 18: Open-JS.Rnw:548-567
###################################################
set.seed(3401)
z2 <- recruitable <- died <- recruited <- matrix(0, M, T)
z2[1:N0,1] <- 1 # First N0 are alive
recruitable[(N0+1):M,1] <- 1
for(t in 2:T) {
    prevN <- sum(z2[,t-1]) # number alive at t-1
    gamma <- nu0*exp(-nu1*prevN) ## Density dependent recruitment rate
    ER <- prevN*gamma # expected number of recruits
    prevA <- sum(recruitable[,t-1]) # Number available to be recruited
    gammaPrime <- ER/prevA
    if(gammaPrime > 1)
        stop("M isn't big enough")
    for(i in 1:M) {
        z2[i,t] <- rbinom(1, 1, z2[i,t-1]*phi + recruitable[i,t-1]*gammaPrime)
        recruitable[i,t] <- 1 - max(z2[i,1:(t)]) # to be recruited
        died[i,t] <- z2[i,t-1]==1 & z2[i,t]==0
        recruited[i,t] <- z2[i,t]==1 & z2[i,t-1]==0
    }
}


###################################################
### code chunk number 19: Open-JS.Rnw:572-576
###################################################
N2 <- colSums(z2) # Population size
Deaths2 <- colSums(died)
Recruits2 <- colSums(recruited)
everAlive2 <- sum(rowSums(z2)>0)


###################################################
### code chunk number 20: Open-JS.Rnw:592-600
###################################################
yall <- array(NA, c(M, J, K, T))
for(i in 1:M) {
    for(t in 1:T) {
        for(j in 1:J) {
            yall[i,j,1:K,t] <- rbinom(K, 1, z2[i,t]*p[i,j])
        }
    }
}


###################################################
### code chunk number 21: Open-JS.Rnw:605-608
###################################################
detected <- rowSums(yall) > 0
y2 <- yall[detected,,,]
str(y2)


###################################################
### code chunk number 22: NDR-DD
###################################################
plot(1:T, N2, ylim=c(0, 50), type="o", pch=16,
     xlab="Year", ylab="")
lines(2:T, Deaths2[-1], col="red", type="o", pch=16)
lines(2:T, Recruits2[-1], col="blue", type="o", pch=16)
legend(1, 50, c("Population size", "Deaths", "Recruits"),
       col=c("black", "red", "blue"), pch=16, lty=1)


###################################################
### code chunk number 23: Open-JS.Rnw:651-654
###################################################
M2 <- nrow(y2) + 75
yz2 <- array(0, c(M2, J, K, T))
yz2[1:nrow(y2),,,] <- y2


###################################################
### code chunk number 24: Open-JS.Rnw:659-663
###################################################
zi <- matrix(0, M2, T)
##zi[1:nrow(y2),] <- 1
zi[1:nrow(y2),] <- z2[detected,] ## cheating
ji2 <- function() list(phi=0.01, z=zi)


###################################################
### code chunk number 25: Open-JS.Rnw:668-675 (eval = FALSE)
###################################################
## jd2 <- list(y=yz2, M=M2, X=X,
##             J=J, K=K, T=T, xlim=xlim, ylim=ylim)
## jp2 <- c("phi", "nu0", "nu1", "p0", "sigma", "N", "Deaths", "Recruits", "Ntot")
## library(rjags)
## jm2 <- jags.model("JS-spatial-DD.jag", jd2, ji2, n.chains=1, n.adapt=200)
## jc2 <- coda.samples(jm2, jp2, 5000)
## ##jc2.2 <- coda.samples(jm2, jp2, 15000)


###################################################
### code chunk number 26: jc2
###################################################
plot(jc2[,c("phi", "nu0", "nu1")])


###################################################
### code chunk number 27: Open-JS.Rnw:705-709
###################################################
Npost <- as.matrix(jc2[,paste("N[", 1:10, "]", sep="")])
Nmed <- apply(Npost, 2, median)
Nupper <- apply(Npost, 2, quantile, prob=0.975)
Nlower <- apply(Npost, 2, quantile, prob=0.025)


###################################################
### code chunk number 28: Npost
###################################################
plot(1:T, N2, type="o", col="blue", ylim=c(0, 100), xlab="Time",
     ylab="Abundance")
points(1:T, Nmed)
arrows(1:T, Nlower, 1:T, Nupper, angle=90, code=3, length=0.05)
legend(1, 100, c("Actual abundance", "Estimated abundance"),
       col=c("blue", "black"), lty=c(1,1), pch=c(1,1))


