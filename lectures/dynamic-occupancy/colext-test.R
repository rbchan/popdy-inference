library(unmarked)

set.seed(54598)
nSites <- 100
nPrimary <- 10
nSecondary <- 3    

beta0.psi <- -1; beta1.psi <- 1
elevation <- rnorm(nSites)
z2 <- psi2 <- matrix(NA, nSites, nPrimary)
psi2[,1] <- plogis(beta0.psi + beta1.psi*elevation)
z2[,1] <- rbinom(n=nSites, size=1, prob=psi2[,1])

epsilon2 <- 0.3 ## Local extinction prob
temperature <- matrix(rnorm(nSites*nPrimary)*elevation, nrow=nSites)    
beta0.gamma <- 1; beta1.gamma <- -1
gamma2 <- plogis(beta0.gamma + beta1.gamma*temperature)
for(k in 2:nPrimary) {
    psi2[,k] <- z2[,k-1]*(1-epsilon2) + (1-z2[,k-1])*gamma2[,k]
    z2[,k] <- rbinom(n=nSites, size=1, prob=psi2[,k])    }

p2 <- 0.2
y2 <- array(NA, c(nSites, nSecondary, nPrimary))
for(i in 1:nSites) {
    for(k in 1:nPrimary) {
        y2[i,,k] <- rbinom(nSecondary, size=1, prob=z2[i,k]*p2)
    } }

y2.wide <- matrix(y2, nrow=nrow(y2)) ## Format as nSites by (nSec*nPrimary)


## This shouldn't throw an error
umf <- unmarkedMultFrame(y=y2.wide, 
                         siteCovs=data.frame(elevation),
                         yearlySiteCovs=list(temp=temperature), ##[,-nPrimary]
                         numPrimary=nPrimary) 


## Try discarding final year
umf <- unmarkedMultFrame(y=y2.wide, 
                         siteCovs=data.frame(elevation),
                         obsCovs=list(temp=temperature[,-nPrimary]),
                         numPrimary=nPrimary) 


## This works
umf <- unmarkedMultFrame(y=y2.wide, 
                         siteCovs=data.frame(elevation),
##                         obsCovs=list(temp=temperature[,-nPrimary]),
                         numPrimary=nPrimary) 





