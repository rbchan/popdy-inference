## Goodness-of-fit for N-mixture model with ZIP

nSites <- 200
nVisits <- 4
psi <- 0.5   # Proportion of extra-Poisson zeros
z <- rbinom(n=nSites, size=1, prob=psi) # Extra zeros
lam <- 5     # expected count when z=1
N <- rpois(n=nSites, lambda=lam*z)      # abundance at each site
p <- 0.3     # detection prob
y <- matrix(NA, nSites, nVisits)
for(i in 1:nSites) {
    y[i,] <- rbinom(n=nVisits, size=N[i], prob=p) # count data
}


library(unmarked)
umf <- unmarkedFramePCount(y=y)

fm <- pcount(~1~1, umf)
fmzip <- pcount(~1~1, umf, mixture="ZIP")

pb <- parboot(fm)
pbzip <- parboot(fmzip)

plot(pb)
plot(pbzip)  ## Should fit but doesn't

pb
pbzip
