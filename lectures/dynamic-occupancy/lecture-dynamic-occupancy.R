#' ---
#' format:
#'   revealjs:
#'     theme: default
#'     pagetitle: Dynamic occupancy models
#'     slide-number: true
#'     echo: true
#' lang: en
#' title: lecture-dynamic-occupancy
#' ---
#' 

#| label: sim-init
set.seed(54598)
nSites <- 100
nPrimary <- 10
z <- psi <- matrix(NA, nSites, nPrimary)
psi[,1] <- 0.5  ## Initial occupancy prob
z[,1] <- rbinom(n=nSites, size=1, prob=psi[,1])

#| label: nsim-dy
epsilon <- 0.3 ## Local extinction prob
gamma <- 0.2   ## Local colonization prob
for(k in 2:nPrimary) {
    psi[,k] <- z[,k-1]*(1-epsilon) + (1-z[,k-1])*gamma
    z[,k] <- rbinom(n=nSites, size=1, prob=psi[,k])    }

#| label: sim-cov3
nSecondary <- 3    
p <- 0.2
y <- array(NA, c(nSites, nSecondary, nPrimary))
for(i in 1:nSites) {
    for(k in 1:nPrimary) {
        y[i,,k] <- rbinom(nSecondary, size=1, prob=z[i,k]*p)
    } }

#| label: sim-nocov-dat
y[1:15,,1]

#| label: sim-nocov-ss1
siteDets <- rowSums(y) # Dets at each site
table(siteDets)        # Frequency

#| label: sim-nocov-ss2
yearDets <- apply(y, 3, sum)
yearDets

#| label: un
#| include: false
library(unmarked)

#| label: un-umf
## Format as nSites by (nSec*nPrimary) matrix
y.wide <- matrix(y, nrow=nrow(y)) 
umf <- unmarkedMultFrame(y=y.wide, numPrimary=nPrimary)

#| label: wfac
#| size: scriptsize
summary(umf)

#| label: un-fit
fm <- colext(~1,~1,~1,~1, umf)    
fm

#| label: un-compare
coef(backTransform(fm, type="psi"))
coef(backTransform(fm, type="col"))
coef(backTransform(fm, type="ext"))
coef(backTransform(fm, type="det"))

#| label: un-actual
psi[1,1]
gamma
epsilon
p

#| label: re
#| fig-height: 5
#| fig-width: 7
#| out-width: 60%
#| fig-align: center
#| code-fold: true
par(mai=c(0.9, 0.9, 0.1, 0.1))  
re <- ranef(fm)
occupied.post <- predict(re, func=colSums, nsim=1000)
plot(1:nPrimary, rowMeans(occupied.post), type="b", pch=16,
     xlab="Time", ylab="Sites occupied", ylim=c(0, 70))
lines(1:nPrimary, colSums(z), type="b", col=4)
segments(1:nPrimary, apply(occupied.post, 1, quantile, prob=0.025),
         1:nPrimary, apply(occupied.post, 1, quantile, prob=0.975))
legend(1, 70, c("Actual", "Estimated"), pch=c(1,16), col=c(4,1), lty=1)

#| label: bugs
#| echo: false
writeLines(readLines("dynocc-model.jag"))

#| label: bugs-data
jags.data <- list(y=y, nSites=nSites,
                  J=nSecondary, K=nPrimary)

#| label: bugs-inits
jags.inits <- function() {
    list(psi1=runif(1), epsilon=runif(1),
         gamma=runif(1), p=runif(1),
         z=matrix(1, nSites, nPrimary))
}

#| label: bugs-pars
jags.pars <- c("psi1", "epsilon", "gamma", "p", "N")

#| label: jagsUI
#| include: false
#| results: hide
library(jagsUI)
library(coda)

#| label: bugs-mcmc
#| results: hide
#| message: false
#| cache: true
library(jagsUI)
jags.post.samples <- jags.basic(data=jags.data, inits=jags.inits,
                                parameters.to.save=jags.pars,
                                model.file="dynocc-model.jag",
                                n.chains=3, n.adapt=100, n.burnin=0,
                                n.iter=2000, parallel=TRUE)

#| label: bugs-sum
summary(jags.post.samples[,jags.pars[1:4]])

#| label: bugs-plot1
#| fig-height: 7
#| fig-width: 7
#| out-width: 80%
#| fig-align: center
plot(jags.post.samples[,jags.pars[1:4]])

#| label: sim-init-cov
set.seed(5098)
beta0.psi <- -1; beta1.psi <- 1
elevation <- rnorm(nSites) ## site covariate
z2 <- psi2 <- matrix(NA, nSites, nPrimary)
psi2[,1] <- plogis(beta0.psi + beta1.psi*elevation)
z2[,1] <- rbinom(n=nSites, size=1, prob=psi2[,1])

#| label: nsim-dy-cov
epsilon2 <- 0.3 ## Local extinction prob
temperature <- matrix(rnorm(nSites*nPrimary)*elevation, ## Primary period covariate
                      nrow=nSites) 
beta0.gamma <- -2; beta1.gamma <- -1
gamma2 <- plogis(beta0.gamma + beta1.gamma*temperature)
for(k in 2:nPrimary) {
    psi2[,k] <- z2[,k-1]*(1-epsilon2) + (1-z2[,k-1])*gamma2[,k-1]
    z2[,k] <- rbinom(n=nSites, size=1, prob=psi2[,k])
}

#| label: sim-cov3-cov
p2 <- 0.2
y2 <- array(NA, c(nSites, nSecondary, nPrimary))
for(i in 1:nSites) {
    for(k in 1:nPrimary) {
        y2[i,,k] <- rbinom(nSecondary, size=1, prob=z2[i,k]*p2)
    } }

#| label: un-umf-cov
y2.wide <- matrix(y2, nrow=nrow(y2)) ## nSites by (nSec*nPrimary) matrix
umf2 <- unmarkedMultFrame(y=y2.wide, numPrimary=nPrimary, siteCovs=data.frame(elevation),
                          yearlySiteCovs=list(temp=temperature))

#| label: wfac-cov
summary(umf2)

#| label: un-fit-cov
fm2 <- colext(~elevation,~temp,~1,~1, umf2)    
fm2

#| label: preddat
pred.data <- data.frame(temp=seq(from=-5, to=4, length=50))

#| label: pred-gamma
gamma.pred <- predict(fm2, newdata=pred.data,
                      type='col', append=TRUE)

#| label: pred-gamma1
#| fig-height: 7
#| fig-width: 9
#| out-width: 60%
#| fig-align: center
#| code-fold: true
plot(Predicted ~ temp, gamma.pred, type="l", ylab="Colonization probability",
     xlab="Standardized temperature", ylim=0:1, lwd=2) 
lines(lower ~ temp, gamma.pred, col=gray(0.6))
lines(upper ~ temp, gamma.pred, col=gray(0.6))

#| label: bugs-cov
#| echo: false
#| comment: ''
#| background: beige
writeLines(readLines("dynocc-model-covars.jag"))

#| label: bugs-data-cov
jags.data.covs <- list(y=y2, nSites=nSites, elevation=elevation,
                       temp=temperature, J=nSecondary, K=nPrimary)

#| label: bugs-inits-cov
jags.inits.covs <- function() {
    list(beta0.psi=rnorm(1), beta1.psi=rnorm(1), epsilon=runif(1),
         beta0.gamma=rnorm(1), beta1.gamma=rnorm(1), p=runif(1),
         z=matrix(1, nSites, nPrimary))
}

#| label: bugs-pars-cov
jags.pars.covs <- c("beta0.psi", "beta1.psi", "epsilon",
                    "beta0.gamma", "beta1.gamma", "p", "N")

#| label: bugs-mcmc-cov
#| results: hide
#| message: false
#| cache: true
library(jagsUI)
jags.post.samples.covs <- jags.basic(
    data=jags.data.covs, inits=jags.inits.covs,
    parameters.to.save=jags.pars.covs,
    model.file="dynocc-model-covars.jag",
    n.chains=3, n.adapt=100, n.burnin=0, n.iter=2000, parallel=TRUE)

#| label: bugs-sum-cov
summary(jags.post.samples.covs[,jags.pars.covs[1:6]])

#| label: bugs-plot1-cov
#| out-width: 70%
#| fig-align: center
plot(jags.post.samples.covs[,jags.pars.covs[1:3]])

#| label: bugs-plot2-cov
#| out-width: 70%
#| fig-align: center
plot(jags.post.samples.covs[,jags.pars.covs[4:6]])

#| label: N-post-covs
library(coda)
N.names <- grep("N\\[", varnames(jags.post.samples.covs)  )
N.post <- as.matrix(jags.post.samples.covs[,N.names])
N.post.mean <- colMeans(N.post)
N.post.lower <- apply(N.post, 2, quantile, prob=0.025)
N.post.upper <- apply(N.post, 2, quantile, prob=0.975)

#| label: N-post-covs-plot
#| fig-height: 7
#| fig-width: 9
#| out-width: 65%
#| fig-align: center
#| code-fold: true
par(mai=c(0.9, 0.9, 0.1, 0.1))  
plot(1:nPrimary, N.post.mean, type="b", ylim=c(0, 70), pch=16,
     xlab="Time", ylab="Sites occupied")
segments(1:nPrimary, N.post.lower, 1:nPrimary, N.post.upper)
lines(1:nPrimary, colSums(z2), type="b", col=4)
legend(1, 70, c("Actual", "Estimated"), pch=c(1,16), col=c(4,1), lty=1)

#| label: beta-post
beta.post <- as.matrix(jags.post.samples.covs[,c(
    "beta0.gamma","beta1.gamma")])
temperature.seq <- seq(-5, 5, length=50)
n.samples <- nrow(beta.post)
gamma.post <- matrix(NA, n.samples, length(temperature.seq))
for(i in 1:n.samples) {
    gamma.post[i,] <- plogis(beta.post[i,1] +
                             beta.post[i,2]*temperature.seq)
}
gamma.post.mean <- colMeans(gamma.post)
gamma.post.lower <- apply(gamma.post, 2, quantile, prob=0.025)
gamma.post.upper <- apply(gamma.post, 2, quantile, prob=0.975)

#| label: gamma-post-covs-plot
#| fig-height: 5
#| fig-width: 7
#| out-width: 65%
#| fig-align: center
#| code-fold: true
par(mai=c(0.9, 0.9, 0.1, 0.1))  
matplot(temperature.seq, t(gamma.post[seq(1,n.samples,10),]), type="l",
        xlab="Temperature (standardized)", ylab="Colonization prob", col=gray(0.9))
lines(temperature.seq, gamma.post.mean, lwd=3, col="blue")
lines(temperature.seq, gamma.post.lower, col="blue")
lines(temperature.seq, gamma.post.upper, col="blue")

