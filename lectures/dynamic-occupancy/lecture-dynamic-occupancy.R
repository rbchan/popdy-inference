## ----sim-init,size='scriptsize',echo=-1---------------------------------------
set.seed(54598)
nSites <- 100
nPrimary <- 10
z <- psi <- matrix(NA, nSites, nPrimary)
psi[,1] <- 0.5  ## Initial occupancy prob
z[,1] <- rbinom(n=nSites, size=1, prob=psi[,1])


## ----nsim-dy,size='scriptsize'------------------------------------------------
epsilon <- 0.3 ## Local extinction prob
gamma <- 0.2   ## Local colonization prob
for(k in 2:nPrimary) {
    psi[,k] <- z[,k-1]*(1-epsilon) + (1-z[,k-1])*gamma
    z[,k] <- rbinom(n=nSites, size=1, prob=psi[,k])    }


## ----sim-cov3,size='scriptsize'-----------------------------------------------
nSecondary <- 3    
p <- 0.2
y <- array(NA, c(nSites, nSecondary, nPrimary))
for(i in 1:nSites) {
    for(k in 1:nPrimary) {
        y[i,,k] <- rbinom(nSecondary, size=1, prob=z[i,k]*p)
    } }


## ----sim-nocov-dat,size='scriptsize'------------------------------------------
y[1:15,,1]


## ----sim-nocov-ss1,size='scriptsize'------------------------------------------
siteDets <- rowSums(y) # Dets at each site
table(siteDets)        # Frequency


## ----sim-nocov-ss2,size='scriptsize'------------------------------------------
yearDets <- apply(y, 3, sum)
yearDets

## ----un,include=FALSE---------------------------------------------------------
library(unmarked)


## ----un-umf,size='scriptsize'-------------------------------------------------
## Format as nSites by (nSec*nPrimary) matrix
y.wide <- matrix(y, nrow=nrow(y)) 
umf <- unmarkedMultFrame(y=y.wide, numPrimary=nPrimary)


## ----wfac,size='scriptsize'---------------------------------------------------
summary(umf)


## ----un-fit,size='tiny'-------------------------------------------------------
fm <- colext(~1,~1,~1,~1, umf)    
fm


## ----un-compare,size='tiny'---------------------------------------------------
## backTransform(fm, type="psi")
coef(backTransform(fm, type="psi"))
coef(backTransform(fm, type="col"))
coef(backTransform(fm, type="ext"))
coef(backTransform(fm, type="det"))


## ----un-actual,size='tiny'----------------------------------------------------
#c(initial=psi[1,1],col=gamma,ext=epsilon,det=p)
psi[1,1]
gamma
epsilon
p


## ----re,size='scriptsize',fig.height=5,out.width='65%',fig.align='center',echo=-1----
par(mai=c(0.9, 0.9, 0.1, 0.1))  
re <- ranef(fm)
occupied.post <- predict(re, func=colSums, nsim=1000)
plot(1:nPrimary, rowMeans(occupied.post), type="b", pch=16,
     xlab="Time", ylab="Sites occupied", ylim=c(0, 70))
lines(1:nPrimary, colSums(z), type="b", col=4)
segments(1:nPrimary, apply(occupied.post, 1, quantile, prob=0.025),
         1:nPrimary, apply(occupied.post, 1, quantile, prob=0.975))
legend(1, 70, c("Actual", "Estimated"), pch=c(1,16), col=c(4,1), lty=1)


## ----bugs0,size='tiny',results='hide'-----------------------------------------
writeLines(readLines("dynocc-model.jag"))

## ----bugs,size='tiny',echo=FALSE,comment='',background='lightblue'------------
writeLines(readLines("dynocc-model.jag"))


## ----bugs-data,size='small'---------------------------------------------------
jags.data <- list(y=y, nSites=nSites,
                  J=nSecondary, K=nPrimary)


## ----bugs-inits,size='small'--------------------------------------------------
jags.inits <- function() {
    list(psi1=runif(1), epsilon=runif(1),
         gamma=runif(1), p=runif(1),
         z=matrix(1, nSites, nPrimary))
}


## ----bugs-pars,size='small'---------------------------------------------------
jags.pars <- c("psi1", "epsilon", "gamma", "p", "N")


## ----jagsUI,include=FALSE,results='hide'--------------------------------------
library(jagsUI)
library(coda)


## ----bugs-mcmc,size='scriptsize',message=FALSE,cache=TRUE,results='hide'------
library(jagsUI)
jags.post.samples <- jags.basic(data=jags.data, inits=jags.inits,
                                parameters.to.save=jags.pars,
                                model.file="dynocc-model.jag",
                                n.chains=3, n.adapt=100, n.burnin=0,
                                n.iter=2000, parallel=TRUE)


## ----bugs-sum,size='tiny'-----------------------------------------------------
summary(jags.post.samples[,jags.pars[1:4]])


## ----bugs-plot1,size='footnotesize',out.width="0.7\\textwidth",fig.align='center'----
plot(jags.post.samples[,jags.pars[1:4]])


## ----sim-init-cov,size='tiny',echo=-1-----------------------------------------
set.seed(5098)
beta0.psi <- -1; beta1.psi <- 1
elevation <- rnorm(nSites) ## site covariate
z2 <- psi2 <- matrix(NA, nSites, nPrimary)
psi2[,1] <- plogis(beta0.psi + beta1.psi*elevation)
z2[,1] <- rbinom(n=nSites, size=1, prob=psi2[,1])


## ----nsim-dy-cov,size='tiny'--------------------------------------------------
epsilon2 <- 0.3 ## Local extinction prob
temperature <- matrix(rnorm(nSites*nPrimary)*elevation, ## Primary period covariate
                      nrow=nSites) 
beta0.gamma <- -2; beta1.gamma <- -1
gamma2 <- plogis(beta0.gamma + beta1.gamma*temperature)
for(k in 2:nPrimary) {
    psi2[,k] <- z2[,k-1]*(1-epsilon2) + (1-z2[,k-1])*gamma2[,k-1]
    z2[,k] <- rbinom(n=nSites, size=1, prob=psi2[,k])
}


## ----sim-cov3-cov,size='tiny'-------------------------------------------------
p2 <- 0.2
y2 <- array(NA, c(nSites, nSecondary, nPrimary))
for(i in 1:nSites) {
    for(k in 1:nPrimary) {
        y2[i,,k] <- rbinom(nSecondary, size=1, prob=z2[i,k]*p2)
    } }


## ----un-umf-cov,size='tiny'---------------------------------------------------
y2.wide <- matrix(y2, nrow=nrow(y2)) ## nSites by (nSec*nPrimary) matrix
umf2 <- unmarkedMultFrame(y=y2.wide, numPrimary=nPrimary, siteCovs=data.frame(elevation),
                          yearlySiteCovs=list(temp=temperature))


## ----wfac-cov,size='tiny'-----------------------------------------------------
summary(umf2)


## ----un-fit-cov,size='tiny'---------------------------------------------------
fm2 <- colext(~elevation,~temp,~1,~1, umf2)    
fm2


## ----preddat,size='footnotesize'----------------------------------------------
pred.data <- data.frame(temp=seq(from=-5, to=4, length=50))


## ----pred-gamma,size='footnotesize'-------------------------------------------
gamma.pred <- predict(fm2, newdata=pred.data,
                      type='col', append=TRUE)


## ----pred-gamma1,fig.height=5.5,size='tiny',out.width='80%',fig.align='center'----
plot(Predicted ~ temp, gamma.pred, type="l", ylab="Occurrence probability",
     xlab="Standardized temperature", ylim=0:1, lwd=2) 
lines(lower ~ temp, gamma.pred, col=gray(0.6))
lines(upper ~ temp, gamma.pred, col=gray(0.6))


## ----bugs-cov0,size='tiny',results='hide'-------------------------------------
writeLines(readLines("dynocc-model-covars.jag"))

## ----bugs-cov,size='tiny',echo=FALSE,comment='',background='lightblue'--------
writeLines(readLines("dynocc-model-covars.jag"))


## ----bugs-data-cov,size='scriptsize'------------------------------------------
jags.data.covs <- list(y=y2, nSites=nSites, elevation=elevation,
                       temp=temperature, J=nSecondary, K=nPrimary)


## ----bugs-inits-cov,size='scriptsize'-----------------------------------------
jags.inits.covs <- function() {
    list(beta0.psi=rnorm(1), beta1.psi=rnorm(1), epsilon=runif(1),
         beta0.gamma=rnorm(1), beta1.gamma=rnorm(1), p=runif(1),
         z=matrix(1, nSites, nPrimary))
}


## ----bugs-pars-cov,size='scriptsize'------------------------------------------
jags.pars.covs <- c("beta0.psi", "beta1.psi", "epsilon",
                    "beta0.gamma", "beta1.gamma", "p", "N")


## ----bugs-mcmc-cov,size='scriptsize',message=FALSE,cache=TRUE,results='hide'----
library(jagsUI)
jags.post.samples.covs <- jags.basic(
    data=jags.data.covs, inits=jags.inits.covs,
    parameters.to.save=jags.pars.covs,
    model.file="dynocc-model-covars.jag",
    n.chains=3, n.adapt=100, n.burnin=0, n.iter=2000, parallel=TRUE)


## ----bugs-sum-cov,size='tiny'-------------------------------------------------
summary(jags.post.samples.covs[,jags.pars.covs[1:6]])


## ----bugs-plot1-cov,size='footnotesize',out.width="0.7\\textwidth",fig.align='center'----
plot(jags.post.samples.covs[,jags.pars.covs[1:3]])


## ----bugs-plot2-cov,size='footnotesize',out.width="0.7\\textwidth",fig.align='center'----
plot(jags.post.samples.covs[,jags.pars.covs[4:6]])


## ----N-post-covs,size='small'-------------------------------------------------
library(coda)
N.names <- grep("N\\[", varnames(jags.post.samples.covs)  )
N.post <- as.matrix(jags.post.samples.covs[,N.names])
N.post.mean <- colMeans(N.post)
N.post.lower <- apply(N.post, 2, quantile, prob=0.025)
N.post.upper <- apply(N.post, 2, quantile, prob=0.975)


## ----N-post-covs-plot,size='scriptsize',fig.height=5,out.width="75%",fig.align="center",echo=-1----
par(mai=c(0.9, 0.9, 0.1, 0.1))  
plot(1:nPrimary, N.post.mean, type="b", ylim=c(0, 70), pch=16,
     xlab="Time", ylab="Sites occupied")
segments(1:nPrimary, N.post.lower, 1:nPrimary, N.post.upper)
lines(1:nPrimary, colSums(z2), type="b", col=4)
legend(1, 70, c("Actual", "Estimated"), pch=c(1,16), col=c(4,1), lty=1)


## ----beta-post,size='footnotesize'--------------------------------------------
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


## ----gamma-post-covs-plot,size='scriptsize',fig.height=5,out.width="75%",fig.align="center",echo=-1----
par(mai=c(0.9, 0.9, 0.1, 0.1))  
matplot(temperature.seq, t(gamma.post[seq(1,n.samples,10),]), type="l",
        xlab="Temperature (standardized)", ylab="Colonization prob", col=gray(0.9))
lines(temperature.seq, gamma.post.mean, lwd=3, col="blue")
lines(temperature.seq, gamma.post.lower, col="blue")
lines(temperature.seq, gamma.post.upper, col="blue")

