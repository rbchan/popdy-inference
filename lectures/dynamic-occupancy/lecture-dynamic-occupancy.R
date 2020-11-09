## ----buildit,include=FALSE,eval=FALSE-----------------------------------------
## rnw2pdf("lecture-dynamic-occupancy")
## rnw2pdf("lecture-dynamic-occupancy", tangle=TRUE)




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
y.wide <- matrix(y, nrow=nrow(y)) ## Format as nSites by (nSec*nPrimary)
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


## ----re,size='scriptsize',fig.height=5,out.width='70%',fig.align='center',echo=-1----
par(mai=c(0.9, 0.9, 0.1, 0.1))  
re <- ranef(fm)
occupied.post <- predict(re, func=colSums, nsim=1000)
plot(1:nPrimary, rowMeans(occupied.post), type="b",
     xlab="Time", ylab="Sites occupied", ylim=c(0, 70))
segments(1:nPrimary, apply(occupied.post, 1, quantile, prob=0.025),
         1:nPrimary, apply(occupied.post, 1, quantile, prob=0.975))


## ----bugs,size='scriptsize',echo=FALSE----------------------------------------
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
set.seed(54598)
beta0.psi <- -1; beta1.psi <- 1
elevation <- rnorm(nSites)
z2 <- psi2 <- matrix(NA, nSites, nPrimary)
psi2[,1] <- plogis(beta0.psi + beta1.psi*elevation)
z2[,1] <- rbinom(n=nSites, size=1, prob=psi2[,1])


## ----nsim-dy-cov,size='scriptsize'--------------------------------------------
epsilon2 <- 0.3 ## Local extinction prob
temperature <- matrix(rnorm(nSites*nPrimary)*elevation, nrow=nSites)    
beta0.gamma <- 1; beta1.gamma <- -1
gamma2 <- plogis(beta0.gamma + beta1.gamma*temperature)
for(k in 2:nPrimary) {
    psi2[,k] <- z2[,k-1]*(1-epsilon2) + (1-z2[,k-1])*gamma2[,k]
    z2[,k] <- rbinom(n=nSites, size=1, prob=psi2[,k])    }


## ----sim-cov3-cov,size='scriptsize'-------------------------------------------
p2 <- 0.2
y2 <- array(NA, c(nSites, nSecondary, nPrimary))
for(i in 1:nSites) {
    for(k in 1:nPrimary) {
        y2[i,,k] <- rbinom(nSecondary, size=1, prob=z2[i,k]*p2)
    } }


## ----un-umf-cov,size='scriptsize'---------------------------------------------
y2.wide <- matrix(y2, nrow=nrow(y2)) ## Format as nSites by (nSec*nPrimary)
umf2 <- unmarkedMultFrame(y=y2.wide, numPrimary=nPrimary,
                          siteCovs=data.frame(elevation))##,
##                          obsCovs=list(temp=temperature[,-nPrimary]))


## ----wfac-cov,size='scriptsize'-----------------------------------------------
summary(umf2)


## ----un-fit-cov,size='tiny'---------------------------------------------------
fm2 <- colext(~elevation,~elevation,~1,~1, umf)    
fm2


## ----bugs-cov,size='scriptsize',echo=FALSE------------------------------------
writeLines(readLines("dynocc-model.jag"))


## ----bugs-data-cov,size='small'-----------------------------------------------
jags.data <- list(y=y, nSites=nSites,
                  J=nSecondary, K=nPrimary)


## ----bugs-inits-cov,size='small'----------------------------------------------
jags.inits <- function() {
    list(psi1=runif(1), epsilon=runif(1),
         gamma=runif(1), p=runif(1),
         z=matrix(1, nSites, nPrimary))
}


## ----bugs-pars-cov,size='small'-----------------------------------------------
jags.pars <- c("psi1", "epsilon", "gamma", "p", "N")


## ----bugs-mcmc-cov,size='scriptsize',message=FALSE,cache=TRUE,results='hide'----
library(jagsUI)
jags.post.samples <- jags.basic(data=jags.data, inits=jags.inits,
                                parameters.to.save=jags.pars,
                                model.file="dynocc-model.jag",
                                n.chains=3, n.adapt=100, n.burnin=0,
                                n.iter=2000, parallel=TRUE)


## ----bugs-sum-cov,size='tiny'-------------------------------------------------
summary(jags.post.samples[,jags.pars[1:4]])


## ----bugs-plot1-cov,size='footnotesize',out.width="0.7\\textwidth",fig.align='center'----
plot(jags.post.samples[,jags.pars[1:4]])

