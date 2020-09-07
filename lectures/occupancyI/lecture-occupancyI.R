## ----buildit,include=FALSE,eval=FALSE-----------------------------------------
## rnw2pdf("lecture-occupancyI")
## rnw2pdf("lecture-occupancyI", tangle=TRUE)




## ----sim-nocov1,size='footnotesize'-------------------------------------------
psi <- 0.5       ## Occurrence probability
p <- 0.2         ## Detection probability
nSites <- 20
nVisits <- 4


## ----sim-nocov2,size='footnotesize'-------------------------------------------
set.seed(3439)    ## Just to make it reproducible
z <- rbinom(nSites, size=1, psi) ## pres/absence


## ----sim-nocov3,size='footnotesize'-------------------------------------------
y <- matrix(NA, nrow=nSites, ncol=nVisits)
for(i in 1:nSites) {
    y[i,] <- rbinom(nVisits, size=1, prob=z[i]*p)
}


## ----sim-nocov-dat,size='scriptsize'------------------------------------------
y


## ----sim-nocov-ss1,size='scriptsize'------------------------------------------
siteDets <- rowSums(y) # Dets at each site
table(siteDets)        # Frequency


## ----sim-nocov-ss2,size='scriptsize'------------------------------------------
naiveOccupancy <- sum(siteDets>0)/nSites
naiveOccupancy 


## ----un-install---------------------------------------------------------------
## install.packages("unmarked")  
library(unmarked)


## ----un-help,eval=FALSE-------------------------------------------------------
## help("unmarked")


## ----un-vig-------------------------------------------------------------------
vignette(package="unmarked")
## vignette("spp-dist")


## ----un-umf,size='small'------------------------------------------------------
umf <- unmarkedFrameOccu(y=y)
summary(umf)


## ----un-fit,size='scriptsize'-------------------------------------------------
fm <- occu(~1 ~1, umf)
summary(fm)


## ----un-back1,size='footnotesize'---------------------------------------------
backTransform(fm, type="state")


## ----un-back2,size='footnotesize'---------------------------------------------
backTransform(fm, type="det")


## ----bup,size='scriptsize'----------------------------------------------------
## Posterior probs Pr(z_i=1 | y_i)       
z.post <- ranef(fm)
## Extract posterior means
psi.conditional <-  bup(z.post, stat="mean") 


## ----psi-zpm,size='scriptsize'------------------------------------------------
round(data.frame(y=y, psi.unconditional=predict(fm, type="state")[,1],
                 psi.conditional), 3)


## ----bup-hist,size='tiny',fig.width=7,fig.height=5,out.width='0.7\\textwidth',fig.align='center'----
nsim <- 1000
sites.occupied.post <- predict(z.post, func=function(x) sum(x), nsim=nsim)
par(mai=c(0.9,0.9,0.1,0.1))
plot(table(sites.occupied.post)/nsim, lwd=5, xlab="Sites occupied",
    ylab="Empirical Bayes posterior probability")
abline(v=sum(z), col="red") ## Actual number occupied


## ----bugs,size='scriptsize'---------------------------------------------------
writeLines(readLines("occupancy-model.jag"))


## ----bugs-data,size='small'---------------------------------------------------
jags.data <- list(y=y, nSites=nSites, nOccasions=nVisits)


## ----bugs-inits,size='small'--------------------------------------------------
jags.inits <- function() {
    list(psi=runif(1), p=runif(1), z=rep(1, nSites))
}


## ----bugs-pars,size='small'---------------------------------------------------
jags.pars <- c("psi", "p", "sitesOccupied")


## ----bugs-mcmc,size='scriptsize',message=FALSE,cache=FALSE--------------------
## install.packages("jagsUI")
library(jagsUI)
jags.post.samples <- jags.basic(data=jags.data, inits=jags.inits,
                                parameters.to.save=jags.pars,
                                model.file="occupancy-model.jag",
                                n.chains=3, n.adapt=100, n.burnin=0,
                                n.iter=2000, parallel=TRUE)


## ----bugs-sum,size='scriptsize'-----------------------------------------------
summary(jags.post.samples)


## ----bugs-plot,size='footnotesize',out.width="0.7\\textwidth",fig.align='center'----
plot(jags.post.samples)

