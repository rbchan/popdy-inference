## ----buildit,include=FALSE,eval=FALSE-----------------------------------------
# rnw2pdf("lecture-occupancyII")
# rnw2pdf("lecture-occupancyII", tangle=TRUE)




## ----sim-cov1,size='scriptsize'-----------------------------------------------
nSites <- 200; nVisits <- 4; set.seed(83) # Make it reproducible
elev <- rnorm(nSites,5000,500) # Site covariate (continuous)
forest <- rnorm(nSites,100,10) # Another continuous site covariate
temp <- matrix(sample(c("Cold", "Hot"), size=nSites*nVisits,
                      replace=TRUE),
               nrow=nSites, ncol=nVisits) # Categorical obs covar
tempHot <- ifelse(temp=="Hot", 1, 0)      # Dummy variable


## ----nsim-cov2,size='scriptsize'----------------------------------------------
beta0 <- 5; beta1 <- -0.001; beta2 <- 0
psi <- plogis(beta0 + beta1*elev + beta2*forest)
alpha0 <- -6.5; alpha1 <- 0.001; alpha2 <- 2
p <- plogis(alpha0 + alpha1*elev + alpha2*tempHot)


## ----sim-cov3,size='scriptsize'-----------------------------------------------
z <- rbinom(nSites, size=1, psi)            ## pres/absence
y <- matrix(NA, nrow=nSites, ncol=nVisits)
for(i in 1:nSites) {
    y[i,] <- rbinom(nVisits, size=1, prob=z[i]*p[i,]) }


## ----sim-sitch,fig.width=7,fig.height=5,size='tiny',out.width='90%',fig.align='center',echo=FALSE----
par(mai=c(0.9,0.9,0.1,0.1))  
plot(function(x) plogis(beta0+beta1*x), from=min(elev), to=max(elev),
     ylab="Probability", col="blue", ylim=0:1, xlab="Elevation", lwd=2)
plot(function(x) plogis(alpha0+alpha1*x), from=min(elev), to=max(elev),
     col="orange", add=TRUE)
legend(5000, 1, c("Occupancy (psi)", "Detection (p)"), lty=c(1, 1),
       col=c("blue", "orange"), lwd=2, cex=1.2)


## ----sim-nocov-dat,size='scriptsize'------------------------------------------
y[1:20,]


## ----sim-nocov-ss1,size='scriptsize'------------------------------------------
siteDets <- rowSums(y) # Dets at each site
table(siteDets)        # Frequency


## ----sim-nocov-ss2,size='scriptsize'------------------------------------------
naiveOccupancy <- sum(siteDets>0)/nSites
naiveOccupancy 

## ----un,include=FALSE---------------------------------------------------------
library(unmarked)


## ----sim-z,size='scriptsize'--------------------------------------------------
sum(z) / nSites


## ----un-umf,size='tiny',warning=FALSE-----------------------------------------
umf <- unmarkedFrameOccu(y=y, siteCovs=data.frame(elev,forest), obsCovs=list(temp=temp))


## ----wfac,size='tiny'---------------------------------------------------------
summary(umf)


## ----umf-zcov,size='footnotesize'---------------------------------------------
mean.elev <- mean(elev); mean.forest <- mean(forest)
sd.elev <- sd(elev); sd.forest <- sd(forest)
elev_s <- (elev-mean.elev)/sd.elev
forest_s <- (forest-mean.forest)/sd.forest
siteCovs(umf)$elev_s <- elev_s
siteCovs(umf)$forest_s <- forest_s


## ----un-fit,size='scriptsize'-------------------------------------------------
fm <- occu(~elev_s+temp ~elev_s+forest_s, umf) # Standardized covariates
fm


## ----preddat,size='footnotesize'----------------------------------------------
pred.data <- data.frame(elev_s=seq(from=-3, to=3, length=50),
                        forest_s=0, temp='Hot') # Standardized
pred.data$elev <- pred.data$elev_s*sd.elev+mean.elev # Orig scale
pred.data$forest <- pred.data$forest_s*sd.forest+mean.forest 


## ----predpsi,size='footnotesize'----------------------------------------------
psi.pred <- predict(fm, newdata=pred.data,
                    type='state', append=TRUE)


## ----predp,size='footnotesize'------------------------------------------------
p.pred <- predict(fm, newdata=pred.data,
                  type='det', append=TRUE)


## ----psi-head,size='footnotesize'---------------------------------------------
print(head(psi.pred), digits=2)


## ----p-head,size='footnotesize'-----------------------------------------------
print(head(p.pred), digits=2)


## ----pred-psi1,fig.width=7,fig.height=5,size='tiny',out.width='80%',fig.align='center',echo=-1----
par(mai=c(0.9,0.9,0.1,0.1))  
plot(Predicted ~ elev_s, psi.pred, type="l", ylab="Occurrence probability", col="blue",
     xlab="Elevation (standardized)", ylim=0:1) 
lines(lower ~ elev_s, psi.pred, col="grey"); lines(upper ~ elev_s, psi.pred, col="grey")


## ----pred-psi1s,fig.width=7,fig.height=5,size='tiny',out.width='80%',fig.align='center',echo=-1----
par(mai=c(0.9,0.9,0.1,0.1))  
plot(Predicted ~ elev, psi.pred, type="l", ylab="Occurrence probability", col="blue",
     xlab="Elevation (original scale)", ylim=0:1)
lines(lower ~ elev, psi.pred, col="grey"); lines(upper ~ elev, psi.pred, col="grey")


## ----psi-actual-est,fig.width=7,fig.height=5,size='tiny',out.width='90%',fig.align='center',echo=FALSE----
par(mai=c(0.9,0.9,0.1,0.1))  
plot(Predicted ~ elev, psi.pred, type="l", ylab="Occurrence probability", col="blue",
     xlab="Elevation (original scale)", ylim=0:1)
lines(lower ~ elev, psi.pred, col="grey"); lines(upper ~ elev, psi.pred, col="grey")
plot(function(x) plogis(beta0 + beta1*x),
     from=min(psi.pred$elev), to=max(psi.pred$elev), add=TRUE,
     lwd=2)
legend(5800, 1, c("Actual", "Estimated"), lwd=c(2,1), col=c("black", "blue"))


## ----bugs,size='scriptsize',echo=FALSE,background='beige',comment=''----------
writeLines(readLines("occupancy-model-covs.jag"))


## ----bugs-data,size='small'---------------------------------------------------
jags.data <- list(y=y, elev=elev_s, forest=forest_s,
                  tempHot=tempHot,
                  nSites=nSites, nOccasions=nVisits)


## ----bugs-inits,size='small'--------------------------------------------------
jags.inits <- function() {
    list(beta0=rnorm(1), alpha0=rnorm(1), z=rep(1, nSites))
}


## ----bugs-pars,size='small'---------------------------------------------------
jags.pars <- c("beta0", "beta1", "beta2",
               "alpha0", "alpha1", "alpha2",
               "sitesOccupied")


## ----bugs-mcmc,size='scriptsize',message=FALSE,cache=FALSE--------------------
library(jagsUI)
jags.post.samples <- jags.basic(data=jags.data, inits=jags.inits,
                                parameters.to.save=jags.pars,
                                model.file="occupancy-model-covs.jag",
                                n.chains=3, n.adapt=100, n.burnin=0,
                                n.iter=2000, parallel=TRUE)


## ----bugs-sum,size='tiny'-----------------------------------------------------
summary(jags.post.samples)


## ----bugs-plot1,size='footnotesize',out.width="0.7\\textwidth",fig.align='center'----
plot(jags.post.samples[,1:3])


## ----bugs-plot2,size='footnotesize',out.width="0.7\\textwidth",fig.align='center'----
plot(jags.post.samples[,c(4:6,8)])


## ----psi-coefs,size='scriptsize'----------------------------------------------
psi.coef.post <- as.matrix(jags.post.samples[,c("beta0","beta1","beta2")])
head(psi.coef.post, n=4)


## ----psi-predmat,size='scriptsize'--------------------------------------------
n.iter <- nrow(psi.coef.post)  
psi.post.pred <- matrix(NA, nrow=n.iter, ncol=nrow(pred.data))


## ----psi-pred-bayes,size='scriptsize'-----------------------------------------
for(i in 1:n.iter) {
    psi.post.pred[i,] <- plogis(psi.coef.post[i,"beta0"] +
                                psi.coef.post[i,"beta1"]*pred.data$elev_s)
}


## ----psi-pred1,size='scriptsize',fig.align='center',out.width='80%',fig.height=5,dev='png',dpi=200----
plot(pred.data$elev, psi.post.pred[1,], type="l", xlab="Elevation",
     ylab="Occurrence probability", ylim=c(0, 1), col=gray(0.8))


## ----psi-pred-post,size='scriptsize',fig.align='center',out.width='80%',fig.height=5,echo=-1,dev='png',cache=FALSE,dpi=200----
plot(pred.data$elev, psi.post.pred[1,], type="l", xlab="elev",
     ylab="Occurrence probability", ylim=c(0, 1), col=gray(0.8))
for(i in 1:n.iter) {
    lines(pred.data$elev, psi.post.pred[i,], col=gray(0.8))
}


## ----psi-pred-post-meanCI,size='tiny',fig.align='center',out.width='80%',fig.height=5,echo=-(1:2),dev='png',cache=FALSE,dpi=200----
plot(pred.data$elev, psi.post.pred[1,], type="l", xlab="elev",
     ylab="Occurrence probability", ylim=c(0, 1), col=gray(0.8))
for(i in 1:n.iter) {
    lines(pred.data$elev, psi.post.pred[i,], col=gray(0.8))
}
pred.post.mean <- colMeans(psi.post.pred)
pred.post.lower <- apply(psi.post.pred, 2, quantile, prob=0.025)
pred.post.upper <- apply(psi.post.pred, 2, quantile, prob=0.975)
lines(pred.data$elev, pred.post.mean, col="blue")
lines(pred.data$elev, pred.post.lower, col="blue", lty=2)
lines(pred.data$elev, pred.post.upper, col="blue", lty=2)

