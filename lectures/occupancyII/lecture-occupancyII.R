## ----buildit,include=FALSE,eval=FALSE-----------------------------------------
## rnw2pdf("lecture-occupancyII")
## rnw2pdf("lecture-occupancyII", tangle=TRUE)




## ----sim-cov1,size='scriptsize'-----------------------------------------------
nSites <- 100; nVisits <- 4; set.seed(83) ## Make it reproducible
x1 <- rnorm(nSites,0,0.5); x2 <- rnorm(nSites,100,10) ## Continuous covs
w <- matrix(sample(c("Cold", "Hot"), size=nSites*nVisits, replace=T),
            nrow=nSites, ncol=nVisits)
wHot <- ifelse(w=="Hot", 1, 0)              ## Dummy variable


## ----nsim-cov2,size='scriptsize'----------------------------------------------
beta0 <- 0; beta1 <- -1; beta2 <- 0
psi <- plogis(beta0 + beta1*x1 + beta2*x2)
alpha0 <- -0.5; alpha1 <- 1; alpha2 <- 2
p <- plogis(alpha0 + alpha1*x1 + alpha2*wHot)


## ----sim-cov3,size='scriptsize'-----------------------------------------------
z <- rbinom(nSites, size=1, psi)            ## pres/absence
y <- matrix(NA, nrow=nSites, ncol=nVisits)
for(i in 1:nSites) {
    y[i,] <- rbinom(nVisits, size=1, prob=z[i]*p[i,])
}


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


## ----un-umf,size='tiny',warning=FALSE-----------------------------------------
umf <- unmarkedFrameOccu(y=y, siteCovs=data.frame(x1,x2), obsCovs=list(w=w))


## ----wfac,size='tiny'---------------------------------------------------------
summary(umf)


## ----umf-zcov,size='footnotesize'---------------------------------------------
mean.x1 <- mean(x1); mean.x2 <- mean(x2)
sd.x1 <- sd(x1); sd.x2 <- sd(x2)
siteCovs(umf)$x1s <- (x1-mean.x1)/sd.x1
siteCovs(umf)$x2s <- (x2-mean.x2)/sd.x2


## ----un-fit,size='tiny'-------------------------------------------------------
fm <- occu(~x1s+w ~x1s+x2s, umf)    ## Notice standardized covariates
fm


## ----un-compare,size='tiny'---------------------------------------------------
rbind('Occupancy coefs (beta)'=c(beta0, beta1, beta2),
      'Detection coefs (alpha)'=c(alpha0, alpha1, alpha2))


## ----preddat,size='footnotesize'----------------------------------------------
pred.data <- data.frame(x1s=seq(from=-3, to=3, length=50),
                        x2s=0, w='Hot')
pred.data$x1 <- pred.data$x1s*sd.x1+mean.x1 ## Back-transform x1s
pred.data$x2 <- pred.data$x2s*sd.x2+mean.x2 ## Back-transform x2s


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
plot(Predicted ~ x1s, psi.pred, type="l", ylab="Occurrence probability", col="blue",
     xlab="Standardized covariate (x1s)", ylim=0:1) 
lines(lower ~ x1s, psi.pred, col="grey"); lines(upper ~ x1s, psi.pred, col="grey")


## ----pred-psi1s,fig.width=7,fig.height=5,size='tiny',out.width='80%',fig.align='center',echo=-1----
par(mai=c(0.9,0.9,0.1,0.1))  
plot(Predicted ~ x1, psi.pred, type="l", ylab="Occurrence probability", col="blue",
     xlab="Original scale covariate (x1)", ylim=0:1)
lines(lower ~ x1, psi.pred, col="grey"); lines(upper ~ x1, psi.pred, col="grey")


## ----pred-plot2,fig.width=7,fig.height=5,size='tiny',out.width='70%',fig.align='center',echo=-1----
par(mai=c(0.9,0.9,0.1,0.1))  
plot(Predicted ~ x1, psi.pred, type="l", ylab="Probability", col="blue", ylim=0:1,
     xlab="x1", lwd=2)
lines(Predicted ~ x1, p.pred, col="orange", lwd=2)
legend(-3, 0.3, c("Occupancy (psi)", "Detection (p)"), lty=c(1, 1), col=c("blue", "orange"), lwd=2)


## ----bugs,size='scriptsize',echo=FALSE,background='lightblue',comment=''------
writeLines(readLines("occupancy-model-covs.jag"))


## ----bugs-data,size='small'---------------------------------------------------
jags.data <- list(y=y, x1=(x1-mean(x1))/sd(x1),
                  x2=(x2-mean(x2))/sd(x2), wHot=wHot,
                  nSites=nSites, nOccasions=nVisits)


## ----bugs-inits,size='small'--------------------------------------------------
jags.inits <- function() {
    list(beta0=rnorm(1), alpha0=rnorm(1), z=rep(1, nSites))
}


## ----bugs-pars,size='small'---------------------------------------------------
jags.pars <- c("beta0", "beta1", "beta2",
               "alpha0", "alpha1", "alpha2", "sitesOccupied")


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
                                psi.coef.post[i,"beta1"]*pred.data$x1s)
}


## ----psi-pred1,size='scriptsize',fig.align='center',out.width='80%',fig.height=5,dev='png',dpi=200----
plot(pred.data$x1, psi.post.pred[1,], type="l", xlab="x1s",
     ylab="Occurrence probability", ylim=c(0, 1), col=gray(0.8))


## ----psi-pred-post,size='scriptsize',fig.align='center',out.width='80%',fig.height=5,echo=-1,dev='png',cache=FALSE,dpi=200----
plot(pred.data$x1, psi.post.pred[1,], type="l", xlab="x1s",
     ylab="Occurrence probability", ylim=c(0, 1), col=gray(0.8))
for(i in 1:n.iter) {
    lines(pred.data$x1, psi.post.pred[i,], col=gray(0.8))
}


## ----psi-pred-post-meanCI,size='tiny',fig.align='center',out.width='80%',fig.height=5,echo=-(1:2),dev='png',cache=FALSE,dpi=200----
plot(pred.data$x1, psi.post.pred[1,], type="l", xlab="x1s",
     ylab="Occurrence probability", ylim=c(0, 1), col=gray(0.8))
for(i in 1:n.iter) {
    lines(pred.data$x1, psi.post.pred[i,], col=gray(0.8))
}
pred.post.mean <- colMeans(psi.post.pred)
pred.post.lower <- apply(psi.post.pred, 2, quantile, prob=0.025)
pred.post.upper <- apply(psi.post.pred, 2, quantile, prob=0.975)
lines(pred.data$x1, pred.post.mean, col="blue")
lines(pred.data$x1, pred.post.lower, col="blue", lty=2)
lines(pred.data$x1, pred.post.upper, col="blue", lty=2)

