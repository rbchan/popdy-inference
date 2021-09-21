## ----buildit,include=FALSE,eval=FALSE---------------------------------------------------
## rnw2pdf("lecture-Nmix-binomial-I")
## rnw2pdf("lecture-Nmix-binomial-I", tangle=TRUE)




## ----sim-nocov1,size='scriptsize'-------------------------------------------------------
nSites <- 100
nVisits <- 4
set.seed(3439)  ## Make it reproducible
lambda1 <- 2.6  ## Expected value of N
N1 <- rpois(n=nSites, lambda=lambda1)


## ----sim-nocov2,size='scriptsize'-------------------------------------------------------
p1 <- 0.3
y1 <- matrix(NA, nrow=nSites, ncol=nVisits)
for(i in 1:nSites) {
    y1[i,] <- rbinom(nVisits, size=N1[i], prob=p1)
}


## ----N1y1,size='scriptsize'-------------------------------------------------------------
cbind(y1, N1)[1:5,]


## ----sim-cov1,size='scriptsize'---------------------------------------------------------
forest <- factor(sample(c("Hardwood", "Mixed", "Pine"), nSites, replace=TRUE))
forestMixed <- ifelse(forest=="Mixed", 1, 0)        ## Dummy
forestPine <- ifelse(forest=="Pine", 1, 0)          ## Dummy
temp <- matrix(rnorm(nSites*nVisits), nrow=nSites, ncol=nVisits)


## ----nsim-cov2,size='scriptsize'--------------------------------------------------------
beta0 <- 0; beta1 <- -1; beta2 <- 1
lambda2 <- exp(beta0 + beta1*forestMixed + beta2*forestPine)
alpha0 <- -2; alpha1 <- 1
p2 <- plogis(alpha0 + alpha1*temp)


## ----sim-cov3,size='scriptsize'---------------------------------------------------------
N2 <- rpois(nSites, lambda=lambda2)         ## local abundance 
y2 <- matrix(NA, nrow=nSites, ncol=nVisits)
for(i in 1:nSites) {
    y2[i,] <- rbinom(nVisits, size=N2[i], prob=p2[i,])
}


## ----sim-nocov-dat,size='scriptsize'----------------------------------------------------
y2[1:20,]


## ----sim-nocov-ss1,size='scriptsize'----------------------------------------------------
# Max count at each site
maxCounts <- apply(y2, 1, max) 
table(maxCounts)              


## ----sim-nocov-ss2,size='scriptsize'----------------------------------------------------
naiveOccupancy <- sum(maxCounts>0)/nSites
naiveOccupancy 

## ----un,include=FALSE-------------------------------------------------------------------
library(unmarked)


## ----un-umf,size='tiny'-----------------------------------------------------------------
umf <- unmarkedFramePCount(y=y2, siteCovs=data.frame(forest), obsCovs=list(temp=temp))


## ----wfac,size='tiny'-------------------------------------------------------------------
summary(umf)


## ----un-fit0,size='tiny'----------------------------------------------------------------
fm0 <- pcount(~1 ~1, umf, K=100)    
fm0


## ----back,size='tiny'-------------------------------------------------------------------
backTransform(fm0, type="state")


## ----un-fit,size='tiny'-----------------------------------------------------------------
fm <- pcount(~temp ~forest, umf, K=100)    
fm


## ----un-compare,size='tiny'-------------------------------------------------------------
c(beta0=beta0, beta1=beta1, beta2=beta2); c(alpha0=alpha0, alpha1=alpha1)


## ----width90,include=FALSE--------------------------------------------------------------
oo <- options(width=90)


## ----coef,size='tiny'-------------------------------------------------------------------
round(coef(fm), digits=4)


## ----upK,size='tiny'--------------------------------------------------------------------
fm.test <- pcount(~temp ~forest, umf, K=150)    
round(coef(fm.test), digits=4)


## ----ranef,size='scriptsize',out.width='80%',fig.align='center',fig.width=9-------------
re <- ranef(fm)
plot(re, layout=c(4,3), subset=site%in%1:12, xlim=c(-1, 11), lwd=5)


## ----Ntotal,size='scriptsize',out.width='40%',fig.align='center'------------------------
N.total.post <- predict(re, func=sum, nsim=1000)
hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")


## ----Ntotal-stats,size='scriptsize'-----------------------------------------------------
c(Estimate=mean(N.total.post), quantile(N.total.post, prob=0.025),
  quantile(N.total.post, prob=0.975))


## ----preddat,size='footnotesize'--------------------------------------------------------
pred.data <- data.frame(forest=c("Hardwood", "Mixed", "Pine"),
                        temp=0) 


## ----predpsi,size='footnotesize'--------------------------------------------------------
lambda.pred <- predict(fm, newdata=pred.data,
                       type='state', append=TRUE)


## ----psi-head,size='footnotesize'-------------------------------------------------------
print(head(lambda.pred), digits=2)


## ----pred-psi1,fig.width=7,fig.height=5.5,size='tiny',out.width='80%',fig.align='center'----
bpx <- barplot(lambda.pred$Predicted, ylab="Expected value of abundance", #col="blue",
               ylim=c(0,3.5), names=lambda.pred$forest, xlab="Forest type"); box()
arrows(bpx, lambda.pred$Predicted, bpx, lambda.pred$Predicted+lambda.pred$SE,
       angle=90, length=0.1)


## ----bugs0,size='scriptsize',echo=FALSE,comment='',background='lightblue'---------------
writeLines(readLines("Nmix-model.jag"))

## ----jagsUI0,include=FALSE--------------------------------------------------------------
library(jagsUI)


## ----bugs-data0,size='scriptsize'-------------------------------------------------------
jags.data0 <- list(y=y2, nSites=nSites, nOccasions=nVisits)


## ----bugs-inits0,size='scriptsize'------------------------------------------------------
jags.inits0 <- function() list(lambda=runif(1), p=runif(1), N=maxCounts)


## ----bugs-pars0,size='scriptsize',results='hide',message=FALSE--------------------------
jags0.post.samples <- jags.basic(data=jags.data0, inits=jags.inits0,
                                 parameters.to.save=c("lambda", "p"), 
                                 model.file="Nmix-model.jag",
                                 n.chains=3, n.adapt=100, n.burnin=0,
                                 n.iter=2000, parallel=TRUE)


## ----bugs,size='scriptsize',echo=FALSE,comment='',background='lightblue'----------------
writeLines(readLines("Nmix-model-covs.jag"))

## ----jagsUI,include=FALSE---------------------------------------------------------------
library(jagsUI)


## ----bugs-data,size='small'-------------------------------------------------------------
jags.data <- list(y=y2, temp=temp,
                  forestMixed=forestMixed,
                  forestPine=forestPine,
                  nSites=nSites, nOccasions=nVisits)


## ----bugs-inits,size='small'------------------------------------------------------------
jags.inits <- function() {
    list(beta0=rnorm(1), alpha0=rnorm(1), N=maxCounts)
}


## ----bugs-pars,size='small'-------------------------------------------------------------
jags.pars <- c("beta0", "beta1", "beta2",
               "alpha0", "alpha1", "totalAbundance")


## ----bugs-mcmc,size='scriptsize',message=FALSE,cache=TRUE-------------------------------
library(jagsUI)
jags.post.samples <- jags.basic(data=jags.data, inits=jags.inits,
                                ## Monitor "N[i]" 
                                parameters.to.save=c(jags.pars, "N"), 
                                model.file="Nmix-model-covs.jag",
                                n.chains=3, n.adapt=100, n.burnin=0,
                                n.iter=2000, parallel=TRUE)


## ----bugs-sum,size='tiny'---------------------------------------------------------------
summary(jags.post.samples[,jags.pars])


## ----localN,out.width='70%',fig.align='center',size='scriptsize'------------------------
plot(jags.post.samples[,paste0("N[", 1:4, "]")])


## ----bugs-plot1,size='footnotesize',out.width="0.7\\textwidth",fig.align='center'-------
plot(jags.post.samples[,jags.pars[1:3]])


## ----bugs-plot2,size='footnotesize',out.width="0.7\\textwidth",fig.align='center'-------
plot(jags.post.samples[,jags.pars[4:5]])


## ----psi-coefs,size='scriptsize'--------------------------------------------------------
p.coef.post <- as.matrix(jags.post.samples[,c("alpha0","alpha1")])
head(p.coef.post, n=4)


## ----p-predmat,size='scriptsize'--------------------------------------------------------
n.iter <- nrow(p.coef.post)
temp.pred <- seq(-3, 3, length=50)
p.post.pred <- matrix(NA, nrow=n.iter, ncol=length(temp.pred))


## ----psi-pred-bayes,size='scriptsize'---------------------------------------------------
for(i in 1:n.iter) {
    p.post.pred[i,] <- plogis(p.coef.post[i,"alpha0"] +
                              p.coef.post[i,"alpha1"]*temp.pred)
}


## ----psi-pred-post-meanCI,size='tiny',fig.align='center',out.width='70%',fig.height=5,echo=-(1),dev='png',cache=TRUE,dpi=200----
par(mai=c(0.9,0.9,0.1,0.1))  
plot(temp.pred, p.post.pred[1,], type="l", xlab="Temperature (standardized)",
     ylab="Detection probability", ylim=c(0, 1), col=gray(0.8))
for(i in seq(1, n.iter, by=10)) {  ## Thin by 10
    lines(temp.pred, p.post.pred[i,], col=gray(0.8))  }
pred.post.mean <- colMeans(p.post.pred)
pred.post.lower <- apply(p.post.pred, 2, quantile, prob=0.025)
pred.post.upper <- apply(p.post.pred, 2, quantile, prob=0.975)
lines(temp.pred, pred.post.mean, col="blue")
lines(temp.pred, pred.post.lower, col="blue", lty=2)
lines(temp.pred, pred.post.upper, col="blue", lty=2)


## ----bugs-data2,size='small'------------------------------------------------------------
jags.data2 <- jags.data
jags.data2$y[] <- NA      ## Replace data with missing values


## ----bugs-mcmc-prior,size='scriptsize',message=FALSE,cache=TRUE-------------------------
jags.prior.samples <- jags.basic(data=jags.data2, inits=jags.inits,
                                 parameters.to.save=jags.pars,
                                 model.file="Nmix-model-covs.jag",
                                 n.chains=3, n.adapt=100, n.burnin=0,
                                 n.iter=2000, parallel=TRUE)


## ----summary-prior,size='tiny'----------------------------------------------------------
summary(jags.prior.samples)


## ----psi-coefs-prior,size='scriptsize'--------------------------------------------------
p.coef.prior <- as.matrix(jags.prior.samples[,c("alpha0","alpha1")])
p.prior.pred <- matrix(NA, nrow=n.iter, ncol=length(temp.pred))
for(i in 1:n.iter) {
    p.prior.pred[i,] <- plogis(p.coef.prior[i,"alpha0"] +
                               p.coef.prior[i,"alpha1"]*temp.pred)
}


## ----prior-post-pred-1,size='scriptsize'------------------------------------------------
pred.prior.mean <- colMeans(p.prior.pred)
pred.prior.lower <- apply(p.prior.pred, 2, quantile, prob=0.025)
pred.prior.upper <- apply(p.prior.pred, 2, quantile, prob=0.975)


## ----prior-post-pred-2,size='scriptsize',fig.height=6,dev='png',dpi=200,fig.show='hide'----
plot(temp.pred, p.post.pred[1,], type="n", xlab="Temperature (standardized)",
     ylab="Detection probability", ylim=c(0, 1.3), col=gray(0.8))
polygon(x=c(temp.pred, rev(temp.pred)),
        y=c(pred.post.lower, rev(pred.post.upper)),
        col=rgb(0,0,1,0.5), border=NA)                   # Post CI
polygon(x=c(temp.pred, rev(temp.pred)),
        y=c(pred.prior.lower, rev(pred.prior.upper)),
        col=rgb(0,1,0,0.5), border=NA)                   # Prior CI
for(i in seq(1, n.iter, by=100)) {  ## Thin by 100
    lines(temp.pred, p.prior.pred[i,], col=gray(0.8))  } # Prior preds
lines(temp.pred, pred.post.mean, col="blue", lwd=2)
lines(temp.pred, pred.prior.mean, col="darkgreen", lwd=2)
legend(-3, 1.3, c("Prior mean", "Prior samples", "Posterior mean"),
       col=c("darkgreen", "grey", "blue"), lwd=2)

