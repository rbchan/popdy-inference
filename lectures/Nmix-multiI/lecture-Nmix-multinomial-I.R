## ----buildit,include=FALSE,eval=FALSE-----------------------------------------
## rnw2pdf("lecture-Nmix-multinomial-I")
## rnw2pdf("lecture-Nmix-multinomial-I", tangle=TRUE)




## ----multi-sim-concept,size='footnotesize'------------------------------------
N <- 20
pi <- c(survived=0.5, depredated=0.3, starved=0.2)
t(rmultinom(n=1, size=N, prob=pi))


## ----include=FALSE,echo=FALSE-------------------------------------------------
set.seed(34889243)


## ----sim-rem-nocov1,size='scriptsize',echo=-1---------------------------------
set.seed(430)
nSites <- 100
lambda1 <- 4.6  ## Expected value of N
N1 <- rpois(n=nSites, lambda=lambda1)


## ----sim-nocov2,size='scriptsize'---------------------------------------------
nPasses <- 3
K <- nPasses+1  # multinomial cells
p1 <- 0.4
pi1 <- c(p1, (1-p1)*p1, (1-p1)*(1-p1)*p1, (1-p1)^3)
y1.all <- matrix(NA, nrow=nSites, ncol=K)
for(i in 1:nSites) {
    y1.all[i,] <- rmultinom(n=1, size=N1[i], prob=pi1)    }


## ----N1y1,size='scriptsize'---------------------------------------------------
y1 <- y1.all[,-K]
head(y1, n=3)


## ----y1-kable,size='normalsize',align='center',echo=FALSE---------------------
colnames(y1) <- c("Pass 1", "Pass 2", "Pass 3")
rownames(y1) <- paste("Site", 1:nrow(y1))
kable(y1[1:10,], format="latex", booktabs=TRUE, table.envir="table")


## ----sim-cov1,size='scriptsize'-----------------------------------------------
streamDepth <- rnorm(nSites)


## ----nsim-cov2,size='scriptsize'----------------------------------------------
beta0 <- 1; beta1 <- 0.5
lambda2 <- exp(beta0 + beta1*streamDepth)
alpha0 <- 0; alpha1 <- -1
p2 <- plogis(alpha0 + alpha1*streamDepth)
pi2 <- t(sapply(p2, function(p) c(p, (1-p)*p, (1-p)^2*p, (1-p)^3)))


## ----sim-cov3,size='scriptsize'-----------------------------------------------
N2 <- rpois(nSites, lambda=lambda2)         ## local abundance 
y2.all <- matrix(NA, nrow=nSites, ncol=K)
for(i in 1:nSites) {
    y2.all[i,] <- rmultinom(n=1, size=N2[i], prob=pi2[i,])
}
y2 <- y2.all[,-K] ## Discard final column... individuals not detected


## ----sim-nocov-dat,size='scriptsize'------------------------------------------
y2[1:19,]


## ----sim-nocov-ss2,size='scriptsize'------------------------------------------
colSums(y2)


## ----sim-nocov-ss3,size='scriptsize'------------------------------------------
sum(y2)

## ----un,include=FALSE---------------------------------------------------------
library(unmarked)


## ----un-umf,size='tiny'-------------------------------------------------------
umf <- unmarkedFrameMPois(y=y2, siteCovs=data.frame(streamDepth), type="removal")


## ----wfac,size='scriptsize'---------------------------------------------------
summary(umf)


## ----un-fit,size='tiny'-------------------------------------------------------
fm <- multinomPois(~streamDepth ~streamDepth, umf)    
fm


## ----un-compare,size='tiny'---------------------------------------------------
c(beta0=beta0, beta1=beta1); c(alpha0=alpha0, alpha1=alpha1)


## ----ranef,size='scriptsize',out.width='80%',fig.align='center',fig.width=9----
re <- ranef(fm, K=15)
plot(re, layout=c(4,3), subset=site%in%1:12, xlim=c(-1, 11), lwd=5)


## ----Ntotal,size='scriptsize',out.width='60%',fig.align='center'--------------
N.total.post <- predict(re, func=sum, nsim=1000)
hist(N.total.post, freq=FALSE, main="", xlab="N total", ylab="Probability")


## ----preddat,size='footnotesize'----------------------------------------------
pred.data <- data.frame(streamDepth=seq(-3, 3, length=20))


## ----predpsi,size='footnotesize'----------------------------------------------
lambda.pred <- predict(fm, newdata=pred.data,
                       type='state', append=TRUE)


## ----psi-head,size='footnotesize'---------------------------------------------
print(head(lambda.pred), digits=2)


## ----pred-lam2,fig.width=7,fig.height=5,size='tiny',out.width='80%',fig.align='center',echo=-1----
par(mai=c(0.9,0.9,0.1,0.1))
plot(Predicted ~ streamDepth, lambda.pred, ylab="Expected value of abundance",
     ylim=c(0,20), xlab="Stream depth", type="l", lwd=2)
lines(lower ~ streamDepth, lambda.pred, col="grey")
lines(upper ~ streamDepth, lambda.pred, col="grey")
points(rowSums(y2)~streamDepth)
lines(lowess(rowSums(y2)~streamDepth), col="blue", lwd=2)  ## Loess line for fun (it's way off)


## ----bugs-removal2,size='scriptsize',echo=FALSE,comment='',background='lightblue'----
writeLines(readLines("removal-mod2.jag"))


## ----bugs-data2,size='small'--------------------------------------------------
jags.data.rem2 <- list(y=y2, n=rowSums(y2),
                       streamDepth=streamDepth,
                       nSites=nSites, nPasses=nPasses)


## ----bugs-inits,size='small'--------------------------------------------------
jags.inits.rem <- function() {
    list(lambda.intercept=runif(1), alpha0=rnorm(1),
         N=rowSums(y2)+rpois(nrow(y2), 2))
}


## ----bugs-pars,size='small'---------------------------------------------------
jags.pars.rem <- c("beta0", "beta1",
                   "alpha0", "alpha1", "totalAbundance")


## ----bugs-mcmc-rem2,size='tiny',message=FALSE,cache=TRUE,results='hide'-------
library(jagsUI); library(coda)
jags.post.rem2 <- jags.basic(data=jags.data.rem2, inits=jags.inits.rem,
                             parameters.to.save=jags.pars.rem, model.file="removal-mod2.jag",
                             n.chains=3, n.adapt=100, n.burnin=0, n.iter=2000, parallel=TRUE)


## ----bugs-sum-rem2,size='tiny'------------------------------------------------
summary(jags.post.rem2[,jags.pars.rem])


## ----bugs-plot1-rem2,size='footnotesize',out.width="0.7\\textwidth",fig.align='center'----
plot(jags.post.rem2[,jags.pars.rem[1:3]])


## ----bugs-plot2-rem2,size='footnotesize',out.width="0.7\\textwidth",fig.align='center'----
plot(jags.post.rem2[,jags.pars.rem[4:5]])


## ----bugs-removal3,size='scriptsize',echo=FALSE,comment='',background='lightblue'----
writeLines(readLines("removal-mod3.jag"))


## ----bugs-data3,size='small'--------------------------------------------------
jags.data.rem3 <- list(y=y2, streamDepth=streamDepth,
                       nSites=nSites)


## ----bugs-mcmc-rem3,size='scriptsize',message=FALSE,cache=TRUE,results='hide'----
jags.post.rem3 <- jags.basic(data=jags.data.rem3, inits=jags.inits.rem,
                             parameters.to.save=jags.pars.rem,
                             model.file="removal-mod3.jag",
                             n.chains=3, n.adapt=100, n.burnin=0,
                             n.iter=2000, parallel=TRUE)


## ----bugs-sum-rem3,size='tiny'------------------------------------------------
summary(jags.post.rem3[,jags.pars.rem])


## ----bugs-plot1-rem3,size='footnotesize',out.width="0.7\\textwidth",fig.align='center'----
plot(jags.post.rem3[,jags.pars.rem[1:3]])


## ----bugs-plot2-rem3,size='footnotesize',out.width="0.7\\textwidth",fig.align='center'----
plot(jags.post.rem3[,jags.pars.rem[4:5]])


## ----bugs-removal4,size='scriptsize',echo=FALSE,comment='',background='lightblue'----
writeLines(readLines("removal-mod4.jag"))


## ----bugs-data4,size='small'--------------------------------------------------
jags.data.rem4 <- list(y=y2, streamDepth=streamDepth,
                       nSites=nSites)


## ----bugs-mcmc-rem4,size='scriptsize',message=FALSE,cache=TRUE,results='hide'----
jags.post.rem4 <- jags.basic(data=jags.data.rem4, inits=jags.inits.rem,
                             parameters.to.save=jags.pars.rem,
                             model.file="removal-mod4.jag",
                             n.chains=3, n.adapt=100, n.burnin=0,
                             n.iter=2000, parallel=TRUE)


## ----bugs-sum-rem4,size='tiny'------------------------------------------------
summary(jags.post.rem4[,jags.pars.rem[1:4]])


## ----bugs-plot1-rem4,size='footnotesize',out.width="0.7\\textwidth",fig.align='center'----
plot(jags.post.rem4[,jags.pars.rem[1:4]])


## ----gel2,size='tiny'---------------------------------------------------------
library(coda)
gelman.diag(jags.post.rem2[,jags.pars.rem[1:4]], multivar=FALSE) ## Conditional-on-N and n


## ----gel3,size='tiny',echo=-(1:2)---------------------------------------------
#Rhat3 <- gelman.diag(jags.post.rem3)$psrf
#Rhat3[jags.pars.rem[1:4],1,drop=FALSE]
gelman.diag(jags.post.rem3[,jags.pars.rem[1:4]], multivar=FALSE) ## Sequential binomial


## ----gel4,size='tiny',echo=-(1:2)---------------------------------------------
#Rhat4 <- gelman.diag(jags.post.rem4)$psrf
#Rhat4[jags.pars.rem[1:4],1,drop=FALSE]
gelman.diag(jags.post.rem4[,jags.pars.rem[1:4]], multivar=FALSE) ## Marginalized N


## ----ess2,size='scriptsize'---------------------------------------------------
effectiveSize(jags.post.rem2[,jags.pars.rem[1:4]])


## ----ess3,size='scriptsize'---------------------------------------------------
effectiveSize(jags.post.rem3[,jags.pars.rem[1:4]])


## ----ess4,size='scriptsize'---------------------------------------------------
effectiveSize(jags.post.rem4[,jags.pars.rem[1:4]])

