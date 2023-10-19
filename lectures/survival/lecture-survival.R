## ----sim0,size='scriptsize',echo=-1,out.width="0.6\\textwidth",fig.align='center',fig.width=9----
set.seed(7829)
n <- 100                          ## Sample size
hazard0 <- 1/100
survival.time0 <- rexp(n, hazard0)
hist(survival.time0, xlab="Survival time (time to mortality)", main="")
abline(v=c(mean(survival.time0), 100), col=c("blue", "black"))
legend("topright", legend=c("Sample mean", "Population mean"),
       col=c("blue", "black"), lwd=2)


## ----sim1,size='scriptsize',echo=-1,out.width="0.5\\textwidth",fig.align='center'----
set.seed(7829)
n <- 200                          ## Sample size
x <- rnorm(n)                     ## Covariate
beta0 <- -5; beta1 <- -0.5        ## Hazard coefficients
hazard <- exp(beta0 + beta1*x) 
survival.time <- rexp(n, hazard)
## summary(survival.time)
hist(survival.time, xlab="Survival time (time to mortality)", main="")


## ----cen,size='footnotesize'--------------------------------------------------
ctime <- 500    ## Censoring time
censored <- ifelse(survival.time>ctime, 1, 0)
table(censored) ## Most individuals died before censoring time


## ----cen1,size='footnotesize'-------------------------------------------------
survival.time.c <- ifelse(censored, ctime, survival.time) 
summary(survival.time.c)


## ----life-lines,size='footnotesize',fig.width=8,out.width='90%',fig.align='center',echo=FALSE----
plot(0, type="n", xlim=c(0, ctime), ylim=c(1,20),
     xlab="Time", ylab="Individual")
abline(h=1:20, col=gray(0.8))
segments(rep(0, 20), 1:20, survival.time.c[1:20], 1:20)
points(rep(0,20), 1:20, pch=16, col="seagreen3", cex=1)
points(survival.time.c[1:20], 1:20, pch=17,
       col=ifelse(survival.time.c[1:20]==ctime, "blue", "red"))
legend(170, 23, c("mortality", "censoring"), pch=17, col=c("red", "blue"),
       xpd=TRUE, horiz=TRUE)


## ----km,out.width="0.6\\textwidth",fig.align='center',size='footnotesize'-----
library(survival)
y <- Surv(survival.time.c, 1-censored)
plot(y, xlab="Time", ylab="Survivorship", main="Kaplan-Meier")


## ----survreg,size='scriptsize'------------------------------------------------
summary(pm1 <- survreg(y ~ x, dist="exponential"))


## ----plothaz,out.width="0.5\\textwidth",fig.align='center',size='scriptsize'----
plot(function(x) exp(beta0 + beta1*x), from=-3, to=3,
     xlab="Covariate", ylab="Hazard", col="black", ylim=c(0, 0.035))
alpha.hat <- coef(pm1)   ## Extract the estimates
plot(function(x) 1/exp(alpha.hat[1] + alpha.hat[2]*x), ## survreg formulation
     from=-3, to=3, col="blue", add=TRUE, lty=2)
legend(0.5, 0.035, c("Estimated hazard", "Hazard"), col=c(4,1), lty=2:1)


## ----coxreg,size='scriptsize'-------------------------------------------------
(cox1 <- coxph(y ~ x))


## ----cox-S,out.width="0.65\\textwidth",fig.width=8,fig.align='center',size='scriptsize'----
plot(survfit(cox1, newdata=data.frame(x=c(-2,0,2))),
     xlab="Time", ylab="Survivorship", col=1:3,
     main="Cox proportional hazards fit")
legend(350, 1, paste("x =", c(-2, 0, 2)), lty=1, col=1:3)


## ----exp-jag,eval=FALSE-------------------------------------------------------
## writeLines(readLines("surv-exp.jag"))

## ----exp-jag2,size='footnotesize',background='lightblue',comment='',echo=FALSE----
writeLines(readLines("surv-exp.jag"))


## ----load,include=FALSE-------------------------------------------------------
library(jagsUI)
library(coda)


## ----jd-surv-exp,size='scriptsize'--------------------------------------------
survival.time.jags <- survival.time.c
survival.time.jags[censored==1] <- NA
jd.exp <- list(survivalTime=survival.time.jags, censored=censored,
               censorTime=rep(ctime, length(censored)),
               x=x, n=length(censored))


## ----ji-jp-exp,size='scriptsize'----------------------------------------------
ji.exp <- function() list(beta0=rnorm(1), beta1=rnorm(1)) 
jp.exp <- c("beta0", "beta1")


## ----js,size='scriptsize',cache=TRUE,results='hide',eval=TRUE-----------------
jags.post.samples.exp <- jags.basic(data=jd.exp, inits=ji.exp,
                                    parameters.to.save=jp.exp,
                                    model.file="surv-exp.jag",
                                    n.chains=3, n.adapt=100, n.burnin=0,
                                    n.iter=2000, parallel=TRUE)


## ----post-samps-exp,fig.align='center',out.width='65%'------------------------
plot(jags.post.samples.exp)


## ----haz-post,fig.width=8,out.width='90%',fig.align='center',echo=FALSE-------
beta.post <- as.matrix(jags.post.samples.exp)[,c("beta0","beta1")]
x.pred <- seq(-3, 3, by=0.2)
lambda.post <- matrix(NA_real_, nrow(beta.post), length(x.pred))
for(i in 1:nrow(beta.post))  lambda.post[i,] <- exp(beta.post[i,1] + beta.post[i,2]*x.pred)
matplot(x.pred, y=t(lambda.post[1:500,]), type="l", col=rgb(0,0,1,0.02),
        xlab="Covariate (x)", ylab="Hazard")
lines(x.pred, colMeans(lambda.post), type="l", col="blue", lwd=2)
plot(function(x) exp(beta0 + beta1*x), from=-3, to=3, col="orange", lwd=2, add=TRUE)
legend(1, 0.06, c("Posterior sample", "Posterior mean", "Actual"),
       col=c(rgb(0,0,1,0.2), "blue", "orange"), lty=1, lwd=c(1,2,2))


## ----sim-norobust,size='footnotesize',echo=-1---------------------------------
set.seed(34918)  
maxTime <- 10           ## Time period
n <- 25                 ## nIndividuals
phi <- 0.7              ## Survival probability over 1 time step


## ----sim-norbust2,size='footnotesize'-----------------------------------------
z <- matrix(NA, n, maxTime)
first <- rpois(n, 1)+1  ## random release dates
for(i in 1:n) {
    z[i,first[i]] <- 1  ## Known alive at release
    for(t in (first[i]+1):maxTime) {
        z[i,t] <- rbinom(1, 1, z[i,t-1]*phi) ## Alive/dead state
    }
}


## ----z,size='small'-----------------------------------------------------------
z[1:15,]


## ----jags-surv-dtime,eval=FALSE-----------------------------------------------
## writeLines(readLines("surv-dtime.jag"))

## ----jags-surv-dtime2,background='lightblue',comment='',echo=FALSE------------
writeLines(readLines("surv-dtime.jag"))


## ----jd-dtime0,size='scriptsize'----------------------------------------------
jd.dtime0 <- list(z=z, n=n, first=first, maxTime=maxTime)


## ----ji-dtime0,size='scriptsize'----------------------------------------------
ji.dtime0 <- function() list(phi=runif(1))
jp.dtime0 <- c("phi")


## ----js-dtime0,size='scriptsize',cache=TRUE,results='hide',eval=TRUE----------
jags.post.samples.dtime0 <- jags.basic(data=jd.dtime0, inits=ji.dtime0,
                                       parameters.to.save=jp.dtime0,
                                       model.file="surv-dtime.jag",
                                       n.chains=3, n.adapt=100, n.burnin=0,
                                       n.iter=2000, parallel=TRUE)


## ----post-samps-dtime0,fig.align='center',out.width='65%'---------------------
plot(jags.post.samples.dtime0)


## ----phi-post,size='tiny',fig.width=8,out.width='60%',fig.align='center',echo=-3----
phi.post <- as.matrix(jags.post.samples.dtime0)[,"phi"]
S.post <- sapply(phi.post, function(x) x^(0:10))
par(mai=c(0.9,0.9,0.1,0.1))
matplot(0:10, S.post[,1:1000], type="l", xlab="Time", ylab="Survivorship",
        col=rgb(0,0.5,1,0.02))
lines(0:10, rowMeans(S.post), col="royalblue", lwd=3)
lines(0:10, phi^(0:10), col="cadetblue", lwd=3)
legend(6, 1, c("Posterior sample", "Posterior mean", "Actual"),
       col=c(rgb(0,0.5,1,0.2), "royalblue", "cadetblue"), lty=1, lwd=c(1,3,3))


## ----jags-surv-dtime-tcovs,eval=FALSE-----------------------------------------
## writeLines(readLines("surv-dtime-tcovs.jag"))

## ----jags-surv-dtime2-tcovs,background='lightblue',comment='',echo=FALSE------
writeLines(readLines("surv-dtime-tcovs.jag"))


## ----comp-risk-sim,size='small'-----------------------------------------------
beta0 <- c(-6, -4, -3) ## log-hazard. Not covariates
lambda <- exp(beta0)   ## Cause-specific hazard
pi <- lambda / (1+sum(lambda))
pi[4] <- 1-sum(pi)     ## Probability of surviving t to t+1


## ----comp-risk-sim2,size='small'----------------------------------------------
nDeer <- 100
nDays <- 100
z <- matrix(NA, nDeer, nDays)
z[,1] <- 4  ## Everyone starts in state 4 (alive)


## ----comp-risk-Phi,size='scriptsize',echo=-4----------------------------------
Phi <- diag(4)
rownames(Phi) <- colnames(Phi) <- c("cougar", "bear", "wolf", "alive")
Phi[4,] <- pi
# kable(Phi, format="latex", digits=3, booktabs=TRUE)


## ----comp-risk-z,size='scriptsize'--------------------------------------------
for(i in 1:nDeer) {
    for(t in 2:nDays) {
        z[i,t] <- which(rmultinom(n=1, size=1, prob=Phi[z[i,t-1],])==1)
    }
}


## ----jags-surv-dtime-comp-risk,size='scriptsize',eval=FALSE-------------------
## writeLines(readLines("surv-dtime-comp-risks.jag"))


## ----jags-surv-dtime2-comp-risk,size='scriptsize',background='lightblue',comment='',echo=FALSE----
writeLines(readLines("surv-dtime-comp-risks.jag"))


## ----jd-comp-risk,size='scriptsize'-------------------------------------------
Phi.data <- diag(4)
Phi.data[4,] <- NA
jd.comp.risk <- list(z=z, n=nDeer, first=rep(1,nrow(z)), maxTime=ncol(z),
                     nRisks=3, Phi=Phi.data)


## ----ji-comp-risk,size='scriptsize'-------------------------------------------
ji.comp.risk <- function() list(beta0=log(runif(3)))
jp.comp.risk <- c("beta0", "pi")


## ----js-comp-risk,size='scriptsize',cache=TRUE,results='hide',eval=TRUE-------
jags.ps.comp.risk <- jags.basic(data=jd.comp.risk, inits=ji.comp.risk,
                                parameters.to.save=jp.comp.risk,
                                model.file="surv-dtime-comp-risks.jag",
                                n.chains=3, n.adapt=100, n.burnin=0,
                                n.iter=2000, parallel=TRUE)


## ----post-samps-comp-risk,fig.align='center',out.width='75%',echo=FALSE-------
plot(jags.ps.comp.risk[,c("beta0[1]", "beta0[2]", "beta0[3]")])


## ----post-samps-comp-risk-pi,fig.align='center',out.width='75%',echo=FALSE----
plot(jags.ps.comp.risk[,c("pi[1]", "pi[2]", "pi[3]", "pi[4]")])

