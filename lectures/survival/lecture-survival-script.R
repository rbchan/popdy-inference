#' ---
#' format:
#'   revealjs:
#'     theme: default
#'     pagetitle: Survival models
#'     slide-number: true
#'     echo: true
#' lang: en
#' title: lecture-survival-script
#' ---
#' 

#| label: sim0
#| out-width: 60%
#| fig-align: center
#| fig-width: 9
set.seed(7829)
n <- 100                          ## Sample size
hazard0 <- 1/100
survival.time0 <- rexp(n, hazard0)
hist(survival.time0, xlab="Survival time (time to mortality)", main="")
abline(v=c(mean(survival.time0), 100), col=c("blue", "black"))
legend("topright", legend=c("Sample mean", "Population mean"),
       col=c("blue", "black"), lwd=2)

#| label: sim1
#| out-width: 70%
#| fig-align: center
set.seed(7829)
n <- 200                          ## Sample size
x <- rnorm(n)                     ## Covariate
beta0 <- -5; beta1 <- -0.5        ## Hazard coefficients
hazard <- exp(beta0 + beta1*x) 
survival.time <- rexp(n, hazard)
## summary(survival.time)
hist(survival.time, xlab="Survival time (time to mortality)", main="")

#| label: cen
ctime <- 500    ## Censoring time
censored <- ifelse(survival.time>ctime, 1, 0)
table(censored) ## Most individuals died before censoring time

#| label: cen1
survival.time.c <- ifelse(censored, ctime, survival.time) 
summary(survival.time.c)

#| label: life-lines
#| out-width: 90%
#| fig-align: center
#| fig-width: 8
#| code-fold: true
plot(0, type="n", xlim=c(0, ctime), ylim=c(1,20),
     xlab="Time", ylab="Individual")
abline(h=1:20, col=gray(0.8))
segments(rep(0, 20), 1:20, survival.time.c[1:20], 1:20)
points(rep(0,20), 1:20, pch=16, col="seagreen3", cex=1)
points(survival.time.c[1:20], 1:20,
       pch=ifelse(survival.time.c[1:20]==ctime, 18, 17),
       col=ifelse(survival.time.c[1:20]==ctime, "blue", "red"))
legend(170, 24, c("mortality", "censoring"), pch=c(17,18), col=c("red", "blue"),
       xpd=TRUE, horiz=TRUE)

#| label: km
#| out-width: 60%
#| fig-align: center
library(survival)
y <- Surv(survival.time.c, 1-censored)
plot(y, xlab="Time", ylab="Survivorship", main="Kaplan-Meier")

#| label: survreg
summary(pm1 <- survreg(y ~ x, dist="exponential"))

#| label: plothaz
#| out-width: 50%
#| fig-align: center
#| code-fold: true
plot(function(x) exp(beta0 + beta1*x), from=-3, to=3,
     xlab="Covariate", ylab="Hazard", col="black", ylim=c(0, 0.035))
alpha.hat <- coef(pm1)   ## Extract the estimates
plot(function(x) 1/exp(alpha.hat[1] + alpha.hat[2]*x), ## survreg formulation
     from=-3, to=3, col="blue", add=TRUE, lty=2)
legend(0.5, 0.035, c("Estimated hazard", "Hazard"), col=c(4,1), lty=2:1)

#| label: coxreg
(cox1 <- coxph(y ~ x))

#| label: cox-S
#| out-width: 65%
#| fig-align: center
#| fig-width: 8
plot(survfit(cox1, newdata=data.frame(x=c(-2,0,2))),
     xlab="Time", ylab="Survivorship", col=1:3,
     main="Cox proportional hazards fit")
legend(350, 1, paste("x =", c(-2, 0, 2)), lty=1, col=1:3)

#| label: exp-jag2
#| size: footnotesize
#| echo: false
writeLines(readLines("surv-exp.jag"))

#| label: load
#| include: false
library(jagsUI)
library(coda)

#| label: jd-surv-exp
survival.time.jags <- survival.time.c
survival.time.jags[censored==1] <- NA
jd.exp <- list(survivalTime=survival.time.jags, censored=censored,
               censorTime=rep(ctime, length(censored)),
               x=x, n=length(censored))

#| label: ji-jp-exp
ji.exp <- function() list(beta0=rnorm(1), beta1=rnorm(1)) 
jp.exp <- c("beta0", "beta1")

#| label: js
#| cache: true
#| results: hide
#| eval: true
jags.post.samples.exp <- jags.basic(data=jd.exp, inits=ji.exp,
                                    parameters.to.save=jp.exp,
                                    model.file="surv-exp.jag",
                                    n.chains=3, n.adapt=100, n.burnin=0,
                                    n.iter=2000, parallel=TRUE)

#| label: post-samps-exp
#| out-width: 65%
#| fig-align: center
plot(jags.post.samples.exp)

#| label: haz-post
#| out-width: 90%
#| fig-align: center
#| fig-width: 8
#| echo: false
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

#| label: sim-norobust
set.seed(34918)  
maxTime <- 10           ## Time period
n <- 25                 ## nIndividuals
phi <- 0.7              ## Survival probability over 1 time step

#| label: sim-norbust2
z <- matrix(NA, n, maxTime)
first <- rpois(n, 1)+1  ## random release dates
for(i in 1:n) {
    z[i,first[i]] <- 1  ## Known alive at release
    for(t in (first[i]+1):maxTime) {
        z[i,t] <- rbinom(1, 1, z[i,t-1]*phi) ## Alive/dead state
    }
}

#| label: z
z[1:15,]

#| label: jags-surv-dtime2
#| echo: false
writeLines(readLines("surv-dtime.jag"))

#| label: jd-dtime0
jd.dtime0 <- list(z=z, n=n, first=first, maxTime=maxTime)

#| label: ji-dtime0
ji.dtime0 <- function() list(phi=runif(1))
jp.dtime0 <- c("phi")

#| label: js-dtime0
#| cache: true
#| results: hide
#| eval: true
jags.post.samples.dtime0 <- jags.basic(data=jd.dtime0, inits=ji.dtime0,
                                       parameters.to.save=jp.dtime0,
                                       model.file="surv-dtime.jag",
                                       n.chains=3, n.adapt=100, n.burnin=0,
                                       n.iter=2000, parallel=TRUE)

#| label: post-samps-dtime0
#| out-width: 65%
#| fig-align: center
plot(jags.post.samples.dtime0)

#| label: phi-post
#| out-width: 70%
#| fig-align: center
#| fig-width: 8
#| code-fold: true
phi.post <- as.matrix(jags.post.samples.dtime0)[,"phi"]
S.post <- sapply(phi.post, function(x) x^(0:10))
par(mai=c(0.9,0.9,0.1,0.1))
matplot(0:10, S.post[,1:1000], type="l", xlab="Time", ylab="Survivorship",
        col=rgb(0,0.5,1,0.02))
lines(0:10, rowMeans(S.post), col="royalblue", lwd=3)
lines(0:10, phi^(0:10), col="cadetblue", lwd=3)
legend(6, 1, c("Posterior sample", "Posterior mean", "Actual"),
       col=c(rgb(0,0.5,1,0.2), "royalblue", "cadetblue"), lty=1, lwd=c(1,3,3))

#| label: jags-surv-dtime2-tcovs
#| echo: false
writeLines(readLines("surv-dtime-tcovs.jag"))

#| label: comp-risk-sim
beta0 <- c(-6, -4, -3) ## log-hazard. No covariates
lambda <- exp(beta0)   ## Cause-specific hazard
pi <- lambda / (1+sum(lambda))
pi[4] <- 1-sum(pi)     ## Probability of surviving t to t+1

#| label: comp-risk-sim2
nDeer <- 100
nDays <- 100
z <- matrix(NA, nDeer, nDays)
z[,1] <- 4  ## Everyone starts in state 4 (alive)

#| label: comp-risk-Phi
Phi <- diag(4)
rownames(Phi) <- colnames(Phi) <- c("cougar", "bear", "wolf", "alive")
Phi[4,] <- pi
# kable(Phi, format="latex", digits=3, booktabs=TRUE)

#| label: comp-risk-z
for(i in 1:nDeer) {
    for(t in 2:nDays) {
        z[i,t] <- which(rmultinom(n=1, size=1, prob=Phi[z[i,t-1],])==1)
    }
}

#| label: jags-surv-dtime2-comp-risk
#| echo: false
writeLines(readLines("surv-dtime-comp-risks.jag"))

#| label: jd-comp-risk
Phi.data <- diag(4)
Phi.data[4,] <- NA
jd.comp.risk <- list(z=z, n=nDeer, first=rep(1,nrow(z)), maxTime=ncol(z),
                     nRisks=3, Phi=Phi.data)

#| label: ji-comp-risk
ji.comp.risk <- function() list(beta0=log(runif(3)))
jp.comp.risk <- c("beta0", "pi")

#| label: js-comp-risk
#| results: hide
#| eval: true
#| message: false
jags.ps.comp.risk <- jags.basic(data=jd.comp.risk, inits=ji.comp.risk,
                                parameters.to.save=jp.comp.risk,
                                model.file="surv-dtime-comp-risks.jag",
                                n.chains=3, n.adapt=100, n.burnin=0,
                                n.iter=2000, parallel=TRUE)

#| label: post-samps-comp-risk
#| out-width: 75%
#| fig-align: center
#| echo: false
plot(jags.ps.comp.risk[,c("beta0[1]", "beta0[2]", "beta0[3]")])

#| label: post-samps-comp-risk-pi
#| out-width: 75%
#| fig-align: center
#| fig-width: 7
#| echo: false
#| fig-height: 9
plot(jags.ps.comp.risk[,c("pi[1]", "pi[2]", "pi[3]", "pi[4]")])

