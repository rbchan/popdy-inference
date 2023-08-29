## ----buildit,include=FALSE,eval=FALSE-----------------------------------------
## rnw2pdf("lecture-stats-basics")
## rnw2pdf("lecture-stats-basics", tangle=TRUE)




## ----slr-data,size='tiny'-----------------------------------------------------
width <- c(0.30, 0.91, 0.89, 0.24, 0.77, 0.56, 0.59, 0.92, 0.81, 0.59)  ## x
mass <- c(0.08, 0.59, 0.18, 0.17, 0.42, 0.71, 0.49, 0.75, 0.46, 0.04)   ## y


## ----slr-viz,echo=FALSE-------------------------------------------------------
plot(width, mass, xlim=c(0,1))
abline(lm(mass~width))


## ----slr-fit,size='scriptsize'------------------------------------------------
lm(mass~width)


## ----lm-jag,size="scriptsize",comment="",echo=FALSE,background='beige'--------
  writeLines(readLines("lm.jag"))


## ----lm-jd,size='footnotesize'------------------------------------------------
jd.lm <- list(x=width, y=mass, n=length(mass))
ji.lm <- function() c(beta0=rnorm(1), beta1=0, sigmaSq=runif(1))
jp.lm <- c("beta0", "beta1", "sigmaSq")


## ----lm-jags,size='scriptsize',results='hide',warning=FALSE,cache=FALSE-------
library(jagsUI)  
js.lm <- jags.basic(data=jd.lm, inits=ji.lm, parameters.to.save=jp.lm,
                    model.file="lm.jag", n.chains=1, n.iter=1000)


## ----lm-jags-sum,size='scriptsize'--------------------------------------------
round(summary(js.lm)$quant, 2)


## ----lm-jags-viz,echo=FALSE---------------------------------------------------
plot(js.lm)


## ----linmod-out,size='tiny'---------------------------------------------------
lm(mass~width)


## ----linmod-xc,include=FALSE--------------------------------------------------
set.seed(3400)  
species <- gl(4, 25) 
mass <- rnorm(100, model.matrix(~species)%*%c(10,1,-1,2), 5)
ym <- tapply(mass, species, mean)
yse <- sqrt(sum(resid(lm(mass~species))^2)/96)/sqrt(25)
bpx <- barplot(ym, ylim=c(0, 15), xlab="Species",
               ylab="Group mean", cex.lab=1.3)
arrows(bpx, ym, bpx, ym+yse, angle=90, code=3, length=0.05)


## ----linmod-xc-out,size='tiny'------------------------------------------------
lm(mass~species)


## ----read-grouse,size='tiny'--------------------------------------------------
grouse.data <- read.csv("grouse_data_glm.csv", row.names=1)
grouse.data$route <- factor(grouse.data$route)
grouse.data$utmZone <- factor(grouse.data$utmZone)
str(grouse.data)


## ----grouse-fm1,size='scriptsize'---------------------------------------------
fm1 <- lm(abundance ~ elevation, grouse.data)
summary(fm1)


## ----grouse-pred-dat,size='scriptsize'----------------------------------------
elev.min <- min(grouse.data$elevation)
elev.max <- max(grouse.data$elevation)
seq.length <- 20 ## Determines how smooth the function looks in GLMs 
elev.seq <- seq(from=elev.min, to=elev.max, length.out=seq.length)
pred.data <- data.frame(elevation=elev.seq)


## ----grouse-pred,size='scriptsize'--------------------------------------------
pred.elev <- predict(fm1, newdata=pred.data, se=TRUE)


## ----grouse-pred-plot,fig.width=7,fig.height=5,out.width="0.85\\textwidth",fig.align='center',size='scriptsize'----
plot(abundance ~ elevation, data=grouse.data, ylim=c(0,2))
lines(elev.seq, pred.elev$fit, col="blue", lwd=2)


## ----logit-p,size='tiny'------------------------------------------------------
beta0 <- 5
beta1 <- -0.08
elevation <- 100
(logit.p <- beta0 + beta1*elevation)


## ----inv-logit,size='tiny'----------------------------------------------------
p <- exp(logit.p)/(1+exp(logit.p))
p


## ----plogis,size='tiny'-------------------------------------------------------
plogis(logit.p)


## ----logit,size='tiny'--------------------------------------------------------
log(p/(1-p))
qlogis(p)


## ----nologit,fig.show='hide',fig.width=6,fig.height=4,size='scriptsize'-------
plot(function(x) 5 + -0.08*x, from=0, to=100,
     xlab="Elevation", ylab="logit(prob of occurrence)")


## ----logit2,fig.show='hide',fig.width=6,fig.height=4,size='scriptsize'--------
plot(function(x) plogis(5 + -0.08*x), from=0, to=100,
     xlab="Elevation", ylab="Probability of occurrence")


## ----simFrogs,echo=FALSE,results='hide'---------------------------------------
set.seed(43340)
n <- 30
elev <- round(runif(n, 0, 500))
habitat <- gl(3, 10, labels=c("Oak", "Maple", "Pine"))
beta0 <- -1
beta1 <- 0.01
mu <- plogis(beta0 + beta1*elev)
summary(mu)
y <- rbinom(n, 1, mu)
frogData <- data.frame(presence=y,
                       abundance=rpois(n, exp(beta0 + beta1*elev)),
                       elevation=elev, habitat)
glm1 <- glm(presence ~ elev+habitat, family=binomial(link="logit"), data=frogData)
summary(glm1)
anova(glm1)


## ----binom1,echo=FALSE,fig.width=7,fig.height=6,out.width="0.9\\textwidth"----
plot(0:5, dbinom(0:5, 5, 0.5), type="h",
     xlab="Number of 'successes'", ylab="Probability",
     lend="butt", lwd=5, col="blue", ylim=c(0,0.6),
     main="Binomial(N=5, p=0.5)")


## ----binom2,echo=FALSE,fig.width=7,fig.height=6,out.width="0.9\\textwidth"----
plot(0:5, dbinom(0:5, 5, 0.9), type="h",
     xlab="Number of 'successes'", ylab="Probability",
     lend="butt", lwd=5, col="blue", ylim=c(0,0.6),
     main="Binomial(N=5, p=0.9)")


## ----logitreg,size='scriptsize'-----------------------------------------------
logitreg1 <- glm(presence ~ elevation, data=grouse.data,
                 family=binomial(link="logit"))
logitreg1


## ----grouse-lrpred-dat,size='scriptsize'--------------------------------------
elev.min <- min(grouse.data$elevation)
elev.max <- max(grouse.data$elevation)
seq.length <- 20 ## Determines how smooth the function looks in GLMs 
elev.seq <- seq(from=elev.min, to=elev.max, length.out=seq.length)
pred.data.lr <- data.frame(elevation=elev.seq)


## ----grouse-lrpred,size='scriptsize'------------------------------------------
pred.elev <- predict(logitreg1, newdata=pred.data.lr, type="link",
                     se=TRUE)


## ----grouse-lrpred-plot,fig.width=7,fig.height=5,out.width="0.75\\textwidth",fig.align='center',size='tiny'----
plot(presence ~ elevation, data=grouse.data, ylim=c(0,1))
lines(elev.seq, plogis(pred.elev$fit), col="blue", lwd=2)
lines(elev.seq, plogis(pred.elev$fit+pred.elev$se.fit), col="blue", lwd=1, lty=2)
lines(elev.seq, plogis(pred.elev$fit-pred.elev$se.fit), col="blue", lwd=1, lty=2)


## ----pois1,fig.show='hide',echo=FALSE-----------------------------------------
x <- 0:25
plot(x, dpois(x, lambda=1), type="h", lwd=5, col="blue", lend="butt",
     xlab="Response variable", ylab="Probability",
     main=expression(paste("Poisson(", lambda, "= 1)", sep="")), cex.lab=1.5 )

## ----pois2,fig.show='hide',echo=FALSE-----------------------------------------
plot(x, dpois(x, lambda=5), type="h", lwd=5, col="blue", lend="butt",
     xlab="Response variable", ylab="Probability",
     main=expression(paste("Poisson(", lambda, "= 5)", sep="")), cex.lab=1.5 )

## ----pois3,fig.show='hide',echo=FALSE-----------------------------------------
plot(x, dpois(x, lambda=10), type="h", lwd=5, col="blue", lend="butt",
     xlab="Response variable", ylab="Probability",
     main=expression(paste("Poisson(", lambda, "= 10)", sep="")), cex.lab=1.5 )


## ----nolog,fig.show='hide',fig.width=7,fig.height=5,size='footnotesize'-------
plot(function(x) 5 + -0.08*x, from=0, to=100,
     xlab="Elevation", ylab="log(Expected abundance)")


## ----log,fig.show='hide',fig.width=7,fig.height=5,size='footnotesize'---------
plot(function(x) exp(5 + -0.08*x), from=0, to=100,
     xlab="Elevation", ylab="Expected abundance")


## ----pois-cov-----------------------------------------------------------------
n <- 100  
x <- rnorm(n)


## ----pois-lam-----------------------------------------------------------------
beta0 <- -1
beta1 <- 1
lam <- exp(beta0 + beta1*x)


## ----pois-y-------------------------------------------------------------------
y <- rpois(n=n, lambda=lam)


## ----pois-fit,size='small',eval=FALSE-----------------------------------------
## (poisreg1 <- glm(y ~ x, family=poisson(link="log")))


## ----glm-jag,size="scriptsize",comment="",echo=FALSE,background='beige'-------
  writeLines(readLines("glm.jag"))


## ----glm-jd-------------------------------------------------------------------
jd.glm <- list(x=x, y=y, n=length(y))
ji.glm <- function() c(beta0=rnorm(1), beta1=rnorm(1))
jp.glm <- c("beta0", "beta1")


## ----glm-jags,size='scriptsize',results='hide',warning=FALSE,cache=FALSE------
library(jagsUI)  
js.glm <- jags.basic(data=jd.glm, inits=ji.glm, parameters.to.save=jp.glm,
                     model.file="glm.jag", n.chains=2, n.iter=1000)


## ----glm-jags-sum,size='tiny'-------------------------------------------------
round(summary(js.glm)$quant, 2)


## ----glm-jags-viz,echo=FALSE--------------------------------------------------
plot(js.glm)

