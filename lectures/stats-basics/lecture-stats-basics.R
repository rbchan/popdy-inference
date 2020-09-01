## ----buildit,include=FALSE,eval=FALSE-----------------------------------------
## rnw2pdf("lecture-stats-basics")
## rnw2pdf("lecture-stats-basics", tangle=TRUE)




## ----linmod,include=FALSE-----------------------------------------------------
set.seed(3400)  
x1 <- runif(100, 0, 50)  
y <- rnorm(100, 10 + 1*x1, 5)
plot(x1, y)
abline(lm(y~x1))


## ----linmod-out,size='tiny'---------------------------------------------------
lm(y~x1)


## ----linmod-xc,include=FALSE--------------------------------------------------
set.seed(3400)  
xc <- gl(4, 25) 
y <- rnorm(100, model.matrix(~xc)%*%c(10,1,-1,2), 5)
ym <- tapply(y, xc, mean)
yse <- sqrt(sum(resid(lm(y~xc))^2)/96)/sqrt(25)
bpx <- barplot(ym, ylim=c(0, 15), xlab="Treatment group",
               ylab="Group mean", cex.lab=1.3)
arrows(bpx, ym, bpx, ym+yse, angle=90, code=3, length=0.05)


## ----linmod-xc-out,size='tiny'------------------------------------------------
lm(y~xc)


## ----read-grouse,size='tiny'--------------------------------------------------
grouse.data <- read.csv("grouse_data_glm.csv", row.names=1)
grouse.data$route <- factor(grouse.data$route)
grouse.data$utmZone <- factor(grouse.data$utmZone)
str(grouse.data)


## ----grouse-fm1,size='scriptsize'---------------------------------------------
fm1 <- lm(abundance ~ elevation + utmZone, grouse.data)
summary(fm1)


## ----grouse-pred-dat,size='scriptsize'----------------------------------------
elev.min <- min(grouse.data$elevation)
elev.max <- max(grouse.data$elevation)
seq.length <- 20 ## Determines how smooth the function looks in GLMs 
elev.seq <- seq(from=elev.min, to=elev.max, length.out=seq.length)
pred.data.west <- data.frame(elevation=elev.seq, utmZone="16S")
pred.data.east <- data.frame(elevation=elev.seq, utmZone="17S")


## ----grouse-pred,size='scriptsize'--------------------------------------------
pred.west <- predict(fm1, newdata=pred.data.west, se=TRUE)
pred.east <- predict(fm1, newdata=pred.data.east, se=TRUE)


## ----grouse-pred-plot,fig.width=7,fig.height=5,out.width="0.85\\textwidth",fig.align='center',size='scriptsize'----
plot(abundance ~ elevation, data=grouse.data, ylim=c(0,2))
lines(elev.seq, pred.west$fit, col="blue", lwd=2)
lines(elev.seq, pred.east$fit, col="grey", lwd=2)
legend(900, 2, c("West", "East"), lty=1, col=c("blue","grey"), lwd=2)


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

