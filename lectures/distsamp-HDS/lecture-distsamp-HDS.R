## ----buildit,include=FALSE,eval=FALSE-----------------------------------------
## rnw2pdf("lecture-distsamp-HDS")
## rnw2pdf("lecture-distsamp-HDS", tangle=TRUE)




## ----hn,size='scriptsize',fig.width=7,fig.height=5,out.width="70%",echo=-1----
par(mai=c(0.9,0.9,0.1,0.1))  
sigma1 <- 25; sigma2 <- 50
plot(function(x) exp(-x^2/(2*sigma1^2)), from=0, to=100, 
     xlab="Distance (x)", ylab="Detection probability (p)")
plot(function(x) exp(-x^2/(2*sigma2^2)), from=0, to=100, add=TRUE, col=4)
legend(70, 1, c("sigma=25", "sigma=50"), lty=1, col=c("black","blue"))


## ----nexp,size='scriptsize',fig.width=7,fig.height=5,out.width="70%",echo=-1----
par(mai=c(0.9,0.9,0.1,0.1))  
sigma1 <- 25; sigma2 <- 50
plot(function(x) exp(-x/sigma1), from=0, to=100, 
     xlab="Distance (x)", ylab="Detection probability (p)")
plot(function(x) exp(-x/sigma2), from=0, to=100, add=TRUE, col=4)
legend(70, 1, c("sigma=25", "sigma=50"), lty=1, col=c("black","blue"))


## ----haz,size='scriptsize',fig.width=7,fig.height=5,out.width="70%",echo=-1----
par(mai=c(0.9,0.9,0.1,0.1))  
a1 <- 25; a2 <- 50; b1 <- 2; b2 <- 10
plot(function(x) 1-exp(-(x/a1)^(-b1)), from=0, to=100, 
     xlab="Distance (x)", ylab="Detection probability (p)", ylim=0:1)
plot(function(x) 1-exp(-(x/a2)^(-b2)), from=0, to=100, add=TRUE, col=4)
legend(70, 1, c("a=25, b=2", "a=50, b=10"), lty=1, col=c("black","blue"))


## ----pbar-hn,size='footnotesize'----------------------------------------------
B <- 100                                         # Transect width
g <- function(x, sigma=25) exp(-x^2/(2*sigma^2)) # g(x)
pdf <- function(x) 1/B                           # p(x), constant


## ----pbar-hn-int,size='footnotesize'------------------------------------------
gp <- function(x) g(x)*pdf(x)
(pbar <- integrate(gp, lower=0, upper=B)$value)


## ----pbar-hn-int2,size='footnotesize'-----------------------------------------
(pbar <- integrate(g, lower=0, upper=B)$value / B)
(pbar <- (pnorm(B,0,25) - pnorm(0,0,25)) / dnorm(0,0,25) / B)


## ----pbar-hn-pt,size='footnotesize'-------------------------------------------
B <- 100                                         # Transect width
g <- function(x, sigma=25) exp(-x^2/(2*sigma^2)) # g(x)
pdf <- function(x) 2*x/B^2                       # p(x)


## ----pbar-hn-int-pt,size='footnotesize'---------------------------------------
gp <- function(x) g(x)*pdf(x)
(pbar <- integrate(gp, lower=0, upper=B)$value)


## ----pbar-hn-int2-pt,size='footnotesize'--------------------------------------
sigma <- 25
(pbar <- (sigma^2*(1-exp(-B^2/(2*sigma^2))) -
          sigma^2*(1-exp(-0^2/(2*sigma^2)))) * 2*pi/(pi*B^2))


## ----include=FALSE,echo=FALSE-------------------------------------------------
set.seed(34889243)


## ----sim-hds-nocov1,size='scriptsize'-----------------------------------------
nSites <- 100; lambda1 <- 2.6  ## Expected value of N
N1 <- rpois(n=nSites, lambda=lambda1)


## ----sim-hds-nocov2,size='scriptsize'-----------------------------------------
J <- 5                     # distance bins
sigma <- 50                # scale parameter
B <- 100; L <- 100         # transect widths (B) and lengths (L)
b <- seq(0, B, length=J+1) # distance break points
psi <- diff(b)/B           # Pr(x is in bin j)
pbar1 <- numeric(J)        # average detection probability
pi1 <- numeric(J+1)        # multinomial cell probs
for(j in 1:J) {
    pbar1[j] <- integrate(function(x) exp(-x^2/(2*sigma^2)),
                          lower=b[j], upper=b[j+1])$value / diff(b)[j]
    pi1[j] <- pbar1[j]*psi[j]
}
pi1[J+1] <- 1-sum(pi1[1:J])


## ----sim-hds-y1,size='scriptsize'---------------------------------------------
y1.all <- matrix(NA, nrow=nSites, ncol=J+1)
for(i in 1:nSites) {
    y1.all[i,] <- rmultinom(n=1, size=N1[i], prob=pi1)    }
y1 <- y1.all[,1:J]  ## Drop final cell


## ----dist-hist2,fig.height=5,out.width="90%",size='scriptsize',echo=-1--------
par(mai=c(0.9,0.9,0.1,0.1))  
plot(b[-(J+1)]+10, colSums(y1), type="h", lwd=80, lend=2, col="skyblue4",
     xlim=c(0,100), ylim=c(0, 70), xlab="Distance", ylab="Detections")


## ----sim-hds-cov1,size='scriptsize'-------------------------------------------
elevation <- rnorm(nSites)  
beta0 <- 2; beta1 <- 1
lambda2 <- exp(beta0 + beta1*elevation)
N2 <- rpois(n=nSites, lambda=lambda2)


## ----sim-hds-cov2,size='scriptsize'-------------------------------------------
noise <- rnorm(nSites)
alpha0 <- 3; alpha1 <- -0.5
sigma2 <- exp(alpha0 + alpha1*noise)
pi2 <- matrix(NA, nSites, J+1) # multinomial cell probs
for(i in 1:nSites) {
  for(j in 1:J) {
      pi2[i,j] <- integrate(function(x) exp(-x^2/(2*sigma2[i]^2)),
          lower=b[j], upper=b[j+1])$value / (b[j+1]-b[j]) * psi[j] }
  pi2[i,J+1] <- 1-sum(pi2[i,1:J]) } 


## ----sim-hds-y2,size='scriptsize'---------------------------------------------
y2.all <- matrix(NA, nrow=nSites, ncol=J+1)
for(i in 1:nSites) {
    y2.all[i,] <- rmultinom(n=1, size=N2[i], prob=pi2[i,])    }
y2 <- y2.all[,1:J]


## ----sim-nocov-dat,size='tiny'------------------------------------------------
y2[1:25,]


## ----sim-nocov-ss1,size='scriptsize'------------------------------------------
# Max count at each site
maxCounts <- apply(y2, 1, max) 
naiveOccupancy <- sum(maxCounts>0)/nSites
naiveOccupancy 


## ----sim-nocov-ss2,size='scriptsize'------------------------------------------
colSums(y2)


## ----sim-nocov-ss3,size='scriptsize'------------------------------------------
sum(y2)

## ----un,include=FALSE---------------------------------------------------------
library(unmarked)


## ----un-umf,size='tiny'-------------------------------------------------------
umf <- unmarkedFrameDS(y=y2, siteCovs=data.frame(elevation,noise), dist.breaks=b,
                       tlength=rep(L, nSites), survey="line", unitsIn="m")


## ----wfac,size='tiny'---------------------------------------------------------
summary(umf)


## ----un-fit,size='tiny'-------------------------------------------------------
## fm <- distsamp(~noise ~elevation, umf, keyfun="exp")     # negative exp
## fm <- distsamp(~noise ~elevation, umf, keyfun="hazard")  # hazard rate
fm <- distsamp(~noise ~elevation, umf, keyfun="halfnorm")   # half-normal
fm


## ----un-compare,size='tiny'---------------------------------------------------
c(beta0=beta0, beta1=beta1); c(alpha0=alpha0, alpha1=alpha1)


## ----preddat,size='footnotesize'----------------------------------------------
pred.data <- data.frame(noise=seq(-3, 3, length=20))


## ----predpsi,size='footnotesize'----------------------------------------------
sigma.pred <- predict(fm, newdata=pred.data,
                      type='det', append=TRUE)


## ----psi-head,size='footnotesize'---------------------------------------------
print(head(sigma.pred), digits=2)


## ----pred-sigma,fig.width=7,fig.height=5,size='tiny',out.width='80%',fig.align='center',echo=-c(1)----
par(mai=c(0.9,0.9,0.1,0.1))
plot(Predicted ~ noise, sigma.pred, ylab="Scale parameter (sigma)",
     ylim=c(0,100), xlab="Noise level", type="l")
lines(lower ~ noise, sigma.pred, col="grey")
lines(upper ~ noise, sigma.pred, col="grey")


## ----pred-sigma2,fig.width=7,fig.height=5,size='tiny',out.width='80%',fig.align='center',echo=-c(1,5:10)----
par(mai=c(0.9,0.9,0.1,0.1))
plot(Predicted ~ noise, sigma.pred, ylab="Scale parameter (sigma)",
     ylim=c(0,100), xlab="Noise level", type="l")
lines(lower ~ noise, sigma.pred, col="grey")
lines(upper ~ noise, sigma.pred, col="grey")
arrows((log(40)-3.137)/-0.369, 40, 0.1, 63, len=0.1, col="blue")
points((log(40)-3.137)/-0.369, 40, col="orange", pch=16, cex=2)
par( fig=c(0.5,0.95,0.55,0.95), new=TRUE, las=1)#, mar=c(0,0,0,0) )
plot(function(x) exp(-x^2/(2*40^2)), 0, 100, xlab="Dist", ylab="p",
     col="orange", lwd=2, ylim=c(0,1), las=1)


## ----pred-sigma2-2,fig.width=7,fig.height=5,size='tiny',out.width='80%',fig.align='center',echo=-c(1,5:10)----
par(mai=c(0.9,0.9,0.1,0.1))
plot(Predicted ~ noise, sigma.pred, ylab="Scale parameter (sigma)",
     ylim=c(0,100), xlab="Noise level", type="l")
lines(lower ~ noise, sigma.pred, col="grey")
lines(upper ~ noise, sigma.pred, col="grey")
arrows((log(15)-3.137)/-0.369, 15, 1.1, 55, len=0.1, col="blue")
points((log(15)-3.137)/-0.369, 15, col="red", pch=16, cex=2)
par( fig=c(0.5,0.95,0.55,0.95), new=TRUE, las=1)#, mar=c(0,0,0,0) )
plot(function(x) exp(-x^2/(2*15^2)), 0, 100, xlab="Dist", ylab="p",
     col="red", lwd=2, ylim=c(0,1), las=1)


## ----bugs-line,size='tiny',echo=FALSE,comment='',background='lightblue'-------
writeLines(readLines("distsamp-line-mod.jag"))


## ----bugs-data2,size='footnotesize'-------------------------------------------
jags.data.line <- list(y=y2, n=rowSums(y2),
                       b=b,           # Distance break points
                       psi=diff(b)/B, # Pr(occuring in bin j)
                       elevation=elevation, noise=noise,
                       nSites=nSites, nBins=J)


## ----bugs-inits,size='footnotesize'-------------------------------------------
jags.inits.line <- function() {
    list(lambda.intercept=runif(1), alpha0=rnorm(1, 5),
         N=rowSums(y2)+rpois(nrow(y2), 2))
}


## ----bugs-pars,size='small'---------------------------------------------------
jags.pars.line <- c("beta0", "beta1",
                    "alpha0", "alpha1", "totalAbundance")


## ----bugs-mcmc-line,size='scriptsize',message=FALSE,cache=TRUE,results='hide'----
library(jagsUI)
jags.post.line <- jags.basic(data=jags.data.line, inits=jags.inits.line,
                             parameters.to.save=jags.pars.line,
                             model.file="distsamp-line-mod.jag",
                             n.chains=3, n.adapt=100, n.burnin=0,
                             n.iter=2000, parallel=TRUE)


## ----jags-sum-line,size='scriptsize',cache=TRUE-------------------------------
round(summary(jags.post.line)$quantile, digits=3)


## ----bugs-plot1-rem2,size='footnotesize',out.width="0.7\\textwidth",fig.align='center',cache=TRUE----
plot(jags.post.line[,jags.pars.line[1:3]])


## ----bugs-plot2-rem2,size='footnotesize',out.width="0.7\\textwidth",fig.align='center',cache=TRUE----
plot(jags.post.line[,jags.pars.line[4:5]])


## ----include=FALSE,echo=FALSE-------------------------------------------------
set.seed(34889243)


## ----sim-hds-nocovN3-pt,size='scriptsize'-------------------------------------
nSites <- 100; lambda1 <- 2.6  ## Expected value of N
N3 <- rpois(n=nSites, lambda=lambda1)


## ----sim-hds-nocov3-pt,size='scriptsize'--------------------------------------
J <- 5                     # distance bins
sigma <- 50                # scale parameter
B <- 100                   # plot radius (B)
b <- seq(0, B, length=J+1) # distance break points
area <- pi*b^2             # area of each circle
psi <- (area[-1]-area[-(J+1)]) / area[J+1]
pbar3 <- numeric(J)        # average detection probability
pi3 <- numeric(J+1)        # multinomial cell probs
for(j in 1:J) {
    pbar3[j] <- integrate(function(x) exp(-x^2/(2*sigma^2))*x,
                          lower=b[j], upper=b[j+1])$value *
                          (2*pi/diff(area)[j])
    pi3[j] <- pbar3[j]*psi[j] }; pi3[J+1] <- 1-sum(pi3[1:J])


## ----sim-hds-y3,size='scriptsize'---------------------------------------------
y3.all <- matrix(NA, nrow=nSites, ncol=J+1)
for(i in 1:nSites) {
    y3.all[i,] <- rmultinom(n=1, size=N3[i], prob=pi3)    }
y3 <- y3.all[,1:J]  ## Drop final cell


## ----dist-hist3,fig.height=5,out.width="90%",size='scriptsize',echo=-1--------
par(mai=c(0.9,0.9,0.1,0.1))  
plot(b[-(J+1)]+10, colSums(y3), type="h", lwd=80, lend=2, col="skyblue4",
     xlim=c(0,100), ylim=c(0, 70), xlab="Distance", ylab="Detections")


## ----sim-hds-covN4-pt,size='scriptsize'---------------------------------------
elevation <- rnorm(nSites)  
beta0 <- 2; beta1 <- 1
lambda4 <- exp(beta0 + beta1*elevation)
N4 <- rpois(n=nSites, lambda=lambda4)


## ----sim-hds-cov4-pt,size='scriptsize'----------------------------------------
noise <- rnorm(nSites)
alpha0 <- 3; alpha1 <- -0.5
sigma4 <- exp(alpha0 + alpha1*noise)
pi4 <- matrix(NA, nSites, J+1) # multinomial cell probs
for(i in 1:nSites) {
  for(j in 1:J) {
      pi4[i,j] <- integrate(function(x) exp(-x^2/(2*sigma4[i]^2))*x,
          lower=b[j], upper=b[j+1])$value*(2*pi/diff(area)[j])*psi[j] }
  pi4[i,J+1] <- 1-sum(pi4[i,1:J]) } 


## ----sim-hds-y4-pt,size='scriptsize'------------------------------------------
y4.all <- matrix(NA, nrow=nSites, ncol=J+1)
for(i in 1:nSites) {
    y4.all[i,] <- rmultinom(n=1, size=N4[i], prob=pi4[i,])    }
y4 <- y4.all[,1:J]


## ----sim-nocov-dat4-pt,size='tiny'--------------------------------------------
y4[1:25,]


## ----sim-nocov-ssO-pt,size='scriptsize'---------------------------------------
# Max count at each site
maxCounts <- apply(y4, 1, max) 
naiveOccupancy <- sum(maxCounts>0)/nSites
naiveOccupancy 


## ----sim-nocov-ss4-pt,size='scriptsize'---------------------------------------
colSums(y4)


## ----sim-nocov-ss4sum-pt,size='scriptsize'------------------------------------
sum(y4)


## ----un-umf-pt,size='tiny'----------------------------------------------------
umf4 <- unmarkedFrameDS(y=y4, siteCovs=data.frame(elevation,noise), dist.breaks=b,
                       survey="point", unitsIn="m")


## ----wfac-pt,size='tiny'------------------------------------------------------
summary(umf4)


## ----un-fit-pt,size='tiny'----------------------------------------------------
## fm4 <- distsamp(~noise ~elevation, umf4, keyfun="exp")     # negative exp
## fm4 <- distsamp(~noise ~elevation, umf4, keyfun="hazard")  # hazard rate
fm4 <- distsamp(~noise ~elevation, umf4, keyfun="halfnorm")   # half-normal
fm4


## ----un-compare-pt,size='tiny'------------------------------------------------
c(beta0=beta0, beta1=beta1); c(alpha0=alpha0, alpha1=alpha1)


## ----bugs-pt,size='tiny',echo=FALSE,comment='',background='lightblue'---------
writeLines(readLines("distsamp-point-mod.jag"))


## ----bugs-data-pt,size='footnotesize'-----------------------------------------
jags.data.pt <- list(y=y4, n=rowSums(y4), area=diff(area),
                     b=b,           # Distance break points
                     psi=psi,       # Pr(occuring in bin j)
                     elevation=elevation, noise=noise,
                     nSites=nSites, nBins=J)


## ----bugs-inits-pt,size='footnotesize'----------------------------------------
jags.inits.pt <- function() {
    list(lambda.intercept=runif(1), alpha0=rnorm(1, 5),
         N=rowSums(y4)+rpois(nrow(y4), 2))
}


## ----bugs-pars-pt,size='small'------------------------------------------------
jags.pars.pt <- c("beta0", "beta1",
                  "alpha0", "alpha1", "totalAbundance")


## ----bugs-mcmc-pt,size='scriptsize',message=FALSE,cache=TRUE,results='hide'----
jags.post.pt <- jags.basic(data=jags.data.pt, inits=jags.inits.pt,
                           parameters.to.save=jags.pars.pt,
                           model.file="distsamp-point-mod.jag",
                           n.chains=3, n.adapt=100, n.burnin=0,
                           n.iter=2000, parallel=TRUE)


## ----jags-sum-pt,size='scriptsize',cache=TRUE---------------------------------
round(summary(jags.post.pt)$quantile, digits=3)

