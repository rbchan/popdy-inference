## ----sim-pars,size='scriptsize'-----------------------------------------------
T <- 10      # years/primary periods
K <- 3       # 3 secondary sampling occasion
lambda <- 25 # Expected value of abundance in year 1
M <- 500     # Easiest way to simulate data is using data augmentation
phi <- 0.7   # Apparent survival
gamma <- 0.3 # Per-capital recruitment rate
p0 <- 0.4    # Baseline capture prob
sigma <- 0.1 # Scale parameter of encounter rate function


## ----sim-p,size='scriptsize'--------------------------------------------------
set.seed(340)
co <- seq(0.25, 0.75, length=5)
x <- cbind(rep(co, each=5), rep(co, times=5))
J <- nrow(x)  ## nTraps
xlim <- ylim <- c(0,1)
## Activity centers, no dispersal
s <- cbind(runif(M, xlim[1], xlim[2]), runif(M, ylim[1], ylim[2]))
d <- p <- matrix(NA, M, J)
for(i in 1:M) {
    d[i,] <- sqrt((s[i,1]-x[,1])^2 + (s[i,2]-x[,2])^2)
    p[i,] <- p0*exp(-d[i,]^2/(2*sigma^2)) } # capture prob 


## ----sim-SR,size='scriptsize'-------------------------------------------------
set.seed(034)
z <- recruitable <- died <- recruited <- matrix(0, M, T)
z[,1] <- rbinom(M, 1, lambda/M) # alive at t=1
recruitable[,1] <- 1-z[,1]
N <- integer(T)
N[1] <- sum(z[,1])
for(t in 2:T) {
    ER <- N[t-1]*gamma # expected number of recruits
    prevA <- sum(recruitable[,t-1]) # Number available to be recruited
    gammaPrime <- ER/prevA
    if(gammaPrime > 1) stop("M isn't big enough")
    z[,t] <- rbinom(M, 1, (1-recruitable[,t-1])*z[,t-1]*phi +
                          recruitable[,t-1]*gammaPrime)
    recruitable[,t] <- recruitable[,t-1]*(1-z[,t])
    N[t] <- sum(z[,t])  }


## ----sim-N,size='scriptsize'--------------------------------------------------
died <- (z[,1:(T-1)]==1) & (z[,2:T]==0)
recruited <- (z[,1:(T-1)]==0) & (z[,2:T]==1)
Deaths <- colSums(died)
Recruits <- colSums(recruited)
everAlive <- sum(rowSums(z)>0)


## ----dynamics1,size='tiny',fig.height=5,fig.align='center',out.width='90%',echo=-1----
par(mai=c(0.9,0.9,0.1,0.1))
plot(1:T, N, type="b", xlab="Year", ylab="", ylim=c(0, 50), pch=16)
lines(2:T, Deaths, type="b", col="red2", pch=16); lines(2:T, Recruits, type="b", col="seagreen2", pch=16)
legend(1, 50, c("Abundance","Mortalities","Recruits"), lty=1, pch=16, col=c("black","red2","seagreen2"))


## ----sim-yall,size='footnotesize'---------------------------------------------
yall.bern <- array(NA, c(M, J, K, T))   ## For Bernoulli
yall <- array(NA, c(M, J, T))           ## For Binomial
for(i in 1:M) {
    for(t in 1:T) {
        for(j in 1:J) {
            yall.bern[i,j,1:K,t] <- rbinom(K, 1, z[i,t]*p[i,j])
            yall[i,j,t] <- rbinom(1, K, z[i,t]*p[i,j])
        }
    }
}


## ----sim-y,size='footnotesize'------------------------------------------------
y.bern <- yall.bern[rowSums(yall.bern)>0,,,]
y <- yall[rowSums(yall)>0,,]
str(y)


## ----mask,size='scriptsize',out.width="60%",fig.align="center",results='hide'----
library(openpopscr)
trap.df <- data.frame(x*1000); colnames(trap.df) <- c("x","y")
traps <- read.traps(data=trap.df, detector="proximity")
mask <- make.mask(traps=traps, buffer=250)
plot(mask); points(traps, pch=3, col="blue", lwd=2)


## ----JS-model-new,size='scriptsize'-------------------------------------------
y.secr <- y.bern
year <- rep(slice.index(y.bern, 4), y.secr)  ## Primary period
day <- rep(slice.index(y.bern, 3), y.secr)   ## Secondary period
caps <- data.frame(session=1,
                   animal=rep(slice.index(y.bern, 1), y.secr),
                   occasion=(year-1)*K+day,
                   trap=rep(slice.index(y.bern, 2), y.secr))
capthist <- make.capthist(captures=caps, traps=traps, noccasions=T*K)


## ----format-openpop,size='scriptsize'-----------------------------------------
js.data <- ScrData$new(capthist, mask, primary=rep(1:T, each=K))


## ----js-mod-fit,size='scriptsize',results='hide',cache=TRUE,eval=FALSE--------
## start <- get_start_values(js.data, model = "JsModel")
## mod <- JsModel$new(list(lambda0~1, sigma~1, D~1, phi~1, beta~1), js.data,
##                    start=start)
##                    start=list(lambda0=0.1, sigma=100, D=100,
##                               phi=0.8, beta=0.2))
## mod$fit()
## ## mod


## ----cjs-mod-est,size='scriptsize',eval=FALSE---------------------------------
## mod$get_par("lambda0", k = 1, j = 1)
## mod$get_par("sigma", k = 1, j = 1)
## mod$get_par("D")                 ## Superpopulation density
## mod$get_par("beta", k = 1, m=1)  ## Proportion of superpop alive in year 1
## mod$get_par("phi", k = 1, m=1)   ## Survival


## ----jagsmod1,size='tiny',comment='',echo=FALSE,background='lightblue'--------
writeLines(readLines("JS-spatial.jag"))


## ----load-jagsUI-coda,include=FALSE-------------------------------------------
library(jagsUI)
library(coda)


## ----sim-aug,size='scriptsize'------------------------------------------------
M <- nrow(y) + 50   ## Trial and error
yz <- array(0, c(M, J, T))
yz[1:nrow(y),,] <- y
jags.data1 <- list(y=yz, M=M, x=x, J=J, K=K, T=T, xlim=xlim, ylim=ylim)


## ----zi,size='scriptsize'-----------------------------------------------------
zi <- matrix(0, M, T)
zi[1:nrow(y),] <- 1
ji1 <- function() list(phi=0.01, gamma=0.001, z=zi)
jp1 <- c("phi", "gamma", "p0", "sigma", "N", "Deaths", "Recruits", "Ntot")


## ----jags-run,size='scriptsize',results='hide',cache=TRUE---------------------
library(jagsUI)
jags.post.samples1 <- jags.basic(data=jags.data1, inits=ji1,
                                 parameters.to.save=jp1,
                                 model.file="JS-spatial.jag",
                                 n.chains=3, n.adapt=100, n.burnin=0,
                                 n.iter=2000, parallel=TRUE)


## ----Ntot,size='tiny',out.width='70%',fig.align='center'----------------------
hist(as.matrix(jags.post.samples1[,"Ntot"]), xlab="Total population size",
     ylab="", main="", freq=FALSE, xlim=c(nrow(y), M))
abline(v=M, lwd=3, col="blue")


## ----jc1,out.width='60%',fig.align='center',size='scriptsize'-----------------
plot(jags.post.samples1[,c("phi", "gamma", "p0", "sigma")])


## ----jcN1-4,include=FALSE,echo=FALSE------------------------------------------
plot(jags.post.samples1[,c("N[1]", "N[2]", "N[3]", "N[4]")])

## ----jcN5-8,include=FALSE,echo=FALSE------------------------------------------
plot(jags.post.samples1[,c("N[5]", "N[6]", "N[7]", "N[8]")])


## ----N-post-samples,size='scriptsize'-----------------------------------------
Npost <- as.matrix(jags.post.samples1[,paste("N[", 1:10, "]", sep="")])
Nmed <- apply(Npost, 2, median)
Nupper <- apply(Npost, 2, quantile, prob=0.975)
Nlower <- apply(Npost, 2, quantile, prob=0.025)


## ----Npost,size='scriptsize',fig.height=5,fig.show='hide',echo=-1-------------
par(mai=c(0.9,0.9,0.1,0.1))  
plot(1:T, N, type="b", ylim=c(0, 60), xlab="Time",
     ylab="Abundance", pch=16)
arrows(1:T, Nlower, 1:T, Nupper, angle=90, code=3,
       length=0.05, col=gray(0.7))
points(1:T, Nmed, pch=16, col=gray(0.7))
legend(1, 60, c("Actual", "Estimated"),
       col=c("black", gray(0.7)), lty=c(NA,1), pch=c(16,16))


## ----R-post-samples,size='scriptsize'-----------------------------------------
Rpost <- as.matrix(jags.post.samples1[,paste("Recruits[",1:9,"]",sep="")])
Rmed <- apply(Rpost, 2, median)
Rupper <- apply(Rpost, 2, quantile, prob=0.975)
Rlower <- apply(Rpost, 2, quantile, prob=0.025)


## ----Rpost,size='scriptsize',fig.height=5,fig.show='hide',echo=-1-------------
par(mai=c(0.9,0.9,0.1,0.1))  
plot(1:(T-1), Recruits, type="b", ylim=c(0, 30), xlab="Time",
     ylab="Recruits", pch=16, col="seagreen2")
arrows(1:(T-1), Rlower, 1:(T-1), Rupper, angle=90, code=3,
       length=0.05, col=gray(0.7))
points(1:(T-1), Rmed, pch=16, col=gray(0.7))
legend(1, 30, c("Actual", "Estimated"),
       col=c("seagreen2", gray(0.7)), lty=c(NA,1), pch=c(16,16))


## ----D-post-samples,size='scriptsize'-----------------------------------------
Dpost <- as.matrix(jags.post.samples1[,paste("Deaths[", 1:9, "]", sep="")])
Dmed <- apply(Dpost, 2, median)
Dupper <- apply(Dpost, 2, quantile, prob=0.975)
Dlower <- apply(Dpost, 2, quantile, prob=0.025)


## ----Dpost,size='scriptsize',fig.height=5,fig.show='hide',echo=-1-------------
par(mai=c(0.9,0.9,0.1,0.1))  
plot(1:(T-1), Deaths, type="b", ylim=c(0, 30), xlab="Time",
     ylab="Mortalities", pch=16, col="red2")
arrows(1:(T-1), Dlower, 1:(T-1), Dupper, angle=90, code=3,
       length=0.05, col=gray(0.7))
points(1:(T-1), Dmed, pch=16, col=gray(0.7))
legend(1, 30, c("Actual", "Estimated"),
       col=c("red2", gray(0.7)), lty=c(NA,1), pch=c(16,16))

