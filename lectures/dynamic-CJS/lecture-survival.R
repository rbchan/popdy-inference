## ----sim-norobust,size='footnotesize',echo=-1---------------------------------
set.seed(34918)  
T <- 10               ## primary periods (eg, years)
K <- 1                ## only 1 secondary sampling occasion
n <- 25
phi <- 0.7
p <- 0.4
z <- matrix(NA, n, T)
y <- matrix(NA, n, T)
first <- rpois(n, 1)+1 ## random release dates


## ----sim-norbust2,size='footnotesize'-----------------------------------------
for(i in 1:n) {
    z[i,first[i]] <- 1  ## Known alive at release
    y[i,first[i]] <- 1
    for(t in (first[i]+1):T) {
        z[i,t] <- rbinom(1, 1, z[i,t-1]*phi) ## Alive/dead state
        y[i,t] <- rbinom(1, 1, z[i,t]*p)     ## Data
    }
}


## ----z------------------------------------------------------------------------
z[,1:7]


## ----y------------------------------------------------------------------------
y[,1:7]


## ----install-marked,include=FALSE---------------------------------------------
if(!require(marked)) install.packages("marked")
library(marked)


## ----marked-data,size='footnotesize'------------------------------------------
library(marked)
y.marked <- ifelse(is.na(y), 0, y)
cap.histories <- data.frame(
    ch=apply(y.marked, 1, paste, collapse=""))
head(cap.histories, n=3)


## ----marked-fit,size='footnotesize',results='hide',message=FALSE--------------
fm0 <- crm(data=cap.histories, model="CJS", hessian=TRUE,
           model.parameters=list(Phi=list(formula=~1),
                                 p=list(formula=~1)))


## ----marked-predict,size='footnotesize'---------------------------------------
predict(fm0)


## ----jags-nonsp,size='scriptsize'---------------------------------------------
writeLines(readLines("CJS-nonspatial.jag"))


## ----jagsUI,include=FALSE,results='hide'--------------------------------------
library(jagsUI)
library(coda)


## ----nonsp-data,size='scriptsize'---------------------------------------------
jags.data.nonsp <- list(y=y, n=n, first=first, T=T)


## ----zi,size='scriptsize'-----------------------------------------------------
zi.nonsp <- matrix(NA, n, T)
for(i in 1:n) {
    zi.nonsp[i,(first[i]+1):T] <- 1
    }
jags.inits.nonsp <- function() list(phi=runif(1), p=runif(1), z=zi.nonsp)
jags.pars.nonsp <- c("phi", "p")


## ----jc1,size='scriptsize',message=FALSE,results='hide',cache=TRUE------------
library(jagsUI)
jags.post.samples.nonsp <- jags.basic(data=jags.data.nonsp,
                                      inits=jags.inits.nonsp,
                                      parameters.to.save=jags.pars.nonsp,
                                      model.file="CJS-nonspatial.jag",
                                      n.chains=3, n.adapt=100, n.burnin=0,
                                      n.iter=2000, parallel=TRUE)


## ----jc1-plot,out.width="65%",fig.align='center',size='scriptsize'------------
plot(jags.post.samples.nonsp)


## ----sim2,size='tiny'---------------------------------------------------------
T <- 10       # primary occasions
K <- 3        # secondary occasions
n <- 25       # nIndividuals
phi <- 0.7
p0 <- 0.4     # baseline capture prob
sigma <- 0.1  # scale parameter


## ----X,size='tiny'------------------------------------------------------------
co <- seq(0.25, 0.75, length=5)
x <- cbind(rep(co, each=5),
           rep(co, times=5))
J <- nrow(x)
xlim <- ylim <- c(0,1)
s <- cbind(runif(n, xlim[1], xlim[2]),
           runif(n, ylim[1], ylim[2]))
d <- p.sp <- matrix(NA, n, J)
for(i in 1:n) {
    d[i,] <- sqrt((s[i,1]-x[,1])^2 +
                  (s[i,2]-x[,2])^2)
    p.sp[i,] <- p0*exp(-d[i,]^2/(2*sigma^2))
}


## ----z2,size='tiny'-----------------------------------------------------------
z <- matrix(NA, n, T)
y.sp <- array(NA, c(n, J, K, T)) ## 4D array
first <- rpois(n, 1)+1 ## random release dates
for(i in 1:n) {
    z[i,first[i]] <- 1
    for(t in (first[i]+1):T) {
        z[i,t] <- rbinom(1, 1,
                         z[i,t-1]*phi)
        for(j in 1:J) {
            y.sp[i,j,1:K,t] <- rbinom(K, 1,
                z[i,t]*p.sp[i,j])
        }
    }
}


## ----install-load,include=FALSE-----------------------------------------------
if(!require(remotes))
    install.packages("remotes")
library(remotes)
if(!require(openpopscr)) {
    install_github("r-glennie/openpopscr", build = TRUE,
                   build_opts = c("--no-resave-data", "--no-manual"),
                   build_vignettes = TRUE)
}
library(openpopscr)
if(!require(secr))
    install.packages("secr")
library(secr)


## ----install-openpopscr,size='scriptsize',message=FALSE,results='hide',eval=FALSE----
## library(remotes)
## install_github("r-glennie/openpopscr", build = TRUE,
##                build_opts = c("--no-resave-data", "--no-manual"),
##                build_vignettes = TRUE)
## library(openpopscr)
## library(secr)


## ----mask,size='scriptsize',out.width="60%",fig.align="center"----------------
trap.df <- data.frame(x*1000); colnames(trap.df) <- c("x","y")
traps <- read.traps(data=trap.df, detector="proximity")
mask <- make.mask(traps=traps, buffer=250)
plot(mask); points(traps, pch=3, col="blue", lwd=2)


## ----CJS-model-new,size='scriptsize'------------------------------------------
y.sp.secr <- y.sp; y.sp.secr[is.na(y.sp)] <- 0
caps <- data.frame(session=1,
                   animal=rep(slice.index(y.sp.secr, 1), y.sp.secr),
                   occasion=rep(slice.index(y.sp.secr, 3:4), y.sp.secr),
                   trap=rep(slice.index(y.sp.secr, 2), y.sp.secr))
capthist <- make.capthist(captures=caps, traps=traps, noccasions=30)


## ----format-openpop,size='scriptsize'-----------------------------------------
cjs.data <- ScrData$new(capthist, mask, primary=rep(1:10, each=3))


## ----cjs-mod-fit,size='scriptsize',results='hide',cache=TRUE------------------
mod <- CjsModel$new(list(lambda0~1, sigma~1, phi~1), cjs.data,
                    start=list(lambda0=1.5, sigma=50, phi=0.5))
mod$fit()
mod


## ----cjs-mod-est,size='scriptsize'--------------------------------------------
mod$get_par("lambda0", k = 1, j = 1)
mod$get_par("sigma", k = 1, j = 1)
mod$get_par("phi", k = 1, m=1)


## ----jags-sp,size='tiny'------------------------------------------------------
writeLines(readLines("CJS-spatial.jag"))


## ----sp-data,size='scriptsize'------------------------------------------------
jags.data.sp <- list(y=y.sp, n=n, first=first, x=x,
                     J=J, K=K, T=T, xlim=xlim, ylim=ylim)


## ----ji2,size='scriptsize'----------------------------------------------------
zi.sp <- matrix(NA, n, T)
for(i in 1:n) {
    zi.sp[i,(first[i]+1):T] <- 1
    }
jags.inits.sp <- function() list(phi=runif(1), z=zi.sp)
jags.pars.sp <- c("phi", "p0", "sigma")


## ----jc2,size='scriptsize',cache=TRUE,results='hide'--------------------------
jags.post.samples.sp <- jags.basic(data=jags.data.sp,
                                   inits=jags.inits.sp,
                                   parameters.to.save=jags.pars.sp,
                                   model.file="CJS-spatial.jag",
                                   n.chains=3, n.adapt=100, n.burnin=0,
                                   n.iter=2000, parallel=TRUE)


## ----jc2-plot,size='scriptsize',out.width="60%",fig.align='center'------------
plot(jags.post.samples.sp)


## ----move,size='tiny'---------------------------------------------------------
T <- 10; K <- 3; n <- 25 ## same as before
phi <- 0.7; p0 <- 0.4; sigma <- 0.1
## dispersal parameter (of random walk)
tau <- 0.05


## ----s,size='tiny'------------------------------------------------------------
## state-space should be bigger
xlim <- ylim <- c(0, 1)
first <- rpois(n, 1)+1
s <- array(NA, c(n, 2, T))
for(i in 1:n) {
    s[i,,first[i]] <- cbind(runif(1, xlim[1],
                                  xlim[2]),
                            runif(1, ylim[1],
                                  ylim[2]))
    for(t in (first[i]+1):T) {
        s[i,1,t] <- rnorm(1, s[i,1,t-1], tau)
        s[i,2,t] <- rnorm(1, s[i,2,t-1], tau)
    }
}
d <- p <- array(NA, c(n, J, T))
for(i in 1:n) {
    for(t in 1:T) {
        d[i,,t] <- sqrt((s[i,1,t]-x[,1])^2 +
                        (s[i,2,t]-x[,2])^2)
        p[i,,t] <- p0*exp(-d[i,,t]^2 /
                          (2*sigma^2))
    }
}


## ----z3,size='tiny'-----------------------------------------------------------
z <- matrix(NA, n, T)
y.move <- array(NA, c(n, J, K, T)) 
for(i in 1:n) {
    z[i,first[i]] <- 1
    for(t in (first[i]+1):T) {
        z[i,t] <- rbinom(1, 1, z[i,t-1]*phi)
        for(j in 1:J) {
            y.move[i,j,1:K,t] <-
                rbinom(K, 1, z[i,t]*p[i,j,t])
        }
    }
}


## ----s1,size='scriptsize',out.width="65%",fig.align="center"------------------
plot(t(s[1,,]), pch=16, type="o", xlab="x", ylab="y",
     xlim=c(0, 1), ylim=c(0, 1), asp=1, col="blue")


## ----jags-sp-move,size='tiny'-------------------------------------------------
writeLines(readLines("CJS-spatial-move.jag"))


## ----sp-move-data,size='scriptsize'-------------------------------------------
jd.sp.move <- list(y=y.move, n=n, first=first, x=x,
                   J=J, K=K, T=T, xlim=xlim, ylim=ylim)


## ----ji3,size='scriptsize'----------------------------------------------------
zi <- matrix(NA, n, T)
for(i in 1:n) {
    zi[i,(first[i]+1):T] <- 1
    }
ji.sp.move <- function() list(phi=runif(1), z=zi)
jp.sp.move <- c("phi", "p0", "sigma", "tau")


## ----jc3,size='scriptsize',results='hide',cache=TRUE--------------------------
jps.sp.move <- jags.basic(data=jd.sp.move,
                          inits=ji.sp.move,
                          parameters.to.save=jp.sp.move,
                          model.file="CJS-spatial-move.jag",
                          n.chains=3, n.adapt=100, n.burnin=0,
                          n.iter=2000, parallel=TRUE)


## ----jc3-plot,size='scriptsize',out.width="60%",fig.align='center'------------
plot(jps.sp.move[,jp.sp.move])

