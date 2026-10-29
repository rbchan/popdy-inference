#' ---
#' format:
#'   revealjs:
#'     theme: default
#'     pagetitle: Cormack-Jolly-Seber Model
#'     slide-number: true
#'     echo: true
#' lang: en
#' title: lecture-CJS
#' ---
#' 

#| label: sim-norobust
set.seed(34918)  
T <- 10               ## primary periods (eg, years)
K <- 1                ## only 1 secondary sampling occasion
n <- 25
phi <- 0.7
p <- 0.4
z <- matrix(NA, n, T)
y <- matrix(NA, n, T)
first <- rpois(n, 1)+1 ## random release dates

#| label: sim-norbust2
for(i in 1:n) {
    z[i,first[i]] <- 1  ## Known alive at release
    y[i,first[i]] <- 1
    for(t in (first[i]+1):T) {
        z[i,t] <- rbinom(1, 1, z[i,t-1]*phi) ## Alive/dead state
        y[i,t] <- rbinom(1, 1, z[i,t]*p)     ## Data
    }
}

#| label: z
z[,1:7]

#| label: y-peak
y[,1:7]

#| label: install-marked
#| include: false
#| eval: true
if(!require(marked)) install.packages("marked")
library(marked)

#| label: marked-data
library(marked)
y.marked <- ifelse(is.na(y), 0, y)
cap.histories <- data.frame(
    ch=apply(y.marked, 1, paste, collapse=""))
head(cap.histories, n=3)

#| label: marked-fit
#| results: hide
#| warning: false
fm0 <- crm(data=cap.histories, model="CJS", hessian=TRUE,
           model.parameters=list(Phi=list(formula=~1),
                                 p=list(formula=~1)))

#| label: marked-predict
predict(fm0)

#| label: jags-nonsp
#| background: beige
#| comment: ''
#| echo: false
writeLines(readLines("CJS-nonspatial.jag"))

#| label: jagsUI
#| include: false
#| results: hide
library(jagsUI)
library(coda)

#| label: nonsp-data
jags.data.nonsp <- list(y=y, n=n, first=first, T=T)

#| label: zi
zi.nonsp <- matrix(NA, n, T)
for(i in 1:n) {
    zi.nonsp[i,(first[i]+1):T] <- 1
    }
jags.inits.nonsp <- function() list(phi=runif(1), p=runif(1), z=zi.nonsp)
jags.pars.nonsp <- c("phi", "p")

#| label: jc1
#| results: hide
#| message: false
#| cache: true
library(jagsUI)
jags.post.samples.nonsp <- jags.basic(data=jags.data.nonsp,
                                      inits=jags.inits.nonsp,
                                      parameters.to.save=jags.pars.nonsp,
                                      model.file="CJS-nonspatial.jag",
                                      n.chains=3, n.adapt=100, n.burnin=0,
                                      n.iter=2000, parallel=TRUE)

#| label: jc1-plot
#| out-width: 65%
#| fig-align: center
plot(jags.post.samples.nonsp)

#| label: sim2
T <- 10       # primary occasions
K <- 3        # secondary occasions
n <- 25       # nIndividuals
phi <- 0.7
p0 <- 0.4     # baseline capture prob
sigma <- 0.1  # scale parameter

#| label: X
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

#| label: z2
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

#| label: install-openpopscr
#| eval: false
#| results: hide
#| message: false
# library(remotes)
# install_github("r-glennie/openpopscr", build = TRUE,
#                build_opts = c("--no-resave-data", "--no-manual"),
#                build_vignettes = TRUE)
# library(openpopscr)
# library(secr)

#| label: load
#| include: false
library(secr)
library(openpopscr)

#| label: mask
#| out-width: 60%
#| fig-align: center
#| fig-height: 7
#| fig-width: 7
trap.df <- data.frame(x*1000); colnames(trap.df) <- c("x","y")
traps <- read.traps(data=trap.df, detector="proximity")
mask <- make.mask(traps=traps, buffer=250)
plot(mask); points(traps, pch=3, col="blue", lwd=2)

#| label: CJS-model-new
y.sp.secr <- y.sp
y.sp.secr[is.na(y.sp)] <- 0
caps <- data.frame(session=1,
                   animal=rep(slice.index(y.sp.secr, 1), y.sp.secr),
                   occasion=rep(slice.index(y.sp.secr, 3:4), y.sp.secr),
                   trap=rep(slice.index(y.sp.secr, 2), y.sp.secr))
capthist <- make.capthist(captures=caps, traps=traps, noccasions=30)

#| label: format-openpop
cjs.data <- ScrData$new(capthist, mask, primary=rep(1:10, each=3))

#| label: cjs-mod-fit
#| results: hide
#| cache: true
mod <- CjsModel$new(list(lambda0~1, sigma~1, phi~1), cjs.data,
                    start=list(lambda0=1.5, sigma=50, phi=0.5))
mod$fit()
mod

#| label: cjs-mod-est
mod$get_par("lambda0", k = 1, j = 1)
mod$get_par("sigma", k = 1, j = 1)
mod$get_par("phi", k = 1, m=1)

#| label: jags-sp
#| echo: false
writeLines(readLines("CJS-spatial.jag"))

#| label: sp-data
jags.data.sp <- list(y=y.sp, n=n, first=first, x=x,
                     J=J, K=K, T=T, xlim=xlim, ylim=ylim)

#| label: ji2
zi.sp <- matrix(NA, n, T)
for(i in 1:n) {
    zi.sp[i,(first[i]+1):T] <- 1
    }
jags.inits.sp <- function() list(phi=runif(1), z=zi.sp)
jags.pars.sp <- c("phi", "p0", "sigma")

#| label: jc2
#| results: hide
#| cache: true
jags.post.samples.sp <- jags.basic(data=jags.data.sp,
                                   inits=jags.inits.sp,
                                   parameters.to.save=jags.pars.sp,
                                   model.file="CJS-spatial.jag",
                                   n.chains=3, n.adapt=100, n.burnin=0,
                                   n.iter=2000, parallel=TRUE)

#| label: jc2-plot
#| out-width: 60%
#| fig-align: center
plot(jags.post.samples.sp)

#| label: move
T <- 10; K <- 5; n <- 25 ## same as before
phi <- 0.9; p0 <- 0.4; sigma <- 0.05
## dispersal parameter (of random walk)
tau <- 0.05

#| label: s
## state-space should be bigger
xlim <- ylim <- c(0, 1)
first <- rpois(n, 1)+1
s <- array(NA, c(n, 2, T))
for(i in 1:n) {
    s[i,,first[i]] <-
        cbind(runif(1, xlim[1], xlim[2]),
              runif(1, ylim[1], ylim[2]))
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

#| label: z3
z <- matrix(NA, n, T)
y.move <- array(NA, c(n, J, T)) 
for(i in 1:n) {
    z[i,first[i]] <- 1
    for(t in (first[i]+1):T) {
        z[i,t] <- rbinom(1, 1, z[i,t-1]*phi)
        for(j in 1:J) {
            y.move[i,j,t] <-
                rbinom(1, K, z[i,t]*p[i,j,t])
        }
    }
}

#| label: s1
#| out-width: 65%
#| fig-align: center
#| fig-height: 7
#| fig-width: 7
#| code-fold: true
plot(t(s[1,,]), pch=16, type="o", xlab="x", ylab="y",
     xlim=c(0, 1), ylim=c(0, 1), asp=1, col="blue")
points(x, pch=3)

#| label: s-all
#| echo: false
#| out-width: 65%
#| fig-align: center
par(mai=c(0.1,0.1,0.1,0.1))  
blue <- seq(0.1, 0.9, length=nrow(s))
plot(t(s[1,,]), pch=1, type="o", xaxt="n", yaxt="n", frame=FALSE, ann=FALSE,
     xlim=c(-0.2, 1.2), ylim=c(-0.2, 1.2), asp=1, col=rgb(0,0,blue[1],0.5))
for(i in 2:nrow(s)) lines(t(s[i,,]), pch=i, type="o", col=rgb(0,0,blue[i],0.5))
points(x, pch=3)

#| label: jags-sp-move
#| background: beige
#| comment: ''
#| echo: false
writeLines(readLines("CJS-spatial-move.jag"))

#| label: sp-move-data
jd.sp.move <- list(y=y.move, n=n, first=first, x=x, J=J, K=K, T=T, xlim=xlim, ylim=ylim)

#| label: ji3
zi <- matrix(NA, n, T)
for(i in 1:n) {
  zi[i,(first[i]+1):T] <- 1
}
ji.sp.move <- function() list(phi=runif(1), z=zi, s=s) ## Cheating with s
jp.sp.move <- c("phi", "p0", "sigma", "tau")
good1 <- which.max(rowSums(z, na.rm=TRUE)+rowSums(y.move, na.rm=TRUE))
time1 <- first[good1]:T
jp.s1x <- paste0("s[", good1, ",1,", time1, "]")
jp.s1y <- paste0("s[", good1, ",2,", time1, "]")
jp.z1 <- paste0("z[", good1, ",", time1, "]") 

#| label: jc3
#| results: hide
#| cache: true
jps.sp.move <- jags.basic(data=jd.sp.move, inits=ji.sp.move,
                          parameters.to.save=c(jp.sp.move, jp.s1x, jp.s1y, jp.z1),
                          model.file="CJS-spatial-move.jag", n.chains=3, n.adapt=100,
                          n.burnin=0, n.iter=2000, parallel=TRUE)

#| label: jc3-plot
#| echo: false
#| out-width: 60%
#| fig-align: center
#| fig-height: 8
#| fig-width: 7
plot(jps.sp.move[,jp.sp.move])

#| label: jc3-s
#| echo: false
#| out-width: 90%
#| fig-align: center
#| fig-height: 7
#| fig-width: 7
s1.post <- as.matrix(window(jps.sp.move[,c(jp.s1x,jp.s1y)], start=1501, thin=10))
z1.post <- as.matrix(window(jps.sp.move[,jp.z1], start=1501, thin=10))
plot(x, pch=3, xlim=c(0,1), ylim=c(0,1), asp=1, ann=FALSE, axes=FALSE)
for(i in 1:nrow(s1.post)) {
    alive1 <- which(z1.post[i,]==1)
    T1 <- length(alive1)
    arrows(s1.post[i,alive1[-T1]], s1.post[i,(alive1+T1)[-T1]],
           s1.post[i,alive1[-1]], s1.post[i,(alive1+T1)[-1]],
           col=rgb(0,0,1,0.1), length=0.05)
}
traps1 <- which(rowSums(y.move[good1,,], na.rm=TRUE)>0)
points(x[traps1,,drop=FALSE], pch=3, cex=2, lwd=2, col="red")
lines(t(s[good1,,]), type="o", pch=16, col="blue")

