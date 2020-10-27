## ----buildit,include=FALSE,eval=FALSE-----------------------------------------
## rnw2pdf("lecture-cap-recap-spatial")
## rnw2pdf("lecture-cap-recap-spatial", tangle=TRUE)




## ----ppp1,size='footnotesize',out.width="50%",fig.align="center"--------------
lambda1 <- 25; A <- 1     ## lambda1=density, A=area
N <- rpois(1, lambda1*A)  
s <- cbind(runif(N, 0, 1), runif(N, 0, 1))
plot(s, pch=16, col="blue", xlab="Easting", ylab="Northing",
     xlim=c(0,1), ylim=c(0,1), cex.lab=1.5, asp=1)


## ----include=FALSE,eval=FALSE-------------------------------------------------
## ## Simulate covariate from normal distribution with exponential covariance
## delta <- 0.01
## grid0 <- seq(0+delta/2, 1-delta/2, delta)
## grid <- cbind(rep(grid0, each=length(grid0)),
##               rep(grid0, times=length(grid0)))
## dist.grid <- as.matrix(dist(grid))
## V <- exp(-dist.grid*10)
## set.seed(4984)
## w <- t(chol(V)) %*% rnorm(nrow(V))
## 
## elevation <- 1000 + w*100
## 
## library(raster)
## elevation.r <- rasterFromXYZ(cbind(grid[,1], grid[,2], elevation))
## plot(elevation.r)
## writeRaster(elevation.r, filename="elevation.tif", overwrite=TRUE)
## 
## 
## 
## ## library(lattice)
## 
## ## levelplot(w ~ grid[,1]+grid[,2], aspect="iso", at=seq(min(w), max(w), len=100))
## 
## ## lambda <- exp(-7 + 1.5*w)
## ## (Lambda <- sum(lambda*delta^2*1e4))
## ## (N <- rpois(1, Lambda))
## 
## ## levelplot(lambda ~ grid[,1]+grid[,2], aspect="iso",
## ##           at=seq(min(lambda), max(lambda), len=100),
## ##           col.regions=topo.colors(100))
## 
## ## pi <- lambda/sum(lambda)
## 
## ## s <- grid[sample(nrow(grid), size=N, replace=TRUE, prob=pi),]
## 
## ## levelplot(lambda ~ grid[,1]+grid[,2], aspect="iso",
## ##           panel=function(...) {
## ##               panel.levelplot(...)
## ##               lpoints(s, pch=16, col="blue")
## ##           },
## ##           at=seq(min(lambda), max(lambda), len=100),
## ##           col.regions=terrain.colors(100))
## 
## 


## ----ippp1,size='footnotesize',fig.width=7.2,out.width="60%",fig.align="center",results="hide"----
library(raster)
elevation <- raster("elevation.tif")
plot(elevation, col=topo.colors(100), main="Elevation")


## ----ippp2,size='footnotesize',fig.width=7.2,out.width="60%",fig.align="center"----
beta0 <- -15
beta1 <- 0.01 #0.005
lambda <- exp(beta0 + beta1*elevation) # Intensity function
plot(lambda, col=terrain.colors(100), main="Density surface")


## ----ippp3,size='footnotesize'------------------------------------------------
set.seed(538)  
ds <- 1                            ## Pixel area is 1 ha
lambda.values <- values(lambda)    ## Convert raster to vector
Lambda <- sum(lambda.values*ds)    ## E(N)
(N <- rpois(1, Lambda))            ## Realized N


## ----ipp4,size='footnotesize'-------------------------------------------------
n.pixels <- length(lambda)
jitter <- 0.005                    ## Half width of pixel 
s.pixels <- sample(n.pixels, size=N, replace=TRUE,
                   prob=lambda.values/Lambda)
elevation.xyz <- as.data.frame(elevation, xy=TRUE)
s <- elevation.xyz[s.pixels,c("x","y")] +
    cbind(runif(N, -jitter, jitter),runif(N, -jitter, jitter))


## ----ippp5,size='scriptsize',fig.width=7.2,out.width="70%",fig.align="center"----
plot(lambda, col=terrain.colors(100),
     main="Density surface with activity centers")
points(s, pch=16, cex=1, col="blue")


## ----traps1,size='scriptsize',fig.width=7.2,out.width="60%",fig.align="center"----
x <- cbind(rep(seq(0.15, 0.85, by=0.1), each=8),
           rep(seq(0.15, 0.85, by=0.1), times=8))  ## Trap locations
plot(lambda, col=terrain.colors(100),
     main="Density surface with activity centers and traps")
points(s, pch=16, col="blue") ## Activity center locations
points(x, pch=3)              ## Trap locations


## ----dist1,size='footnotesize'------------------------------------------------
J <- nrow(x)                 ## nTraps
dist.sx <- matrix(NA, N, J)  
for(i in 1:N) {
    dist.sx[i,] <- sqrt((s[i,1]-x[,1])^2 + (s[i,2]-x[,2])^2)
}


## ----dist2,size='footnotesize'------------------------------------------------
dist.sx[1:4,1:5]


## ----p1,size='footnotesize'---------------------------------------------------
g0 <- 0.2
sigma <- 0.05
p <- g0*exp(-dist.sx^2/(2*sigma^2))


## ----p2,size='footnotesize'---------------------------------------------------
print(p[1:4,1:5], digits=3)


## ----y1,size='footnotesize'---------------------------------------------------
K <- 5                          # nOccasions
y.all <- array(NA, c(N, J, K))
for(i in 1:N) {
    for(j in 1:J) {
        y.all[i,j,] <- rbinom(K, 1, prob=p[i,j])
    }
}


## ----y2,size='footnotesize'---------------------------------------------------
captured <- rowSums(y.all)>0
y <- y.all[captured,,]


## ----y3,size='footnotesize'---------------------------------------------------
y[1:2,1:5,1]


## ----n,size='footnotesize'----------------------------------------------------
(n <- nrow(y))


## ----cap-freq,size='footnotesize'---------------------------------------------
y.tilde <- rowSums(y)
table(y.tilde)


## ----sp-cap-freq,size='footnotesize'------------------------------------------
y.nok <- apply(y, c(1, 2), sum)
y.nojk <- apply(y.nok>0, 1, sum)
table(y.nojk)


## ----spider,size='scriptsize',fig.width=7.15,fig.show='hide'------------------
plot(lambda, col=terrain.colors(100),
     main="Density surface, activity centers, traps, and capture locs")
s.cap <- s[captured,]
for(i in 1:n) {
    traps.i <- which(rowSums(y[i,,])>0)
    for(j in 1:length(traps.i)) {
        segments(s.cap[i,1], s.cap[i,2],
                 x[traps.i[j],1], x[traps.i[j],2], col=gray(0.9))
    }
}
points(s[captured,], pch=16, col="blue") ## Activity center locations
points(s[!captured,], pch=1, col="blue") ## Activity center locations
points(x, pch=3)              ## Trap locations


## ----write,include=FALSE,results="hide"---------------------------------------
ch.out <- data.frame(session=1,
                     individual=rep(slice.index(y, 1), y),
                     occasion=rep(slice.index(y, 3), y),
                     trap=rep(slice.index(y, 2), y))
write.table(ch.out, file="encounter_data_file.csv",
            row.names=FALSE, col.names=FALSE, sep=",")
traps.out <- data.frame(trap=1:nrow(x), x*1000)
write.table(traps.out, file="trap_data_file.csv",
            row.names=FALSE, col.names=FALSE, sep=",")
library(secr)


## ----secr-in,warning=FALSE,size='tiny'----------------------------------------
library(secr)  
sch <- read.capthist(captfile="encounter_data_file.csv",
                     trapfile="trap_data_file.csv",
                     detector="multi", fmt="trapID")
summary(sch)


## ----secr-plot,out.width="70%",fig.align="center"-----------------------------
plot(sch)


## ----secr-M0,size='scriptsize',cache=TRUE-------------------------------------
fm.M0 <- secr.fit(sch, model=list(D=~1, g0=~1, sigma=~1),
                  buffer=150, trace=FALSE)
coef(fm.M0)


## ----secr-M0-real,size='scriptsize'-------------------------------------------
predict(fm.M0)


## ----secr-Mt,size='scriptsize',cache=TRUE-------------------------------------
fm.Mt <- secr.fit(sch, model=list(D=~1, g0=~t, sigma=~1),
                  buffer=150, trace=FALSE, ncores=3)
coef(fm.Mt)


## ----secr-Mt-real,size='scriptsize'-------------------------------------------
predict(fm.Mt)


## ----secr-Mb,size='scriptsize',cache=TRUE-------------------------------------
fm.Mb <- secr.fit(sch, model=list(D=~1, g0=~b, sigma=~1),
                  buffer=150, trace=FALSE, ncores=3)
coef(fm.Mb)


## ----secr-Mb-real,size='scriptsize'-------------------------------------------
predict(fm.Mb)


## ----regionN-M0,size='footnotesize'-------------------------------------------
region.N(fm.M0)


## ----regionN-Mt,size='footnotesize'-------------------------------------------
region.N(fm.Mt)


## ----regionN-Mb,size='footnotesize'-------------------------------------------
region.N(fm.Mb)


## ----buffer1,size='scriptsize',warning=FALSE,cache=TRUE-----------------------
predict(update(fm.M0, buffer=100))[1,]


## ----buffer2,size='scriptsize',warning=FALSE,cache=TRUE-----------------------
predict(update(fm.M0, buffer=150))[1,]


## ----buffer3,size='scriptsize',warning=FALSE,cache=TRUE-----------------------
predict(update(fm.M0, buffer=200))[1,]


## ----buffer4,size='scriptsize',warning=FALSE,cache=TRUE-----------------------
predict(update(fm.M0, buffer=250))[1,]


## ----buffer5,size='scriptsize',warning=FALSE,cache=TRUE-----------------------
predict(update(fm.M0, buffer=300))[1,]


## ----aic,size='tiny'----------------------------------------------------------
AIC(fm.M0, fm.Mt, fm.Mb)


## ----bugs-SC0,size='scriptsize'-----------------------------------------------
writeLines(readLines("SCR0.jag"))


## ----rjags,include=FALSE,results="hide"---------------------------------------
library(rjags)


## ----jd-SCR0-aug,size='scriptsize'--------------------------------------------
M <- 150
y.aug <- array(0, c(M, J, K))
y.aug[1:nrow(y),,] <- y
jags.data.SCR0 <- list(y=y.aug, M=M, J=J, K=K,
                       x=x, xlim=c(0,1), ylim=c(0,1))


## ----ji-M0-aug,size='scriptsize'----------------------------------------------
ji.SCR0 <- function() {
    list(z=rep(1,M), psi=runif(1),
         s=cbind(runif(M), runif(M)),
         g0=runif(1), sigma=runif(1, 0.05, 0.1)) }
jp.SCR0 <- c("g0", "sigma", "EN", "N")
library(jagsUI)


## ----mcmc-M0-aug,size='scriptsize',results='hide',cache=TRUE------------------
jags.post.SCR0 <- jags.basic(data=jags.data.SCR0, inits=ji.SCR0,
                             parameters.to.save=jp.SCR0,
                             model.file="SCR0.jag",
                             n.chains=3, n.adapt=100, n.burnin=0,
                             n.iter=2000, parallel=TRUE)


## ----summary-mcmc-SCR0,size='tiny'--------------------------------------------
summary(jags.post.SCR0)


## ----plot-mcmc-SCR0,size='footnotesize',out.width="0.7\\textwidth",fig.align='center'----
plot(jags.post.SCR0[,jp.SCR0])


## ----bugs-SC0-faster,size='tiny'----------------------------------------------
writeLines(readLines("SCR0-faster.jag"))


## ----jd-SCR0-aug-faster,size='scriptsize'-------------------------------------
y.tilde <- apply(y, c(1,2), sum)
n <- nrow(y)
jags.data.SCR0.faster <- list(y.tilde=y.tilde, n=n, M=M, J=J,
                              z=c(rep(1, n), rep(NA, M-n)),
                              K=K, zero=rep(0, M), x=x,
                              xlim=c(0,1), ylim=c(0,1))


## ----ji-M0-aug-faster,size='scriptsize'---------------------------------------
ji.SCR0.faster <- function() {
    list(z=c(rep(NA, n), rep(0,M-n)), psi=runif(1),
         s=cbind(runif(M), runif(M)),
         g0=runif(1), sigma=runif(1, 0.05, 0.1)) }


## ----mcmc-M0-aug-faster,size='scriptsize',results='hide',cache=TRUE-----------
jags.post.SCR0.faster <- jags.basic(data=jags.data.SCR0.faster,
                                    inits=ji.SCR0.faster,
                                    parameters.to.save=jp.SCR0,
                                    model.file="SCR0-faster.jag",
                                    n.chains=3, n.adapt=100, n.burnin=0,
                                    n.iter=2000, parallel=TRUE)


## ----summary-mcmc-SCR0-faster,size='tiny'-------------------------------------
summary(jags.post.SCR0.faster)


## ----plot-mcmc-SCR0-faster,size='footnotesize',out.width="0.7\\textwidth",fig.align='center'----
plot(jags.post.SCR0.faster[,jp.SCR0])

