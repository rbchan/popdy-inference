## ----buildit,include=FALSE,eval=FALSE-----------------------------------------
# rnw2pdf("lecture-cap-recap-spatial-II")
# rnw2pdf("lecture-cap-recap-spatial-II", tangle=TRUE)




## ----ippp1,size='scriptsize',fig.width=7.2,out.width="60%",fig.align="center",results="hide",message=FALSE----
library(raster)
elevation <- raster("elevation.tif")
delta <- res(elevation)[1]  ## resolution
plot(elevation, col=topo.colors(100), main="Elevation")


## ----ippp2,size='tiny',fig.width=7.2,out.width="60%",fig.align="center"-------
beta0 <- -15
beta1 <- 0.01 #0.005
lambda <- exp(beta0 + beta1*elevation) # Intensity function
plot(lambda, col=terrain.colors(100), main="Density surface", zlim=c(0,0.18))


## ----ippp3,size='footnotesize'------------------------------------------------
set.seed(538)  
ds <- 1                          ## Pixel area is 1 ha
lambda.values <- values(lambda)  ## Convert raster to vector
Lambda <- sum(lambda.values*ds)  ## E(N)
(N <- rpois(1, Lambda))          ## Realized N


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




## ----secr-in,warning=FALSE,size='scriptsize',results='hide',message=FALSE-----
library(secr)  
sch <- read.capthist(captfile="encounter_data_file.csv",
                   trapfile="trap_data_file.csv",
                   detector="proximity", fmt="trapID")


## ----format-elevation,size='scriptsize'---------------------------------------
elevation.xyz <- as.data.frame(elevation, xy=TRUE)
elevation.xyz$y <- round(elevation.xyz$y, 3) ## Fix numerical fuzz
elevation.xyz.m <- elevation.xyz
elevation.xyz.m$x <- elevation.xyz$x*1000  ## Convert units to meters
elevation.xyz.m$y <- elevation.xyz$y*1000
elevation.xyz.m$elevation <- scale(elevation.xyz$elevation) ## Standardize
elevation.m <- rasterFromXYZ(elevation.xyz.m)


## ----make-mask,size='scriptsize'----------------------------------------------
library(sp)
elev.sp <- as(elevation.m, "SpatialGridDataFrame")
trp <- traps(sch)
mask <- make.mask(trp, buffer=150, spacing=10)
mask <- addCovariates(mask, spatialdata=elev.sp)


## ----plot-mask,fig.height=5,out.width="95%",fig.align='center',echo=-1,size='scriptsize'----
par(mai=c(0.1, 0.1, 0.1, 0.1))
plot(mask, covariate="elevation")
plot(trp, add=TRUE)


## ----secr-elev,size='scriptsize',cache=TRUE-----------------------------------
fm.elev <- secr.fit(sch, model=list(D=~elevation, g0=~1, sigma=~1),
                    mask=mask, trace=FALSE)  ## Don't use 'buffer'
coef(fm.elev)


## ----secr-M0-real,size='scriptsize'-------------------------------------------
predict(fm.elev)


## ----secr-dsurf,size='scriptsize',fig.height=6.8,echo=-1,size='scriptsize',fig.show='hide',warning=FALSE----
dsurf <- predictDsurface(fm.elev)
dsurf.r <- raster(dsurf, covariate="D.0")
pix.area <- (delta*1000)^2
plot(dsurf.r/pix.area, col=terrain.colors(100), zlim=c(0,0.18),
     main="Estimated density surface")


## ----fxi,eval=TRUE,size='scriptsize',out.width='65%',fig.align='center'-------
fxiContour(fm.elev, i=1)


## ----regionN-M0,size='small'--------------------------------------------------
region.N(fm.elev)


## ----bugs-SCR-elev-R,size='tiny',eval=FALSE-----------------------------------
# writeLines(readLines("SCR-elev-v1.jag"))

## ----bugs-SCR-elev,size='tiny',echo=FALSE,background='beige',comment=''-------
writeLines(readLines("SCR-elev-v1.jag"))


## ----bugs-SCR-elev-v2-R,size='tiny',eval=FALSE--------------------------------
# writeLines(readLines("SCR-elev-v2.jag"))

## ----bugs-SCR-elev-v2,size='tiny',echo=FALSE,comment='',background='beige'----
writeLines(readLines("SCR-elev-v2.jag"))


## ----rjags,include=FALSE,results="hide"---------------------------------------
library(rjags)


## ----lookup,size='scriptsize',fig.show='hide',echo=-1-------------------------
#par(mai=c(0.1,0.1,0.3,0.1))  
G <- nrow(elevation.xyz)    ## nPixels  
lookup <- matrix(1:G, nrow=nrow(elevation), ncol=ncol(elevation),
                 byrow=TRUE) ## To be consistent with raster package
delta <- res(elevation)[1]  ## resolution
image(seq(delta/2, 1-delta/2, length=100),
      seq(delta/2, 1-delta/2, length=100), lookup, asp=1,
      xlab="x", ylab="y", frame=FALSE, col=0, main="Pixel ID")
text(elevation.xyz$x[seq(5, G, by=10)],
     elevation.xyz$y[seq(5, G, by=10)], seq(5, G, by=10), cex=0.4)


## ----jd-SCR-elev-faster,size='scriptsize'-------------------------------------
M <- 150  
y.tilde <- apply(y, c(1,2), sum)
n <- nrow(y)
jags.data.SCR.elev <- list(
    y.tilde=y.tilde, n=n, M=M, J=J, G=G, delta=delta,
    lookup=lookup, pixelArea=delta^2*1e4, ## Convert to hectares
    elevation=elevation.xyz$elevation, z=c(rep(1, n), rep(NA, M-n)),
    zeros=rep(0, M), K=K, zero.cap=rep(0, M), x=x, xlim=c(0,1), ylim=c(0,1))


## ----jp-ji-SCR-elev,size='scriptsize'-----------------------------------------
jp.SCR.elev <- c("beta0", "beta1", "g0", "sigma", "EN", "N")
ji.SCR.elev <- function() {
    list(z=c(rep(NA, n), rep(0,M-n)),#psi=runif(1),beta0=runif(1,-20,-15),
         beta1=0.01, s=cbind(runif(M), runif(M)),
         g0=runif(1), sigma=runif(1, 0.05, 0.1)) }


## ----mcmc-SCR-elev,size='scriptsize',results='hide',cache=TRUE,message=FALSE----
library(jagsUI)  
jags.post.SCR.elev <- jags.basic(
    data=jags.data.SCR.elev, inits=ji.SCR.elev,
    parameters.to.save=jp.SCR.elev, model.file="SCR-elev-v2.jag",
    n.chains=3, n.adapt=100, n.iter=1000, parallel=TRUE)


## ----summary-mcmc-SCR0-faster,size='tiny'-------------------------------------
summary(jags.post.SCR.elev)


## ----plot-mcmc-SCR-elev1,size='footnotesize',out.width="0.7\\textwidth",fig.align='center'----
plot(jags.post.SCR.elev[,jp.SCR.elev[1:3]])


## ----plot-mcmc-SCR-elev2,size='footnotesize',out.width="0.7\\textwidth",fig.align='center'----
plot(jags.post.SCR.elev[,jp.SCR.elev[4:6]])


## ----extract-beta,size='scriptsize'-------------------------------------------
beta.post <- as.matrix(jags.post.SCR.elev[,c("beta0","beta1")])
n.samples <- nrow(beta.post)


## ----lambda-samples,size='scriptsize'-----------------------------------------
lambda.post <- matrix(NA, n.pixels, n.samples)
for(i in 1:n.samples) {
    lambda.post[,i] <- exp(
        beta.post[i,"beta0"] +
        beta.post[i,"beta1"]*elevation.xyz$elevation)
}
lambda.post.mean <- rowMeans(lambda.post)
lambda.post.lower <- apply(lambda.post, 1, quantile, prob=0.025)
lambda.post.upper <- apply(lambda.post, 1, quantile, prob=0.975)


## ----lambda-post,size='scriptsize',fig.height=5,fig.show="hide",dev='png',dpi=400----
library(latticeExtra)
trellis.par.set(regions=list(col=hcl.colors(100))) ## Colors
panel1 <- levelplot(lambda.post.lower ~ x+y, elevation.xyz, aspect="iso",
                    at=seq(0,0.20,0.005), ## Resolution of color key
                    colorkey=list(space="bottom"),
                    xlab="", ylab="")
panel2 <- levelplot(lambda.post.mean ~x+y, elevation.xyz, aspect="iso",
                    at=seq(0,0.20,0.005))
panel3 <- levelplot(lambda.post.upper ~x+y, elevation.xyz, aspect="iso",
                    at=seq(0,0.20,0.005))
panels <- c(panel1, panel2, panel3)
panels <- update( 
    panels, scales=list(draw=FALSE), layout=c(3,1), 
    xlab="", ylab="", 
    strip=strip.custom(bg=gray(0.8),
        factor.levels=c("Lower CI","Posterior mean","Upper CI")))
plot(panels)


## ----mcmc-SCR-elev-sz,size='scriptsize',results='hide',cache=TRUE-------------
jags.post.SCR.elev.sz <- jags.basic(
    data=jags.data.SCR.elev, inits=ji.SCR.elev,
    parameters.to.save=c(jp.SCR.elev,"s","z"),
    model.file="SCR-elev-v2.jag",
    n.chains=3, n.adapt=100, n.iter=1000, parallel=TRUE)


## ----extract-sz,size='scriptsize'---------------------------------------------
s1.post <- as.matrix(jags.post.SCR.elev.sz[,paste0("s[", 1:M, ",1]")])
s2.post <- as.matrix(jags.post.SCR.elev.sz[,paste0("s[", 1:M, ",2]")])
z.post <- as.matrix(jags.post.SCR.elev.sz[,paste0("z[", 1:M, "]")])


## ----cut-s,size='scriptsize'--------------------------------------------------
lambda.r.post <- array(NA, c(sqrt(n.pixels), sqrt(n.pixels), n.samples))
for(i in 1:n.samples) {
    si1.post.d <- cut(s1.post[i,z.post[i,]==1], breaks=seq(0, 1, delta))
    si2.post.d <- cut(s2.post[i,z.post[i,]==1], breaks=seq(0, 1, delta))
    pixel.counts <- table(si2.post.d, si1.post.d)
    lambda.r.post[,,i] <- pixel.counts[sqrt(n.pixels):1,]
}


## ----lambda-r,size='scriptsize',fig.height=3,out.width='99%',echo=-1----------
par(omi=c(0,0,0,0.2))  
lambda.r.mean <- raster(apply(lambda.r.post, c(1,2), mean))
lambda.e.mean <- raster(matrix(lambda.post.mean, sqrt(n.pixels), byrow=TRUE))
lambda.re <- stack(lambda.r.mean, lambda.e.mean)
names(lambda.re) <- c("Realized", "Expected")
plot(lambda.re, zlim=c(0,0.18))


## ----si-post,size='scriptsize',echo=-(1),fig.show='hide'----------------------
par(mai=c(0.0,0.0,0.0,0.0))
plot(x, asp=1, xlim=0:1, ylim=0:1, pch=3, axes=FALSE, ann=FALSE)
points(s1.post[,1], s2.post[,1], pch=16, col=rgb(0,0,1,0.1))
contour(MASS::kde2d(s1.post[,1], s2.post[,1]), add=TRUE)
points(x[y.tilde[1,]>0,,drop=FALSE], pch=3, col="red", lwd=2, cex=2)


## ----si-post-i,size='scriptsize',echo=FALSE,fig.show='hide'-------------------
for(i in 1:n) {
    par(mai=c(0.0,0.0,0.0,0.0))
    plot(x, asp=1, xlim=0:1, ylim=0:1, pch=3, axes=FALSE, ann=FALSE)
    points(s1.post[,i], s2.post[,i], pch=16, col=rgb(0,0,1,0.1))
    contour(MASS::kde2d(s1.post[,i], s2.post[,i]), add=TRUE)
    points(x[y.tilde[i,]>0,,drop=FALSE], pch=3, col="red", lwd=2, cex=2)
}


## ----si-post-y0,size='scriptsize',echo=-(1:2),fig.show='hide'-----------------
par(mai=c(0.0,0.0,0.0,0.0))
## plot(x, asp=1, xlim=0:1, ylim=0:1, pch=3, axes=FALSE, ann=FALSE)
image(MASS::kde2d(s1.post[z.post[,150]==1,150],
                  s2.post[z.post[,150]==1,150]))
points(x, asp=1, pch=3)
##points(s1.post[z.post[,150]==1,150],
##       s2.post[z.post[,150]==1,150], pch=16, col=rgb(0,0,1,0.8))

