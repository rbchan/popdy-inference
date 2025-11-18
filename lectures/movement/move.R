## Movement models




## Movement models option 1
## Random walk movement models

T <- 50
s <- matrix(NA, nrow=T, ncol=2)
s[1,] <- c(0.5,0.5)

for(t in 2:T) {
    ## The two lines below are equivalent
##    s[t,] <- rnorm(n=2, mean=s[t-1,], sd=0.1)
    s[t,] <- s[t-1,] + rnorm(n=2, mean=c(0,0), sd=0.05)
}

s


plot(s, asp=1, xlim=c(0,1), ylim=c(0,1))
arrows(s[1:49,1], s[1:49,2],
       s[2:50,1], s[2:50,2], length=0.1, col="blue")



## Correlated random walk
T <- 50
s <- matrix(NA, nrow=T, ncol=2)
s[1,] <- c(0.5,0.5)

gamma <- 0.9
beta <- 0 #pi/2
R <- gamma*matrix(c(cos(beta), sin(beta),
                    -sin(beta), cos(beta)),
            2, 2)

s[2,] <- rnorm(2, s[1,], 0.05)
for(t in 3:T) {
    diff <- s[t-1,]-s[t-2,]
    s[t,] <- s[t-1,] + R%*%diff +
        (rnorm(n=2, mean=c(0,0), sd=0.11))
}

s


plot(s, asp=1) #, xlim=c(0,1), ylim=c(0,1))
arrows(s[1:49,1], s[1:49,2],
       s[2:50,1], s[2:50,2], length=0.1, col="blue")











set.seed(34038)
nSteps <- 500
uRW <- uCRW <- uBRW <- uBCRW <- matrix(0, nSteps, 2)
sigma <- 0.01     ## SD of white noise (epsilon)
gamma <- 0.8      ## Damping of correlated walk
theta <- 0         ## Correlation (-pi,pi)
R <- gamma*matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)), nrow=2, ncol=2)
s <- c(0,0)          ## Activity center
rho <- 0.1           ## Bias towards activity center
B <- diag(rho,2)
du2 <- du <- c(0,0)  ## Differences
for(t in 2:nSteps) {
    uRW[t,] <- uRW[t-1,]+rnorm(n=2,mean=0,sd=sigma)
    dsu <- s-uBRW[t-1,]; if(t>2) du <- uCRW[t-1,]-uCRW[t-2,]
    uCRW[t,] <- uCRW[t-1,] + R%*%du + rnorm(n=2, mean=0, sd=sigma)
    uBRW[t,] <- uBRW[t-1,] + B%*%dsu + rnorm(n=2, mean=0, sd=sigma)
    dsu2 <- s-uBCRW[t-1,]; if(t>2) du2 <- uBCRW[t-1,]-uBCRW[t-2,]
    uBCRW[t,] <- uBCRW[t-1,] + B%*%dsu2 + R%*%du2 + rnorm(n=2, mean=0, sd=sigma)
}
uall <- data.frame(Walk=c(rep("Vanilla",nSteps), rep("Correlated", nSteps),
                          rep("Biased",nSteps), rep("Biased correlated",nSteps)),
                   rbind(uRW, uCRW, uBRW, uBCRW))
uall$Walk <- factor(uall$Walk, levels=c("Vanilla", "Biased", "Correlated", "Biased correlated"))
names(uall) <- c("walk", "u1", "u2")
xyplot(u2~u1|walk, data=uall, type='l', aspect='iso',
       scales=list(relation='free',draw=FALSE), xlab="", ylab="",
       layout=c(4,1), 
       panel=function(...) {
           panel.xyplot(...)
           if(panel.number()%in%c(2,4))
               lpoints(0, 0, pch=16, col='black', ##'lightcoral',
                       cex=1)
       }, strip=strip.custom(bg=gray(0.7)), as.table=TRUE)









set.seed(430)
delta <- 0.025                          ## Resolution
grid0 <- seq(delta/2, 1-delta/2, delta)    
grid <- cbind(rep(grid0, each=length(grid0)), rep(grid0, times=length(grid0)))
distmat <- as.matrix(dist(grid))        ## Distance b/w grid points
npixels <- nrow(distmat)
V <- exp(-0.1*distmat)                  ## Covariance
R <- chol(V)
X <- t(R) %*% rnorm(npixels)            ## Gaussian random field
Xr <- rasterFromXYZ(cbind(grid, X),     ## Convert to raster
    crs="+proj=utm +zone=17 +datum=NAD83 +units=m +no_defs +type=crs")
terr <- terrain(Xr, opt=c("slope","aspect"))  ## Slope and aspect
vfieldU <- terr$slope * sin(terr$aspect)      ## x-coord of arrow head
vfieldV <- terr$slope * cos(terr$aspect)      ## y-coord of arrow head
vfieldData <- cbind(as.data.frame(Xr, xy=TRUE), v1=values(vfieldU), v2=values(vfieldV))
names(vfieldData) <- c("s1","s2","X","v1","v2")





set.seed(311409) 
uERW <- matrix(0.5, nSteps, 2)                   ## Initial location
alpha <- 0.02; D <- diag(alpha,2)                ## Effect of envi gradient
xmin <- ymin <- 0+delta; xmax <- ymax <- 1-delta ## Bounds
enviro <- numeric(nSteps)                        ## Environment along path
enviro[1] <- extract(Xr, uERW[1,,drop=FALSE])
for(t in 2:nSteps) {
    ulast <- uERW[t-1,,drop=FALSE]
    grad <- c(extract(vfieldU, ulast), extract(vfieldV, ulast))
    unew <- ulast + t(D%*%grad) + rnorm(2, 0, sigma)
    if(unew[1]<xmin) unew[1] <- xmin+(xmin-unew[1]) ## Reflection at boundary
    if(unew[1]>xmax) unew[1] <- xmax-(unew[1]-xmax)
    if(unew[2]<ymin) unew[2] <- ymin+(ymin-unew[2])
    if(unew[2]>ymax) unew[2] <- ymax-(unew[2]-ymax)
    uERW[t,] <- unew                   ## New location
    enviro[t] <- extract(Xr, unew)     ## Environment at new location
}



colrs <- colorRampPalette(c("#333a7b","#4b6982","#70c6c7","#b4ffd8","#ffffff"), alpha=0.9)
levPlot <- levelplot(X~s1+s2, vfieldData, aspect="iso", col.regions=colrs(100),
                     colorkey=list(space='top', at=seq(min(X),max(X),len=100), width=0.5),
                     xlab="", ylab="", scales=list(draw=FALSE),
##                     par.settings=list(axis.line = list(col = "transparent")),
                     panel=function(...) {
                         panel.levelplot(...)
                         with(vfieldData, larrows(s1, s2, s1+v1*0.02, s2+v2*0.02,
                                                  length=0.02, col=rgb(0,0,0,0.5)))
                         llines(uERW[,1], uERW[,2], col=rgb(1,0,0,0.9), lwd=2)
                     })
wirePlot <- wireframe(X~s1+s2, vfieldData, shade=FALSE, drape=TRUE, col.regions=colrs(100),
                      xlab="",ylab="",zlab="", scales=list(draw=FALSE),
                      x2=uERW[,1], y2=uERW[,2], z2=enviro,
                      col=rgb(0,0,0,0.1),
                      panel.3d.wireframe=function(x,y,z, xlim, ylim, zlim, col,
                                                  xlim.scaled,ylim.scaled,zlim.scaled,
                                                  x2,y2,z2,...) {
                          panel.3dwire(x=x, y=y, z=z, xlim=xlim, ylim=ylim, zlim=zlim,
                                       xlim.scaled=xlim.scaled, ylim.scaled=ylim.scaled,
                                       zlim.scaled=zlim.scaled, col=col, ...)
                          xx <- xlim.scaled[1]+diff(xlim.scaled)*(x2-xlim[1])/diff(xlim)
                          yy <- ylim.scaled[1]+diff(ylim.scaled)*(y2-ylim[1])/diff(ylim)
                          zz <- zlim.scaled[1]+diff(zlim.scaled)*(z2-zlim[1])/diff(zlim)
                          panel.3dscatter(x=xx, y=yy, z=zz, type='l', col=rgb(1,0,0,0.9),
                                          lwd=2, ...) 
                      }, screen=list(z=20,x=-20), zoom=1.01)
plot(c(levPlot,wirePlot))













## Make a raster to use as a spatial covariate
delta <- 0.025 ## Resolution
grid0 <- seq(delta/2, 1-delta/2, delta)
grid <- cbind(rep(grid0, each=length(grid0)),
              rep(grid0, times=length(grid0)))

plot(grid, asp=1, pch=3)

distmat <- as.matrix(dist(grid))

npixels <- nrow(distmat)

## Covariance matrix
V <- exp(-0.1*distmat)
R <- chol(V)

## This is a draw from a multivariate normal distribution
X <- t(R) %*% rnorm(npixels)

library(lattice)

levelplot(X ~ grid[,1]+grid[,2], aspect="iso")















## Movement model option 1, aka:
## Step-selection movement model
## Point process movement model






set.seed(3490)
envData <- as.data.frame(Xr, xy=TRUE); colnames(envData) <- c("u1","u2","x")
J <- nrow(envData); nSteps <- 100
distmatX <- as.matrix(dist(envData[,c("u1","u2")]))
beta1 <- 2; beta2 <- -15
s <- 821      ## Activity center is pixel 821
dist2s <- distmatX[s,]
alpha1 <- -50
rcat <- function(pvec) which(rmultinom(n=1, size=1, prob=pvec)==1)
f <- exp(beta1*envData$x + beta2*dist2s)   ## Relative probability of use
loc <- integer(nSteps); loc[1] <- 128      ## Initial location
for(t in 2:nSteps) { ## 
    g <- exp(alpha1*distmat[loc[t-1],])
    fg <- f*g             ## Unnormalized RSF
    rsf <- fg/sum(f*g)    ## Pr(selecting pixel j | u_t)
    loc[t] <- rcat(rsf)
}
ucoords <- envData[loc,c("u1","u2")]
levelplot(x~u1+u2, envData, aspect="iso", scales=list(draw=FALSE),
          region=FALSE, contour=TRUE, xlab="", ylab="", col=rgb(0,0,0,0.2),
          panel=function(...) {
              panel.levelplot(...)
              larrows(ucoords[-nSteps,1], ucoords[-nSteps,2],
                      ucoords[-1,1], ucoords[-1,2], length=0.05,
                      col=rgb(0,0,0,0.9))
              lpoints(envData[s,c("u1","u2"),drop=FALSE], pch=16,
                      col='blue', cex=1.5)
          })









T <- 50
pixelID <- rep(NA, T)

pixelID[1] <- 256

beta1 <- -10
beta2 <- 5

for(t in 2:T) {
    h <- exp(beta1*distmat[,pixelID[t-1]] + beta2*X)
    theta <- h/sum(h)
    pixelID[t] <- which(rmultinom(n=1, size=1, prob=theta)==1)
}

pixelID

path <- grid[pixelID,]

Xmat <- matrix(X, sqrt(npixels), sqrt(npixels), byrow=TRUE)

plot(path, asp=1, xlim=c(0,1), ylim=c(0,1))
image(Xmat, add=TRUE)
arrows(path[1:49,1], path[1:49,2],
       path[2:50,1], path[2:50,2], length=0.1, col="yellow")


