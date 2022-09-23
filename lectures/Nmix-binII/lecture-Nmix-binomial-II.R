## ----buildit,include=FALSE,eval=FALSE-----------------------------------------
## rnw2pdf("lecture-Nmix-binomial-II")
## rnw2pdf("lecture-Nmix-binomial-II", tangle=TRUE)




## ----grouse-in,size='footnotesize',cache=FALSE,results='hide',message=FALSE----
library(unmarked)
grouse.data <- read.csv("grouse_data_Nmix.csv", row.names=1)
grouse.umf <- unmarkedFramePCount(
    y=grouse.data[,paste0("grouse",1:3)],
    siteCovs=grouse.data[,c("utmE","utmN","elevation")],
    obsCovs=list(temp=grouse.data[,paste0("Temperature.",1:3)]))


## ----grouse-stand,size='footnotesize',cache=FALSE-----------------------------
site.covs.s <- scale(siteCovs(grouse.umf))
colnames(site.covs.s) <- paste0(colnames(site.covs.s), ".s")
siteCovs(grouse.umf) <- cbind(siteCovs(grouse.umf), site.covs.s)
obsCovs(grouse.umf) <- scale(obsCovs(grouse.umf))


## ----grouse-mods,size='footnotesize',warning=FALSE,cache=FALSE----------------
fm1 <- pcount(~temp ~ elevation.s+utmE.s+utmN.s, grouse.umf, K=50)
fm2 <- pcount(~temp ~ elevation.s+utmN.s, grouse.umf, K=50)
fm3 <- pcount(~temp ~ elevation.s, grouse.umf, K=50)
fm4 <- pcount(~1 ~ elevation.s+utmN.s, grouse.umf, K=50)
fm5 <- pcount(~1 ~ elevation.s, grouse.umf, K=50)
fm6 <- pcount(~1 ~ 1, grouse.umf, K=50)


## ----grouse-fitlist,size='footnotesize',warning=FALSE-------------------------
grouse.models <- fitList('lam(elev+utmE+utmN)p(temp)'=fm1,
                         'lam(elev+utmN)p(temp)'=fm2,
                         'lam(elev)p(ptemp)'=fm3,
                         'lam(elev+utmN)p(.)'=fm4,
                         'lam(elev)p(.)'=fm5,
                         'lam(.)p(.)'=fm6)



## ----addNA,size='footnotesize'------------------------------------------------
na.sites <- apply(is.na(site.covs.s), 1, any)
grouse.counts <- getY(grouse.umf)
grouse.counts[na.sites,] <- NA
grouse.umf@y <- grouse.counts


## ----grouse-mods2,size='scriptsize',warning=FALSE,cache=TRUE------------------
fm1 <- pcount(~temp ~ elevation.s+utmE.s+utmN.s, grouse.umf, K=50)
fm2 <- pcount(~temp ~ elevation.s+utmN.s, grouse.umf, K=50)
fm3 <- pcount(~temp ~ elevation.s, grouse.umf, K=50)
fm4 <- pcount(~1 ~ elevation.s+utmN.s, grouse.umf, K=50)
fm5 <- pcount(~1 ~ elevation.s, grouse.umf, K=50)
fm6 <- pcount(~1 ~ 1, grouse.umf, K=50)


## ----grouse-fitlist2,size='scriptsize',warning=FALSE--------------------------
grouse.models <- fitList('lam(elev+utmE+utmN)p(temp)'=fm1,
                         'lam(elev+utmN)p(temp)'=fm2,
                         'lam(elev)p(ptemp)'=fm3,
                         'lam(elev+utmN)p(.)'=fm4,
                         'lam(elev)p(.)'=fm5,
                         'lam(.)p(.)'=fm6)



## ----aic-table,size='scriptsize'----------------------------------------------
modSel(grouse.models)


## ----elev-sp,fig.width=9.7,out.width='80%',fig.align='center',size='scriptsize',results='hide',message=FALSE,warning=FALSE----
library(raster); library(sf)
load("state_boundaries.gzip")
elev <- raster("elev_utm16.tif")
plot(elev, main="Elevation")
plot(ga.nc.sc.tn.utm, add=TRUE)


## ----utmN-sp,fig.width=9.7,out.width='80%',fig.align='center',size='scriptsize'----
utmN <- elev
utmN[] <- yFromCell(elev, cell=1:length(elev))
plot(utmN, main="UTM northing")
plot(ga.nc.sc.tn.utm, add=TRUE)


## ----means-sds,size='scriptsize'----------------------------------------------
attributes(site.covs.s)[3] ## Means
attributes(site.covs.s)[4] ## SDs


## ----std-rast,size='scriptsize'-----------------------------------------------
elev.s <- (elev-666.1)/186.6        
utmN.s <- (utmN-3852681.2)/18030.2  


## ----stack,size='scriptsize'--------------------------------------------------
elev.utm <- stack(elev.s,utmN.s)
names(elev.utm) <- c("elevation.s", "utmN.s")


## ----lam-pred-map,size='scriptsize',results='hide',fig.width=10.5,out.width='85%',fig.align='center',cache=TRUE----
grouse.pred.map <- predict(fm2, newdata=elev.utm, type="state")
plot(grouse.pred.map)


## ----Elam-pred-map,size='scriptsize',results='hide',fig.width=9.7,out.width='90%',fig.align='center'----
plot(grouse.pred.map$Predicted, main="Grouse distribution")
plot(ga.nc.sc.tn.utm, add=TRUE)


## ----bugs-data,size='scriptsize'----------------------------------------------
jags.data <- list(
    y=grouse.counts,
    elevation=site.covs.s[,"elevation.s"],
    utmE=site.covs.s[,"utmE.s"],
    utmN=site.covs.s[,"utmN.s"],
    temp=as.matrix(grouse.data[,paste0("Temperature.",1:3)]),
    nSites=nrow(grouse.counts),
    nOccasions=ncol(grouse.counts))
jags.data$temp <- (jags.data$temp-mean(jags.data$temp, na.rm=TRUE))/
    sd(jags.data$temp, na.rm=TRUE) # standardize temperature


## ----bugs-inits,size='scriptsize'---------------------------------------------
jags.inits <- function() {
    list(lambda.intercept=runif(1), alpha0=rnorm(1),
         N=rep(2, jags.data$nSites))
}


## ----bugs-pars,size='scriptsize'----------------------------------------------
jags.pars <- c("beta0", "beta1", "beta2", "beta3",
               "alpha0", "alpha1", "totalAbundance",
               "ld.y.dot", "ld.ydot.N")


## ----bugs,size='tiny',echo=FALSE,comment='',background='lightblue'------------
writeLines(readLines("Nmix-model-grouse.jag"))

## ----jagsUI,include=FALSE-----------------------------------------------------
library(jagsUI)


## ----bugs-mcmc,size='footnotesize',message=FALSE,cache=TRUE,results='hide'----
library(jagsUI)
jags.data1 <- jags.data
jags.data1$modswitch <- c(1,1,1,1) ## Include all covariates
jm1 <- jags.basic(data=jags.data1, inits=jags.inits,
                  parameters.to.save=jags.pars, 
                  model.file="Nmix-model-grouse.jag",
                  n.chains=3, n.adapt=100, n.burnin=0,
                  n.iter=2000, parallel=TRUE)


## ----bugs-mcmc2,size='footnotesize',message=FALSE,cache=TRUE,results='hide'----
jags.data2 <- jags.data; jags.data2$modswitch <- c(1,0,1,1) 
jm2 <- jags.basic(data=jags.data2, inits=jags.inits,
                  parameters.to.save=jags.pars, 
                  model.file="Nmix-model-grouse.jag",
                  n.chains=3, n.adapt=100, n.burnin=0,
                  n.iter=2000, parallel=TRUE)


## ----bugs-mcmc3,size='tiny',message=FALSE,cache=TRUE,results='hide'-----------
jags.data3 <- jags.data; jags.data3$modswitch <- c(1,0,0,1)
jm3 <- jags.basic(data=jags.data3, inits=jags.inits, parameters.to.save=jags.pars,
                  model.file="Nmix-model-grouse.jag", n.chains=3, n.adapt=100, n.burnin=0,
                  n.iter=2000, parallel=TRUE)


## ----bugs-mcmc4,size='tiny',message=FALSE,cache=TRUE,results='hide'-----------
jags.data4 <- jags.data; jags.data4$modswitch <- c(1,0,1,0)
jm4 <- jags.basic(data=jags.data4, inits=jags.inits, parameters.to.save=jags.pars,
                  model.file="Nmix-model-grouse.jag", n.chains=3, n.adapt=100, n.burnin=0,
                  n.iter=2000, parallel=TRUE)


## ----bugs-mcmc5,size='tiny',message=FALSE,cache=TRUE,results='hide'-----------
jags.data5 <- jags.data; jags.data5$modswitch <- c(1,0,0,0)
jm5 <- jags.basic(data=jags.data5, inits=jags.inits, parameters.to.save=jags.pars,
                  model.file="Nmix-model-grouse.jag", n.chains=3, n.adapt=100, n.burnin=0,
                  n.iter=2000, parallel=TRUE)


## ----bugs-mcmc6,size='tiny',message=FALSE,cache=TRUE,results='hide'-----------
jags.data6 <- jags.data; jags.data6$modswitch <- c(0,0,0,0)
jm6 <- jags.basic(data=jags.data6, inits=jags.inits, parameters.to.save=jags.pars,
                  model.file="Nmix-model-grouse.jag", n.chains=3, n.adapt=100, n.burnin=0,
                  n.iter=2000, parallel=TRUE)


## ----waic-fn-ld-op,size='small',include=FALSE---------------------------------
waic <- function(x, focus=c("y", "yN")) {
    vars <- coda::varnames(x)
    if(focus[1]=="y") {
        ld.samples <- as.matrix(x[,grep("ld.y.dot", vars)])
    } else if(focus[1]=="yN") {
        ld.samples <- as.matrix(x[,grep("ld.ydot.N", vars)])
    } else stop("focus should be either 'y' or 'yN'")
    lppd <- sum(log(colMeans(exp(ld.samples))))
    penalty <- sum(apply(ld.samples, 2, var))
    return(-2*(lppd-penalty))
}

## ----waic-fn,size='footnotesize'----------------------------------------------
waic <- function(x) {
    ## Parameter names
    vars <- coda::varnames(x)
    ## Extract log-density of y at each site for each post sample
    ld.samples <- as.matrix(x[,grep("ld.y.dot", vars)])
    ## Compute log-pointwise-predictive-density
    lppd <- sum(log(colMeans(exp(ld.samples))))
    ## Compute penalty
    penalty <- sum(apply(ld.samples, 2, var))
    ## Return WAIC
    return(-2*(lppd-penalty))
}


## ----waic1,size='scriptsize'--------------------------------------------------
(waic1 <- waic(jm1))

## ----waic2,size='scriptsize'--------------------------------------------------
(waic2 <- waic(jm2))

## ----waic3,size='scriptsize'--------------------------------------------------
(waic3 <- waic(jm3))


## ----waic4,size='scriptsize'--------------------------------------------------
(waic4 <- waic(jm4))

## ----waic5,size='scriptsize'--------------------------------------------------
(waic5 <- waic(jm5))

## ----waic6,size='scriptsize'--------------------------------------------------
(waic6 <- waic(jm6))


## ----waic-table,size='scriptsize'---------------------------------------------
waic.table <- data.frame(WAIC=c(waic1,waic2,waic3,waic4,waic5,waic5))
rownames(waic.table) <- paste("Model", 1:6)
waic.table$delta <- waic.table$WAIC-min(waic.table$WAIC)
waic.table <- waic.table[order(waic.table$WAIC),]
knitr::kable(waic.table, format="latex", digits=2, booktabs=TRUE)


## ----zip-sim,size='scriptsize'------------------------------------------------
nSites <- 20
nVisits <- 3
psi <- 0.5   # Proportion of extra-Poisson zeros
z <- rbinom(n=nSites, size=1, prob=psi) # Extra zeros
lam <- 5     # expected count when z=1
N <- rpois(n=nSites, lambda=lam*z)      # abundance at each site
p <- 0.5     # detection prob
y <- matrix(NA, nSites, nVisits)
for(i in 1:nSites) {
    y[i,] <- rbinom(n=nVisits, size=N[i], prob=p) # count data
}


## ----resid2,fig.width=9,out.width='60%',fig.align='center'--------------------
plot(fm2)


## ----parboot,size='scriptsize',warning=FALSE,cache=TRUE,fig.width=9,out.width='70%',fig.align='center'----
fitstat <- function(fm) { ## Sum of squared-residuals
    return(c(SSE=sum(residuals(fm)^2, na.rm=TRUE))) 
}
pb <- parboot(fm2, statistic=fitstat, nsim=200, ncores=3); plot(pb)


## ----parboot2,size='scriptsize'-----------------------------------------------
pb


## ----negbin,size='tiny',warning=FALSE,cache=TRUE------------------------------
(fm.nb <- pcount(~temp ~ elevation.s+utmN.s, data=grouse.umf, K=50, mixture="NB"))


## ----zip,size='tiny',warning=FALSE,cache=TRUE---------------------------------
(fm.zip <- pcount(~temp ~ elevation.s+utmN.s, data=grouse.umf, K=50, mixture="ZIP"))


## ----bugs2,size='tiny',echo=FALSE,comment='',background='lightblue'-----------
writeLines(readLines("Nmix-model-grouse2.jag"))


## ----bugs-mcmc-gof,size='footnotesize',message=FALSE,cache=TRUE,results='hide'----
jm <- jags.basic(data=jags.data, inits=jags.inits,
                 parameters.to.save=c(jags.pars, "SSE", "SSE.new"), 
                 model.file="Nmix-model-grouse2.jag",
                 n.chains=3, n.adapt=100, n.burnin=0,
                 n.iter=2000, parallel=TRUE)


## ----bugs-mcmc-gof-plot,size='footnotesize',fig.width=8,fig.align='center',out.width='55%'----
plot(as.matrix(jm[,c("SSE", "SSE.new")]))
abline(a=0, b=1, col="red")

