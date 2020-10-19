## ----buildit,include=FALSE,eval=FALSE-----------------------------------------
## rnw2pdf("lecture-cap-recap-nonspatial")
## rnw2pdf("lecture-cap-recap-nonspatial", tangle=TRUE)




## ----include=FALSE,echo=FALSE-------------------------------------------------
set.seed(34889243)


## ----sim-M0-pars,size='scriptsize'--------------------------------------------
N <- 100
p <- 0.2


## ----sim-M0-ch,size='scriptsize'----------------------------------------------
J <- 4  ## Occasions
y.all <- matrix(NA, N, J)
for(i in 1:N) {
    y.all[i,] <- rbinom(J, 1, p)
}


## ----sim-M0-y1,size='scriptsize'----------------------------------------------
captured <- rowSums(y.all)>0
(n <- sum(captured))
y <- y.all[captured,]
y[1:3,]


## ----M0-hist,size='scriptsize'------------------------------------------------
histories <- apply(y, 1, paste, collapse="")
sort(table(histories))


## ----M0-freq,size='scriptsize'------------------------------------------------
y.tilde <- rowSums(y)
sort(table(y.tilde))


## ----sim-Mt-pars,size='scriptsize'--------------------------------------------
p.t <- c(0.3, 0.5, 0.2, 0.4)


## ----sim-Mt-ch,size='scriptsize'----------------------------------------------
y.all.Mt <- matrix(NA, N, J)
for(i in 1:N) {
    y.all.Mt[i,] <- rbinom(J, 1, p.t) }


## ----sim-Mt-y1,size='scriptsize'----------------------------------------------
captured.Mt <- rowSums(y.all.Mt)>0
n.Mt <- sum(captured.Mt)
y.Mt <- y.all.Mt[captured.Mt,]
y.Mt[1:3,]
colSums(y.Mt)


## ----sim-Mb-pars,size='scriptsize'--------------------------------------------
p.b <- 0.3
c <- 0.5  ## Trap happy


## ----sim-Mb-ch,size='scriptsize'----------------------------------------------
y.all.Mb <- matrix(NA, N, J)
prevcap <- matrix(FALSE, N, J)
for(i in 1:N) {
    y.all.Mb[i,1] <- rbinom(1, 1, p.b)
    for(j in 2:J) {
        prevcap[i,j] <- any(y.all.Mb[i,1:(j-1)]>0)
        prob <- ifelse(prevcap[i,j], c, p.b)
        y.all.Mb[i,j] <- rbinom(1, 1, prob)
    }
}


## ----sim-Mb-y1,size='scriptsize'----------------------------------------------
captured.Mb <- rowSums(y.all.Mb)>0
n.Mb <- sum(captured.Mb)
y.Mb <- y.all.Mb[captured.Mb,]
## prevcap[1:3,]


## ----sim-Mh-pars,size='scriptsize'--------------------------------------------
logit.p.bar <- -1  ## Mean p on logit scale
logit.p.var <- 1   ## SD of p on logit scale
p.h <- plogis(rnorm(N, logit.p.bar, logit.p.var))


## ----sim-Mh-ch,size='scriptsize'----------------------------------------------
y.all.Mh <- matrix(NA, N, J)
for(i in 1:N) {
    y.all.Mh[i,] <- rbinom(J, 1, p.h[i])
}


## ----sim-Mh-y1,size='scriptsize'----------------------------------------------
captured.Mh <- rowSums(y.all.Mh)>0
n.Mh <- sum(captured.Mh)
y.Mh <- y.all.Mh[captured.Mh,]
#y.Mh[1:3,]


## ----nll-M0,echo=TRUE,size='scriptsize'---------------------------------------
nll.M0 <- function(pars, y) {           ## Negative log-likelihood
    n <- nrow(y);       J <- ncol(y)
    N <- exp(pars[1])
    n0 <- N-n
    if(n0<0) return(NA)
    p <- plogis(pars[2])
    ld.y1 <- sum(dbinom(y, 1, p, log=TRUE))
    p0 <- (1-p)^J
    ld.n0 <- lgamma(N+1)-lgamma(n0+1)+n0*log(p0)
    nll <- -(ld.y1+ld.n0)
    return(nll)
}


## ----opt-nll-M0, size='scriptsize'--------------------------------------------
fm.M0 <- optim(c(log.N=4,logit.p=0), nll.M0, y=y, hessian=TRUE)
fm.M0.est <- data.frame(Estimate=c(fm.M0$par[1], fm.M0$par[2]),
                        SE=sqrt(diag(solve(fm.M0$hessian))))
fm.M0.est


## ----opt-nll-M0-back, size='scriptsize'---------------------------------------
c(N.hat=exp(fm.M0$par[1]), p.hat=plogis(fm.M0$par[2]))


## ----dg,size='scriptsize'-----------------------------------------------------
c(N=N, p=p)


## ----eval=FALSE,include=FALSE,echo=FALSE--------------------------------------
## nll.M0.2 <- function(pars, y) {
##     n <- nrow(y);       J <- ncol(y)
##     N <- exp(pars[1])
##     n0 <- N-n
##     if(n0<0) return(NA)
##     p <- plogis(pars[2])
##     J <- ncol(y)
##     p.star <- 1-(1-p)^J
##     ydot <- rowSums(y)
##     lcd.y <- sum(dbinom(ydot, J, p, log=TRUE)-log(p.star))
##     ld.n <- lchoose(N, n)+n*log(p.star)+n0*log(1-p.star)
##     nll <- -(lcd.y+ld.n)
##     return(nll)
## }
## 
## fm.M0.2 <- optim(c(5,0), nll.M0.2, y=y, hessian=TRUE)
## c(N=exp(fm.M0.2$par[1]), p=plogis(fm.M0.2$par[2]))


## ----eval=FALSE,include=FALSE,echo=FALSE--------------------------------------
## nll.M0.cn <- function(pars, y) {
##     p <- plogis(pars[1])
##     J <- ncol(y)
##     p.star <- 1-(1-p)^J
##     ydot <- rowSums(y)
##     nll <- -sum(dbinom(ydot, J, p, log=TRUE)-log(p.star))
##     return(nll)
## }
## 
## fm.M0.cn <- optim(0, nll.M0.cn, y=y, hessian=TRUE, method="Brent", lower=-10, upper=10)
## c(p=plogis(fm.M0.cn$par[1]))


## ----Rcapture,size='scriptsize'-----------------------------------------------
## install.packages("Rcapture")
library(Rcapture)
closedp(y)


## ----mra,size='scriptsize'----------------------------------------------------
## install.packages("mra")
library(mra)
F.huggins.estim(~1, histories=y)


## ----RMark,size='scriptsize',warning=FALSE,results='hide',cache=TRUE----------
## install.packages("RMark") ## Must install MARK TOO!!
library(RMark)
y.ch <- data.frame(ch=apply(y, 1, paste, collapse=""))
mark.M0 <- mark(data=y.ch, model="Closed", silent=TRUE,
                model.parameters=list(p=list(formula=~1,share=TRUE)))


## ----RMark-pn0,size='scriptsize'----------------------------------------------
mark.M0$results$real      


## ----RMark-N,size='scriptsize'------------------------------------------------
mark.M0$results$derived   


## ----RMark-Mt,size='scriptsize',warning=FALSE,results='hide',cache=TRUE-------
yt.ch <- data.frame(ch=apply(y.Mt, 1, paste, collapse=""))
mark.Mt <- mark(data=yt.ch, model="Closed", silent=TRUE,
                model.parameters=list(p=list(formula=~time,share=TRUE)))


## ----RMark-pn0-Mt,size='scriptsize'-------------------------------------------
mark.Mt$results$real      


## ----RMark-N-Mt,size='scriptsize'---------------------------------------------
mark.Mt$results$derived   


## ----RMark-Mb,size='scriptsize',warning=FALSE,results='hide',cache=TRUE-------
yb.ch <- data.frame(ch=apply(y.Mb, 1, paste, collapse=""))
mark.Mb <- mark(data=yb.ch, model="Closed", silent=TRUE) ## Default model


## ----RMark-pn0-Mb,size='scriptsize'-------------------------------------------
mark.Mb$results$real      


## ----RMark-N-Mb,size='scriptsize'---------------------------------------------
mark.Mb$results$derived   


## ----aug-y,size='scriptsize'--------------------------------------------------
M <- 200
y.aug <- matrix(0, M, J)
y.aug[1:nrow(y.Mb),] <- y.Mb


## ----un-prevcap,size='scriptsize'---------------------------------------------
prevcap.aug <- matrix(FALSE, M, J)
prevcap.aug[1:nrow(y.Mb),] <- prevcap[captured.Mb,]
prevcap.aug <- ifelse(prevcap.aug, "Yes", "No")


## ----un-occ,size='scriptsize'-------------------------------------------------
occasion.aug <- matrix(as.character(1:J), M, J, byrow=TRUE)


## ----un-umf,results='hide',size='scriptsize',warning=FALSE--------------------
library(unmarked)
umf <- unmarkedFrameOccu(y=y.aug,
                         obsCovs=list(prevcap=prevcap.aug,
                                      occasion=occasion.aug))


## ----un-M0,size='scriptsize'--------------------------------------------------
fm.M0 <- occu(~1~1, umf)
fm.M0


## ----un-Mt,size='scriptsize'--------------------------------------------------
fm.Mt <- occu(~occasion~1, umf)
fm.Mt


## ----un-Mb,size='scriptsize'--------------------------------------------------
fm.Mb <- occu(~prevcap~1, umf)
fm.Mb


## ----un-re,size='footnotesize'------------------------------------------------
re.M0 <- ranef(fm.M0)
re.Mt <- ranef(fm.Mt)
re.Mb <- ranef(fm.Mb)


## ----un-N,size='footnotesize'-------------------------------------------------
N.post.M0 <- predict(re.M0, func=sum, nsim=1000)
N.post.Mt <- predict(re.Mt, func=sum, nsim=1000)
N.post.Mb <- predict(re.Mb, func=sum, nsim=1000)


## ----un-N-post, size='tiny',fig.height=3,out.width="100%",fig.show='hide', echo=-1----
par(mfrow=c(1,3), mai=c(0.6,0.6,0.5,0.1))
plot(table(N.post.M0), xlab="N", ylab="Posterior probability", main="M0",
     xlim=c(75, 130), ylim=c(0,150)); abline(v=N, col="red")
plot(table(N.post.Mt), xlab="N", ylab="Posterior probability", main="Mt",
     xlim=c(75, 130), ylim=c(0,150)); abline(v=N, col="red")
plot(table(N.post.Mb), xlab="N", ylab="Posterior probability", main="Mb",
     xlim=c(75, 130), ylim=c(0,150)); abline(v=N, col="red")


## ----bugs-M0-aug,size='small'-------------------------------------------------
writeLines(readLines("M0-aug.jag"))


## ----jd-M0-aug,size='scriptsize'----------------------------------------------
y.aug <- matrix(0, M, J)
y.aug[1:n,] <- y
jags.data.M0 <- list(y=y.aug, M=M, J=J)


## ----ji-M0-aug,size='scriptsize'----------------------------------------------
ji.M0 <- function() list(z=rep(1,M), psi=runif(1), p=runif(1))
jp.M0 <- c("p", "psi", "N")


## ----mcmc-M0-aug,size='scriptsize',results='hide'-----------------------------
library(jagsUI)
jags.post.M0 <- jags.basic(data=jags.data.M0, inits=ji.M0,
                           parameters.to.save=jp.M0,
                           model.file="M0-aug.jag",
                           n.chains=3, n.adapt=100, n.burnin=0,
                           n.iter=2000, parallel=TRUE)


## ----plot-mcmc-M0,size='footnotesize',out.width="0.7\\textwidth",fig.align='center',cache=TRUE----
plot(jags.post.M0[,jp.M0])


## ----bugs-M0,size='footnotesize'----------------------------------------------
writeLines(readLines("M0.jag"))


## ----jd-M0-noaug,size='scriptsize'--------------------------------------------
n0max <- 1000  ## Upper limit of prior on n0
jags.data.M0.noaug <- list(y=y.aug, n=n, J=J,
                           ## Prior probs for n0
                           n0probs=rep(1/n0max, n0max), zero=0)


## ----ji-M0-noaug,size='scriptsize'--------------------------------------------
ji.M0.noaug <- function() list(n0=rpois(1, 5), p=runif(1, 0, 0.1))
jp.M0.noaug <- c("p", "N")


## ----mcmc-M0-noaug,size='scriptsize',results='hide'---------------------------
jags.post.M0.noaug <- jags.basic(data=jags.data.M0.noaug,
                                 inits=ji.M0.noaug,
                                 parameters.to.save=jp.M0.noaug,
                                 model.file="M0.jag",
                                 n.chains=3, n.adapt=100, n.burnin=0,
                                 n.iter=2000, parallel=TRUE)


## ----plot-mcmc-M0-noaug,size='footnotesize',out.width="0.7\\textwidth",fig.align='center',cache=TRUE----
plot(jags.post.M0.noaug[,jp.M0.noaug])


## ----bugs-Mt-aug,size='small'-------------------------------------------------
writeLines(readLines("Mt-aug.jag"))


## ----ji-Mt-aug,size='scriptsize'----------------------------------------------
jags.data.Mt <- jags.data.M0
ji.Mt <- function() list(z=rep(1,M), psi=runif(1), p=runif(4))
jp.Mt <- c("p", "psi", "N")


## ----mcmc-Mt-aug,size='scriptsize',results='hide'-----------------------------
jags.post.Mt <- jags.basic(data=jags.data.Mt, inits=ji.Mt,
                           parameters.to.save=jp.Mt,
                           model.file="Mt-aug.jag",
                           n.chains=3, n.adapt=100, n.burnin=0,
                           n.iter=2000, parallel=TRUE)


## ----plot-mcmc-Mt1,size='footnotesize',out.width="0.7\\textwidth",fig.align='center',cache=TRUE----
plot(jags.post.Mt[,paste0("p[", 1:4, "]")])

