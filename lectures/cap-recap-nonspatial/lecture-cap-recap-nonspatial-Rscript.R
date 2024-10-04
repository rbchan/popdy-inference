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


## ----sim-Mhm-pars,size='scriptsize'-------------------------------------------
mixture.prob <- 0.6
group <- rbinom(N, 1, mixture.prob)    ## Two groups
p.Mh.mix <- ifelse(group==0, 0.2, 0.7) ## Two-point mixture


## ----sim-Mhm-ch,size='scriptsize'---------------------------------------------
y.all.Mh.mix <- matrix(NA, N, J)
for(i in 1:N) {
    y.all.Mh.mix[i,] <- rbinom(J, 1, p.Mh.mix[i])
}


## ----sim-Mhm-y1,size='scriptsize'---------------------------------------------
captured.Mh.mix <- rowSums(y.all.Mh.mix)>0
n.Mh.mix <- sum(captured.Mh.mix)
y.Mh.mix <- y.all.Mh.mix[captured.Mh.mix,]


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


## ----RMark,size='scriptsize',warning=FALSE,results='hide',cache=FALSE,warning=FALSE,message=FALSE----
## install.packages("RMark") ## Must install MARK too!!
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


## ----RMark-Mh,size='scriptsize',warning=FALSE,results='hide',cache=FALSE------
yh.ch <- data.frame(ch=apply(y.Mh.mix, 1, paste, collapse=""))
mark.Mh <- mark(data=yh.ch, silent=TRUE, model="HetClosed", # "HugHet", 
                model.parameters=list(p=list(formula=~mixture,share=TRUE)))


## ----RMark-pn0-Mh,size='scriptsize'-------------------------------------------
mark.Mh$results$real      


## ----RMark-N-Mh,size='scriptsize'---------------------------------------------
mark.Mh$results$derived   


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


## ----bugs-M0-aug,size='footnotesize',comment='',echo=-1,background='beige'----
writeLines(readLines("M0-aug.jag"))


## ----jd-M0-aug,size='scriptsize'----------------------------------------------
y.aug <- matrix(0, M, J)
y.aug[1:n,] <- y
jags.data.M0 <- list(y=y.aug, M=M, J=J)


## ----ji-M0-aug,size='scriptsize'----------------------------------------------
ji.M0 <- function() list(z=rep(1,M), psi=runif(1), p=runif(1))
jp.M0 <- c("p", "psi", "N")


## ----mcmc-M0-aug,size='scriptsize',results='hide',warning=FALSE,message=FALSE----
library(jagsUI)
jags.post.M0 <- jags.basic(data=jags.data.M0, inits=ji.M0,
                           parameters.to.save=jp.M0,
                           model.file="M0-aug.jag",
                           n.chains=3, n.adapt=100, n.burnin=0,
                           n.iter=2000, parallel=TRUE)


## ----plot-mcmc-M0,size='footnotesize',out.width="0.7\\textwidth",fig.align='center',cache=TRUE----
plot(jags.post.M0[,jp.M0])


## ----bugs-Mb-aug,size='footnotesize',comment='',echo=-1,background='beige'----
writeLines(readLines("Mb-aug.jag"))


## ----ji-Mb-aug,size='scriptsize'----------------------------------------------
jags.data.Mb <- jags.data.M0
jags.data.Mb$y[1:nrow(y.Mb),] <- y.Mb
jags.data.Mb$prevcap <- ifelse(prevcap.aug=="Yes", 1, 0)
ji.Mb <- function() list(z=rep(1,M), psi=runif(1), p=runif(1), c=runif(1))
jp.Mb <- c("p", "c", "psi", "N")


## ----mcmc-Mb-aug,size='scriptsize',results='hide'-----------------------------
jags.post.Mb <- jags.basic(data=jags.data.Mb, inits=ji.Mb,
                           parameters.to.save=jp.Mb,
                           model.file="Mb-aug.jag",
                           n.chains=3, n.adapt=100, n.burnin=0,
                           n.iter=2000, parallel=TRUE)


## ----plot-mcmc-Mb,size='footnotesize',out.width="0.7\\textwidth",fig.align='center',cache=TRUE----
plot(jags.post.Mb[,jp.Mb])


## ----bugs-Mhm-aug,size='footnotesize',comment='',echo=-1,background='beige'----
writeLines(readLines("Mhm-aug.jag"))


## ----bugs-Mh-aug,size='footnotesize',comment='',echo=-1,background='beige'----
writeLines(readLines("Mh-aug.jag"))

