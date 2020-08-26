## ----buildit,include=FALSE,eval=FALSE-----------------------------------------
## rnw2pdf("lecture-popdy-basics")
## rnw2pdf("lecture-popdy-basics", tangle=TRUE)




## ----geo1,size='scriptsize',out.width='0.9\\textwidth',fig.width=8,fig.height=5,fig.align="center"----
time <- 0:100         
T <- length(time)     ## number of time points
r <- 0.01             ## growth rate
N1 <- 10*(1+r)^time   ## abundance at time t. N(0)=10
plot(time, N1, xlab="Time", ylab="Abundance", type="l", cex.lab=1.5)


## ----geo2,size='tiny',out.width='0.8\\textwidth',fig.width=8,fig.height=5,fig.align="center"----
time <- 0:100         
T <- length(time)     ## number of time points
r <- 0.01             ## growth rate
N2 <- c(10, rep(NA, T-1))  ## initial abundance = 10
for(t in 2:T) {
    N2[t] <- N2[t-1] + N2[t-1]*r
}
plot(time, N1, xlab="Time", ylab="Abundance", type="l", cex.lab=1.5)
lines(time, N2, col=rgb(0,0,1,0.5), lwd=10)
legend(0, 27, c("Geometric growth - option 1", "Geometric growth - option 2"),
       lwd=c(1,10), col=c(1, rgb(0,0,1,0.5)))


## ----geo3code,size='small'----------------------------------------------------
time <- 0:100         
T <- length(time)     ## number of time points
N3 <- c(10, rep(NA, T-1))  ## initial abundance = 10
r <- rep(NA, T-1)
rbar <- 0.01
sigma <- 0.2
for(t in 2:T) {
    r[t-1] <- rnorm(n=1, mean=rbar, sd=sigma)
    N3[t] <- N3[t-1] + N3[t-1]*r[t-1]
}


## ----geo3,size='scriptsize',out.width='0.99\\textwidth',fig.width=8,fig.height=5,fig.align="center"----
plot(time, N3, xlab="Time", ylab="Abundance", type="l", cex.lab=1.5)


## ----Nl0,echo=FALSE-----------------------------------------------------------
Time <- 0:100
T <- length(Time)
rmax <- .1
K <- 50
Nl <- rep(0, T)
Nl[1] <- 2
for(t in 2:T) {
    Nl[t] <- Nl[t-1] + Nl[t-1]*rmax*(1 - Nl[t-1]/K)
}
Nl4 <- Nl3 <- Nl2 <- Nl
Nl4[1] <- Nl3[1] <- Nl2[1] <- 2
for(t in 2:T) {
    Nl2[t] <- Nl2[t-1] + Nl2[t-1]*0.5*(1-Nl2[t-1]/K)
    Nl3[t] <- Nl3[t-1] + Nl3[t-1]*2.0*(1-Nl3[t-1]/K)
    Nl4[t] <- Nl4[t-1] + Nl4[t-1]*3.0*(1-Nl4[t-1]/K)
}

## ----Nl,include=FALSE,echo=FALSE,fig.width=8,fig.height=6---------------------
plot(Time, Nl, lwd=4, type="l", ylim=c(0, 100), cex.lab=1.3,
     xlab="Time (t)", ylab="Population size (N)", col="purple")
legend(0, 100, c("r=0.1", "", "", ""),
       lwd=4, col=c("purple", NA, NA, NA))

## ----Nl2,include=FALSE,echo=FALSE,fig.width=8,fig.height=6--------------------
plot(Time, Nl, lwd=4, type="l", ylim=c(0, 100), cex.lab=1.3,
     xlab="Time (t)", ylab="Population size (N)", col="purple")
lines(Time, Nl2, lwd=4, col="blue")
legend(0, 100, c("r=0.1", "r=0.5", "", ""),
       lwd=4, col=c("purple", "blue", NA, NA))

## ----Nl3,include=FALSE,echo=FALSE,fig.width=8,fig.height=6--------------------
plot(Time, Nl, lwd=4, type="l", ylim=c(0, 100), cex.lab=1.3,
     xlab="Time (t)", ylab="Population size (N)", col="purple",
     main="Damped oscillation")
lines(Time, Nl2, lwd=4, col="blue")
lines(Time, Nl3, lwd=4, col="orange")
legend(0, 100, c("r=0.1", "r=0.5", "r=2.0", ""),
       lwd=4, col=c("purple", "blue", "orange", NA))

## ----Nl4,include=FALSE,echo=FALSE,fig.width=8,fig.height=6--------------------
plot(Time, Nl, lwd=4, type="l", ylim=c(0, 100), cex.lab=1.3,
     xlab="Time (t)", ylab="Population size (N)", col="purple",
     main="Chaos")
lines(Time, Nl2, lwd=4, col="blue")
lines(Time, Nl3, lwd=4, col="orange")
lines(Time, Nl4, lwd=4, col="gray")
legend(0, 100, c("r=0.1", "r=0.5", "r=2.0", "r=3.0"),
       lwd=4, col=c("purple", "blue", "orange", "gray"))


## ----proj1code,size='footnotesize'--------------------------------------------
T <- 30                         ## time steps
n <- matrix(NA, nrow=3, ncol=T) ## age-class abundance matrix    
n[,1] <- c(50, 40, 10)          ## abundance at t=1
s <- c(0.4, 0.5, 0.3)           ## survival rates
f <- c(0, 0.8, 1.7)             ## fecundities

A <- matrix(c(f, s[1], 0, 0, 0, s[2], s[3]), 
            nrow=3, ncol=3, byrow=TRUE)  
A                               ## Leslie matrix

for(t in 2:T) {
    n[,t] <- A %*% n[,t-1]      ## matrix multiplication
}

lambda <- n[,-1] / n[,-T]       ## growth rates


## ----proj1,size='scriptsize',out.width='0.99\\textwidth',fig.width=8,fig.height=5,fig.align="center"----
matplot(1:T, t(n), type="o", pch=16, xlab="Time", ylab="Population size",
        cex.lab=1.3, ylim=c(0, 60), col=c("black", "orange", "purple"))
legend(20, 60, c("Age class 1", "Age class 2", "Age class 3"),
       col=c("black", "orange", "purple"), pch=16, lty=1:3)


## ----lambda1,include=FALSE,echo=FALSE,fig.width=8,fig.height=6----------------
matplot(t(lambda), type="o", pch=16,
        xlab="Time", ylab="Population growth rate (lambda)",
        cex.lab=1.3,
        col=c("black", "orange", "purple"))
legend(20, 2.4, c("Age class 1", "Age class 2", "Age class 3"),
       col=c("black", "orange", "purple"), pch=16, lty=1:3)


## ----prop1,include=FALSE,echo=FALSE,fig.width=8,fig.height=6------------------
N <- colSums(n)
C <- sweep(n, 2, N, "/")
matplot(1:T, t(C), type="o", pch=16, xlab="Time",
        cex.lab=1.3,
        ylab="Proportion in age class",
        col=c("black", "orange", "purple"), ylim=c(0, 1))
legend(20, 1, c("Age class 1", "Age class 2", "Age class 3"),
       col=c("black", "orange", "purple"), pch=16, lty=1:3)


## ----proj1recode,size='footnotesize',eval=FALSE-------------------------------
## for(t in 2:T) {
##     n[,t] <- A %*% n[,t-1]
## }


## ----eig,size='scriptsize'----------------------------------------------------
vw <- eigen(A)
Re(vw$values[1])                          ## lambda
Re(vw$vectors[,1]/sum(vw$vectors[,1]))    ## stable age distribution
vw2 <- eigen(t(A))
Re(vw2$vectors[,1]/sum(vw2$vectors[,1]))  ## reproductive values

