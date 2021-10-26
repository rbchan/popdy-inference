## Import bear data from 'ga-bear-pva' github repo


load("ga_bear_data_females.gzip")

ls()

str(ch4d.f)

apply(ch4d.f, 4, sum)

dimnames(ch4d.f)[[3]] <- paste0("week", 1:8)

dimnames(ch4d.f)


## Extract data from 2016

caps.2016 <- ch4d.f[,,,"2016"]
caps.2016 <- caps.2016[rowSums(caps.2016)>0,,]

str(caps.2016)





## Format for secr

caps.2016.secr <- data.frame(session=1,
                             individual=rep(slice.index(caps.2016, 1), caps.2016),
                             occasion=rep(slice.index(caps.2016, 3), caps.2016),
                             trap=rep(slice.index(caps.2016, 2), caps.2016))


## Export secr data
write.table(caps.2016.secr, file="secr_bear_encounter_data_file.csv",
            row.names=FALSE, col.names=FALSE, sep=",")
traps.secr <- data.frame(trap=1:nrow(traps), traps)
write.table(traps.secr, file="secr_bear_trap_data_file.csv",
            row.names=FALSE, col.names=FALSE, sep=",")

library(secr)


sch <- read.capthist(captfile="secr_bear_encounter_data_file.csv",
                     trapfile="secr_bear_trap_data_file.csv",
                     detector="proximity", fmt="trapID")
summary(sch)

fm.M0 <- secr.fit(sch, buffer=2000)
fm.Mt <- secr.fit(sch, model=list(D=~1, g0=~t, sigma=~1))

region.N(fm.M0)
region.N(fm.Mt)


## Export data for JAGS
save(caps.2016, traps, file="jags_bear_data.gzip")






## Heather's code

caps <- dget("GAbears_caps.txt")
traps <- dget("GABears_traps.txt")
oper <- read.csv("GABears_trapeffort.csv")[,-1]

y.tilde <- apply(caps, c(1,2), sum)

n <- nrow(caps)
library(raster)
statespace <- raster("state-space360.tif")
plot(statespace)
points(traps)
buff <- 4000
xlim <- c(min(traps[,1] -buff), max(traps[,1]+buff))
ylim <- c(min(traps[,2] -buff), max(traps[,2]+buff))
lines(c(xlim, rev(xlim), xlim[1]), c(ylim[1], ylim[1], ylim[2], ylim[2], ylim[1]))

#augment data

M <- 100
J <- nrow(traps)
K <- ncol(oper)

writeLines(
"model{
    psi ~ dunif(0, 1)
    g0 ~ dunif(0, 1)
    sigma ~ dunif(0, 5000)
  for(i in 1:M) {
      s[i,1] ~ dunif(xlim[1], xlim[2])
      s[i,2] ~ dunif(ylim[1], ylim[2])
      z[i] ~ dbern(psi)
      for(j in 1:J) {
        dist[i,j] <- sqrt((s[i,1]-x[j,1])^2 + (s[i,2]-x[j,2])^2)
        p[i,j] <- g0*exp(-dist[i,j]^2/(2*sigma^2))
      }
    }
    for(i in 1:n) { ## Model for observed capture histories
      for(j in 1:J) {
        y.tilde[i,j] ~ dbinom(p[i,j], K)
      }
    }
    for(i in (n+1):M) { ## Model for augmented guys
      PrAtLeastOneCap[i] <- 1-prod(1-p[i,])^K
      zero[i] ~ dbern(PrAtLeastOneCap[i]*z[i])
    }
    EN <- M*psi
    N <- sum(z)
}"
                      
           ,"SCR0_bears.txt")

jd <- list(y.tilde = y.tilde, n =n, M=M, J=J, z = c(rep(1,n), rep(NA, M-n)), K = K, zero = rep(0,M), x=traps, xlim = xlim, ylim = ylim)
s.init <- array(NA, dim = c(n, 2))
for(i in 1:nrow(caps)){
  t <- which(caps[i,,] > 0, arr.ind = T)[,1]
  s.init[i,1] <- mean(traps[t,1])
  s.init[i,2] <- mean(traps[t,2])
}



ji <- function(){
  list(z = c(rep(NA,n), rep(0, M-n)), psi = n/M, s = rbind(s.init, cbind(runif(M-n, xlim[1], xlim[2]), runif(M-n, ylim[1], ylim[2]))), g0 = runif(1), sigma = 4000)
}

library(jagsUI)
bears.j <- jags.basic(data = jd, inits = ji, parameters.to.save = c("g0", "sigma", "EN", "N", "psi"), model.file = "SCR0_bears.txt", n.chains = 3, n.adapt = 100, n.burnin = 0, n.iter = 2000, parallel = T )






