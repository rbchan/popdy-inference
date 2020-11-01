## Import data
caps.in <- read.csv("captures.csv")
nets.in <- read.csv("Constant_Effort_Nets.csv")



## Format date
caps <- caps.in
caps$NetDate <- as.Date(caps.in$NetDate,
                        format="%d-%b-%y")
caps$datetime <- as.POSIXct(
    paste(as.character(caps$NetDate), caps$CapTime),
    format="%Y-%m-%d %H:%M:%S")
caps$year <- format(caps$NetDate, "%Y")
caps$month <- format(caps$NetDate, "%m")
caps$day <- format(caps$NetDate, "%d")

head(caps.in)
head(caps)




## Format net data
## Only keep nets with 2 occasions
open.2020 <- nets.in[,c("Net", "site",
                        grep("open2020", colnames(nets.in), value=TRUE))]
nets.2020 <- nets.in[rowSums(open.2020[,-(1:2)])==2 & nets.in$UTM_N>1e6,
                     c("Net", "site", "UTM_W", "UTM_N")]
colnames(nets.2020) <- c("net", "site", "utmW", "utmN")

nrow(nets.2020)

plot(utmN~utmW, nets.2020, asp=1, pch=3)



## Format BTBW data
btbw <- subset(caps, Species=="BTBW")

table(btbw$year)

btbw.2020 <- subset(btbw, year=="2020" & Net%in%nets.2020$net)

## Add occasion
btbw.2020$occasion <- NA
site.btbw.2020 <- unique(btbw.2020$Site)


## for(i in 1:length(site.btbw.2020)) {
##     btbw.2020.site <- subset(btbw.2020, Site==site.btbw.2020[i])
## }

for(i in 1:nrow(btbw.2020)) {
    all.dates.i <- unique(sort(
        subset(btbw.2020, Site==btbw.2020[i,"Site"],
               select="NetDate", drop=TRUE)))
    date.i <- btbw.2020[i,"NetDate"]
    occasion.i <- which(all.dates.i %in% date.i)
    if(length(occasion.i)<1) stop("bad")
    btbw.2020$occasion[i] <- occasion.i
}

btbw.2020

btbw.2020.out <- data.frame(session=1,
                            btbw.2020[,c("BandNumber", "occasion", "Net")])
colnames(btbw.2020.out) <- c("session", "band", "occasion", "net")
btbw.2020.out$band <- as.character(btbw.2020.out$band)


nets.2020.out <- nets.2020[,c("net", "utmW", "utmN")]


library(secr)


all(btbw.2020.out$net %in% nets.2020$net)

nets.2020.secr <- nets.2020.out[,-1]
colnames(nets.2020.secr) <- c("x","y")
rownames(nets.2020.secr) <- nets.2020.out$net


nets.2020.traps <- read.traps(data=nets.2020.secr, detector="proximity")

ch <- make.capthist(captures=btbw.2020.out,
                    traps=nets.2020.traps,
                    fmt="trapID")

summary(ch)


secr.fit(ch, model=list(D=~1, g0=~1, sigma=~1), buffer=1500,
         start=c(-1, -4, 4), trace=FALSE)


secr.fit(ch, model=list(D=~1, g0=~bk, sigma=~1), buffer=10000,
         start=c(-1, 0, 0, 2), trace=FALSE)
