list.files()


## Import data

sites.in <- read.csv("Site_Data_WILD8390.csv")
dets.in <- read.csv("Detection_Data_WILD8390.csv")
surveys.in <- read.csv("Survey_Covs_WILD8390.csv")

str(sites.in)
str(dets.in)
str(surveys.in)


## Subset data by sites surveyed during *standard* survey

site.names <- c(
    grep("PCCow", sites.in$PointName, value=TRUE),
    grep("PCDry", sites.in$PointName, value=TRUE),
    grep("PCWS0", sites.in$PointName, value=TRUE))
site.names

## Sites
sites <- subset(sites.in, PointName %in% site.names,
                select=-1)
rownames(sites) <- sites$PointName
sites$PointName <- NULL
str(sites)

## Surveys
names(surveys.in)
surveys <- subset(surveys.in, PointName %in% site.names,
                  select=!(names(surveys.in) %in%
                           c("X", "min_time",
                             "stand_time", "Year")))
str(surveys)


## Detections
names(dets.in)
dets <- subset(dets.in, Site %in% site.names, select=-1)
dets$PointName <- factor(dets$Site)
dets$species <- factor(dets$species)
dets$Site <- NULL
dets$datetime <- as.POSIXct(paste(dets$Date, dets$Time),
                            format="%Y-%m-%d %H:%M:%OS")
dets$year <- factor(format(dets$datetime, format="%Y"))
dets$Date <- dets$Time <- NULL
dets$interval <- factor(as.character(dets$TimeInterval))
dets$TimeInterval <- NULL


str(dets)





## Put detection data in a site x survey x occasion x year species array

## TODO: Make sure a site wasn't surveyed more than once during each year

det.arr <- with(dets, table(PointName, interval, year, species))
str(det.arr)


## TODO: Write a function to do this for any species
dets.cawa.14 <- det.arr[,,"2014","CAWA"]
dets.cawa.15 <- det.arr[,,"2015","CAWA"]
dets.cawa.16 <- det.arr[,,"2016","CAWA"]
dets.cawa.17 <- det.arr[,,"2017","CAWA"]
dets.cawa.18 <- det.arr[,,"2018","CAWA"]
dets.cawa.19 <- det.arr[,,"2019","CAWA"]
dets.cawa.20 <- det.arr[,,"2020","CAWA"]

cawa.14 <- data.frame(y=ifelse(rowSums(dets.cawa.14)>0, 1L, 0L),
                      Site=rownames(dets.cawa.14),
                      sites[rownames(dets.cawa.14),],
                      Year="2014")
str(cawa.14)
cawa.15 <- data.frame(y=ifelse(rowSums(dets.cawa.15)>0, 1L, 0L),
                      Site=rownames(dets.cawa.15),
                      sites[rownames(dets.cawa.15),],
                      Year="2015")
cawa.16 <- data.frame(y=ifelse(rowSums(dets.cawa.16)>0, 1L, 0L),
                      Site=rownames(dets.cawa.16),
                      sites[rownames(dets.cawa.16),],
                      Year="2016")
cawa.17 <- data.frame(y=ifelse(rowSums(dets.cawa.17)>0, 1L, 0L),
                      Site=rownames(dets.cawa.17),
                      sites[rownames(dets.cawa.17),],
                      Year="2017")
cawa.18 <- data.frame(y=ifelse(rowSums(dets.cawa.18)>0, 1L, 0L),
                      Site=rownames(dets.cawa.18),
                      sites[rownames(dets.cawa.18),],
                      Year="2018")
cawa.19 <- data.frame(y=ifelse(rowSums(dets.cawa.19)>0, 1L, 0L),
                      Site=rownames(dets.cawa.19),
                      sites[rownames(dets.cawa.19),],
                      Year="2019")
cawa.20 <- data.frame(y=ifelse(rowSums(dets.cawa.20)>0, 1L, 0L),
                      Site=rownames(dets.cawa.20),
                      sites[rownames(dets.cawa.20),],
                      Year="2020")
str(cawa.20)

cawa <- rbind(cawa.14, cawa.15, cawa.16, cawa.17, cawa.18, cawa.19, cawa.20)


str(cawa)




fm1 <- glm(y ~ Elevation+I(Elevation^2), binomial, cawa)
fm2 <- glm(y ~ Year, binomial, cawa)
fm3 <- glm(y ~ Elevation+I(Elevation^2)+Year, binomial, cawa)
fm4 <- glm(y ~ (Elevation+I(Elevation^2))*Year, binomial, cawa)

AIC(fm1, fm2, fm3, fm4)

summary(fm3)

elev.min <- min(cawa$Elevation)
elev.max <- max(cawa$Elevation)

yrs <- 2014:2020

pred.data <- data.frame(Year=factor(rep(as.character(yrs), each=50)),
                        Elevation=rep(seq(elev.min, elev.max, length=50),
                                      times=length(yrs)))

pred3 <- predict(fm3, newdata=pred.data, type="response", se.fit=TRUE)

pred3dat <- data.frame(pred=pred3$fit, pred.data)

plot(pred ~ Elevation, pred3dat, subset=Year=="2014",
     type="l", ylim=c(0, 1))
lines(pred ~ Elevation, pred3dat, subset=Year=="2015", col=2)
lines(pred ~ Elevation, pred3dat, subset=Year=="2016", col=3)
lines(pred ~ Elevation, pred3dat, subset=Year=="2017", col=4)
lines(pred ~ Elevation, pred3dat, subset=Year=="2018", col=5)
lines(pred ~ Elevation, pred3dat, subset=Year=="2019", col=6)
lines(pred ~ Elevation, pred3dat, subset=Year=="2020", col=7)
