list.files()


## Import Coweeta data

sites.in <- read.csv("Site_Data_WILD8390.csv")
dets.in <- read.csv("Detection_Data_WILD8390.csv")
surveys.in <- read.csv("Survey_Covs_WILD8390.csv")

str(sites.in)
str(dets.in)
str(surveys.in)


## Import grouse data

grouse.dets.in <- read.csv("GrouseDetectionsQuery.csv")
grouse.surveys.in <- read.csv("RuffedGrouseDrummingSurveys Query.csv")

str(grouse.dets.in)
str(grouse.surveys.in)


grouse.surveys <- grouse.surveys.in
grouse.surveys$routePoint <- factor(paste(grouse.surveys$Route,
                                          grouse.surveys$Point.ID, sep="_"))

with(grouse.surveys,
     ave(seq_along(routePoint), routePoint, FUN=seq_along))

## Subset Coweeta data by sites surveyed during *standard* survey

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
formatLogisticRegData <- function(detArray, sites, species) {
    dim.dets <- dim(detArray)
    nYears <- dim.dets[3]
    yrs <- dimnames(detArray)[[3]]
    dets <- detArray[,,1,species]
    counts.yr1 <- rowSums(dets)
    dat <- data.frame(abundance=counts.yr1,
                      presence=ifelse(counts.yr1>0, 1L, 0L),
                      sites[rownames(dets),],
                      year=yrs[1])
    for(t in 2:nYears) {
        dets.t <- detArray[,,t,species]
        counts.yr.t <- rowSums(dets.t)
        dat.t <- data.frame(abundance=counts.yr.t,
                            presence=ifelse(counts.yr.t>0, 1L, 0L),
                            sites[rownames(dets.t),],
                            year=yrs[t])
        dat <- rbind(dat, dat.t)
    }
    colnames(dat) <- c("abundance", "presence", "utmN", "utmW", "elevation", "year")
    return(dat)
}
        


cawaData <- formatLogisticRegData(det.arr, sites, "CAWA")
howaData <- formatLogisticRegData(det.arr, sites, "HOWA")
veerData <- formatLogisticRegData(det.arr, sites, "VEER")
wothData <- formatLogisticRegData(det.arr, sites, "WOTH")

str(cawaData)


write.csv(cawaData, file="cawa_data_glm.csv")



fm1 <- glm(presence ~ elevation+I(elevation^2), binomial, cawaData)
fm2 <- glm(presence ~ year, binomial, cawaData)
fm3 <- glm(presence ~ elevation+I(elevation^2)+year, binomial, cawaData)
fm4 <- glm(presence ~ (elevation+I(elevation^2))*year, binomial, cawaData)
fm5 <- glm(presence ~ elevation*year, binomial, cawaData)
fm6 <- glm(presence ~ elevation+year, binomial, cawaData)

AIC(fm1, fm2, fm3, fm4, fm5, fm6)

summary(fm3)

elev.min <- min(cawaData$elevation)
elev.max <- max(cawaData$elevation)

yrs <- 2014:2020

pred.data <- data.frame(year=factor(rep(as.character(yrs), each=50)),
                        elevation=rep(seq(elev.min, elev.max, length=50),
                                      times=length(yrs)))

pred3 <- predict(fm3, newdata=pred.data, type="response", se.fit=TRUE)

pred3dat <- data.frame(pred=pred3$fit, pred.data)

plot(pred ~ elevation, pred3dat, subset=year=="2014",
     type="l", ylim=c(0, 1))
lines(pred ~ elevation, pred3dat, subset=year=="2015", col=2)
lines(pred ~ elevation, pred3dat, subset=year=="2016", col=3)
lines(pred ~ elevation, pred3dat, subset=year=="2017", col=4)
lines(pred ~ elevation, pred3dat, subset=year=="2018", col=5)
lines(pred ~ elevation, pred3dat, subset=year=="2019", col=6)
lines(pred ~ elevation, pred3dat, subset=year=="2020", col=7)
