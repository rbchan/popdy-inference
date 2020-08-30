list.files()


## ------------------- Import Coweeta data ---------------------------

sites.in <- read.csv("Site_Data_WILD8390.csv")
dets.in <- read.csv("Detection_Data_WILD8390.csv")
surveys.in <- read.csv("Survey_Covs_WILD8390.csv")

str(sites.in)
str(dets.in)
str(surveys.in)


## -------------------- Import grouse data ------------------------

grouse.dets.in <- read.csv("GrouseDetectionsQuery.csv")
grouse.surveys.in <- read.csv("RuffedGrouseDrummingSurveys Query.csv")

str(grouse.dets.in)
str(grouse.surveys.in)



## -------------------- Format grouse data ------------------------


## Create a unique ID for each stop on each route
grouse.surveys <- grouse.surveys.in
grouse.surveys$routePoint <- factor(paste(grouse.surveys$Route,
                                          grouse.surveys$Point.ID, sep="_"))

## Ditto, but use factor levels from above for routePoint
grouse.dets <- grouse.dets.in
grouse.dets$routePoint <- factor(paste(grouse.dets.in$Route,
                                       grouse.dets.in$Point.ID, sep="_"),
                                 levels=levels(grouse.surveys$routePoint))



## Create occasion index
grouse.surveys$occasion <- with(grouse.surveys,
     ave(seq_along(routePoint), routePoint, FUN=seq_along))


grouse.surveys.wide <- reshape(
    grouse.surveys, idvar="routePoint", timevar="occasion",
    v.names=c("Date", "Time", "Observer", "Cloud.Cover", "Temperature",
              "Wind.Speed", "Disturbance", "Precipitation",
              "Ericaceous.Cover"),
    drop=c("ID", "Sunrise", "Surveyed.Songbirds", "Notes", "Proofed."),
    direction="wide")
                               

str(grouse.surveys.wide)


## Combine detection and survey data
grouse.counts <- unclass(table(grouse.dets$routePoint))

grouse.data <- data.frame(abundance=grouse.counts,
                         presence=ifelse(grouse.counts>0, 1L, 0L))
rownames(grouse.data) <- names(grouse.counts)

reorder.surveys <- match(rownames(grouse.data),
                         grouse.surveys.wide$routePoint)

keepvars <- c("Route..", "Coordinates..easting.", "Coordinates..northing.",
              "UTM.Zone")

grouse.data <- cbind(grouse.data,
                     grouse.surveys.wide[reorder.surveys,keepvars])
colnames(grouse.data) <- c("abundance", "presence", "route", "utmE", "utmN", "utmZone")

str(grouse.data)


with(grouse.data, plot(utmN ~ utmE, asp=1, pch=3))


## Spatial data
## install.packages("USAboundaries")
## install.packages("elevatr")

library(USAboundaries)
library(elevatr)
library(sf)
library(rgdal)

cs <- make_EPSG()

grep("zone=16", cs$prj4, value=TRUE)

utm.z16 <- "+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs +type=crs"
utm.z17 <- "+proj=utm +zone=17 +datum=WGS84 +units=m +no_defs +type=crs"

grouse.elev <- get_elev_point(grouse.data[,c("utmE","utmN")],
                              prj=utm.z17, src="aws")

grouse.data$elevation <- grouse.elev@data$elevation
grouse.data$zone

str(grouse.data)



## TODO project to UTM
us.states <- us_states()
ga.nc.sc.tn <- st_geometry(us.states)[us.states$stusps %in%
                                      c("GA", "NC", "SC", "TN")]

plot(ga.nc.sc.tn)
with(grouse.data, points(utmN ~ utmE, pch=3, col=4))



write.csv(grouse.data, file="grouse_data_glm.csv")


with(grouse.data, plot(abundance ~ elevation))
with(grouse.data, boxplot(abundance ~ utmZone))


gm1 <- glm(presence ~ elevation, binomial, grouse.data)
gm2 <- glm(presence ~ utmZone, binomial, grouse.data)
gm3 <- glm(presence ~ elevation*utmZone, binomial, grouse.data)
gm4 <- glm(presence ~ elevation+I(elevation^2), binomial, grouse.data)
gm5 <- glm(presence ~ elevation+I(elevation^2)+utmZone, binomial,
           grouse.data)

summary(gm1)
summary(gm2)
summary(gm3)
summary(gm4)
summary(gm5)


library(lme4)

grouse.data$route <- factor(grouse.data$route)


gmm1 <- glmer(presence ~ scale(elevation)+I(scale(elevation)^2)+(1|route),
              data=grouse.data,
              family=binomial)
summary(gmm1)




## -------------------- Format Coweeta data ------------------------


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


## Function to format logistic regression data for selected species
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
        


cawa.data <- formatLogisticRegData(det.arr, sites, "CAWA")
howa.data <- formatLogisticRegData(det.arr, sites, "HOWA")
veer.data <- formatLogisticRegData(det.arr, sites, "VEER")
woth.data <- formatLogisticRegData(det.arr, sites, "WOTH")

str(cawa.data)


write.csv(cawa.data, file="cawa_data_glm.csv")



fm1 <- glm(presence ~ elevation+I(elevation^2), binomial, cawa.data)
fm2 <- glm(presence ~ year, binomial, cawa.data)
fm3 <- glm(presence ~ elevation+I(elevation^2)+year, binomial, cawa.data)
fm4 <- glm(presence ~ (elevation+I(elevation^2))*year, binomial, cawa.data)
fm5 <- glm(presence ~ elevation*year, binomial, cawa.data)
fm6 <- glm(presence ~ elevation+year, binomial, cawa.data)

AIC(fm1, fm2, fm3, fm4, fm5, fm6)

summary(fm3)

elev.min <- min(cawa.data$elevation)
elev.max <- max(cawa.data$elevation)

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
