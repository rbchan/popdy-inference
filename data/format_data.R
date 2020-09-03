list.files()


## ------------------- Import Coweeta data ---------------------------

sites.in <- read.csv("Site_Data_WILD8390.csv")
dets.in <- read.csv("Detection_Data_WILD8390.csv")
surveys.in <- read.csv("Survey_Covs_WILD8390.csv")

str(sites.in)
str(dets.in)
str(surveys.in)


## -------------------- Import grouse data ------------------------


## For grouse.locs.in, the proj information is:
## NAD83 UTM Zone 16N
## GCS83 


grouse.dets.in <- read.csv("GrouseDetectionsQuery.csv")
grouse.surveys.in <- read.csv("RuffedGrouseDrummingSurveys Query.csv")
grouse.locs.in <- read.csv("RUGR_Survey_Locations.csv")

str(grouse.dets.in)
str(grouse.surveys.in)
str(grouse.locs.in)


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

grouse.locs <- grouse.locs.in[,1:5]
grouse.locs$routePoint <- factor(sub("-", "_", grouse.locs$ident))
rownames(grouse.locs) <- grouse.locs$routePoint

setdiff(grouse.locs$routePoint, grouse.surveys$routePoint)
setdiff(grouse.surveys$routePoint, grouse.locs$routePoint)


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

keepvars <- c("Route..", ##"Coordinates..easting.", "Coordinates..northing.",
              "UTM.Zone")

grouse.data <- cbind(grouse.data,
                     grouse.surveys.wide[reorder.surveys,keepvars])
colnames(grouse.data) <- c("abundance", "presence", "route", "utmE", "utmN",
                           "utmZone")

str(grouse.data)


grouse.data$utmE <- grouse.locs[match(rownames(grouse.data),
                                      grouse.locs$routePoint), "x_proj"]
grouse.data$utmN <- grouse.locs[match(rownames(grouse.data),
                                      grouse.locs$routePoint), "y_proj"]

str(grouse.data)




plot(grouse.data[,c("utmE", "utmN")], asp=1)
points(grouse.locs[,c("x_proj", "y_proj")], pch=16, col=2)
points(grouse.data[,c("utmE", "utmN")])




## Spatial data
## install.packages("USAboundaries")
## install.packages("elevatr")

library(USAboundaries)
library(elevatr)
library(sf)
library(rgdal)
library(sp)
library(raster)
library(lme4)


## The grouse data are a bit tricky because the study area spans two UTM Zones
## cs <- make_EPSG()
## grep("zone=16", cs$prj4, value=TRUE)
grep("zone=16", cs$prj4, value=TRUE)

## proj4strings
## utm.z16 <- "+proj=utm +zone=16S +datum=WGS84 +units=m +no_defs +type=crs"
## utm.z17 <- "+proj=utm +zone=17S +datum=WGS84 +units=m +no_defs +type=crs"
utm.z16 <- "+proj=utm +zone=16S +datum=NAD83 +units=m +no_defs +type=crs"
utm.z17 <- "+proj=utm +zone=17S +datum=NAD83 +units=m +no_defs +type=crs"
longlat <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"

us.states <- us_states()
ga.nc.sc.tn <- st_geometry(us.states)[us.states$stusps %in%
                                      c("GA", "NC", "SC", "TN")]

ga.nc.sc.tn.utm <- st_transform(ga.nc.sc.tn, crs=CRS(utm.z16))


## We have a problem
## plot(Coordinates..northing. ~ Coordinates..easting., grouse.surveys.in,
##      asp=1, pch=3)
## plot(ga.nc.sc.tn.utm, add=TRUE)




## st_crs(ga.nc.sc.tn)
## cs[cs$code==4326,]

## Convert to SpatialPoints class
grouse.coords.z16 <- SpatialPoints(grouse.data[!is.na(grouse.data$utmE),
                                               c("utmE", "utmN")],
                                   proj4string=CRS(utm.z16))
grouse.coords.z17 <- SpatialPoints(grouse.data[,c("utmE", "utmN")],
                                   proj4string=CRS(utm.z17))

plot(ga.nc.sc.tn.utm, axes=TRUE)
points(grouse.coords.z16)



## Project to longlat
grouse.coords.longlat <- coordinates(spTransform(grouse.coords.z16, CRS(longlat)))
grouse.coords.longlat[grouse.data$utmZone=="17S",] <- 
    coordinates(spTransform(grouse.coords.z17, CRS(longlat)))[grouse.data$utmZone=="17S",]

all(rownames(grouse.data)==rownames(grouse.coords.longlat)) ## Must be TRUE


grouse.coords.utm17 <- spTransform(SpatialPoints(grouse.coords.longlat,
                                                 proj4string=CRS(longlat)),
                                   CRS(utm.z17))

grouse.coords.utm16 <- spTransform(SpatialPoints(grouse.coords.longlat,
                                                 proj4string=CRS(longlat)),
                                   CRS(utm.z16))



coordinates(grouse.coords.z17)-coordinates(grouse.coords.utm17)
coordinates(grouse.coords.z16)-coordinates(grouse.coords.utm16)

plot(ga.nc.sc.tn.utm, axes=TRUE)
points(grouse.coords.utm17)

plot(ga.nc.sc.tn.utm, xlim=c(2e5,4e5), ylim=c(37e5,4e6))
points(grouse.coords.utm17)


plot(ga.nc.sc.tn)
points(grouse.coords.longlat)


plot(grouse.coords.longlat, asp=1)
points(grouse.coords.longlat[grouse.data$utmZone=="16S",], pch=16, col=2)
points(grouse.coords.longlat[grouse.data$utmZone=="17S",], pch=16, col=4)
plot(ga.nc.sc.tn, add=TRUE)


plot(grouse.coords.longlat, asp=1)
points(grouse.coords.longlat, pch=16, cex=grouse.data$abundance*2, col=rgb(0,0,1,0.5))
plot(ga.nc.sc.tn, add=TRUE)

plot(grouse.coords.utm17, asp=1)
points(grouse.coords.utm17, pch=16, cex=grouse.data$abundance*2, col=rgb(0,0,1,0.5))
plot(ga.nc.sc.tn.utm, add=TRUE)


plot(grouse.coords.utm16, asp=1)
points(grouse.coords.utm16, pch=16, cex=grouse.data$abundance*2, col=rgb(0,0,1,0.5))
plot(ga.nc.sc.tn.utm, add=TRUE)


with(grouse.data, {
    plot(utmN ~ utmE, asp=1, pch=3)
    points(utmN ~ utmE, pch=16, cex=abundance, col=rgb(0,0,1,0.5))
##    points(utmN ~ utmE, subset=utmZone=="16S", pch=16, cex=2, col=2)
    plot(ga.nc.sc.tn.utm, add=TRUE)
})





dir.create("../lectures/stats-basics/figs")

pdf("../lectures/stats-basics/figs/grouse_map_locs_longlat.pdf", width=6, height=5)
par(mai=c(0.7,0.8,0.1,0.1))
plot(ga.nc.sc.tn, axes=TRUE, xlim=c(-85.6,-82.9), ylim=c(34.5, 35.5), las=1)
points(grouse.coords.longlat, pch=3, cex=0.5, col=gray(0.7))
text(-83.5, 34.4, "Georgia", pos=1)
text(-83.5, 35.2, "North Carolina", pos=1)
text(-85, 35.2, "Tennessee", pos=1)
dev.off()
system("gopen ../lectures/stats-basics/figs/grouse_map_locs_longlat.pdf")


pdf("../lectures/stats-basics/figs/grouse_map_locs.pdf", width=6, height=5)
par(mai=c(0.7,0.9,0.1,0.1))
plot(ga.nc.sc.tn.utm, axes=TRUE, xlim=c(71e4, 85e4), ylim=c(378e4, 39e5), las=1)
points(grouse.data[,c("utmE", "utmN")], pch=3, cex=0.5, col=gray(0.7))
text(8e5, 381e4, "Georgia", pos=1)
text(8e5, 389e4, "North Carolina", pos=1)
text(72e4, 389e4, "Tennessee", pos=1)
dev.off()
system("gopen ../lectures/stats-basics/figs/grouse_map_locs.pdf")


pdf("../lectures/stats-basics/figs/grouse_map_locs_dets_longlat.pdf", width=6, height=5)
par(mai=c(0.7,0.8,0.1,0.1))
plot(ga.nc.sc.tn, axes=TRUE, xlim=c(-85,-83), ylim=c(34.6, 35.2), las=1)
points(grouse.coords.longlat, pch=3, cex=0.5, col=gray(0.7))
points(grouse.coords.longlat, pch=16, cex=grouse.data$abundance*2, col=rgb(0,0,1,0.5))
text(-83.5, 34.5, "Georgia", pos=1)
text(-83.5, 35.2, "North Carolina", pos=1)
text(-84.7, 35.2, "Tennessee", pos=1)
dev.off()
system("gopen ../lectures/stats-basics/figs/grouse_map_locs_dets_longlat.pdf")


pdf("../lectures/stats-basics/figs/grouse_map_locs_dets.pdf", width=6, height=5)
par(mai=c(0.7,0.9,0.1,0.1))
plot(ga.nc.sc.tn.utm, axes=TRUE, xlim=c(71e4, 85e4), ylim=c(378e4, 39e5), las=1)
points(grouse.data[,c("utmE", "utmN")], pch=3, cex=0.5, col=gray(0.7))
points(grouse.data[,c("utmE", "utmN")], pch=16, cex=grouse.data$abundance*2,
       col=rgb(0,0,1,0.5))
text(8e5, 381e4, "Georgia", pos=1)
text(8e5, 389e4, "North Carolina", pos=1)
text(72e4, 389e4, "Tennessee", pos=1)
dev.off()
system("gopen ../lectures/stats-basics/figs/grouse_map_locs_dets.pdf")


##rar <- raster(extent(-85, -83, 34.6, 35.2), crs=longlat, res=0.001)
rar <- raster(extent(70e4, 88e4, 378e4, 39e5), crs=utm.z16, res=100)
## rar[] <- 0


region.elev <- get_elev_raster(rar, z=7)#,
##prj="+proj=longlat +datum=WGS84 +no_defs", ##src="aws",
#                               clip="bbox")

region.elev <- crop(region.elev, rar)

plot(region.elev)




pdf("../lectures/stats-basics/figs/grouse_map_elev_locs_dets_longlat.pdf", width=12.4, height=5)
par(mai=c(0.7,0.8,0.2,0.1))
plot(region.elev, xlim=c(-84.9,-83), ylim=c(34.4, 35.2))
plot(ga.nc.sc.tn, add=TRUE)#, las=1)
points(grouse.coords.longlat, pch=3, cex=0.5, col=2)
points(grouse.coords.longlat, pch=16, cex=grouse.data$abundance*2, col=rgb(0,0,1,0.5))
text(-83.5, 34.7, "Georgia", pos=1)
text(-83.5, 35.2, "North Carolina", pos=1)
text(-84.7, 35.2, "Tennessee", pos=1)
dev.off()
system("gopen ../lectures/stats-basics/figs/grouse_map_elev_locs_dets_longlat.pdf")


pdf("../lectures/stats-basics/figs/grouse_map_elev_locs_dets.pdf", width=8.3, height=5)
par(mai=c(0.7,0.9,0.1,0.1))
plot(region.elev, xlim=c(70e4, 87e4), ylim=c(379e4, 39e5), las=1)
plot(ga.nc.sc.tn.utm, add=TRUE)
points(grouse.data[,c("utmE", "utmN")], pch=3, cex=0.5, col=gray(0.7))
points(grouse.data[,c("utmE", "utmN")], pch=16, cex=grouse.data$abundance*2,
       col=rgb(0,0,1,0.5))
text(8e5, 381e4, "Georgia", pos=1)
text(8e5, 389e4, "North Carolina", pos=1)
text(73e4, 389e4, "Tennessee", pos=1)
dev.off()
system("gopen ../lectures/stats-basics/figs/grouse_map_elev_locs_dets.pdf")




## grouse.elev <- get_elev_point(as.data.frame(grouse.coords.longlat),
##                               prj=longlat, src="aws")
grouse.elev <- get_elev_point(as.data.frame(grouse.data[!is.na(grouse.data$utmE),
                                                        c("utmE","utmN")]),
                              prj=utm.z16, src="aws")

grouse.data$elevation <- NA_real_
grouse.data$elevation[!is.na(grouse.data$utmE)] <- grouse.elev@data$elevation

str(grouse.data)






write.csv(grouse.data, file="grouse_data_glm.csv")

plot(grouse.data$elev, grouse.data$abund)

with(grouse.data, plot(abundance ~ elevation))
with(grouse.data, boxplot(abundance ~ utmZone))


gm1 <- glm(presence ~ elevation, binomial, grouse.data)
gm2 <- glm(presence ~ utmZone, binomial, grouse.data)
gm3 <- glm(presence ~ elevation*utmZone, binomial, grouse.data)
gm4 <- glm(presence ~ elevation+I(elevation^2), binomial, grouse.data)
gm5 <- glm(presence ~ elevation+I(elevation^2)+utmZone, binomial,
           grouse.data)

AIC(gm1,gm2,gm3,gm4,gm5)


summary(gm1)
summary(gm2)
summary(gm3)
summary(gm4)
summary(gm5)



grouse.data$route <- factor(grouse.data$route)


gmm1 <- glmer(presence ~ scale(elevation)+utmZone+(1|route),
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
formatGLMData <- function(detArray, sites, species) {
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
        


cawa.data <- formatGLMData(det.arr, sites, "CAWA")
howa.data <- formatGLMData(det.arr, sites, "HOWA")
veer.data <- formatGLMData(det.arr, sites, "VEER")
woth.data <- formatGLMData(det.arr, sites, "WOTH")

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
