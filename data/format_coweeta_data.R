list.files()


## ------------------- Import Coweeta data ---------------------------

sites.in <- read.csv("Site_Data_WILD8390.csv")
dets.in <- read.csv("Detection_Data_WILD8390.csv")
surveys.in <- read.csv("Survey_Covs_WILD8390.csv")

str(sites.in)
str(dets.in)
str(surveys.in)




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

surveys$Date <- as.Date(surveys$SurvDate, format="%Y-%m-%d")
surveys <- surveys[order(surveys$PointName, surveys$Date),]
str(surveys)



## Create occasion index
## surveys.2017 <- subset(surveys, format(Date, format="%Y")=="2017")

## surveys.2017$occasion <- with(surveys.2017,
##      ave(seq_along(PointName), PointName, FUN=seq_along))
## str(surveys.2017)
## surveys.2017.wide <- reshape(surveys.2017, idvar="PointName",
##                              timevar="occasion", direction="wide")

## str(surveys.2017.wide)

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



## Export 2017 CAWA data
cawa.data.2017 <- ifelse(unclass(det.arr[,,"2017","CAWA"])>1, 1L, 0L)

surveys.2017 <- subset(surveys, format(Date, format="%Y")=="2017" &
                                PointName %in% rownames(cawa.data.2017))
surveys.2017$occasion <- with(surveys.2017,
                              ave(seq_along(PointName), PointName, FUN=seq_along))
surveys.2017 <- surveys.2017[surveys.2017$occasion==1,]


## surveys.2017$occasion <- with(surveys.2017,
##      ave(seq_along(PointName), PointName, FUN=seq_along))


    
match(surveys.2017$PointName, rownames(cawa.data.2017))

cawa.data.2017 <- cbind(cawa.data.2017,
                        surveys.2017[,c("SurvDate","SurvTime","Precipitation",
                                        "Wind","Noise")])
cawa.data.2017$Elevation <- sites.in[match(rownames(cawa.data.2017),
                                           sites.in$PointName),"Elevation"]
colnames(cawa.data.2017)[1:4] <- c("cawa1","cawa2","cawa3","cawa4")

str(cawa.data.2017)

write.csv(cawa.data.2017, file="cawa_data_2017_occu.csv")








## Export 2017 CAWA count data
cawa.counts.2017 <- unclass(det.arr[,,"2017","CAWA"])


cawa.counts.2017 <- cbind(cawa.counts.2017,
                        surveys.2017[,c("SurvDate","SurvTime","Precipitation",
                                        "Wind","Noise")])
cawa.counts.2017$Elevation <- sites.in[match(rownames(cawa.counts.2017),
                                             sites.in$PointName),"Elevation"]
colnames(cawa.counts.2017)[1:4] <- c("cawa1","cawa2","cawa3","cawa4")

str(cawa.counts.2017)

write.csv(cawa.counts.2017, file="cawa_data_2017_binNmix.csv")







## Format distance sampling data

library(unmarked)

str(dets)

dets.btbw.2019 <- subset(dets, species=="BTBW" & year=="2019")

ds.btbw <- formatDistData(dets.btbw.2019, distCol="Distance",
                          transectNameCol="PointName",
                          dist.breaks=seq(0,100,20))


all(rownames(sites) == rownames(ds.btbw)) ## Must be TRUE

surveys.2019 <- surveys[format(surveys$Date, "%Y")=="2019",]

all(rownames(sites) == surveys.2019$PointName) ## Must be TRUE



dsdat.btbw <- cbind(sites, ds.btbw, surveys.2019)
dsdat.btbw <- dsdat.btbw[,!(colnames(dsdat.btbw) %in% c("PointName", "Date"))]

str(dsdat.btbw)


write.csv(dsdat.btbw, file="btbw_data_distsamp.csv")



umf <- unmarkedFrameDS(y=data.matrix(dsdat.btbw[,4:8]),
                       survey="point", dist.breaks=seq(0, 100, 20),
                       unitsIn="m")

distsamp(~1~1, umf)

