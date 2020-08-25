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

str(dets)
