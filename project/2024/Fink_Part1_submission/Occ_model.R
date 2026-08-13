###############################
# Fink Estimation Project
# Deer Dynamic Occupancy Model
###############################

#Read in data

d<- read.csv("./deer_occ_combind.csv")
covs<- read.csv('./covariates_std.csv')
obsCovs<- read.csv('./doy.csv')

#drop site name for occ matrix
d<- d[,-1]

#change obscovs to correct format for unmarked
obs<- list(doy = as.matrix(obsCovs))
print(dim(obs)) 

# set primary period and sampling days
library(unmarked)
nPrimary<- 3
n  <- 30

#create unmarked frame
umf <- unmarkedMultFrame(y=d, numPrimary=nPrimary, siteCovs=data.frame(covs),
                          obsCovs=obs)
summary(umf)

#fit NULL model
fm <- colext(~1,~1,~1,~1, umf)    
fm


#fit detection model
m_doy <- colext(~1,~1,~1,~doy, umf)   
m_doy

#Univariate Models with doy for detection
#income
m_inc <- colext(~Ave_med_incom_std,~1,~1,~doy, umf)    
m_inc

#race
m_race <- colext(~prop_white_std,~1,~1,~doy, umf)    
m_race

#cover
m_cover <- colext(~cover_std,~1,~1,~doy, umf)    
m_cover

#canopy
m_canopy <- colext(~canopy_std,~1,~1,~doy, umf)    
m_canopy

#imperv
m_imperv <- colext(~imprev_std,~1,~1,~doy, umf)    
m_imperv

#roads
m_roads <- colext(~roads_std,~1,~1,~doy, umf)    
m_roads
  

#best model was included roads so built on rds and
#only included covariates that were not correlated

#roads_cover
m_roads_cover <- colext(~roads_std+cover_std,~1,~1,~doy, umf)    
m_roads_cover

#roads_inc
m_roads_inc <- colext(~roads_std+Ave_med_incom_std,~1,~1,~doy, umf)    
m_roads_inc

#roads_race
m_roads_race <- colext(~roads_std+prop_white_std,~1,~1,~doy, umf)    
m_roads_race

################
# so deer was not the best species to make my point 
# all they care about are roads...

# will probably change doy to temperature 