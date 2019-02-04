################### Geospatial Analysis ###########################

#### Clean the working directory

rm(list=ls())


##### Load libraries

library(sp)
library(rgdal)
library(rgeos)
library(maptools)
library(dplyr)
library(broom)
library(data.table)
library(ggplot2)

#### Load location data ##################

# Set work directory----

#work.dir=("C:/Users/00097191/Google Drive/Analysis_GlobalArchive_Bodysize_Nestor")
work.dir=("G:/My drive/Analysis_GlobalArchive_Bodysize_Nestor") #Nestor
setwd("~/workspace")
dir()

# Set sub directories----

data.dir=paste(work.dir,"Data",sep="/")
plots.dir=paste(work.dir,"Plots",sep="/")
model.out=paste(work.dir,"ModelOut",sep="/")
functions.dir=paste(work.dir,"Functions",sep="/")

# Read in data----

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.csv")%>%
  dplyr::select(-X,-site,-location,-sst1961_20)%>%
  #filter(dst2ramp<150000)%>%
  filter(ausbath<60)%>%
  filter(depth<60)%>%
  #filter(dst2water<80000)%>%
  filter(relief<40)%>%
  filter(dst2mland<100000)%>%
  #filter(no3_sd<2.5)%>%
  #filter(t_sd<2.5)%>%
  #filter(t_m<26)%>%
  rename(response=number)%>%
  filter(Taxa=="Large")%>%
  mutate(status=ifelse(status=="Fished",0,1))%>%
  na.omit()%>%
  glimpse()

### Transform predictors

dat$gravity.200<-log1p(dat$gravity.200)
dat$pop.den.200<-sqrt(dat$pop.den.200)
dat$dst2ramp<-log1p(dat$dst2ramp)
dat$dst2coast<-log1p(dat$dst2coast)
dat$dst2water<-log1p(dat$dst2water)
dat$no3_m<-log1p(dat$no3_m)
dat$relief<-log1p(dat$relief)

## Check depth distribution

dat$ausbath<-(-dat$ausbath)
hist(dat$ausbath) ## Highly skewed - problably need to cut to the 60m isobar

## Code year as a factor - to include as a random effect

dat$year<-as.factor(dat$year)
dat$status<-as.factor(dat$status)
levels(dat$status)
plot(dat$status)

## Convert into an spatial object

coordinates(dat)<-~longitude+latitude
class(dat)

### Now we can set the reference system to the widely used WGS84

proj4string(dat)<-CRS("+init=epsg:4326")
dat@proj4string

#### Fit a top model
#### First I will do marginal predictions (no random effects) and not accounting for the offset (generates rates - really low numbers)

library(mgcv)

dat$sampling.deployment.duration.min.<-as.factor(dat$sampling.deployment.duration.min.)
levels(dat$sampling.deployment.duration.min.)

# Large
glimpse(dat)
Ra.gamm=gam(response ~ s(relief, k = 3, bs = "cr") + te(gravity.200, ausbath, k = 3, bs = "cr") + status + s(campaignid,bs="re") + s(sampling.deployment.duration.min.,bs="re"),
            family=tw(),data=dat)
plot(Ra.gamm,all.terms = TRUE,pages=1)
summary(Ra.gamm)
AIC(Ra.gamm)


## Get Australia coastline
setwd("C:/Users/22373243/Dropbox/Nestor_ascii")
dir()

australia<-readOGR('.','Australiaboundary67')
proj4string(australia)
australia <- spTransform(australia, CRS("+proj=longlat +datum=WGS84"))
australia<-fortify(australia)

############################# Rasters #####################

library(raster)

#### Read in rasters for the top large models ## Relief ## Depth ## Dst2ramp

## Large
dir()
relief<-raster("relief.asc")
ausbath<-raster("ausbath.asc")
gravity.200<-raster("r.gravity.asc")

relief
ausbath
gravity.200

### We need to fix ausbath raster because it gives weird values

ausbath[ausbath>0]<-NA
ausbath[ausbath<(-60)]<-NA

### Take a look at the layers

plot(relief)
plot(ausbath)
plot(gravity.200)

### Transform raster layers

relief<-log1p(relief)
gravity.200<-log1p(gravity.200)

plot(relief)
plot(gravity.200)


### Set coordinate reference system for gravity to WGS84

projection(gravity.200) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0 "

### At this stage you need to rasterize your polygons - to obtain grid cells stauts - 2 levels - Fished vs. No-take

## 1. Get Australia Marine Parks polygon shapefile
dir()
mpa<-readOGR('.','mpa_corrected')
class(mpa)
bbox(mpa)

### Now we need to rasterize
## We can rasterize to any of the rasters that we already have
## Check extent is the same - they are at rasterize them using that
## Taking the mask from relief raster

shp2raster <- function(shp, mask.raster, label, value, transform = FALSE, proj.from = NA,
                       proj.to = NA, map = TRUE) {
  require(raster, rgdal)
  
  # use transform==TRUE if the polygon is not in the same coordinate system as
  # the output raster, setting proj.from & proj.to to the appropriate
  # projections
  if (transform == TRUE) {
    proj4string(shp) <- proj.from
    shp <- spTransform(shp, proj.to)
  }
  
  # convert the shapefile to a raster based on a standardised background
  # raster
  r <- rasterize(shp, mask.raster)
  # set the cells associated with the shapfile to the specified value
  r[!is.na(r)] <- value
  # merge the new raster with the mask raster and export to the working
  # directory as a tif file
  r <- mask(merge(r, mask.raster), mask.raster, filename = label, format = "GTiff",
            overwrite = T)
  
  # plot map of new raster
  if (map == TRUE) {
    plot(r, main = label, axes = F, box = F)
  }
  
  names(r) <- label
  return(r)
}

## Set the background cells in the raster to 0

relief[!is.na(relief)]<-0

## MPAs to a raster

status <- shp2raster(shp = mpa,
                         mask.raster = relief, label = "status", value = 1,
                         transform = FALSE)

## Load relief again
setwd("C:/Users/22373243/Dropbox/Nestor_ascii")
dir()
relief<-raster("relief.asc")
relief<-log1p(relief)
plot(relief)

### Modify the values of the rasters
### Change status to original categories

# relief[relief>40]<-NA

names(gravity.200) <- c('gravity.200')
names(gravity.200)

names(relief) <- c('relief')
names(relief)

### Combine the 4 layers -- ausbath, relief, gravity.200 and status into a RasterStack object

predict.large<-stack(ausbath,relief,gravity.200,status)
predict.large@layers
levels(dat$campaignid)

## Predict to stack raster layer for each year

prediction.large<-predict(predict.large,Ra.gamm, type = 'response', se.fit=T, index=1:2, progress='text',const=(data.frame(campaignid='2015-Summer-PSGLMP-StereoBRUVs')))

## Export predictions as a raster layer

writeRaster(prediction.large,'prediction.large.asc')

## Convert raster to a data frame

prediction.large_df<-as.data.frame(prediction.large,xy=TRUE,na.rm=TRUE)
glimpse(prediction.large_df)
colnames(prediction.large_df)<-c("Longitude","Latitude","Relative abundance (MaxN)","SE")
glimpse(prediction.large_df)

### Plot the predictions and CI

min(prediction.large_df$`Relative abundance (MaxN)`)
max(prediction.large_df$`Relative abundance (MaxN)`)


prediction.large_df$`Relative abundance (MaxN)`<-log1p(prediction.large_df$`Relative abundance (MaxN)`)

### Set the colour Rampalette we want for plotting

library(viridis)
library(ggplot2)

ggplot()+
  geom_tile(data=prediction.large_df,aes(x=Longitude,y=Latitude,fill=`Relative abundance (MaxN)`),alpha=0.8)+
  geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="lightgrey")+
  scale_fill_viridis(option = "magma",direction = -1)+
  coord_equal()+
  xlim(110,160)+
  ylim(-45,-9)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA,size = 1))


ggplot()+
  geom_tile(data=prediction.large_df,aes(x=Longitude,y=Latitude,fill=CI),alpha=0.8)+
  geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="lightgrey")+
  scale_fill_viridis(option = "magma",direction = -1)+
  coord_equal()+
  xlim(110,160)+
  ylim(-45,-9)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA,size = 1))

large<-as.data.frame(large)


ggplot()+
  geom_point(data=large,aes(x=longitude,y=latitude,colour=status),alpha=0.8)+
  geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="lightgrey")+
  geom_polygon(data=mpa,aes(x=long,y=lat,group=group),fill="green")+
  coord_equal()+
  xlim(110,160)+
  ylim(-45,-9)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA,size = 1))
