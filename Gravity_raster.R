########### Creating a gravity layer for Prediction Across Australia #############
###### Created by Nestor E. Bosch ###############
###### Any bugs report to nbosch1989@gmail.com ###########

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
library(raster)


##### Create a raster layer #########
##### Use the initial predictive space - from dst2ramp raster

setwd("C:/Users/22373243/Dropbox/Nestor_ascii")
setwd("~/workspace")
dir()

### Load dst2ramp raster

dst2ramp<-raster("dst2ramp.asc")
plot(dst2ramp)
dst2ramp

### Create an empty raster with the same extent as dst2ramp

pop.den<-raster(dst2ramp)
pop.den

### Convert dst2ramp to dataframe

pop.den_df<-as.data.frame(dst2ramp,xy=TRUE,na.rm=T)%>%
  dplyr::select(-dst2ramp)%>%
  glimpse()


### Bring in LandsCan 2011 population density grid - in WGS84

setwd("C:/Users/22373243/Dropbox/Nestor_ascii/LandScan Global 2011/lspop2011.folder")
dir()
pop.den.data<-raster("dblbnd.adf")
pop.den.data<-crop(pop.den.data,extent(109.2207,156.6676,-44.14824,-8.179126))
projection(pop.den.data)
plot(pop.den.data)

### Aggregate cell values in a 200 km buffer - fact = 200

pop.den.data<-aggregate(pop.den.data,fact=200,fun=sum,na.rm=T)
plot(pop.den.data)
pop.den.data


### Extract population density values

points.200<-extract(pop.den.data,pop.den_df,df=T)
glimpse(points.200)

### Include values in pop.density_df

pop.den_df$pop.den.200<-points.200$dblbnd
glimpse(pop.den_df)

### Rasterize
### First need to convert to a spatialpoints object

coordinates(pop.den_df)<-~x+y
class(pop.den_df)

r.pop.den.200<-rasterize(pop.den_df, dst2ramp, pop.den_df$pop.den.200,fun=mean)
plot(r.pop.den.200)

### Export population density raster

writeRaster(r.pop.den.200,'r.pop.den.200.asc')

### Clean working directory

rm(dst2ramp,points.200,pop.den,pop.den_df,pop.den.data)

### Bring in dst2townc raster

dir()
dst2townc<-raster("dst2townc.asc")
dst2townc
plot(dst2townc)

dst2townc<-dst2townc/1000
plot(dst2townc)

### Create gravity layer - population density/distance to town centre^2

r.gravity<-r.pop.den.200/dst2townc^2
plot(r.gravity)

writeRaster(r.gravity,'r.gravity.asc')


