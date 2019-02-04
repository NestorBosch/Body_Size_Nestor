################### Model Performance ####################################
### Refers to methods for measuring the performance of a given predictive model on new test datasets
### The basic idea consist of dividing the dataset into two sets
### 1. The training set - used to train the model
### 2. Testing or validation set - used to validate the model - by estimating predictive error

### Basic strategy consist on -
### 1. Build the model on a training dataset
### 2. Apply the model on a new test data set to make predictions
### 3. Compute the prediction errors

### Several statistical metrics to quantify the overall quality of regression models
### 1. R2 - squared correlation between the observed outcome values and the predicted values by the model - the higher the better the model
### 2. Root Mean Squared Error (RMSE) - average prediction error made by the model in predicting the outcome for an observation - the lower the better
### 3. Mean Absolute Error (MAE) - less sensitive to outliers - the lower the better


### Method - Repeated k-fold cross-validation
## The final model error is taken as the mean error from the number of repeats
## We will use 5-fold cross validation with 50 repeats



### Load libraries

library(tidyverse)
library(caret)
library(tidyr)
library(dplyr)
library(forcats)
library(mgcv)
library(MuMIn)
library(car)
library(doBy)
library(gplots)
library(RColorBrewer)
library(doParallel)
library(gamm4)
library(RCurl)#needed to download data from GitHub
library(doSNOW)

### Read in the data

rm(list=ls())

# GAMM Size models----

name<-"Size"


# Set work directory----

#work.dir=("G:/My drive/Analysis_GlobalArchive_Bodysize_Nestor") #Nestor
work.dir=("~/workspace/Body_Size_Nestor") ### Ecocloud platform

# Set sub directories----

data.dir=paste(work.dir,"Data",sep="/")
plots.dir=paste(work.dir,"Plots",sep="/")
model.out=paste(work.dir,"ModelOut",sep="/")
functions.dir=paste(work.dir,"Functions",sep="/")
layers.dir=paste(work.dir,"Layers",sep="/")

# Read in data----

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.csv")%>%
  select(-X,-site,-location,-sst1961_20)%>%
  #filter(dst2ramp<150000)%>%
  filter(ausbath<100)%>%
  filter(depth<100)%>%
  #filter(dst2water<80000)%>%
  #filter(relief<40)%>%
  #filter(dst2mland<100000)%>%
  #filter(no3_sd<2.5)%>%
  #filter(t_sd<2.5)%>%
  #filter(t_m<26)%>%
  rename(response=number)%>%
  #filter(Taxa=="Large")%>%
  na.omit()%>%
  glimpse()


## Code year and soak time as a factor - to include as a random effect

dat$year<-as.factor(dat$year)
dat$sampling.deployment.duration.min.<-as.numeric(dat$sampling.deployment.duration.min.)
plot(dat$sampling.deployment.duration.min.)

## Transform predictors

dat$gravity.200<-log1p(dat$gravity.200)
dat$pop.den.200<-sqrt(dat$pop.den.200)
dat$dst2ramp<-log1p(dat$dst2ramp)
dat$dst2coast<-log1p(dat$dst2coast)
dat$dst2water<-log1p(dat$dst2water)
dat$no3_m<-log1p(dat$no3_m)
dat$relief<-log1p(dat$relief)

## Create Large, Medium and Small data sets

large<-dat%>%
  filter(Taxa=="Large")%>%
  glimpse()

medium<-dat%>%
  filter(Taxa=="Medium")%>%
  glimpse()

small<-dat%>%
  filter(Taxa=="Small")%>%
  glimpse()

#### Investigating patterns in the residuals

###### Large - relief + gravity.te.ausbath + status

### Fit the GAM model to the data set

La.gamm=gam(response ~ s(relief, k = 3, bs = "cr") + te(gravity.200, ausbath, k = 3, bs = "cr") + status + s(campaignid,bs="re"),
            family=tw(),offset = log(sampling.deployment.duration.min.) ,data=large)
plot(La.gamm,all.terms = TRUE,pages=1)
summary(La.gamm)
AIC(La.gamm)
par(mfrow=c(1,1))
gam.check(La.gamm)

## Set working directory

setwd(plots.dir)

jpeg("Residuals_large_no_filter.jpg", width = 300, height = 300)
par(mfrow=c(3,3))

hist(residuals(La.gamm))

plot(large$response~ fitted(La.gamm))
abline(h=0, lty="dotted")
lines(lowess(fitted(La.gamm), large$response), col="red")

plot(residuals(La.gamm)~ fitted(La.gamm))
abline(h=0, lty="dotted")
lines(lowess(fitted(La.gamm), residuals(La.gamm)), col="red")

plot(large$relief,residuals(La.gamm))
abline(h=0, lty="dotted")
lines(lowess(large$relief, residuals(La.gamm)), col="red")

plot(large$gravity.200,residuals(La.gamm))
abline(h=0, lty="dotted")
lines(lowess(large$gravity.200, residuals(La.gamm)), col="red")

plot(large$ausbath,residuals(La.gamm))
abline(h=0, lty="dotted")
lines(lowess(large$ausbath, residuals(La.gamm)), col="red")

plot(large$status,residuals(La.gamm))

dev.off()

###### Medium - relief + gravity.te.ausbath + status

### Fit the GAM model to the data set

Ma.gamm=gam(response ~ s(relief, k = 3, bs = "cr") + te(gravity.200, ausbath, k = 3, bs = "cr") + status + s(campaignid,bs="re"),
            family=tw(),offset = log(sampling.deployment.duration.min.) ,data=medium)
plot(Ma.gamm,all.terms = TRUE,pages=1)
summary(Ma.gamm)
AIC(Ma.gamm)
par(mfrow=c(1,1))
gam.check(Ma.gamm)

## Set working directory

setwd(plots.dir)

jpeg("Residuals_medium_no_filter.jpg", width = 300, height = 300)
par(mfrow=c(3,3))

hist(residuals(Ma.gamm))

plot(medium$response~ fitted(Ma.gamm))
abline(h=0, lty="dotted")
lines(lowess(fitted(Ma.gamm), medium$response), col="red")

plot(residuals(Ma.gamm)~ fitted(Ma.gamm))
abline(h=0, lty="dotted")
lines(lowess(fitted(Ma.gamm), residuals(Ma.gamm)), col="red")

plot(medium$relief,residuals(Ma.gamm))
abline(h=0, lty="dotted")
lines(lowess(medium$relief, residuals(Ma.gamm)), col="red")

plot(medium$gravity.200,residuals(Ma.gamm))
abline(h=0, lty="dotted")
lines(lowess(medium$gravity.200, residuals(Ma.gamm)), col="red")

plot(medium$ausbath,residuals(Ma.gamm))
abline(h=0, lty="dotted")
lines(lowess(medium$ausbath, residuals(Ma.gamm)), col="red")

plot(medium$status,residuals(Ma.gamm))

dev.off()


###### Small

### Fit the GAM model to the data set

Sa.gamm=gam(response ~ s(t_sd, k = 3, bs = "cr") + s(ausbath,by=status, k = 3, bs = "cr") + s(dst2ramp,by=status, k = 3, bs = "cr") + status + s(campaignid, bs = "re"),
            family=tw(),offset = log(sampling.deployment.duration.min.),data=small)
plot(Sa.gamm,all.terms = T,pages = 1)
summary(Sa.gamm)
gam.check(Sa.gamm)

## Set working directory

setwd(plots.dir)

jpeg("Residuals_small_no_filter.jpg", width = 300, height = 300)
par(mfrow=c(3,3))

hist(residuals(Sa.gamm))

plot(small$response~ fitted(Sa.gamm))
abline(h=0, lty="dotted")
lines(lowess(fitted(Sa.gamm), small$response), col="red")


plot(residuals(Sa.gamm)~ fitted(Sa.gamm))
abline(h=0, lty="dotted")
lines(lowess(fitted(Sa.gamm), residuals(Sa.gamm)), col="red")

plot(small$t_sd,residuals(Sa.gamm))
abline(h=0, lty="dotted")
lines(lowess(small$t_sd, residuals(Sa.gamm)), col="red")

plot(small$ausbath,residuals(Sa.gamm))
abline(h=0, lty="dotted")
lines(lowess(small$ausbath, residuals(Sa.gamm)), col="red")

plot(small$dst2ramp,residuals(Sa.gamm))
abline(h=0, lty="dotted")
lines(lowess(small$dst2ramp, residuals(Sa.gamm)), col="red")

plot(small$status,residuals(Sa.gamm))

dev.off()

## Define predictive ability using 5-fold cross-validation

####### Large MaxN models

## Randomly shuffle the data

large<-large[sample(nrow(large)),]

#Create 5 equally size folds

folds <- cut(seq(1,nrow(large)),breaks=5,labels=FALSE)
#Perform 10 fold cross validation
for(i in 1:5){
  #Segement your data by fold using the which() function 
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- large[testIndexes, ]
  trainData <- large[-testIndexes, ]
  #Use the test and train data partitions however you desire...
}

### Fit the GAM model to the training data set

Ra.gamm=gam(response ~ s(relief, k = 3, bs = "cr") + te(gravity.200, ausbath, k = 3, bs = "cr") + status + te(longitude,latitude, k = 60, bs = "cr"),
            family=tw(),offset = sampling.deployment.duration.min. ,data=large)
plot(Ra.gamm, all.terms=TRUE,pages=1)
summary(Ra.gamm)

# Calculate Root Mean Squared Error and Normalized Root Mean Squared Error

predicts<-predict(Ra.gamm,testData,type="response")

rmse <- sqrt(mean((testData$response-predicts)^2))
rmse

nrmse <- sqrt(mean((testData$response-predicts)^2))/( max(testData$response)-min(testData$response) )
nrmse
nrmse*100 ## In percentage

############ Medium MaxN models

## Randomly shuffle the data

medium<-medium[sample(nrow(medium)),]

#Create 5 equally size folds

folds <- cut(seq(1,nrow(medium)),breaks=5,labels=FALSE)

#Perform 10 fold cross validation

for(i in 1:5){
  #Segement your data by fold using the which() function 
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- medium[testIndexes, ]
  trainData <- medium[-testIndexes, ]
  #Use the test and train data partitions however you desire...
}

### Fit the GAM model to the training data set

Ma.gamm=gam(response ~ s(ausbath, k = 3, bs = "cr") + s(t_m, by=status, k = 3, bs = "cr") + s(t_sd,by=status, k = 3, bs = "cr") + status + s(year, bs = "re") + s(day,k=3,bs="cr"),
            family=tw(),data=trainData)

# Calculate Root Mean Squared Error and Normalized Root Mean Squared Error

predicts<-predict(Ma.gamm,testData,type="response")

rmse <- sqrt(mean((testData$response-predicts)^2))
rmse

nrmse <- sqrt(mean((testData$response-predicts)^2))/( max(testData$response)-min(testData$response) )
nrmse
nrmse*100 ## In percentage

############ Small MaxN models

## Randomly shuffle the data

small<-small[sample(nrow(small)),]

#Create 5 equally size folds

folds <- cut(seq(1,nrow(small)),breaks=5,labels=FALSE)

#Perform 5 fold cross validation

for(i in 1:5){
  #Segement your data by fold using the which() function 
  testIndexes <- which(folds==i,arr.ind=TRUE)
  testData <- small[testIndexes, ]
  trainData <- small[-testIndexes, ]
  #Use the test and train data partitions however you desire...
}

### Fit the GAM model to the training data set

Sa.gamm=gam(response ~ s(relief, k = 3, bs = "cr") + s(t_m, k = 3, bs = "cr") + s(t_sd, k = 3, bs = "cr") + s(dst2shelf, k = 3, bs = "cr") + s(year, bs = "re") + s(day,k=3,bs="cr"),
            family=tw(),data=trainData)

# Calculate Root Mean Squared Error and Normalized Root Mean Squared Error

predicts<-predict(Sa.gamm,testData,type="response")

rmse <- sqrt(mean((testData$response-predicts)^2))
rmse

nrmse <- sqrt(mean((testData$response-predicts)^2))/( max(testData$response)-min(testData$response) )
nrmse
nrmse*100 ## In percentage

plot(La.gamm,all.terms = T,pages = 1)
summary(La.gamm)
plot(Ma.gamm,all.terms = T)

######### Spatial Plot of residuals ###########################################

### Large

res.large<-residuals(La.gamm)
large$residuals<-res.large
glimpse(large)

hist(large$residuals)

## Find breaks in the residuals

library(classInt)
library(viridis)

## Generate Fisher natural breaks

break.points<-classIntervals(large$residuals,5,"fisher")$brks
break.points

## Get Australia coastline

library(sp)
library(rgdal)
library(rgeos)
library(maptools)
library(broom)

setwd("G:/My Drive/Map layers")
dir()

australia<-readOGR('.','Australiaboundary67')
proj4string(australia)
australia <- spTransform(australia, CRS("+proj=longlat +datum=WGS84"))
australia<-fortify(australia)

large.residuals<-ggplot()+
  geom_point(data=large,aes(x=longitude,y=latitude,colour=residuals),size=0.01,alpha=0.8)+
  geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="lightgrey")+
  scale_colour_gradientn(colours=rainbow(6))+
  #scale_color_viridis(option = "magma",direction=-1,breaks=break.points)+
  coord_equal()+
  xlim(110,160)+
  ylim(-45,-9)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA,size = 1))

large.residuals


##### Convert into raster

setwd("C:/Users/22373243/Dropbox/Nestor_ascii")
dir()
relief<-raster("relief.asc")

### Create spatial Points dataframe

coordinates(large)<-~longitude+latitude
class(large)
proj4string(large)<-CRS("+init=epsg:4326")
large@proj4string

## Rasterize it

raster.large<-rasterize(large,relief,field="residuals")

## Plot it

raster.large_df<-as.data.frame(raster.large,xy=TRUE,na.rm=TRUE)
glimpse(raster.large_df)
colnames(raster.large_df) <- c("Longitude","Latitude","residuals")
glimpse(raster.large_df)

ggplot()+
  geom_tile(data=raster.large_df,aes(x=Longitude,y=Latitude,fill=residuals),alpha=0.8)+
  geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="lightgrey")+
  scale_fill_gradientn(colours=rainbow(6))+
  #scale_fill_viridis(option = "magma",direction = -1)+
  coord_equal()+
  xlim(110,160)+
  ylim(-45,-9)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA,size = 1))



### Medium

res.medium<-residuals(Ma.gamm)
medium$residuals<-res.medium
glimpse(medium)

hist(medium$residuals)

## Find breaks in the residuals

library(classInt)
library(viridis)

## Generate Fisher natural breaks

break.points<-classIntervals(medium$residuals,5,"fisher")$brks
break.points

## Get Australia coastline

library(sp)
library(rgdal)
library(rgeos)
library(maptools)
library(broom)

setwd("G:/My Drive/Map layers")
dir()

australia<-readOGR('.','Australiaboundary67')
proj4string(australia)
australia <- spTransform(australia, CRS("+proj=longlat +datum=WGS84"))
australia<-fortify(australia)

medium.residuals<-ggplot()+
  geom_point(data=medium,aes(x=longitude,y=latitude,colour=residuals),alpha=0.8)+
  geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="lightgrey")+
  scale_colour_gradientn(colours=rainbow(6))+
  #scale_color_viridis(option = "magma",direction=-1,breaks=break.points)+
  coord_equal()+
  xlim(110,160)+
  ylim(-45,-9)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA,size = 1))

medium.residuals

### Small

res.small<-residuals(Sa.gamm)
small$residuals<-res.small
glimpse(small)

hist(small$residuals)

## Find breaks in the residuals

library(classInt)
library(viridis)

## Generate Fisher natural breaks

break.points<-classIntervals(small$residuals,5,"fisher")$brks
break.points

## Get Australia coastline

library(sp)
library(rgdal)
library(rgeos)
library(maptools)
library(broom)

setwd("G:/My Drive/Map layers")
dir()

australia<-readOGR('.','Australiaboundary67')
proj4string(australia)
australia <- spTransform(australia, CRS("+proj=longlat +datum=WGS84"))
australia<-fortify(australia)

small.residuals<-ggplot()+
  geom_point(data=small,aes(x=longitude,y=latitude,colour=residuals),alpha=0.8)+
  geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="lightgrey")+
  scale_colour_gradientn(colours=rainbow(6))+
  #scale_color_viridis(option = "magma",direction=-1,breaks=break.points)+
  coord_equal()+
  xlim(110,160)+
  ylim(-45,-9)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA,size = 1))

small.residuals





