# A simple function for full subsets multiple regression in ecology with R
# 
# R. Fisher
# S.K. Wilson
# S.M. Sin
# A.C. Lee
# T.J. Langlois


# Script information----
# This script is designed to work with long format data - where response variables are stacked one upon each other (see http://tidyr.tidyverse.org/)
# There are random factors
# We have used a Tweedie error distribution to account for the high occurence of zero values in the dataset.
# We have implemented the ramdom effects and Tweedie error distribution using the mgcv() package

# Install packages

install.packages("Rcpp")
install.packages("car")
install.packages("doBy")
install.packages("gplots")
install.packages("doParallel")
install.packages("gamm4")
install.packages("doSNOW")


# librarys----

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
library(ggplot2)
library(devtools)

rm(list=ls())

# GAMM Size models----

name<-"Size"


# Set work directory----

work.dir=("C:/Users/22373243/Dropbox/Projects/Analysis/Analysis_GlobalArchive_Bodysize_Nestor") # Nestor
setwd("~/workspace") # Ecocloud

# Set sub directories----

data.dir=paste(work.dir,"Data",sep="/")
plots.dir=paste(work.dir,"Plots",sep="/")
model.out=paste(work.dir,"ModelOut",sep="/")
functions.dir=paste(work.dir,"Functions",sep="/")

# Install functions from package ----

devtools::install_github("beckyfisher/FSSgam_package")

library(FSSgam)

?full.subsets.gam

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

### If you want to check how the model performs run with a subset of the data - In this case I will use ca. 1/4 of the data = 2500 samples

dat<-dat[sample(nrow(dat), 2000), ]

## Code year and soak time as a factor - to include as a random effect
            
dat$year<-as.factor(dat$year)
dat$sampling.deployment.duration.min.<-as.numeric(dat$sampling.deployment.duration.min.)
plot(dat$sampling.deployment.duration.min.)

## Rescale biomass - to kg

## Transform soak time to the natural log - to use as an offset 

dat$sampling.deployment.duration.min.<-log(dat$sampling.deployment.duration.min.)
plot(dat$sampling.deployment.duration.min.)

#### Check distribution of the response

## Large

par(mfrow=c(1,1))

large<-dat%>%
  filter(Taxa=="Large")%>%
  glimpse()
hist(dat$response)

### Medium
  
medium<-dat%>%
  filter(Taxa=="Medium")%>%
  glimpse()
hist(medium$response)

### Small

small<-dat%>%
  filter(Taxa=="Small")%>%
  glimpse()
hist(small$response)

levels(dat$campaignid)
levels(dat$method)
levels(dat$status)

#### Check distribution of predictors dst2ramp and ausbath by status
#### We need to make sure that there is fairly even distribution of points across the range of the predictors 
#### For both to allow interactions to be valid

large.fished<-large%>%
  filter(status=="Fished")%>%
  glimpse()

large.notake<-large%>%
  filter(status=="No-take")%>%
  glimpse()

plot.new()
par(mfrow=c(2,2))

## Dst2ramp

hist(large.fished$dst2ramp)
hist(large.notake$dst2ramp)
plot(large.fished$dst2ramp)
plot(large.notake$dst2ramp)

## Depth

plot.new()
par(mfrow=c(2,2))

hist(large.fished$depth)
hist(large.notake$depth)
plot(large.fished$depth)
plot(large.notake$depth)

## Ausbath

plot.new()
par(mfrow=c(2,2))

hist(large.fished$ausbath)
hist(large.notake$ausbath)
plot(large.fished$ausbath)
plot(large.notake$ausbath)

# Set predictor (i.e. continous) variables----
names(dat)
pred.vars=c("sand","relief","t_m","t_sd","no3_m","ausbath",
            "dst2water","dst2coast","dst2ramp","pop.den.200","Eastness","Northness","gravity.200","longitude")

test.dat=dat

# Check for correlation of predictor variables- remove anything highly correlated (>0.90)---
round(cor(test.dat[,pred.vars]),2)  ## There is some error here when applied to Bio-oracle covariates

# nothing is highly correlated

# Plot of likely transformations
plot.new()
 
par(mfrow=c(3,2))
for (i in pred.vars) {
  x<-test.dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

# Review of individual predictors - we have to make sure they have an even distribution---
#If the data are squewed to low numbers try sqrt>log or if squewed to high numbers try ^2 of ^3
# Is data is clumped into distances - of a feature of interest - try as a linear predictor

## Remove potential outlier
## This is done when I bring in the data

## Transform predictors if necessary
## However - I don't like to transfor predictor as the transformation is tricky to interpret
## Yet - when plotting predicts you can predict in the response scale using type="response"

dat$gravity.200<-log1p(dat$gravity.200)
dat$pop.den.200<-sqrt(dat$pop.den.200)
dat$dst2ramp<-log1p(dat$dst2ramp)
dat$dst2coast<-log1p(dat$dst2coast)
dat$dst2water<-log1p(dat$dst2water)
dat$no3_m<-log1p(dat$no3_m)
dat$relief<-log1p(dat$relief)


# Run the full subset model selection----
# Set directory for the model outputs----

setwd(model.out)

# Check to make sure Response vector has not more than 80% zeros---

unique.vars=unique(as.character(dat$Taxa))
unique.vars.use=character()
for(i in 1:length(unique.vars)){
  temp.dat=dat[which(dat$Taxa==unique.vars[i]),]
  if(length(which(temp.dat$response==0))/nrow(temp.dat)<0.85){
    unique.vars.use=c(unique.vars.use,unique.vars[i])}
}
unique.vars.use
write.csv(unique.vars.use,file=paste(name,"unique.vars.use.csv",sep = "_"))

# Set variables needed for the FSS function-

names(dat)
resp.vars=unique.vars
use.dat=dat
factor.vars=c("status") # Status as a Factor with two levels
linear.vars=c("longitude") 
#cyclic.vars=c("aspect")
null.vars=c("campaignid","sampling.deployment.duration.min.")
out.all=list()
var.imp=list()

# Loop through the FSS function for each Taxa----
  
for(i in 1:length(resp.vars)){use.dat=dat[which(dat$Taxa==resp.vars[i]),]
  
  Model1=gam(response ~ s(relief,k=3,bs="cr") + s(campaignid,bs="re") + s(sampling.deployment.duration.min.,bs="re"),
             family=tw(), data=use.dat)
  
  out.list=full.subsets.gam(use.dat=use.dat,
                            test.fit=Model1,
                            pred.vars.cont=pred.vars,
                            pred.vars.fact=factor.vars,
                            linear.vars = linear.vars,
                            factor.smooth.interactions = list(
                              fact.vars=c("status"),
                              cont.vars=c("ausbath","dst2ramp","pop.den.200","gravity.200")),
                              #linear.vars=c("ausbath")),
                            #cyclic.vars = cyclic.vars,
                            smooth.smooth.interactions = c("gravity.200","ausbath"),
                            k=3,
                            max.predictors=4,
                            null.terms="s(campaignid,bs='re') + s(sampling.deployment.duration.min.,bs='re')",
                            max.models=10000,
                            #save.model.fits = F,
                            parallel=T)
  names(out.list)
  
  out.list$failed.models # examine the list of failed models
  mod.table=out.list$mod.data.out  # look at the model selection table
  mod.table=mod.table[order(mod.table$AICc),]
  mod.table$cumsum.wi=cumsum(mod.table$wi.AICc)
  out.i=mod.table[which(mod.table$delta.AICc<=2),]
  out.all=c(out.all,list(out.i))
  #var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Either raw importance score
  var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Or importance score weighted by r2
  
  # plot the best models
  for(m in 1:nrow(out.i)){
    best.model.name=as.character(out.i$modname[m])
    
    png(file=paste(name,m,resp.vars[i],"mod_fits.png",sep="_"))
    if(best.model.name!="null"){
      par(mfrow=c(3,1),mar=c(9,4,3,1))
      best.model=out.list$success.models[[best.model.name]]
      plot(best.model,all.terms=T,pages=1,residuals=T,pch=16)
      mtext(side=2,text=resp.vars[i],outer=F)}  
    dev.off()
  }
}

# Model fits and importance scores---
names(out.all)=resp.vars
names(var.imp)=resp.vars
all.mod.fits=do.call("rbind",out.all)
all.var.imp=do.call("rbind",var.imp)
write.csv(all.mod.fits[,-2],file=paste(name,"all.mod.fits.csv",sep="_"))
write.csv(all.var.imp,file=paste(name,"all.var.imp.csv",sep="_"))

