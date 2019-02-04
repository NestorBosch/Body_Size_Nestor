# A simple function for full subsets multiple regression in ecology with R
# 
# R. Fisher
# S.K. Wilson
# S.M. Sin
# A.C. Lee
# T.J. Langlois

# Script information----
# Plotting the most parsimonious GAMM mixed models - where there are interactions between continuous and factor predictor variables --

# librarys----
detach("package:plyr", unload=TRUE)#will error - no worries
detach("package:Rcpp", unload=TRUE)#will error - no worries
library(tidyr)
library(dplyr)
options(dplyr.width = Inf) #enables head() to display all coloums
library(ggplot2)
library(mgcv)
library(gridExtra)
library(grid)
library(RCurl) #needed to download data from GitHub
library(forcats)
# library(ggpubr)
# library(cowplot)
# library(readxl)
citation("ggpubr")
citation("cowplot")

rm(list=ls())

study<-"body.size"

# Set work directory----

work.dir=("C:/Users/22373243/Dropbox/Projects/Analysis/Analysis_GlobalArchive_Bodysize_Nestor") ### Desktop
setwd("~/workspace") ### Ecocloud platform

# Set sub directories----

data.dir=paste(work.dir,"Data",sep="/")
plots.dir=paste(work.dir,"Plots",sep="/")
model.out=paste(work.dir,"ModelOut",sep="/")
functions.dir=paste(work.dir,"Functions",sep="/")


# Plotting defaults----

# Theme-
Theme1 <-
  theme( # use theme_get() to see available options
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    # legend.background = element_rect(fill="white"),
    legend.background = element_blank(),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.text = element_text(size=15),
    legend.title = element_blank(),
    legend.position = c(0.2, 0.8),
    text=element_text(size=15),
    strip.text.y = element_text(size = 15,angle = 0),
    axis.title.x=element_text(vjust=0.3, size=15),
    axis.title.y=element_text(vjust=0.6, angle=90, size=15),
    axis.text.x=element_text(size=15),
    axis.text.y=element_text(size=15),
    axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
    axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
    strip.background = element_blank())

Theme.blank <-  # From Tims lobster plotting script 
  # For a blank plot to fill gaps in multyplot layout
  theme( # use theme_get() to see available options
    
    legend.background = element_rect(fill="transparent"),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    #text=element_text(size=15, family="Times New Roman"),
    legend.text = element_text(size=18),
    legend.title = element_text(size=18),
    #strip.text.x = element_text(size = 18,angle = 0),
    panel.grid.minor=element_blank(), panel.grid.major=element_blank(),
    #axis.title.x=element_blank(),
    axis.title.y=element_text(hjust= 0.5, angle=90, size=18),
    axis.text.x=element_blank(),
    axis.text.y=element_blank(),
    panel.background = element_blank(),
    #panel.border = element_rect(colour = "black", fill=NA),
    legend.position = c(0.2, 0.9),
    #plot.title=element_text(size=16, hjust=0.5, angle=0, face="bold"),
    plot.title=element_text(colour="black",face="bold", hjust=0.5, size=18),
    #plot.title = element_custom(),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_blank(),
    strip.background = element_blank())



# # functions for summarising data on plots----
se <- function(x) sd(x) / sqrt(length(x))
se.min <- function(x) (mean(x)) - se(x)
se.max <- function(x) (mean(x)) + se(x)


# Bring in and format the data----

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

## Create natural log of soak time

dat$sampling.deployment.duration.min.<-log(dat$sampling.deployment.duration.min.)
plot(dat$sampling.deployment.duration.min.)

## Transform predictors

dat$gravity.200<-log1p(dat$gravity.200)
dat$pop.den.200<-sqrt(dat$pop.den.200)
dat$dst2ramp<-log1p(dat$dst2ramp)
dat$dst2coast<-log1p(dat$dst2coast)
dat$dst2water<-log1p(dat$dst2water)
dat$no3_m<-log1p(dat$no3_m)
dat$relief<-log1p(dat$relief)

## Create data for each size class

large<-dat%>%
  filter(Taxa=='Large')%>%
  glimpse()

medium<-dat%>%
  filter(Taxa=='Medium')%>%
  glimpse()

small<-dat%>%
  filter(Taxa=='Small')%>%
  glimpse()


# Manually make the most parsimonious GAM models for each taxa ---- Large, Medium and Small
# pred.vars=c()

# MODEL # Large # relief + status +gravity.by.status + ausbath.t.status

### Plot of top models

La.gamm=gam(response~s(relief,k=3,bs="cr") + s(gravity.200,by=status,k=3,bs="cr")+ s(ausbath,by=status,k=3,bs="cr") + status,
            family=tw(),offset = log(sampling.deployment.duration.min.) ,data=large)
plot(La.gamm, all.terms=TRUE,pages=1)
summary(La.gamm)
AIC(La.gamm)

#+ s(campaignid,bs="re")

Ra.gamm=gam(response ~ s(relief, k = 3, bs = "cr") + te(gravity.200, ausbath, k = 3, bs = "cr") + status + s(campaignid,bs="re"),
            family=tw(),offset = sampling.deployment.duration.min. ,data=large)
plot(Ra.gamm,all.terms = TRUE,pages=1)
summary(Ra.gamm)
AIC(Ra.gamm)

Na.gamm=gam(response ~ s(relief, k = 3, bs = "cr") + te(gravity.200, ausbath, k = 3, bs = "cr") + status,
            family=tw(),offset = sampling.deployment.duration.min. ,data=large)
plot(Na.gamm,all.terms = TRUE,pages=1)
summary(Na.gamm)
AIC(Na.gamm)


#### Test correlation

testdata.1 <-expand.grid(relief=seq(min(large$relief),max(large$relief),length.out = 20),
                       gravity.200=seq(min(large$gravity.200),max(large$gravity.200),length.out = 20),
                       ausbath=seq(min(large$ausbath),max(large$ausbath),length.out = 20),
                       campaignid=(Ra.gamm$model$campaignid),
                       #sampling.deployment.duration.min.=1,
                       #year=(La.gamm$model$year),
                       status=c("Fished","No-take"))%>%
  glimpse()

fits.1 <- predict(Ra.gamm, newdata=testdata.1, type='response', se.fit=T)

testdata.2 <-expand.grid(relief=seq(min(large$relief),max(large$relief),length.out = 20),
                         gravity.200=seq(min(large$gravity.200),max(large$gravity.200),length.out = 20),
                         ausbath=seq(min(large$ausbath),max(large$ausbath),length.out = 20),
                         #campaignid=(Ra.gamm$model$campaignid),
                         #sampling.deployment.duration.min.=1,
                         #year=(La.gamm$model$year),
                         status=c("Fished","No-take"))%>%
  glimpse()

fits.2 <- predict(Na.gamm, newdata=testdata.2, type='response', se.fit=T)

plot(fit.values.1 ~ fit.values.2)
cor(fit.values.1,fit.values.2)

glimpse(fit.values.1)

fit.values.1<-fits.1$fit
fit.values.2<-fits.2$fit

fit.values.1<-head(fit.values.1,16000)

rm(La.predicts.relief,La.predicts.ausbath,La.predicts.status)

#pages = s(relief, k = 3, bs = "cr") + te(gravity.200, ausbath, k = 3, bs = "cr")

# predict - relief from MODEL ----

testdata <-expand.grid(relief=seq(min(large$relief),max(large$relief),length.out = 20),
                        gravity.200=mean(La.gamm$model$gravity.200),
                        ausbath=mean(La.gamm$model$ausbath),
                        #campaignid=(La.gamm$model$campaignid),
                        #sampling.deployment.duration.min.=1,
                        #year=(La.gamm$model$year),
                        status=c("Fished","No-take"))%>%
  glimpse()

fits <- predict(La.gamm, newdata=testdata, type='response', se.fit=T)

#head(fits,2)
La.predicts.relief = testdata%>%data.frame(fits)%>%
  group_by(relief)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()%>%
  glimpse()
write.csv(La.predicts.relief,"predicts.csv") #there is some BUG in dplyr - that this fixes
La.predicts.relief<-read.csv("predicts.csv")%>%
  glimpse()

# predict - gravity by status from MODEL ----

testdata <-expand.grid(gravity.200=seq(min(large$gravity.200),max(large$gravity.200),length.out = 20),
                       relief=mean(La.gamm$model$relief),
                       ausbath=mean(La.gamm$model$ausbath),
                       #campaignid=(La.gamm$model$campaignid),
                       #sampling.deployment.duration.min.=1,
                       status=c("Fished","No-take"))%>%
  glimpse()

fits <- predict(La.gamm, newdata=testdata, type='response', se.fit=T)

#head(fits,2)
La.predicts.gravity = testdata%>%data.frame(fits)%>%
  group_by(gravity.200,status)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()%>%
  glimpse()
write.csv(La.predicts.gravity,"predicts.csv") #there is some BUG in dplyr - that this fixes
La.predicts.gravity<-read.csv("predicts.csv")%>%
  glimpse()


# predict - ausbath.by.status from MODEL ----

testdata <-expand.grid(ausbath=seq(min(large$ausbath),max(large$ausbath),length.out = 20),
                       relief=mean(La.gamm$model$relief),
                       gravity.200=mean(La.gamm$model$gravity.200),
                       #campaignid=(La.gamm$model$campaignid),
                       #sampling.deployment.duration.min.=1,
                       status=c("Fished","No-take"))%>%
  glimpse()

fits <- predict(La.gamm, newdata=testdata, type='response', se.fit=T)

#head(fits,2)
La.predicts.ausbath = testdata%>%data.frame(fits)%>%
  group_by(ausbath,status)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()%>%
  glimpse()
write.csv(La.predicts.ausbath,"predicts.csv") #there is some BUG in dplyr - that this fixes
La.predicts.ausbath<-read.csv("predicts.csv")%>%
  glimpse()

## Predict - status - from model

testdata <- expand.grid(relief=mean(La.gamm$model$relief),
                        ausbath=mean(La.gamm$model$ausbath),
                        gravity.200=mean(La.gamm$model$gravity.200),
                        #campaignid=(La.gamm$model$campaignid),
                        #sampling.deployment.duration.min.=1,
                        status=c("Fished","No-take"))%>%
  glimpse()

fits <- predict(La.gamm, newdata=testdata, type='response', se.fit=T)

#head(fits,2)
La.predicts.status = testdata%>%data.frame(fits)%>%
  group_by(status)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()%>%
  glimpse()
write.csv(La.predicts.status,"predicts.csv") #there is some BUG in dplyr - that this fixes
La.predicts.status<-read.csv("predicts.csv")%>%
  glimpse()

## Plot large

ggmod.La.relief<-ggplot()+
  ylab("")+
  xlab("Relief")+
  #ylim(0,4)+
  #xlim(0,5)+
  #geom_point(data=large,aes(x=relief,y=response),alpha=0.2,size=2,show.legend = TRUE)+
  geom_line(data=La.predicts.relief,aes(x=relief,y=response),alpha=0.9)+
  geom_line(data=La.predicts.relief,aes(x=relief,y=response - se.fit),linetype="dashed",alpha=0.9)+
  geom_line(data=La.predicts.relief,aes(x=relief,y=response + se.fit),linetype="dashed",alpha=0.9)+
  theme_classic()+
  Theme1+
  annotate("text", x =-Inf, y=Inf, label = "(A) Large (> 30 cm)",vjust = 1, hjust = -.1,size=5)

ggmod.La.relief

ggmod.La.gravity.by.status<-ggplot(aes(x=gravity.200,y=response,colour=status),data=La.predicts.gravity) +
  scale_color_manual(labels = c("Fished", "No-take"),values=c("brown1", "blue3"))+
  ylab("")+
  xlab('Gravity')+
  #ylim(0,4)+
  #geom_point(data=large,aes(x=gravity.200,y=response),alpha=0.2,size=2,show.legend = TRUE)+
  geom_line(data=La.predicts.gravity,aes(x=gravity.200,y=response,colour=status),alpha=0.9,show.legend = F)+
  geom_line(data=La.predicts.gravity,aes(x=gravity.200,y=response - se.fit),linetype="dashed",alpha=0.9,show.legend = F)+
  geom_line(data=La.predicts.gravity,aes(x=gravity.200,y=response + se.fit),linetype="dashed",alpha=0.9,show.legend = F)+
  theme_classic()+
  Theme1+
  theme(legend.position="none")+
  annotate("text", x = -Inf, y=Inf, label = "(C)",vjust = 1, hjust = -.1,size=5)

ggmod.La.gravity.by.status

ggmod.ausbath.by.Status<-ggplot(aes(x=ausbath,y=response,colour=status),data=La.predicts.ausbath) +
  scale_color_manual(labels = c("Fished", "No-take"),values=c("brown1", "blue3"))+
  ylab("")+
  xlab('Depth (m)')+
  #ylim(0,4)+
  #geom_point(alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=La.predicts.ausbath,aes(x=ausbath,y=response,colour=status),alpha=0.9,show.legend = F)+
  geom_line(data=La.predicts.ausbath,aes(x=ausbath,y=response - se.fit),linetype="dashed",alpha=0.9,show.legend = F)+
  geom_line(data=La.predicts.ausbath,aes(x=ausbath,y=response + se.fit),linetype="dashed",alpha=0.9,show.legend = F)+
  theme_classic()+
  Theme1+
  theme(legend.position="none")+
  annotate("text", x = -Inf, y=Inf, label = "(C)",vjust = 1, hjust = -.1,size=5)

ggmod.ausbath.by.Status

ggmod.La.status<- ggplot(aes(x=status,y=response,fill=status), data=La.predicts.status) +
  ylab("")+
  xlab('Status')+
  #ylim(0,4)+
  scale_colour_manual(labels = c("Fished", "No-Take"),values=c("brown1", "blue3"))+
  scale_fill_manual(labels = c("Fished", "No-Take"),values=c("brown1", "blue3"))+
  geom_bar(stat = "identity",show.legend=TRUE)+
  geom_errorbar(aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5,show.legend=FALSE) +
  theme_classic()+
  Theme1+
  #annotate("Text",x="No-Take",y=3.5,label="***",size=10)+
  annotate("text", x = -Inf, y=Inf, label = "(D)",vjust = 1, hjust = -.1,size=5)+
  annotate("text", x = -Inf, y=Inf, label = "",vjust = 1, hjust = -.1,size=5,fontface="italic")

ggmod.La.status

## Plot together in 3 by 4 plot

setwd(plots.dir)

# MAke blank plot to fill gaps in multyplot layout

blank <- ggplot(aes(x=relief,y=response,colour=status), data=large) +
  ylab(NULL)+
  xlab(NULL)+
  Theme.blank
blank

# To see what they will look like use grid.arrange()
## Large

grid.arrange(ggmod.La.relief,ggmod.La.gravity.by.status,ggmod.ausbath.by.Status,ggmod.La.status,nrow=2)

# MODEL # Medium - t_m + t_sd + dst2coast + status

Ma.gamm=gam(response~s(t_m, k = 3, bs = "cr") + s(t_sd, k = 3, bs = "cr") + s(ausbath,by=status,k=3,bs="cr") + status + s(year,bs="re")+ s(day,k=3,bs="cr"),
            family=tw(),data=medium)
plot(Ma.gamm, all.terms=TRUE,pages=1)
summary(Ma.gamm)


#+ s(year, bs = "re") + s(day, k = 3, bs = "cr")

# predict - t_m from MODEL ----

testdata <- expand.grid(t_m=seq(min(medium$t_m),max(medium$t_m),length.out = 80),
                        t_sd=mean(Ma.gamm$model$t_sd),
                        ausbath=mean(Ma.gamm$model$ausbath),
                        day=mean(Ma.gamm$model$day),
                        year=(Ma.gamm$model$year),
                        status=c("Fished","No-take"))%>%
  glimpse()

fits <- predict(Ma.gamm, newdata=testdata, type='response', se.fit=T)

#head(fits,2)
Ma.predicts.t_m = testdata%>%data.frame(fits)%>%
  group_by(t_m)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()%>%
  glimpse()
write.csv(Ma.predicts.t_m,"predicts.csv") #there is some BUG in dplyr - that this fixes
Ma.predicts.t_m<-read.csv("predicts.csv")%>%
  glimpse()

# predict - t_sd from MODEL ----

testdata <- expand.grid(t_sd=seq(min(medium$t_sd),max(medium$t_sd),length.out = 80),
                        t_m=mean(Ma.gamm$model$t_m),
                        ausbath=mean(Ma.gamm$model$ausbath),
                        day=mean(Ma.gamm$model$day),
                        year=(Ma.gamm$model$year),
                        status=c("Fished","No-take"))%>%
  glimpse()

fits <- predict(Ma.gamm, newdata=testdata, type='response', se.fit=T)

#head(fits,2)
Ma.predicts.t_sd = testdata%>%data.frame(fits)%>%
  group_by(t_sd)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()%>%
  glimpse()
write.csv(Ma.predicts.t_sd,"predicts.csv") #there is some BUG in dplyr - that this fixes
Ma.predicts.t_sd<-read.csv("predicts.csv")%>%
  glimpse()

## Predict - depth by status - from model ---

testdata <- expand.grid(ausbath=seq(min(medium$ausbath),max(medium$ausbath),length.out = 80),
                        t_m=mean(Ma.gamm$model$t_m),
                        t_sd=mean(Ma.gamm$model$t_sd),
                        day=mean(Ma.gamm$model$day),
                        year=(Ma.gamm$model$year),
                        status=c("Fished","No-take"))%>%
  glimpse()

fits <- predict(Ma.gamm, newdata=testdata, type='response', se.fit=T)

#head(fits,2)
Ma.predicts.depth.by.status = testdata%>%data.frame(fits)%>%
  group_by(ausbath,status)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()%>%
  glimpse()
write.csv(Ma.predicts.depth.by.status,"predicts.csv") #there is some BUG in dplyr - that this fixes
Ma.predicts.depth.by.status<-read.csv("predicts.csv")%>%
  glimpse()

## Predict - status - from model

testdata <- expand.grid(t_m=mean(Ma.gamm$model$t_m),
                        t_sd=mean(Ma.gamm$model$t_sd),
                        ausbath=mean(Ma.gamm$model$ausbath),
                        day=mean(Ma.gamm$model$day),
                        year=(Ma.gamm$model$year),
                        status=c("Fished","No-take"))%>%
  glimpse()

fits <- predict(Ma.gamm, newdata=testdata, type='response', se.fit=T)

#head(fits,2)
Ma.predicts.status = testdata%>%data.frame(fits)%>%
  group_by(status)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()%>%
  glimpse()
write.csv(Ma.predicts.status,"predicts.csv") #there is some BUG in dplyr - that this fixes
Ma.predicts.status<-read.csv("predicts.csv")%>%
  glimpse()

## Plot Medium

ggmod.Ma.t_m<-ggplot()+
  ylab("Relative abundance (MaxN)")+
  xlab("Temperature (mean)")+
  ylim(0,11)+
  #xlim(0,41)+
  #geom_point(data=large,aes(x=relief,y=response),alpha=0.7,size=2,show.legend = TRUE)+
  geom_line(data=Ma.predicts.t_m,aes(x=t_m,y=response),alpha=0.9)+
  geom_line(data=Ma.predicts.t_m,aes(x=t_m,y=response - se.fit),linetype="dashed",alpha=0.9)+
  geom_line(data=Ma.predicts.t_m,aes(x=t_m,y=response + se.fit),linetype="dashed",alpha=0.9)+
  theme_classic()+
  Theme1+
  annotate("text", x =-Inf, y=Inf, label ="(E) Medium (20 - 30 cm)",vjust = 1, hjust = -.1,size=5)

ggmod.Ma.t_m

ggmod.Ma.t_sd<-ggplot()+
  ylab("")+
  xlab("Temperature (sd)")+
  ylim(0,11)+
  #xlim(0,41)+
  #geom_point(data=large,aes(x=t_sd,y=response),alpha=0.7,size=2,show.legend = TRUE)+
  geom_line(data=Ma.predicts.t_sd,aes(x=t_sd,y=response),alpha=0.9)+
  geom_line(data=Ma.predicts.t_sd,aes(x=t_sd,y=response - se.fit),linetype="dashed",alpha=0.9)+
  geom_line(data=Ma.predicts.t_sd,aes(x=t_sd,y=response + se.fit),linetype="dashed",alpha=0.9)+
  theme_classic()+
  Theme1+
  annotate("text", x =-Inf, y=Inf, label = "(F)",vjust = 1, hjust = -.1,size=5)

ggmod.Ma.t_sd

ggmod.Ma.depth.by.status<-ggplot(aes(x=ausbath,y=response,colour=status),data=Ma.predicts.depth.by.status) +
  scale_color_manual(labels = c("Fished", "No-take"),values=c("brown1", "blue3"))+
  ylab("")+
  xlab('Depth (m)')+
  ylim(0,5.5)+
  #geom_point(alpha=0.75, size=2,show.legend=FALSE)+
  geom_line(data=Ma.predicts.depth.by.status,aes(x=ausbath,y=response,colour=status),alpha=0.9,show.legend = F)+
  geom_line(data=Ma.predicts.depth.by.status,aes(x=ausbath,y=response - se.fit),linetype="dashed",alpha=0.9,show.legend = F)+
  geom_line(data=Ma.predicts.depth.by.status,aes(x=ausbath,y=response + se.fit),linetype="dashed",alpha=0.9,show.legend = F)+
  theme_classic()+
  Theme1+
  theme(legend.position="none")+
  annotate("text", x = -Inf, y=Inf, label = "(C)",vjust = 1, hjust = -.1,size=5)

ggmod.Ma.depth.by.status

ggmod.Ma.status<- ggplot(aes(x=status,y=response,fill=status), data=Ma.predicts.status) +
  ylab("")+
  xlab('Status')+
  ylim(0,5)+
  scale_colour_manual(labels = c("Fished", "No-Take"),values=c("brown1", "blue3"))+
  scale_fill_manual(labels = c("Fished", "No-Take"),values=c("brown1", "blue3"))+
  geom_bar(stat = "identity",show.legend=FALSE)+
  geom_errorbar(aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5,show.legend=FALSE) +
  theme_classic()+
  Theme1+
  #annotate("Text",x="No-Take",y=3.5,label="***",size=10)+
  annotate("text", x = -Inf, y=Inf, label = "(H)",vjust = 1, hjust = -.1,size=5)+
  annotate("text", x = -Inf, y=Inf, label = "",vjust = 1, hjust = -.1,size=5,fontface="italic")

ggmod.Ma.status

## Plot together in 3 by 4 plot

setwd(plots.dir)

# MAke blank plot to fill gaps in multyplot layout

blank <- ggplot(aes(x=relief,y=response,colour=status), data=large) +
  ylab(NULL)+
  xlab(NULL)+
  Theme.blank
blank

# To see what they will look like use grid.arrange()
## Medium

grid.arrange(ggmod.Ma.t_m,ggmod.Ma.t_sd,ggmod.Ma.depth.by.status,ggmod.Ma.status,nrow=1)

# MODEL # Small - relief + t_m + t_sd + dst2shelf

Sa.gamm=gam(response~s(relief, k = 3, bs = "cr") + s(t_m, k = 3, bs = "cr") + s(t_sd, k = 3, bs = "cr") + s(year,bs="re")+ s(day,k=3,bs="cr"),
            family=tw(),data=small)
plot(Sa.gamm, all.terms=TRUE,pages=1)
summary(Sa.gamm)

# predict - relief from MODEL ----

testdata <- expand.grid(relief=seq(min(small$relief),max(small$relief),length.out = 80),
                        t_m=mean(Sa.gamm$model$t_m),
                        t_sd=mean(Sa.gamm$model$t_sd),
                        #aspect=mean(Sa.gamm$model$aspect),
                        day=mean(Sa.gamm$model$day),
                        year=(Sa.gamm$model$year))%>%
  glimpse()

fits <- predict(Sa.gamm, newdata=testdata, type='response', se.fit=T)

#head(fits,2)
Sa.predicts.relief = testdata%>%data.frame(fits)%>%
  group_by(relief)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()%>%
  glimpse()
write.csv(Sa.predicts.relief,"predicts.csv") #there is some BUG in dplyr - that this fixes
Sa.predicts.relief<-read.csv("predicts.csv")%>%
  glimpse()

# predict - t_m from MODEL ----

testdata <- expand.grid(t_m=seq(min(small$t_m),max(small$t_m),length.out = 80),
                        relief=mean(Sa.gamm$model$relief),
                        t_sd=mean(Sa.gamm$model$t_sd),
                        #aspect=mean(Sa.gamm$model$aspect),
                        day=mean(Sa.gamm$model$day),
                        year=(Sa.gamm$model$year))%>%
  glimpse()

fits <- predict(Sa.gamm, newdata=testdata, type='response', se.fit=T)

#head(fits,2)
Sa.predicts.t_m = testdata%>%data.frame(fits)%>%
  group_by(t_m)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()%>%
  glimpse()
write.csv(Sa.predicts.t_m,"predicts.csv") #there is some BUG in dplyr - that this fixes
Sa.predicts.t_m<-read.csv("predicts.csv")%>%
  glimpse()

# predict - t_sd from MODEL ----

testdata <- expand.grid(t_sd=seq(min(small$t_sd),max(small$t_sd),length.out = 80),
                        relief=mean(Sa.gamm$model$relief),
                        t_m=mean(Sa.gamm$model$t_m),
                        #aspect=mean(Sa.gamm$model$aspect),
                        day=mean(Sa.gamm$model$day),
                        year=(Sa.gamm$model$year))%>%
  glimpse()

fits <- predict(Sa.gamm, newdata=testdata, type='response', se.fit=T)

#head(fits,2)
Sa.predicts.t_sd = testdata%>%data.frame(fits)%>%
  group_by(t_sd)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()%>%
  glimpse()
write.csv(Sa.predicts.t_sd,"predicts.csv") #there is some BUG in dplyr - that this fixes
Sa.predicts.t_sd<-read.csv("predicts.csv")%>%
  glimpse()

# predict - aspect from MODEL ----

testdata <- expand.grid(aspect=seq(min(small$aspect),max(small$aspect),length.out = 80),
                        relief=mean(Sa.gamm$model$relief),
                        t_m=mean(Sa.gamm$model$t_m),
                        t_sd=mean(Sa.gamm$model$t_sd),
                        aspect=mean(Sa.gamm$model$aspect),
                        day=mean(Sa.gamm$model$day),
                        year=(Sa.gamm$model$year))%>%
  glimpse()

fits <- predict(Sa.gamm, newdata=testdata, type='response', se.fit=T)

#head(fits,2)
Sa.predicts.aspect = testdata%>%data.frame(fits)%>%
  group_by(aspect)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()%>%
  glimpse()
write.csv(Sa.predicts.aspect,"predicts.csv") #there is some BUG in dplyr - that this fixes
Sa.predicts.aspect<-read.csv("predicts.csv")%>%
  glimpse()

## Plot Small

ggmod.Sa.relief<-ggplot()+
  ylab("")+
  xlab("Relief")+
  ylim(0,3.6)+
  #xlim(0,41)+
  #geom_point(data=large,aes(x=relief,y=response),alpha=0.7,size=2,show.legend = TRUE)+
  geom_line(data=Sa.predicts.relief,aes(x=relief,y=response),alpha=0.9)+
  geom_line(data=Sa.predicts.relief,aes(x=relief,y=response - se.fit),linetype="dashed",alpha=0.9)+
  geom_line(data=Sa.predicts.relief,aes(x=relief,y=response + se.fit),linetype="dashed",alpha=0.9)+
  theme_classic()+
  Theme1+
  annotate("text", x =-Inf, y=Inf, label = "(I) Small (< 20 cm)",vjust = 1, hjust = -.1,size=5)

ggmod.Sa.relief

ggmod.Sa.t_m<-ggplot()+
  ylab("")+
  xlab("Temperature (mean)")+
  ylim(0,3.6)+
  #xlim(0,41)+
  #geom_point(data=large,aes(x=t_sd,y=response),alpha=0.7,size=2,show.legend = TRUE)+
  geom_line(data=Sa.predicts.t_m,aes(x=t_m,y=response),alpha=0.9)+
  geom_line(data=Sa.predicts.t_m,aes(x=t_m,y=response - se.fit),linetype="dashed",alpha=0.9)+
  geom_line(data=Sa.predicts.t_m,aes(x=t_m,y=response + se.fit),linetype="dashed",alpha=0.9)+
  theme_classic()+
  Theme1+
  annotate("text", x =-Inf, y=Inf, label = "(J)",vjust = 1, hjust = -.1,size=5)

ggmod.Sa.t_m

ggmod.Sa.t_sd<-ggplot()+
  ylab("")+
  xlab("Temperature (sd)")+
  ylim(0,3.6)+
  #xlim(0,41)+
  #geom_point(data=large,aes(x=t_sd,y=response),alpha=0.7,size=2,show.legend = TRUE)+
  geom_line(data=Sa.predicts.t_sd,aes(x=t_sd,y=response),alpha=0.9)+
  geom_line(data=Sa.predicts.t_sd,aes(x=t_sd,y=response - se.fit),linetype="dashed",alpha=0.9)+
  geom_line(data=Sa.predicts.t_sd,aes(x=t_sd,y=response + se.fit),linetype="dashed",alpha=0.9)+
  theme_classic()+
  Theme1+
  annotate("text", x =-Inf, y=Inf, label = "(K)",vjust = 1, hjust = -.1,size=5)

ggmod.Sa.t_sd

ggmod.Sa.aspect<-ggplot()+
  ylab("")+
  xlab("Seabed aspect")+
  #ylim(0,3)+
  #xlim(0,41)+
  #geom_point(data=large,aes(x=t_sd,y=response),alpha=0.7,size=2,show.legend = TRUE)+
  geom_line(data=Sa.predicts.aspect,aes(x=aspect,y=response),alpha=0.9)+
  geom_line(data=Sa.predicts.aspect,aes(x=aspect,y=response - se.fit),linetype="dashed",alpha=0.9)+
  geom_line(data=Sa.predicts.aspect,aes(x=aspect,y=response + se.fit),linetype="dashed",alpha=0.9)+
  theme_classic()+
  Theme1+
  annotate("text", x =-Inf, y=Inf, label = "(K)",vjust = 1, hjust = -.1,size=5)

ggmod.Sa.aspect



## Plot together in 3 by 4 plot

setwd(plots.dir)

# MAke blank plot to fill gaps in multyplot layout

blank <- ggplot(aes(x=relief,y=response,colour=status), data=large) +
  ylab(NULL)+
  xlab(NULL)+
  Theme.blank
blank

# To see what they will look like use grid.arrange()
## Large

grid.arrange(ggmod.Sa.relief,ggmod.Sa.t_m,ggmod.Sa.t_sd,nrow=1)


# Use arrangeGrob ONLY - as we can pass this to ggsave! Note use of raw ggplot's
setwd(plots.dir)

name<-"Body_size"
arrangeGrob(ggmod.La.relief,ggmod.La.depth.by.status,ggmod.La.dst2ramp.by.Status,ggmod.La.status,
            ggmod.Ma.t_m,ggmod.Ma.t_sd,ggmod.Ma.depth.by.status,ggmod.Ma.status,
            ggmod.Sa.relief,ggmod.Sa.t_m,ggmod.Sa.t_sd,nrow=3)%>%
  ggsave(file= paste("gamm",study,name,"png",sep="."), width = 27, height = 30,units = "cm")

