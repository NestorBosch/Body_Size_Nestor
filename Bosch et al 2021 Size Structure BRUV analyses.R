#############################################################################################
############## Bosch et al. 2021. Body size - Australia national s-BRUVs synthesis ########## 

### Index 
## (1) Data collation 
## (2) Generating a snapshot of the analyses dataframes 
## (3) GAMM models - Assemblage-level
## (4) Plot of top models - Model Validation - GAMMs - Assemblage-level
## (5) Importane scores - Assemblage-level models
## (6) GLMM models - Regional species groups
## (7) Model averaged coefficients - Regional species groups
## (8) Summed AICc weights for Anthropogenic, habitat, and environmental predictors
## (9) Extended Data Tables and Figures

## Clean the working directory

rm(list=ls())

## Load libraries 

library(tidyr)
library(dplyr)
library(forcats)
library(ggplot2)
library(stringr)
library(classInt)
library(rgdal)
library(viridis)
library(tibble)
library(spdep)
library(geosphere)
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
library(devtools)
library(FSSgam)
library(ggpubr)

### Set plotting defaults

Theme1 <-
  theme( # use theme_get() to see available options
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill="white"),
    legend.key = element_blank(), # switch off the rectangle around symbols in the legend
    legend.text = element_text(size=10),
    legend.title = element_text(size=14, face="bold"),
    legend.position = "right",
    legend.direction="vertical",
    text=element_text(size=14),
    strip.text.y = element_text(size = 14,angle = 0),
    axis.title.x=element_text(vjust=0.3, size=14),
    axis.title.y=element_text(vjust=0.6, angle=90, size=14),
    axis.text.x=element_text(size=12,angle = 0, hjust=1,vjust=0.5),
    axis.text.y=element_text(size=12),
    axis.line.x=element_line(colour="black", size=0.5,linetype='solid'),
    axis.line.y=element_line(colour="black", size=0.5,linetype='solid'),
    strip.background = element_blank())

# Study and directories----

name<-"AU_BRUVs_Fish_Size"

# Set work directory---

work.dir=setwd("C:/Users/22373243/Dropbox/Projects/Analysis/Analysis_GlobalArchive_Bodysize_Nestor")  ### Nestor

# Set sub directories----

data.dir=paste(work.dir,"Data",sep="/")
plots.dir=paste(work.dir,"Plots",sep="/")
model.out=paste(work.dir,"ModelOut",sep="/")
func.dir=paste(work.dir,"Functions",sep="/")

############################################################################################
############################ (1) Data collation ############################################

## Bring in Au Length Metadata ----

setwd(data.dir)
dir()

length.metadata<-read.csv("Au.Length.Metadata.2019-11-14.csv")%>%
  dplyr::select(-site,-management.and.monitoring.areas,-location,-marine.region,-observer,-successful,-successful.count,-successful.length,-custodian)%>%
  dplyr::mutate(id=paste(campaignid,sample,sep="."))%>%
  dplyr::filter(!campaignid%in%c('2016-09_Scott & Rowleys_FinPrint_stereoBRUVS'))%>%
  dplyr::filter(!grepl('Christmas.Island|Cocos.Keeling.Islands|LHI',id))%>% ## Filter regions outside the continental platform
  dplyr::filter(!grepl('QLD',state))%>% ## Filter deployments in Queensland where abundance is not reported as MaxN
  glimpse()

## Convert time into decimal hours 

length.metadata$time<-as.character(length.metadata$time)

hhmmss2dec <- function(x) {
  xlist <- strsplit(x,split=":")
  h <- as.numeric(sapply(xlist,"[",1))
  m <- as.numeric(sapply(xlist,"[",2))
  s <- as.numeric(sapply(xlist,"[",3))
  xdec <- h+(m/60)+(s/3600)
  return(xdec)
}

length.metadata$time<-hhmmss2dec(length.metadata$time)
length.metadata$time<-as.numeric(length.metadata$time)
hist(length.metadata$time)

##### Filter final lenghtmetadata - daylight 
##### Also - select only those variables of interest

names(length.metadata)
length.metadata<-length.metadata%>%
  filter(time>=9&time<=17)%>%
  glimpse()
hist(length.metadata$time)

### Make sure status is coded as a factor

length.metadata$status<-as.factor(length.metadata$status)
levels(length.metadata$status)
plot(length.metadata$status)

#### Check depth distribution of samples

hist(length.metadata$depth)

## Plot depth in space

ggplot(length.metadata,aes(x=longitude,y=latitude,size=depth,colour=depth))+
  geom_point()+
  scale_colour_viridis()+
  theme_classic()

## Explore depth across states

ggplot(length.metadata,aes(x=depth,fill=state))+
  geom_density()+
  scale_x_continuous(breaks=c(0,50,100,200,400))+
  facet_wrap(~state)+
  theme_classic()

## Explore depth across management status

ggplot(length.metadata,aes(x=depth,fill=status))+
  geom_density()+
  scale_x_continuous(breaks=c(0,50,100,200,400))+
  facet_wrap(~status)+
  theme_classic()

## Explore depth across states
## Everything comparable - although shallow areas are underrepresented in Tasmania and NSW

ggplot(length.metadata,aes(x=depth,fill=state))+
  geom_density()+
  scale_x_continuous(breaks=c(0,50,100,200,400))+
  facet_wrap(~state)+
  theme_classic()

## Sparsity of data beyond 50 m depth - we retain stereo-BRUVs deployments within the 50 m isobath

length.metadata<-length.metadata%>%
  filter(depth<=50)%>%
  glimpse()

### Create temporal variables = Year, Month, day of the year

glimpse(length.metadata$date)

length.metadata<-length.metadata%>%
  separate(col=date,into=c("month","day","year"),sep="/")%>%
  dplyr::select(-day)%>% ## drop day - as we will generate a variable of day of the year - might be use as cyclical smooth in models
  glimpse()

### Convert into factors

length.metadata$year<-as.factor(length.metadata$year)
unique(length.metadata$year) ### 2004 - 2017

length.metadata$month<-as.factor(length.metadata$month)
unique(length.metadata$month<-as.factor(length.metadata$month)) ## January - December

#### Bring in covariates metadata ----

dir()
covariates.metadata<-read.csv("legal.analysis.20200213.csv")%>% ### To bring the corrected covariates from previous datasets
  na_if(-9999)%>%
  dplyr::select(-X)%>%
  dplyr::mutate(id=paste(campaignid,sample,sep="."))%>%
  filter(Taxa=='Legal')%>%
  glimpse()

#### Select variables of interest

names(covariates.metadata)

covariates.metadata<-covariates.metadata%>%
  dplyr::select(id,
                slope,
                relief,
                ausbath,
                t_m,
                t_sd,
                no3_m,
                no3_sd,
                po4_m,
                po4_sd,
                dst2ramp,
                pop.den.50,
                gravity.50,
                method,
                soak.time,
                Realm)%>%
  glimpse()

##### Join length metadata and covariates

names(length.metadata)
names(covariates.metadata)

length.metadata.join<-left_join(length.metadata,covariates.metadata,by="id")%>%
  glimpse()

##### There is an error with left_join in dplyr, and duplicated rows are created

length.metadata.join<-length.metadata.join[!duplicated(length.metadata.join$id), ]

## Check amount of NAs

sum(is.na(length.metadata.join))/prod(dim(length.metadata.join))*100
apply(length.metadata.join,2,function(col)sum(is.na(col))/length(col))*100

##### Check final distribution of samples in the map

map.length.metadata.join<-ggplot(data=length.metadata.join,aes(x=longitude,y=latitude))+
  geom_point()+
  ggtitle("length.metadata.join")+
  theme_classic()

map.length.metadata.join

#### Bring in length data ----

setwd(data.dir)
dir()

length<-read.csv("Au.length.and.biomass.2019-11-14.csv")%>%
  dplyr::filter(!grepl('Christmas.Island|Cocos.Keeling.Islands|LHI',id))%>%
  mutate(scientific=paste(genus,species,sep=" "))%>%
  glimpse()
length$length<-length$length/10
hist(length$length)
sum(length$number) ## 530,271 individuals measured

## Join with Auhstralia life history and filter ----
# Actinopterigii
# Fishing mortality (Commercial and Recreational)
## Demersal or benthic (no pelagic species)
## We will retain only species with minimum legal size regulations

setwd(data.dir)
dir()

lifehistory<-read.csv("life.history.csv")%>%
  dplyr::select(class,scientific,fishing.mortality,fishing.type,rls.water.column,rls.trophic.level,
                minlegal.wa,minlegal.nsw,minlegal.vic,minlegal.tas,minlegal.sa)%>%
  glimpse()

### Include average legal size to life.history

lifehistory$minlegal<-rowMeans(lifehistory[,7:11],na.rm=TRUE)%>%
  glimpse()

## Create filters and use species list to filter out length data

lifehistory<-lifehistory%>%
  dplyr::filter(class%in%c('Actinopterygii'))%>% ## Only ray-finned fishes
  dplyr::filter(grepl('Y',fishing.mortality))%>% ## Only species with fishing mortality
  dplyr::filter(fishing.type!='B')%>% ## Only species recreationally and commercially targeted (no by-catch)
  dplyr::filter(!rls.water.column%in%c('pelagic non-site attached','pelagic site attached'))%>% ## Only demersal and benthic
  dplyr::filter(minlegal!="NaN")%>% ## Exclude species with no minimun legal size
  glimpse()

## Clean dataset

lifehistory<-lifehistory%>%
  dplyr::select(scientific,minlegal)%>%
  glimpse()

## Create a vector of species to filter the length data

species<-as.vector(unique(lifehistory$scientific))

### Join and filter 
### Use this to retain species with an associated minimum legal size

length<-length%>%
  right_join(lifehistory,by=c("scientific"))%>%
  dplyr::filter(!is.na(length))%>%
  mutate(minlegal=minlegal/10)%>% ## Convert minimum legal size to cm
  mutate(mass.g=mass.g/1000)%>% ## Convert minimum legal size to cm
  mutate(Taxa=ifelse(length>=minlegal,'Legal','Sub-legal'))%>% ###For classifying only in legal and sublegal
  glimpse()
hist(length$length.cm)
hist(length$mass.g)
sum(is.na(length$mass.g)) ## ~ 1% of the observations - biomass could not be estimated

## Length - weight relationship for reviewer

ggplot(length,aes(x=length.cm,y=mass.g,colour=scientific))+
  geom_point(show.legend = F)+
  stat_smooth(method = "gam",show.legend = F)+
  theme_classic()

### We will filter out species from the family Monacanthidae - as this have little recreational and commercial interst across states
### Also - some species with little fisheries value - that might be problematic for the analyses - e.g. forming large schools

length<-length%>%
  dplyr::filter(!scientific%in%c("Ophthalmolepis lineolatus","Meuschenia hippocrepis","Meuschenia freycineti","Meuschenia galii",
                                 "Meuschenia flavolineata","Acanthaluteres vittiger",
                                 "Eubalichthys mosaicus","Acanthaluteres brownii",
                                 "Eubalichthys gunnii","Meuschenia australis",
                                 "Meuschenia trachylepis","Meuschenia venusta",
                                 "Nelusetta ayraud","Thamnaconus degeni","Dinolestes lewini"))%>%
  glimpse()

## Create a summary of species ordered by decreasing abundances (MaxN)

test<-length%>%
  dplyr::group_by(scientific)%>%
  dplyr::summarise(N=sum(number))%>%
  glimpse()
unique(length$scientific) ## 82 species
n_distinct(length$scientific) 
unique(length$family) # 18 families
n_distinct(length$family) 
sum(length$number)  # 62,237

## Create a vector with species 90th of the body size distribution

species<-as.vector(unique(length$scientific))

dataframe=data.frame()

for (i in species) {
  
  dat.sum<-length%>%
    filter(scientific%in%i)%>%
    glimpse()
  
  hist(dat.sum$length.cm)
  x<-as.data.frame(quantile(dat.sum$length.cm,probs = c(0.90))) ## 90th percetile of body size distribtion (cm)
  rownames(x)<-i
  colnames(x)<-'nighty.percentile'
  
  dataframe<-rbind(dataframe,x)
    
}

dataframe<-rownames_to_column(dataframe, var = "scientific")

## Join with length data

length<-left_join(length,dataframe,by="scientific")%>%
  glimpse()

## Create another level of Taxa for analyses of the largest fishes

length<-length%>%
  mutate(Taxa.2=ifelse(length.cm>nighty.percentile,"Large",Taxa))%>%
  glimpse()

########## Generate Assemblage-level data analyses frames ----

## Large

length.large<-length%>%
  dplyr::filter(Taxa.2=='Large')%>%
  dplyr::select(campaignid,sample,number,mass.g)%>%
  dplyr::right_join(length.metadata.join,by=c("campaignid","sample"))%>% ## To bring in true zeros
  replace_na(list(number = 0))%>%## Fill in true zeros
  replace_na(list(mass.g = 0))%>%## Fill in true zeros
  dplyr::group_by(campaignid,sample)%>%
  dplyr::summarise(number=sum(number,na.rm=TRUE),
                   biomass=sum(mass.g,na.rm = TRUE))%>%
  dplyr::right_join(length.metadata.join,by=c("campaignid","sample"))%>%
  dplyr::mutate(Taxa="Large")%>%
  glimpse()

##### Check proportion of NAs in the final dataset

sum(is.na(length.large))/prod(dim(length.large))*100
apply(length.large,2,function(col)sum(is.na(col))/length(col))*100

## Export csv

write.csv(length.large,"large.analysis.assemblage-level.csv")

## Legal

length.legal<-length%>%
  dplyr::filter(Taxa=='Legal')%>%
  dplyr::select(campaignid,sample,number,mass.g)%>%
  dplyr::right_join(length.metadata.join,by=c("campaignid","sample"))%>% ## To bring in true zeros
  replace_na(list(number = 0))%>%## Fill in true zeros
  replace_na(list(mass.g = 0))%>%## Fill in true zeros
  dplyr::group_by(campaignid,sample)%>%
  dplyr::summarise(number=sum(number),
                   biomass=sum(mass.g,na.rm = TRUE))%>%
  dplyr::right_join(length.metadata.join,by=c("campaignid","sample"))%>%
  dplyr::mutate(Taxa="Legal")%>%
  glimpse()

## Sub-legal

length.sublegal<-length%>%
  dplyr::filter(Taxa=='Sub-legal')%>%
  dplyr::select(campaignid,sample,number,mass.g)%>%
  dplyr::right_join(length.metadata.join,by=c("campaignid","sample"))%>% ## To bring in true zeros
  replace_na(list(number = 0))%>%## Fill in true zeros
  replace_na(list(mass.g = 0))%>%## Fill in true zeros
  dplyr::group_by(campaignid,sample)%>%
  dplyr::summarise(number=sum(number),
                   biomass=sum(mass.g,na.rm = TRUE))%>%
  dplyr::right_join(length.metadata.join,by=c("campaignid","sample"))%>%
  dplyr::mutate(Taxa="Sub-legal")%>%
  glimpse()

### Join them back together

length.analysis<-bind_rows(length.legal,length.sublegal)%>%
  glimpse()

##### Check proportion of NAs in the final dataset

sum(is.na(length.analysis))/prod(dim(length.analysis))*100
apply(length.analysis,2,function(col)sum(is.na(col))/length(col))*100

## Export csv

write.csv(length.analysis,"length.analysis.assemblage-level.csv")

## Explore the correlation between the metrics

cor(length.legal$number,length.legal$biomass)
cor(length.sublegal$number,length.sublegal$biomass)

## Clean working directory

rm(length.analysis,length.large,length.legal,length.sublegal)

####### Generate regional species groups analyses dataframes -----

## Create a species groups variable

length<-length%>%
  mutate(species.group=ifelse(scientific%in%c('Lethrinus nebulosus','Lethrinus spp','Lethrinus atkinsoni','Lethrinus punctulatus',
                                'Lethrinus miniatus','Lethrinus laticaudis','Lethrinus microdon','Lethrinus olivaceus'),"Lethrinus spp.",
                ifelse(scientific%in%c('Chrysophrys auratus'),'Chrysophrys auratus',
                       ifelse(scientific%in%c('Lutjanus sebae','Lutjanus carponotatus','Lutjanus argentimaculatus',
                                              'Lutjanus lemniscatus','Lutjanus malabaricus','Lutjanus russellii'),'Lutjanus spp.',
                              ifelse(scientific%in%c('Choerodon cyanodus','Choerodon rubescens','Choerodon schoenleinii'),'Choerodon spp.',
                                     ifelse(scientific%in%c('Plectropomus maculatus', 'Plectropomus areolatus', 'Plectropomus leopardus','Plectropomus spp'),'Plectropomus spp.',
                                            ifelse(scientific%in%c('Nemadactylus douglasii','Nemadactylus macropterus','Nemadactylus valenciennesi'),'Nemadactylus spp.',
                                                   ifelse(scientific%in%c('Notolabrus fucicola','Notolabrus tetricus'),'Notolabrus spp.','Others'))))))))%>%
  glimpse()
unique(length$species.group)

########## Generate analyses data analyses

## Large

length.large<-length%>%
  filter(!species.group%in%c('Others'))%>%
  dplyr::filter(Taxa.2=='Large')%>%
  dplyr::select(campaignid,sample,species.group,number,mass.g)%>%
  dplyr::right_join(length.metadata.join,by=c("campaignid","sample"))%>% ## To bring in true zeros
  complete(species.group,nesting(campaignid,sample),fill = list(number = 0,
                                                       mass.g = 0))%>%
  dplyr::group_by(campaignid,sample,species.group)%>%
  dplyr::summarise(number=sum(number,na.rm=TRUE),
                   biomass=sum(mass.g,na.rm = TRUE))%>%
  dplyr::right_join(length.metadata.join,by=c("campaignid","sample"))%>%
  dplyr::mutate(Taxa="Large")%>%
  filter(!is.na(species.group))%>%
  glimpse()
unique(length.large$species.group)

##### Check proportion of NAs in the final dataset

sum(is.na(length.large))/prod(dim(length.large))*100
apply(length.large,2,function(col)sum(is.na(col))/length(col))*100

## Export csv

write.csv(length.large,"large.analysis.species.groups.csv")

## Legal

length.legal<-length%>%
  filter(!species.group%in%c('Others'))%>%
  dplyr::filter(Taxa=='Legal')%>%
  dplyr::select(campaignid,sample,species.group,number,mass.g)%>%
  dplyr::right_join(length.metadata.join,by=c("campaignid","sample"))%>% ## To bring in true zeros
  complete(species.group,nesting(campaignid,sample),fill = list(number = 0,
                                                                mass.g = 0))%>%
  dplyr::group_by(campaignid,sample,species.group)%>%
  dplyr::summarise(number=sum(number,na.rm=TRUE),
                   biomass=sum(mass.g,na.rm = TRUE))%>%
  dplyr::right_join(length.metadata.join,by=c("campaignid","sample"))%>%
  dplyr::mutate(Taxa="Legal")%>%
  filter(!is.na(species.group))%>%
  glimpse()
unique(length.legal$species.group)

## Sub-legal

length.sublegal<-length%>%
  filter(!species.group%in%c('Others'))%>%
  dplyr::filter(Taxa=='Sub-legal')%>%
  dplyr::select(campaignid,sample,species.group,number,mass.g)%>%
  dplyr::right_join(length.metadata.join,by=c("campaignid","sample"))%>% ## To bring in true zeros
  complete(species.group,nesting(campaignid,sample),fill = list(number = 0,
                                                                mass.g = 0))%>%
  dplyr::group_by(campaignid,sample,species.group)%>%
  dplyr::summarise(number=sum(number,na.rm=TRUE),
                   biomass=sum(mass.g,na.rm = TRUE))%>%
  dplyr::right_join(length.metadata.join,by=c("campaignid","sample"))%>%
  dplyr::mutate(Taxa='Sub-legal')%>%
  filter(!is.na(species.group))%>%
  glimpse()
unique(length.sublegal$species.group)

### Join them back together

length.analysis<-bind_rows(length.legal,length.sublegal)%>%
  glimpse()

##### Check proportion of NAs in the final dataset

sum(is.na(length.analysis))/prod(dim(length.analysis))*100
apply(length.analysis,2,function(col)sum(is.na(col))/length(col))*100

## Export csv

write.csv(length.analysis,"length.analysis.species.group.csv")

## Explore the correlation between the metrics

cor(length.legal$number,length.legal$biomass)
cor(length.sublegal$number,length.sublegal$biomass)

## Supporting Information Table S2 ----
## Species-specific MLS, size, and depth distribution

### Generate a species list

species.list<-as.data.frame(unique(length$scientific))%>%
  glimpse()
colnames(species.list)<-'scientific'
species<-as.vector(unique(species.list$scientific))

### Join Minimun legal sizes by state

setwd(data.dir)
dir()

lifehistory<-read.csv("life.history.csv")%>%
  dplyr::select(family,scientific, minlegal.wa,minlegal.nsw,minlegal.vic,minlegal.tas,minlegal.sa)%>%
  glimpse()

species.list<-left_join(lifehistory,species.list,by=c("scientific"))%>%
  dplyr::filter(scientific%in%species)%>%
  glimpse()

## Obtain symmaries of lengths

glimpse(length)

sampled.depths<-length.legal%>%
  dplyr::select(id,depth)%>%
  ungroup()%>%
  glimpse()

species.summaries<-length%>%
  right_join(sampled.depths,by="id")%>%
  dplyr::group_by(scientific)%>%
  dplyr::summarise(length_min=min(length.cm),
                   length_max=max(length.cm),
                   length_mean=mean(length.cm,na.rm=T),
                   length_sd=sd(length.cm,na.rm=T),
                   depth_min=min(depth),
                   depth_max=max(depth),
                   depth_mean=mean(depth,na.rm=T),
                   depth_sd=sd(depth,na.rm = T))%>%
  glimpse()

## with life history

species.summaries<-left_join(species.summaries,species.list,by="scientific")%>%
  glimpse()

### Convert MLS to cm

species.summaries$minlegal.wa<-species.summaries$minlegal.wa/10
species.summaries$minlegal.nsw<-species.summaries$minlegal.nsw/10
species.summaries$minlegal.vic<-species.summaries$minlegal.vic/10
species.summaries$minlegal.tas<-species.summaries$minlegal.tas/10
species.summaries$minlegal.sa<-species.summaries$minlegal.sa/10

## Round decimals

species.summaries<-species.summaries%>%
  mutate_at(vars(2:9), round, 0)%>%
  glimpse()

## Reorder columns

names(species.summaries)

col_order<-c("family","scientific","minlegal.wa","minlegal.nsw","minlegal.vic","minlegal.sa","minlegal.tas",
             "length_min","length_max","length_mean","length_sd","depth_min","depth_max","depth_mean","depth_sd")

species.summaries<-species.summaries[,col_order]

## Join columns

species.summaries<-species.summaries%>%
  dplyr::mutate(Size=paste(length_mean,length_sd,sep=" ± "),
                Size_range=paste(length_min,length_max,sep=" - "),
                Depth=paste(depth_mean,depth_sd,sep=" ± "),
                Depth_range=paste(depth_min,depth_max,sep=" - "))%>%
  dplyr::select(-length_mean,-length_sd,-length_min,-length_max,-depth_mean,-depth_sd,
                -depth_min,-depth_max)%>%
  glimpse()

species.summaries<-species.summaries%>%
  dplyr::mutate(Size=paste(Size_range,Size,sep=" ("),
                Depth=paste(Depth_range,Depth,sep=" ("))%>%
  dplyr::select(-Size_range,-Depth_range)%>%
  dplyr::mutate(p=")")%>%
  dplyr::mutate(Size=paste(Size,p,sep=""),
                Depth=paste(Depth,p,sep = ""))%>%
  dplyr::select(-p)%>%
  arrange(family)%>%
  dplyr::filter(!is.na(family))%>%
  glimpse()

## Write Table S2 to the data folder

setwd(data.dir)
dir()

write.csv(species.summaries,"Table.S2.Species_summaries.csv",quote=TRUE)

### Create species contributions dataset

sum(length$number)

species.contributions<-length%>%
  dplyr::group_by(family,scientific)%>%
  dplyr::summarise(N=sum(number),
                   prop=(N/62237)*100)%>%
  arrange(desc(prop))%>%
  glimpse()

species.contributions$cum<-cumsum(species.contributions$prop)%>%
  glimpse()

## Export as a csv file

setwd(data.dir)
dir()

write.csv(species.contributions,"Table.SX.species.contributions.csv")

##################################################################################################
############################# (2) Generating a snapshot of the data analyses #####################

#### (A) Assemblage-level dataframes ----
#### Note - this is carried out for each dataframe - not repeated in the code

## Bring dataset

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.assemblage-level.csv")%>% ## Replace your dataframe here!
  dplyr::select(-X)%>%
  na.omit()%>%
  glimpse()

## Calculate linear distance between samples

coordinates(dat)<-~longitude+latitude
class(dat)

proj4string(dat)<-CRS("+init=epsg:4326")
dat@proj4string

## For distance based on sphere through the geosphere package

## Distance based in projected coordinates - 1994 GDA Lambert

?spDists

## Transform dat to Lambert

dat<-spTransform(dat,CRS("+init=epsg:3112"))
proj4string(dat)
coordinates(dat)

d.matrix <- spDists(dat)
min(d.matrix) ## 0 duplicated samples
max(d.matrix) ## 3,988,200

## Sample clusters - within 250 m

hc <- hclust(as.dist(d.matrix), method="complete")

d=500

dat$clust.500m <- cutree(hc, h=d)
glimpse(dat@data)
unique(dat$clust.500m)
n_distinct(dat$clust.500m) ## 3,154 unique samples within 500 m

## Reconvert to WGS84 prior to converting to a dataframe

dat<-spTransform(dat,CRS("+init=epsg:4326"))
proj4string(dat)
coordinates(dat)

## Reconvert into a dataframe

dat<-as.data.frame(dat)%>%
  glimpse()
dat$clust.500m<-as.factor(dat$clust.500m)

### Filter closest year to 2011

ggplot(dat,aes(x=year,fill=Realm))+
  geom_histogram(colour="black")+
  facet_wrap(~Realm)+
  theme_bw()

dat<-dat%>%
  group_by(clust.500m,status,Taxa)%>%
  slice(which.min(abs(year-2011)))%>% ## Samples are not excluded because these correspond to status
  ungroup()%>%
  glimpse()
unique(dat$campaignid)
n_distinct(dat$campaignid) ## 95 CampaignIDs
n_distinct(dat$id) ## 3,172 stereo-BRUVs deployments

# Find if clusters have more than one year

test<-dat%>%
  group_by(clust.500m)%>%
  summarise(N=n_distinct(year))%>%
  glimpse()
sum(test$N>1) ## 8 clusters that had more than one deployment - but this correspond to Fished and NTR

## Find if clusters have more than one year

test<-dat%>%
  group_by(clust.500m)%>%
  summarise(N=n_distinct(month))%>%
  glimpse()
sum(test$N>1) ## 8 clusters that had more than one deployment - but this correspond to Fished and NTR

## Find clusters of deployments within 4-km of each other - in case we need to use as a random effect in the models

## Calculate linear distance between cluster samples

coordinates(dat)<-~longitude+latitude
class(dat)

proj4string(dat)<-CRS("+init=epsg:4326")
dat@proj4string

?spDists

## Transform dat to Lambert

dat<-spTransform(dat,CRS("+init=epsg:3112"))
proj4string(dat)
coordinates(dat)

d.matrix <- spDists(dat)
min(d.matrix) ## 0 duplicated samples
max(d.matrix) ## 4,074.416

## Find clusters based on distance apart - hierarchical clustering with complete linkage

hc <- hclust(as.dist(d.matrix), method="complete")

## Reef cluster 4 km

d=4000

dat$clust.4km <- cutree(hc, h=d)
glimpse(dat)
unique(dat$clust.4km)
n_distinct(dat$clust.4km) ## 687 reef clusters (4 km)

## Reconvert to WGS84 prior to converting to a dataframe

dat<-spTransform(dat,CRS("+init=epsg:4326"))
proj4string(dat)
coordinates(dat)

## Reconvert into a dataframe

dat<-as.data.frame(dat)%>%
  glimpse()

## Convert cluster 4 km into a factor

dat$clust.4km<-as.factor(dat$clust.4km)

## Export coordinates to download NPP (min,max,mean,sd) from MSEC

coordinates<-dat%>%
  select(id,longitude,latitude)%>%
  glimpse()

write.csv(coordinates,"coordinates.csv")

####Bring in Global legal dataset - extract SEASYNC variables

setwd(data.dir)
dir()

seasync<-read.csv("MSEC.csv")%>%
  #dplyr::select(id,npp_mean,npp_min,npp_max,npp_sd)%>%
  distinct(id,.keep_all = T)%>%
  glimpse()

dat<-left_join(dat,seasync,by="id")%>%
  filter(!is.na(npp_mean))%>%
  glimpse()

##### Check proportion of NAs in the final dataset

sum(is.na(dat))/prod(dim(dat))*100
apply(dat,2,function(col)sum(is.na(col))/length(col))*100

## Export final dataset

setwd(data.dir)
write.csv(dat,"length.analysis.assemblage-level.csv")

## (B) Regional-species groups ----
## Note - this is done repetitively for each dataframe (i.e. large vs. legal and sub-legal)

## (B.1) Lethrinus spp. ----

## Bring in samples id from snapshot of the assemblage-level model

setwd(data.dir)
dir()

sample.id<-read.csv("length.analysis.assemblage-level.csv")%>%
  dplyr::select(id)%>%
  glimpse()
x<-as.vector(unique(sample.id$id))

## Bring species group data - filter deployments only present in the temporal snapshot

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.species.group.csv")%>%
  dplyr::select(-X)%>%
  filter(species.group%in%c('Lethrinus spp.'))%>%
  filter(id%in%x)%>%
  na.omit()%>%
  glimpse()
unique(dat$species.group)
n_distinct(dat$id)

#### Trimmed each individual dataframe to the geographic range of occurrence 
## Filter <1% Probability of occurrence to filter out records of vagrant indivuals

## First we will need to group legal and sub-legal taxa 

test.dat<-dat%>%
  group_by(state,campaignid,id,longitude,latitude,species.group)%>%
  summarise(response=sum(number))%>%
  glimpse()

## Create presence/absence dataset

test.dat.incidence<-test.dat%>%
  mutate(response=ifelse(response>0,1,response))%>%
  glimpse()

## Create a dataframe with presence only data

test.dat.presence<-test.dat.incidence%>%
  filter(response>0)%>%
  glimpse()

## Check distributional of presence/absence data

map<-ggplot()+
  geom_point(data=test.dat.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=test.dat.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,160))+
  scale_y_continuous(expand = c(0,0),limits = c(-45,-9))+
  ggtitle('(a) Lethrinus spp.')+
  theme_classic()+
  Theme1
map

## Subset dataset

x<-test.dat.incidence%>%
  arrange(desc(latitude))%>%
  glimpse()

quantile(test.dat.presence$latitude, c(0.01, .50, .99)) ## Find percentiles

test.dat.incidence<-test.dat.incidence%>%
  dplyr::filter(latitude>=(-29))%>%
  dplyr::filter(latitude<=(-13.92173))%>%
  dplyr::filter(longitude<130)%>%
  glimpse()
n_distinct(test.dat.incidence$id) ## 1,616 samples

test.dat.presence<-test.dat.presence%>%
  dplyr::filter(latitude>=(-29))%>%
  dplyr::filter(latitude<=(-13.92173))%>%
  dplyr::filter(longitude<130)%>%
  glimpse()
n_distinct(test.dat.presence$id) ## 185 samples

## Plot again to cross-check

map<-ggplot()+
  geom_point(data=test.dat.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=test.dat.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,130))+
  scale_y_continuous(expand = c(0,0),limits=c(-30,-13))+
  ggtitle('(a) Letrhinus spp.')+
  theme_classic()+
  Theme1
map

## Get a vector of unique Campaignids were the species was present (i.e. samples)

samples<-as.vector(unique(test.dat.incidence$id))

## Subset original dataset to stereo-BRUVs deployments within the biogeographic range of occurrence

dat<-dat%>%
  filter(id%in%samples)%>%
  glimpse()
sum(dat$number == 0 ) / length(dat$number) ## 88 % zeroes

## Store dataset to join latter

lethrinus.dat<-dat

## (B.2) Chrysophrys auratus ----

sample.id<-read.csv("length.analysis.assemblage-level.csv")%>%
  dplyr::select(id)%>%
  glimpse()
x<-as.vector(unique(sample.id$id))

## Bring species group data - filter deployments only present in the temporal snapshot

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.species.group.csv")%>%
  dplyr::select(-X)%>%
  filter(species.group%in%c('Chrysophrys auratus'))%>%
  filter(id%in%x)%>%
  na.omit()%>%
  glimpse()
unique(dat$species.group)
n_distinct(dat$id)

#### Trimmed each individual dataframe to the geographic range of occurrence 
## Filter <1% Probability of occurrence to filter out records of vagrant indivuals

## First we will need to group legal and sub-legal taxa 

test.dat<-dat%>%
  group_by(state,campaignid,id,longitude,latitude,species.group)%>%
  summarise(response=sum(number))%>%
  glimpse()

## Create presence/absence dataset

test.dat.incidence<-test.dat%>%
  mutate(response=ifelse(response>0,1,response))%>%
  glimpse()

## Create a dataframe with presence only data

test.dat.presence<-test.dat.incidence%>%
  filter(response>0)%>%
  glimpse()

## Check distributional of presence/absence data

map<-ggplot()+
  geom_point(data=test.dat.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=test.dat.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,160))+
  scale_y_continuous(expand = c(0,0),limits = c(-45,-9))+
  ggtitle('(b) Chrysophrys auratus')+
  theme_classic()+
  Theme1
map

## Subset dataset

x<-test.dat.incidence%>%
  arrange(desc(latitude))%>%
  glimpse()

quantile(test.dat.presence$latitude, c(0.01, .50, .99)) ## Find percentiles

test.dat.incidence<-test.dat.incidence%>%
  filter(latitude<=(-25))%>%
  glimpse()
n_distinct(test.dat.incidence$id) ## 1,808 samples

test.dat.presence<-test.dat.presence%>%
  filter(latitude<=(-25))%>%
  glimpse()
n_distinct(test.dat.presence$id) ## 136 samples

## Plot again to cross-check

map<-ggplot()+
  geom_point(data=test.dat.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=test.dat.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,160))+
  scale_y_continuous(expand = c(0,0),limits=c(-45,-25))+
  ggtitle('(b) Chrysophurus auratus')+
  theme_classic()+
  Theme1
map

## Get a vector of unique Campaignids were the species was present (i.e. samples)

samples<-as.vector(unique(test.dat.incidence$id))

## Subset original dataset to stereo-BRUVs deployments within the biogeographic range of occurrence

dat<-dat%>%
  filter(id%in%samples)%>%
  glimpse()
sum(dat$number == 0 ) / length(dat$number) ## 81 % zeroes

## Store dataset to join latter

snapper.dat<-dat

## (B.3) Lutjanus spp. ----

sample.id<-read.csv("length.analysis.assemblage-level.csv")%>%
  dplyr::select(id)%>%
  glimpse()
x<-as.vector(unique(sample.id$id))

## Bring species group data - filter deployments only present in the temporal snapshot

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.species.group.csv")%>%
  dplyr::select(-X)%>%
  filter(species.group%in%c('Lutjanus spp.'))%>%
  filter(id%in%x)%>%
  na.omit()%>%
  glimpse()
unique(dat$species.group)
n_distinct(dat$id)

#### Trimmed each individual dataframe to the geographic range of occurrence 
## Filter <1% Probability of occurrence to filter out records of vagrant indivuals

## First we will need to group legal and sub-legal taxa 

test.dat<-dat%>%
  group_by(state,campaignid,id,longitude,latitude,species.group)%>%
  summarise(response=sum(number))%>%
  glimpse()

## Create presence/absence dataset

test.dat.incidence<-test.dat%>%
  mutate(response=ifelse(response>0,1,response))%>%
  glimpse()

## Create a dataframe with presence only data

test.dat.presence<-test.dat.incidence%>%
  filter(response>0)%>%
  glimpse()

## Check distributional of presence/absence data

map<-ggplot()+
  geom_point(data=test.dat.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=test.dat.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,160))+
  scale_y_continuous(expand = c(0,0),limits = c(-45,-9))+
  ggtitle('(c) Lutjanus spp.')+
  theme_classic()+
  Theme1
map

## Subset dataset

x<-test.dat.incidence%>%
  arrange(desc(latitude))%>%
  glimpse()

quantile(test.dat.presence$latitude, c(0.01, .50, .99)) ## Find percentiles

test.dat.incidence<-test.dat.incidence%>%
  dplyr::filter(latitude>=(-24))%>%
  dplyr::filter(longitude<130)%>%
  glimpse()
n_distinct(test.dat.incidence$id) ## 1,335 samples

test.dat.presence<-test.dat.presence%>%
  dplyr::filter(latitude>=(-24))%>%
  dplyr::filter(longitude<130)%>%
  glimpse()
n_distinct(test.dat.presence$id) ## 124 samples

## Plot again to cross-check

map<-ggplot()+
  geom_point(data=test.dat.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=test.dat.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,130))+
  scale_y_continuous(expand = c(0,0),limits=c(-30,-13))+
  ggtitle('(c) Lutjanus spp.')+
  theme_classic()+
  Theme1
map

## Get a vector of unique Campaignids were the species was present (i.e. samples)

samples<-as.vector(unique(test.dat.incidence$id))

## Subset original dataset to stereo-BRUVs deployments within the biogeographic range of occurrence

dat<-dat%>%
  filter(id%in%samples)%>%
  glimpse()
sum(dat$number == 0 ) / length(dat$number) ## 90 % zeroes

## Store dataset to join latter

Lutjanus.dat<-dat

## (B.4) Choerodon spp. ----

sample.id<-read.csv("length.analysis.assemblage-level.csv")%>%
  dplyr::select(id)%>%
  glimpse()
x<-as.vector(unique(sample.id$id))

## Bring species group data - filter deployments only present in the temporal snapshot

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.species.group.csv")%>%
  dplyr::select(-X)%>%
  filter(species.group%in%c('Choerodon spp.'))%>%
  filter(id%in%x)%>%
  na.omit()%>%
  glimpse()
unique(dat$species.group)
n_distinct(dat$id)

#### Trimmed each individual dataframe to the geographic range of occurrence 
## Filter <1% Probability of occurrence to filter out records of vagrant indivuals

## First we will need to group legal and sub-legal taxa 

test.dat<-dat%>%
  group_by(state,campaignid,id,longitude,latitude,species.group)%>%
  summarise(response=sum(number))%>%
  glimpse()

## Create presence/absence dataset

test.dat.incidence<-test.dat%>%
  mutate(response=ifelse(response>0,1,response))%>%
  glimpse()

## Create a dataframe with presence only data

test.dat.presence<-test.dat.incidence%>%
  filter(response>0)%>%
  glimpse()

## Check distributional of presence/absence data

map<-ggplot()+
  geom_point(data=test.dat.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=test.dat.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,160))+
  scale_y_continuous(expand = c(0,0),limits = c(-45,-9))+
  ggtitle('(d) Choerodon spp.')+
  theme_classic()+
  Theme1
map

## Subset dataset

x<-test.dat.incidence%>%
  arrange(desc(latitude))%>%
  glimpse()

quantile(test.dat.presence$latitude, c(0.01, .50, .99)) ## Find percentiles

test.dat.incidence<-test.dat.incidence%>%
  dplyr::filter(latitude>=(-33.53508))%>%
  dplyr::filter(longitude<130)%>%
  glimpse()
n_distinct(test.dat.incidence$id) ## 2171 samples

test.dat.presence<-test.dat.presence%>%
  dplyr::filter(latitude>=(-33.53508))%>%
  dplyr::filter(longitude<130)%>%
  glimpse()
n_distinct(test.dat.presence$id) ## 170 samples

## Plot again to cross-check

map<-ggplot()+
  geom_point(data=test.dat.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=test.dat.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,130))+
  scale_y_continuous(expand = c(0,0),limits=c(-35,-13))+
  ggtitle('(d) Choerodon spp.')+
  theme_classic()+
  Theme1
map

## Get a vector of unique Campaignids were the species was present (i.e. samples)

samples<-as.vector(unique(test.dat.incidence$id))

## Subset original dataset to stereo-BRUVs deployments within the biogeographic range of occurrence

dat<-dat%>%
  filter(id%in%samples)%>%
  glimpse()
sum(dat$number == 0 ) / length(dat$number) ## 92 % zeroes

## Store dataset to join latter

Choerodon.dat<-dat

## (B.5) Plectropomus spp. ----

sample.id<-read.csv("length.analysis.assemblage-level.csv")%>%
  dplyr::select(id)%>%
  glimpse()
x<-as.vector(unique(sample.id$id))

## Bring species group data - filter deployments only present in the temporal snapshot

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.species.group.csv")%>%
  dplyr::select(-X)%>%
  filter(species.group%in%c('Plectropomus spp.'))%>%
  filter(id%in%x)%>%
  na.omit()%>%
  glimpse()
unique(dat$species.group)
n_distinct(dat$id)

#### Trimmed each individual dataframe to the geographic range of occurrence 
## Filter <1% Probability of occurrence to filter out records of vagrant indivuals

## First we will need to group legal and sub-legal taxa 

test.dat<-dat%>%
  group_by(state,campaignid,id,longitude,latitude,species.group)%>%
  summarise(response=sum(number))%>%
  glimpse()

## Create presence/absence dataset

test.dat.incidence<-test.dat%>%
  mutate(response=ifelse(response>0,1,response))%>%
  glimpse()

## Create a dataframe with presence only data

test.dat.presence<-test.dat.incidence%>%
  filter(response>0)%>%
  glimpse()

## Check distributional of presence/absence data

map<-ggplot()+
  geom_point(data=test.dat.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=test.dat.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,160))+
  scale_y_continuous(expand = c(0,0),limits = c(-45,-9))+
  ggtitle('(e) Plectropomus spp.')+
  theme_classic()+
  Theme1
map

## Subset dataset

x<-test.dat.incidence%>%
  arrange(desc(latitude))%>%
  glimpse()

quantile(test.dat.presence$latitude, c(0.01, .50, .99)) ## Find percentiles

test.dat.incidence<-test.dat.incidence%>%
  dplyr::filter(latitude>=(-29))%>%
  dplyr::filter(latitude<=(-13.92173))%>%
  dplyr::filter(longitude<130)%>%
  glimpse()
n_distinct(test.dat.incidence$id) ## 1,616 samples

test.dat.presence<-test.dat.presence%>%
  dplyr::filter(latitude>=(-29))%>%
  dplyr::filter(latitude<=(-13.92173))%>%
  dplyr::filter(longitude<130)%>%
  glimpse()
n_distinct(test.dat.presence$id) ## 73 samples

## Plot again to cross-check

map<-ggplot()+
  geom_point(data=test.dat.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=test.dat.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,130))+
  scale_y_continuous(expand = c(0,0),limits=c(-30,-13))+
  ggtitle('(e) Plectropomus spp.')+
  theme_classic()+
  Theme1
map

## Get a vector of unique Campaignids were the species was present (i.e. samples)

samples<-as.vector(unique(test.dat.incidence$id))

## Subset original dataset to stereo-BRUVs deployments within the biogeographic range of occurrence

dat<-dat%>%
  filter(id%in%samples)%>%
  glimpse()
sum(dat$number == 0 ) / length(dat$number) ## 95 % zeroes

## Store dataset to join latter

Plectropomus.dat<-dat

## (B.6) Nemadactylus spp. ----

sample.id<-read.csv("length.analysis.assemblage-level.csv")%>%
  dplyr::select(id)%>%
  glimpse()
x<-as.vector(unique(sample.id$id))

## Bring species group data - filter deployments only present in the temporal snapshot

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.species.group.csv")%>%
  dplyr::select(-X)%>%
  filter(species.group%in%c('Nemadactylus spp.'))%>%
  filter(id%in%x)%>%
  na.omit()%>%
  glimpse()
unique(dat$species.group)
n_distinct(dat$id)

#### Trimmed each individual dataframe to the geographic range of occurrence 
## Filter <1% Probability of occurrence to filter out records of vagrant indivuals

## First we will need to group legal and sub-legal taxa 

test.dat<-dat%>%
  group_by(state,campaignid,id,longitude,latitude,species.group)%>%
  summarise(response=sum(number))%>%
  glimpse()

## Create presence/absence dataset

test.dat.incidence<-test.dat%>%
  mutate(response=ifelse(response>0,1,response))%>%
  glimpse()

## Create a dataframe with presence only data

test.dat.presence<-test.dat.incidence%>%
  filter(response>0)%>%
  glimpse()

## Check distributional of presence/absence data

map<-ggplot()+
  geom_point(data=test.dat.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=test.dat.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,160))+
  scale_y_continuous(expand = c(0,0),limits = c(-45,-9))+
  ggtitle('(f) Nemadactylus spp.')+
  theme_classic()+
  Theme1
map

## Subset dataset

x<-test.dat.incidence%>%
  arrange(desc(latitude))%>%
  glimpse()

quantile(test.dat.presence$latitude, c(0.01, .50, .99)) ## Find percentiles

test.dat.incidence<-test.dat.incidence%>%
  filter(latitude<=-29.77585)%>%
  glimpse()
n_distinct(test.dat.incidence$id) ## 1,493 samples

test.dat.presence<-test.dat.presence%>%
  filter(latitude<=-29.77585)%>%
  glimpse()
n_distinct(test.dat.presence$id) ## 60 samples

## Plot again to cross-check

map<-ggplot()+
  geom_point(data=test.dat.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=test.dat.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,160))+
  scale_y_continuous(expand = c(0,0),limits=c(-45,-13))+
  ggtitle('(f) Nemadactylus spp.')+
  theme_classic()+
  Theme1
map

## Get a vector of unique Campaignids were the species was present (i.e. samples)

samples<-as.vector(unique(test.dat.incidence$id))

## Subset original dataset to stereo-BRUVs deployments within the biogeographic range of occurrence

dat<-dat%>%
  filter(id%in%samples)%>%
  glimpse()
sum(dat$number == 0 ) / length(dat$number) ## 95 % zeroes

## Store dataset to join latter

Nemadactylus.dat<-dat

## (B.7) Notolabrus spp. ----

sample.id<-read.csv("length.analysis.assemblage-level.csv")%>%
  dplyr::select(id)%>%
  glimpse()
x<-as.vector(unique(sample.id$id))

## Bring species group data - filter deployments only present in the temporal snapshot

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.species.group.csv")%>%
  dplyr::select(-X)%>%
  filter(species.group%in%c('Notolabrus spp.'))%>%
  filter(id%in%x)%>%
  na.omit()%>%
  glimpse()
unique(dat$species.group)
n_distinct(dat$id)

#### Trimmed each individual dataframe to the geographic range of occurrence 
## Filter <1% Probability of occurrence to filter out records of vagrant indivuals

## First we will need to group legal and sub-legal taxa 

test.dat<-dat%>%
  group_by(state,campaignid,id,longitude,latitude,species.group)%>%
  summarise(response=sum(number))%>%
  glimpse()

## Create presence/absence dataset

test.dat.incidence<-test.dat%>%
  mutate(response=ifelse(response>0,1,response))%>%
  glimpse()

## Create a dataframe with presence only data

test.dat.presence<-test.dat.incidence%>%
  filter(response>0)%>%
  glimpse()

## Check distributional of presence/absence data

map<-ggplot()+
  geom_point(data=test.dat.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=test.dat.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,160))+
  scale_y_continuous(expand = c(0,0),limits = c(-45,-9))+
  ggtitle('(g) Notolabrus spp.')+
  theme_classic()+
  Theme1
map

## Subset dataset

x<-test.dat.incidence%>%
  arrange(desc(latitude))%>%
  glimpse()

quantile(test.dat.presence$latitude, c(0.01, .50, .99)) ## Find percentiles

test.dat.incidence<-test.dat.incidence%>%
  filter(state%in%c('VIC','SA','TAS'))%>%
  glimpse()
n_distinct(test.dat.incidence$id) ## 439 samples

test.dat.presence<-test.dat.presence%>%
  filter(state%in%c('VIC','SA','TAS'))%>%
  glimpse()
n_distinct(test.dat.presence$id) ## 62 samples

## Plot again to cross-check

map<-ggplot()+
  geom_point(data=test.dat.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=test.dat.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,160))+
  scale_y_continuous(expand = c(0,0),limits=c(-45,-13))+
  ggtitle('(g) Notolabrus spp.')+
  theme_classic()+
  Theme1
map

## Get a vector of unique Campaignids were the species was present (i.e. samples)

samples<-as.vector(unique(test.dat.incidence$id))

## Subset original dataset to stereo-BRUVs deployments within the biogeographic range of occurrence

dat<-dat%>%
  filter(id%in%samples)%>%
  glimpse()
sum(dat$number == 0 ) / length(dat$number) ## 85 % zeroes

## Store dataset to join latter

Notolabrus.dat<-dat

## Join all datasets

dat.final<-bind_rows(lethrinus.dat,snapper.dat,Lutjanus.dat,
                     Choerodon.dat,Plectropomus.dat,Nemadactylus.dat,Notolabrus.dat)%>%
  glimpse()
unique(dat.final$species.group)

## Export dataset

setwd(data.dir)
dir()

write.csv(dat.final,"length.analysis.species.group.csv")

##################################################################################################
########################### (3) GAMM models - Assemblage-level ###################################

## (A) Large abundance ----

?full.subsets.gam

# Read in data 

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.assemblage-level.csv")%>%
  filter(relief<200)%>% ## Some problems with few observations close to coast - with extreme values of relief
  rename(response=number)%>%
  filter(Taxa%in%c('Large'))%>%
  na.omit()%>%
  glimpse()
unique(dat$Taxa)

## Explore the temporal distribution of the data

min(dat$year)
max(dat$year)

year.plot<-ggplot(dat,aes(x=year,fill=Realm))+
  geom_histogram(colour="black")+
  facet_wrap(~Realm)+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  Theme1
year.plot

min(dat$month)
max(dat$month)

month.plot<-ggplot(dat,aes(x=month,fill=Realm))+
  geom_histogram(colour="black")+
  facet_wrap(~Realm)+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  Theme1
month.plot

## Arrange all plots in a grob

ggarrange(year.plot,month.plot,
          nrow = 2,ncol = 1,align = 'hv',common.legend = T,legend = "none")

## Export plot

setwd("C:/Users/22373243/Dropbox/Projects/Analysis/Analysis_GlobalArchive_Bodysize_Nestor/Plots") 
dir()

name<-'Temporal_Distribution'

ggsave(paste(name,"Top_ranked_20201002.png",sep="."),width = 21, height = 15,units = "cm",dpi = 600)
ggsave(paste(name,"Top_ranked_20201002.pdf",sep="."),width = 21, height = 15,units = "cm",dpi = 600,useDingbats=FALSE)

## Explore spatial variation in gravity values

ggplot(dat,aes(x=longitude,y=latitude,size=log1p(gravity.50),colour=log1p(gravity.50)))+
  geom_point()+
  scale_colour_viridis()+
  theme_classic()

## Check the distribution of depth by state

ggplot(dat,aes(x=depth,fill=state))+
  geom_density()+
  scale_x_continuous(breaks=c(0,30,60,100))+
  facet_wrap(~state,ncol = 3,nrow=2)+
  theme_classic()

### Code covariates

dat$slope<-as.numeric(dat$slope)
dat$relief<-as.numeric(dat$relief)
dat$ausbath<-as.numeric(dat$ausbath)
dat$depth<-as.numeric(dat$depth)
dat$t_m<-as.numeric(dat$t_m)
dat$t_sd<-as.numeric(dat$t_sd)
dat$no3_m<-as.numeric(dat$no3_m)
dat$no3_sd<-as.numeric(dat$no3_sd)
dat$po4_m<-as.numeric(dat$po4_m)
dat$po4_sd<-as.numeric(dat$po4_sd)
dat$npp_mean<-as.numeric(dat$npp_mean)
dat$npp_min<-as.numeric(dat$npp_min)
dat$npp_max<-as.numeric(dat$npp_max)
dat$npp_sd<-as.numeric(dat$npp_sd)
dat$dst2ramp<-as.numeric(dat$dst2ramp)
dat$pop.den.50<-as.numeric(dat$pop.den.50)
dat$gravity.50<-as.numeric(dat$gravity.50)
dat$status<-as.factor(dat$status)
dat$campaignid<-as.factor(dat$campaignid)
dat$state<-as.factor(dat$state)
dat$clust.4km<-as.factor(dat$clust.4km)

glimpse(dat)

# Explore samples corresponding to Fished and NTR

table(dat$status)
plot(dat$status)

## Set predictor variables to include in the models
## Note - only predictors identified during exploratory stages - correlations < 0.8 and VIF < 5

names(dat)
pred.vars.cont=c("relief","depth","t_m","t_sd","no3_m","npp_mean","npp_sd","gravity.50")

# Plot of likely transformations 

plot.new()

par(mfrow=c(3,2))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

## Transform predictors 

dat$gravity.50<-log1p(dat$gravity.50)
dat$npp_sd<-log1p(dat$npp_sd)
dat$npp_mean<-log1p(dat$npp_mean)
dat$no3_m<-log1p(dat$no3_m)
dat$relief<-log1p(dat$relief)

## Identify outliers 

plot.new()
par(mfrow=c(1,3))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  boxplot(x)
}

## Explore outliers in NPPmean
## Inshore areas in the Exmouth Guld and Areas of the Pilbara - which high tidal regimes - and high turbidity
## Hence - we retain 

ggplot(dat,aes(x=longitude,y=latitude,size=npp_mean,colour=npp_mean))+
  geom_point()+
  ylim(-25,-20)+
  xlim(110,130)+
  scale_colour_viridis()+
  theme_classic()

## Check Final predictor correlations 

library(corrplot)
library(RColorBrewer)

## Visualize

M<-dat[,pred.vars.cont]
M<-round(cor(M),2)

plot.new()
par(mfrow=c(1,1))
corrplot(M,method="number")

# Check for correlation of predictor variables- remove anything highly correlated (>0.70) 

round(cor(dat[,pred.vars.cont]),2)

## Correlation plots and distributions with ggally

#install.packages('GGally')
library(GGally)

ggpairs(dat,columns = pred.vars.cont, title="correlogram with ggpairs()") 

### Check distributions to fit individual models

plot.new()
par(mfrow=c(1,2))

hist(dat$response)
plot(dat$response)

sum(dat$response == 0 ) / length(dat$response) ## 73% zeroes - We can try to model abundance here

## Quick plots on raw abundances

library(viridis)

ggplot(dat,aes(x=longitude,y=latitude,size=response,colour=response,alpha=response))+
  geom_jitter(width = 1)+
  scale_colour_viridis()+
  ggtitle('(a) Large')+
  geom_point()+
  theme_classic()

## Due to being heavily zero-inflated we will convert the data to presence/absence
## And fit with a binomial distribution with a logit link function

dat<-dat%>%
  mutate(response=ifelse(response>0,1,0))%>%
  glimpse()

# Check to make sure Response vector has not more than 80% zeros---

setwd(model.out)
dir()

unique.vars=unique(as.character(dat$Taxa))
unique.vars.use=character()
for(i in 1:length(unique.vars)){
  temp.dat=dat[which(dat$Taxa==unique.vars[i]),]
  if(length(which(temp.dat$response==0))/nrow(temp.dat)<1){
    unique.vars.use=c(unique.vars.use,unique.vars[i])}
}
unique.vars.use
write.csv(unique.vars.use,file=paste(name,"unique.vars.use.csv",sep = "_"))

# Set variables needed for the FSS function-

pred.vars.cont=pred.vars.cont
resp.vars=unique.vars
use.dat=dat
factor.vars=c("status")
null.vars=c()
out.all=list()
var.imp=list()

# Loop through the FSS function for each Taxa----

for(i in 1:length(resp.vars)){use.dat=dat[which(dat$Taxa==resp.vars[i]),]

Model1=uGamm(response~s(t_m,k=4,bs='cr'),
             family=binomial(link = "logit"), random=~(1|state/clust.4km) + (1|campaignid),
             data=use.dat,
             lme4=TRUE)

model.set=generate.model.set(use.dat=use.dat,
                             test.fit=Model1,
                             pred.vars.cont=pred.vars.cont,
                             pred.vars.fact=factor.vars,
                             #factor.smooth.interactions = factor.vars,
                             #linear.vars = linear.vars,
                             factor.smooth.interactions = list(
                               fact.vars=c("status"),
                               cont.vars=c("gravity.50")),
                             cov.cutoff = 0.28,
                             #cyclic.vars = cyclic.vars,
                             # smooth.smooth.interactions = c("gravity.50","depth","gravity.50"),
                             k=4,
                             max.predictors=4)

out.list=fit.model.set(model.set,
                       max.models=100000,
                       save.model.fits = F,
                       report.unique.r2 = T,
                       parallel=T,
                       VI.mods="all")

names(out.list)

out.list$failed.models # examine the list of failed models
mod.table=out.list$mod.data.out  # look at the model selection table
mod.table=mod.table[order(mod.table$AIC),]
mod.table$cumsum.wi=cumsum(mod.table$wi.AIC)  
out.i=mod.table[which(mod.table$delta.AIC<=2),]
out.all=c(out.all,list(out.i))
var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Either raw importance score
#var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Or importance score weighted by r2

# plot the best models
for(m in 1:nrow(out.i)){
  best.model.name=as.character(out.i$modname[m])
  
  png(file=paste(name,m,resp.vars[i],"mod_fits.png",sep="_"))
  if(best.model.name!="null"){
    par(mfrow=c(3,1),mar=c(9,4,3,1))
    best.model=update(Model1,out.list$success.models[[best.model.name]])
    plot(best.model$gam,all.terms=T,pages=1,residuals=T,pch=16)
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
write.csv(mod.table,file=paste(name,"mod.table.csv",sep="_"))

## (B) Legal vs. Sublegal abundance ----

?full.subsets.gam

# Read in data 

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.assemblage-level.csv")%>%
  filter(relief<200)%>% ## Some problems with few observations close to coast - with extreme values of relief
  rename(response=number)%>%
  filter(!Taxa%in%c('Large'))%>%
  na.omit()%>%
  glimpse()
unique(dat$Taxa)

## Explore spatial variation in gravity values

ggplot(dat,aes(x=longitude,y=latitude,size=log1p(gravity.50),colour=log1p(gravity.50)))+
  geom_point()+
  scale_colour_viridis()+
  theme_classic()

## Check the distribution of depth by state

ggplot(dat,aes(x=depth,fill=state))+
  geom_density()+
  scale_x_continuous(breaks=c(0,30,60,100))+
  facet_wrap(~state,ncol = 3,nrow=2)+
  theme_classic()

### Code covariates

dat$slope<-as.numeric(dat$slope)
dat$relief<-as.numeric(dat$relief)
dat$ausbath<-as.numeric(dat$ausbath)
dat$depth<-as.numeric(dat$depth)
dat$t_m<-as.numeric(dat$t_m)
dat$t_sd<-as.numeric(dat$t_sd)
dat$no3_m<-as.numeric(dat$no3_m)
dat$no3_sd<-as.numeric(dat$no3_sd)
dat$po4_m<-as.numeric(dat$po4_m)
dat$po4_sd<-as.numeric(dat$po4_sd)
dat$npp_mean<-as.numeric(dat$npp_mean)
dat$npp_min<-as.numeric(dat$npp_min)
dat$npp_max<-as.numeric(dat$npp_max)
dat$npp_sd<-as.numeric(dat$npp_sd)
dat$dst2ramp<-as.numeric(dat$dst2ramp)
dat$pop.den.50<-as.numeric(dat$pop.den.50)
dat$gravity.50<-as.numeric(dat$gravity.50)
dat$status<-as.factor(dat$status)
dat$campaignid<-as.factor(dat$campaignid)
dat$state<-as.factor(dat$state)
dat$clust.4km<-as.factor(dat$clust.4km)

glimpse(dat)

# Explore samples corresponding to Fished and NTR

table(dat$status)
plot(dat$status)

## Set predictor variables to include in the models
## Note - only predictors identified during exploratory stages - correlations < 0.8 and VIF < 5

names(dat)
pred.vars.cont=c("relief","depth","t_m","t_sd","no3_m","npp_mean","npp_sd","gravity.50")

# Plot of likely transformations 

plot.new()

par(mfrow=c(3,2))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

## Transform predictors 

dat$gravity.50<-log1p(dat$gravity.50)
dat$npp_sd<-log1p(dat$npp_sd)
dat$npp_mean<-log1p(dat$npp_mean)
dat$no3_m<-log1p(dat$no3_m)
dat$relief<-log1p(dat$relief)

## Identify outliers 

plot.new()
par(mfrow=c(1,3))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  boxplot(x)
}

## Explore outliers in NPPmean
## Inshore areas in the Exmouth Guld and Areas of the Pilbara - which high tidal regimes - and high turbidity
## Hence - we retain 

ggplot(dat,aes(x=longitude,y=latitude,size=npp_mean,colour=npp_mean))+
  geom_point()+
  ylim(-25,-20)+
  xlim(110,130)+
  scale_colour_viridis()+
  theme_classic()

## Check Final predictor correlations 

library(corrplot)
library(RColorBrewer)

## Visualize

M<-dat[,pred.vars.cont]
M<-round(cor(M),2)

plot.new()
par(mfrow=c(1,1))
corrplot(M,method="number")

# Check for correlation of predictor variables- remove anything highly correlated (>0.70) 

round(cor(dat[,pred.vars.cont]),2)

## Correlation plots and distributions with ggally

#install.packages('GGally')
library(GGally)

ggpairs(dat,columns = pred.vars.cont, title="correlogram with ggpairs()") 

### Check distributions to fit individual models

plot.new()
par(mfrow=c(1,2))

## Legal

legal<-dat%>%
  dplyr::filter(Taxa=="Legal")%>%
  glimpse()
hist(legal$response)
plot(legal$response)

sum(legal$response == 0 ) / length(legal$response) ## 46% zeroes

## Sub-legal

sublegal<-dat%>%
  dplyr::filter(Taxa=="Sub-legal")%>%
  glimpse()
hist(sublegal$response)
plot(sublegal$response)

sum(sublegal$response == 0 ) / length(sublegal$response) ## 37% zeroes

## Quick plots on raw abundances

library(viridis)

legal.plot<-ggplot(legal,aes(x=longitude,y=latitude,size=response,colour=response,alpha=response))+
  geom_jitter(width = 1)+
  scale_colour_viridis()+
  ggtitle('(a) Legal')+
  geom_point()+
  theme_classic()
legal.plot

sublegal.plot<-ggplot(sublegal,aes(x=longitude,y=latitude,size=response,colour=response,alpha=response))+
  geom_jitter(width = 1)+
  scale_colour_viridis()+
  ggtitle('(b) Sub-legal')+
  geom_point()+
  theme_classic()
sublegal.plot

## Arrange and export

ggarrange(legal.plot,sublegal.plot,ncol=2,align = 'hv',common.legend = T,legend = "bottom")

# Check to make sure Response vector has not more than 80% zeros---

setwd(model.out)
dir()

unique.vars=unique(as.character(dat$Taxa))
unique.vars.use=character()
for(i in 1:length(unique.vars)){
  temp.dat=dat[which(dat$Taxa==unique.vars[i]),]
  if(length(which(temp.dat$response==0))/nrow(temp.dat)<1){
    unique.vars.use=c(unique.vars.use,unique.vars[i])}
}
unique.vars.use
write.csv(unique.vars.use,file=paste(name,"unique.vars.use.csv",sep = "_"))

# Set variables needed for the FSS function-

pred.vars.cont=pred.vars.cont
resp.vars=unique.vars
use.dat=dat
factor.vars=c("status")
null.vars=c()
out.all=list()
var.imp=list()

# Loop through the FSS function for each Taxa----

for(i in 1:length(resp.vars)){use.dat=dat[which(dat$Taxa==resp.vars[i]),]

Model1=uGamm(response~s(t_m,k=4,bs='cr'),
             family=negbin(1), random=~(1|state/clust.4km) + (1|campaignid),
             data=use.dat,
             lme4=TRUE)

model.set=generate.model.set(use.dat=use.dat,
                             test.fit=Model1,
                             pred.vars.cont=pred.vars.cont,
                             pred.vars.fact=factor.vars,
                             #factor.smooth.interactions = factor.vars,
                             #linear.vars = linear.vars,
                             factor.smooth.interactions = list(
                               fact.vars=c("status"),
                               cont.vars=c("gravity.50")),
                             cov.cutoff = 0.28,
                             #cyclic.vars = cyclic.vars,
                             # smooth.smooth.interactions = c("gravity.50","depth","gravity.50"),
                             k=4,
                             max.predictors=4)

out.list=fit.model.set(model.set,
                       max.models=100000,
                       save.model.fits = F,
                       report.unique.r2 = T,
                       parallel=T,
                       VI.mods="all")

names(out.list)

out.list$failed.models # examine the list of failed models
mod.table=out.list$mod.data.out  # look at the model selection table
mod.table=mod.table[order(mod.table$AIC),]
mod.table$cumsum.wi=cumsum(mod.table$wi.AIC)  
out.i=mod.table[which(mod.table$delta.AIC<=2),]
out.all=c(out.all,list(out.i))
var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Either raw importance score
#var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Or importance score weighted by r2

# plot the best models
for(m in 1:nrow(out.i)){
  best.model.name=as.character(out.i$modname[m])
  
  png(file=paste(name,m,resp.vars[i],"mod_fits.png",sep="_"))
  if(best.model.name!="null"){
    par(mfrow=c(3,1),mar=c(9,4,3,1))
    best.model=update(Model1,out.list$success.models[[best.model.name]])
    plot(best.model$gam,all.terms=T,pages=1,residuals=T,pch=16)
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
write.csv(mod.table,file=paste(name,"mod.table.csv",sep="_"))

## (B) Legal vs. Sublegal biomass ----

?full.subsets.gam

# Read in data 

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.assemblage-level.csv")%>%
  filter(relief<200)%>% ## Some problems with few observations close to coast - with extreme values of relief
  rename(response=biomass)%>%
  filter(!Taxa%in%c('Large'))%>%
  na.omit()%>%
  glimpse()
unique(dat$Taxa)

## Explore spatial variation in gravity values

ggplot(dat,aes(x=longitude,y=latitude,size=log1p(gravity.50),colour=log1p(gravity.50)))+
  geom_point()+
  scale_colour_viridis()+
  theme_classic()

## Check the distribution of depth by state

ggplot(dat,aes(x=depth,fill=state))+
  geom_density()+
  scale_x_continuous(breaks=c(0,30,60,100))+
  facet_wrap(~state,ncol = 3,nrow=2)+
  theme_classic()

### Code covariates

dat$slope<-as.numeric(dat$slope)
dat$relief<-as.numeric(dat$relief)
dat$ausbath<-as.numeric(dat$ausbath)
dat$depth<-as.numeric(dat$depth)
dat$t_m<-as.numeric(dat$t_m)
dat$t_sd<-as.numeric(dat$t_sd)
dat$no3_m<-as.numeric(dat$no3_m)
dat$no3_sd<-as.numeric(dat$no3_sd)
dat$po4_m<-as.numeric(dat$po4_m)
dat$po4_sd<-as.numeric(dat$po4_sd)
dat$npp_mean<-as.numeric(dat$npp_mean)
dat$npp_min<-as.numeric(dat$npp_min)
dat$npp_max<-as.numeric(dat$npp_max)
dat$npp_sd<-as.numeric(dat$npp_sd)
dat$dst2ramp<-as.numeric(dat$dst2ramp)
dat$pop.den.50<-as.numeric(dat$pop.den.50)
dat$gravity.50<-as.numeric(dat$gravity.50)
dat$status<-as.factor(dat$status)
dat$campaignid<-as.factor(dat$campaignid)
dat$state<-as.factor(dat$state)
dat$clust.4km<-as.factor(dat$clust.4km)

glimpse(dat)

# Explore samples corresponding to Fished and NTR

table(dat$status)
plot(dat$status)

## Set predictor variables to include in the models
## Note - only predictors identified during exploratory stages - correlations < 0.8 and VIF < 5

names(dat)
pred.vars.cont=c("relief","depth","t_m","t_sd","no3_m","npp_mean","npp_sd","gravity.50")

# Plot of likely transformations 

plot.new()

par(mfrow=c(3,2))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

## Transform predictors 

dat$gravity.50<-log1p(dat$gravity.50)
dat$npp_sd<-log1p(dat$npp_sd)
dat$npp_mean<-log1p(dat$npp_mean)
dat$no3_m<-log1p(dat$no3_m)
dat$relief<-log1p(dat$relief)

## Identify outliers 

plot.new()
par(mfrow=c(1,3))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  boxplot(x)
}

## Explore outliers in NPPmean
## Inshore areas in the Exmouth Guld and Areas of the Pilbara - which high tidal regimes - and high turbidity
## Hence - we retain 

ggplot(dat,aes(x=longitude,y=latitude,size=npp_mean,colour=npp_mean))+
  geom_point()+
  ylim(-25,-20)+
  xlim(110,130)+
  scale_colour_viridis()+
  theme_classic()

## Check Final predictor correlations 

library(corrplot)
library(RColorBrewer)

## Visualize

M<-dat[,pred.vars.cont]
M<-round(cor(M),2)

plot.new()
par(mfrow=c(1,1))
corrplot(M,method="number")

# Check for correlation of predictor variables- remove anything highly correlated (>0.70) 

round(cor(dat[,pred.vars.cont]),2)

## Correlation plots and distributions with ggally

#install.packages('GGally')
library(GGally)

ggpairs(dat,columns = pred.vars.cont, title="correlogram with ggpairs()") 

### Check distributions to fit individual models

plot.new()
par(mfrow=c(1,2))

## Legal

legal<-dat%>%
  dplyr::filter(Taxa=="Legal")%>%
  glimpse()
hist(legal$response)
plot(legal$response)

sum(legal$response == 0 ) / length(legal$response) ## 48% zeroes

## Sub-legal

sublegal<-dat%>%
  dplyr::filter(Taxa=="Sub-legal")%>%
  glimpse()
hist(sublegal$response)
plot(sublegal$response)

sum(sublegal$response == 0 ) / length(sublegal$response) ## 38% zeroes

## Quick plots on raw abundances

library(viridis)

legal.plot<-ggplot(legal,aes(x=longitude,y=latitude,size=response,colour=response,alpha=response))+
  geom_jitter(width = 1)+
  scale_colour_viridis()+
  ggtitle('(a) Legal')+
  geom_point()+
  theme_classic()
legal.plot

sublegal.plot<-ggplot(sublegal,aes(x=longitude,y=latitude,size=response,colour=response,alpha=response))+
  geom_jitter(width = 1)+
  scale_colour_viridis()+
  ggtitle('(b) Sub-legal')+
  geom_point()+
  theme_classic()
sublegal.plot

## Arrange and export

library(ggpubr)

ggarrange(legal.plot,sublegal.plot,ncol=2,align = 'hv',common.legend = T,legend = "bottom")

# Check to make sure Response vector has not more than 80% zeros---

setwd(model.out)
dir()

unique.vars=unique(as.character(dat$Taxa))
unique.vars.use=character()
for(i in 1:length(unique.vars)){
  temp.dat=dat[which(dat$Taxa==unique.vars[i]),]
  if(length(which(temp.dat$response==0))/nrow(temp.dat)<1){
    unique.vars.use=c(unique.vars.use,unique.vars[i])}
}
unique.vars.use
write.csv(unique.vars.use,file=paste(name,"unique.vars.use.csv",sep = "_"))

# Set variables needed for the FSS function-

pred.vars.cont=pred.vars.cont
resp.vars=unique.vars
use.dat=dat
factor.vars=c("status")
null.vars=c()
out.all=list()
var.imp=list()

# Loop through the FSS function for each Taxa----

for(i in 1:length(resp.vars)){use.dat=dat[which(dat$Taxa==resp.vars[i]),]

Model1=uGamm(response~s(t_m,k=4,bs='cr'),
             family=Tweedie(1.2), random=~(1|state/clust.4km) + (1|campaignid),
             data=use.dat,
             lme4=TRUE)

model.set=generate.model.set(use.dat=use.dat,
                             test.fit=Model1,
                             pred.vars.cont=pred.vars.cont,
                             pred.vars.fact=factor.vars,
                             #factor.smooth.interactions = factor.vars,
                             #linear.vars = linear.vars,
                             factor.smooth.interactions = list(
                               fact.vars=c("status"),
                               cont.vars=c("gravity.50")),
                             cov.cutoff = 0.28,
                             #cyclic.vars = cyclic.vars,
                             # smooth.smooth.interactions = c("gravity.50","depth","gravity.50"),
                             k=4,
                             max.predictors=4)

out.list=fit.model.set(model.set,
                       max.models=100000,
                       save.model.fits = F,
                       report.unique.r2 = T,
                       parallel=T,
                       VI.mods="all")

names(out.list)

out.list$failed.models # examine the list of failed models
mod.table=out.list$mod.data.out  # look at the model selection table
mod.table=mod.table[order(mod.table$AIC),]
mod.table$cumsum.wi=cumsum(mod.table$wi.AIC)  
out.i=mod.table[which(mod.table$delta.AIC<=2),]
out.all=c(out.all,list(out.i))
var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Either raw importance score
#var.imp=c(var.imp,list(out.list$variable.importance$aic$variable.weights.raw)) #Or importance score weighted by r2

# plot the best models
for(m in 1:nrow(out.i)){
  best.model.name=as.character(out.i$modname[m])
  
  png(file=paste(name,m,resp.vars[i],"mod_fits.png",sep="_"))
  if(best.model.name!="null"){
    par(mfrow=c(3,1),mar=c(9,4,3,1))
    best.model=update(Model1,out.list$success.models[[best.model.name]])
    plot(best.model$gam,all.terms=T,pages=1,residuals=T,pch=16)
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
write.csv(mod.table,file=paste(name,"mod.table.csv",sep="_"))

################################################################################################
########################### (4) Plot of top models GAMM ########################################

## (A) Large Presence-absence ----

?full.subsets.gam

# Read in data 

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.assemblage-level.csv")%>%
  filter(relief<200)%>% ## Some problems with few observations close to coast - with extreme values of relief
  rename(response=number)%>%
  filter(Taxa%in%c('Large'))%>%
  na.omit()%>%
  glimpse()
unique(dat$Taxa)

### Code covariates

dat$slope<-as.numeric(dat$slope)
dat$relief<-as.numeric(dat$relief)
dat$ausbath<-as.numeric(dat$ausbath)
dat$depth<-as.numeric(dat$depth)
dat$t_m<-as.numeric(dat$t_m)
dat$t_sd<-as.numeric(dat$t_sd)
dat$no3_m<-as.numeric(dat$no3_m)
dat$no3_sd<-as.numeric(dat$no3_sd)
dat$po4_m<-as.numeric(dat$po4_m)
dat$po4_sd<-as.numeric(dat$po4_sd)
dat$npp_mean<-as.numeric(dat$npp_mean)
dat$npp_min<-as.numeric(dat$npp_min)
dat$npp_max<-as.numeric(dat$npp_max)
dat$npp_sd<-as.numeric(dat$npp_sd)
dat$dst2ramp<-as.numeric(dat$dst2ramp)
dat$pop.den.50<-as.numeric(dat$pop.den.50)
dat$gravity.50<-as.numeric(dat$gravity.50)
dat$status<-as.factor(dat$status)
dat$campaignid<-as.factor(dat$campaignid)
dat$state<-as.factor(dat$state)
dat$clust.4km<-as.factor(dat$clust.4km)

glimpse(dat)

# Explore samples corresponding to Fished and NTR

table(dat$status)
plot(dat$status)

## Set predictor variables to include in the models
## Note - only predictors identified during exploratory stages - correlations < 0.8 and VIF < 5

names(dat)
pred.vars.cont=c("relief","depth","t_m","t_sd","no3_m","npp_mean","npp_sd","gravity.50")

## Transform predictors 

dat$gravity.50<-log1p(dat$gravity.50)
dat$npp_sd<-log1p(dat$npp_sd)
dat$npp_mean<-log1p(dat$npp_mean)
dat$no3_m<-log1p(dat$no3_m)
dat$relief<-log1p(dat$relief)

### Check distributions to fit individual models

plot.new()
par(mfrow=c(1,2))

hist(dat$response)
plot(dat$response)

sum(dat$response == 0 ) / length(dat$response) ## 73% zeroes - Presence/absence model

## Convert to presence/absence

dat<-dat%>%
  mutate(response=ifelse(response>0,1,0))%>%
  glimpse()

## Top-ranked 

### M1 - depth + gravity.50 + relief ----

M1=gamm4(response ~ s(depth, k = 4, bs = "cr") + s(gravity.50, k = 4, bs = "cr") + 
         s(relief, k = 4, bs = "cr") + status,
         family=binomial(link = "logit"), random=~(1|state/clust.4km) + (1|campaignid),
         data=dat)
plot(M1$gam,all.terms = TRUE,pages=1,scale=0,residuals=F)
summary(M1$gam)
summary(M1$mer)
AICc(M1$mer)
BIC(M1$mer)

cor(fitted(M1$mer), dat$response)^2 #R^2 for lmer model
cor(fitted(M1$gam), dat$response)^2 #R^2 for gam model

## Plot of top-ranked

library(sjPlot)

## Legal ----

gravity.large<-plot_model(M1,type ="pred",terms="gravity.50",alpha=0.3)  + 
  ggtitle('(c)')  + theme_classic() + Theme1 + xlab("Gravity") +  ylab('Probability of occurrence (%)') + geom_rug(sides = "b")
gravity.large

depth.large<-plot_model(M1,type ="pred",terms="depth",alpha=0.3)  + xlab ("Depth") +
  ggtitle('(d)')+ theme_classic() + Theme1 + ylab("") +  geom_rug(sides = "b")
depth.large

relief.large<-plot_model(M1,type ="pred",terms="relief",alpha=0.3)  + 
  ggtitle('(e)')  + theme_classic() + Theme1 +  xlab ("Relief") + geom_rug(sides = "b")
relief.large

## Arrange all plots in a grob

ggarrange(relief.large,depth.large,gravity.large,nrow = 1,ncol = 4,align = 'hv',common.legend = T,legend = "bottom")

## Plot the effect of management status manually

## Legal ---

unique(dat$state)
unique(dat$campaignid)
unique(dat$clust.4km)
unique(dat$status)

testdata <- expand.grid(depth=mean(dat$depth),
                        relief=mean(dat$relief),
                        gravity.50=mean(dat$gravity.50),
                        state='WA.Temperate Australasia',
                        campaignid='2016-10_North.Kimberley_stereoBRUVs',
                        clust.4km='1',
                        status=c("Fished","No-take"))%>%
  glimpse()

fits <- predict(M1$gam, newdata=testdata, type='response', se.fit=T,exclude = c("s(state)","s(campaignid)","s(clust.4km)"))

La.predicts = testdata%>%data.frame(fits)%>%
  group_by(status)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()%>%
  glimpse()
write.csv(La.predicts,"predicts.csv") #there is some BUG in dplyr - that this fixes
La.predicts<-read.csv("predicts.csv")%>%
  glimpse()

Status.large<- ggplot(aes(x=status,y=response,fill=status), data=La.predicts) +
  ylab("")+
  xlab('Status')+
  #ylim(0,6)+
  scale_fill_manual(labels = c("Fished", "No-take"),values=c("brown1", "blue3"))+
  geom_bar(stat = "identity",show.legend=TRUE,colour="black")+
  geom_errorbar(aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5,show.legend=FALSE) +
  scale_x_discrete(expand = c(0.3,0.3))+
  scale_y_continuous(expand=c(0,0))+
  theme_classic()+
  Theme1
Status.large

### Check residuals of models selected for prediction

## Extract residuals from our modes

res1<-resid(M1$mer, type = "deviance")

## Plot histogram of residuals

plot.new()
par(mfrow=c(1,1))

hist(res1,main = '(a)')

## Residuals vs. fitted values ----

## Large

setwd(plots.dir)
dir()

pdf("Large_Assemblage.pdf",width = 3,height = 3,useDingbats=FALSE) 

F1<-fitted(M1$mer)

plot(x=F1,
     y=res1,
     xlab="Fitted values",
     ylab="Deviance residuals",
     cex.lab=1.5)
abline(h=0,lty=2)

dev.off()

## Put the residuals in the dataset

dat$res<-res1

### Check for residual auto-correlation ----

## Convert into an spatial object

coordinates(dat)<-~longitude+latitude
class(dat)

### Now we can set the reference system to the widely used WGS84

proj4string(dat)<-CRS("+init=epsg:4326")
dat@proj4string

proj4string(dat)<-CRS("+init=epsg:4326")
dat@proj4string

## Compute semivariogram

library(gstat)

dat.response<-variogram(dat$response~1,data=dat,cutoff=20,width=0.1)
dat.res<-variogram(dat$res~1,data=dat,cutoff=20,width=0.1)

## Export semivariograms

## Large response

pdf("Global_Legal_Semivariogram_response.pdf",useDingbats=FALSE) 

plot(dat.response,main="(a) Large",ylab="Semivariance",xlab="Distance (km)")

dev.off()

## Large residuals 

pdf("Global_Legal_Semivariogram_residuals.pdf",useDingbats=FALSE) 

plot(dat.res,main="(b) Large",ylab="Semivariance",xlab="Distance (km)")

dev.off()

### Calculate Moran's I based on distance via the spdep package
### We would calculated manually for each basin, diversity index combination and store it in a csv file
### However - for future reference it would be good to loop all of this together

## Global Moran's I ---- 

library(geosphere)

## Legal

dists <- distm(dat, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(dat$response,dists.inv)

## Legal (res 1)

Moran.I(dat$res,dists.inv) 

## (B) Legal vs. sub-legal Abundance ----

## Bring in the data 

# Read in data 

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.assemblage-level.csv")%>%
  filter(relief<200)%>% ## Some problems with few observations close to coast - with extreme values of relief
  rename(response=number)%>%
  filter(!Taxa%in%c('Large'))%>%
  na.omit()%>%
  glimpse()
unique(dat$Taxa)

### Code covariates

dat$slope<-as.numeric(dat$slope)
dat$relief<-as.numeric(dat$relief)
dat$ausbath<-as.numeric(dat$ausbath)
dat$depth<-as.numeric(dat$depth)
dat$t_m<-as.numeric(dat$t_m)
dat$t_sd<-as.numeric(dat$t_sd)
dat$no3_m<-as.numeric(dat$no3_m)
dat$no3_sd<-as.numeric(dat$no3_sd)
dat$po4_m<-as.numeric(dat$po4_m)
dat$po4_sd<-as.numeric(dat$po4_sd)
dat$npp_mean<-as.numeric(dat$npp_mean)
dat$npp_min<-as.numeric(dat$npp_min)
dat$npp_max<-as.numeric(dat$npp_max)
dat$npp_sd<-as.numeric(dat$npp_sd)
dat$dst2ramp<-as.numeric(dat$dst2ramp)
dat$pop.den.50<-as.numeric(dat$pop.den.50)
dat$gravity.50<-as.numeric(dat$gravity.50)
dat$status<-as.factor(dat$status)
dat$campaignid<-as.factor(dat$campaignid)
dat$state<-as.factor(dat$state)
dat$clust.4km<-as.factor(dat$clust.4km)

glimpse(dat)

## Transform predictors 

dat$gravity.50<-log1p(dat$gravity.50)
dat$npp_sd<-log1p(dat$npp_sd)
dat$npp_mean<-log1p(dat$npp_mean)
dat$no3_m<-log1p(dat$no3_m)
dat$relief<-log1p(dat$relief)

### Check distributions to fit individual models

plot.new()
par(mfrow=c(1,2))

## Legal

legal<-dat%>%
  dplyr::filter(Taxa=="Legal")%>%
  glimpse()
hist(legal$response)
plot(legal$response)

sum(legal$response == 0 ) / length(legal$response) ## 48% zeroes

## Sub-legal

sublegal<-dat%>%
  dplyr::filter(Taxa=="Sub-legal")%>%
  glimpse()
hist(sublegal$response)
plot(sublegal$response)

sum(sublegal$response == 0 ) / length(sublegal$response) ## 38% zeroes

## Top-ranked 

## (A) Legal

### M1 - depth + gravity + relief + status ----

M1=gamm4(response ~ s(depth, k = 4, bs = "cr") + s(relief, k = 4, bs = "cr") + 
           s(gravity.50, k = 4, bs = "cr") + status,
         family=negbin(1), random=~(1|state/clust.4km) + (1|campaignid),
         data=legal)
plot(M1$gam,all.terms = TRUE,pages=1,scale=0,residuals=F)
summary(M1$gam)
summary(M1$mer)
AICc(M1$mer)
BIC(M1$mer)

cor(fitted(M1$mer), legal$response)^2 #R^2 for lmer model
cor(fitted(M1$gam), legal$response)^2 #R^2 for gam model

### Sub-legal 
## M8 - t_m + t_sd + status ----

M2=gamm4(response ~ s(t_m, k = 4, bs = "cr") + s(t_sd, k = 4, bs = "cr") + status,
         family=negbin(1), random=~(1|state/clust.4km) + (1|campaignid),
         data=sublegal)
plot(M2$gam,all.terms = TRUE,pages=1,scale=0,residuals=F)
summary(M2$gam)
summary(M2$mer)
AICc(M2$mer)
BIC(M2$mer)

cor(fitted(M2$mer), sublegal$response)^2 #R^2 for lmer model
cor(fitted(M2$gam), sublegal$response)^2 #R^2 for gam model

## Plot of top-ranked

library(sjPlot)

## Legal ----

gravity.legal<-plot_model(M1,type ="pred",terms="gravity.50",alpha=0.3)  + 
  ggtitle('(c)')  + theme_classic() + Theme1 + xlab("Gravity") +  geom_rug(sides = "b")
gravity.legal

depth.legal<-plot_model(M1,type ="pred",terms="depth",alpha=0.3)  + xlab ("Depth") +
  ggtitle('(d)')+ theme_classic() + Theme1 + ylab("") +  geom_rug(sides = "b")
depth.legal

relief.legal<-plot_model(M1,type ="pred",terms="relief",alpha=0.3)  + 
  ggtitle('(e)')  + theme_classic() + Theme1 +  ylab('Relative abundance (MaxN)') + 
  xlab ("Relief") + geom_rug(sides = "b")
relief.legal

## Arrange all plots in a grob

ggarrange(relief.legal,depth.legal,gravity.legal,nrow = 1,ncol = 4,align = 'hv',common.legend = T,legend = "bottom")

## Sub-legal ----

SSTmean.sublegal<-plot_model(M2,type ="pred",terms="t_m",alpha=0.3)  + xlab ("SSTmean") +
  ggtitle('(g)')+ theme_classic() + Theme1 + ylab("Relative abundance (MaxN)") + geom_rug(sides = "b")
SSTmean.sublegal

SSTsd.sublegal<-plot_model(M2,type ="pred",terms="t_sd",alpha=0.3)  + 
  ggtitle('(h)')  + theme_classic() + Theme1 + xlab("SSTsd") + ylab("") + geom_rug(sides = "b")
SSTsd.sublegal

## Plot the effect of management status manually

## Legal ---

unique(dat$state)
unique(dat$campaignid)
unique(dat$clust.4km)
unique(dat$status)

testdata <- expand.grid(depth=mean(legal$depth),
                        relief=mean(legal$relief),
                        gravity.50=mean(legal$gravity.50),
                        state='WA.Temperate Australasia',
                        campaignid='2016-10_North.Kimberley_stereoBRUVs',
                        clust.4km='1',
                        status=c("Fished","No-take"))%>%
  glimpse()

fits <- predict(M1$gam, newdata=testdata, type='response', se.fit=T,exclude = c("s(state)","s(campaignid)","s(clust.4km)"))

La.predicts = testdata%>%data.frame(fits)%>%
  group_by(status)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()%>%
  glimpse()
write.csv(La.predicts,"predicts.csv") #there is some BUG in dplyr - that this fixes
La.predicts<-read.csv("predicts.csv")%>%
  glimpse()

Status.legal<- ggplot(aes(x=status,y=response,fill=status), data=La.predicts) +
  ylab("")+
  xlab('Status')+
  #ylim(0,6)+
  scale_fill_manual(labels = c("Fished", "No-take"),values=c("brown1", "blue3"))+
  geom_bar(stat = "identity",show.legend=TRUE,colour="black")+
  geom_errorbar(aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5,show.legend=FALSE) +
  scale_x_discrete(expand = c(0.3,0.3))+
  scale_y_continuous(expand=c(0,0))+
  theme_classic()+
  Theme1
Status.legal

## Sub-legal ---

unique(dat$state)
unique(dat$campaignid)
unique(dat$clust.4km)
unique(dat$status)

testdata <- expand.grid(t_m=mean(sublegal$t_m),
                        t_sd=mean(sublegal$t_sd),
                        state='WA.Temperate Australasia',
                        campaignid='2016-10_North.Kimberley_stereoBRUVs',
                        clust.4km='1',
                        status=c("Fished","No-take"))%>%
  glimpse()

fits <- predict(M2$gam, newdata=testdata, type='response', se.fit=T,exclude = c("s(state)","s(campaignid)","s(clust.4km)"))

La.predicts = testdata%>%data.frame(fits)%>%
  group_by(status)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()%>%
  glimpse()
write.csv(La.predicts,"predicts.csv") #there is some BUG in dplyr - that this fixes
La.predicts<-read.csv("predicts.csv")%>%
  glimpse()

Status.sublegal<- ggplot(aes(x=status,y=response,fill=status), data=La.predicts) +
  ylab("")+
  xlab('Status')+
  #ylim(0,6)+
  scale_fill_manual(labels = c("Fished", "No-take"),values=c("brown1", "blue3"))+
  geom_bar(stat = "identity",show.legend=TRUE,colour="black")+
  geom_errorbar(aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5,show.legend=FALSE) +
  scale_x_discrete(expand = c(0.3,0.3))+
  scale_y_continuous(expand=c(0,0))+
  theme_classic()+
  Theme1
Status.sublegal

## Arrange all plots in a grob

ggarrange(Status.large,gravity.large,depth.large,relief.large,
          Status.legal,gravity.legal,depth.legal,relief.legal,
          Status.sublegal,SSTmean.sublegal,SSTsd.sublegal,
          nrow = 3,ncol = 4,align = 'hv',common.legend = T,legend = "none")

## Export plot

setwd("C:/Users/22373243/Dropbox/Projects/Analysis/Analysis_GlobalArchive_Bodysize_Nestor/Plots") 

name<-'Top_ranked_20210401'

ggsave(paste(name,".png",sep="."),width = 19.5, height = 15,units = "cm",dpi = 600)
ggsave(paste(name,".pdf",sep="."),width = 19.5, height = 15,units = "cm",dpi = 600,useDingbats=FALSE)

### Check residuals of models selected for prediction

## Extract residuals from our modes

res1<-resid(M1$mer, type = "deviance")
res2<-resid(M2$mer, type = "deviance")

## Plot histogram of residuals

plot.new()
par(mfrow=c(2,2))

hist(res1,main = '(a)')
hist(res2,main='(b)')

## Residuals vs. fitted values ----

## Legal

setwd(plots.dir)
dir()

plot.new()
par(mfrow=c(1,1))

pdf("Legal_Assemblage.pdf",width = 3,height = 3,useDingbats=FALSE) 

F1<-fitted(M1$mer)

plot(x=F1,
     y=res1,
     xlab="Fitted values",
     ylab="Deviance residuals",
     cex.lab=1.5)
abline(h=0,lty=2)

dev.off()


## Sub-legal

plot.new()
par(mfrow=c(1,1))

pdf("Sub-legal_Assemblage.pdf",width = 3,height = 3,useDingbats=FALSE) 

F2<-fitted(M2$mer)

plot(x=F2,
     y=res2,
     xlab="Fitted values",
     ylab="Deviance residuals",
     cex.lab=1.5)
abline(h=0,lty=2)


dev.off()

## Put the residuals in the dataset

legal$res<-res1
sublegal$res<-res2

### Check for residual auto-correlation ----

## Convert into an spatial object

coordinates(legal)<-~longitude+latitude
class(legal)

coordinates(sublegal)<-~longitude+latitude
class(sublegal)

### Now we can set the reference system to the widely used WGS84

proj4string(legal)<-CRS("+init=epsg:4326")
legal@proj4string

proj4string(sublegal)<-CRS("+init=epsg:4326")
sublegal@proj4string

## Compute semivariogram

library(gstat)

legal.var.response<-variogram(legal$response~1,data=legal,cutoff=20,width=0.1)
sublegal.var.response<-variogram(sublegal$response~1,data=legal,cutoff=20,width=0.1)
legal.var.res<-variogram(legal$res~1,data=legal,cutoff=20,width=0.1)
sublegal.var.res<-variogram(sublegal$res~1,data=legal,cutoff=20,width=0.1)

## Export semivariograms

## Legal response

pdf("Global_Legal_Semivariogram_response.pdf",useDingbats=FALSE) 

plot(legal.var.response,main="(a) Legal (response)",ylab="Semivariance",xlab="Distance (km)")

dev.off()

## Legal residuals 

pdf("Global_Legal_Semivariogram_residuals.pdf",useDingbats=FALSE) 

plot(legal.var.res,main="(b) Legal (residuals)",ylab="Semivariance",xlab="Distance (km)")

dev.off()

## Sub-legal response

pdf("Global_Sublegal_Semivariogram_response.pdf",useDingbats=FALSE) 

plot(sublegal.var.response,main="(c) Sub-legal (residuals-nested)",ylab="Semivariance",xlab="Distance (km)")

dev.off()

## Legal residual (3)

pdf("Global_Sublegal_Semivariogram_residuals.pdf",useDingbats=FALSE) 

plot(sublegal.var.res,main="(c) Sub-legal (residuals)",ylab="Semivariance",xlab="Distance (km)",ylim=c(0.5,1.2))

dev.off()

### Calculate Moran's I based on distance via the spdep package
### We would calculated manually for each basin, diversity index combination and store it in a csv file
### However - for future reference it would be good to loop all of this together

## Global Moran's I ---- 

library(geosphere)

## Legal

dists <- distm(legal, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(legal$response,dists.inv)

## Legal (res 1)

Moran.I(legal$res,dists.inv) 

## Sub-legal

dists <- distm(sublegal, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(sublegal$response,dists.inv)

## Legal (res)

Moran.I(sublegal$res,dists.inv) 

## (C) Legal vs. Sub-legal Biomass ----

## Bring in the data 

# Read in data 

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.assemblage-level.csv")%>%
  filter(relief<200)%>% ## Some problems with few observations close to coast - with extreme values of relief
  rename(response=biomass)%>%
  filter(!Taxa%in%c('Large'))%>%
  na.omit()%>%
  glimpse()
unique(dat$Taxa)

### Code covariates

dat$slope<-as.numeric(dat$slope)
dat$relief<-as.numeric(dat$relief)
dat$ausbath<-as.numeric(dat$ausbath)
dat$depth<-as.numeric(dat$depth)
dat$t_m<-as.numeric(dat$t_m)
dat$t_sd<-as.numeric(dat$t_sd)
dat$no3_m<-as.numeric(dat$no3_m)
dat$no3_sd<-as.numeric(dat$no3_sd)
dat$po4_m<-as.numeric(dat$po4_m)
dat$po4_sd<-as.numeric(dat$po4_sd)
dat$npp_mean<-as.numeric(dat$npp_mean)
dat$npp_min<-as.numeric(dat$npp_min)
dat$npp_max<-as.numeric(dat$npp_max)
dat$npp_sd<-as.numeric(dat$npp_sd)
dat$dst2ramp<-as.numeric(dat$dst2ramp)
dat$pop.den.50<-as.numeric(dat$pop.den.50)
dat$gravity.50<-as.numeric(dat$gravity.50)
dat$status<-as.factor(dat$status)
dat$campaignid<-as.factor(dat$campaignid)
dat$state<-as.factor(dat$state)
dat$clust.4km<-as.factor(dat$clust.4km)

glimpse(dat)

## Transform predictors 

dat$gravity.50<-log1p(dat$gravity.50)
dat$npp_sd<-log1p(dat$npp_sd)
dat$npp_mean<-log1p(dat$npp_mean)
dat$no3_m<-log1p(dat$no3_m)
dat$relief<-log1p(dat$relief)

### Check distributions to fit individual models

plot.new()
par(mfrow=c(1,2))

## Legal

legal<-dat%>%
  dplyr::filter(Taxa=="Legal")%>%
  glimpse()
hist(legal$response)
plot(legal$response)

sum(legal$response == 0 ) / length(legal$response) ## 48% zeroes

## Sub-legal

sublegal<-dat%>%
  dplyr::filter(Taxa=="Sub-legal")%>%
  glimpse()
hist(sublegal$response)
plot(sublegal$response)

sum(sublegal$response == 0 ) / length(sublegal$response) ## 38% zeroes

## Top-ranked 

## (A) Legal

M1=gamm4(response ~ s(depth, k = 4, bs = "cr") + s(relief, k = 4, bs = "cr") + 
           s(gravity.50, k = 4, bs = "cr") + status,
         family=Tweedie(1.2), random=~(1|state/clust.4km) + (1|campaignid),
         data=legal)
plot(M1$gam,all.terms = TRUE,pages=1,scale=0,residuals=F)
summary(M1$gam)
summary(M1$mer)
AICc(M1$mer)
BIC(M1$mer)

cor(fitted(M1$mer), legal$response)^2 #R^2 for lmer model
cor(fitted(M1$gam), legal$response)^2 #R^2 for gam model

## Plot of top-ranked

library(sjPlot)

## Legal ----

gravity.legal<-plot_model(M1,type ="pred",terms="gravity.50",alpha=0.3)  + 
  ggtitle('(c)')  + theme_classic() + Theme1 + xlab("Gravity") +  geom_rug(sides = "b")
gravity.legal

depth.legal<-plot_model(M1,type ="pred",terms="depth",alpha=0.3)  + xlab ("Depth") +
  ggtitle('(d)')+ theme_classic() + Theme1 + ylab("") +  geom_rug(sides = "b")
depth.legal

relief.legal<-plot_model(M1,type ="pred",terms="relief",alpha=0.3)  + 
  ggtitle('(e)')  + theme_classic() + Theme1 +  ylab('Relative abundance (MaxN)') + 
  xlab ("Relief") + geom_rug(sides = "b")
relief.legal

## Plot the effect of management status manually

## Legal ---

unique(dat$state)
unique(dat$campaignid)
unique(dat$clust.4km)
unique(dat$status)

testdata <- expand.grid(depth=mean(legal$depth),
                        relief=mean(legal$relief),
                        gravity.50=mean(legal$gravity.50),
                        state='WA.Temperate Australasia',
                        campaignid='2016-10_North.Kimberley_stereoBRUVs',
                        clust.4km='1',
                        status=c("Fished","No-take"))%>%
  glimpse()

fits <- predict(M1$gam, newdata=testdata, type='response', se.fit=T,exclude = c("s(state)","s(campaignid)","s(clust.4km)"))

La.predicts = testdata%>%data.frame(fits)%>%
  group_by(status)%>% #only change here
  dplyr::summarise(response=mean(fit),se.fit=mean(se.fit))%>%
  ungroup()%>%
  glimpse()
write.csv(La.predicts,"predicts.csv") #there is some BUG in dplyr - that this fixes
La.predicts<-read.csv("predicts.csv")%>%
  glimpse()

Status.legal<- ggplot(aes(x=status,y=response,fill=status), data=La.predicts) +
  ylab("")+
  xlab('Status')+
  #ylim(0,6)+
  scale_fill_manual(labels = c("Fished", "No-take"),values=c("brown1", "blue3"))+
  geom_bar(stat = "identity",show.legend=TRUE,colour="black")+
  geom_errorbar(aes(ymin = response-se.fit,ymax = response+se.fit),width = 0.5,show.legend=FALSE) +
  scale_x_discrete(expand = c(0.3,0.3))+
  scale_y_continuous(expand=c(0,0))+
  theme_classic()+
  Theme1
Status.legal

## Arrange plots in a grob

ggarrange(Status.legal,gravity.legal,depth.legal,relief.legal,
          nrow = 1,ncol = 4,align = 'hv',common.legend = T,legend = "bottom")

#########################################################################################################
################################### (5) Importance scores - Assemblage-level ############################

# colour ramps-
re <- colorRampPalette(c("mistyrose", "red2","darkred"))(200)
bl <- colorRampPalette(c("grey14","grey55","grey62"))(200)  

# Labels-
legend_title<-"Importance"

## Read in importance scores

setwd('C:/Users/22373243/Dropbox/Projects/Analysis/Analysis_GlobalArchive_Bodysize_Nestor/ModelOut/Assemblage-level models/Legal-Abundance')
dir()

dat.taxa<-read.csv("Size_all.var.imp.csv",check.names = F)%>%
  dplyr::rename(Response=Taxa)%>%
  gather(key=predictor,value=importance,2:ncol(.))%>%
  dplyr::rename(Importance=importance)%>%
  dplyr::mutate(Response=fct_relevel(Response,'Sub-legal','Legal','Large'))%>%
  dplyr::mutate(predictor=fct_relevel(predictor,'Status','Gravity','Depth','Relief','SSTmean',
                                      'SSTsd','NPPmean','NPPsd','Nitrate'))%>%
  glimpse()

#gg.importance.scores

ggheat.body.size <- ggplot(dat.taxa, aes(x=predictor,y=Response,fill=Importance))+ 
  geom_tile(show.legend=T) + 
  scale_fill_gradient(low="white",high=re,labels=c('0','0.25','0.5','0.75','1'))+
  xlab(NULL)+
  ylab(NULL)+
  scale_y_discrete(expand=c(0,0))+
  scale_x_discrete(expand=c(0,0))+
  Theme1+
  # ggtitle('(b)')+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust=1,size=18,colour = "black"),
        axis.text.y = element_text(angle = 0,hjust = 1,vjust=1,size=18,colour = "black"),
        legend.spacing.x = unit(1, 'cm'),
        legend.title=element_text(size = 14),
        legend.text=element_text(size=14),
        legend.key.size = unit(1, "cm"),
        legend.position = "none")
# geom_text(aes(label=label),size=10)
ggheat.body.size 

####### Export as a png file

setwd("C:/Users/22373243/Dropbox/Projects/Analysis/Analysis_GlobalArchive_Bodysize_Nestor/Plots") 

name<-'Importance.scores_20200401'

ggsave(paste(name,".png",sep="."),width = 21, height = 6,units = "cm",dpi = 600)
ggsave(paste(name,".pdf",sep="."),width = 21, height = 6,units = "cm",dpi = 600,useDingbats=FALSE)

##########################################################################################################
############################# (6) GLMM models - Regional Species Groups ##################################

#### (A) Chrysophrys auratus (Pink Snapper) ----

# Read in data 

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.species.group.csv")%>%
  filter(species.group%in%c('Chrysophrys auratus'))%>%
  filter(relief<200)%>%
  rename(response=number)%>%
  na.omit()%>%
  glimpse()
unique(dat$species.group)
unique(dat$Taxa)

## Check sampling effort per levels of random effect

test<-dat%>%
  group_by(Ecoregion)%>%
  summarise(campaign_N=n_distinct(campaignid),
            clust_N=n_distinct(clust.4km),
            N=n_distinct(id))%>%
  glimpse()

## Check distribution of gravity values

library(viridis)

ggplot(dat,aes(x=longitude,y=latitude,colour=log1p(gravity.50),size=log1p(gravity.50)))+
  scale_color_viridis()+
  geom_point()+
  theme_classic()

## Code covariates 

dat$slope<-as.numeric(dat$slope)
dat$relief<-as.numeric(dat$relief)
dat$ausbath<-as.numeric(dat$ausbath)
dat$depth<-as.numeric(dat$depth)
dat$t_m<-as.numeric(dat$t_m)
dat$t_sd<-as.numeric(dat$t_sd)
dat$no3_m<-as.numeric(dat$no3_m)
dat$no3_sd<-as.numeric(dat$no3_sd)
dat$po4_m<-as.numeric(dat$po4_m)
dat$po4_sd<-as.numeric(dat$po4_sd)
dat$npp_mean<-as.numeric(dat$npp_mean)
dat$npp_min<-as.numeric(dat$npp_min)
dat$npp_max<-as.numeric(dat$npp_max)
dat$npp_sd<-as.numeric(dat$npp_sd)
dat$dst2ramp<-as.numeric(dat$dst2ramp)
dat$pop.den.50<-as.numeric(dat$pop.den.50)
dat$gravity.50<-as.numeric(dat$gravity.50)
dat$status<-as.factor(dat$status)
dat$campaignid<-as.factor(dat$campaignid)
dat$state<-as.factor(dat$state)
dat$clust.4km<-as.factor(dat$clust.4km)

glimpse(dat)

## Set predictor variables

names(dat)
pred.vars.cont=c("relief","depth","t_m","t_sd","no3_m","npp_mean","npp_sd","gravity.50")

# Plot of likely transformations 

plot.new()

par(mfrow=c(3,2))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

## Transform predictors 

dat$gravity.50<-log1p(dat$gravity.50)
dat$npp_mean<-log1p(dat$npp_mean)
dat$no3_m<-log1p(dat$no3_m)
dat$relief<-log1p(dat$relief)

## Identify outliers 

plot.new()
par(mfrow=c(1,3))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  boxplot(x)
}

## Spatial plots of outliers

ggplot(dat,aes(x=longitude,y=latitude,colour=no3_m))+
  geom_point()+
  scale_color_viridis()+
  theme_classic()+
  Theme1

## Explore initial correlations between covariates

library(corrplot)
library(RColorBrewer)

## Visualize

M<-dat[,pred.vars.cont]
M<-round(cor(M),2)

plot.new()
par(mfrow=c(1,1))
corrplot(M,method="number")

## Scale (mean 0 and SD 1) predictors 

pred.vars.cont

dat$relief<-scale(dat$relief)
dat$depth<-scale(dat$depth)
dat$t_m<-scale(dat$t_m)
dat$t_sd<-scale(dat$t_sd)
dat$no3_m<-scale(dat$no3_m)
dat$npp_mean<-scale(dat$npp_mean)
dat$npp_sd<-scale(dat$npp_sd)
dat$gravity.50<-scale(dat$gravity.50)

### Check distributions to fit individual models
## For each Taxa - species group combination

taxa.name<-unique(dat$Taxa)

plot.new()
par(mfrow=c(3,2))
for (i in taxa.name) {
  hist(dat$response[dat$Taxa==i],main=i)#Looks best
  plot(dat$response[dat$Taxa==i])
}

## Create individual dataframes for modelling

## Large

large<-dat%>%
  dplyr::filter(Taxa=="Large")%>%
  glimpse()
hist(large$response)
plot(large$response)

sum(large$response == 0 ) / length(large$response) ## 92% zeroes

## Convert to presence/absence

large<-large%>%
  mutate(response=ifelse(response>0,1,0))%>%
  glimpse()

## Legal

legal<-dat%>%
  dplyr::filter(Taxa=="Legal")%>%
  glimpse()
hist(legal$response)
plot(legal$response)

sum(legal$response == 0 ) / length(legal$response) ## 83% zeroes

## Sub-legal

sublegal<-dat%>%
  dplyr::filter(Taxa=="Sub-legal")%>%
  glimpse()
hist(sublegal$response)
plot(sublegal$response)

sum(sublegal$response == 0 ) / length(sublegal$response) ## 79% zeroes

## Model selection in MuMIn - via the function dredge
## Fit candidate models based on different combinations of covariates that are ecologically sensible
## Exclude predictors with correlations > 0.5 or VIF > 5
## Limit the complexity to a maximum of 4 predictor variables

library(glmmTMB)

## (A1) Large ----

# First create a full model with all potential predictors and interactions

pred.vars.cont

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = binomial(link = "logit"),
                     data = large,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.5(Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column

sel.table$Model<-rownames(sel.table)

test<-get.models(results,subset = delta <= 2)

test$`173`$call$formula
test$`301`$call$formula
test$`45`$call$formula
test$`47`$call$formula
test$`297`$call$formula
test$`171`$call$formula
test$`43`$call$formula
test$`169`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Chrysophurus_auratus_Large_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"GLM_Chrysophurus_auratus_Importance_Large_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## (A2) Legal ----

# Create a full model

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth  + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = legal,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.5(Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column
sel.table$Model<-rownames(sel.table)

test<-get.models(results,subset = delta <= 2)

test$`297`$call$formula
test$`43`$call$formula
test$`301`$call$formula
test$`293`$call$formula
test$`299`$call$formula
test$`45`$call$formula
test$`47`$call$formula
test$`291`$call$formula
test$`39`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Chrysophurus_auratus_Legal_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"GLM_Chrysophurus_auratus_Importance_Legal_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## (A3) Sub-legal ----

## Create a full model 

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth  + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = sublegal,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.28 (Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column
sel.table$Model<-rownames(sel.table)
sel.table

## Get models 

test<-get.models(results,subset = delta <= 2)

## Write variable names

test$`107`$call$formula
test$`75`$call$formula
test$`331`$call$formula
test$`79`$call$formula
test$`43`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Chrysophurus_auratus_Sublegal_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"GLM_Chrysophurus_auratus_Importance_Sublegal_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## Plot of top models and model validation ----

library(glmmTMB)
library(performance)
library(DHARMa)

## Large 

large.model<- glmmTMB(response ~ gravity.50 + no3_m + npp_mean + relief +
                    (1|Ecoregion/clust.4km) + (1|campaignid),
                  family = binomial(link = "logit"),
                  data = large,REML = TRUE)
r2_nakagawa(large.model)

## Quick check via Dharma

setwd(plots.dir)
dir()

pdf("Large_C.Auratus.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(large.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(large.model)

## Extract Pearson residuals

res1<-resid(large.model, type = "pearson")

## Put the residuals in the dataset

large$res<-res1

## Legal  - (1) NPPmean + relief + SSTsd

legal.model<- glmmTMB(response ~ npp_mean + relief + t_sd +
                    (1|Ecoregion/clust.4km) + (1|campaignid),
                  family = nbinom2(link = "log"),
                  data = legal,REML = TRUE)
r2_nakagawa(legal.model)

## Quick check via Dharma

setwd(plots.dir)
dir()

pdf("Legal_C.Auratus.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(legal.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(legal.model)

## Extract Pearson residuals

res1<-resid(legal.model, type = "pearson")

## Put the residuals in the dataset

legal$res<-res1

## Sublegal - (1) gravity + NPPmean + relief + status

sublegal.model<- glmmTMB(response ~ gravity.50 + npp_mean + relief + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = sublegal,REML = TRUE)

r2_nakagawa(sublegal.model)

## Quick check via Dharma 

setwd(plots.dir)
dir()

pdf("Sub-legal_C.Auratus.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(sublegal.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(sublegal.model)

## Extract Pearson residuals

res1<-resid(sublegal.model, type = "pearson")

## Put the residuals in the dataset

sublegal$res<-res1

### Check for residual auto-correlation ----

## Convert into an spatial object

coordinates(large)<-~longitude+latitude
class(large)

coordinates(legal)<-~longitude+latitude
class(legal)

coordinates(sublegal)<-~longitude+latitude
class(sublegal)

### Now we can set the reference system to the widely used WGS84

proj4string(large)<-CRS("+init=epsg:4326")
large@proj4string

proj4string(legal)<-CRS("+init=epsg:4326")
legal@proj4string

proj4string(sublegal)<-CRS("+init=epsg:4326")
sublegal@proj4string

## Compute semivariogram

library(gstat)

large.var.response<-variogram(large$response~1,data=large,cutoff=20,width=0.1)
legal.var.response<-variogram(legal$response~1,data=legal,cutoff=20,width=0.1)
sublegal.var.response<-variogram(sublegal$response~1,data=legal,cutoff=20,width=0.1)
large.var.res<-variogram(large$res~1,data=large,cutoff=20,width=0.1)
legal.var.res<-variogram(legal$res~1,data=legal,cutoff=20,width=0.1)
sublegal.var.res<-variogram(sublegal$res~1,data=legal,cutoff=20,width=0.1)

## Export semivariograms

## Large response

plot(large.var.response,main="(a) Legal (response)",ylab="Semivariance",xlab="Distance (km)")
plot(large.var.res,main="(b) Legal (residuals)",ylab="Semivariance",xlab="Distance (km)")

## Legal response

plot(legal.var.response,main="(a) Legal (response)",ylab="Semivariance",xlab="Distance (km)")
plot(legal.var.res,main="(b) Legal (residuals)",ylab="Semivariance",xlab="Distance (km)")

## Sub-legal response

plot(sublegal.var.response,main="(c) Sub-legal (residuals-nested)",ylab="Semivariance",xlab="Distance (km)")
plot(sublegal.var.res,main="(c) Sub-legal (residuals)",ylab="Semivariance",xlab="Distance (km)",ylim=c(0.5,1.2))

### Calculate Moran's I based on distance via the spdep package
### We would calculated manually for each basin, diversity index combination and store it in a csv file
### However - for future reference it would be good to loop all of this together

## Global Moran's I ---- 

library(geosphere)

## Large

dists <- distm(large, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(large$response,dists.inv)

## Legal (res 1)

Moran.I(large$res,dists.inv) 

## Legal

dists <- distm(legal, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(legal$response,dists.inv)

## Legal (res 1)

Moran.I(legal$res,dists.inv) 

## Sub-legal

dists <- distm(sublegal, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(sublegal$response,dists.inv)

## Legal (res)

Moran.I(sublegal$res,dists.inv) 

#### (B) Lethrinus spp. ----

# Read in data 

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.species.group.csv")%>%
  filter(species.group%in%c('Lethrinus spp.'))%>%
  filter(relief<200)%>%
  rename(response=number)%>%
  na.omit()%>%
  glimpse()
unique(dat$species.group)
unique(dat$Taxa)

## Check sampling effort per levels of random effect

test<-dat%>%
  group_by(Ecoregion)%>%
  summarise(campaign_N=n_distinct(campaignid),
            clust_N=n_distinct(clust.4km),
            N=n_distinct(id))%>%
  glimpse()

## Check distribution of gravity values

library(viridis)

ggplot(dat,aes(x=longitude,y=latitude,colour=log1p(gravity.50),size=log1p(gravity.50)))+
  scale_color_viridis()+
  geom_point()+
  theme_classic()

## Code covariates 

dat$slope<-as.numeric(dat$slope)
dat$relief<-as.numeric(dat$relief)
dat$ausbath<-as.numeric(dat$ausbath)
dat$depth<-as.numeric(dat$depth)
dat$t_m<-as.numeric(dat$t_m)
dat$t_sd<-as.numeric(dat$t_sd)
dat$no3_m<-as.numeric(dat$no3_m)
dat$no3_sd<-as.numeric(dat$no3_sd)
dat$po4_m<-as.numeric(dat$po4_m)
dat$po4_sd<-as.numeric(dat$po4_sd)
dat$npp_mean<-as.numeric(dat$npp_mean)
dat$npp_min<-as.numeric(dat$npp_min)
dat$npp_max<-as.numeric(dat$npp_max)
dat$npp_sd<-as.numeric(dat$npp_sd)
dat$dst2ramp<-as.numeric(dat$dst2ramp)
dat$pop.den.50<-as.numeric(dat$pop.den.50)
dat$gravity.50<-as.numeric(dat$gravity.50)
dat$status<-as.factor(dat$status)
dat$campaignid<-as.factor(dat$campaignid)
dat$state<-as.factor(dat$state)
dat$clust.4km<-as.factor(dat$clust.4km)

glimpse(dat)

## Set predictor variables

names(dat)
pred.vars.cont=c("relief","depth","t_m","t_sd","no3_m","npp_mean","npp_sd","gravity.50")

# Plot of likely transformations 

plot.new()

par(mfrow=c(3,2))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

## Transform predictors 

dat$gravity.50<-log1p(dat$gravity.50)
dat$npp_sd<-log1p(dat$npp_sd)
dat$npp_mean<-log1p(dat$npp_mean)
dat$no3_m<-log1p(dat$no3_m)
dat$relief<-log1p(dat$relief)

## Identify outliers 

plot.new()
par(mfrow=c(1,3))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  boxplot(x)
}

## Spatial plots of outliers

ggplot(dat,aes(x=longitude,y=latitude,colour=no3_m))+
  geom_point()+
  scale_color_viridis()+
  theme_classic()+
  Theme1

## Explore initial correlations between covariates

library(corrplot)
library(RColorBrewer)

## Visualize

M<-dat[,pred.vars.cont]
M<-round(cor(M),2)

plot.new()
par(mfrow=c(1,1))
corrplot(M,method="number")

## Scale (mean 0 and SD 1) predictors 

pred.vars.cont

dat$relief<-scale(dat$relief)
dat$depth<-scale(dat$depth)
dat$t_m<-scale(dat$t_m)
dat$t_sd<-scale(dat$t_sd)
dat$no3_m<-scale(dat$no3_m)
dat$npp_mean<-scale(dat$npp_mean)
dat$npp_sd<-scale(dat$npp_sd)
dat$gravity.50<-scale(dat$gravity.50)

## Create individual dataframes for modelling

plot.new()
par(mfrow=c(3,2))

## Large

large<-dat%>%
  dplyr::filter(Taxa=="Large")%>%
  glimpse()
hist(large$response)
plot(large$response)

sum(large$response == 0 ) / length(large$response) ## 88% zeroes

## Convert to presence/absence

large<-large%>%
  mutate(response=ifelse(response>0,1,0))%>%
  glimpse()

## Legal

legal<-dat%>%
  dplyr::filter(Taxa=="Legal")%>%
  glimpse()
hist(legal$response)
plot(legal$response)

sum(legal$response == 0 ) / length(legal$response) ## 64% zeroes

## Sub-legal

sublegal<-dat%>%
  dplyr::filter(Taxa=="Sub-legal")%>%
  glimpse()
hist(sublegal$response)
plot(sublegal$response)

sum(sublegal$response == 0 ) / length(sublegal$response) ## 60% zeroes

## Run models 

## (B1) Large ----

# that's better

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = binomial(link = "logit"),
                     data = large,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.5(Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column

sel.table$Model<-rownames(sel.table)

test<-get.models(results,subset = delta <= 2)

test$`331`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Letrhinus_spp._Large_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"GLM_Lethrinidae_Importance_Large_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## (B2) Legal ----

## Fit a full model 

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth  + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = legal,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.5(Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column
sel.table$Model<-rownames(sel.table)

test<-get.models(results,subset = delta <= 2)

test$`269`$call$formula
test$`77`$call$formula
test$`329`$call$formula
test$`267`$call$formula
test$`15`$call$formula
test$`281`$call$formula
test$`75`$call$formula
test$`29`$call$formula
test$`333`$call$formula
test$`27`$call$formula
test$`89`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Lethrinidae_Legal_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"GLM_Lethrinidae_Importance_Legal_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## (B3) Sub-legal ----

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth  + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = sublegal,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.28 (Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column
sel.table$Model<-rownames(sel.table)
sel.table

## Get models 

test<-get.models(results,subset = delta <= 2)

## Write variable names

test$`138`$call$formula
test$`142`$call$formula
test$`134`$call$formula
test$`154`$call$formula
test$`386`$call$formula
test$`394`$call$formula
test$`390`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Lethrinidae_Sublegal_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"Lethrinus_Nebulosus_Lethrinidae_Sublegal_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## Top-ranked models and model validation----

## Large 

large.model<- glmmTMB(response ~  gravity.50 + npp_mean + status + t_sd +
                      (1|Ecoregion/clust.4km) + (1|campaignid),
                      family = binomial(link = "logit"),
                      data = large,REML = TRUE)
r2_nakagawa(large.model)

## Quick check via Dharma

setwd(plots.dir)
dir()

pdf("Large_Lethrinus_spp..pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(large.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(large.model)

## Extract Pearson residuals

res1<-resid(large.model, type = "pearson")

## Put the residuals in the dataset

large$res<-res1

## Legal - (1) Nitrate + NPPmean + SSTsd

legal.model<- glmmTMB(response ~ no3_m + npp_mean + t_sd +
                    (1|Ecoregion/clust.4km) + (1|campaignid),
                  family = nbinom2(link = "log"),
                  data = legal,REML = TRUE)
r2_nakagawa(legal.model)

## Quick check via Dharma

setwd(plots.dir)
dir()

pdf("Legal_Lethrinus_spp.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(legal.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(legal.model)

## Extract Pearson residuals

res1<-resid(legal.model, type = "pearson")

## Put the residuals in the dataset

legal$res<-res1

## Sublegal - (1) depth + NPPmean + SSTmean

sublegal.model<- glmmTMB(response ~ depth + npp_mean + t_m +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = sublegal,REML = TRUE)
r2_nakagawa(sublegal.model)

## Quick check via Dharma 

setwd(plots.dir)
dir()

pdf("Sub-legal_Letrhinus_spp.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(sublegal.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(sublegal.model)

## Extract Pearson residuals

res1<-resid(sublegal.model, type = "pearson")

## Put the residuals in the dataset

sublegal$res<-res1

### Check for residual auto-correlation ----

## Convert into an spatial object

coordinates(large)<-~longitude+latitude
class(large)

coordinates(legal)<-~longitude+latitude
class(legal)

coordinates(sublegal)<-~longitude+latitude
class(sublegal)

### Now we can set the reference system to the widely used WGS84

proj4string(large)<-CRS("+init=epsg:4326")
large@proj4string

proj4string(legal)<-CRS("+init=epsg:4326")
legal@proj4string

proj4string(sublegal)<-CRS("+init=epsg:4326")
sublegal@proj4string

## Compute semivariogram

library(gstat)

large.var.response<-variogram(large$response~1,data=large,cutoff=20,width=0.1)
legal.var.response<-variogram(legal$response~1,data=legal,cutoff=20,width=0.1)
sublegal.var.response<-variogram(sublegal$response~1,data=legal,cutoff=20,width=0.1)
large.var.res<-variogram(large$res~1,data=large,cutoff=20,width=0.1)
legal.var.res<-variogram(legal$res~1,data=legal,cutoff=20,width=0.1)
sublegal.var.res<-variogram(sublegal$res~1,data=legal,cutoff=20,width=0.1)

## Export semivariograms

## Large response

plot(large.var.response,main="(a) Legal (response)",ylab="Semivariance",xlab="Distance (km)")
plot(large.var.res,main="(b) Legal (residuals)",ylab="Semivariance",xlab="Distance (km)")

## Legal response

plot(legal.var.response,main="(a) Legal (response)",ylab="Semivariance",xlab="Distance (km)")
plot(legal.var.res,main="(b) Legal (residuals)",ylab="Semivariance",xlab="Distance (km)")

## Sub-legal response

plot(sublegal.var.response,main="(c) Sub-legal (residuals-nested)",ylab="Semivariance",xlab="Distance (km)")
plot(sublegal.var.res,main="(c) Sub-legal (residuals)",ylab="Semivariance",xlab="Distance (km)",ylim=c(0.5,1.2))

### Calculate Moran's I based on distance via the spdep package
### We would calculated manually for each basin, diversity index combination and store it in a csv file
### However - for future reference it would be good to loop all of this together

## Global Moran's I ---- 

library(geosphere)

## Large

dists <- distm(large, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(large$response,dists.inv)

## Legal (res 1)

Moran.I(large$res,dists.inv) 

## Legal

dists <- distm(legal, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(legal$response,dists.inv)

## Legal (res 1)

Moran.I(legal$res,dists.inv) 

## Sub-legal

dists <- distm(sublegal, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(sublegal$response,dists.inv)

## Legal (res)

Moran.I(sublegal$res,dists.inv) 

#### (C) Plectropomus spp ----

# Read in data 

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.species.group.csv")%>%
  filter(species.group%in%c('Plectropomus spp.'))%>%
  filter(relief<200)%>%
  rename(response=number)%>%
  na.omit()%>%
  glimpse()
unique(dat$species.group)
unique(dat$Taxa)

## Check sampling effort per levels of random effect

test<-dat%>%
  group_by(Ecoregion)%>%
  summarise(campaign_N=n_distinct(campaignid),
            clust_N=n_distinct(clust.4km),
            N=n_distinct(id))%>%
  glimpse()

## Check distribution of gravity values

library(viridis)

ggplot(dat,aes(x=longitude,y=latitude,colour=log1p(gravity.50),size=log1p(gravity.50)))+
  scale_color_viridis()+
  geom_point()+
  theme_classic()

## Code covariates 

dat$slope<-as.numeric(dat$slope)
dat$relief<-as.numeric(dat$relief)
dat$ausbath<-as.numeric(dat$ausbath)
dat$depth<-as.numeric(dat$depth)
dat$t_m<-as.numeric(dat$t_m)
dat$t_sd<-as.numeric(dat$t_sd)
dat$no3_m<-as.numeric(dat$no3_m)
dat$no3_sd<-as.numeric(dat$no3_sd)
dat$po4_m<-as.numeric(dat$po4_m)
dat$po4_sd<-as.numeric(dat$po4_sd)
dat$npp_mean<-as.numeric(dat$npp_mean)
dat$npp_min<-as.numeric(dat$npp_min)
dat$npp_max<-as.numeric(dat$npp_max)
dat$npp_sd<-as.numeric(dat$npp_sd)
dat$dst2ramp<-as.numeric(dat$dst2ramp)
dat$pop.den.50<-as.numeric(dat$pop.den.50)
dat$gravity.50<-as.numeric(dat$gravity.50)
dat$status<-as.factor(dat$status)
dat$campaignid<-as.factor(dat$campaignid)
dat$state<-as.factor(dat$state)
dat$clust.4km<-as.factor(dat$clust.4km)

glimpse(dat)

## Set predictor variables

names(dat)
pred.vars.cont=c("relief","depth","t_m","t_sd","no3_m","npp_mean","npp_sd","gravity.50")

# Plot of likely transformations 

plot.new()

par(mfrow=c(3,2))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

## Transform predictors 

dat$gravity.50<-log1p(dat$gravity.50)
dat$npp_sd<-log1p(dat$npp_sd)
dat$npp_mean<-log1p(dat$npp_mean)
dat$no3_m<-log1p(dat$no3_m)
dat$relief<-log1p(dat$relief)

## Identify outliers 

plot.new()
par(mfrow=c(1,3))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  boxplot(x)
}

## Spatial plots of outliers

ggplot(dat,aes(x=longitude,y=latitude,colour=no3_m))+
  geom_point()+
  scale_color_viridis()+
  theme_classic()+
  Theme1

## Explore initial correlations between covariates

library(corrplot)
library(RColorBrewer)

## Visualize

M<-dat[,pred.vars.cont]
M<-round(cor(M),2)

plot.new()
par(mfrow=c(1,1))
corrplot(M,method="number")

## Scale (mean 0 and SD 1) predictors 

pred.vars.cont

dat$relief<-scale(dat$relief)
dat$depth<-scale(dat$depth)
dat$t_m<-scale(dat$t_m)
dat$t_sd<-scale(dat$t_sd)
dat$no3_m<-scale(dat$no3_m)
dat$npp_mean<-scale(dat$npp_mean)
dat$npp_sd<-scale(dat$npp_sd)
dat$gravity.50<-scale(dat$gravity.50)

## Create individual dataframes for modelling

plot.new()
par(mfrow=c(3,2))

## Large

large<-dat%>%
  dplyr::filter(Taxa=="Large")%>%
  glimpse()
hist(large$response)
plot(large$response)

sum(large$response == 0 ) / length(large$response) ## 95% zeroes

## Convert to presence/absence

large<-large%>%
  mutate(response=ifelse(response>0,1,0))%>%
  glimpse()

## Legal

legal<-dat%>%
  dplyr::filter(Taxa=="Legal")%>%
  glimpse()
hist(legal$response)
plot(legal$response)

sum(legal$response == 0 ) / length(legal$response) ## 92% zeroes

## Sub-legal

sublegal<-dat%>%
  dplyr::filter(Taxa=="Sub-legal")%>%
  glimpse()
hist(sublegal$response)
plot(sublegal$response)

sum(sublegal$response == 0 ) / length(sublegal$response) ## 82% zeroes

## Run models 

library(glmmTMB)

## (C1) Large ----

# that's better

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = binomial(link = "logit"),
                     data = large,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.5(Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column

sel.table$Model<-rownames(sel.table)

test<-get.models(results,subset = delta <= 2)

test$`77`$call$formula
test$`93`$call$formula
test$`89`$call$formula
test$`333`$call$formula
test$`109`$call$formula
test$`329`$call$formula
test$`105`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Plectropomus_spp._Large_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"GLM_Plectropomus_Importance_Large_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## (C2) Legal ----

## Fit a full model 

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth  + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = legal,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.5(Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column
sel.table$Model<-rownames(sel.table)

test<-get.models(results,subset = delta <= 2)

test$`45`$call$formula
test$`105`$call$formula
test$`109`$call$formula
test$`101`$call$formula
test$`77`$call$formula
test$`57`$call$formula
test$`297`$call$formula
test$`293`$call$formula
test$`301`$call$formula
test$`61`$call$formula
test$`587`$call$formula
test$`121`$call$formula
test$`353`$call$formula
test$`165`$call$formula
test$`269`$call$formula
test$`357`$call$formula
test$`361`$call$formula
test$`169`$call$formula
test$`173`$call$formula
test$`15`$call$formula
test$`47`$call$formula
test$`53`$call$formula
test$`329`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Plectropomus_Legal_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"GLM_Plectropomus_Importance_Legal_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## (C3) Sub-legal ----

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth  + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = sublegal,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.28 (Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column
sel.table$Model<-rownames(sel.table)
sel.table

## Get models 

test<-get.models(results,subset = delta <= 2)

## Write variable names

test$`54`$call$formula
test$`114`$call$formula
test$`50`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Plectropomus_Sublegal_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"Plectropomus_Sublegal_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## Top-ranked models and model validation----

## Large 

large.model<- glmmTMB(response ~ no3_m + npp_mean + status +
                        (1|Ecoregion/clust.4km) + (1|campaignid),
                      family = binomial(link = "logit"),
                      data = large,REML = TRUE)
r2_nakagawa(large.model)

## Quick check via Dharma

setwd(plots.dir)
dir()

pdf("Large_Plectropomus_spp..pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(large.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(large.model)

## Extract Pearson residuals

res1<-resid(large.model, type = "pearson")

## Put the residuals in the dataset

large$res<-res1

## Legal - (1) Nitrate + NPPmean + SSTsd

legal.model<- glmmTMB(response ~ no3_m + npp_mean + t_sd +
                    (1|Ecoregion/clust.4km) + (1|campaignid),
                  family = nbinom2(link = "log"),
                  data = legal,REML = TRUE)
r2_nakagawa(legal.model)

## Quick check via Dharma

setwd(plots.dir)
dir()

pdf("Legal_Plectropomus_spp.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(legal.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(legal.model)

## Extract Pearson residuals

res1<-resid(legal.model, type = "pearson")

## Put the residuals in the dataset

legal$res<-res1

## Sublegal - (1) depth + NPPmean + SSTmean

sublegal.model<- glmmTMB(response ~ depth + npp_sd + no3_m + relief +
                     (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = sublegal,REML = TRUE)
r2_nakagawa(sublegal.model)


## Quick check via Dharma 

setwd(plots.dir)
dir()

pdf("Sub-legal_Plectropomus_spp.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(sublegal.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(sublegal.model)

## Extract Pearson residuals

res1<-resid(sublegal.model, type = "pearson")

## Put the residuals in the dataset

sublegal$res<-res1

### Check for residual auto-correlation ----

## Convert into an spatial object

coordinates(large)<-~longitude+latitude
class(large)

coordinates(legal)<-~longitude+latitude
class(legal)

coordinates(sublegal)<-~longitude+latitude
class(sublegal)

### Now we can set the reference system to the widely used WGS84

proj4string(large)<-CRS("+init=epsg:4326")
large@proj4string

proj4string(legal)<-CRS("+init=epsg:4326")
legal@proj4string

proj4string(sublegal)<-CRS("+init=epsg:4326")
sublegal@proj4string

## Compute semivariogram

library(gstat)

large.var.response<-variogram(large$response~1,data=large,cutoff=20,width=0.1)
legal.var.response<-variogram(legal$response~1,data=legal,cutoff=20,width=0.1)
sublegal.var.response<-variogram(sublegal$response~1,data=legal,cutoff=20,width=0.1)
large.var.res<-variogram(large$res~1,data=large,cutoff=20,width=0.1)
legal.var.res<-variogram(legal$res~1,data=legal,cutoff=20,width=0.1)
sublegal.var.res<-variogram(sublegal$res~1,data=legal,cutoff=20,width=0.1)

## Export semivariograms

## Large response

plot(large.var.response,main="(a) Legal (response)",ylab="Semivariance",xlab="Distance (km)")
plot(large.var.res,main="(b) Legal (residuals)",ylab="Semivariance",xlab="Distance (km)")

## Legal response

plot(legal.var.response,main="(a) Legal (response)",ylab="Semivariance",xlab="Distance (km)")
plot(legal.var.res,main="(b) Legal (residuals)",ylab="Semivariance",xlab="Distance (km)")

## Sub-legal response

plot(sublegal.var.response,main="(c) Sub-legal (residuals-nested)",ylab="Semivariance",xlab="Distance (km)")
plot(sublegal.var.res,main="(c) Sub-legal (residuals)",ylab="Semivariance",xlab="Distance (km)",ylim=c(0.5,1.2))

### Calculate Moran's I based on distance via the spdep package
### We would calculated manually for each basin, diversity index combination and store it in a csv file
### However - for future reference it would be good to loop all of this together

## Global Moran's I ---- 

library(geosphere)

## Large

dists <- distm(large, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(large$response,dists.inv)

## Legal (res 1)

Moran.I(large$res,dists.inv) 

## Legal

dists <- distm(legal, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(legal$response,dists.inv)

## Legal (res 1)

Moran.I(legal$res,dists.inv) 

## Sub-legal

dists <- distm(sublegal, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(sublegal$response,dists.inv)

## Legal (res)

Moran.I(sublegal$res,dists.inv) 

#### (D) Lutjanus spp ----

# Read in data 

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.species.group.csv")%>%
  filter(species.group%in%c('Lutjanus spp.'))%>%
  filter(relief<200)%>%
  rename(response=number)%>%
  na.omit()%>%
  glimpse()
unique(dat$species.group)
unique(dat$Taxa)

## Check sampling effort per levels of random effect

test<-dat%>%
  group_by(Ecoregion)%>%
  summarise(campaign_N=n_distinct(campaignid),
            clust_N=n_distinct(clust.4km),
            N=n_distinct(id))%>%
  glimpse()

## Check distribution of gravity values

library(viridis)

ggplot(dat,aes(x=longitude,y=latitude,colour=log1p(gravity.50),size=log1p(gravity.50)))+
  scale_color_viridis()+
  geom_point()+
  theme_classic()

## Code covariates 

dat$slope<-as.numeric(dat$slope)
dat$relief<-as.numeric(dat$relief)
dat$ausbath<-as.numeric(dat$ausbath)
dat$depth<-as.numeric(dat$depth)
dat$t_m<-as.numeric(dat$t_m)
dat$t_sd<-as.numeric(dat$t_sd)
dat$no3_m<-as.numeric(dat$no3_m)
dat$no3_sd<-as.numeric(dat$no3_sd)
dat$po4_m<-as.numeric(dat$po4_m)
dat$po4_sd<-as.numeric(dat$po4_sd)
dat$npp_mean<-as.numeric(dat$npp_mean)
dat$npp_min<-as.numeric(dat$npp_min)
dat$npp_max<-as.numeric(dat$npp_max)
dat$npp_sd<-as.numeric(dat$npp_sd)
dat$dst2ramp<-as.numeric(dat$dst2ramp)
dat$pop.den.50<-as.numeric(dat$pop.den.50)
dat$gravity.50<-as.numeric(dat$gravity.50)
dat$status<-as.factor(dat$status)
dat$campaignid<-as.factor(dat$campaignid)
dat$state<-as.factor(dat$state)
dat$clust.4km<-as.factor(dat$clust.4km)

glimpse(dat)

## Set predictor variables

names(dat)
pred.vars.cont=c("relief","depth","t_m","t_sd","no3_m","npp_mean","npp_sd","gravity.50")

# Plot of likely transformations 

plot.new()

par(mfrow=c(3,2))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

## Transform predictors 

dat$gravity.50<-log1p(dat$gravity.50)
dat$npp_sd<-log1p(dat$npp_sd)
dat$npp_mean<-log1p(dat$npp_mean)
dat$no3_m<-log1p(dat$no3_m)
dat$relief<-log1p(dat$relief)

## Identify outliers 

plot.new()
par(mfrow=c(1,3))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  boxplot(x)
}

## Spatial plots of outliers

ggplot(dat,aes(x=longitude,y=latitude,colour=no3_m))+
  geom_point()+
  scale_color_viridis()+
  theme_classic()+
  Theme1

## Explore initial correlations between covariates

library(corrplot)
library(RColorBrewer)

## Visualize

M<-dat[,pred.vars.cont]
M<-round(cor(M),2)

plot.new()
par(mfrow=c(1,1))
corrplot(M,method="number")

## Scale (mean 0 and SD 1) predictors 

pred.vars.cont

dat$relief<-scale(dat$relief)
dat$depth<-scale(dat$depth)
dat$t_m<-scale(dat$t_m)
dat$t_sd<-scale(dat$t_sd)
dat$no3_m<-scale(dat$no3_m)
dat$npp_mean<-scale(dat$npp_mean)
dat$npp_sd<-scale(dat$npp_sd)
dat$gravity.50<-scale(dat$gravity.50)

## Create individual dataframes for modelling

plot.new()
par(mfrow=c(3,2))

## Large

large<-dat%>%
  dplyr::filter(Taxa=="Large")%>%
  glimpse()
hist(large$response)
plot(large$response)

sum(large$response == 0 ) / length(large$response) ## 90% zeroes

## Convert to presence/absence

large<-large%>%
  mutate(response=ifelse(response>0,1,0))%>%
  glimpse()

## Legal

legal<-dat%>%
  dplyr::filter(Taxa=="Legal")%>%
  glimpse()
hist(legal$response)
plot(legal$response)

sum(legal$response == 0 ) / length(legal$response) ## 87% zeroes

## Sub-legal

sublegal<-dat%>%
  dplyr::filter(Taxa=="Sub-legal")%>%
  glimpse()
hist(sublegal$response)
plot(sublegal$response)

sum(sublegal$response == 0 ) / length(sublegal$response) ## 70% zeroes

## Run models 

library(glmmTMB)

## (D1) Large ----

# that's better

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = binomial(link = "logit"),
                     data = large,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.5(Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column

sel.table$Model<-rownames(sel.table)

test<-get.models(results,subset = delta <= 2)

test$`197`$call$formula
test$`389`$call$formula
test$`141`$call$formula
test$`201`$call$formula
test$`149`$call$formula
test$`205`$call$formula
test$`453`$call$formula
test$`449`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Lutjanus_spp._Large_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"GLM_Lutjanus_Importance_Large_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## (D2) Legal ----

## Fit a full model 

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth  + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = legal,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.5(Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column
sel.table$Model<-rownames(sel.table)

test<-get.models(results,subset = delta <= 2)

test$`325`$call$formula
test$`269`$call$formula
test$`277`$call$formula
test$`389`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Plectropomus_Legal_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"GLM_Plectropomus_Importance_Legal_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## (D3) Sub-legal ----

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth  + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = sublegal,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.28 (Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column
sel.table$Model<-rownames(sel.table)
sel.table

## Get models 

test<-get.models(results,subset = delta <= 2)

## Write variable names

test$`262`$call$formula
test$`70`$call$formula
test$`22`$call$formula
test$`134`$call$formula
test$`326`$call$formula
test$`278`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Plectropomus_Sublegal_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"Plectropomus_Sublegal_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## Top-ranked models and model validation ----

## Large 

large.model<- glmmTMB(response ~ no3_m + status + t_m +
                      (1|Ecoregion/clust.4km) + (1|campaignid),
                      family = binomial(link = "logit"),
                      data = large,REML = TRUE)
r2_nakagawa(large.model)

## Quick check via Dharma

setwd(plots.dir)
dir()

pdf("Large_Lutjanus_spp.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(large.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(large.model)

## Extract Pearson residuals

res1<-resid(large.model, type = "pearson")

## Put the residuals in the dataset

large$res<-res1

## Legal - (1) Nitrate + NPPmean + SSTsd

legal.model<- glmmTMB(response ~ no3_m + npp_mean + t_sd +
                    (1|Ecoregion/clust.4km) + (1|campaignid),
                  family = nbinom2(link = "log"),
                  data = legal,REML = TRUE)
r2_nakagawa(legal.model)

## Quick check via Dharma

setwd(plots.dir)
dir()

pdf("Legal_Lutjanus_spp.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(legal.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(legal.model)

## Extract Pearson residuals

res1<-resid(legal.model, type = "pearson")

## Put the residuals in the dataset

legal$res<-res1

## Sublegal - (1) depth + NPPmean + SSTmean

sublegal.model<-glmmTMB(response ~ depth + npp_mean + t_m +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = sublegal,REML = TRUE)
r2_nakagawa(sublegal.model)

## Quick check via Dharma 

setwd(plots.dir)
dir()

pdf("Sub-legal_Lutjanus_spp.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(sublegal.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(sublegal.model)

## Extract Pearson residuals

res1<-resid(sublegal.model, type = "pearson")

## Put the residuals in the dataset

sublegal$res<-res1

### Check for residual auto-correlation ----

## Convert into an spatial object

coordinates(large)<-~longitude+latitude
class(large)

coordinates(legal)<-~longitude+latitude
class(legal)

coordinates(sublegal)<-~longitude+latitude
class(sublegal)

### Now we can set the reference system to the widely used WGS84

proj4string(large)<-CRS("+init=epsg:4326")
large@proj4string

proj4string(legal)<-CRS("+init=epsg:4326")
legal@proj4string

proj4string(sublegal)<-CRS("+init=epsg:4326")
sublegal@proj4string

## Compute semivariogram

library(gstat)

large.var.response<-variogram(large$response~1,data=large,cutoff=20,width=0.1)
legal.var.response<-variogram(legal$response~1,data=legal,cutoff=20,width=0.1)
sublegal.var.response<-variogram(sublegal$response~1,data=legal,cutoff=20,width=0.1)
large.var.res<-variogram(large$res~1,data=large,cutoff=20,width=0.1)
legal.var.res<-variogram(legal$res~1,data=legal,cutoff=20,width=0.1)
sublegal.var.res<-variogram(sublegal$res~1,data=legal,cutoff=20,width=0.1)

## Export semivariograms

## Large response

plot(large.var.response,main="(a) Legal (response)",ylab="Semivariance",xlab="Distance (km)")
plot(large.var.res,main="(b) Legal (residuals)",ylab="Semivariance",xlab="Distance (km)")

## Legal response

plot(legal.var.response,main="(a) Legal (response)",ylab="Semivariance",xlab="Distance (km)")
plot(legal.var.res,main="(b) Legal (residuals)",ylab="Semivariance",xlab="Distance (km)")

## Sub-legal response

plot(sublegal.var.response,main="(c) Sub-legal (residuals-nested)",ylab="Semivariance",xlab="Distance (km)")
plot(sublegal.var.res,main="(c) Sub-legal (residuals)",ylab="Semivariance",xlab="Distance (km)",ylim=c(0.5,1.2))

### Calculate Moran's I based on distance via the spdep package
### We would calculated manually for each basin, diversity index combination and store it in a csv file
### However - for future reference it would be good to loop all of this together

## Global Moran's I ---- 

library(geosphere)

## Large

dists <- distm(large, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(large$response,dists.inv)

## Legal (res 1)

Moran.I(large$res,dists.inv) 

## Legal

dists <- distm(legal, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(legal$response,dists.inv)

## Legal (res 1)

Moran.I(legal$res,dists.inv) 

## Sub-legal

dists <- distm(sublegal, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(sublegal$response,dists.inv)

## Legal (res)

Moran.I(sublegal$res,dists.inv) 

#### (E) Choerodon spp ----

# Read in data 

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.species.group.csv")%>%
  filter(species.group%in%c('Choerodon spp.'))%>%
  filter(relief<200)%>%
  rename(response=number)%>%
  na.omit()%>%
  glimpse()
unique(dat$species.group)
unique(dat$Taxa)

## Check sampling effort per levels of random effect

test<-dat%>%
  group_by(Ecoregion)%>%
  summarise(campaign_N=n_distinct(campaignid),
            clust_N=n_distinct(clust.4km),
            N=n_distinct(id))%>%
  glimpse()

## Check distribution of gravity values

library(viridis)

ggplot(dat,aes(x=longitude,y=latitude,colour=log1p(gravity.50),size=log1p(gravity.50)))+
  scale_color_viridis()+
  geom_point()+
  theme_classic()

## Code covariates 

dat$slope<-as.numeric(dat$slope)
dat$relief<-as.numeric(dat$relief)
dat$ausbath<-as.numeric(dat$ausbath)
dat$depth<-as.numeric(dat$depth)
dat$t_m<-as.numeric(dat$t_m)
dat$t_sd<-as.numeric(dat$t_sd)
dat$no3_m<-as.numeric(dat$no3_m)
dat$no3_sd<-as.numeric(dat$no3_sd)
dat$po4_m<-as.numeric(dat$po4_m)
dat$po4_sd<-as.numeric(dat$po4_sd)
dat$npp_mean<-as.numeric(dat$npp_mean)
dat$npp_min<-as.numeric(dat$npp_min)
dat$npp_max<-as.numeric(dat$npp_max)
dat$npp_sd<-as.numeric(dat$npp_sd)
dat$dst2ramp<-as.numeric(dat$dst2ramp)
dat$pop.den.50<-as.numeric(dat$pop.den.50)
dat$gravity.50<-as.numeric(dat$gravity.50)
dat$status<-as.factor(dat$status)
dat$campaignid<-as.factor(dat$campaignid)
dat$state<-as.factor(dat$state)
dat$clust.4km<-as.factor(dat$clust.4km)

glimpse(dat)

## Set predictor variables

names(dat)
pred.vars.cont=c("relief","depth","t_m","t_sd","no3_m","npp_mean","npp_sd","gravity.50")

# Plot of likely transformations 

plot.new()

par(mfrow=c(3,2))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

## Transform predictors 

dat$gravity.50<-log1p(dat$gravity.50)
dat$npp_sd<-log1p(dat$npp_sd)
dat$npp_mean<-log1p(dat$npp_mean)
dat$no3_m<-log1p(dat$no3_m)
dat$relief<-log1p(dat$relief)

## Identify outliers 

plot.new()
par(mfrow=c(1,3))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  boxplot(x)
}

## Spatial plots of outliers

ggplot(dat,aes(x=longitude,y=latitude,colour=no3_m))+
  geom_point()+
  scale_color_viridis()+
  theme_classic()+
  Theme1

## Explore initial correlations between covariates

library(corrplot)
library(RColorBrewer)

## Visualize

M<-dat[,pred.vars.cont]
M<-round(cor(M),2)

plot.new()
par(mfrow=c(1,1))
corrplot(M,method="number")

## Scale (mean 0 and SD 1) predictors 

pred.vars.cont

dat$relief<-scale(dat$relief)
dat$depth<-scale(dat$depth)
dat$t_m<-scale(dat$t_m)
dat$t_sd<-scale(dat$t_sd)
dat$no3_m<-scale(dat$no3_m)
dat$npp_mean<-scale(dat$npp_mean)
dat$npp_sd<-scale(dat$npp_sd)
dat$gravity.50<-scale(dat$gravity.50)

## Create individual dataframes for modelling

plot.new()
par(mfrow=c(3,2))

## Large

large<-dat%>%
  dplyr::filter(Taxa=="Large")%>%
  glimpse()
hist(large$response)
plot(large$response)

sum(large$response == 0 ) / length(large$response) ## 92% zeroes

## Convert to presence/absence

large<-large%>%
  mutate(response=ifelse(response>0,1,0))%>%
  glimpse()

## Legal

legal<-dat%>%
  dplyr::filter(Taxa=="Legal")%>%
  glimpse()
hist(legal$response)
plot(legal$response)

sum(legal$response == 0 ) / length(legal$response) ## 86% zeroes

## Sub-legal

sublegal<-dat%>%
  dplyr::filter(Taxa=="Sub-legal")%>%
  glimpse()
hist(sublegal$response)
plot(sublegal$response)

sum(sublegal$response == 0 ) / length(sublegal$response) ## 70% zeroes

## Run models 

library(glmmTMB)

## (E1) Large ----

# that's better

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = binomial(link = "logit"),
                     data = large,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.5(Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column

sel.table$Model<-rownames(sel.table)

test<-get.models(results,subset = delta <= 2)

test$`293`$call$formula
test$`309`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Choerodon_spp._Large_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"GLM_Choerodon_Importance_Large_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## (E2) Legal ----

## Fit a full model 

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth  + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = legal,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.5(Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column
sel.table$Model<-rownames(sel.table)

test<-get.models(results,subset = delta <= 2)

test$`39`$call$formula
test$`47`$call$formula
test$`299`$call$formula
test$`43`$call$formula
test$`295`$call$formula
test$`291`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Plectropomus_Legal_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"GLM_Plectropomus_Importance_Legal_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## (E3) Sub-legal ----

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth  + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = sublegal,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.28 (Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column
sel.table$Model<-rownames(sel.table)
sel.table

## Get models 

test<-get.models(results,subset = delta <= 2)

## Write variable names

test$`390`$call$formula
test$`386`$call$formula
test$`394`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Plectropomus_Sublegal_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"Plectropomus_Sublegal_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## Top-ranked models and model validation ----

## Large 

large.model<- glmmTMB(response ~   no3_m + relief + t_sd +
                        (1|Ecoregion/clust.4km) + (1|campaignid),
                      family = binomial(link = "logit"),
                      data = large,REML = TRUE)
r2_nakagawa(large.model)

## Quick check via Dharma

setwd(plots.dir)
dir()

pdf("Large_Choerodon_spp.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(large.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(large.model)

## Extract Pearson residuals

res1<-resid(large.model, type = "pearson")

## Put the residuals in the dataset

large$res<-res1

## Legal - (1) Nitrate + NPPmean + SSTsd

legal.model<- glmmTMB(response ~ no3_m + npp_mean + t_sd +
                    (1|Ecoregion/clust.4km) + (1|campaignid),
                  family = nbinom2(link = "log"),
                  data = legal,REML = TRUE)
r2_nakagawa(legal.model)

## Quick check via Dharma

setwd(plots.dir)
dir()

pdf("Legal_Choerodon_spp.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(legal.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(legal.model)

## Extract Pearson residuals

res1<-resid(legal.model, type = "pearson")

## Put the residuals in the dataset

legal$res<-res1

## Sublegal - (1) depth + NPPmean + SSTmean

sublegal.model<- glmmTMB(response ~ depth + npp_mean + t_m +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = sublegal,REML = TRUE)
r2_nakagawa(sublegal.model)

## Quick check via Dharma 

setwd(plots.dir)
dir()

pdf("Sub-legal_Choerodon_spp.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(sublegal.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(sublegal.model)

## Extract Pearson residuals

res1<-resid(sublegal.model, type = "pearson")

## Put the residuals in the dataset

sublegal$res<-res1

### Check for residual auto-correlation ----

## Convert into an spatial object

coordinates(large)<-~longitude+latitude
class(large)

coordinates(legal)<-~longitude+latitude
class(legal)

coordinates(sublegal)<-~longitude+latitude
class(sublegal)

### Now we can set the reference system to the widely used WGS84

proj4string(large)<-CRS("+init=epsg:4326")
large@proj4string

proj4string(legal)<-CRS("+init=epsg:4326")
legal@proj4string

proj4string(sublegal)<-CRS("+init=epsg:4326")
sublegal@proj4string

## Compute semivariogram

library(gstat)

large.var.response<-variogram(large$response~1,data=large,cutoff=20,width=0.1)
legal.var.response<-variogram(legal$response~1,data=legal,cutoff=20,width=0.1)
sublegal.var.response<-variogram(sublegal$response~1,data=legal,cutoff=20,width=0.1)
large.var.res<-variogram(large$res~1,data=large,cutoff=20,width=0.1)
legal.var.res<-variogram(legal$res~1,data=legal,cutoff=20,width=0.1)
sublegal.var.res<-variogram(sublegal$res~1,data=legal,cutoff=20,width=0.1)

## Export semivariograms

## Large response

plot(large.var.response,main="(a) Legal (response)",ylab="Semivariance",xlab="Distance (km)")
plot(large.var.res,main="(b) Legal (residuals)",ylab="Semivariance",xlab="Distance (km)")

## Legal response

plot(legal.var.response,main="(a) Legal (response)",ylab="Semivariance",xlab="Distance (km)")
plot(legal.var.res,main="(b) Legal (residuals)",ylab="Semivariance",xlab="Distance (km)")

## Sub-legal response

plot(sublegal.var.response,main="(c) Sub-legal (residuals-nested)",ylab="Semivariance",xlab="Distance (km)")
plot(sublegal.var.res,main="(c) Sub-legal (residuals)",ylab="Semivariance",xlab="Distance (km)",ylim=c(0.5,1.2))

### Calculate Moran's I based on distance via the spdep package
### We would calculated manually for each basin, diversity index combination and store it in a csv file
### However - for future reference it would be good to loop all of this together

## Global Moran's I ---- 

library(geosphere)

## Large

dists <- distm(large, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(large$response,dists.inv)

## Legal (res 1)

Moran.I(large$res,dists.inv) 

## Legal

dists <- distm(legal, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(legal$response,dists.inv)

## Legal (res 1)

Moran.I(legal$res,dists.inv) 

## Sub-legal

dists <- distm(sublegal, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(sublegal$response,dists.inv)

## Legal (res)

Moran.I(sublegal$res,dists.inv) 

#### (F) Nemadactylus spp ----

# Read in data 

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.species.group.csv")%>%
  filter(species.group%in%c('Nemadactylus spp.'))%>%
  filter(relief<200)%>%
  rename(response=number)%>%
  na.omit()%>%
  glimpse()
unique(dat$species.group)
unique(dat$Taxa)

## Check sampling effort per levels of random effect

test<-dat%>%
  group_by(Ecoregion)%>%
  summarise(campaign_N=n_distinct(campaignid),
            clust_N=n_distinct(clust.4km),
            N=n_distinct(id))%>%
  glimpse()

## Check distribution of gravity values

library(viridis)

ggplot(dat,aes(x=longitude,y=latitude,colour=log1p(gravity.50),size=log1p(gravity.50)))+
  scale_color_viridis()+
  geom_point()+
  theme_classic()

## Code covariates 

dat$slope<-as.numeric(dat$slope)
dat$relief<-as.numeric(dat$relief)
dat$ausbath<-as.numeric(dat$ausbath)
dat$depth<-as.numeric(dat$depth)
dat$t_m<-as.numeric(dat$t_m)
dat$t_sd<-as.numeric(dat$t_sd)
dat$no3_m<-as.numeric(dat$no3_m)
dat$no3_sd<-as.numeric(dat$no3_sd)
dat$po4_m<-as.numeric(dat$po4_m)
dat$po4_sd<-as.numeric(dat$po4_sd)
dat$npp_mean<-as.numeric(dat$npp_mean)
dat$npp_min<-as.numeric(dat$npp_min)
dat$npp_max<-as.numeric(dat$npp_max)
dat$npp_sd<-as.numeric(dat$npp_sd)
dat$dst2ramp<-as.numeric(dat$dst2ramp)
dat$pop.den.50<-as.numeric(dat$pop.den.50)
dat$gravity.50<-as.numeric(dat$gravity.50)
dat$status<-as.factor(dat$status)
dat$campaignid<-as.factor(dat$campaignid)
dat$state<-as.factor(dat$state)
dat$clust.4km<-as.factor(dat$clust.4km)

glimpse(dat)

## Set predictor variables

names(dat)
pred.vars.cont=c("relief","depth","t_m","t_sd","no3_m","npp_mean","npp_sd","gravity.50")

# Plot of likely transformations 

plot.new()

par(mfrow=c(3,2))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

## Transform predictors 

dat$gravity.50<-log1p(dat$gravity.50)
dat$no3_m<-log1p(dat$no3_m)
dat$relief<-log1p(dat$relief)

## Identify outliers 

plot.new()
par(mfrow=c(1,3))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  boxplot(x)
}

## Spatial plots of outliers

ggplot(dat,aes(x=longitude,y=latitude,colour=no3_m))+
  geom_point()+
  scale_color_viridis()+
  theme_classic()+
  Theme1

## Explore initial correlations between covariates

library(corrplot)
library(RColorBrewer)

## Visualize

M<-dat[,pred.vars.cont]
M<-round(cor(M),2)

plot.new()
par(mfrow=c(1,1))
corrplot(M,method="number")

## Scale (mean 0 and SD 1) predictors 

pred.vars.cont

dat$relief<-scale(dat$relief)
dat$depth<-scale(dat$depth)
dat$t_m<-scale(dat$t_m)
dat$t_sd<-scale(dat$t_sd)
dat$no3_m<-scale(dat$no3_m)
dat$npp_mean<-scale(dat$npp_mean)
dat$npp_sd<-scale(dat$npp_sd)
dat$gravity.50<-scale(dat$gravity.50)

## Create individual dataframes for modelling

plot.new()
par(mfrow=c(3,2))

## Large

large<-dat%>%
  dplyr::filter(Taxa=="Large")%>%
  glimpse()
hist(large$response)
plot(large$response)

sum(large$response == 0 ) / length(large$response) ## 96% zeroes

## Convert to presence/absence

large<-large%>%
  mutate(response=ifelse(response>0,1,0))%>%
  glimpse()

## Legal

legal<-dat%>%
  dplyr::filter(Taxa=="Legal")%>%
  glimpse()
hist(legal$response)
plot(legal$response)

sum(legal$response == 0 ) / length(legal$response) ## 82% zeroes

## Sub-legal

sublegal<-dat%>%
  dplyr::filter(Taxa=="Sub-legal")%>%
  glimpse()
hist(sublegal$response)
plot(sublegal$response)

sum(sublegal$response == 0 ) / length(sublegal$response) ## 92% zeroes

## Run models 

library(glmmTMB)

## (F1) Large ----

# that's better

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = binomial(link = "logit"),
                     data = large,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.5(Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column

sel.table$Model<-rownames(sel.table)

test<-get.models(results,subset = delta <= 2)

test$`325`$call$formula
test$`101`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Nemadactylus_spp._Large_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"GLM_Nemadactylus_Importance_Large_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## (F2) Legal ----

## Fit a full model 

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth  + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = legal,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.5(Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column
sel.table$Model<-rownames(sel.table)

test<-get.models(results,subset = delta <= 2)

test$`38`$call$formula
test$`294`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Plectropomus_Legal_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"GLM_Plectropomus_Importance_Legal_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## (F3) Sub-legal ----

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth  + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = sublegal,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.28 (Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column
sel.table$Model<-rownames(sel.table)
sel.table

## Get models 

test<-get.models(results,subset = delta <= 2)

## Write variable names

test$`262`$call$formula
test$`322`$call$formula
test$`386`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Plectropomus_Sublegal_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"Plectropomus_Sublegal_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## Top-ranked models and model validation ----

## Large 

large.model<- glmmTMB(response ~ no3_m + status + t_sd +
                        (1|Ecoregion/clust.4km) + (1|campaignid),
                      family = binomial(link = "logit"),
                      data = large,REML = TRUE)
r2_nakagawa(large.model)

## Quick check via Dharma

setwd(plots.dir)
dir()

pdf("Large_Nemadactylus_spp.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(large.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(large.model)

## Extract Pearson residuals

res1<-resid(large.model, type = "pearson")

## Put the residuals in the dataset

large$res<-res1

## Legal - (1) Nitrate + NPPmean + SSTsd

legal.model<- glmmTMB(response ~depth + no3_m + relief +
                    (1|Ecoregion/clust.4km) + (1|campaignid),
                  family = nbinom2(link = "log"),
                  data = legal,REML = TRUE)
r2_nakagawa(legal.model)

## Quick check via Dharma

setwd(plots.dir)
dir()

pdf("Legal_Nemadactylus_spp.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(legal.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(legal.model)

## Extract Pearson residuals

res1<-resid(legal.model, type = "pearson")

## Put the residuals in the dataset

legal$res<-res1

## Sublegal - (1) depth + SSTsd + Nitrate

sublegal.model<- glmmTMB(response ~ depth + no3_m + t_sd +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = sublegal,REML = TRUE)
r2_nakagawa(sublegal.model)

## Quick check via Dharma 

setwd(plots.dir)
dir()

pdf("Sub-legal_Nemadactylus_spp.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(sublegal.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(sublegal.model)

## Extract Pearson residuals

res1<-resid(sublegal.model, type = "pearson")

## Put the residuals in the dataset

sublegal$res<-res1

### Check for residual auto-correlation ----

## Convert into an spatial object

coordinates(large)<-~longitude+latitude
class(large)

coordinates(legal)<-~longitude+latitude
class(legal)

coordinates(sublegal)<-~longitude+latitude
class(sublegal)

### Now we can set the reference system to the widely used WGS84

proj4string(large)<-CRS("+init=epsg:4326")
large@proj4string

proj4string(legal)<-CRS("+init=epsg:4326")
legal@proj4string

proj4string(sublegal)<-CRS("+init=epsg:4326")
sublegal@proj4string

## Compute semivariogram

library(gstat)

large.var.response<-variogram(large$response~1,data=large,cutoff=20,width=0.1)
legal.var.response<-variogram(legal$response~1,data=legal,cutoff=20,width=0.1)
sublegal.var.response<-variogram(sublegal$response~1,data=legal,cutoff=20,width=0.1)
large.var.res<-variogram(large$res~1,data=large,cutoff=20,width=0.1)
legal.var.res<-variogram(legal$res~1,data=legal,cutoff=20,width=0.1)
sublegal.var.res<-variogram(sublegal$res~1,data=legal,cutoff=20,width=0.1)

## Export semivariograms

## Large response

plot(large.var.response,main="(a) Legal (response)",ylab="Semivariance",xlab="Distance (km)")
plot(large.var.res,main="(b) Legal (residuals)",ylab="Semivariance",xlab="Distance (km)")

## Legal response

plot(legal.var.response,main="(a) Legal (response)",ylab="Semivariance",xlab="Distance (km)")
plot(legal.var.res,main="(b) Legal (residuals)",ylab="Semivariance",xlab="Distance (km)")

## Sub-legal response

plot(sublegal.var.response,main="(c) Sub-legal (residuals-nested)",ylab="Semivariance",xlab="Distance (km)")
plot(sublegal.var.res,main="(c) Sub-legal (residuals)",ylab="Semivariance",xlab="Distance (km)",ylim=c(0.5,1.2))

### Calculate Moran's I based on distance via the spdep package
### We would calculated manually for each basin, diversity index combination and store it in a csv file
### However - for future reference it would be good to loop all of this together

## Global Moran's I ---- 

library(geosphere)

## Large

dists <- distm(large, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(large$response,dists.inv)

## Legal (res 1)

Moran.I(large$res,dists.inv) 

## Legal

dists <- distm(legal, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(legal$response,dists.inv)

## Legal (res 1)

Moran.I(legal$res,dists.inv) 

## Sub-legal

dists <- distm(sublegal, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(sublegal$response,dists.inv)

## Legal (res)

Moran.I(sublegal$res,dists.inv) 

#### (G) Notolabrus spp ----

# Read in data 

setwd(data.dir)
dir()

dat<-read.csv("length.analysis.species.group.csv")%>%
  filter(species.group%in%c('Notolabrus spp.'))%>%
  filter(relief<200)%>%
  rename(response=number)%>%
  na.omit()%>%
  glimpse()
unique(dat$species.group)
unique(dat$Taxa)

## Check sampling effort per levels of random effect

test<-dat%>%
  group_by(Ecoregion)%>%
  summarise(campaign_N=n_distinct(campaignid),
            clust_N=n_distinct(clust.4km),
            N=n_distinct(id))%>%
  glimpse()

## Check distribution of gravity values

library(viridis)

ggplot(dat,aes(x=longitude,y=latitude,colour=log1p(gravity.50),size=log1p(gravity.50)))+
  scale_color_viridis()+
  geom_point()+
  theme_classic()

## Code covariates 

dat$slope<-as.numeric(dat$slope)
dat$relief<-as.numeric(dat$relief)
dat$ausbath<-as.numeric(dat$ausbath)
dat$depth<-as.numeric(dat$depth)
dat$t_m<-as.numeric(dat$t_m)
dat$t_sd<-as.numeric(dat$t_sd)
dat$no3_m<-as.numeric(dat$no3_m)
dat$no3_sd<-as.numeric(dat$no3_sd)
dat$po4_m<-as.numeric(dat$po4_m)
dat$po4_sd<-as.numeric(dat$po4_sd)
dat$npp_mean<-as.numeric(dat$npp_mean)
dat$npp_min<-as.numeric(dat$npp_min)
dat$npp_max<-as.numeric(dat$npp_max)
dat$npp_sd<-as.numeric(dat$npp_sd)
dat$dst2ramp<-as.numeric(dat$dst2ramp)
dat$pop.den.50<-as.numeric(dat$pop.den.50)
dat$gravity.50<-as.numeric(dat$gravity.50)
dat$status<-as.factor(dat$status)
dat$campaignid<-as.factor(dat$campaignid)
dat$state<-as.factor(dat$state)
dat$clust.4km<-as.factor(dat$clust.4km)

glimpse(dat)

## Set predictor variables

names(dat)
pred.vars.cont=c("relief","depth","t_m","t_sd","no3_m","npp_mean","npp_sd","gravity.50")

# Plot of likely transformations 

plot.new()

par(mfrow=c(3,2))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  hist(sqrt(x))
  plot(sqrt(x))
  hist(log(x+1))
  plot(log(x+1))
}

## Transform predictors 

dat$gravity.50<-log1p(dat$gravity.50)
dat$relief<-log1p(dat$relief)

## Identify outliers 

plot.new()
par(mfrow=c(1,3))
for (i in pred.vars.cont) {
  x<-dat[ ,i]
  x = as.numeric(unlist(x))  
  hist((x))#Looks best
  plot((x),main = paste(i))
  boxplot(x)
}

## Spatial plots of outliers

ggplot(dat,aes(x=longitude,y=latitude,colour=no3_m))+
  geom_point()+
  scale_color_viridis()+
  theme_classic()+
  Theme1

## Explore initial correlations between covariates

library(corrplot)
library(RColorBrewer)

## Visualize

M<-dat[,pred.vars.cont]
M<-round(cor(M),2)

plot.new()
par(mfrow=c(1,1))
corrplot(M,method="number")

## Scale (mean 0 and SD 1) predictors 

pred.vars.cont

dat$relief<-scale(dat$relief)
dat$depth<-scale(dat$depth)
dat$t_m<-scale(dat$t_m)
dat$t_sd<-scale(dat$t_sd)
dat$no3_m<-scale(dat$no3_m)
dat$npp_mean<-scale(dat$npp_mean)
dat$npp_sd<-scale(dat$npp_sd)
dat$gravity.50<-scale(dat$gravity.50)

## Create individual dataframes for modelling

plot.new()
par(mfrow=c(3,2))

## Large

large<-dat%>%
  dplyr::filter(Taxa=="Large")%>%
  glimpse()
hist(large$response)
plot(large$response)

sum(large$response == 0 ) / length(large$response) ## 85% zeroes

## Convert to presence/absence

large<-large%>%
  mutate(response=ifelse(response>0,1,0))%>%
  glimpse()

## Legal

legal<-dat%>%
  dplyr::filter(Taxa=="Legal")%>%
  glimpse()
hist(legal$response)
plot(legal$response)

sum(legal$response == 0 ) / length(legal$response) ## 50% zeroes

## Sub-legal

sublegal<-dat%>%
  dplyr::filter(Taxa=="Sub-legal")%>%
  glimpse()
hist(sublegal$response)
plot(sublegal$response)

sum(sublegal$response == 0 ) / length(sublegal$response) ## 58% zeroes

## Run models 

library(glmmTMB)

## (G1) Large ----

# that's better

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = binomial(link = "logit"),
                     data = large,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.5(Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column

sel.table$Model<-rownames(sel.table)

test<-get.models(results,subset = delta <= 2)

test$`449`$call$formula
test$`453`$call$formula
test$`481`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Notolabrus_spp._Large_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"GLM_Notolabrus_Importance_Large_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## (G2) Legal ----

## Fit a full model 

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth  + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = legal,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.5(Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column
sel.table$Model<-rownames(sel.table)

test<-get.models(results,subset = delta <= 2)

test$`99`$call$formula
test$`353`$call$formula
test$`101`$call$formula
test$`225`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Plectropomus_Legal_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"GLM_Plectropomus_Importance_Legal_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## (G3) Sub-legal ----

full.model<- glmmTMB(response ~ relief + t_m + t_sd + no3_m + depth  + npp_mean + npp_sd +
                       gravity.50 + gravity.50:status + status +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = sublegal,REML = TRUE)
summary(full.model)

# Exclude Variables from the same model that have correlations > 0.28 (Graham, 2003)

smat <- abs(cor(legal[,pred.vars.cont])) <= .5
smat[!lower.tri(smat)] <- NA

# only allow a maximum of 4 and minimum of 1 parameters in each model

results<-dredge(full.model,m.lim = c(3,4),subset = smat,rank = AIC)

# what's it look like, hmm AIC with small sample bias adjustment AICc
# delta AICc, and the model weights
results

# create a confidence set of models using the subset function
# select models with delta AICc less than 5
# IMPORTANT: Weights have been renormalized!!

subset(results, delta <=2)

# coerce the object out.put into a data frame
# elements 6-10 in out.put have what we want

sel.table<-as.data.frame(subset(results, delta <=2))[13:17]
sel.table

# a little clean-up, lets round things a bit

sel.table[,2:3]<- round(sel.table[,2:3],2)
sel.table[,4:5]<- round(sel.table[,4:5],3)

## lets be sure to put the model names in a column
sel.table$Model<-rownames(sel.table)
sel.table

## Get models 

test<-get.models(results,subset = delta <= 2)

## Write variable names

test$`166`$call$formula
test$`162`$call$formula
test$`226`$call$formula

# write to a file, here a comma separated values format
# make sure your working directory is properly specified

setwd(model.out)
dir()

write.csv(sel.table,"GLM_Plectropomus_Sublegal_AIC.csv", row.names = F)

# Importance weights for individual predictor variables
# calculated using the importance function

MuMIn::importance(results) ## Error

importance.scores<-as.data.frame(MuMIn::importance(results))
importance.scores[,1]<- round(importance.scores[,1],2)

importance.scores$Predictor<-rownames(importance.scores)

## Export importance scores

setwd(model.out)
dir()

write.csv(importance.scores,"Plectropomus_Sublegal_AIC.csv", row.names = F)

## Obtain Model averaged coefficients and SE

# Model average using all candidate models, always use revised.var = TRUE

MA.ests<-model.avg(results, revised.var = TRUE)
summary(MA.ests)

## Top-ranked models and model validation ----

## Large 

large.model<- glmmTMB(response ~ status + t_m + t_sd +
                        (1|Ecoregion/clust.4km) + (1|campaignid),
                      family = binomial(link = "logit"),
                      data = large,REML = TRUE)
r2_nakagawa(large.model)

## Quick check via Dharma

setwd(plots.dir)
dir()

pdf("Large_Notolabrus_spp.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(large.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(large.model)

## Extract Pearson residuals

res1<-resid(large.model, type = "pearson")

## Put the residuals in the dataset

large$res<-res1

## Legal - (1) gravity + relief + status

legal.model<- glmmTMB(response ~ gravity.50 + relief + status +
                  (1|Ecoregion/clust.4km) + (1|campaignid),
                  family = nbinom2(link = "log"),
                  data = legal,REML = TRUE)
r2_nakagawa(legal.model)

## Quick check via Dharma

setwd(plots.dir)
dir()

pdf("Legal_Notolabrus_spp.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(legal.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(legal.model)

## Extract Pearson residuals

res1<-resid(legal.model, type = "pearson")

## Put the residuals in the dataset

legal$res<-res1

## Sublegal - (1) depth + nitrate + relief + SSTmean

sublegal.model<- glmmTMB(response ~ depth + no3_m + relief + t_m +
                       (1|Ecoregion/clust.4km) + (1|campaignid),
                     family = nbinom2(link = "log"),
                     data = sublegal,REML = TRUE)
r2_nakagawa(sublegal.model)

## Quick check via Dharma 

setwd(plots.dir)
dir()

pdf("Sub-legal_Notolabrus_spp.pdf",width = 7,height = 7,useDingbats=FALSE) 

M1_res <- simulateResiduals(sublegal.model)
plot(M1_res)

dev.off()

## Test zero-inflation

testZeroInflation(sublegal.model)

## Extract Pearson residuals

res1<-resid(sublegal.model, type = "pearson")

## Put the residuals in the dataset

sublegal$res<-res1

### Check for residual auto-correlation ----

## Convert into an spatial object

coordinates(large)<-~longitude+latitude
class(large)

coordinates(legal)<-~longitude+latitude
class(legal)

coordinates(sublegal)<-~longitude+latitude
class(sublegal)

### Now we can set the reference system to the widely used WGS84

proj4string(large)<-CRS("+init=epsg:4326")
large@proj4string

proj4string(legal)<-CRS("+init=epsg:4326")
legal@proj4string

proj4string(sublegal)<-CRS("+init=epsg:4326")
sublegal@proj4string

## Compute semivariogram

library(gstat)

large.var.response<-variogram(large$response~1,data=large,cutoff=20,width=0.1)
legal.var.response<-variogram(legal$response~1,data=legal,cutoff=20,width=0.1)
sublegal.var.response<-variogram(sublegal$response~1,data=legal,cutoff=20,width=0.1)
large.var.res<-variogram(large$res~1,data=large,cutoff=20,width=0.1)
legal.var.res<-variogram(legal$res~1,data=legal,cutoff=20,width=0.1)
sublegal.var.res<-variogram(sublegal$res~1,data=legal,cutoff=20,width=0.1)

## Export semivariograms

## Large response

plot(large.var.response,main="(a) Legal (response)",ylab="Semivariance",xlab="Distance (km)")
plot(large.var.res,main="(b) Legal (residuals)",ylab="Semivariance",xlab="Distance (km)")

## Legal response

plot(legal.var.response,main="(a) Legal (response)",ylab="Semivariance",xlab="Distance (km)")
plot(legal.var.res,main="(b) Legal (residuals)",ylab="Semivariance",xlab="Distance (km)")

## Sub-legal response

plot(sublegal.var.response,main="(c) Sub-legal (residuals-nested)",ylab="Semivariance",xlab="Distance (km)")
plot(sublegal.var.res,main="(c) Sub-legal (residuals)",ylab="Semivariance",xlab="Distance (km)",ylim=c(0.5,1.2))

### Calculate Moran's I based on distance via the spdep package
### We would calculated manually for each basin, diversity index combination and store it in a csv file
### However - for future reference it would be good to loop all of this together

## Global Moran's I ---- 

library(geosphere)

## Large

dists <- distm(large, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(large$response,dists.inv)

## Legal (res 1)

Moran.I(large$res,dists.inv) 

## Legal

dists <- distm(legal, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(legal$response,dists.inv)

## Legal (res 1)

Moran.I(legal$res,dists.inv) 

## Sub-legal

dists <- distm(sublegal, fun = distGeo)
min(dists) ## Distance of 0 correspond to temporally auto-correlated data
max(dists) ## 5398.58

## We need to generate a matrix of inverse distance weights
### Replace the diagonal with 0

dists.inv<- 1/dists
diag(dists.inv) <- 0

## Remove infinite values in the matrix

# dists.inv[is.infinite(dists.inv)] <- 0 ## I am not too sure this step in correct

## Calculate Moran's I test statistic - ape package

library(ape)

## Legal (response)

Moran.I(sublegal$response,dists.inv)

## Legal (res)

Moran.I(sublegal$res,dists.inv) 

#########################################################################################################
##################### (7) Model averaged coefficients - Regional Species Groups ########################

library(forcats)

setwd('C:/Users/22373243/Dropbox/Projects/Analysis/Analysis_GlobalArchive_Bodysize_Nestor/ModelOut/Species_models_20200811')
dir()

dat<-read.csv("Model_Averaged_Coefficients_All_1.csv")%>%
  dplyr::select(-X)%>%
  dplyr::mutate(Size=fct_relevel(Size,'Sub-legal','Legal','Large'))%>%
  glimpse()

## Gravity

coefficients<-ggplot()+
  xlab('')+
  ggtitle('(a)')+
  geom_point(data=dat,aes(x=Taxa,y=Estimate,shape=Size,colour=Size),size=3,show.legend = T,
             position=position_dodge(width=0.5))+
  geom_errorbar(data=dat,aes(x=Taxa,ymin=Estimate-SE, ymax=Estimate+SE,colour=Size), width=.2,show.legend = F,
                position=position_dodge(width=0.5))+
  scale_color_manual(labels = c("Sub-legal","Legal","Large"),values=c("lightgrey","darkgrey","black"))+
  coord_flip()+
  facet_wrap(~Predictor,scales = "free_x")+
  geom_hline(yintercept = 0,linetype="dashed")+
  theme_classic()+
  Theme1+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
coefficients

## Export plot

setwd(plots.dir)
dir()

name<-'Fig.3.Model average coefficients_20200404'

ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm")
ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",useDingbats=FALSE)

########################################################################################################
############################### (8) Summed Importance scores (AICweights) ##############################

## Bring in summed AICc weights for anthropogenic, habitat, and environmental predictors

setwd('C:/Users/22373243/Dropbox/Projects/Analysis/Analysis_GlobalArchive_Bodysize_Nestor/ModelOut/Datasets_Plotting')
dir()

dat<-read.csv("Summed_AICc_All.csv")%>%
  mutate(Taxa.2=paste(Model,Taxa,sep = "."))%>%
  glimpse()

## Plot summed AICweights

unique(dat$Model)
unique(dat$Taxa.2)

dat<-dat%>%
  dplyr::mutate(Taxa.2=fct_relevel(Taxa.2,"Choerodon spp..Sub-legal","Choerodon spp..Legal","Choerodon spp..Large",
                                   "Chrysophrys auratus.Sub-legal","Chrysophrys auratus.Legal","Chrysophrys auratus.Large",
                                   "Lethrinus spp..Sub-legal","Lethrinus spp..Legal","Lethrinus spp..Large",
                                   "Lutjanus spp..Sub-legal","Lutjanus spp..Legal","Lutjanus spp..Large",
                                   "Nemadactylus spp..Sub-legal","Nemadactylus spp..Legal","Nemadactylus spp..Large",
                                   "Notolabrus spp..Sub-legal","Notolabrus spp..Legal","Notolabrus spp..Large",
                                   "Plectropomus spp..Sub-legal","Plectropomus spp..Legal","Plectropomus spp..Large",
                                   "Assemblage.Sub-legal","Assemblage.Legal","Assemblage.Large"))%>%
  glimpse()


importance<-ggplot(data=dat,aes(x=Taxa.2,y=Importance,fill=Taxa))+
  geom_bar(stat = "identity",colour="black")+
  scale_fill_manual(labels = c("Large","Legal","Sub-legal"),values=c("black","darkgrey","lightgrey"))+
  facet_wrap(~Predictor)+
  coord_flip()+
  theme_classic()+
  Theme1+
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
importance

## Export plot

setwd(plots.dir)
dir()

name<-'Fig.4.Summed_AICweights_20200406'

ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm")
ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",useDingbats=FALSE)

##########################################################################################################
############################# (9) Extended Data and Figures ##############################################

### Set plotting defaults

# Read in data----

setwd(data.dir)
dir()

dat<-read.csv("legal.analysis.20200722_Snapshot.csv")%>%
  dplyr::select(-X)%>%
  dplyr::filter(!campaignid%in%c('2016-09_Scott & Rowleys_FinPrint_stereoBRUVS'))%>%
  dplyr::filter(depth<=50)%>%
  dplyr::rename(response=number)%>%
  na.omit()%>%
  glimpse()

n_distinct(dat$Taxa)
unique(dat$Taxa)
n_distinct(dat$campaignid)
unique(dat$campaignid)
table(dat$status)
n_distinct(dat$id)
unique(dat$state)
hist(dat$depth)
hist(dat$ausbath)

### Check NAs in the data

sum(is.na(dat))/prod(dim(dat))*100
apply(dat,2,function(col)sum(is.na(col))/length(col))*100

## Filter extreme value in the relief covariate

dat<-dat%>%
  filter(relief<200)%>%
  filter(dst2shelf<400000)%>%
  glimpse()

ggplot(dat,aes(x=longitude,y=latitude))+
  geom_point()+
  theme_classic()

## Transform predictors 

dat$gravity.census<-log1p(dat$gravity.census)
dat$wave_mean<-log1p(dat$wave_mean)
dat$npp_mean<-log1p(dat$npp_mean) ## There is outliers we will need to deal with later
dat$no3_m<-log1p(dat$no3_m)
dat$relief<-log1p(dat$relief)

#### Map of spatial context of the study ----

### Bring in Marine Ecoregions of the world

setwd("C:/Users/22373243/Dropbox/Nestor_ascii_2/MEOW")
dir()

marine.ecoregions<-readOGR('.','meow_ecos')
proj4string(marine.ecoregions)
marine.ecoregions<-crop(marine.ecoregions,extent(105,160,-48,-9))
plot(marine.ecoregions)
marine.ecoregions<-fortify(marine.ecoregions,region = "Lat_Zone")

## Filter out Ecoregions not wihtin limits of the Continental Platform of Australia

marine.ecoregions<-marine.ecoregions%>%
  filter(!grepl('Java Transitional|Tropical Southwestern Pacific|Lord Howe and Norfolk Islands|Eastern Coral Triangle|
|Western Coral Triangle',id))%>%
  glimpse()
unique(marine.ecoregions$id)

### Bring in state shapefile

setwd('C:/Users/22373243/Dropbox/Nestor_ascii_2')
dir()

states<-readOGR('.','Australia_Polygon')
proj4string(states)
states<-crop(states,extent(105,160,-48,-9))
bbox(states)
plot(states)

### Explore the polygon data

test<-states@data%>%
  glimpse()
unique(test$name)

### Convert to a dataframe

states<-fortify(states,region='name')
unique(states$id)

### Retain only main states with coastline

states<-states%>%
  filter(grepl('New South Wales|Northern Territory|Queensland|South Australia|
               |Tasmania|Victoria|Western Australia',id))%>%
  glimpse()

### Quick Ggplot on spatial distribution of samples at a national scale

col_list<-c("blue4","red")

map<-ggplot()+
  geom_polygon(data=marine.ecoregions,aes(x=long,y=lat,group=group,fill=id),color="gray",size=0.5,alpha=0.1)+
  scale_fill_manual(values=col_list)+
  geom_polygon(data=states,aes(x=long,y=lat,group=group),color="black",fill="gray20",size=0.5)+
  #geom_polygon(data=mpa,aes(x=long,y=lat,group=group),fill="green")+
  geom_point(data=dat%>%filter(status%in%c('Fished')),aes(x=longitude,y=latitude),colour="red",show.legend = F,size=2,alpha=0.5)+
  geom_point(data=dat%>%filter(status%in%c('No-take')),aes(x=longitude,y=latitude),colour="blue3",show.legend = F,size=2,alpha=1)+
  geom_jitter(width = 1)+
  xlab("Latitude")+
  ylab("Longitude")+
  xlim(105,160)+
  ylim(-48,-8)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  #labs(fill="Ecoregions")+
  annotate("text",x=122,y=-25,label="WA",colour = "white",size=7)+
  annotate("text",x=133.5,y=-20,label="NT",colour = "white",size=7)+
  annotate("text",x=144,y=-22,label="QLD",colour = "white",size=7)+
  annotate("text",x=135,y=-29,label="SA",colour = "white",size=7)+
  annotate("text",x=147,y=-32,label="NSW",colour = "white",size=7)+
  annotate("text",x=144,y=-37,label="VIC",colour = "white",size=6)+
  annotate("text",x=146.7,y=-42,label="TAS",colour = "white",size=5)+
  #annotate("text",x=139,y=-13,label="1",colour = "black",size=7)+
  #annotate("text",x=119,y=-18,label="2",colour = "black",size=7)+
  #annotate("text",x=112,y=-28,label="3",colour = "black",size=7)+
  #annotate("text",x=129,y=-33.5,label="4",colour = "black",size=7)+
  #annotate("text",x=141,y=-40,label="5",colour = "black",size=7)+
  #annotate("text",x=155,y=-30,label="6",colour = "black",size=7)+
  #annotate("text",x=152,y=-22,label="7",colour = "black",size=7)+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black",fill = NA,size = 1),
        legend.title = element_blank(),
        legend.text=element_text(size=14),
        text = element_text(size=14),
        legend.position = c(0.3,0.1),
        legend.direction = "horizontal",
        legend.spacing.x = unit(0.5, 'cm'))
map

### Export plot

setwd(plots.dir)
dir()

name<-'Spatial_context'

ggsave(paste(name,".png",sep="."),width = 21, height = 15,units = "cm")
ggsave(paste(name,".pdf",sep="."),width = 21, height = 15,units = "cm",useDingbats=FALSE)


#### Table summary of Dataset ----

table(dat$status)
n_distinct(dat$id) ## 3,613 samples
sum(dat$response) ## 20,091
unique(dat$state)

Realm.summary<-dat%>%
  dplyr::group_by(Realm,Taxa)%>%
  dplyr::summarise(N=n_distinct(id),
                   prop=N/3613,
                   Ind=sum(response),
                   Ind_prop=Ind/20091,
                   Depth_mean=mean(depth,na.rm=T),
                   Depth_sd=sd(depth,na.rm = T),
                   Depth_min=min(depth),
                   Depth_max=max(depth))%>%
  glimpse()

State.summary<-dat%>%
  dplyr::group_by(state,Taxa)%>%
  dplyr::summarise(N=n_distinct(id),
                   prop=N/3613,
                   Ind=sum(response),
                   Ind_prop=Ind/20091,
                   Depth_mean=mean(depth,na.rm=T),
                   Depth_sd=sd(depth,na.rm = T),
                   Depth_min=min(depth),
                   Depth_max=max(depth))%>%
  glimpse()


#### Gravity + Reserve size distribution ----

## Bring in Australia coastline

setwd('C:/Users/22373243/Dropbox/Nestor_ascii_2')
dir()

australia<-readOGR('.','Australiaboundary67')
proj4string(australia)
australia <- spTransform(australia, CRS("+proj=longlat +datum=WGS84"))

## Crop australia to the extent

e<-extent(112,155,-45,-7)

australia<-crop(australia,e)
plot(australia)

australia.continental<-fortify(australia)

## Bring in Gravity raster layer and convert into dataframe to plot with geom_tile (ggplot2)

plot.new()
par(mfrow=(c(1,1)))

gravity<-raster("r.gravity.50.log.tif")
gravity
proj4string(gravity)
coordinates(gravity)
plot(gravity)

gravity<-as.data.frame(gravity,xy=T,na.rm=TRUE)%>%
  glimpse()
colnames(gravity)<-c('longitude','latitude','Gravity')
glimpse(gravity)

## Bring in MPA dataset from CAPAD

setwd('C:/Users/22373243/Dropbox/Nestor_ascii_2/Nestor_ascii')
dir()

mpa<-readOGR('.','capad_2018_marine_edits_WA_FINAL_v2_WGS84')
proj4string(mpa)
mpa <- spTransform(mpa, CRS("+proj=longlat +datum=WGS84"))

mpa<-fortify(mpa)

#### Plot distribution of gravity and MPA across Australia continental shelf systems

gravity.plot<-ggplot()+
  geom_tile(data=gravity,aes(x=longitude,y=latitude,fill=Gravity),alpha=0.8)+
  scale_fill_viridis(option = "magma",direction = -1)+
  geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="gray12",alpha=0.8)+
  geom_polygon(data=mpa,aes(x=long,y=lat,group=group),fill=NA,colour="darkgreen",alpha=0.8)+
  scale_color_gradient()+
  scale_x_continuous(expand = c(0,0),limits = c(112,155))+
  scale_y_continuous(expand=c(0,0),limits=c(-45,-7))+
  ggtitle('(a)')+
  theme_classic()+
  Theme1+
  theme(legend.position = "none",
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.direction = "horizontal",
        legend.title = element_text(size=10),
        legend.text = element_text(size=10))
gravity.plot

### Export plot

setwd(plots.dir)
dir()

name<-'National'

ggsave(paste(name,"gravity.50.png",sep="."),width = 21, height = 15,units = "cm")
ggsave(paste(name,"gravity.50.pdf",sep="."),width = 21, height = 15,units = "cm",useDingbats=FALSE)

#### Create factor for reserve typ -
### < 1 km2, 1 - 10 (km2), 10 - 100 (km2), > 100 km2

ntr<-dat%>%
  dplyr::filter(Taxa=='Legal')%>%
  dplyr::mutate(reserve.type=cut(reserve.si,breaks = c(-Inf,0,10,100,Inf),labels = c("Fished","<10km2","10-100km2",">100km2")))%>%
  glimpse()
range(ntr$reserve.si)
table(ntr$reserve.type)

ntr$gravity.census<-log1p(ntr$gravity.census)
hist(ntr$gravity.census)
ntr$gravity.50<-log1p(ntr$gravity.50)
hist(ntr$gravity.50)

### Plot reserve size

## Desity plots ----

# unique(ntr$state)
# range(ntr$gravity.census)
# unique(ntr$reserve.type)
# 
# WA.plot<-ggplot(data=ntr%>%filter(state=='WA'),aes(x=gravity.census,fill=reserve.type))+
#   geom_density(alpha=0.3)+
#   scale_fill_manual(labels = c("Fished", "<10km2","10-100km2",">100km2"),
#                     values=c("red","darkmagenta","orange","darkgreen"))+
#   theme_classic()+
#   xlab('Gravity')+
#   ylab('Density of data')+
#   ggtitle('(b) WA')+
#   scale_x_continuous(expand = c(0,0),limits=c(0,13))+
#   scale_y_continuous(expand=c(0,0))+
#   Theme1+
#   theme(legend.direction = "horizontal",
#         legend.position = "bottom")
# WA.plot
# 
# VIC.plot<-ggplot(data=ntr%>%filter(state=='VIC'),aes(x=gravity.census,fill=reserve.type))+
#   geom_density(alpha=0.3)+
#   scale_fill_manual(labels = c("Fished", "<10km2","10-100km2",">100km2"),
#                     values=c("red","darkmagenta","orange","darkgreen"))+
#   theme_classic()+
#   xlab('Gravity')+
#   ylab('')+
#   ggtitle('(c) VIC')+
#   scale_x_continuous(expand = c(0,0),limits=c(0,13))+
#   scale_y_continuous(expand=c(0,0))+
#   Theme1+
#   theme(legend.direction = "horizontal",
#         legend.position = "bottom")
# VIC.plot
# 
# NSW.plot<-ggplot(data=ntr%>%filter(state=='NSW'),aes(x=gravity.census,fill=reserve.type))+
#   geom_density(alpha=0.3)+
#   scale_fill_manual(labels = c("Fished", "<10km2","10-100km2",">100km2"),
#                     values=c("red","darkmagenta","orange","darkgreen"))+
#   theme_classic()+
#   xlab('Gravity')+
#   ylab('')+
#   ggtitle('(d) NSW')+
#   scale_x_continuous(expand = c(0,0),limits=c(0,13))+
#   scale_y_continuous(expand=c(0,0))+
#   Theme1+
#   theme(legend.direction = "horizontal",
#         legend.position = "bottom")
# NSW.plot
# 
# SA.plot<-ggplot(data=ntr%>%filter(state=='SA'),aes(x=gravity.census,fill=reserve.type))+
#   geom_density(alpha=0.3)+
#   scale_fill_manual(labels = c("Fished", "<10km2","10-100km2",">100km2"),
#                     values=c("red","darkmagenta","orange","darkgreen"))+
#   theme_classic()+
#   xlab('Gravity')+
#   ylab('')+
#   ggtitle('(e) SA')+
#   scale_x_continuous(expand = c(0,0),limits=c(0,13))+
#   scale_y_continuous(expand=c(0,0))+
#   Theme1+
#   theme(legend.direction = "horizontal",
#         legend.position = "bottom")
# SA.plot
# 
# ggplot(data=ntr%>%filter(state=='TAS')) + aes(x = gravity.census, y = ..count../sum(..count..),fill=reserve.type) + 
#   geom_histogram(binwidth = 0.2)+
#   scale_fill_manual(labels = c("Fished", "<10km2","10-100km2",">100km2"),
#                     values=c("red","darkmagenta","orange","darkgreen"))+
#   theme_classic()+
#   xlab('Gravity')+
#   ylab('')+
#   ggtitle('(f) TAS')+
#   scale_x_continuous(expand = c(0,0),limits=c(0,13))+
#   scale_y_continuous(expand=c(0,0),labels = scales::percent)+
#   Theme1+
#   theme(legend.direction = "horizontal",
#         legend.position = "bottom")
# 
# 
# 
# TAS.plot<-ggplot(data=ntr%>%filter(state=='TAS'),aes(x=gravity.census,fill=reserve.type))+
#   geom_density(alpha=0.3)+
#   scale_fill_manual(labels = c("Fished", "<10km2","10-100km2",">100km2"),
#                     values=c("red","darkmagenta","orange","darkgreen"))+
#   theme_classic()+
#   xlab('Gravity')+
#   ylab('')+
#   ggtitle('(f) TAS')+
#   scale_x_continuous(expand = c(0,0),limits=c(0,13))+
#   scale_y_continuous(expand=c(0,0),labels = scales::percent)+
#   Theme1+
#   theme(legend.direction = "horizontal",
#         legend.position = "bottom")
# TAS.plot
# 
# ### Arrage plots
# 
# ggarrange(WA.plot,VIC.plot,NSW.plot,SA.plot,TAS.plot,
#           ncol=5,nrow = 1,align = 'hv',
#           common.legend = T,legend = "bottom")

## Histograms ----

unique(ntr$state)
range(ntr$gravity.census)
unique(ntr$reserve.type)
max(ntr$gravity.50)

WA.plot<-ggplot(data=ntr%>%filter(state=='WA')) + aes(x = gravity.50, y = ..count../sum(..count..),fill=reserve.type) +
  geom_histogram(binwidth = 1,alpha=0.3,colour="black")+
  scale_fill_manual(labels = c("Fished", "<10km2","10-100km2",">100km2"),
                    values=c("red","darkmagenta","orange","darkgreen"))+
  theme_classic()+
  xlab('Gravity')+
  ylab('Frequency (%)')+
  ggtitle('(b) WA')+
  scale_x_continuous(expand = c(0,0),limits=c(0,14))+
  scale_y_continuous(expand=c(0,0),labels = scales::percent)+
  Theme1+
  theme(legend.direction = "horizontal",
        legend.position = "bottom")
WA.plot

VIC.plot<-ggplot(data=ntr%>%filter(state=='VIC')) + aes(x = gravity.50, y = ..count../sum(..count..),fill=reserve.type) +
  geom_histogram(binwidth = 1,alpha=0.3,colour="black")+
  scale_fill_manual(labels = c("Fished", "<10km2","10-100km2",">100km2"),
                    values=c("red","darkmagenta","orange","darkgreen"))+
  theme_classic()+
  xlab('Gravity')+
  ylab('')+
  ggtitle('(c) VIC')+
  scale_x_continuous(expand = c(0,0),limits=c(0,14))+
  scale_y_continuous(expand=c(0,0),labels = scales::percent)+
  Theme1+
  theme(legend.direction = "horizontal",
        legend.position = "bottom")
VIC.plot

NSW.plot<-ggplot(data=ntr%>%filter(state=='NSW')) + aes(x = gravity.50, y = ..count../sum(..count..),fill=reserve.type) +
  geom_histogram(binwidth = 1,alpha=0.3,colour="black")+
  scale_fill_manual(labels = c("Fished", "<10km2","10-100km2",">100km2"),
                    values=c("red","darkmagenta","orange","darkgreen"))+
  theme_classic()+
  xlab('Gravity')+
  ylab('')+
  ggtitle('(d) NSW')+
  scale_x_continuous(expand = c(0,0),limits=c(0,14))+
  scale_y_continuous(expand=c(0,0),labels = scales::percent)+
  Theme1+
  theme(legend.direction = "horizontal",
        legend.position = "bottom")
NSW.plot

SA.plot<-ggplot(data=ntr%>%filter(state=='SA')) + aes(x = gravity.50, y = ..count../sum(..count..),fill=reserve.type) +
  geom_histogram(binwidth = 1,alpha=0.3,colour="black")+
  scale_fill_manual(labels = c("Fished", "<10km2","10-100km2",">100km2"),
                    values=c("red","darkmagenta","orange","darkgreen"))+
  theme_classic()+
  xlab('Gravity')+
  ylab('')+
  ggtitle('(e) SA')+
  scale_x_continuous(expand = c(0,0),limits=c(0,14))+
  scale_y_continuous(expand=c(0,0),labels = scales::percent)+
  Theme1+
  theme(legend.direction = "horizontal",
        legend.position = "bottom")
SA.plot

TAS.plot<-ggplot(data=ntr%>%filter(state=='TAS')) + aes(x = gravity.50, y = ..count../sum(..count..),fill=reserve.type) +
  geom_histogram(binwidth = 1,alpha=0.3,colour="black")+
  scale_fill_manual(labels = c("Fished", "<10km2","10-100km2",">100km2"),
                    values=c("red","darkmagenta","orange","darkgreen"))+
  theme_classic()+
  xlab('Gravity')+
  ylab('')+
  ggtitle('(f) TAS')+
  scale_x_continuous(expand = c(0,0),limits=c(0,14))+
  scale_y_continuous(expand=c(0,0),labels = scales::percent)+
  Theme1+
  theme(legend.direction = "horizontal",
        legend.position = "bottom")
TAS.plot

### Arrage plots

ggarrange(WA.plot,VIC.plot,NSW.plot,SA.plot,TAS.plot,
          ncol=5,nrow = 1,align = 'hv',
          common.legend = T,legend = "bottom")

### Export plot

setwd(plots.dir)
dir()

name<-'Reserve_size'

ggsave(paste(name,"gravity.50.png",sep="."),width = 21, height = 6,units = "cm")
ggsave(paste(name,"gravity.50.pdf",sep="."),width = 21, height = 6,units = "cm",useDingbats=FALSE)

### Arrange with cowplot

library(cowplot)

bottom_row=ggarrange(WA.plot,VIC.plot,NSW.plot,SA.plot,TAS.plot,
                     ncol=5,nrow = 1,common.legend=T,legend="bottom",align='hv')
bottom_row

ggdraw()+
  draw_plot(gravity.plot,x=0.1,y=0.4,width = 0.8,height=0.6)+
  draw_plot(bottom_row,x=0,y=0,width = 1 ,height=0.4)

### Export plot

setwd(plots.dir)
dir()

name<-'Legal'

ggsave(paste(name,"Gravity.50_NTR.png",sep="."),width = 21, height = 20,units = "cm")
ggsave(paste(name,"Gravity.50_NTR.pdf",sep="."),width = 21, height = 20,units = "cm",useDingbats=FALSE)

#### Geographic extent of occurrence species groups ----

## Bring in Australia coastline

setwd('C:/Users/22373243/Dropbox/Nestor_ascii_2')
dir()

australia<-readOGR('.','Australiaboundary67')
proj4string(australia)
australia <- spTransform(australia, CRS("+proj=longlat +datum=WGS84"))

## Crop australia to the extent

e<-extent(112,155,-45,-7)

australia<-crop(australia,e)
plot(australia)

australia.continental<-fortify(australia)

## Bring in individual species groups dataset and trimmed to their biogeographic extent of occurrence (excl. vagrant records)

## (1) Pink Snapper

setwd(data.dir)
dir()

snapper.dat<-read.csv("length_analysis_Chrysophurus_auratus.csv")%>%
  dplyr::group_by(state,id,longitude,latitude,scientific)%>%
  dplyr::summarise(response=sum(response))%>%
  na.omit()%>%
  glimpse()

## Create presence/absence dataset

snapper.incidence<-snapper.dat%>%
  mutate(response=ifelse(response>0,1,response))%>%
  glimpse()

## Create a dataframe with presence only data

snapper.presence<-snapper.incidence%>%
  filter(response>0)%>%
  glimpse()

## Check distributional of presence/absence data

snapper.map<-ggplot()+
  geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="gray12",alpha=0.8)+
  geom_point(data=snapper.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=snapper.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,160))+
  scale_y_continuous(expand = c(0,0),limits = c(-45,-9))+
  ggtitle('(a) Chrysophrys auratus')+
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(size=10))
snapper.map

## (2) Emperors

setwd(data.dir)
dir()

emperors.dat<-read.csv("length_analysis_Lethrinidae.csv")%>%
  dplyr::group_by(state,id,longitude,latitude,scientific)%>%
  dplyr::summarise(response=sum(response))%>%
  dplyr::filter(!grepl('Scott & Rowleys',id))%>%
  na.omit()%>%
  glimpse()

## Create presence/absence dataset

emperors.incidence<-emperors.dat%>%
  mutate(response=ifelse(response>0,1,response))%>%
  glimpse()

## Create a dataframe with presence only data

emperors.presence<-emperors.incidence%>%
  filter(response>0)%>%
  glimpse()

## Check distributional of presence/absence data

emperors.map<-ggplot()+
  geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="gray12",alpha=0.8)+
  geom_point(data=emperors.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=emperors.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,160))+
  scale_y_continuous(expand = c(0,0),limits = c(-45,-9))+
  ggtitle('(b) Lethrinus spp.')+
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(size=10))
emperors.map

## (3) Coral Trouts

setwd(data.dir)
dir()

coral.trouts.dat<-read.csv("length_analysis_Plectropomus_spp.csv")%>%
  dplyr::group_by(state,id,longitude,latitude,scientific)%>%
  dplyr::summarise(response=sum(response))%>%
  dplyr::filter(!grepl('Scott & Rowleys',id))%>%
  na.omit()%>%
  glimpse()

## Create presence/absence dataset

coral.trouts.incidence<-coral.trouts.dat%>%
  mutate(response=ifelse(response>0,1,response))%>%
  glimpse()

## Create a dataframe with presence only data

coral.trouts.presence<-coral.trouts.incidence%>%
  filter(response>0)%>%
  glimpse()

## Check distributional of presence/absence data

coral.trouts.map<-ggplot()+
  geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="gray12",alpha=0.8)+
  geom_point(data=coral.trouts.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=coral.trouts.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,160))+
  scale_y_continuous(expand = c(0,0),limits = c(-45,-9))+
  ggtitle('(c) Plectropomus spp.')+
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(size=10))
coral.trouts.map

## (4) Baldchin groper

setwd(data.dir)
dir()

baldchin.dat<-read.csv("length_analysis_Choerodon_spp.csv")%>%
  dplyr::group_by(state,id,longitude,latitude,scientific)%>%
  dplyr::summarise(response=sum(response))%>%
  dplyr::filter(!grepl('Scott & Rowleys',id))%>%
  na.omit()%>%
  glimpse()

## Create presence/absence dataset

baldchin.incidence<-baldchin.dat%>%
  mutate(response=ifelse(response>0,1,response))%>%
  glimpse()

## Create a dataframe with presence only data

baldchin.presence<-baldchin.incidence%>%
  filter(response>0)%>%
  glimpse()

## Check distributional of presence/absence data

baldchin.map<-ggplot()+
  geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="gray12",alpha=0.8)+
  geom_point(data=baldchin.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=baldchin.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,160))+
  scale_y_continuous(expand = c(0,0),limits = c(-45,-9))+
  ggtitle('(d) Choerodon spp.')+
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(size=10))
baldchin.map

## (5) Morwongs

setwd(data.dir)
dir()

morwongs.dat<-read.csv("length_analysis_Nemadactylus_spp.csv")%>%
  dplyr::group_by(state,id,longitude,latitude,scientific)%>%
  dplyr::summarise(response=sum(response))%>%
  dplyr::filter(!grepl('Scott & Rowleys',id))%>%
  na.omit()%>%
  glimpse()

## Create presence/absence dataset

morwongs.incidence<-morwongs.dat%>%
  mutate(response=ifelse(response>0,1,response))%>%
  glimpse()

## Create a dataframe with presence only data

morwongs.presence<-morwongs.incidence%>%
  filter(response>0)%>%
  glimpse()

## Check distributional of presence/absence data

morwongs.map<-ggplot()+
  geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="gray12",alpha=0.8)+
  geom_point(data=morwongs.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=morwongs.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,160))+
  scale_y_continuous(expand = c(0,0),limits = c(-45,-9))+
  ggtitle('(e) Nemadactylus spp.')+
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(size=10))
morwongs.map

## (6) Wrasses

setwd(data.dir)
dir()

wrasses.dat<-read.csv("length_analysis_Notolabrus_spp.csv")%>%
  dplyr::group_by(state,id,longitude,latitude,scientific)%>%
  dplyr::summarise(response=sum(response))%>%
  dplyr::filter(!grepl('Scott & Rowleys',id))%>%
  na.omit()%>%
  glimpse()

## Create presence/absence dataset

wrasses.incidence<-wrasses.dat%>%
  mutate(response=ifelse(response>0,1,response))%>%
  glimpse()

## Create a dataframe with presence only data

wrasses.presence<-wrasses.incidence%>%
  filter(response>0)%>%
  glimpse()

## Check distributional of presence/absence data

wrasses.map<-ggplot()+
  geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="gray12",alpha=0.8)+
  geom_point(data=wrasses.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=wrasses.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,160))+
  scale_y_continuous(expand = c(0,0),limits = c(-45,-9))+
  ggtitle('(f) Notolabrus spp.')+
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(size=10))
wrasses.map

## (7) Lutjanids

setwd(data.dir)
dir()

Lutjanids.dat<-read.csv("length_analysis_Lutjanidae.csv")%>%
  dplyr::group_by(state,id,longitude,latitude,scientific)%>%
  dplyr::summarise(response=sum(response))%>%
  dplyr::filter(!grepl('Scott & Rowleys',id))%>%
  na.omit()%>%
  glimpse()

## Create presence/absence dataset

Lutjanids.incidence<-Lutjanids.dat%>%
  mutate(response=ifelse(response>0,1,response))%>%
  glimpse()

## Create a dataframe with presence only data

Lutjanids.presence<-Lutjanids.dat%>%
  filter(response>0)%>%
  glimpse()


## Check distributional of presence/absence data

Lutjanids.map<-ggplot()+
  geom_polygon(data=australia,aes(x=long,y=lat,group=group),fill="gray12",alpha=0.8)+
  geom_point(data=Lutjanids.incidence,aes(x=longitude,y=latitude),colour="red",size=3)+
  geom_point(data=Lutjanids.presence,aes(x=longitude,y=latitude),colour="blue",size=3)+
  scale_x_continuous(expand=c(0,0),limits = c(110,160))+
  scale_y_continuous(expand = c(0,0),limits = c(-45,-9))+
  ggtitle('(g) Lutjanus spp.')+
  theme_classic()+
  Theme1+
  theme(plot.title = element_text(size=10))
Lutjanids.map

### Arrage plots

ggarrange(snapper.map,emperors.map,coral.trouts.map,
          baldchin.map,morwongs.map,wrasses.map,Lutjanids.map,
          ncol=3,nrow = 3,align = 'hv',
          common.legend = T,legend = "bottom")

### Export plot

setwd(plots.dir)
dir()

name<-'EDF2_Species_groups_map'

ggsave(paste(name,".png",sep="."),width = 21, height = 20,units = "cm")
ggsave(paste(name,".pdf",sep="."),width = 21, height = 20,units = "cm",useDingbats=FALSE)


