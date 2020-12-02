##########################################################################
#
# MACGILLIVRAY PRION ADULT SURVIVAL FOR GOUGH ISLAND
#
##########################################################################

# written by Steffen.oppel@rspb.org.uk in December 2020
# adapted from Montserrat Oriole code from 2014

library(tidyverse)
library(data.table)
library(lubridate)
library(RMark)
filter<-dplyr::filter
select<-dplyr::select


#########################################################################
# LOAD AND PREPARE DATA FOR ADULT ANNUAL SURVIVAL ESTIMATION
#########################################################################
## completely revised on 19 Oct 2020

##### LOAD FORMATTED RINGING DATA ###########
setwd("C:/STEFFEN/RSPB/UKOT/Gough/ANALYSIS/SeabirdSurvival")

## run the RODBC import of CMR data in a 32-bit version of R
#system(paste0("C:/PROGRA~1/R/R-35~1.1/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival\\RODBC_CMR_import.R")), wait = TRUE, invisible = FALSE, intern = T)
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival"), silent=T)
load("GOUGH_seabird_CMR_data.RData")

###### FILTER ADULT DATA FROM RAW CONTACTS ########
contacts<-contacts %>% filter(SpeciesCode %in% c("MAPR","BBPR","PRIO")) %>% filter(Location=="Prion Cave") %>%
  mutate(Age=ifelse(is.na(Age),"Adult",as.character(Age))) %>%
  mutate(Contact_Season=ifelse(is.na(Contact_Season),"2017-18",as.character(Contact_Season))) %>%
  mutate(Contact_Season=ifelse(Contact_Season=="2020-21","2019-20",as.character(Contact_Season))) %>%
  filter(!Age=="Chick")

head(contacts)  
unique(contacts$Age)

EncHist<-contacts %>% group_by(BirdID,Contact_Season) %>%
  summarise(n=length(Date_Time)) %>%
  spread(key=Contact_Season,value=n,fill=0)  ### 0 for 'not seen'
dim(EncHist)

#### FORMAT FOR SIMPLE CJS MODEL ############
CH<-as.matrix(EncHist[,2:ncol(EncHist)], dimnames=F)
CH<-ifelse(CH>0,1,2)  ##1=seen, 2=not seen



########################################################################
############### SIMPLE ANALYSIS WITH TRANSIENTS ##################
########################################################################

MARKcapt.hist<-contacts %>% group_by(BirdID,Contact_Season) %>%
  summarise(n=length(Date_Time)) %>%
  mutate(n=ifelse(n>0,1,0)) %>%
  spread(key=Contact_Season,value=n,fill=0) %>%
  unite("ch", 2:tail(names(.),1), sep = "")


#creating the MARK project in R - defining what model and what groups, and the year when first individual was marked
MAPR.proc<-process.data(MARKcapt.hist, model="CJS", begin.time=2014) 

#specifying the design data (PIMs)
MAPR.ddl<-make.design.data(MAPR.proc)

# adding the transients
MAPR.ddl=add.design.data(MAPR.proc,MAPR.ddl,parameter="Phi",type="age",bins=c(0,1,6),name="transients", replace=TRUE)

### Automated workflow to run many models in a function wrapper ###

# Create a function called do.analysis that will create all the models that I want
do.analysis<-function()
{
# Within the function define the set of formulae that I want for Phi and p
Phi.dot<-list(formula=~1)
Phi.trans<-list(formula=~transients)
Phi.time<-list(formula=~time)
p.time<-list(formula=~time)
p.dot<-list(formula=~1)

cml<-create.model.list("CJS")

model.list<-mark.wrapper(cml,data=MAPR.proc,ddl=MAPR.ddl)

# Return the list of model results as the value of the function
return(model.list)
}
myresults<-do.analysis()

### show the results and calculate Evidence ratios and AIC weights for model table
myresults
model.table<-model.table(myresults,use.lnl=TRUE)
model.table$EvidenceRatio<-2.71828182845904523536^(-0.5*model.table$DeltaAICc) 
model.table
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

write.table(model.table, "MAPR_surv_model_selection_table_RMark.csv", sep=",")




########################################################################
############### CLEAN UP WORKSPACE ##################
########################################################################

rm(list=ls(all=TRUE))
cleanup(ask=FALSE)

