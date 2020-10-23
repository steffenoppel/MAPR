##########################################################################
#
# MACGILLIVRAY PRION SURVIVAL ANALYSIS 2014-2018
#
##########################################################################
# based on Kery and Schaub 2012, Chapter 7
# modified by Steffen oppel, October 2018

## UPDATED ON 18 JUNE 2020 to account for transients
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0222241

## MODIFIED ON 19 OCTOBER 2020 to change input data
## adjusted temporal variability of p and phi
## model converged and yielded sensible estimates -> adopt for population model!

library(tidyverse)
library(data.table)
library(lubridate)
library(jagsUI)
filter<-dplyr::filter
select<-dplyr::select


##### LOAD FORMATTED RINGING DATA ###########
#setwd("A:/RSPB/UKOT/Gough/ANALYSIS/SeabirdSurvival")
setwd("C:/STEFFEN/RSPB/UKOT/Gough/ANALYSIS/SeabirdSurvival")
#raw<- as.data.frame(fread("MAPR_enc_hist.csv", header=T))
## run the RODBC import of CMR data in a 32-bit version of R
#system(paste0("C:/PROGRA~1/R/R-35~1.1/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival\\RODBC_CMR_import.R")), wait = TRUE, invisible = FALSE)
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

#### FORMAT FOR MULTIEVENT CJS MODEL ############
CH<-as.matrix(EncHist[,2:ncol(EncHist)], dimnames=F)
CH<-ifelse(CH>0,1,2)  ##1=seen, 2=not seen



#########################################################################
# SPECIFY MULTI-EVENT MODEL ACCOUNTING FOR TRANSIENTS
#########################################################################

sink("MAPR_Transience_MultiEvent_ptime_constemig.jags")
cat("
    model {
    
    # -------------------------------------------------
    # Parameters:
    # phi: survival probability for adults
    # p: recapture probability when breeding
    # -------------------------------------------------

    # Priors and constraints
    for (t in 1:(n.occasions-1)){
      p[t] ~ dunif(0, 1)
    }

    mean.phi ~ dunif(0, 1)             # Prior for mean survival
    emigrate ~ dunif(0,0.25)

    # States (S):
    # 1 dead
    # 2 alive in Prion Cave
    # 3 alive as transient
    
    # Observations (O):
    # 1 observed
    # 2 not observed
    # -------------------------------------------------
    

    
    # -------------------------------------------------
    # Define state-transition and observation matrices 
    # -------------------------------------------------
    
    for (i in 1:nind){
    
    for (t in f[i]:(n.occasions-1)){
    
    # Define probabilities of state S(t+1) [last dim] given S(t) [first dim]
    
    ps[1,i,t,1]<-1    ## dead birds stay dead
    ps[1,i,t,2]<-0
    ps[1,i,t,3]<-0
    
    ps[2,i,t,1]<-(1-mean.phi)
    ps[2,i,t,2]<-mean.phi*(1-emigrate)
    ps[2,i,t,3]<-mean.phi*emigrate
    
    ps[3,i,t,1]<-(1-mean.phi)
    ps[3,i,t,2]<-0
    ps[3,i,t,3]<-mean.phi
    
    # Define probabilities of O(t) [last dim] given S(t)  [first dim]
    
    po[1,i,t,1]<-0
    po[1,i,t,2]<-1
    
    po[2,i,t,1]<-p[t]
    po[2,i,t,2]<-(1-p[t])
    
    po[3,i,t,1]<-0
    po[3,i,t,2]<-1
    
    } #t
    } #i


    # Likelihood 
    for (i in 1:nind){
      # Define latent state at first capture
      z[i,f[i]] <- 2 ## alive when first marked
      for (t in (f[i]+1):n.occasions){
        # State process: draw S(t) given S(t-1)
        z[i,t] ~ dcat(ps[z[i,t-1], i, t-1,])
        # Observation process: draw O(t) given S(t)
        y[i,t] ~ dcat(po[z[i,t], i, t-1,])
      } #t
    } #i


    }
    ",fill = TRUE)
sink()



#########################################################################
# PREPARE DATA FOR MODEL
#########################################################################
head(CH)


# Compute vector with occasion of first capture
get.first <- function(x) min(which(x==1))
f <- apply(CH, 1, get.first)

# Bundle data
jags.data <- list(y = CH, f = f, n.occasions = dim(CH)[2], nind = dim(CH)[1])

# Initial values 
inits <- function(){list(mean.phi = runif(1, 0.7, 1),
                         p = runif(dim(CH)[2]-1, 0, 1))}
 

# Parameters monitored
parameters <- c("mean.phi","phi", "mean.p","p","emigrate")

# MCMC settings
ni <- 30000
nt <- 5
nb <- 5000
nc <- 3



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########  FIT AND COMPARE MODELS WITH DIFFERENT TIME PARAMETERIATION    ###########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
inits <- function(){list(mean.phi = runif(1, 0.7, 1),
                         p = runif(dim(CH)[2]-1, 0, 1))}
MAPRall <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\MAPR\\MAPR_Transience_MultiEvent.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T) #  
out1<-as.data.frame(MAPRall$summary)
out1$parameter<-row.names(MAPRall$summary)
out1$model<-"fulltime"


inits <- function(){list(mean.phi = runif(1, 0.7, 1),
                         mean.p = runif(1, 0, 1))}
MAPRphitime <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\MAPR\\\\MAPR_Transience_MultiEvent_phitime.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T) # 
out2<-as.data.frame(MAPRphitime$summary)
out2$parameter<-row.names(MAPRphitime$summary)
out2$model<-"phitime"


inits <- function(){list(mean.phi = runif(1, 0.7, 1),
                         p = runif(dim(CH)[2]-1, 0, 1))}
MAPRptime <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\MAPR\\MAPR_Transience_MultiEvent_ptime.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T) # 
out3<-as.data.frame(MAPRptime$summary)
out3$parameter<-row.names(MAPRptime$summary)
out3$model<-"ptime"



inits <- function(){list(mean.phi = runif(1, 0.7, 1),
                         mean.p = runif(1, 0, 1))}
MAPRconstant <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\MAPR\\MAPR_Transience_MultiEvent_constant.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T) # 
out4<-as.data.frame(MAPRconstant$summary)
out4$parameter<-row.names(MAPRconstant$summary)
out4$model<-"constant"



inits <- function(){list(mean.phi = runif(1, 0.7, 1),
                         p = runif(dim(CH)[2]-1, 0, 1))}
MAPRptimeconstemig <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\MAPR\\MAPR_Transience_MultiEvent_ptime_constemig.jags",
                  n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T) # 
out5<-as.data.frame(MAPRptimeconstemig$summary)
out5$parameter<-row.names(MAPRptimeconstemig$summary)
out5$model<-"ptime_constemig"



MAPRconstemig <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\MAPR\\MAPR_Transience_MultiEvent_constemig.jags",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T) #  
out6<-as.data.frame(MAPRconstemig$summary)
out6$parameter<-row.names(MAPRconstemig$summary)
out6$model<-"constemig"


### COMBINE OUTPUT FROM ALL 6 MODELS
out<-bind_rows(out1,out2,out3,out4,out5,out6)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# EXTRACT DIC TO COMPARE MODELS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### CALCULATE AICc:
## AICc is the Akaike Information Criterion corrected for small sample sizes calculated as:
## (2 * K - deviance) + (2 * K) * (K+1) / (n - K - 1)
## DIC is essential for highly complex hierarchical models - https://rss.onlinelibrary.wiley.com/doi/full/10.1111/1467-9868.00353
## infos about DIC: https://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-dic/#q12

pd_dic <- function(x) {
  data.frame(n.eff.parameters=x$pD, n.parameters=dim(x$summary)[1]-3,DIC=x$DIC)
}
DIC_tab<-bind_rows(pd_dic(MAPRall),pd_dic(MAPRphitime),pd_dic(MAPRptime),
                   pd_dic(MAPRconstant),pd_dic(MAPRconstemig),pd_dic(MAPRptimeconstemig)) %>%
  mutate(model=c("fulltime","phitime","ptime","constant", "constemig","ptime_constemig")) %>%
  arrange(DIC) %>%
  mutate(deltaDIC=DIC-DIC[1])

ModSelTab<-out %>% dplyr::select(model, parameter,mean) %>%
  filter(parameter=="deviance") %>%
  mutate(mean=round(mean,3)) %>%
  spread(key=parameter, value=mean, fill="not included") %>%
  left_join(DIC_tab, by="model") %>%
  mutate(n.parameters=ifelse(model=="fulltime",n.parameters-1,n.parameters)) %>% ### remove the derived 'mean.phi' from parameter count
  mutate(n.parameters=ifelse(model=="ran_time",n.parameters-4,n.parameters)) %>% ### remove the 'phi[t]' from parameter count
  mutate(AICc=(2*n.parameters)+as.numeric(deviance) + (2*n.parameters) * (n.parameters+1) / (dim(rCH)[1]-n.parameters-1)) %>%
  arrange(AICc) %>%
  mutate(deltaAIC=AICc-AICc[1])
ModSelTab

fwrite(ModSelTab,"TableS1_ModelSelection.csv")

