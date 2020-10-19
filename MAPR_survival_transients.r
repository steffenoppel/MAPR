##########################################################################
#
# MACGILLIVRAY PRION SURVIVAL ANALYSIS 2014-2018
#
##########################################################################
# based on Kery and Schaub 2012, Chapter 7
# modified by Steffen oppel, October 2018

## UPDATED ON 18 JUNE 2020 to account for transients
## leads to completely unrealistic estimate of survival
#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0222241

## MODIFIED ON 19 OCTOBER 2020 to change input data
## adjusted temporal variability of p and phi
## model converged and yielded sensible estimates -> adopt for population model!


library(tidyverse)
library(data.table)
library(lubridate)
library(R2WinBUGS)
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
## find the birds marked very late in 2015:
latemarked<-contacts %>% filter(SpeciesCode %in% c("MAPR","BBPR","PRIO")) %>% filter(Location=="Prion Cave") %>% filter(!Age=="Chick") %>%
  filter(month(Date_Time) %in% c(12,1,2)) %>% group_by(BirdID) %>% summarise(n=length(Date_Time)) ### %>% filter(day(Date_Time)>15) 

transients<-contacts %>% filter(SpeciesCode %in% c("MAPR","BBPR","PRIO")) %>% filter(Location=="Prion Cave") %>% filter(!Age=="Chick") %>%
  group_by(BirdID) %>% summarise(n=length(Date_Time)) %>% filter(n==1) 


contacts<-contacts %>% filter(SpeciesCode %in% c("MAPR","BBPR","PRIO")) %>% filter(Location=="Prion Cave") %>%
  #filter(is.na(Age))
  mutate(Age=ifelse(is.na(Age),"Adult",as.character(Age))) %>%
  mutate(Contact_Season=ifelse(is.na(Contact_Season),"2017-18",as.character(Contact_Season))) %>%
  mutate(Contact_Season=ifelse(Contact_Season=="2020-21","2019-20",as.character(Contact_Season))) %>%
  filter(!Age=="Chick")
#filter(Contact_Year==2015) %>% filter(month(Date_Time)==12) %>% filter(day(Date_Time)>15)
head(contacts)  ## seabird count data
unique(contacts$Age)

EncHist<-contacts %>% group_by(BirdID,Contact_Season) %>%
  summarise(n=length(Date_Time)) %>%
  spread(key=Contact_Season,value=n,fill=0)  ### 0 for 'not seen'
dim(EncHist)

#### FORMAT FOR SIMPLE CJS MODEL ############
CH<-as.matrix(EncHist[,2:ncol(EncHist)], dimnames=F)
CH<-ifelse(CH>0,1,2)  ##1=seen, 2=not seen



#########################################################################
# SPECIFY MULTI-EVENT MODEL ACCOUNTING FOR TRANSIENTS
#########################################################################

sink("MAPR_Transience_MultiEvent_ptime.jags")
cat("
    model {
    
    # -------------------------------------------------
    # Parameters:
    # phi: survival probability for adults
    # p: recapture probability when breeding
    # -------------------------------------------------

    # Priors and constraints
    for (t in 1:(n.occasions-1)){
      #logit(phi[t]) <- mu.phi + surv.raneff[t]
      #surv.raneff[t] ~ dnorm(0, tau.phi)
      p[t] ~ dunif(0, 1)
      emigrate[t] ~ dunif(0,0.25)
    }

    mean.phi ~ dunif(0, 1)             # Prior for mean survival
    #mu.phi <- log(mean.phi / (1-mean.phi)) # Logit transformation
    #sigma.phi ~ dunif(0, 5)                # Prior for standard deviation
    #tau.phi <- pow(sigma.phi, -2)
    #emigrate ~ dunif(0,0.25)

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
    ps[2,i,t,2]<-mean.phi*(1-emigrate[t])
    ps[2,i,t,3]<-mean.phi*emigrate[t]
    
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
parameters <- c("mean.phi","phi", "p","emigrate")

# MCMC settings
ni <- 10000
nt <- 2
nb <- 5000
nc <- 3

# Call JAGS from R
MAPRsurv <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival\\MAPR_Transience_MultiEvent_ptime.jags",
                 n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb,parallel=T) #  




#########################################################################
# PRODUCE OUTPUT TABLE
#########################################################################

out<-as.data.frame(MAPRsurv$summary)
out$parameter<-row.names(MAPRsurv$summary)
out
write.table(out,"MAPR_Gough_Survival_estimates_transients.csv", sep=",", row.names=F)



#########################################################################
# PRODUCE OUTPUT GRAPH
#########################################################################
## only makes sense if survival varies by year

# out[1:11,] %>% select(c(1,5,2,3,7)) %>%
#   setNames(c('Mean', 'Median','SD','lcl', 'ucl')) %>%
#   mutate(Year=colnames(rCH)[1:11]) %>%
#   
#   ggplot(aes(y=Median, x=Year)) + geom_point(size=2.5)+
#   geom_errorbar(aes(ymin=lcl, ymax=ucl), width=.1)+
#   ylab("Annual adult survival probability") +
#   scale_y_continuous(breaks=seq(0.5,1,0.1), limits=c(0.5,1))+
#   #scale_x_continuous(breaks=seq(2006,2017,1))+
#   theme(panel.background=element_rect(fill="white", colour="black"), 
#         axis.text=element_text(size=18, color="black"), 
#         axis.title=element_text(size=20),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(), 
#         panel.border = element_blank())
# 
# dev.off()



