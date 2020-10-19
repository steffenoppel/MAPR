##########################################################################
#
# MACGILLIVRAY PRION POPULATION MODEL FOR GOUGH ISLAND
#
##########################################################################
# based on Kery and Schaub 2012, Chapter 7
# modified by Steffen oppel, October 2018

library(tidyverse)
library(jagsUI)
library(data.table)
library(lubridate)
library(popbio)
library(doParallel)
library(foreach)



#########################################################################
# 1. LOAD AND PREPARE DATA FOR ANNUAL SURVIVAL ESTIMATION
#########################################################################


##### LOAD FORMATTED RINGING DATA ###########
setwd("A:/RSPB/UKOT/Gough/ANALYSIS/SeabirdSurvival")
setwd("C:/STEFFEN/RSPB/UKOT/Gough/ANALYSIS/SeabirdSurvival")
#raw<- as.data.frame(fread("MAPR_enc_hist.csv", header=T))
## run the RODBC import of CMR data in a 32-bit version of R
system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival\\RODBC_CMR_import.R")), wait = TRUE, invisible = FALSE)
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival"), silent=T)
load("GOUGH_seabird_CMR_data.RData")

###### FILTER DATA FROM RAW CONTACTS ########
contacts<-contacts %>% filter(SpeciesCode=="MAPR")
head(contacts)  ## seabird count data

EncHist<-contacts %>% group_by(BirdID,Contact_Year) %>%
  summarise(n=length(Date_Time)) %>%
  spread(key=Contact_Year,value=n,fill=0)

#### FORMAT FOR SIMPLE CJS MODEL ############
CH<-as.matrix(EncHist[,2:ncol(EncHist)], dimnames=F)
CH<-ifelse(CH>0,1,0)
head(CH)

# ELIMINATE TRANSIENTS ONLY OBSERVED IN A SINGLE YEAR
del <- apply(CH[,1:ncol(CH)], 1, sum)
dim(CH)
rCH<-CH[!(del==1),]
dim(rCH)

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x==1))
f <- apply(rCH, 1, get.first)

# Bundle data
jags.data <- list(y = rCH, f = f, n.occasions = dim(rCH)[2], nind = dim(rCH)[1])

# Initial values 
inits <- function(){list(phi = runif(1, 0.7, 1),
                         pp = runif(1, 0, 1))}
 

# Parameters monitored
parameters <- c("phi", "p")

# MCMC settings
ni <- 10000
nt <- 2
nb <- 5000
nc <- 4

# Call JAGS from R (model created in C:\STEFFEN\RSPB\UKOT\Gough\ANALYSIS\SeabirdSurvival\MAPR_survival.r)
MAPRsurv <- autojags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival\\MAPR_CJS.jags",
                 n.chains = nc, n.thin = nt, n.burnin = nb,parallel=T) # n.iter = ni, 

## compile output
out<-as.data.frame(MAPRsurv$summary)
out$parameter<-row.names(MAPRsurv$summary)
out



#########################################################################
# 2. LOAD AND PREPARE DATA FOR BREEDING SUCCESS SUMMARY
#########################################################################

## run the RODBC import in a 32-bit version of R
system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdBreedingSuccess\\RODBC_nest_import.r")), wait = TRUE, invisible = FALSE)
#system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\Users\\Gough Conservation\\Documents\\Gough Birders\\2018-2019\\12.Monthly reports 2018-19\\RODBC_imports.r")), wait = FALSE, invisible = FALSE)

try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdBreedingSuccess"), silent=T)
#try(setwd("C:\\Users\\Gough Conservation\\Documents\\Gough Birders\\2018-2019\\12.Monthly reports 2018-19"), silent=T)

load("GOUGH_nest_data.RData")
head(nestsDB)  ## nest data
head(visDB)  ## nest visit data


##  SELECT DATA FOR TARGET SPECIES AND SUMMARISE NEST SUCCESS ####

head(nestsDB)
export<-nestsDB %>% filter(Species==SP) %>% filter(Year>2013) %>%
  mutate(count=1) %>%
  group_by(Species,Year) %>%
  summarise(SampSize=sum(count),BreedSucc=mean(SUCCESS))





#########################################################################
# 3. SPECIFY SIMPLE POPULATION MODEL AND ESTIMATE POPULATION GROWTH RATE
#########################################################################

### SPECIFY RANGE OF PARAMETERS ###

#pop.size<-seq(5000,10000,1000)				### population size in individuals
Sb<-seq(out[1,3],out[1,7],0.01)					### survival of adult breeders from model above
S1<-seq(0.4,0.9,0.1)		            ### survival of first year birds from ...
F<-max(export$BreedSucc)						### breeding success, take maximum as this is anyway very low


### CREATING THE POPULATION MATRIX ###

seabird.matrix<-expression(
  0,0,(F*0.5*S1),
  S1,0,0,
  0,Sb,Sb)


stable.stage(seabird.matrix)


### COMPREHENSIVE TABLE OF ALL COMBINATIONS OF DEMOGRAPHIC PARAMETERS ###

simul_in<-expand.grid(Sb, S1,F)
dim(simul_in)
names(simul_in)<-c('Sb','S1','F')
SIM_OUT<-data.frame()


### setup parallel backend to use 8 processors

cl<-makeCluster(8)
registerDoParallel(cl, cores=8)



### CALCULATING STABLE AGE DISTRIBUTION

SIM_OUT<-foreach(s=c(1:dim(simul_in)[1]), .packages='popbio',.combine=rbind) %dopar% {
  
  
### CREATE LESLIE MATRIX WITH SUBSET OF VITAL RATES
  
  seabird.vr<-list(F=simul_in[s,3],S1=simul_in[s,2],Sb=simul_in[s,1])
  A<-matrix(sapply(seabird.matrix, eval,seabird.vr , NULL), nrow=sqrt(length(seabird.matrix)), byrow=TRUE)
  
  
  ### CALCULATING POPULATION GROWTH RATE FOR COMBINATION OF DEMOGRAPHIC PARAMETERS
  
  out<-simul_in[s,]
  out$lambda<-lambda(A)
  return(out)
}
stopCluster(cl)






#########################################################################
# PRODUCE POPULATION PROJECTIONS
#########################################################################

cl<-makeCluster(8)
registerDoParallel(cl, cores=8)

POP_Project<-foreach(s=c(1:dim(simul_in)[1]), .packages='tidyverse',.combine=rbind) %dopar% {
  
  traj<-expand.grid(year=seq(2020,2050), pop.size=10000,Sb=SIM_OUT[s,1]) %>% left_join(SIM_OUT[s,], by="Sb")
    for (y in seq(2021,2050)) {
      traj$pop.size[traj$year==y]<-traj$pop.size[traj$year==(y-1)]*SIM_OUT[s,4]
    }
  traj$Scenario<-s
  return(traj)
  
}
stopCluster(cl)

head(POP_Project)  
  
  
#########################################################################
# PRODUCE OUTPUT GRAPH OF POPULATION SIMULATIONS
#########################################################################
POP_Project %>% mutate(Scenario=as.factor(Scenario)) %>%
  
ggplot(aes(y=pop.size, x=year, colour=Scenario)) + geom_line(size=0.8)+
  ylab("Population size") +

  theme(panel.background=element_rect(fill="white", colour="black"),
        axis.text=element_text(size=18, color="black"),
        axis.title=element_text(size=20),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank())




