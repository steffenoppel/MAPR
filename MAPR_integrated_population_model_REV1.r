##########################################################################
#
# MACGILLIVRAY PRION POPULATION MODEL FOR GOUGH ISLAND
#
##########################################################################

# population model adapted from Jiguet, F., Robert, A., Micol, T. and Barbraud, C. (2007), Quantifying stochastic and deterministic threats to island seabirds: last endemic prions face extinction from falcon peregrinations. Animal Conservation, 10: 245-253. doi:10.1111/j.1469-1795.2007.00100.x
# implemented in JAGS based on Kery and Schaub 2012
# written by Steffen.oppel@rspb.org.uk in May 2020

# REVISION in OCTOBER 2020:
# updated CMR data input - removed chicks and changed encounter occasion (from year to season)
# switched to m-array to allow GoF test for survival model - but used this only for GoF test
# included temporal variation in phi and p in survival model
# included all transients and switched to a multi-event survival model
# re-inserted fecundity from previous version v3 that distinguishes between good and bad years

library(tidyverse)
library(jagsUI)
library(data.table)
library(lubridate)
library(popbio)
library(janitor)
filter<-dplyr::filter
select<-dplyr::select


#########################################################################
# 1. LOAD AND PREPARE DATA FOR ADULT ANNUAL SURVIVAL ESTIMATION
#########################################################################
## completely revised on 19 Oct 2020

##### LOAD FORMATTED RINGING DATA ###########
setwd("C:/STEFFEN/RSPB/UKOT/Gough/ANALYSIS/SeabirdSurvival")

## run the RODBC import of CMR data in a 32-bit version of R
system(paste0("C:/PROGRA~1/R/R-35~1.1/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival\\RODBC_CMR_import.R")), wait = TRUE, invisible = FALSE, intern = T)
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

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x==1))
f <- apply(CH, 1, get.first)


#### SUMMARISE NUMBER OF DETECTIONS
EncHist %>% ungroup() %>% gather(key='Season', value='capt',-BirdID) %>% group_by(BirdID) %>%
  summarise(n_capt=sum(capt)) %>% tabyl(n_capt)
  


#########################################################################
# 2. LOAD AND PREPARE DATA FOR BREEDING SUCCESS SUMMARY
#########################################################################

## run the RODBC import in a 32-bit version of R
system(paste0("C:/PROGRA~1/R/R-35~1.1/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdBreedingSuccess\\RODBC_nest_import.r")), wait = TRUE, invisible = FALSE, intern=T)
#system(paste0(Sys.getenv("R_HOME"), "/bin/i386/Rscript.exe ", shQuote("C:\\Users\\Gough Conservation\\Documents\\Gough Birders\\2018-2019\\12.Monthly reports 2018-19\\RODBC_imports.r")), wait = FALSE, invisible = FALSE)

try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdBreedingSuccess"), silent=T)
#try(setwd("C:\\Users\\Gough Conservation\\Documents\\Gough Birders\\2018-2019\\12.Monthly reports 2018-19"), silent=T)

load("GOUGH_nest_data.RData")
head(nestsDB)  ## nest data
head(visDB)  ## nest visit data


##  SELECT DATA FOR TARGET SPECIES AND SUMMARISE NEST SUCCESS ####

head(nestsDB)
succ<-nestsDB %>% filter(Species=="MAPR") %>% filter(Year>2013) %>%
  mutate(count=1) %>%
  group_by(Species,Year) %>%
  summarise(R=sum(count),J=sum(SUCCESS))

rm(contacts,nestsDB,visDB)
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\MAPR")
# save.image("MAPR_IPM_input_data.RData")
# load("MAPR_IPM_input_data.RData")
sum(succ$R)


#########################################################################
# 3. Specify BASIC POPULATION MODEL WITH TWO SCENARIOS
#########################################################################

### DEMOGRAPHIC PARAMETERS 

#Juvenile survival: 	0.728 	from Barbraud & Weimerskirch (2003), Oro et al. (2004)
#Immature survival: 	0.894 	from Barbraud & Weimerskirch (2003)
#Adult survival: 	0.894 	from Barbraud & Weimerskirch (2003)
#Age at maturity: 	4 	from	Warham (1990), Oro et al. (2004)
#Female breeding success: 	0.519 	from Nevoux & Barbraud (2005)

### Calculation of stable age distribution 
### CREATING THE POPULATION MATRIX ###
# 
# seabird.matrix<-matrix(c(
#   0,0,0,0,0.519*0.5,
#   0.728,0,0,0,0,
#   0,0.894,0,0,0,
#   0,0,0.894,0,0,
#   0,0,0,0.894,0.894),ncol=5, byrow=T)
# stable.stage(seabird.matrix)
# 


setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\MAPR")
sink("MAPR_IPM_REV1_varfec.jags")
cat("
  
  
  model {
    #-------------------------------------------------
    # - population model for the MacGillivray's Prion population
    # - age structured model with 4 age classes 
    # - adult survival based on CMR ringing data with m-array and temporal variation
    # - productivity based on Prion Cave nest monitoring data
    # - TWO future scenarios to project population growth with and without eradication
    # -------------------------------------------------
    
    #-------------------------------------------------  
    # 1. PRIORS FOR ALL DATA SETS
    #-------------------------------------------------
    
    
    # -------------------------------------------------        
    # 1.1. Priors and constraints FOR FECUNDITY
    # -------------------------------------------------
    
    mean.fec[1] ~ dunif(0,1)         ## uninformative prior for BAD YEARS
    mean.fec[2] ~ dunif(0,1)         ## uninformative prior for GOOD YEARS
    prop.good ~ dunif(0,0.3)        ## proportion of years that is good or bad (to allow past variation when good years were more common)
    orig.fec ~ dunif(0.88,0.94)        ## uninformative prior for ORIGINAL FECUNDITY in proportion of years with good (similar to 2016) fecundity
    full.fec ~ dnorm(0.519,100) T(0.1,1)     ## prior for full fecundity without predation from Nevoux & Barbraud (2005) - very high precision
    fec.decrease <- (prop.good-orig.fec)/(58-0)   ## 58 years elapsed between original pop size data in 1957 and start of productivity time series in 2014
    
    # -------------------------------------------------        
    # 1.2. Priors and constraints FOR SURVIVAL
    # -------------------------------------------------
    
    # -------------------------------------------------
    # Parameters:
    # phi: survival probability for adults
    # p: recapture probability when breeding
    # emigrate: probability to emigrate into inaccessible part of Prion Cave
    # -------------------------------------------------
    
    # Priors and constraints
    for (t in 1:(n.occasions-1)){
      #logit(phi[t]) <- mu.phi + surv.raneff[t]
      #surv.raneff[t] ~ dnorm(0, tau.phi)
      p[t] ~ dunif(0, 1)
      emigrate[t] ~ dunif(0,0.25)
    }
    
    mean.phi ~ dunif(0, 1)             # Prior for mean survival
    #juv.surv.prop ~ dnorm(mean.juv.surv.prop,1000) T(0,1)
    mean.juv.surv ~ dunif(0.63,0.78)    ## based on juvenile survival for Balearic shearwaters in the Med.
    #mean.juv.surv ~ dunif(0.70,0.85)    ## based on juvenile survival for Balearic shearwaters in the Med, and Grey-faced Petrels

    
    #-------------------------------------------------  
    # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
    #-------------------------------------------------
    
    # -------------------------------------------------        
    # 2.1. System process: female based matrix model
    # -------------------------------------------------
    for (scen in 1:2){

      ### INITIAL VALUES FOR COMPONENTS FOR YEAR 1 - based on stable stage distribution from previous model
      fec.proj[scen,1]<-mean.fec[year.prop.good[scen,1]+1]    ## takes good or bad year fecundity
      year.prop.good[scen,1] ~ dbern(orig.fec)
      
      ### INITIAL VALUES FOR COMPONENTS FOR YEAR 1 - based on stable stage distribution from previous model
      JUV[1,scen]<-round(Ntot.breed[1,scen]*0.5*(mean.fec[year.prop.good[scen,1]+1]))
      N1[1,scen]<-round(Ntot.breed[1,scen]*0.5*(mean.fec[year.prop.good[scen,1]+1])*mean.juv.surv)
      N2[1,scen]<-round(Ntot.breed[1,scen]*0.5*(mean.fec[year.prop.good[scen,1]+1])*mean.juv.surv*mean.phi)
      N3[1,scen]<-round(Ntot.breed[1,scen]*0.5*(mean.fec[year.prop.good[scen,1]+1])*mean.juv.surv*mean.phi*mean.phi)
      Ntot.breed[1,scen] ~ dunif(2000000,5000000)         # initial value of population size
      
      for (tt in 2:65){

        ## LINEARLY DECREASING PROBABILITY OF A GOOD YEAR FROM 1956 to 2014
        year.fec.prop[scen,tt]<- max(0,min(1,(orig.fec + fec.decrease*tt))) ## calculate yearly proportion of good breeding year, but constrain to 0-1 to avoid invalid parent value
        year.prop.good[scen,tt] ~ dbern(year.fec.prop[scen,tt])
        fec.proj[scen,tt]<-mean.fec[year.prop.good[scen,tt]+1]    ## takes good or bad year fecundity
        breed.prop[scen,tt] ~ dunif(0.85,0.95)    ## breeding propensity
        
        ## THE PRE-BREEDERS ##
        JUV[tt,scen] ~ dbin(fec.proj[scen,tt],round(0.5 * Ntot.breed[tt,scen]*breed.prop[scen,tt]))                                   ### number of locally produced FEMALE chicks
        N1[tt,scen]  ~ dbin(mean.juv.surv, max(2,round(JUV[tt-1,scen])))                                               ### number of 1-year old survivors 
        N2[tt,scen] ~ dbin(mean.phi, max(2,round(N1[tt-1,scen])))                                                      ### number of 2-year old survivors
        N3[tt,scen] ~ dbin(mean.phi, max(2,round(N2[tt-1,scen])))                                                      ### number of 3-year old survivors
        
        ## THE BREEDERS ##
        Ntot.breed[tt,scen] ~ dbin(mean.phi, max(2,round(N3[tt-1,scen]+Ntot.breed[tt-1,scen])))                            ### the annual number of breeding birds is the sum of old breeders and recent recruits
        
      } # tt
      
      for (tt in 66:PROJ){
        
        ## SELECT GOOD OR BAD OR RODENT FREE FECUNDITY FOR FUTURE
        year.fec.prop[scen,tt]<- min(1,max(0,(orig.fec + fec.decrease*tt))) ## calculate yearly proportion of good breeding year, but constrain to 0-1 to avoid invalid parent value
        year.prop.good[scen,tt] ~ dbern(year.fec.prop[scen,tt])
        fec.proj[scen,tt]<-max(mean.fec[year.prop.good[scen,tt]+1],(scen-1)*full.fec)    ## takes current fecundity for scenario 1 and full fecundity for scenario 2
        breed.prop[scen,tt] ~ dunif(0.85,0.95)    ## breeding propensity
        
        ## THE PRE-BREEDERS ##
        JUV[tt,scen] ~ dbin(fec.proj[scen,tt],round(0.5 * Ntot.breed[tt,scen]*breed.prop[scen,tt]))                                   ### need a discrete number otherwise dbin will fail, dpois must be >0
        N1[tt,scen]  ~ dbin(mean.juv.surv, max(2,round(JUV[tt-1,scen])))                                               ### number of 1-year old survivors 
        N2[tt,scen] ~ dbin(mean.phi, max(2,round(N1[tt-1,scen])))                                                      ### number of 2-year old survivors
        N3[tt,scen] ~ dbin(mean.phi, max(2,round(N2[tt-1,scen])))                                                      ### number of 3-year old survivors
        
        ## THE BREEDERS ##
        Ntot.breed[tt,scen] ~ dbin(mean.phi, max(2,round(N3[tt-1,scen]+Ntot.breed[tt-1,scen])))                            ### the annual number of breeding birds is the sum of old breeders and recent recruits
        
      } # tt
      
    } # scen    
    


    
    # -------------------------------------------------        
    # 2.2. Likelihood for fecundity: Poisson regression from the number of surveyed broods
    # -------------------------------------------------
    for (t in 1:(T.fec)){      ### T-1 or not
      J[t] ~ dpois(rho.fec[t])
      rho.fec[t] <- R[t]*mean.fec[goodyear[t]+1]
      goodyear[t] ~ dbern(prop.good)
    } #	close loop over every year in which we have fecundity data
    
    
    # -------------------------------------------------        
    # 2.3. Likelihood for adult survival from multi-event model
    # -------------------------------------------------
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
    
    
    # -------------------------------------------------        
    # 4. DERIVED PARAMETERS
    # -------------------------------------------------

    ## DERIVED POPULATION GROWTH RATE 
    for (scen in 1:2){
      for (tt in 1:33){
        lambda[tt,scen]<-Ntot.breed[tt+67,scen]/max(1,Ntot.breed[tt+66,scen])
        loglam[tt,scen]<-log(lambda[tt,scen])
      } ## end of tt
      
      growth.rate[scen] <- exp((1/(33))*sum(loglam[1:(33),scen]))  ### geometric mean growth rate
      
    } ## end of scen
    
  }  ## END MODEL LOOP
    
    
    
    ",fill = TRUE)
sink()





#########################################################################
# 4. SET UP AND RUN INTEGRATED POPULATION MODEL
#########################################################################

# Bundle data
jags.data <- list(## survival
  y = CH,
  f = f,
  n.occasions = dim(CH)[2],
  nind = dim(CH)[1],
  #mean.juv.surv.prop= 0.728/0.894,  ## juvenile survival based on proportion of adult survival from Jiguet 2007
  
  ## fecundity
  R =succ$R,
  J=succ$J,
  T.fec=length(succ$J),
  goodyear=c(0,0,1,0,0,0),
  
  ## population process
  Ntot.obs=matrix(c(3500000,800000,3500000,800000), ncol=2), ### adjusted for v6 to 2-5 million
  PROJ=66+36)  ### adjusted for v6 to 103 years

# Initial values 
inits <- function(){list(mean.phi = runif(1, 0.7, 1),
                         p = runif(dim(CH)[2]-1, 0, 1),
                         orig.fec= runif(1, 0.75, 0.85))}  ### adjusted for REV1 as frequency of good years


# Parameters monitored
#parameters <- c("orig.fec","mean.fec","fec.decrease","fec.drop","mean.juv.surv","mean.phi","growth.rate","lambda","Ntot.breed")
parameters <- c("orig.fec","mean.fec","fec.decrease","prop.good","mean.juv.surv","mean.phi","growth.rate","lambda","Ntot.breed")

# MCMC settings
ni <- 150000
nt <- 10
nb <- 50000
nc <- 3

# Call JAGS from R (model created below)
MAPR_IPM <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\MAPR\\MAPR_IPM_REV1_varfec.jags",  ## changed from v4 to v6 on 10 Aug
                     n.chains = nc, n.thin = nt, n.burnin = nb,parallel=T,n.iter = ni)


### save model workspace
setwd("C:\\STEFFEN\\MANUSCRIPTS\\in_press\\MAPR_pop_model")
#save.image("MAPR_IPM_REV1.RData")
load("MAPR_IPM_REV1.RData")







#########################################################################
# 5. SUMMARISE OUTPUT AND PLOT POPULATION TRAJECTORY
#########################################################################
## compile output
out<-as.data.frame(MAPR_IPM$summary)
out$parameter<-row.names(MAPR_IPM$summary)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE OUTPUT TABLE FOR MANUSCRIPT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(out)

TABLE1<-out %>% filter(parameter %in% c('mean.fec[1]','mean.fec[2]','orig.fec','prop.good','fec.decrease','mean.phi','mean.juv.surv','growth.rate[1]','growth.rate[2]')) %>%
  select(parameter,c(5,3,7))

names(TABLE1)<-c("Parameter","Median","lowerCL","upperCL")
TABLE1$Parameter<-c("original proportion of good breeding year (1956)","productivity (poor year)","productivity (good year)","annual decline in frequency of good breeding year","current proportion of good breeding year (2014-2019)","first year survival probability","annual adult survival probability","annual population growth rate (no eradication)","annual population growth rate (with eradication)")
TABLE1
#fwrite(TABLE1,"MAPR_demographic_parameter_estimates_REV1.csv")

## FORMAT TABLE FOR MANUSCRIPT

TABLE1<-TABLE1 %>% mutate(MED=paste(round(Median,3)," (",round(lowerCL,3)," - ", round(upperCL,3),")", sep="")) %>%
  select(Parameter,MED) %>%
  rename(`Median (95% credible interval)`=MED)
setwd("C:\\STEFFEN\\MANUSCRIPTS\\submitted\\MAPR_pop_model")
#fwrite(TABLE1,"TABLE1.csv")

## REPORT QUANTITIES FOR RESULTS SECTION
sum(succ$R)
dim(CH)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GRAPH 1: POPULATION TRAJECTORY UNDER BOTH SCENARIOS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## retrieve the past population estimates (2006-2019)
MAPRpop<-out[(grep("Ntot.breed\\[",out$parameter)),c(12,5,4,6)] %>%
  mutate(Year=rep(seq(1956,2057),2)) %>%
  mutate(scenario=as.numeric(str_extract_all(parameter,"\\(?[0-9]+\\)?", simplify=TRUE)[,2])) %>%
  mutate(Scenario=ifelse(scenario==1,"no eradication","with eradication")) %>%
  filter(!(Scenario=="with eradication" & Year<2024)) %>%
  #filter((Scenario=="with eradication")) %>%
  #rename(parm=parameter,median=`50%`,lcl=`2.5%`,ucl=`97.5%`) %>%
  rename(parm=parameter,median=`50%`,lcl=`25%`,ucl=`75%`) %>%
  dplyr::select(parm,Scenario,Year,median,lcl,ucl)


### summary for manuscript
MAPRpop %>% filter(Year==2020)
MAPRpop %>% filter(Year==1956)
174684/3502083

### CREATE PLOT FOR BASELINE TRAJECTORY
MAPRpop$ucl[MAPRpop$ucl>5000000]<-4999999

ggplot()+
  geom_line(data=MAPRpop, aes(x=Year, y=median, color=Scenario), size=1)+
  geom_ribbon(data=MAPRpop,aes(x=Year, ymin=lcl,ymax=ucl, fill=Scenario),alpha=0.2)+
  
  ## format axis ticks
  scale_y_continuous(name="MacGillivray's Prion pairs (millions)", limits=c(0,5000000),breaks=seq(0,5000000,500000),labels=seq(0,5,0.5))+
  scale_x_continuous(name="Year", limits=c(1956,2057), breaks=seq(1956,2056,20), labels=as.character(seq(1956,2056,20)))+
  
  ## add count data
  geom_segment(aes(x=1956, xend=1956,y=0.4*5000000,yend=0.5*10000000),lineend = "round", size=2, colour="darkblue") +
  geom_segment(aes(x=2000, xend=2000,y=0.4*1500000,yend=0.5*2000000),lineend = "round", size=2, colour="darkblue") +
  
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=14, color="black"), 
        axis.title=element_text(size=18),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"),
        legend.position = c(0.73,0.89),
        legend.key = element_rect(fill = NA),
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))
ggsave("Figure1.tif", width = 9, height = 6, device='tiff', dpi=300)
ggsave("MAPR_population_projection_REV1_CI50.jpg", width=9, height=6)


### CREATE INSET PLOT FOR FUTURE PROJECTION
MAPRpop$ucl[MAPRpop$ucl>1000000]<-999999

ggplot()+
  geom_line(data=MAPRpop[MAPRpop$Year>2019,], aes(x=Year, y=median, color=Scenario), size=1)+
  geom_ribbon(data=MAPRpop[MAPRpop$Year>2019,],aes(x=Year, ymin=lcl,ymax=ucl, fill=Scenario),alpha=0.2)+
  
  ## format axis ticks
  scale_y_continuous(name="MacGillivray's Prion pairs (millions)", limits=c(0,1000000),breaks=seq(0,1000000,100000),labels=seq(0,1,0.1))+
  scale_x_continuous(name="Year", limits=c(2020,2057), breaks=seq(2021,2056,5), labels=as.character(seq(2021,2056,5)))+
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=14, color="black"), 
        axis.title=element_text(size=18),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"),
        legend.position = c(0.2,0.85),
        legend.key = element_rect(fill = NA),
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))
ggsave("MAPR_population_projection_REV1_CI50_INSET.jpg", width=9, height=6)






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GRAPH 2: EXTINCTION PROBABILITY OVER TIME
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### EXTRACT AND COMBINE DATA FROM ALL CHAINS 
selcol<-grep("Ntot.breed",dimnames(MAPR_IPM$samples[[1]])[[2]])   ## FIND COLUMS WE NEED
allchainsamples <- data.frame()
for(chain in 1:3) {
  samplesout<-as.data.frame(MAPR_IPM$samples[[chain]][,selcol]) %>% gather(key="parm", value="value")
  allchainsamples <- rbind(allchainsamples,as.data.frame(samplesout))
}

  
### CALCULATE EXTINCTION PROBABILITY
extprop <- allchainsamples %>%
    mutate(scen.index=as.numeric(str_extract_all(parm,"\\(?[0-9]+\\)?", simplify=TRUE)[,2])) %>%
    mutate(Scenario=ifelse(scen.index==1,"no eradication","with eradication")) %>%
    mutate(Year=as.numeric(str_extract_all(parm,"\\(?[0-9]+\\)?", simplify=TRUE)[,1])+1955) %>%
    
    mutate(n=1, inc=ifelse(value<500,1,0)) %>%   ### DEFINE EXTINCTION PROBABILITY HERE
    group_by(Scenario,Year) %>%
    summarise(ext.prob=sum(inc)/sum(n)) %>%
    filter(Year>2019)
  
head(allchainsamples)
head(extprop)
dim(extprop)


## CREATE A COLOUR PALETTE FOR THE NUMBER OF CHICKS RELEASED
colfunc <- colorRampPalette(c("firebrick","cornflowerblue"))


ggplot(data=extprop)+
  geom_line(aes(x=Year, y=ext.prob, color=Scenario), size=1)+
  
  ## format axis ticks
  scale_y_continuous(name="Probability of extinction (%)", limits=c(0,0.35),breaks=seq(0,0.35,0.05), labels=as.character(seq(0,35,5)))+
  scale_x_continuous(name="Year", breaks=seq(2020,2055,5), labels=as.character(seq(2020,2055,5)))+
  guides(color=guide_legend(title="Scenario"),fill=guide_legend(title="Scenario"))+
  scale_colour_manual(palette=colfunc)+
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=14, color="black"),
        axis.text.x=element_text(size=14, color="black",angle=45, vjust = 1, hjust=1),
        axis.title=element_text(size=18),
        legend.text=element_text(size=14, color="black"),
        legend.title=element_text(size=14, color="black"),
        legend.key = element_rect(fill = NA),
        legend.position = c(0.15,0.90),
        strip.text.x=element_text(size=14, color="black"),
        strip.text.y=element_text(size=14, color="black"),
        strip.background=element_rect(fill="white", colour="black"))
ggsave("Figure2.tif", width = 9, height = 6, device='tiff', dpi=300)
ggsave("MAPR_extinction_probability_REV1_250.jpg", width=9, height=6)



### CREATE TABLE 2 FOR MANUSCRIPT ###

head(extprop)
TABLE2<- extprop %>%
  filter(Year==2057) %>%
  select(ext.prob,Scenario) %>%
  mutate(ext.prob=ext.prob*100) %>%
  mutate(ext.prob=ifelse(ext.prob<1,"< 1%",paste0(ext.prob,"%")))
fwrite(TABLE2,"TABLE2.csv")






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PROPORTION OF SIMULATIONS WITH NEGATIVE GROWTH RATE AFTER MOUSE ERADICATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

### EXTRACT AND COMBINE DATA FROM ALL CHAINS 
futlamsamples <- data.frame()
for(chain in 1:3) {
  samplesout<-as.data.frame(MAPR_IPM$samples[[chain]][,9]) %>% gather(key="parm", value="value")
  futlamsamples <- rbind(futlamsamples,as.data.frame(samplesout))
}
head(futlamsamples)

futlamsamples %>% mutate(decline=ifelse(value<0.95,1,0)) %>%
  tabyl(decline)






##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## SUMMARISE WIND AND RAIN FOR MAPR BREED SUCCESS YEARS
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## count number of precipitation-free days

try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\DATA\\Weather"), silent=T)
weather<-fread("Gough_weather_data_2000_2020.csv")
head(weather)
weathersummary<-weather %>% group_by(Year) %>%
  summarise(mean_wind=mean(Windspeed, na.rm=T),max_wind=max(Windspeed, na.rm=T),rain=sum(Rain)) %>%
  filter(Year %in% succ$Year)
weathersummary

### McClelland et al. 2018 used 'rain free days' but that metric is suspicious in 2014!
rainsummary<-weather %>% group_by(Year, Month, Day) %>%
  summarise(rain=sum(Rain)) %>%
  mutate(rainfree=ifelse(rain==0,1,0)) %>%
  ungroup() %>%
  group_by(Year) %>%
  summarise(rainfreedays=sum(rainfree)) %>%
  filter(Year %in% succ$Year)
rainsummary



##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## PRODUCE WEATHER CHART
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

ggplot(weathersummary, mapping = aes(x = Year, y = mean_wind, group = 1)) + 
  geom_bar(mapping = aes(y = rain/100), stat = "identity", color="cornflowerblue", fill="cornflowerblue", width = 0.5) + 
  geom_line(color="indianred", size=1.5) +
  

  ## format axis ticks
  scale_y_continuous("Mean wind speed (m/s)", sec.axis = sec_axis(~.*100, name = "Precipitation (mm)")) +
  scale_x_continuous(name="Year", limits=c(2013.6,2019.4), breaks=seq(2014,2019,1)) +
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text=element_text(size=14, color="black"),
        axis.title=element_text(size=18))
