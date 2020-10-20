##########################################################################
#
# MACGILLIVRAY PRION POPULATION MODEL FOR GOUGH ISLAND
#
##########################################################################

# population model adapted from Jiguet, F., Robert, A., Micol, T. and Barbraud, C. (2007), Quantifying stochastic and deterministic threats to island seabirds: last endemic prions face extinction from falcon peregrinations. Animal Conservation, 10: 245-253. doi:10.1111/j.1469-1795.2007.00100.x
# implemented in JAGS based on Kery and Schaub 2012
# written by Steffen.oppel@rspb.org.uk in May 2020

# revised in June 2020: start projection in 1955 with 5 million pairs, 2001 had 1 million pairs, project to 2050 (100 years)
# attempted multi-event formulation for survival to account for transients, but yielded very low adult survival - not included (see MAPR_survival_transients.r in C:\STEFFEN\RSPB\UKOT\Gough\ANALYSIS\SeabirdSurvival)
# retained capture history for only non-transients (birds observed at least twice)


# updated on 22 June 2020 after receiving advice from Adam Butler on how to include the count data
# added stochastic count node and adjusted other priors as suggested
# DID NOT CONVERGE - reverted to very narrow prior on orig.fec

# finalised on 15 July 2020 by including 2014 productivity data

# revised on 10 August 2020 to include Peter Ryan's comments
# adjusted input data to 2-5 million and expanded projection to 36 years (=3 generations)

# REVISION in OCTOBER 2020:
# updated CMR data input - removed chicks and changed encounter occasion (from year to season)
# switched to m-array to allow GoF test for survival model - but used this only for GoF test
# included temporal variation in phi and p in survival model
# included all transients and switched to a multi-event survival model

library(tidyverse)
library(jagsUI)
library(data.table)
library(lubridate)
library(popbio)
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

# Compute vector with occasion of first capture
get.first <- function(x) min(which(x==1))
f <- apply(CH, 1, get.first)



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
save.image("MAPR_IPM_input_data.RData")
load("MAPR_IPM_input_data.RData")

succ<-succ[1,] %>% mutate(Year=2014,R=63,J=0) %>% bind_rows(succ) %>% arrange(Year)
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
sink("MAPR_IPM_REV1.jags")
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
    
    mean.fec ~ dunif(0,0.5)         ## uninformative prior for CURRENT FECUNDITY
    orig.fec ~ dunif(0.3,0.35)         ## uninformative prior for ORIGINAL FECUNDITY
    full.fec ~ dnorm(0.519,100) T(0.1,1)     ## prior for full fecundity without predation from Nevoux & Barbraud (2005) - very high precision
    fec.decrease <- (mean.fec-orig.fec)/(66-1) ## 66 years elapsed between original pop size data in 1957 and start of productivity time series in 2014
    
    
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
    mean.juv.surv ~ dunif(0.65, 0.80) 
    
    #-------------------------------------------------  
    # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
    #-------------------------------------------------
    
    # -------------------------------------------------        
    # 2.1. System process: female based matrix model
    # -------------------------------------------------
    for (scen in 1:2){
      
      ### INITIAL VALUES FOR COMPONENTS FOR YEAR 1 - based on stable stage distribution from previous model
      JUV[1,scen]<-round(Ntot.breed[1,scen]*0.5*(orig.fec))
      N1[1,scen]<-round(Ntot.breed[1,scen]*0.5*(orig.fec)*mean.juv.surv)
      N2[1,scen]<-round(Ntot.breed[1,scen]*0.5*(orig.fec)*mean.juv.surv*mean.phi)
      N3[1,scen]<-round(Ntot.breed[1,scen]*0.5*(orig.fec)*mean.juv.surv*mean.phi*mean.phi)
      Ntot.breed[1,scen] ~ dunif(2000000,5000000)         # initial value of population size
      
      for (tt in 2:65){

        ## DECREASING FECUNDITY FROM UNKNOWN ORIGINAL in 1956 to known current fecundity in 2015
        fec.proj[scen,tt]<-orig.fec + fec.decrease*tt    ## models linear decrease of fecundity
        
        ## THE PRE-BREEDERS ##
        JUV[tt,scen] ~ dbin(fec.proj[scen,tt],round(0.5 * Ntot.breed[tt,scen]))                                   ### number of locally produced FEMALE chicks
        N1[tt,scen]  ~ dbin(mean.juv.surv, max(2,round(JUV[tt-1,scen])))                                               ### number of 1-year old survivors 
        N2[tt,scen] ~ dbin(mean.phi, max(2,round(N1[tt-1,scen])))                                                      ### number of 2-year old survivors
        N3[tt,scen] ~ dbin(mean.phi, max(2,round(N2[tt-1,scen])))                                                      ### number of 3-year old survivors
        
        ## THE BREEDERS ##
        Ntot.breed[tt,scen] ~ dbin(mean.phi, max(2,round(N3[tt-1,scen]+Ntot.breed[tt-1,scen])))                            ### the annual number of breeding birds is the sum of old breeders and recent recruits
        
      } # tt
      
      for (tt in 66:PROJ){
        
        ## SELECT CURRENT OR RODENT FREE FECUNDITY FOR FUTURE
        fec.proj[scen,tt]<-max(mean.fec,(scen-1)*full.fec)    ## takes current fecundity for scenario 1 and full fecundity for scenario 2 
        
        ## THE PRE-BREEDERS ##
        JUV[tt,scen] ~ dbin(fec.proj[scen,tt],round(0.5 * Ntot.breed[tt,scen]))                                   ### need a discrete number otherwise dbin will fail, dpois must be >0
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
      rho.fec[t] <- R[t]*mean.fec
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
  
  ## population process
  Ntot.obs=matrix(c(3500000,800000,3500000,800000), ncol=2), ### adjusted for v6 to 2-5 million
  PROJ=66+36)  ### adjusted for v6 to 103 years

# Initial values 
inits <- function(){list(mean.phi = runif(1, 0.7, 1),
                         p = runif(dim(CH)[2]-1, 0, 1),
                         orig.fec= runif(1, 0.38, 0.40))}  ### adjusted for v6 to 0.25-0.35


# Parameters monitored
parameters <- c("orig.fec","mean.fec","fec.decrease","fec.drop","mean.juv.surv","mean.phi","growth.rate","lambda","Ntot.breed")

# MCMC settings
ni <- 1500000
nt <- 10
nb <- 50000
nc <- 3

# Call JAGS from R (model created below)
MAPR_IPM <- jags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\MAPR\\MAPR_IPM_REV1.jags",  ## changed from v4 to v6 on 10 Aug
                     n.chains = nc, n.thin = nt, n.burnin = nb,parallel=T,n.iter = ni)


### save model workspace
setwd("C:\\STEFFEN\\MANUSCRIPTS\\submitted\\MAPR_pop_model")
save.image("MAPR_IPM_REV1.RData")
#load("MAPR_IPM_REV1.RData")







#########################################################################
# 5. SUMMARISE OUTPUT AND PLOT POPULATION TRAJECTORY
#########################################################################
## compile output
out<-as.data.frame(MAPR_IPM$summary)
out$parameter<-row.names(MAPR_IPM$summary)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE OUTPUT TABLE FOR REPORT /MANUSCRIPT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(out)

TABLE1<-out %>% filter(parameter %in% c('orig.fec','mean.fec','fec.decrease','fec.drop','mean.phi','mean.juv.surv','growth.rate[1]','growth.rate[2]')) %>%
  select(parameter,c(5,3,7))

names(TABLE1)<-c("Parameter","Median","lowerCL","upperCL")
TABLE1$Parameter<-c("original productivity (1956)","current productivity (2014-2019)","productivity decrease (1956-2014)","first year survival probability","annual adult survival probability","annual population growth rate (no eradication)","annual population growth rate (with eradication)")
TABLE1
setwd("C:\\STEFFEN\\MANUSCRIPTS\\submitted\\MAPR_pop_model")
fwrite(TABLE1,"MAPR_demographic_parameter_estimates_v4b.csv")


## REPORT QUANTITIES FOR RESULTS SECTION
sum(succ$R)
dim(CH)
length(del[which(del==1)])
mean(del[which(del>1)])


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
171867/3498443

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
colfunc <- colorRampPalette(c("cornflowerblue", "firebrick"))


ggplot(data=extprop)+
  geom_line(aes(x=Year, y=ext.prob, color=Scenario), size=1)+
  
  ## format axis ticks
  scale_y_continuous(name="Probability of extinction (%)", limits=c(0,0.25),breaks=seq(0,0.25,0.05), labels=as.character(seq(0,25,5)))+
  scale_x_continuous(name="Year", breaks=seq(2020,2055,5), labels=as.character(seq(2020,2055,5)))+
  guides(color=guide_legend(title="Scenario"),fill=guide_legend(title="Scenario"))+
  scale_colour_manual(palette=colfunc)+
  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=14, color="black"),
        axis.text.x=element_text(size=12, color="black",angle=45, vjust = 1, hjust=1),
        axis.title=element_text(size=18),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"),
        legend.key = element_rect(fill = NA),
        strip.text.x=element_text(size=14, color="black"),
        strip.text.y=element_text(size=14, color="black"),
        strip.background=element_rect(fill="white", colour="black"))

ggsave("MAPR_extinction_probability_REV1_250.jpg", width=9, height=6)



### CREATE TABLE 2 FOR MANUSCRIPT ###

head(extprop)
TABLE2<- extprop %>%
  filter(Year==2057) %>%
  select(ext.prob,Scenario) %>%
  mutate(ext.prob=ext.prob*100) %>%
  mutate(ext.prob=ifelse(ext.prob<1,"< 1%",paste0(ext.prob,"%")))
fwrite(TABLE2,"TABLE2.csv")





##### PLAYING WITH PRIORS
orig.fec<- 0.35
mean.fec=0.08
fec.drop<-mean.fec/orig.fec

fec.drop ~ dunif(0,1)      ## uninformative prior for DECREASE IN FECUNDITY
orig.fec ~ dunif(0,1)         ## uninformative prior for ORIGINAL FECUNDITY
full.fec ~ dnorm(0.519,100) T(0.1,1)     ## prior for full fecundity without predation from Nevoux & Barbraud (2005) - very high precision
mean.fec <- orig.fec * (1 - fec.drop)
fec.decrease <- fec.drop * orig.fec / (66 - 1)