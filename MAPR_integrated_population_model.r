##########################################################################
#
# MACGILLIVRAY PRION POPULATION MODEL FOR GOUGH ISLAND
#
##########################################################################

# population model adapted from Jiguet, F., Robert, A., Micol, T. and Barbraud, C. (2007), Quantifying stochastic and deterministic threats to island seabirds: last endemic prions face extinction from falcon peregrinations. Animal Conservation, 10: 245-253. doi:10.1111/j.1469-1795.2007.00100.x
# implemented in JAGS based on Kery and Schaub 2012
# written by Steffen.oppel@rspb.org.uk in May 2020

# revised in June 2020: start projection in 1955 with 5 million pairs, 2001 had 1 million pairs, project to 2050 (100 years)


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

##### LOAD FORMATTED RINGING DATA ###########
setwd("C:/STEFFEN/RSPB/UKOT/Gough/ANALYSIS/SeabirdSurvival")

## run the RODBC import of CMR data in a 32-bit version of R
system(paste0("C:/PROGRA~1/R/R-35~1.1/bin/i386/Rscript.exe ", shQuote("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival\\RODBC_CMR_import.R")), wait = TRUE, invisible = FALSE, intern = T)
try(setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival"), silent=T)
load("GOUGH_seabird_CMR_data.RData")

###### FILTER DATA FROM RAW CONTACTS ########
contacts<-contacts %>% filter(SpeciesCode %in% c("MAPR","BBPR","PRIO")) %>% filter(Location=="Prion Cave")
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





#########################################################################
# 3. SET UP AND RUN INTEGRATED POPULATION MODEL
#########################################################################

# Bundle data
jags.data <- list(## survival
                  y = rCH,
                  f = f,
                  n.occasions = dim(rCH)[2],
                  nind = dim(rCH)[1],
                  mean.juv.surv.prop= 0.728/0.894,  ## juvenile psurvival based on proportion of adult survival from Jiguet 2007
                  
                  ## fecundity
                  R =succ$R,
                  J=succ$J,
                  T.fec=length(succ$J),
                  
                  ## population process
                  PROJ=100,
                  POP.SIZE=4500000
                  )

# Initial values 
inits <- function(){list(phi = runif(1, 0.7, 1),
                         pp = runif(1, 0, 1))}


# Parameters monitored
parameters <- c("mean.fec","full.fec","juv.surv","phi","p","growth.rate","lambda","Ntot.breed")

# MCMC settings
ni <- 10000
nt <- 2
nb <- 5000
nc <- 4

# Call JAGS from R (model created below)
MAPR_IPM <- autojags(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\MAPR\\MAPR_IPM_v2.jags",
                     n.chains = nc, n.thin = nt, n.burnin = nb,parallel=T) # n.iter = ni, 




 
  
#########################################################################
# 4. SUMMARISE OUTPUT AND PLOT POPULATION TRAJECTORY
#########################################################################
## compile output
out<-as.data.frame(MAPR_IPM$summary)
out$parameter<-row.names(MAPR_IPM$summary)



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CREATE OUTPUT TABLE FOR REPORT /MANUSCRIPT
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
head(out)

TABLE1<-out %>% filter(parameter %in% c('mean.fec','full.fec','phi','juv.surv','growth.rate[1]','growth.rate[2]')) %>%
  select(parameter,c(5,3,7))

names(TABLE1)<-c("Parameter","Median","lowerCL","upperCL")
TABLE1$Parameter<-c("current fecundity","mouse-free fecundity","first year survival","adult survival","population growth rate (no eradication)","population growth rate (with eradication)")
TABLE1
#fwrite(TABLE1,"MAPR_demographic_parameter_estimates_v1.csv")






#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GRAPH 1: POPULATION TRAJECTORY UNDER BOTH SCENARIOS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## retrieve the past population estimates (2006-2019)
MAPRpop<-out[(grep("Ntot.breed\\[",out$parameter)),c(12,5,3,7)] %>%
  mutate(Year=rep(seq(1956,2055),2)) %>%
  mutate(scenario=as.numeric(str_extract_all(parameter,"\\(?[0-9]+\\)?", simplify=TRUE)[,2])) %>%
  mutate(Scenario=ifelse(scenario==1,"no eradication","with eradication")) %>%
  filter(!(Scenario=="with eradication" & Year<2025)) %>%
  #filter((Scenario=="with eradication")) %>%
  rename(parm=parameter,median=`50%`,lcl=`2.5%`,ucl=`97.5%`) %>%
  dplyr::select(parm,Scenario,Year,median,lcl,ucl)



### CREATE PLOT FOR BASELINE TRAJECTORY
MAPRpop$ucl[MAPRpop$ucl>6000000]<-6000000

ggplot()+
  geom_line(data=MAPRpop, aes(x=Year, y=median, color=Scenario), size=1)+
  geom_ribbon(data=MAPRpop,aes(x=Year, ymin=lcl,ymax=ucl, fill=Scenario),alpha=0.2)+

  ## format axis ticks
  scale_y_continuous(name="N MacGillivray's Prion on Gough (1000s)", limits=c(0,6000000),breaks=seq(0,6000000,500000),labels=seq(0,6000,500))+
  scale_x_continuous(name="Year", limits=c(1956,2056), breaks=seq(1956,2056,20), labels=as.character(seq(1956,2056,20)))+
  
  ## add count data
  geom_segment(aes(x=1956, xend=1956,y=0.4*10000000,yend=0.5*10000000),lineend = "round", size=2, colour="darkblue") +
  geom_segment(aes(x=2000, xend=2000,y=0.4*1500000,yend=0.5*2000000),lineend = "round", size=2, colour="darkblue") +

  
  ## beautification of the axes
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=14, color="black"), 
        axis.title=element_text(size=18),
        legend.text=element_text(size=12, color="black"),
        legend.title=element_text(size=14, color="black"),
        legend.key = element_rect(fill = NA),
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))

ggsave("MAPR_population_projection_v2.jpg", width=9, height=6)








#########################################################################
# 5. Specify BASIC POPULATION MODEL WITH TWO SCENARIOS
#########################################################################

### DEMOGRAPHIC PARAMETERS 

#Juvenile survival: 	0.728 	from Barbraud & Weimerskirch (2003), Oro et al. (2004)
#Immature survival: 	0.894 	fromBarbraud & Weimerskirch (2003)
#Adult survival: 	0.894 	from Barbraud & Weimerskirch (2003)
#Age at maturity: 	4 	from	Warham (1990), Oro et al. (2004)
#Female breeding success: 	0.519 	from Nevoux & Barbraud (2005)


### Calculation of stable age distribution 
### CREATING THE POPULATION MATRIX ###

seabird.matrix<-matrix(c(
  0,0,0,0,0.519*0.5,
  0.728,0,0,0,0,
  0,0.894,0,0,0,
  0,0,0.894,0,0,
  0,0,0,0.894,0.894),ncol=5, byrow=T)
stable.stage(seabird.matrix)





setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\MAPR")
sink("MAPR_IPM_v2.jags")
cat("
    
    
model {
  #-------------------------------------------------
  # integrated population model for the MacGillivray's Prion population
  # - age structured model with 4 age classes 
  # - adult survival based on CMR ringing data
  # - productivity based on Prion Cave nest monitoring data
  # - simplified population process with informed prior for adults skipping breeding and uninformed immatures recruiting
  # - TWO future scenarios to project population growth with and without eradication
  # -------------------------------------------------
    
  #-------------------------------------------------  
  # 1. PRIORS FOR ALL DATA SETS
  #-------------------------------------------------
    
    
    # -------------------------------------------------        
    # 1.1. Priors and constraints FOR FECUNDITY
    # -------------------------------------------------
    
    mean.fec ~ dunif(0,1)         ## uninformative prior with upper bound from Nevoux & Barbraud (2005)
    full.fec ~ dnorm(0.519,1000)     ## prior for full fecundity without predation from Nevoux & Barbraud (2005) - very high precision

    
    # -------------------------------------------------        
    # 1.2. Priors and constraints FOR SURVIVAL
    # -------------------------------------------------
    
    phi ~ dunif(0.7, 1) 
    p ~ dunif(0, 1)
    juv.surv <- juv.surv.prop*phi
    juv.surv.prop ~ dnorm(mean.juv.surv.prop,1000) T(0,1)
    
  #-------------------------------------------------  
  # 2. LIKELIHOODS AND ECOLOGICAL STATE MODEL
  #-------------------------------------------------
    
    # -------------------------------------------------        
    # 2.1. System process: female based matrix model
    # -------------------------------------------------
    for (scen in 1:2){

      fec.proj[scen]<-max(mean.fec,(scen-1)*full.fec)    ## takes current fecundity for scenario 1 and full fecundity for scenario 2 

      ### INITIAL VALUES FOR COMPONENTS FOR YEAR 1 - based on stable stage distribution from previous model
    
      JUV[1,scen]<-max(2,round(Ntot.breed[1,scen]*0.5*(mean.fec+0.16)))
      N1[1,scen]<-round(Ntot.breed[1,scen]*0.5*(mean.fec+0.17)*juv.surv)
      N2[1,scen]<-round(Ntot.breed[1,scen]*0.5*(mean.fec+0.18)*juv.surv*phi)
      N3[1,scen]<-round(Ntot.breed[1,scen]*0.5*(mean.fec+0.19)*juv.surv*phi*phi)
      Ntot.breed[1,scen] ~ dnorm(POP.SIZE,10)         # initial value of population size

      for (tt in 2:66){
    
        ## THE PRE-BREEDERS ##
    
        nestlings[tt,scen] <- round((mean.fec+0.30/tt) * 0.5 * Ntot.breed[tt,scen])                                                    ### number of locally produced FEMALE chicks
        JUV[tt,scen] ~ dpois(nestlings[tt,scen])                                                                     ### need a discrete number otherwise dbin will fail, dpois must be >0
        N1[tt,scen]  ~ dbin(juv.surv, max(2,round(JUV[tt-1,scen])))                                                    ### number of 1-year old survivors 
        N2[tt,scen] ~ dbin(phi, max(2,round(N1[tt-1,scen])))                                                      ### number of 2-year old survivors
        N3[tt,scen] ~ dbin(phi, max(2,round(N2[tt-1,scen])))                                                       ### number of 3-year old survivors

    
        ## THE BREEDERS ##
    
        Ntot.breed[tt,scen] ~ dbin(phi, max(2,round(N3[tt-1,scen]+Ntot.breed[tt-1,scen])))                            ### the annual number of breeding birds is the sum of old breeders and recent recruits

      } # tt

      for (tt in 67:PROJ){
    
        ## THE PRE-BREEDERS ##
    
        nestlings[tt,scen] <- round(fec.proj[scen] * 0.5 * Ntot.breed[tt,scen])                                                    ### number of locally produced FEMALE chicks
        JUV[tt,scen] ~ dpois(nestlings[tt,scen])                                                                     ### need a discrete number otherwise dbin will fail, dpois must be >0
        N1[tt,scen]  ~ dbin(juv.surv, max(2,round(JUV[tt-1,scen])))                                                    ### number of 1-year old survivors 
        N2[tt,scen] ~ dbin(phi, max(2,round(N1[tt-1,scen])))                                                      ### number of 2-year old survivors
        N3[tt,scen] ~ dbin(phi, max(2,round(N2[tt-1,scen])))                                                       ### number of 3-year old survivors

    
        ## THE BREEDERS ##
    
        Ntot.breed[tt,scen] ~ dbin(phi, max(2,round(N3[tt-1,scen]+Ntot.breed[tt-1,scen])))                            ### the annual number of breeding birds is the sum of old breeders and recent recruits

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
    # 2.3. Likelihood for adult and juvenile survival from CMR
    # -------------------------------------------------
    
    for (i in 1:nind){
      # Define latent state at first capture
      z[i,f[i]] <- 1

        for (t in (f[i]+1):n.occasions){
          # State process
          z[i,t] ~ dbern(mu1[i,t])
          mu1[i,t] <- phi * z[i,t-1]

          # Observation process
          y[i,t] ~ dbern(mu2[i,t])
          mu2[i,t] <- p * z[i,t]
        } #t
    } #i


  # -------------------------------------------------        
  # 4. DERIVED POPULATION GROWTH RATE
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