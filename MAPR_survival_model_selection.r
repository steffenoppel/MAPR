##########################################################################
#
# MACGILLIVRAY PRION SURVIVAL ANALYSIS 2014-2018 - MODEL SELECTION
#
##########################################################################
# based on Kery and Schaub 2012, Chapter 7
# modified by Steffen oppel, October 2018
# REVISION 1, 16 OCT 2020: goodness-of-fit test, model selection
# switched to MARRAY because GoF test not reliable in state-space formulation
# model does not run in JAGS, so run in WinBUGS

# multi-event model including all transients now yields sensible estimates after revision on 19 Oct 2020
# model comparison now with full data and GoF test with full data to justify this model

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
  spread(key=Contact_Season,value=n,fill=0)
dim(EncHist)

#### FORMAT FOR SIMPLE CJS MODEL ############
CH<-as.matrix(EncHist[,2:ncol(EncHist)], dimnames=F)
CH<-ifelse(CH>0,1,0)  ##1=seen, 2=not seen


#### ELIMINATE TRANSIENTS ONLY OBSERVED IN A SINGLE YEAR
del <- apply(CH[,1:ncol(CH)], 1, sum)   ### adjusted to use only first year to exclude transients when ringing was set up - many loafers were ringed
dim(CH)
rCH<-CH[!(del==1),]
dim(rCH)


#### ELIMINATE ONLY TRANSIENTS MARKED AFTER MID DECEMBER AND ONLY OBSERVED IN A SINGLE YEAR
del <- apply(CH[,1:ncol(CH)], 1, sum)   ### adjusted to use only first year to exclude transients when ringing was set up - many loafers were ringed
del2<-EncHist$BirdID %in% latemarked$BirdID
dim(CH)
rCH<-CH[!(del==1 & del2==TRUE),]
dim(rCH)

#### ELIMINATE TRANSIENTS THAT WERE ONLY EVER CAUGHT ONCE
del3<-EncHist$BirdID %in% transients$BirdID
dim(CH)
rCH<-CH[!(del3==TRUE),]
dim(rCH)


#### Create Function to create a m-array based on capture-histories (CH)
marray <- function(CH){
  nind <- dim(CH)[1]
  n.occasions <- dim(CH)[2]
  m.array <- matrix(data = 0, ncol = n.occasions+1, nrow = n.occasions)
  # Calculate the number of released individuals at each time period
  for (t in 1:n.occasions){
    m.array[t,1] <- sum(CH[,t])
  }
  for (i in 1:nind){
    pos <- which(CH[i,]!=0)
    g <- length(pos)
    for (z in 1:(g-1)){
      m.array[pos[z],pos[z+1]] <- m.array[pos[z],pos[z+1]] + 1
    } #z
  } #i
  # Calculate the number of individuals that is never recaptured
  for (t in 1:n.occasions){
    m.array[t,n.occasions+1] <- m.array[t,1] - sum(m.array[t,2:n.occasions])
  }
  out <- m.array[1:(n.occasions-1),2:(n.occasions+1)]
  return(out)
}

# Create the m-array from the capture-histories
marr <- marray(CH)
rmarr <- marray(rCH)

# Parameters monitored
parameters <- c("phi", "p", "fit", "fit.new","mean.phi","sd.phi")

# MCMC settings
ni <- 15000
nt <- 3
nb <- 5000
nc <- 3

# MODEL WITH TIME VARIATION
inits <- function(){list(phi = runif(dim(marr)[2]-1, 0.7, 1), p = runif(dim(marr)[2]-1, 0, 1))}  
jags.data <- list(marr = marr, n.occasions = dim(marr)[2])
MAPRall <- bugs(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival\\MAPR_CJS_marray.bugs",
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = FALSE, bugs.directory = "C:\\STEFFEN\\Software\\WinBUGS14", working.directory = getwd())
MAPRall$summary


# MODEL WITH RANDOM TIME VARIATION
inits <- function(){list(mean.phi = runif(1, 0.7, 1), p = runif(dim(marr)[2]-1, 0, 1))}
parameters <- c("phi", "p", "mean.phi", "fit", "fit.new")
MAPRraneff <- bugs(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival\\MAPR_CJS_marray_temp_raneff.bugs",
                n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = FALSE, bugs.directory = "C:\\STEFFEN\\Software\\WinBUGS14", working.directory = getwd())
MAPRraneff$summary


# MODEL WITH CONSTANT PARAMETERS
parameters <- c("phi", "p", "fit", "fit.new")
inits <- function(){list(phi = runif(1, 0.7, 1), p = runif(1, 0, 1))}  
MAPRconstant <- bugs(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival\\MAPR_CJS_marray_constant.bugs",
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = FALSE, bugs.directory = "C:\\STEFFEN\\Software\\WinBUGS14", working.directory = getwd())

MAPRconstant$summary


# MODEL WITH CONSTANT SURVIVAL
inits <- function(){list(phi = runif(1, 0.7, 1), p = runif(dim(marr)[2]-1, 0, 1))}  
MAPRptime <- bugs(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival\\MAPR_CJS_marray_ptime.bugs",
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = FALSE, bugs.directory = "C:\\STEFFEN\\Software\\WinBUGS14", working.directory = getwd())

MAPRptime$summary


# MODEL WITH CONSTANT RECAPTURE
inits <- function(){list(phi = runif(dim(marr)[2]-1, 0.7, 1), p = runif(1, 0, 1))}  
MAPRphitime <- bugs(jags.data, inits, parameters, "C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\SeabirdSurvival\\MAPR_CJS_marray_phitime.bugs",
                   n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb, debug = FALSE, bugs.directory = "C:\\STEFFEN\\Software\\WinBUGS14", working.directory = getwd())

MAPRphitime$summary


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
########     LOAD ALL PREVIOUS MODEL RESULTS AND COMPARE MODELS    ###########################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

out1<-as.data.frame(MAPRall$summary)
out1$parameter<-row.names(MAPRall$summary)
out1$model<-"fulltime"

out2<-as.data.frame(MAPRphitime$summary)
out2$parameter<-row.names(MAPRphitime$summary)
out2$model<-"phitime"

out3<-as.data.frame(MAPRptime$summary)
out3$parameter<-row.names(MAPRptime$summary)
out3$model<-"ptime"

out4<-as.data.frame(MAPRconstant$summary)
out4$parameter<-row.names(MAPRconstant$summary)
out4$model<-"constant"

out5<-as.data.frame(MAPRraneff$summary)
out5$parameter<-row.names(MAPRraneff$summary)
out5$model<-"ran_time"


### COMBINE OUTPUT FROM ALL 4 MODELS
out<-bind_rows(out1,out2,out3,out4,out5)

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
DIC_tab<-bind_rows(pd_dic(MAPRall),pd_dic(MAPRphitime),pd_dic(MAPRptime),pd_dic(MAPRconstant),pd_dic(MAPRraneff)) %>%
  mutate(model=c("fulltime","phitime","ptime","constant", "ran_time")) %>%
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

fwrite(ModSelTab[,c(1,2,4:8)],"MAPR_surv_model_selection_table.csv")




#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# GOF TESTS of top models
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("C:\\STEFFEN\\RSPB\\UKOT\\Gough\\ANALYSIS\\PopulationModel\\MAPR")
load("MAPR_model_selection_GoF.RData")

# CONSTANT SURVIVAL
MASS::eqscplot(MAPRptime$sims.list$fit, MAPRptime$sims.list$fit.new, xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", las = 1,  
     xlim=range(MAPRptime$sims.list$fit, MAPRptime$sims.list$fit.new), ylim=range(MAPRptime$sims.list$fit, MAPRptime$sims.list$fit.new),bty ="n") 
abline(0, 1, lwd=2, col='red')
mean(MAPRptime$sims.list$fit.new > MAPRptime$sims.list$fit)

# FULL TIME
MASS::eqscplot(MAPRall$sims.list$fit, MAPRall$sims.list$fit.new, xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", las = 1,  
     xlim=range(MAPRall$sims.list$fit, MAPRall$sims.list$fit.new), ylim=range(MAPRall$sims.list$fit, MAPRall$sims.list$fit.new),bty ="n") 
abline(0, 1, lwd=2, col='red')
mean(MAPRall$sims.list$fit.new > MAPRall$sims.list$fit)

# RANDOM TIME
MASS::eqscplot(MAPRraneff$sims.list$fit, MAPRraneff$sims.list$fit.new, xlab = "Discrepancy actual data", ylab = "Discrepancy replicate data", las = 1,  
               xlim=range(MAPRraneff$sims.list$fit, MAPRraneff$sims.list$fit.new), ylim=range(MAPRraneff$sims.list$fit, MAPRraneff$sims.list$fit.new),bty ="n") 
abline(0, 1, lwd=2, col='red')
mean(MAPRraneff$sims.list$fit.new > MAPRraneff$sims.list$fit)




### MAKE FIGURE S1 FOR MANUSCRIPT
plotdata1<-data.frame(obs=MAPRptime$sims.list$fit, sim=MAPRptime$sims.list$fit.new, model="with transients")
plotdata2<-data.frame(obs=MAPRall$sims.list$fit, sim=MAPRall$sims.list$fit.new, model="no transients")

bind_rows(plotdata1,plotdata2) %>%
  ggplot() + geom_point(aes(x=obs,y=sim), size=1,colour="darkgrey") +
  facet_wrap(~model) +
  geom_abline(intercept=0, slope=1,colour="darkred", size=1.5)+
  scale_x_continuous(limits=c(0,15),breaks=seq(0,15,5), labels=seq(0,15,5))+
  scale_y_continuous(limits=c(0,15),breaks=seq(0,15,5), labels=seq(0,15,5))+
  xlab("Discrepancy observed data") +
  ylab("Discrepancy simulated data") +
  theme(panel.background=element_rect(fill="white", colour="black"), 
        axis.text=element_text(size=18, color="black"), 
        axis.title=element_text(size=20), 
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank())



#########################################################################
# Specify M-ARRAY CJS model with RANDOM TIME EFFECTS
#########################################################################

sink("MAPR_CJS_marray_temp_raneff.bugs")
cat("
    model {

    # Priors and constraints
    for (t in 1:(n.occasions-1)){
      phi[t] <- mu.phi+ surv.raneff[t]      # Priors for survival
      surv.raneff[t] ~ dnorm(0,tau.phi)
      p[t]  <- mu.p + capt.raneff[t]           # Priors for recapture
      capt.raneff[t] ~ dnorm(0,tau.p)
    }

    ## PRIORS FOR MEAN PARAMETERS
    mean.phi ~ dunif(0.7,1)
    mean.p ~ dunif(0,1)

    ## LOGIT-TRANSFORMATION
    mu.phi<-log(mean.phi/(1-mean.phi))
    mu.p<-log(mean.p/(1-mean.p))

    ## ERROR FOR RANDOM EFFECTS
    sigma.phi~dunif(0,2)
    tau.phi<-pow(sigma.phi,-2)
    sigma.p~dunif(0,10)
    tau.p<-pow(sigma.p,-2)

    
    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
      marr[t,1:n.occasions] ~ dmulti(pr[t, ], r[t])
    }
    
    # Calculate the number of birds released each year
    for (t in 1:(n.occasions-1)){
      r[t] <- sum(marr[t, ])
    }
    # Define the cell probabilities of the m-array
    # Main diagonal
    for (t in 1:(n.occasions-1)){
      q[t] <- 1-p[t]                # Probability of non-recapture
      pr[t,t] <- phi[t]*p[t]
    
      # Above main diagonal
      for (j in (t+1):(n.occasions-1)){
        pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
      } #j
    
      # Below main diagonal
      for (j in 1:(t-1)){
        pr[t,j] <- 0
      } #j
    } #t
    
    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
      pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
    } #t
    
    # Assess model fit using Freeman-Tukey statistic
    # Compute fit statistics for observed data
    for (t in 1:(n.occasions-1)){
      for (j in 1:n.occasions){
        expmarr[t,j] <- r[t]*pr[t,j]
        E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
    } #t
    
    # Generate replicate data and compute fit stats from them
    for (t in 1:(n.occasions-1)){
      marr.new[t,1:n.occasions] ~ dmulti(pr[t,], r[t])
    
      for (j in 1:n.occasions){
        E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
    } #t
    
    fit <- sum(E.org[,])
    fit.new <- sum(E.new[,])
    }
    ",fill = TRUE)
sink()






#########################################################################
# Specify M-ARRAY CJS model with TIME DEPENDENT PARAMETERS 
#########################################################################

sink("MAPR_CJS_marray.bugs")
cat("
    model {
    # Priors and constraints
    for (t in 1:(n.occasions-1)){
      phi[t] ~ dunif(0.7, 1)         # Priors for survival
      p[t] ~ dunif(0, 1)           # Priors for recapture
    }

    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
      marr[t,1:n.occasions] ~ dmulti(pr[t, ], r[t])
    }

    # Calculate the number of birds released each year
    for (t in 1:(n.occasions-1)){
      r[t] <- sum(marr[t, ])
    }
    # Define the cell probabilities of the m-array
    # Main diagonal
    for (t in 1:(n.occasions-1)){
      q[t] <- 1-p[t]                # Probability of non-recapture
      pr[t,t] <- phi[t]*p[t]
    
      # Above main diagonal
      for (j in (t+1):(n.occasions-1)){
        pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p[j]
      } #j

      # Below main diagonal
      for (j in 1:(t-1)){
        pr[t,j] <- 0
      } #j
    } #t

    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
      pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
    } #t
    
    # Assess model fit using Freeman-Tukey statistic
    # Compute fit statistics for observed data
    for (t in 1:(n.occasions-1)){
      for (j in 1:n.occasions){
        expmarr[t,j] <- r[t]*pr[t,j]
        E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
    } #t

    # Generate replicate data and compute fit stats from them
    for (t in 1:(n.occasions-1)){
      marr.new[t,1:n.occasions] ~ dmulti(pr[t,], r[t])

      for (j in 1:n.occasions){
        E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
    } #t

    fit <- sum(E.org[,])
    fit.new <- sum(E.new[,])
    mean.phi<-mean(phi[])
    sd.phi<-sd(phi[])
  }
    ",fill = TRUE)
sink()



#########################################################################
# Specify M-ARRAY CJS model with CONSTANT RECAPTURE 
#########################################################################

sink("MAPR_CJS_marray_phitime.bugs")
cat("
    model {
    # Priors and constraints
    for (t in 1:(n.occasions-1)){
      phi[t] ~ dunif(0.7, 1)         # Priors for survival
    }
    p ~ dunif(0, 1)           # Priors for recapture

    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
      marr[t,1:n.occasions] ~ dmulti(pr[t, ], r[t])
    }

    # Calculate the number of birds released each year
    for (t in 1:(n.occasions-1)){
      r[t] <- sum(marr[t, ])
    }
    # Define the cell probabilities of the m-array
    # Main diagonal
    for (t in 1:(n.occasions-1)){
      q[t] <- 1-p                # Probability of non-recapture
      pr[t,t] <- phi[t]*p
    
      # Above main diagonal
      for (j in (t+1):(n.occasions-1)){
        pr[t,j] <- prod(phi[t:j])*prod(q[t:(j-1)])*p
      } #j

      # Below main diagonal
      for (j in 1:(t-1)){
        pr[t,j] <- 0
      } #j
    } #t

    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
      pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
    } #t
    
    # Assess model fit using Freeman-Tukey statistic
    # Compute fit statistics for observed data
    for (t in 1:(n.occasions-1)){
      for (j in 1:n.occasions){
        expmarr[t,j] <- r[t]*pr[t,j]
        E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
    } #t

    # Generate replicate data and compute fit stats from them
    for (t in 1:(n.occasions-1)){
      marr.new[t,1:n.occasions] ~ dmulti(pr[t,], r[t])

      for (j in 1:n.occasions){
        E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
    } #t

    fit <- sum(E.org[,])
    fit.new <- sum(E.new[,])
  }
    ",fill = TRUE)
sink()




#########################################################################
# Specify M-ARRAY CJS model with CONSTANT PARAMETERS 
#########################################################################

sink("MAPR_CJS_marray_constant.bugs")
cat("
    model {
    # Priors and constraints
      phi ~ dunif(0.7, 1)         # Priors for survival
      p ~ dunif(0, 1)           # Priors for recapture

    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
      marr[t,1:n.occasions] ~ dmulti(pr[t, ], r[t])
    }

    # Calculate the number of birds released each year
    for (t in 1:(n.occasions-1)){
      r[t] <- sum(marr[t, ])
    }
    # Define the cell probabilities of the m-array
    # Main diagonal
    for (t in 1:(n.occasions-1)){
      q[t] <- 1-p                # Probability of non-recapture
      pr[t,t] <- phi*p
    
      # Above main diagonal
      for (j in (t+1):(n.occasions-1)){
        pr[t,j] <- pow(phi,j-t)*prod(q[t:(j-1)])*p
      } #j

      # Below main diagonal
      for (j in 1:(t-1)){
        pr[t,j] <- 0
      } #j
    } #t

    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
      pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
    } #t
    
    # Assess model fit using Freeman-Tukey statistic
    # Compute fit statistics for observed data
    for (t in 1:(n.occasions-1)){
      for (j in 1:n.occasions){
        expmarr[t,j] <- r[t]*pr[t,j]
        E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
    } #t

    # Generate replicate data and compute fit stats from them
    for (t in 1:(n.occasions-1)){
      marr.new[t,1:n.occasions] ~ dmulti(pr[t,], r[t])

      for (j in 1:n.occasions){
        E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
    } #t

    fit <- sum(E.org[,])
    fit.new <- sum(E.new[,])
  }
    ",fill = TRUE)
sink()



#########################################################################
# Specify M-ARRAY CJS model with CONSTANT SURVIVAL 
#########################################################################

sink("MAPR_CJS_marray_ptime.bugs")
cat("
    model {
    # Priors and constraints
    for (t in 1:(n.occasions-1)){
      p[t] ~ dunif(0, 1)           # Priors for recapture
    }
    phi ~ dunif(0.7, 1)         # Priors for survival

    # Define the multinomial likelihood
    for (t in 1:(n.occasions-1)){
      marr[t,1:n.occasions] ~ dmulti(pr[t, ], r[t])
    }

    # Calculate the number of birds released each year
    for (t in 1:(n.occasions-1)){
      r[t] <- sum(marr[t, ])
    }
    # Define the cell probabilities of the m-array
    # Main diagonal
    for (t in 1:(n.occasions-1)){
      q[t] <- 1-p[t]                # Probability of non-recapture
      pr[t,t] <- phi*p[t]
    
      # Above main diagonal
      for (j in (t+1):(n.occasions-1)){
        pr[t,j] <- pow(phi,j-t)*prod(q[t:(j-1)])*p[j]
      } #j

      # Below main diagonal
      for (j in 1:(t-1)){
        pr[t,j] <- 0
      } #j
    } #t

    # Last column: probability of non-recapture
    for (t in 1:(n.occasions-1)){
      pr[t,n.occasions] <- 1-sum(pr[t,1:(n.occasions-1)])
    } #t
    
    # Assess model fit using Freeman-Tukey statistic
    # Compute fit statistics for observed data
    for (t in 1:(n.occasions-1)){
      for (j in 1:n.occasions){
        expmarr[t,j] <- r[t]*pr[t,j]
        E.org[t,j] <- pow((pow(marr[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
    } #t

    # Generate replicate data and compute fit stats from them
    for (t in 1:(n.occasions-1)){
      marr.new[t,1:n.occasions] ~ dmulti(pr[t,], r[t])

      for (j in 1:n.occasions){
        E.new[t,j] <- pow((pow(marr.new[t,j], 0.5)-pow(expmarr[t,j], 0.5)), 2)
      } #j
    } #t

    fit <- sum(E.org[,])
    fit.new <- sum(E.new[,])
  }
    ",fill = TRUE)
sink()


