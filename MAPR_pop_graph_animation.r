##########################################################################
#
# MACGILLIVRAY PRION POPULATION MODEL GRAPH ANIMATION
#
##########################################################################

# based on Jones et al 2021: Mouse eradication is required to prevent local extinction of an endangered seabird on an oceanic island
# written by steffen.oppel@rspb.org.uk
# 22 January 2021

library(tidyverse)
library(data.table)
library(lubridate)
library(janitor)
library(gganimate)
filter<-dplyr::filter
select<-dplyr::select


#########################################################################
# 1. LOAD MODEL OUTPUT
#########################################################################

setwd("C:\\STEFFEN\\MANUSCRIPTS\\in_press\\MAPR_pop_model")
load("MAPR_IPM_REV1.RData")







#########################################################################
# 2. SUMMARISE DATA FOR PLOTTING
#########################################################################
## compile output
out<-as.data.frame(MAPR_IPM$summary)
out$parameter<-row.names(MAPR_IPM$summary)

## retrieve the past population estimates (2006-2019)
MAPRpop<-out[(grep("Ntot.breed\\[",out$parameter)),c(12,5,4,6)] %>%
  mutate(Year=rep(seq(1956,2057),2)) %>%
  mutate(scenario=as.numeric(str_extract_all(parameter,"\\(?[0-9]+\\)?", simplify=TRUE)[,2])) %>%
  mutate(Scenario=ifelse(scenario==1,"no eradication","with eradication")) %>%
  filter(!(Scenario=="with eradication" & Year<2024)) %>%
  rename(parm=parameter,median=`50%`,lcl=`25%`,ucl=`75%`) %>%
  dplyr::select(parm,Scenario,Year,median,lcl,ucl)




#########################################################################
# 3. CREATE ANIMATED PLOT
#########################################################################

MAPRanim<-ggplot()+
  ## add count data
  geom_segment(aes(x=1956, xend=1956,y=0.4*5000000,yend=0.5*10000000),lineend = "round", size=2, colour="darkblue") +
  geom_segment(aes(x=2000, xend=2000,y=0.4*1500000,yend=0.5*2000000),lineend = "round", size=2, colour="darkblue") +
  geom_text(aes(x=1956, y=5100000, label="Estimate 1956"),colour="darkblue", size=6,hjust=0,vjust=0) +
  geom_text(aes(x=2000, y=3700000, label="Estimate 2000"),colour="darkblue", size=6,hjust=0.5,vjust=1) +
  
  ## add TODAYS DATE
  geom_segment(aes(x=2021.2, xend=2021.2,y=0,yend=2500000),lineend = "round", size=1, colour="darkgrey", linetype='dotted') +
  geom_text(aes(x=2021.2, y=2700000, label="Today"),colour="darkgrey", size=6,hjust=0.5,vjust=1) +
  
  ## add the modelled population trajectory
  geom_line(data=MAPRpop, aes(x=Year, y=median, color=Scenario), size=1)+
  geom_ribbon(data=MAPRpop,aes(x=Year, ymin=lcl,ymax=ucl, fill=Scenario),alpha=0.2)+
  transition_reveal(along=Year)+
  
  ## format axis ticks
  scale_y_continuous(name="breeding pairs (millions)", limits=c(0,5250009),breaks=seq(0,5000000,500000),labels=seq(0,5,0.5))+
  scale_x_continuous(name="Year", limits=c(1956,2057), breaks=seq(1956,2056,10), labels=as.character(seq(1956,2056,10)))+

  
  ## beautification of the axes
  ggtitle("MacGillivray's Prion population on Gough Island") +
  theme(panel.background=element_rect(fill="white", colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.y=element_text(size=18, color="black"),
        axis.text.x=element_text(size=14, color="black"), 
        axis.title=element_text(size=18),
        plot.title = element_text(size = 32, face = "bold",hjust = 0.5),
        legend.text=element_text(size=16, color="darkblue"),
        legend.title=element_text(size=18, color="darkblue"),
        legend.position = c(0.73,0.89),
        legend.key = element_rect(fill = NA),
        strip.text.x=element_text(size=18, color="black"), 
        strip.background=element_rect(fill="white", colour="black"))


# animate in a two step process:
animate(MAPRanim, height = 600, width =900, fps=20, duration=7)
anim_save("MAPR_test.gif")
