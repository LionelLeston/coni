library(lubridate)
library(tidyverse)
library(intrval)
library(dplyr)
library(tidyr)
library(intervals)
#Read in the temperature data you will use
Temperature<-read.csv("0_data/raw/Temperature.csv")
str(Temperature)
Temperature$month<-ifelse(Temperature$month<10,paste0("0",Temperature$month),Temperature$month)
Temperature$day<-ifelse(Temperature$day<10,paste0("0",Temperature$day),Temperature$day)
Temperature$date<-as.Date(paste0(Temperature$year,"-",Temperature$month,"-",Temperature$day))
MeanTemp<-Temperature%>%
  group_by(sm_id, date)%>%
  summarise(meantemp=mean(temp_c))
MeanTemp$date<-as.factor(MeanTemp$date)

interval.length<-c("1min","2min","3min","4min","5min","6min","7min",
                   "8min","9min","10min","11min","12min","13min","14min",
                   "15min","16min","17min","18min","19min","20min","1hour")

#interval.length<-c("1hour")
for (H in interval.length){
  #First off, the start and end times in the period DF should be converted to lubridate intervals: 
  Times<-read.csv(paste0("0_data/processed/3_BoomsMapped_SunAndMoon/SunMoonBoomDetections.",H,".int.csv"),header=TRUE)
  Times.temp<-merge(Times, MeanTemp, by=c("sm_id","date"), all.x=TRUE)
  print(paste0("Temperatures added to ",H," ARU intervals"))
  write.csv(Times.temp, file=paste0("0_data/processed/3_BoomsMapped_SunAndMoon/TempSunMoonBoomDetections.",H,".int.csv"))
  #there will still be some NA values for some sites on some nights when temperature wasn't measured
}
