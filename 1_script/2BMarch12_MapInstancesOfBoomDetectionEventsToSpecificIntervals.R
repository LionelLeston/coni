library(lubridate)
library(intrval)
library(dplyr)
library(tidyr)

#Read in the bird data you will use 
BoomEvents<-read.csv("0_data/raw/BoomEvents.csv", header=TRUE)
str(BoomEvents)#12873 observations out of 13194 observation event-days 
#(including days when no Booms were detected at a site)
BoomEvents.NARM<-BoomEvents[!is.na(BoomEvents$boomevent_sec),]

#You will need to merge the Time data to the detection events to get the 
#start time of each recording to calculate the event time of each detection.
Times<-read.csv("0_data/raw/CONIlat_Times.csv", header=TRUE)
str(Times)
#PreSS contains the start-time of each recording
#   Convert Date-PreSS to a POSIXct object
Times$month<-ifelse(Times$month<10,paste0("0",Times$month),Times$month)
Times$day<-ifelse(Times$day<10,paste0("0",Times$day),Times$day)
Times$date<-as.Date(paste0(Times$year,"-",Times$month,"-",Times$day))
Times$startime<-paste0(Times$date," ",Times$PreSS)
Times$starttime.posix<-as.POSIXct(Times$startime, tz = "", format = "%Y-%m-%d %H:%M:%S")# %p

BoomEvents.NARM$month<-ifelse(BoomEvents.NARM$month<10,paste0("0",BoomEvents.NARM$month),BoomEvents.NARM$month)
BoomEvents.NARM$day<-ifelse(BoomEvents.NARM$day<10,paste0("0",BoomEvents.NARM$day),BoomEvents.NARM$day)
BoomEventTimes<-merge(BoomEvents.NARM, Times, by=c("sm_id","year","month","day"))
str(BoomEventTimes)

BoomEventTimes$eventID <- paste0(BoomEventTimes$sm_id, BoomEventTimes$startime, BoomEventTimes$boomevent_sec)
BoomEventTimes$eventTime <- BoomEventTimes$starttime.posix+BoomEventTimes$boomevent_sec
BoomEventDF <- BoomEventTimes[,c("sm_id","date","eventID","eventTime")]
write.csv(BoomEventDF, file="0_data/processed/BoomEventDForiginal.csv")

#example of mapping boom events to one file
#First off, the start and end times in the period DF should be converted to lubridate intervals: 
Times.14min.int<-read.csv("0_data/processed/1_IntervalUsed/Times.14min.interval.csv",header=TRUE)
Times.14min.int$start_time_new<-as.POSIXct(Times.14min.int$start_time_new, tz = "", format = "%Y-%m-%d %H:%M:%S")# %p
Times.14min.int$end_time_new<-as.POSIXct(Times.14min.int$end_time_new, tz = "", format = "%Y-%m-%d %H:%M:%S")# %p

Times.14min.int$intervalsP <- as.interval(Times.14min.int$start_time_new, Times.14min.int$end_time_new)


#make sure both BoomEventDF and Times.14min.int are sorted by site and date
Boom.sortedbysite<- BoomEventDF[order(BoomEventDF$sm_id),] 
Period.sortedbysite<- Times.14min.int[order(Times.14min.int$sm_id),] 

Boom.sortedbysite$sm_id<-as.factor(Boom.sortedbysite$sm_id)
Period.sortedbysite$sm_id<-as.factor(Period.sortedbysite$sm_id)
Boom.sortedbysite$date<-as.factor(Boom.sortedbysite$date)
Period.sortedbysite$date<-as.factor(Period.sortedbysite$date)

sm_id_list<-levels(Boom.sortedbysite$sm_id)
# [1] "7524"  "7620"  "7689"  "7698"  "7699"  "7720"  "7734"  "7755" 
#[9] "7758"  "8245"  "8246"  "12121" "12122" "13482" "13512" "13523"
#[17] "13525" "13577" "13588"

sm_id_listB<-levels(Period.sortedbysite$sm_id)
# [1] "7524"  "7620"  "7689"  "7698"  "7699"  "7720"  "7734"  "7755" 
#[9] "7758"  "8245"  "8246"  "9990"  "9991"  "12011" "12121" "12122"
#[17] "13482" "13512" "13523" "13525" "13577" "13585" "13588"

#Booms weren't detected at all sites surveyed
Period.sortedbysiteB<-Period.sortedbysite[Period.sortedbysite$sm_id %in% c("7524","7620","7689","7698","7699","7720","7734","7755","7758","8245","8246","12121","12122","13482","13512","13523","13525","13577","13588"),]
sm_id_listB<-levels(Period.sortedbysiteB$sm_id)
nrow(Period.sortedbysite)#20697
nrow(Period.sortedbysiteB)#17147
 
sm_id_listB<-levels(droplevels(Period.sortedbysiteB$sm_id))

Boomassignment<-list()#start with empty list for storing results

#list of all possible sites with Boom detection events
for (i in sm_id_list){
  BE1<-Boom.sortedbysite[Boom.sortedbysite$sm_id==i,]
  print(paste0("Boom detection events sorted at  site ",i)) 
  INT1<-Period.sortedbysiteB[Period.sortedbysiteB$sm_id==i,]
  print(paste0("Intervals sorted at  site ",i)) 
  date_list<-levels(droplevels(BE1$date))
  #list of dates from a given site i (with Boom detections, not necessarily all sites)
 for (j in date_list){
  BE2<-BE1[BE1$date==j,]
  if (nrow(BE2)>0){
    print(paste0("Boom detections from site ",i," on a given date ",j))
  }
  
  INT2<-INT1[INT1$date==j,]#SHOULD ONLY BE USING DATES WHERE DETECTIONS OCCURRED AT SITE I
 
  for (k in seq_len(nrow(BE2))){
    row_k<-BE2[k,]
    L<-row_k$eventTime  %[)% INT2[,c("start_time_new","end_time_new")]
      intervals<-as.character(INT2$intervalsP)[L]
      #NOTE: I had to use the interval to character coercion, because otherwise intervals were coerced to their length in seconds by the apply function and as such being not really useful for matching purposes - i.e. all four intervals in this example are the same length
      IDEventInterval<-paste0(row_k$eventID,"####",row_k$eventTime,"####",intervals)
      
      Boomassignment<-append(Boomassignment, IDEventInterval)
  }
  
  print(paste0("detections assigned from site ",i," on date ",j," to intervals in site ",i," on date ",j))
  #assign detections from site i on date j to appropriate interval from 
  #site i on date j
 }
}
 
IDEventInterval.df<-data.frame(unlist(Boomassignment))#
IDEventInterval.finaldf<-IDEventInterval.df %>% separate(unlist.Boomassignment., c("eventID", "eventTime", "interval"), sep="####")
str(IDEventInterval.finaldf)

Boom.sortedbysite$intervals<-IDEventInterval.finaldf$interval

write.csv(Boom.sortedbysite, file="0_data/processed/2_BoomDetectionsMapped/BoomEvent.sorted14min.csv")

Times.14min.int$intervals <- as.character(Times.14min.int$intervalsP)
mergedDF <- merge(Times.14min.int, Boom.sortedbysite, by = c("sm_id","date","intervals"), all.x=TRUE)
write.csv(mergedDF, file="0_data/processed/2_BoomDetectionsMapped/BoomDetections.14min.int.csv")

interval.length<-c("1min","2min","3min","4min","5min","6min","7min",
                   "8min","9min","10min","11min","12min","13min","14min",
                   "15min","16min","17min","18min","19min","20min","1hour")

for (H in interval.length){
  #First off, the start and end times in the period DF should be converted to lubridate intervals: 
  Times.temp<-read.csv(paste0("0_data/processed/1_IntervalUsed/Times.",H,".interval.csv"),header=TRUE)
  Times.temp$start_time_new<-as.POSIXct(Times.temp$start_time_new, tz = "", format = "%Y-%m-%d %H:%M:%S")# %p
  Times.temp$end_time_new<-as.POSIXct(Times.temp$end_time_new, tz = "", format = "%Y-%m-%d %H:%M:%S")# %p
  Times.temp$intervalsP <- as.interval(Times.temp$start_time_new, Times.temp$end_time_new)
  

  #make sure both BoomEventDF and Times.temp are sorted by site and date
  Boom.sortedbysite<- BoomEventDF[order(BoomEventDF$sm_id),] 
  Period.sortedbysite<- Times.temp[order(Times.temp$sm_id),] 
  
  Boom.sortedbysite$sm_id<-as.factor(Boom.sortedbysite$sm_id)
  sm_id_list<-levels(Boom.sortedbysite$sm_id)
  Period.sortedbysite$sm_id<-as.factor(Period.sortedbysite$sm_id)
  Boom.sortedbysite$date<-as.factor(Boom.sortedbysite$date)
  Period.sortedbysite$date<-as.factor(Period.sortedbysite$date)

  Period.sortedbysiteB<-Period.sortedbysite[Period.sortedbysite$sm_id %in% sm_id_list,]
  sm_id_listB<-levels(droplevels(Period.sortedbysiteB$sm_id))
  
  Boomassignment<-list()#start with empty list for storing results
  sm_id_list<-levels(Boom.sortedbysite$sm_id)
  sm_id_listB<-levels(Period.sortedbysite$sm_id)
  
  #list of all possible sites with Boom detection events
  for (i in sm_id_list){
    BE1<-Boom.sortedbysite[Boom.sortedbysite$sm_id==i,]
    print(paste0("Boom detection events sorted at  site ",i)) 
    INT1<-Period.sortedbysiteB[Period.sortedbysiteB$sm_id==i,]
    print(paste0("Intervals sorted at  site ",i)) 
    date_list<-levels(droplevels(BE1$date))
    #list of dates from a given site i (with Boom detections, not necessarily all sites)
    for (j in date_list){
      BE2<-BE1[BE1$date==j,]
      if (nrow(BE2)>0){
        print(paste0("Boom detections from site ",i," on a given date ",j))
      }
      
      INT2<-INT1[INT1$date==j,]#SHOULD ONLY BE USING DATES WHERE DETECTIONS OCCURRED AT SITE I
      
      for (k in seq_len(nrow(BE2))){
        row_k<-BE2[k,]
        L<-row_k$eventTime  %[)% INT2[,c("start_time_new","end_time_new")]
        intervals<-as.character(INT2$intervalsP)[L]
        #NOTE: I had to use the interval to character coercion, because otherwise intervals were coerced to their length in seconds by the apply function and as such being not really useful for matching purposes - i.e. all four intervals in this example are the same length
        IDEventInterval<-paste0(row_k$eventID,"####",row_k$eventTime,"####",intervals)
        
        Boomassignment<-append(Boomassignment, IDEventInterval)
      }
      
      print(paste0("detections assigned from site ",i," on date ",j," to intervals in site ",i," on date ",j))
      #assign detections from site i on date j to appropriate interval from 
      #site i on date j
    }
  }
  IDEventInterval.df<-data.frame(unlist(Boomassignment))#
  IDEventInterval.finaldf<-IDEventInterval.df %>% separate(unlist.Boomassignment., c("eventID", "eventTime", "interval"), sep="####")
  str(IDEventInterval.finaldf)
  
  Boom.sortedbysite$intervals<-IDEventInterval.finaldf$interval
  
  Times.temp$intervals <- as.character(Times.temp$intervalsP)
  mergedDF <- merge(Times.temp, Boom.sortedbysite, by = c("sm_id","date","intervals"), all.x=TRUE)
  write.csv(mergedDF, file=paste0("0_data/processed/2_BoomDetectionsMapped/BoomDetections.",H,".int.csv"))
  
}
