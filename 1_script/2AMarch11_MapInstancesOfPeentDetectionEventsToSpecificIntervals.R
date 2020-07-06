library(lubridate)
library(intrval)
library(dplyr)
library(tidyr)
#Read in the bird data you will use 
PeentEvents<-read.csv("0_data/raw/PeentEvents.csv", header=TRUE)
str(PeentEvents)#61597 observations out of 61933 observation event-days 
#(including days when no Peents were detected at a site)
PeentEvents.NARM<-PeentEvents[!is.na(PeentEvents$peentevent_sec),]

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

PeentEvents.NARM$month<-ifelse(PeentEvents.NARM$month<10,paste0("0",PeentEvents.NARM$month),PeentEvents.NARM$month)
PeentEvents.NARM$day<-ifelse(PeentEvents.NARM$day<10,paste0("0",PeentEvents.NARM$day),PeentEvents.NARM$day)
PeentEventTimes<-merge(PeentEvents.NARM, Times, by=c("sm_id","year","month","day"))
str(PeentEventTimes)

PeentEventTimes$eventID <- paste0(PeentEventTimes$sm_id, PeentEventTimes$startime, PeentEventTimes$peentevent_sec)
PeentEventTimes$eventTime <- PeentEventTimes$starttime.posix+PeentEventTimes$peentevent_sec
PeentEventDF <- PeentEventTimes[,c("sm_id","date","eventID","eventTime")]
write.csv(PeentEventDF, file="0_data/processed/PeentEventDForiginal.csv")

#First off, the start and end times in the period DF should be converted to lubridate intervals: 
Times.10min.int<-read.csv("0_data/processed/1_IntervalUsed/Times.10min.interval.csv",header=TRUE)
Times.10min.int$start_time_new<-as.POSIXct(Times.10min.int$start_time_new, tz = "", format = "%Y-%m-%d %H:%M:%S")# %p
Times.10min.int$end_time_new<-as.POSIXct(Times.10min.int$end_time_new, tz = "", format = "%Y-%m-%d %H:%M:%S")# %p

Times.10min.int$intervalsP <- as.interval(Times.10min.int$start_time_new, Times.10min.int$end_time_new)


#Step 3: The function can the be used on the event DF and both DFs can then be merged to produce the DF I was looking for:
#PeentEventDF$intervals <- lapply(PeentEventDF$eventTime, PeriodAssign, Times.10min.int$intervalsP)

#make sure both PeentEventDF and Times.10min.int are sorted by site and date
Peent.sortedbysite<- PeentEventDF[order(PeentEventDF$sm_id),] 
Period.sortedbysite<- Times.10min.int[order(Times.10min.int$sm_id),] 

Peent.sortedbysite$sm_id<-as.factor(Peent.sortedbysite$sm_id)
Period.sortedbysite$sm_id<-as.factor(Period.sortedbysite$sm_id)
Peent.sortedbysite$date<-as.factor(Peent.sortedbysite$date)
Period.sortedbysite$date<-as.factor(Period.sortedbysite$date)

sm_id_list<-levels(Peent.sortedbysite$sm_id)
# [1] "7524"  "7620"  "7698"  "7699"  "7720"  "7734"  "7755"  "7758" 
#[9] "8245"  "8246"  "9991"  "12121" "13482" "13512" "13523" "13525"
#[17] "13577" "13585" :Note that PeentEventDF starts with sm_id 7524 and Times.10min.int

sm_id_listB<-levels(Period.sortedbysite$sm_id)
# [1] "7524"  "7620"  "7689"  "7698"  "7699"  "7720"  "7734"  "7755" 
#[9] "7758"  "8245"  "8246"  "9990"  "9991"  "12011" "12121" "12122"
#[17] "13482" "13512" "13523" "13525" "13577" "13585" "13588"

#Calls weren't detected at all sites surveyed
Period.sortedbysiteB<-Period.sortedbysite[Period.sortedbysite$sm_id %in% c("7524","7620","7698","7699","7720","7734","7755","7758","8245","8246","9991","12121","13482","13512","13523","13525","13577","13585"),]
sm_id_listB<-levels(Period.sortedbysiteB$sm_id)
nrow(Period.sortedbysite)#28863
nrow(Period.sortedbysiteB)#22582
write.csv(Period.sortedbysiteB, file="0_data/processed/2_Period.sortedbysiteB.csv")
Period.sortedbysiteB<-read.csv("0_data/processed/2_Period.sortedbysiteB.csv", header=TRUE)
Period.sortedbysiteB$sm_id<-as.factor(Period.sortedbysiteB$sm_id)
sm_id_listB<-levels(Period.sortedbysiteB$sm_id)

Peentassignment<-list()#start with empty list for storing results
Period.sortedbysiteB$start_time_new<-as.POSIXct(Period.sortedbysiteB$start_time_new)
Period.sortedbysiteB$end_time_new<-as.POSIXct(Period.sortedbysiteB$end_time_new)

#list of all possible sites with Peent detection events
for (i in sm_id_list){
  BE1<-Peent.sortedbysite[Peent.sortedbysite$sm_id==i,]
  print(paste0("Peent detection events sorted at  site ",i)) 
  INT1<-Period.sortedbysiteB[Period.sortedbysiteB$sm_id==i,]
  print(paste0("Intervals sorted at  site ",i)) 
  date_list<-levels(droplevels(BE1$date))
  #list of dates from a given site i (with Peent detections, not necessarily all sites)
 for (j in date_list){
  BE2<-BE1[BE1$date==j,]
  if (nrow(BE2)>0){
    print(paste0("Peent detections from site ",i," on a given date ",j))
  }
  
  INT2<-INT1[INT1$date==j,]#SHOULD ONLY BE USING DATES WHERE DETECTIONS OCCURRED AT SITE I

  for (k in seq_len(nrow(BE2))){
    row_k<-BE2[k,]
    L<-row_k$eventTime  %[)% INT2[,c("start_time_new","end_time_new")]
    intervals<-as.character(INT2$intervalsP)[L]
    #NOTE: I had to use the interval to character coercion, because otherwise intervals were coerced to their length in seconds by the apply function and as such being not really useful for matching purposes - i.e. all four intervals in this example are the same length
    IDEventInterval<-paste0(row_k$eventID,"####",row_k$eventTime,"####",intervals)
      
    Peentassignment<-append(Peentassignment, IDEventInterval)
  }
  
  print(paste0("detections assigned from site ",i," on date ",j," to intervals in site ",i," on date ",j))
  #assign detections from site i on date j to appropriate interval from 
  #site i on date j
 }
}
IDEventInterval.df<-data.frame(unlist(Peentassignment))#
IDEventInterval.finaldf<-IDEventInterval.df %>% separate(unlist.Peentassignment., c("eventID", "eventTime", "interval"), sep="####")
str(IDEventInterval.finaldf)

Peent.sortedbysite$intervals<-IDEventInterval.finaldf$interval

Times.10min.int$intervals <- as.character(Times.10min.int$intervalsP)
mergedDF <- merge(Times.10min.int, Peent.sortedbysite, by = c("sm_id","date","intervals"), all.x=TRUE)
write.csv(mergedDF, file="0_data/processed/2_PeentDetectionsMapped/PeentDetections.10min.int.csv")

interval.length<-c("1min","2min","3min","4min","5min","6min","7min",
                   "8min","9min","10min","11min","12min","13min","14min",
                   "15min","16min","17min","18min","19min","20min","1hour")

for (H in interval.length){
  #First off, the start and end times in the period DF should be converted to lubridate intervals: 
  Times.temp<-read.csv(paste0("0_data/processed/1_IntervalUsed/Times.",H,".interval.csv"),header=TRUE)
  Times.temp$start_time_new<-as.POSIXct(Times.temp$start_time_new, tz = "", format = "%Y-%m-%d %H:%M:%S")# %p
  Times.temp$end_time_new<-as.POSIXct(Times.temp$end_time_new, tz = "", format = "%Y-%m-%d %H:%M:%S")# %p
  Times.temp$intervalsP <- as.interval(Times.temp$start_time_new, Times.temp$end_time_new)
  
  #make sure both PeentEventDF and Times.temp are sorted by site and date
  Peent.sortedbysite<- PeentEventDF[order(PeentEventDF$sm_id),] 
  Period.sortedbysite<- Times.temp[order(Times.temp$sm_id),] 
  
  Peent.sortedbysite$sm_id<-as.factor(Peent.sortedbysite$sm_id)
  Period.sortedbysite$sm_id<-as.factor(Period.sortedbysite$sm_id)
  Peent.sortedbysite$date<-as.factor(Peent.sortedbysite$date)
  Period.sortedbysite$date<-as.factor(Period.sortedbysite$date)
  
  Period.sortedbysiteB<-Period.sortedbysite[Period.sortedbysite$sm_id %in% sm_id_list,]
  sm_id_listB<-levels(droplevels(Period.sortedbysiteB$sm_id))

  Peentassignment<-list()#start with empty list for storing results
  sm_id_list<-levels(Peent.sortedbysite$sm_id)
  sm_id_listB<-levels(Period.sortedbysite$sm_id)
  Period.sortedbysiteB$start_time_new<-as.POSIXct(Period.sortedbysiteB$start_time_new)
  Period.sortedbysiteB$end_time_new<-as.POSIXct(Period.sortedbysiteB$end_time_new)
  
  #list of all possible sites with Peent detection events
  for (i in sm_id_list){
    BE1<-Peent.sortedbysite[Peent.sortedbysite$sm_id==i,]
    print(paste0("Peent detection events sorted at  site ",i)) 
    INT1<-Period.sortedbysiteB[Period.sortedbysiteB$sm_id==i,]
    print(paste0("Intervals sorted at  site ",i)) 
    date_list<-levels(droplevels(BE1$date))
    #list of dates from a given site i (with Peent detections, not necessarily all sites)
    for (j in date_list){
      BE2<-BE1[BE1$date==j,]
      if (nrow(BE2)>0){
        print(paste0("Peent detections from site ",i," on a given date ",j))
      }
      
      INT2<-INT1[INT1$date==j,]#SHOULD ONLY BE USING DATES WHERE DETECTIONS OCCURRED AT SITE I
      
      for (k in seq_len(nrow(BE2))){
        row_k<-BE2[k,]
        L<-row_k$eventTime  %[)% INT2[,c("start_time_new","end_time_new")]
        intervals<-as.character(INT2$intervalsP)[L]
        #NOTE: I had to use the interval to character coercion, because otherwise intervals were coerced to their length in seconds by the apply function and as such being not really useful for matching purposes - i.e. all four intervals in this example are the same length
        IDEventInterval<-paste0(row_k$eventID,"####",row_k$eventTime,"####",intervals)
        
        Peentassignment<-append(Peentassignment, IDEventInterval)
      }
      
      print(paste0("detections assigned from site ",i," on date ",j," to intervals in site ",i," on date ",j))
      #assign detections from site i on date j to appropriate interval from 
      #site i on date j
    }
  }
  IDEventInterval.df<-data.frame(unlist(Peentassignment))#
  IDEventInterval.finaldf<-IDEventInterval.df %>% separate(unlist.Peentassignment., c("eventID", "eventTime", "interval"), sep="####")
  str(IDEventInterval.finaldf)
  
  Peent.sortedbysite$intervals<-IDEventInterval.finaldf$interval
  
  Times.temp$intervals <- as.character(Times.temp$intervalsP)
  mergedDF <- merge(Times.temp, Peent.sortedbysite, by = c("sm_id","date","intervals"), all.x=TRUE)
  write.csv(mergedDF, file=paste0("0_data/processed/2_PeentDetectionsMapped/PeentDetections.",H,".int.csv"))
  
}
