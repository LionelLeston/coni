library(lubridate)
library(tidyverse)
#Read in the time data you will use 
Times<-read.csv("0_data/raw/CONIlat_Times.csv", header=TRUE)
str(Times)
#PreSS contains the start-time of each recording
#   Convert Date-PreSS to a POSIXct object
Times$month<-ifelse(Times$month<10,paste0("0",Times$month),Times$month)
Times$day<-ifelse(Times$day<10,paste0("0",Times$day),Times$day)
Times$date<-as.Date(paste0(Times$year,"-",Times$month,"-",Times$day))
Times$startime<-paste0(Times$date," ",Times$PreSS)
Times$starttime.posix<-as.POSIXct(Times$startime, tz = "", format = "%Y-%m-%d %H:%M:%S")# %p
Times$endtime.posix<-Times$starttime.posix+3600+Times$SStoSRend+3600
#SS contains the time of sunset
#SStoSR end contains the time in seconds from sunset to sunrise
#To get end-time of each recording add 3600 seconds to SStoSRend
#   Equals POSIXct(Date-PreSS)+3600+SStoSRend+3600

#Create an interval, convert it to an as.duration object, then split into
#smaller sub-duration objects
Times$fullrec.interval <- Times$starttime.posix %--% Times$endtime.posix
Times$fullrec.interval

Times$fullrec.duration <- as.duration(Times$fullrec.interval)
Times$fullrec.duration 
Times$fullrec.seconds<-Times$fullrec.duration@.Data  

#https://stackoverflow.com/questions/51407177/r-lubridate-split-durations-into-sub-durations
#sm_id replacing site_id
Times$sample_id<-paste0(Times$sm_id," ",Times$date)

Times.1minute.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=60) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.1minute.int), file="0_data/processed/1_IntervalUsed/Times.1min.interval.csv")

Times.2minute.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=120) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.2minute.int), file="0_data/processed/IntervalUsed/Times.2min.interval.csv")

Times.3minute.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=180) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.3minute.int), file="0_data/processed/IntervalUsed/Times.3min.interval.csv")

Times.4minute.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=240) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.4minute.int), file="0_data/processed/IntervalUsed/Times.4min.interval.csv")

Times.5minute.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=300) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.5minute.int), file="0_data/processed/IntervalUsed/Times.5min.interval.csv")

Times.6min.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=360) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.6min.int), file="0_data/processed/IntervalUsed/Times.6min.interval.csv")

Times.7min.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=420) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.7min.int), file="0_data/processed/IntervalUsed/Times.7min.interval.csv")

Times.8min.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=480) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.8min.int), file="0_data/processed/IntervalUsed/Times.8min.interval.csv")

Times.9min.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=540) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.9min.int), file="0_data/processed/IntervalUsed/Times.9min.interval.csv")

Times.10minute.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=600) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.10minute.int), file="0_data/processed/IntervalUsed/Times.10min.interval.csv")

Times.11minute.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=660) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.11minute.int), file="0_data/processed/IntervalUsed/Times.11min.interval.csv")

Times.12minute.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=720) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.12minute.int), file="0_data/processed/IntervalUsed/Times.12min.interval.csv")

Times.13minute.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=780) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.13minute.int), file="0_data/processed/IntervalUsed/Times.13min.interval.csv")

Times.14minute.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=840) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.14minute.int), file="0_data/processed/IntervalUsed/Times.14min.interval.csv")

Times.15min.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=900) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.15min.int), file="0_data/processed/IntervalUsed/Times.15min.interval.csv")

Times.16min.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=960) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.16min.int), file="0_data/processed/IntervalUsed/Times.16min.interval.csv")

Times.17min.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=1020) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.17min.int), file="0_data/processed/IntervalUsed/Times.17min.interval.csv")

Times.18min.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=1080) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.18min.int), file="0_data/processed/IntervalUsed/Times.18min.interval.csv")

Times.19min.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=1140) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.19min.int), file="0_data/processed/IntervalUsed/Times.19min.interval.csv")

Times.20min.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=1200) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.20min.int), file="0_data/processed/IntervalUsed/Times.20min.interval.csv")

Times.1hour.int<-Times %>% group_by(sample_id,sm_id) %>% 
  mutate(duration_new = (as.numeric(fullrec.duration)-1) %>% seq(0,.,by=3600) %>% c(fullrec.duration) %>% diff %>% list,
         start_time_new = list(starttime.posix + seconds(c(0,cumsum(head(duration_new[[1]],-1))))),
         end_time_new = list(starttime.posix + seconds(cumsum(duration_new[[1]]))),
         segment_id = list(seq_along(duration_new[[1]]))) %>% 
  unnest %>%
  ungroup
write.csv(as.data.frame(Times.1hour.int), file="0_data/processed/1_IntervalUsed/Times.1hour.interval.csv")

