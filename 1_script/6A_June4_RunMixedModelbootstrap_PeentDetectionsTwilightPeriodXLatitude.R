library(boot)
library(dplyr)
library(gamm4)
library(ggplot2)
library(grid)
library(gridExtra)
library(jpeg)
library(lme4)
library(lubridate)
library(mgcv)
library(MuMIn)
library(plotly)
library(pROC)
library(tidyr)
library(unmarked)


my.theme <- theme_classic() +
  theme(text=element_text(size=10, family="Arial"),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        axis.line.x=element_line(linetype=1),
        axis.line.y=element_line(linetype=1))

# function to obtain random samples stratified by site 
#cannot stratify by 2-week periods because after excluding sites
#without temperature readings one site (13577)
#only had samples from 2 nights within 1 2-week period
fSampleNint <- function(dat, intervals, n) {
  intervals <- enquo(intervals)
  dat %>%
    group_by(sm_id, twilightperiod)%>%
    filter(UQ(intervals) %in% sample(unique(UQ(intervals)), n)) %>%
    slice(sample(row_number()))
}#,season

#You need to run the mixed model first for later, when generating a
#model matrix for the prediction data sets
peentsint<-read.csv("0_data/processed/3_PeentsMapped_SunAndMoon/TempSunMoonPeentDetections.10min.int.csv", header=TRUE)
str(peentsint)
peentsint$PEENTS<-ifelse(is.na(peentsint$eventID),0,1)
#peentsint$halfmonth<-ifelse(peentsint$day<16,"1","2")
#peentsint$season<-as.factor(paste0(peentsint$month,"_",peentsint$halfmonth))

peentsint$start_time_new<- ymd_hms(as.character(peentsint$start_time_new))
peentsint$julian<-yday(as.Date(peentsint$date))
peentsint$hour<-hour(peentsint$start_time_new)
peentsint$minute<-minute(peentsint$start_time_new)
peentsint$start_time_numeric<-peentsint$hour+(peentsint$minute/60)

peentsint$TSSS.s<-scale(peentsint$TSSS, center=TRUE, scale=TRUE)#scaled to reduce rank deficiency of GAMMs
peentsint$ordinal.s<-scale(peentsint$julian, center=TRUE, scale=TRUE)#scaled to reduce rank deficiency of GAMMs
peentsint$TSSS.s2<-scale(peentsint$TSSS^2, center=FALSE, scale=TRUE)#scaled to reduce rank deficiency of GAMMs
peentsint$ordinal.s2<-scale(peentsint$julian^2, center=FALSE, scale=TRUE)#scaled to reduce rank deficiency of GAMMs
peentsint$latitude.s<-scale(peentsint$latitude, center=FALSE, scale=TRUE)#scaled to reduce rank deficiency of GAMMs
peentsint$meantemp.s<-scale(peentsint$meantemp, center=FALSE, scale=TRUE)#not currently being used in these GAMMs
#peentsint$moon.s<-scale(as.vector(peentsint$moon.fraction), center=TRUE, scale=FALSE)#not currently being used in these GAMMs

peentsint$sm_id<-as.factor(peentsint$sm_id)
peentsint$intervals<-as.factor(peentsint$intervals)

peentsintB<-peentsint[which(peentsint$julian<190&peentsint$julian>170),]

#peentsintB<-peentsint[!peentsint$month==8,]#drop August


#peentsintC<-peentsintB[!is.na(peentsintB$meantemp),]#drop intervals missing temperature peentsint

peentsintC<-peentsintB[!peentsintB$sm_id %in% c("7524","7698","8245","13523","13585","13588"),]
nrow(peentsintC)#69469
maxdur<-max(peentsintC$duration)
peentsintD<-peentsintC[peentsintC$duration==maxdur,]#gets rid of any "remainder" intervals
nrow(peentsintD)#69223

peentsintD$twilightperiodB<-peentsintD$twilightperiod
peentsintD$twilightperiodB<-ifelse(peentsintD$start_time_numeric>12,
                                   paste0(peentsintD$twilightperiodB,".1"),
                                   paste0(peentsintD$twilightperiodB,".2"))
peentsintD$twilightperiodB<-as.factor(peentsintD$twilightperiodB)
levels(peentsintD$twilightperiodB)

levels(peentsintD$twilightperiodB)[levels(peentsintD$twilightperiodB)=="AFTER.2"] <- "AFTER"
levels(peentsintD$twilightperiodB)[levels(peentsintD$twilightperiodB)=="BEFORE.1"] <- "BEFORE"
levels(peentsintD$twilightperiodB)[levels(peentsintD$twilightperiodB)=="NIGHT.2"] <- "NIGHT"


peentsintD$twilightperiodB<-factor(peentsintD$twilightperiodB, levels=c("BEFORE","CIVIL.1","NAUTICAL.1","ASTRONOMICAL.1","NIGHT","ASTRONOMICAL.2","NAUTICAL.2","CIVIL.2","AFTER"))
levels(peentsintD$twilightperiodB)

peentsintD$sm_id.interval<-as.factor(paste0(peentsintD$sm_id,peentsintD$intervals))
peentsintD.test <- fSampleNint(peentsintD, intervals, 10)
nrow(peentsintD.test)#4412;5335;5764;5955;3486;4278;4512
peentsintD.train  <- peentsintD[!peentsintD$sm_id.interval %in% peentsintD.test$sm_id.interval,]
nrow(peentsintD.train)#70745;69822;69393;69202;65737;8156;16702

levels(peentsintD.train$twilightperiodB)


#Mixed-effects model run once before bootstrap to get model matrix
d<-fSampleNint(peentsintD.train, intervals, 10)#13 or more intervals is too many sample without replacement

peents.pervisit<-d %>%
  group_by(sm_id, intervalsP) %>% 
  summarize(numdetect= sum(PEENTS), 
            duration=mean(duration_new),
            TSSS=mean(TSSScorr),
            TSSS.s=mean(TSSS.s),
            TSSS.s2=mean(TSSS.s2),  
            meantemp=mean(meantemp),
            meantemp.s=mean(meantemp.s),
            twilightperiod=names(which.max(table(twilightperiod))),
            twilightperiodB=names(which.max(table(twilightperiodB))),
            latitude.s=mean(latitude.s),
            longitude=mean(longitude),
            #moon.fraction=mean(fraction),
            date=names(which.max(table(date))),
            ordinal.day=mean(julian),
            ordinal.s=mean(ordinal.s),
            ordinal.s2=mean(ordinal.s2),
            time=mean(start_time_numeric))
write.csv(peents.pervisit, file="0_data/processed/6A_PeentActivityRates_MixedModels/temp/peentspervisit.csv")
#confirms that after summarizing I have 20 sample observations per site
#at this point, I can run mixed-effects models or GAMMs


#basic twilight period mixed model for 10-minute peent intervals

m.ordtwilightB <- glmer.nb(numdetect ~ twilightperiodB +(1|sm_id), data=peents.pervisit, verbose=TRUE)
#AIC=1208.2 

#TWILIGHT PERIOD-ONLY MODEL

#adding latitude makes model less likely to converge but latitude may no longer be as important
#after excluding sites lacking peent detections, latitudinal range is just over half (45-54) what it was (45-63)
#Now create the bootstrap
bs <- function(data, indices){
  d<-fSampleNint(data, intervals, 12)
  
  peents.pervisit<-d %>%
    group_by(sm_id, intervalsP) %>% 
    summarize(numdetect= sum(PEENTS), 
              duration=mean(duration_new),
              TSSS=mean(TSSScorr), 
              TSSS.s=mean(TSSS.s), 
              meantemp=mean(meantemp),
              meantemp.s=mean(meantemp.s),
              twilightperiod=names(which.max(table(twilightperiod))),
              twilightperiodB=names(which.max(table(twilightperiodB))),
              latitude=mean(latitude),
              latitude.s=mean(latitude.s),
              longitude=mean(longitude),
              #moon.fraction=mean(fraction),
              date=names(which.max(table(date))),
              ordinal.day=mean(julian),
              ordinal.s=mean(ordinal.s),
              ordinal.s2=mean(ordinal.s2),
              time=mean(start_time_numeric))
  
  #basic model for 10-minute peent intervals
  m.ordtwilight <- try(glmer.nb(numdetect ~ twilightperiodB  +(1|sm_id), data=peents.pervisit, verbose=TRUE))
  if (!inherits(m.ordtwilight, "try-error")){
    coef.mem<-try(fixef(m.ordtwilight))
    #coef.mem<-try(data.frame(summary(m.ordtwilight)["coefficients"])$Estimate)
    coef.df<-try(data.frame(coef.mem))
    try(return(coef.df[,1]))#problems trying to return a single column of values
  }
  else {
    coef.df<-data.frame(rep(NA,9))#same as number of model coefficients
    try(return(coef.df[,1]))
  }
}

#peent.all<-read.csv("0_data/processed/3_PeentsMapped_SunAndMoon/TempSunMoonPeentDetections.10min.int.csv", header=TRUE)

peentresults2 <- boot(data=peentsintD.train,
                    statistic=bs, parallel="multicore",
                    R=100)
#print(paste0("bootstraps successfully run for ",fn))


CoefB <-t(peentresults2$t)

CoefB<-data.frame(CoefB)
CoefB[] <- lapply(CoefB, function(x) as.numeric(as.character(x)))
#delete columns with NA values
all_na <- function(x) any(!is.na(x))
CoefB<-CoefB %>% select_if(all_na)
CoefB<-as.matrix(CoefB)
ncol(CoefB)#38
CoefBCI50 <- t(apply(CoefB, 1, quantile, c(0.5, 0.25, 0.75), na.rm=TRUE))
CoefBCI90 <- t(apply(CoefB, 1, quantile, c(0.5, 0.05, 0.95), na.rm=TRUE))
save(peentresults2, CoefB, CoefBCI50, CoefBCI90, file="3_output/data/6A_PeentMixedEffectsTwilightPeriodModels/boot12samples10minuteintPeents_Twilight.RData")
#save as csv
df<-data.frame(CoefB)
row.names(df) = c("(Intercept)",
                  "CIVIL.1",
                  "NAUTICAL.1",
                  "ASTRONOMICAL.1",
                  "NIGHT",
                  "ASTRONOMICAL.2",
                  "NAUTICAL.2",
                  "CIVIL.2",
                  "AFTER")
write.csv(df, file="3_output/data/6A_PeentMixedEffectsTwilightPeriodModels/boot12samples10minuteintPeents_Twilight.csv")
#print(paste0("bootstrap coefficients saved for ",fn))

## box plot of variables with confidence intervals
prednames = c("(Intercept)",
              "CIVIL.1",
              "NAUTICAL.1",
              "ASTRONOMICAL.1",
              "NIGHT",
              "ASTRONOMICAL.2",
              "NAUTICAL.2",
              "CIVIL.2",
              "AFTER")
var.bci<-CoefBCI90
var.bci<-data.frame(var.bci)
var.bci$median<-as.numeric(var.bci$X50.)
var.bci$lcl<-as.numeric(var.bci$X5.)
var.bci$ucl<-as.numeric(var.bci$X95.)
var.bci$prednames<-prednames

GG<-ggplot(var.bci, aes(x=prednames, y=median))+
  geom_point(aes(x=prednames, y=median))+
  geom_errorbar(aes(ymin=lcl,ymax=ucl))+
  geom_hline(yintercept=0)+
  xlab("Predictor")+
  ylab("Effect on peents counted")+coord_flip()+my.theme
ggsave("3_output/figures/Twilight Period Mixed Models/boot12visits10minuteintervalPeent_NoLat.jpeg", plot=GG, width=12, height=6, units=c("in"), dpi=300)


#Validating the mixed effects models, using the test data from earlier
#i.e. boomsintD.test, which still must be summarized by interval
peents.pervisit.test<-peentsintD.test %>%
  group_by(sm_id, intervalsP) %>% 
  summarize(numdetect= sum(PEENTS), 
            duration=mean(duration_new),
            TSSS=mean(TSSScorr), 
            TSSS.s=mean(TSSS.s), 
            meantemp=mean(meantemp),
            meantemp.s=mean(meantemp.s),
            twilightperiod=names(which.max(table(twilightperiod))),
            twilightperiodB=names(which.max(table(twilightperiodB))),
            latitude=mean(latitude),
            latitude.s=mean(latitude.s),
            longitude=mean(longitude),
            #moon.fraction=mean(fraction),
            date=names(which.max(table(date))),
            ordinal.day=mean(julian),
            ordinal.s=mean(ordinal.s),
            ordinal.s2=mean(ordinal.s2),
            time=mean(start_time_numeric))

#model.matrix function doesn't work with glmer.nb?
peents.pervisit.test$Intercept<-1
peents.pervisit.test$twilightperiodBCIVIL.1<-ifelse(peents.pervisit.test$twilightperiodB=="CIVIL.1",1,0)
peents.pervisit.test$twilightperiodBNAUTICAL.1<-ifelse(peents.pervisit.test$twilightperiodB=="NAUTICAL.1",1,0)
peents.pervisit.test$twilightperiodBASTRONOMICAL.1<-ifelse(peents.pervisit.test$twilightperiodB=="ASTRONOMICAL.1",1,0)
peents.pervisit.test$twilightperiodBNIGHT<-ifelse(peents.pervisit.test$twilightperiodB=="NIGHT",1,0)
peents.pervisit.test$twilightperiodBASTRONOMICAL.2<-ifelse(peents.pervisit.test$twilightperiodB=="ASTRONOMICAL.2",1,0)
peents.pervisit.test$twilightperiodBNAUTICAL.2<-ifelse(peents.pervisit.test$twilightperiodB=="NAUTICAL.2",1,0)
peents.pervisit.test$twilightperiodBCIVIL.2<-ifelse(peents.pervisit.test$twilightperiodB=="CIVIL.2",1,0)
peents.pervisit.test$twilightperiodBAFTER<-ifelse(peents.pervisit.test$twilightperiodB=="AFTER",1,0)

Xnew.Test.mm<-as.matrix(peents.pervisit.test[,c("Intercept","twilightperiodBCIVIL.1","twilightperiodBNAUTICAL.1","twilightperiodBASTRONOMICAL.1","twilightperiodBNIGHT","twilightperiodBASTRONOMICAL.2","twilightperiodBNAUTICAL.2","twilightperiodBCIVIL.2","twilightperiodBAFTER")])

## predict and apply inverse link function: this gives an n_new x B+1 matrix
df<-read.csv("3_output/data/6A_PeentMixedEffectsTwilightPeriodModels/boot12samples10minuteintPeents_Twilight.csv",header=TRUE)
CoefB<-as.matrix(df[,c("X7","X9","X12","X19","X21","X22","X24",
                       "X25","X31","X33","X37","X38","X44","X46",
                       "X49","X51","X60","X65","X68","X70","X71",
                       "X72","X73","X74","X80","X84","X85","X87",
                       "X88","X89","X92","X93","X94","X96","X97",
                       "X98","X99","X100")])
Preds <- (exp(Xnew.Test.mm %*% CoefB))

peents.pervisit.test.df<-data.frame(peents.pervisit.test)
Preds.df<-data.frame(Preds)
TestDatapreds<-cbind(peents.pervisit.test.df, Preds.df)

str(TestDatapreds)
write.csv(TestDatapreds, file="3_output/data/6A_PeentMixedEffectsTwilightPeriodModels/TestDatapredsNoLat.csv")

spearmancorlist<-list()
names<-c("X7","X9","X12","X19","X21","X22","X24",
         "X25","X31","X33","X37","X38","X44","X46",
         "X49","X51","X60","X65","X68","X70","X71",
         "X72","X73","X74","X80","X84","X85","X87",
         "X88","X89","X92","X93","X94","X96","X97",
         "X98","X99","X100")
for (i in names){
  TestDatapreds$Predicted<-TestDatapreds[,i]
  spearmancorlist[[i]]<-cor(TestDatapreds$numdetect, TestDatapreds$Predicted, method="spearman")
}
spearmancorvec<-unlist(spearmancorlist)
spearmandf<-data.frame(spearmancorvec)
## calculate median and 90% CIs
RhoN <-quantile(spearmancorvec, c(0.5, 0.05, 0.95), na.rm=TRUE)
#    50%         5%        95% 
#-0.06776611 -0.10461987 -0.04655391 

#basically no correlation at all

#Graph predicted boom counts by site and twilight period
TestDatapreds<-read.csv("3_output/data/6A_PeentMixedEffectsTwilightPeriodModels/TestDatapredsNolat.csv", header=TRUE)
levels(TestDatapreds$twilightperiodB)[levels(TestDatapreds$twilightperiodB)=="ASTRONOMICAL.1"] <- "ASTRO.1"
levels(TestDatapreds$twilightperiodB)[levels(TestDatapreds$twilightperiodB)=="ASTRONOMICAL.2"] <- "ASTRO.2"

TestDatapreds.g<-TestDatapreds%>%
  gather(X7:X100, key="bootstrap", value="predictedcount")

TestDatapreds.g$twilightperiodB<-as.factor(TestDatapreds.g$twilightperiodB)

TestDatapreds.g$twilightperiodB<-factor(TestDatapreds.g$twilightperiodB, levels=c("BEFORE","CIVIL.1","NAUTICAL.1","ASTRO.1","NIGHT","ASTRO.2","NAUTICAL.2","CIVIL.2","AFTER"))
levels(TestDatapreds.g$twilightperiodB)

TestDatapreds.g$predictedcount.r<-round(TestDatapreds.g$predictedcount, digits=1)

# Basic violin plot
p <- ggplot(TestDatapreds.g, aes(x=twilightperiodB, y=predictedcount.r)) + 
  geom_violin()+
  my.theme+
  ylab("Predicted # peents per 10 min")+xlab("Twilight Period")

#Now loop through sites
sitenames<-levels(as.factor(TestDatapreds.g$sm_id))
for (i in sitenames){
  site<-TestDatapreds.g[TestDatapreds.g$sm_id==i,]
  GG <- ggplot(site, aes(x=twilightperiodB, y=predictedcount.r)) + 
    geom_violin()+
    my.theme+
    ylab("Predicted # peents per 10 min")+
    xlab("Twilight Period")+
    ggtitle(paste0("Site ",i))
  
  ggsave(paste0("3_output/figures/Twilight Period Mixed Models/Violin Plots Peent Twilight Only Model/nolatpredictedpeentcountsXtwilightSite",i,".jpeg"), plot=GG, width=12, height=6, units=c("in"), dpi=300)
  print(paste0("plot for ",i," printed"))
}


#TWILIGHT*LATITUDE INTERACTION MODEL
bs2 <- function(data, indices){
  d<-fSampleNint(data, intervals, 12)
  
  peents.pervisit<-d %>%
    group_by(sm_id, intervalsP) %>% 
    summarize(numdetect= sum(PEENTS), 
              duration=mean(duration_new),
              TSSS=mean(TSSScorr), 
              TSSS.s=mean(TSSS.s), 
              meantemp=mean(meantemp),
              meantemp.s=mean(meantemp.s),
              twilightperiod=names(which.max(table(twilightperiod))),
              twilightperiodB=names(which.max(table(twilightperiodB))),
              latitude=mean(latitude),
              latitude.s=mean(latitude.s),
              longitude=mean(longitude),
              #moon.fraction=mean(fraction),
              date=names(which.max(table(date))),
              ordinal.day=mean(julian),
              ordinal.s=mean(ordinal.s),
              ordinal.s2=mean(ordinal.s2),
              time=mean(start_time_numeric))
  
  #basic model for 10-minute peent intervals
  m.ordtwilight <- try(glmer.nb(numdetect ~ twilightperiodB+latitude.s+twilightperiodB*latitude.s  +(1|sm_id), data=peents.pervisit, verbose=TRUE))
  if (!inherits(m.ordtwilight, "try-error")){
    coef.mem<-try(fixef(m.ordtwilight))
    #coef.mem<-try(data.frame(summary(m.ordtwilight)["coefficients"])$Estimate)
    coef.df<-try(data.frame(coef.mem))
    try(return(coef.df[,1]))#problems trying to return a single column of values
  }
  else {
    coef.df<-data.frame(rep(NA,18))#same as number of model coefficients
    try(return(coef.df[,1]))
  }
}

#peent.all<-read.csv("0_data/processed/3_PeentsMapped_SunAndMoon/TempSunMoonPeentDetections.10min.int.csv", header=TRUE)

peentresults3 <- boot(data=peentsintD.train,
                      statistic=bs2, parallel="multicore",
                      R=100)


CoefB2 <-t(peentresults3$t)

CoefB2<-data.frame(CoefB2)
CoefB2[] <- lapply(CoefB2, function(x) as.numeric(as.character(x)))
#delete columns with NA values
all_na <- function(x) any(!is.na(x))
CoefB2<-CoefB2 %>% select_if(all_na)
CoefB2<-as.matrix(CoefB2)
ncol(CoefB2)#38
CoefBCI50 <- t(apply(CoefB2, 1, quantile, c(0.5, 0.25, 0.75), na.rm=TRUE))
CoefBCI90 <- t(apply(CoefB2, 1, quantile, c(0.5, 0.05, 0.95), na.rm=TRUE))
save(peentresults3, CoefB2, CoefBCI50, CoefBCI90, file="3_output/data/6A_PeentMixedEffectsTwilightPeriodModels/boot12samples10minuteintPeents_TwilightLatInt.RData")
#save as csv
df<-data.frame(CoefB2)
row.names(df) = c("(Intercept)",
                  "CIVIL.1",
                  "NAUTICAL.1",
                  "ASTRONOMICAL.1",
                  "NIGHT",
                  "ASTRONOMICAL.2",
                  "NAUTICAL.2",
                  "CIVIL.2",
                  "AFTER",
                  "LAT",
                  "LAT*CIVIL.1",
                  "LAT*NAUTICAL.1",
                  "LAT*ASTRONOMICAL.1",
                  "LAT*NIGHT",
                  "LAT*ASTRONOMICAL.2",
                  "LAT*NAUTICAL.2",
                  "LAT*CIVIL.2",
                  "LAT*AFTER")
write.csv(df, file="3_output/data/6A_PeentMixedEffectsTwilightPeriodModels/boot12samples10minuteintPeents_TwilightLatInt.csv")
#print(paste0("bootstrap coefficients saved for ",fn))

## box plot of variables with confidence intervals
prednames = c("(Intercept)",
              "CIVIL.1",
              "NAUTICAL.1",
              "ASTRONOMICAL.1",
              "NIGHT",
              "ASTRONOMICAL.2",
              "NAUTICAL.2",
              "CIVIL.2",
              "AFTER",
              "LAT",
              "LAT*CIVIL.1",
              "LAT*NAUTICAL.1",
              "LAT*ASTRONOMICAL.1",
              "LAT*NIGHT",
              "LAT*ASTRONOMICAL.2",
              "LAT*NAUTICAL.2",
              "LAT*CIVIL.2",
              "LAT*AFTER")
var.bci<-CoefBCI90
var.bci<-data.frame(var.bci)
var.bci$median<-as.numeric(var.bci$X50.)
var.bci$lcl<-as.numeric(var.bci$X5.)
var.bci$ucl<-as.numeric(var.bci$X95.)
var.bci$prednames<-prednames

GG<-ggplot(var.bci, aes(x=prednames, y=median))+
  geom_point(aes(x=prednames, y=median))+
  geom_errorbar(aes(ymin=lcl,ymax=ucl))+
  geom_hline(yintercept=0)+
  xlab("Predictor")+
  ylab("Effect on peents counted")+coord_flip()+my.theme
ggsave("3_output/figures/Twilight Period Mixed Models/boot12visits10minuteintervalPeent_LatInt.jpeg", plot=GG, width=12, height=6, units=c("in"), dpi=300)
