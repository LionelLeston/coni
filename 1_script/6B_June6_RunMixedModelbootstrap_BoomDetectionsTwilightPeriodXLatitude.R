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
        axis.text.x=element_text(size=8),
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
}#sampled by twilightperiod rather than twilightperiodB to get more samples for test data

#You need to run the mixed model first for later, when generating a
#model matrix for the prediction data sets
boomsint<-read.csv("0_data/processed/3_BoomsMapped_SunAndMoon/TempSunMoonBoomDetections.10min.int.csv", header=TRUE)
str(boomsint)
boomsint$BOOMS<-ifelse(is.na(boomsint$eventID),0,1)
#boomsint$halfmonth<-ifelse(boomsint$day<16,"1","2")
#boomsint$season<-as.factor(paste0(boomsint$month,"_",boomsint$halfmonth))

boomsint$start_time_new<- ymd_hms(as.character(boomsint$start_time_new))
boomsint$julian<-yday(as.Date(boomsint$date))
boomsint$hour<-hour(boomsint$start_time_new)
boomsint$minute<-minute(boomsint$start_time_new)
boomsint$start_time_numeric<-boomsint$hour+(boomsint$minute/60)

boomsint$TSSS.s<-scale(boomsint$TSSS, center=TRUE, scale=TRUE)#scaled to reduce rank deficiency of GAMMs
boomsint$ordinal.s<-scale(boomsint$julian, center=TRUE, scale=TRUE)#scaled to reduce rank deficiency of GAMMs
boomsint$TSSS.s2<-scale(boomsint$TSSS^2, center=FALSE, scale=TRUE)#scaled to reduce rank deficiency of GAMMs
boomsint$ordinal.s2<-scale(boomsint$julian^2, center=FALSE, scale=TRUE)#scaled to reduce rank deficiency of GAMMs
boomsint$latitude.s<-scale(boomsint$latitude, center=FALSE, scale=TRUE)#scaled to reduce rank deficiency of GAMMs
boomsint$meantemp.s<-scale(boomsint$meantemp, center=FALSE, scale=TRUE)#not currently being used in these GAMMs
#boomsint$moon.s<-scale(as.vector(boomsint$moon.fraction), center=TRUE, scale=FALSE)#not currently being used in these GAMMs

boomsint$sm_id<-as.factor(boomsint$sm_id)
boomsint$intervals<-as.factor(boomsint$intervals)

boomsintB<-boomsint[which(boomsint$julian<190&boomsint$julian>170),]

#boomsintB<-boomsint[!boomsint$month==8,]#drop August


#boomsintC<-boomsintB[!is.na(boomsintB$meantemp),]#drop intervals missing temperature boomsint

boomsintC<-boomsintB[!boomsintB$sm_id %in% c("7524","7698","8245","13523","13585","13588"),]
#remove sites that had no "Night" twilight period because they were too far north
nrow(boomsintC)#8877
maxdur<-max(boomsintC$duration)
boomsintD<-boomsintC[boomsintC$duration==maxdur,]#gets rid of any "remainder" intervals
nrow(boomsintD)#8799

boomsintD$twilightperiodB<-boomsintD$twilightperiod
boomsintD$twilightperiodB<-ifelse(boomsintD$start_time_numeric>12,
                                   paste0(boomsintD$twilightperiodB,".1"),
                                   paste0(boomsintD$twilightperiodB,".2"))
boomsintD$twilightperiodB<-as.factor(boomsintD$twilightperiodB)
levels(boomsintD$twilightperiodB)

levels(boomsintD$twilightperiodB)[levels(boomsintD$twilightperiodB)=="AFTER.2"] <- "AFTER"
levels(boomsintD$twilightperiodB)[levels(boomsintD$twilightperiodB)=="BEFORE.1"] <- "BEFORE"
levels(boomsintD$twilightperiodB)[levels(boomsintD$twilightperiodB)=="NIGHT.2"] <- "NIGHT"


boomsintD$twilightperiodB<-factor(boomsintD$twilightperiodB, levels=c("BEFORE","CIVIL.1","NAUTICAL.1","ASTRONOMICAL.1","NIGHT","ASTRONOMICAL.2","NAUTICAL.2","CIVIL.2","AFTER"))
levels(boomsintD$twilightperiodB)

boomsintD$sm_id.interval<-as.factor(paste0(boomsintD$sm_id,boomsintD$intervals))
boomsintD.test <- fSampleNint(boomsintD, intervals, 10)
nrow(boomsintD.test)#
boomsintD.train  <- boomsintD[!boomsintD$sm_id.interval %in% boomsintD.test$sm_id.interval,]
nrow(boomsintD.train)#

levels(boomsintD.train$twilightperiodB)


#Mixed-effects model run once before bootstrap to get model matrix
#Run from "d<-fSampleNint..." to "m.ordtwilightXLatB<-glmer.nb..." until model converges
d<-fSampleNint(boomsintD.train, intervals, 12)#13 or more intervals is too many sample without replacement

booms.pervisit<-d %>%
  group_by(sm_id, intervalsP) %>% 
  summarize(numdetect= sum(BOOMS), 
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
write.csv(booms.pervisit, file="0_data/processed/9_boomActivityRates_MixedModels/temp/boomspervisit.csv")
#confirms that after summarizing I have 20 sample observations per site
#at this point, I can run mixed-effects models or GAMMs

booms.pervisit$twilightperiodB<-as.factor(booms.pervisit$twilightperiodB)
levels(booms.pervisit$twilightperiodB)
booms.pervisit$twilightperiodB<-factor(booms.pervisit$twilightperiodB, levels=c("BEFORE","CIVIL.1","NAUTICAL.1","ASTRONOMICAL.1","NIGHT","ASTRONOMICAL.2","NAUTICAL.2","CIVIL.2","AFTER"))
levels(booms.pervisit$twilightperiodB)

#basic twilight period mixed model for 10-minute boom intervals

#m.ordtwilightB <- glmer.nb(numdetect ~ twilightperiodB +(1|sm_id), data=booms.pervisit, verbose=TRUE)
#AIC=1461.1 
m.ordtwilightXLatB <- glmer.nb(numdetect ~ twilightperiodB+latitude.s+twilightperiodB*latitude.s +(1|sm_id), data=booms.pervisit, verbose=TRUE)
#AIC=1448
 
#Now create the bootstrap
bs <- function(data, indices){
  d<-fSampleNint(data, intervals, 12)
  
  booms.pervisit<-d %>%
    group_by(sm_id, intervalsP) %>% 
    summarize(numdetect= sum(BOOMS), 
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
  
  booms.pervisit$twilightperiodB<-as.factor(booms.pervisit$twilightperiodB)
  levels(booms.pervisit$twilightperiodB)
  booms.pervisit$twilightperiodB<-factor(booms.pervisit$twilightperiodB, levels=c("BEFORE","CIVIL.1","NAUTICAL.1","ASTRONOMICAL.1","NIGHT","ASTRONOMICAL.2","NAUTICAL.2","CIVIL.2","AFTER"))
  levels(booms.pervisit$twilightperiodB)
  
  #basic model for 10-minute boom intervals
  m.ordtwilight <- try(glmer.nb(numdetect ~ twilightperiodB  +(1|sm_id), data=booms.pervisit, verbose=TRUE))
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

#boom.all<-read.csv("0_data/processed/3_boomsMapped_SunAndMoon/TempSunMoonboomDetections.10min.int.csv", header=TRUE)

boomresults2 <- boot(data=boomsintD.train,
                    statistic=bs, parallel="multicore",
                    R=100)
#print(paste0("bootstraps successfully run for ",fn))


CoefB <-t(boomresults2$t)

CoefB<-data.frame(CoefB)
CoefB[] <- lapply(CoefB, function(x) as.numeric(as.character(x)))
#delete columns with NA values
all_na <- function(x) any(!is.na(x))
CoefB<-CoefB %>% select_if(all_na)
CoefB<-as.matrix(CoefB)
ncol(CoefB)#38
CoefBCI50 <- t(apply(CoefB, 1, quantile, c(0.5, 0.25, 0.75), na.rm=TRUE))
CoefBCI90 <- t(apply(CoefB, 1, quantile, c(0.5, 0.05, 0.95), na.rm=TRUE))
save(boomresults2, CoefB, CoefBCI50, CoefBCI90, file="3_output/data/6B_BoomMixedEffectsTwilightPeriodModels/boot12samples10minuteintbooms_TwilightJune20.RData")
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
write.csv(df, file="3_output/data/6B_BoomMixedEffectsTwilightPeriodModels/boot12samples10minuteintbooms_TwilightJune21.csv")
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
  ylab("Effect on booms counted")+coord_flip()+my.theme
ggsave("3_output/figures/Twilight Period Mixed Models/boot12visits10minuteintervalboom_NoLatJune20.jpeg", plot=GG, width=12, height=6, units=c("in"), dpi=300)

#Validating the mixed effects models, using the test data from earlier
#i.e. boomsintD.test, which still must be summarized by interval
booms.pervisit.test<-boomsintD.test %>%
  group_by(sm_id, intervalsP) %>% 
  summarize(numdetect= sum(BOOMS), 
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
booms.pervisit.test$Intercept<-1
booms.pervisit.test$twilightperiodBCIVIL.1<-ifelse(booms.pervisit.test$twilightperiodB=="CIVIL.1",1,0)
booms.pervisit.test$twilightperiodBNAUTICAL.1<-ifelse(booms.pervisit.test$twilightperiodB=="NAUTICAL.1",1,0)
booms.pervisit.test$twilightperiodBASTRONOMICAL.1<-ifelse(booms.pervisit.test$twilightperiodB=="ASTRONOMICAL.1",1,0)
booms.pervisit.test$twilightperiodBNIGHT<-ifelse(booms.pervisit.test$twilightperiodB=="NIGHT",1,0)
booms.pervisit.test$twilightperiodBASTRONOMICAL.2<-ifelse(booms.pervisit.test$twilightperiodB=="ASTRONOMICAL.2",1,0)
booms.pervisit.test$twilightperiodBNAUTICAL.2<-ifelse(booms.pervisit.test$twilightperiodB=="NAUTICAL.2",1,0)
booms.pervisit.test$twilightperiodBCIVIL.2<-ifelse(booms.pervisit.test$twilightperiodB=="CIVIL.2",1,0)
booms.pervisit.test$twilightperiodBAFTER<-ifelse(booms.pervisit.test$twilightperiodB=="AFTER",1,0)

Xnew.Test.mm<-as.matrix(booms.pervisit.test[,c("Intercept","twilightperiodBCIVIL.1","twilightperiodBNAUTICAL.1","twilightperiodBASTRONOMICAL.1","twilightperiodBNIGHT","twilightperiodBASTRONOMICAL.2","twilightperiodBNAUTICAL.2","twilightperiodBCIVIL.2","twilightperiodBAFTER")])
## predict and apply inverse link function: this gives an n_new x B+1 matrix
Preds <- (exp(Xnew.Test.mm %*% CoefB))

booms.pervisit.test.df<-data.frame(booms.pervisit.test)
Preds.df<-data.frame(Preds)
TestDatapreds<-cbind(booms.pervisit.test.df, Preds.df)

str(TestDatapreds)
write.csv(TestDatapreds, file="3_output/data/6B_BoomMixedEffectsTwilightPeriodModels/TestDatapredsNoLat.csv")

spearmancorlist<-list()
names<-c("X1","X2","X3","X4","X5","X6","X7","X8","X9","X10",
         "X11","X12","X13","X14","X15","X16","X17","X18","X19","X20",
         "X21","X22","X23","X24","X25","X26","X27","X28","X29","X30",
         "X31","X32","X33","X34","X35","X36","X37","X38","X39","X40",
         "X41","X42","X43","X44","X45","X46","X47","X48","X49","X50",
         "X51","X52","X53","X54","X55","X56","X57","X58","X59","X60",
         "X61","X62","X63","X64","X65","X66","X67","X68","X69","X70",
         "X71","X72","X73","X74","X75","X76","X77","X78","X79","X80",
         "X81","X82","X83","X84","X85","X86","X87","X88","X89","X90",
         "X91","X92","X93","X94","X95","X96","X97","X98","X99","X100")
for (i in names){
  TestDatapreds$Predicted<-TestDatapreds[,i]
  spearmancorlist[[i]]<-cor(TestDatapreds$numdetect, TestDatapreds$Predicted, method="spearman")
}
spearmancorvec<-unlist(spearmancorlist)
spearmandf<-data.frame(spearmancorvec)
## calculate median and 90% CIs
RhoN <-quantile(spearmancorvec, c(0.5, 0.05, 0.95), na.rm=TRUE)
#    50%         5%        95% 
#0.2527481 0.2413596 0.2607584 

#not great but arguably better predictions than GAMM models

#Graph predicted boom counts by site and twilight period
TestDatapreds<-read.csv("3_output/data/6B_BoomMixedEffectsTwilightPeriodModels/TestDatapredsNolat.csv", header=TRUE)
levels(TestDatapreds$twilightperiodB)[levels(TestDatapreds$twilightperiodB)=="ASTRONOMICAL.1"] <- "ASTRO.1"
levels(TestDatapreds$twilightperiodB)[levels(TestDatapreds$twilightperiodB)=="ASTRONOMICAL.2"] <- "ASTRO.2"

TestDatapreds.g<-TestDatapreds%>%
  gather(X1:X100, key="bootstrap", value="predictedcount")

TestDatapreds.g$twilightperiodB<-as.factor(TestDatapreds.g$twilightperiodB)

TestDatapreds.g$twilightperiodB<-factor(TestDatapreds.g$twilightperiodB, levels=c("BEFORE","CIVIL.1","NAUTICAL.1","ASTRO.1","NIGHT","ASTRO.2","NAUTICAL.2","CIVIL.2","AFTER"))
levels(TestDatapreds.g$twilightperiodB)

TestDatapreds.g$predictedcount.r<-round(TestDatapreds.g$predictedcount, digits=1)

# Basic violin plot
p <- ggplot(TestDatapreds.g, aes(x=twilightperiodB, y=predictedcount.r)) + 
  geom_violin()+
  my.theme+
  ylab("Predicted # booms per 10 min")+xlab("Twilight Period")

#Now loop through sites
sitenames<-levels(as.factor(TestDatapreds.g$sm_id))
for (i in sitenames){
  site<-TestDatapreds.g[TestDatapreds.g$sm_id==i,]
  GG <- ggplot(site, aes(x=twilightperiodB, y=predictedcount.r)) + 
    geom_violin()+
    my.theme+
    ylab("Predicted # booms per 10 min")+
    xlab("Twilight Period")+
    ggtitle(paste0("Site ",i))
  
  ggsave(paste0("3_output/figures/Twilight Period Mixed Models/nolatpredictedboomcountsXtwilightSite",i,".jpeg"), plot=GG, width=12, height=6, units=c("in"), dpi=300)
  print(paste0("plot for ",i," printed"))
}





#TWILIGHT*LATITUDE INTERACTION MODEL
bs2 <- function(data, indices){
  d<-fSampleNint(data, intervals, 12)
  
  booms.pervisit<-d %>%
    group_by(sm_id, intervalsP) %>% 
    summarize(numdetect= sum(BOOMS), 
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
  
  booms.pervisit$twilightperiodB<-as.factor(booms.pervisit$twilightperiodB)
  levels(booms.pervisit$twilightperiodB)
  booms.pervisit$twilightperiodB<-factor(booms.pervisit$twilightperiodB, levels=c("BEFORE","CIVIL.1","NAUTICAL.1","ASTRONOMICAL.1","NIGHT","ASTRONOMICAL.2","NAUTICAL.2","CIVIL.2","AFTER"))
  levels(booms.pervisit$twilightperiodB)
  
  #basic model for 10-minute boom intervals
  m.ordtwilight <- try(glmer.nb(numdetect ~ twilightperiodB+latitude.s+twilightperiodB*latitude.s  +(1|sm_id), data=booms.pervisit, verbose=TRUE))
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

#boom.all<-read.csv("0_data/processed/3_boomsMapped_SunAndMoon/TempSunMoonboomDetections.10min.int.csv", header=TRUE)

boomresults3 <- boot(data=boomsintD.train,
                      statistic=bs2, parallel="multicore",
                      R=100)


CoefB2 <-t(boomresults3$t)

CoefB2<-data.frame(CoefB2)
CoefB2[] <- lapply(CoefB2, function(x) as.numeric(as.character(x)))
#delete columns with NA values
all_na <- function(x) any(!is.na(x))
CoefB2<-CoefB2 %>% select_if(all_na)
CoefB2<-as.matrix(CoefB2)
ncol(CoefB2)#27
CoefBCI50 <- t(apply(CoefB2, 1, quantile, c(0.5, 0.25, 0.75), na.rm=TRUE))
CoefBCI90 <- t(apply(CoefB2, 1, quantile, c(0.5, 0.05, 0.95), na.rm=TRUE))
save(boomresults3, CoefB2, CoefBCI50, CoefBCI90, file="3_output/data/6B_BoomMixedEffectsTwilightPeriodModels/boot12samples10minuteintbooms_TwilightLatIntJune20.RData")
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
write.csv(df, file="3_output/data/6B_BoomMixedEffectsTwilightPeriodModels/boot12samples10minuteintbooms_TwilightLatIntJune20.csv")
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
  ylab("Effect on booms counted")+coord_flip()+my.theme
ggsave("3_output/figures/Twilight Period Mixed Models/boot12visits10minuteintervalboom_LatIntJune20.jpeg", plot=GG, width=12, height=6, units=c("in"), dpi=300)


#Validating the mixed effects models, using the test data from earlier
#i.e. boomsintD.test, which still must be summarized by interval
booms.pervisit.test<-boomsintD.test %>%
  group_by(sm_id, intervalsP) %>% 
  summarize(numdetect= sum(BOOMS), 
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

## model matrix that'll match the coefficients
Xnew.Test.mm <- model.matrix(m.ordtwilightXLatB, booms.pervisit.test)
#model.matrix function doesn't work with glmer.nb?
booms.pervisit.test$Intercept<-1
booms.pervisit.test$twilightperiodBCIVIL.1<-ifelse(booms.pervisit.test$twilightperiodB=="CIVIL.1",1,0)
booms.pervisit.test$twilightperiodBNAUTICAL.1<-ifelse(booms.pervisit.test$twilightperiodB=="NAUTICAL.1",1,0)
booms.pervisit.test$twilightperiodBASTRONOMICAL.1<-ifelse(booms.pervisit.test$twilightperiodB=="ASTRONOMICAL.1",1,0)
booms.pervisit.test$twilightperiodBNIGHT<-ifelse(booms.pervisit.test$twilightperiodB=="NIGHT",1,0)
booms.pervisit.test$twilightperiodBASTRONOMICAL.2<-ifelse(booms.pervisit.test$twilightperiodB=="ASTRONOMICAL.2",1,0)
booms.pervisit.test$twilightperiodBNAUTICAL.2<-ifelse(booms.pervisit.test$twilightperiodB=="NAUTICAL.2",1,0)
booms.pervisit.test$twilightperiodBCIVIL.2<-ifelse(booms.pervisit.test$twilightperiodB=="CIVIL.2",1,0)
booms.pervisit.test$twilightperiodBAFTER<-ifelse(booms.pervisit.test$twilightperiodB=="AFTER",1,0)
booms.pervisit.test$twilightperiodBCIVIL.1XLAT<-ifelse(booms.pervisit.test$twilightperiodB=="CIVIL.1",booms.pervisit.test$latitude.s,0)
booms.pervisit.test$twilightperiodBNAUTICAL.1XLAT<-ifelse(booms.pervisit.test$twilightperiodB=="NAUTICAL.1",booms.pervisit.test$latitude.s,0)
booms.pervisit.test$twilightperiodBASTRONOMICAL.1XLAT<-ifelse(booms.pervisit.test$twilightperiodB=="ASTRONOMICAL.1",booms.pervisit.test$latitude.s,0)
booms.pervisit.test$twilightperiodBNIGHTXLAT<-ifelse(booms.pervisit.test$twilightperiodB=="NIGHT",booms.pervisit.test$latitude.s,0)
booms.pervisit.test$twilightperiodBASTRONOMICAL.2XLAT<-ifelse(booms.pervisit.test$twilightperiodB=="ASTRONOMICAL.2",booms.pervisit.test$latitude.s,0)
booms.pervisit.test$twilightperiodBNAUTICAL.2XLAT<-ifelse(booms.pervisit.test$twilightperiodB=="NAUTICAL.2",booms.pervisit.test$latitude.s,0)
booms.pervisit.test$twilightperiodBCIVIL.2XLAT<-ifelse(booms.pervisit.test$twilightperiodB=="CIVIL.2",booms.pervisit.test$latitude.s,0)
booms.pervisit.test$twilightperiodBAFTERXLAT<-ifelse(booms.pervisit.test$twilightperiodB=="AFTER",booms.pervisit.test$latitude.s,0)

Xnew.Test.mm<-as.matrix(booms.pervisit.test[,c("Intercept","twilightperiodBCIVIL.1","twilightperiodBNAUTICAL.1","twilightperiodBASTRONOMICAL.1","twilightperiodBNIGHT","twilightperiodBASTRONOMICAL.2","twilightperiodBNAUTICAL.2","twilightperiodBCIVIL.2","twilightperiodBAFTER","latitude.s","twilightperiodBCIVIL.1XLAT","twilightperiodBNAUTICAL.1XLAT","twilightperiodBASTRONOMICAL.1XLAT","twilightperiodBNIGHTXLAT","twilightperiodBASTRONOMICAL.2XLAT","twilightperiodBNAUTICAL.2XLAT","twilightperiodBCIVIL.2XLAT","twilightperiodBAFTERXLAT")])
## predict and apply inverse link function: this gives an n_new x B+1 matrix
Preds <- (exp(Xnew.Test.mm %*% CoefB2))

booms.pervisit.test.df<-data.frame(booms.pervisit.test)
Preds.df<-data.frame(Preds)
TestDatapreds<-cbind(booms.pervisit.test.df, Preds.df)

str(TestDatapreds)
write.csv(TestDatapreds, file="3_output/data/6B_BoomMixedEffectsTwilightPeriodModels/TestDatapreds.csv")

spearmancorlist<-list()
names<-c("X5","X9","X12","X22","X24","X26","X27","X32","X43","X48","X49","X56","X61","X62","X66","X71","X73","X74","X78","X80","X82","X84","X85","X92","X97","X99","X100")
for (i in names){
  TestDatapreds$Predicted<-TestDatapreds[,i]
  spearmancorlist[[i]]<-cor(TestDatapreds$numdetect, TestDatapreds$Predicted, method="spearman")
}
spearmancorvec<-unlist(spearmancorlist)
spearmandf<-data.frame(spearmancorvec)
## calculate median and 90% CIs
RhoN <-quantile(spearmancorvec, c(0.5, 0.05, 0.95), na.rm=TRUE)
#    50%         5%        95% 
#0.2893885 0.2483996 0.3076131 

#not great but arguably better predictions than GAMM models

#Graph predicted boom counts by site and twilight period
TestDatapreds<-read.csv("3_output/data/6B_BoomMixedEffectsTwilightPeriodModels/TestDatapreds.csv", header=TRUE)
levels(TestDatapreds$twilightperiodB)[levels(TestDatapreds$twilightperiodB)=="ASTRONOMICAL.1"] <- "ASTRO.1"
levels(TestDatapreds$twilightperiodB)[levels(TestDatapreds$twilightperiodB)=="ASTRONOMICAL.2"] <- "ASTRO.2"

TestDatapreds.g<-TestDatapreds%>%
  gather(X5:X100, key="bootstrap", value="predictedcount")

TestDatapreds.g$twilightperiodB<-as.factor(TestDatapreds.g$twilightperiodB)

TestDatapreds.g$twilightperiodB<-factor(TestDatapreds.g$twilightperiodB, levels=c("BEFORE","CIVIL.1","NAUTICAL.1","ASTRO.1","NIGHT","ASTRO.2","NAUTICAL.2","CIVIL.2","AFTER"))
levels(TestDatapreds.g$twilightperiodB)

TestDatapreds.g$predictedcount.r<-round(TestDatapreds.g$predictedcount, digits=1)

# Basic violin plot
p <- ggplot(TestDatapreds.g, aes(x=twilightperiodB, y=predictedcount.r)) + 
  geom_violin()+
  my.theme+
  ylab("Predicted # booms per 10 min")+xlab("Twilight Period")

#Now loop through sites
sitenames<-levels(as.factor(TestDatapreds.g$sm_id))
for (i in sitenames){
  site<-TestDatapreds.g[TestDatapreds.g$sm_id==i,]
  GG <- ggplot(site, aes(x=twilightperiodB, y=predictedcount.r)) + 
    geom_violin()+
    my.theme+
    ylab("Predicted # booms per 10 min")+
    xlab("Twilight Period")+
    ggtitle(paste0("Site ",i))
  
  ggsave(paste0("3_output/figures/Twilight Period Mixed Models/predictedboomcountsXtwilightSite",i,".jpeg"), plot=GG, width=12, height=6, units=c("in"), dpi=300)
  print(paste0("plot for ",i," printed"))
}

