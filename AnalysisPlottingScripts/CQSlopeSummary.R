#===================================================================================#
# NOTES: This script creates a figure visualizing the concentration discharge slopes
# for the different watersheds during stormflow and baseflow and during different
# seasons
#-----------------------------------------------------------------------------------#
# Galen Gorski                                                                      #
# galengorski@berkeley.edu                                                          #
# 2021-04-20                                                                        #  
#-----------------------------------------------------------------------------------#
#===================================================================================#

#===================================================================================#
#####INSTALL PACKAGES#####
library(tidyverse)
#install.packages('magrittr')
library(magrittr)
#install.packages('Hmisc')
library(Hmisc)
#install.packages('ggpubr')
library(ggpubr)
#install.packages('vioplot')
library(vioplot)
#install.packages('data.table')
library(data.table)

#####
#===================================================================================#

#===================================================================================#
#####READ IN THE DATA#####

#site records list
#site.records <- readRDS('L2/siterecords.rds')

#site records list without low flow period in jefferson
site.records <- readRDS('L2/siterecordsnolowflow.rds')

#events as a list of lists
site.events.ls <- readRDS('L2/siteeventslist.rds')

#events as a list of dataframes
site.events.df <- readRDS('L2/siteeventsdataframe.rds')

#non events a as list of dataframes
#site.non.events <- readRDS('L2/sitenonevents.rds')

#non events as a list of dataframes without the low flow periods in jefferson
site.non.events <- readRDS('L2/sitenoneventsnolowflow.rds')

wshed.names <- c('UPN','USC','MRF','MJF','DVM')
#####
#===================================================================================#

#===================================================================================#
#####MAKING A SUMMARY TABLE OF THE EVENTS#####
#make a list to store the event summary

summary.df <- data.frame(event.start = NA, event.end = NA, event.length = NA, season = NA, cum.q = NA, mean.c = NA, cum.load = NA, cq.slp = NA, cvc.cvq = NA, r2 = NA, bsig = NA, max.flow.percentile = NA)

max.flow.sites <- lapply(site.records, function(x) max(x$X_00060_00003, na.rm = T)) %>% unlist()

site.events.summary.ls <- list(pan = summary.df, sc = summary.df, rf = summary.df, jeff = summary.df, vm = summary.df)

#for each site
for(i in 1:5){
  for(j in 1:length(site.events.ls[[i]])){
    temp.event <- site.events.ls[[i]][[j]]
    if(length(temp.event) == 1){
      next
    }else{
      if(nrow(temp.event[!is.na(temp.event$X_99133_00003),]) < 3){
        next
      }else{
        #event start
        site.events.summary.ls[[i]][j,'event.start'] <- as.character(temp.event$Date[1])
        #event end
        site.events.summary.ls[[i]][j,'event.end'] <- as.character(temp.event$Date[nrow(temp.event)])
        #event.duration
        event.dur <- difftime(temp.event$Date[nrow(temp.event)],temp.event$Date[1], units = 'days') %>% as.numeric()
        site.events.summary.ls[[i]][j,'event.length'] <- event.dur
        #get the month for the max flow during the event
        month <- temp.event[temp.event$X_00060_00003 == max(temp.event$X_00060_00003),]$Date[1] %>%
          strptime(.,format = '%Y-%m-%d', tz = 'UTC') %>%
          format(.,format = '%m') %>%
          as.numeric()
        #that's what we'll use to classify the season
        if(month %in% c(1,2,3)){
          site.events.summary.ls[[i]][j,'season'] <- '2-JFM'
        }else if(month %in% c(4,5,6)){
          site.events.summary.ls[[i]][j,'season'] <- '3-AMJ'
        }else if(month %in% c(7,8,9)){
          site.events.summary.ls[[i]][j,'season'] <- '4-JAS'
        }else{
          site.events.summary.ls[[i]][j,'season'] <- '1-OND'
        }
        #event cumulative discharge
        site.events.summary.ls[[i]][j,'cum.q'] <- sum(temp.event$X_00060_00003)
        #event cumulative nitrate
        site.events.summary.ls[[i]][j,'mean.c'] <- mean(temp.event$X_99133_00003, na.rm = T)
        #number of days with both no3 and discharge observations
        fulldays <- nrow(temp.event[!is.na(temp.event$X_00060_00003) & !is.na(temp.event$X_99133_00003),])
        #event cumulative load
        site.events.summary.ls[[i]][j,'cum.load'] <- sum(temp.event$X_99133_00003*temp.event$X_00060_00003, na.rm = T)/(fulldays/event.dur)
        #event cQ slope
        site.events.summary.ls[[i]][j,'cq.slp'] <- coef(lm(log(temp.event$X_99133_00003)~log(temp.event$X_00060_00003)))[2]
        #cvc.cvq
        sigma.q <- sd(temp.event$X_00060_00003)
        mu.q <- mean(temp.event$X_00060_00003)
        sigma.c <- sd(temp.event$X_99133_00003)
        mu.c <- mean(temp.event$X_99133_00003)
        site.events.summary.ls[[i]][j,'cvc.cvq'] <- (sigma.c/mu.c)/(sigma.q/mu.q)
        #event r2
        site.events.summary.ls[[i]][j,'r2'] <- summary(lm(log(temp.event$X_99133_00003)~log(temp.event$X_00060_00003)))$r.square
        #event b significance
        site.events.summary.ls[[i]][j,'bsig'] <- summary(lm(log(temp.event$X_99133_00003)~log(temp.event$X_00060_00003)))$coefficients[2,'Pr(>|t|)']
        #max flow percentile
        site.events.summary.ls[[i]][j,'max.flow.percentile'] <- max(temp.event$X_00060_00003, na.rm = T)/max.flow.sites[i]
      }
    }
  }
  site.events.summary.ls[[i]] <- site.events.summary.ls[[i]][!is.na(site.events.summary.ls[[i]]$event.start),]
}
#####
#===================================================================================#

#===================================================================================#
#####CLASSIFYING THE SEASON FOR THE NON-EVENT PERIODS#####

#make a month column
for(i in 1:5){
  site.non.events[[i]]$month <- strptime(site.non.events[[i]]$dateTime,format = '%Y-%m-%d', tz = 'UTC') %>% 
    format(.,format = '%m') %>% 
    as.numeric()
}

#make a year column
for(i in 1:5){
  site.non.events[[i]]$year <- strptime(site.non.events[[i]]$dateTime,format = '%Y-%m-%d', tz = 'UTC') %>% 
    format(.,format = '%Y') %>% 
    as.numeric()
}



non.event.summary.df <- data.frame(jfm.cq = NA, amj.cq = NA, jas.cq = NA, ond.cq = NA, jfm.cvccvq = NA, amj.cvccvq = NA, jas.cvccvq = NA, ond.cvccvq = NA, 
                                   jfm.r2 = NA, amj.r2 = NA, jas.r2 = NA, ond.r2 = NA, jfm.bsig = NA, amj.bsig = NA, jas.bsig = NA, ond.bsig = NA)
site.non.event.summary.ls <- list(pan = non.event.summary.df, sc = non.event.summary.df, rf = non.event.summary.df, jeff = non.event.summary.df, vm = non.event.summary.df)

#calculate cq and cvc.cvq by season for non-event periods
for(i in 1:5){
  #jfm
  jfm <- site.non.events[[i]][site.non.events[[i]]$month %in% c(1,2,3),]
  #cq jfm
  site.non.event.summary.ls[[i]]$jfm.cq <- coef(lm(log(jfm$X_99133_00003)~log(jfm$X_00060_00003)))[2] 
  #cvc.cvq jfm
  sigma.q <- sd(jfm$X_00060_00003, na.rm = T)
  mu.q <- mean(jfm$X_00060_00003, na.rm = T)
  sigma.c <- sd(jfm$X_99133_00003, na.rm = T)
  mu.c <- mean(jfm$X_99133_00003, na.rm = T)
  site.non.event.summary.ls[[i]]$jfm.cvccvq <- (sigma.c/mu.c)/(sigma.q/mu.q)
  #r2 jfm
  site.non.event.summary.ls[[i]]$jfm.r2 <- summary(lm(log(jfm$X_99133_00003)~log(jfm$X_00060_00003)))$r.square
  #bsig jfm
  site.non.event.summary.ls[[i]]$jfm.bsig <- summary(lm(log(jfm$X_99133_00003)~log(jfm$X_00060_00003)))$coefficients[2,'Pr(>|t|)']
  
  #amj
  amj <- site.non.events[[i]][site.non.events[[i]]$month %in% c(4,5,6),]
  #cq amj
  site.non.event.summary.ls[[i]]$amj.cq <- coef(lm(log(amj$X_99133_00003)~log(amj$X_00060_00003)))[2] 
  #cvc.cvq amj
  sigma.q <- sd(amj$X_00060_00003, na.rm = T)
  mu.q <- mean(amj$X_00060_00003, na.rm = T)
  sigma.c <- sd(amj$X_99133_00003, na.rm = T)
  mu.c <- mean(amj$X_99133_00003, na.rm = T)
  site.non.event.summary.ls[[i]]$amj.cvccvq <- (sigma.c/mu.c)/(sigma.q/mu.q)
  #r2 amj
  site.non.event.summary.ls[[i]]$amj.r2 <- summary(lm(log(amj$X_99133_00003)~log(amj$X_00060_00003)))$r.square
  #bsig amj
  site.non.event.summary.ls[[i]]$amj.bsig <- summary(lm(log(amj$X_99133_00003)~log(amj$X_00060_00003)))$coefficients[2,'Pr(>|t|)']
  
  #jas
  jas <- site.non.events[[i]][site.non.events[[i]]$month %in% c(7,8,9),]
  #cq jas
  site.non.event.summary.ls[[i]]$jas.cq <- coef(lm(log(jas$X_99133_00003)~log(jas$X_00060_00003)))[2] 
  #cvc.cvq jas
  sigma.q <- sd(jas$X_00060_00003, na.rm = T)
  mu.q <- mean(jas$X_00060_00003, na.rm = T)
  sigma.c <- sd(jas$X_99133_00003, na.rm = T)
  mu.c <- mean(jas$X_99133_00003, na.rm = T)
  site.non.event.summary.ls[[i]]$jas.cvccvq <- (sigma.c/mu.c)/(sigma.q/mu.q)
  #r2 jas
  site.non.event.summary.ls[[i]]$jas.r2 <- summary(lm(log(jas$X_99133_00003)~log(jas$X_00060_00003)))$r.square
  #bsig jas
  site.non.event.summary.ls[[i]]$jas.bsig <- summary(lm(log(jas$X_99133_00003)~log(jas$X_00060_00003)))$coefficients[2,'Pr(>|t|)']
  
  
  #ond
  ond <- site.non.events[[i]][site.non.events[[i]]$month %in% c(10,11,12),]
  #cq ond
  site.non.event.summary.ls[[i]]$ond.cq <- coef(lm(log(ond$X_99133_00003)~log(ond$X_00060_00003)))[2] 
  #cvc.cvq ond
  sigma.q <- sd(ond$X_00060_00003, na.rm = T)
  mu.q <- mean(ond$X_00060_00003, na.rm = T)
  sigma.c <- sd(ond$X_99133_00003, na.rm = T)
  mu.c <- mean(ond$X_99133_00003, na.rm = T)
  site.non.event.summary.ls[[i]]$ond.cvccvq <- (sigma.c/mu.c)/(sigma.q/mu.q)
  #r2 ond
  site.non.event.summary.ls[[i]]$ond.r2 <- summary(lm(log(ond$X_99133_00003)~log(ond$X_00060_00003)))$r.square
  #bsig ond
  site.non.event.summary.ls[[i]]$ond.bsig <- summary(lm(log(ond$X_99133_00003)~log(ond$X_00060_00003)))$coefficients[2,'Pr(>|t|)']
  
}

site.non.event.summary.ls



#make season a column for all site.non.events.ls

for(i in 1:5){
  #site.non.events[[i]]$season <- NA
  for(j in 1:nrow(site.non.events[[i]])){
    if(site.non.events[[i]][j,'month'] %in% c(1,2,3)){
      site.non.events[[i]][j,'season'] <- '2-JFM'
    }else if(site.non.events[[i]][j,'month'] %in% c(4,5,6)){
      site.non.events[[i]][j,'season'] <- '3-AMJ'
    }else if(site.non.events[[i]][j,'month'] %in% c(7,8,9)){
      site.non.events[[i]][j,'season'] <- '4-JAS'
    }else{
      site.non.events[[i]][j,'season'] <- '1-OND'
    }
  }
}

#make a season-year column identifier
for(i in 1:5){
  site.non.events[[i]]$season.yr <- paste(site.non.events[[i]]$season, site.non.events[[i]]$year, sep = '-')
}
#####
#===================================================================================#

#===================================================================================#
#####CLASSIFYING THE ENTIRE RECORDS BY SEASON####

for(i in 1:5){
  site.records[[i]]$month <- strptime(site.records[[i]]$dateTime,format = '%Y-%m-%d', tz = 'UTC') %>% 
    format(.,format = '%m') %>% 
    as.numeric()
}

#make a year column
for(i in 1:5){
  site.records[[i]]$year <- strptime(site.records[[i]]$dateTime,format = '%Y-%m-%d', tz = 'UTC') %>% 
    format(.,format = '%Y') %>% 
    as.numeric()
}

record.summary.df <- data.frame(jfm.cq = NA, amj.cq = NA, jas.cq = NA, ond.cq = NA, jfm.cvccvq = NA, amj.cvccvq = NA, jas.cvccvq = NA, ond.cvccvq = NA, 
                                jfm.r2 = NA, amj.r2 = NA, jas.r2 = NA, ond.r2 = NA, jfm.bsig = NA, amj.bsig = NA, jas.bsig = NA, ond.bsig = NA)
site.record.summary.ls <- list(pan = record.summary.df, sc = record.summary.df, rf = record.summary.df, jeff = record.summary.df, vm = record.summary.df)

#calculate cq and cvc.cvq by season for non-event periods
for(i in 1:5){
  #jfm
  jfm <- site.records[[i]][site.records[[i]]$month %in% c(1,2,3),]
  #cq jfm
  site.record.summary.ls[[i]]$jfm.cq <- coef(lm(log(jfm$X_99133_00003)~log(jfm$X_00060_00003)))[2] 
  #cvc.cvq jfm
  sigma.q <- sd(jfm$X_00060_00003, na.rm = T)
  mu.q <- mean(jfm$X_00060_00003, na.rm = T)
  sigma.c <- sd(jfm$X_99133_00003, na.rm = T)
  mu.c <- mean(jfm$X_99133_00003, na.rm = T)
  site.record.summary.ls[[i]]$jfm.cvccvq <- (sigma.c/mu.c)/(sigma.q/mu.q)
  #r2 jfm
  site.record.summary.ls[[i]]$jfm.r2 <- summary(lm(log(jfm$X_99133_00003)~log(jfm$X_00060_00003)))$r.square
  #bsig jfm
  site.record.summary.ls[[i]]$jfm.bsig <- summary(lm(log(jfm$X_99133_00003)~log(jfm$X_00060_00003)))$coefficients[2,'Pr(>|t|)']
  
  #amj
  amj <- site.records[[i]][site.records[[i]]$month %in% c(4,5,6),]
  #cq amj
  site.record.summary.ls[[i]]$amj.cq <- coef(lm(log(amj$X_99133_00003)~log(amj$X_00060_00003)))[2] 
  #cvc.cvq amj
  sigma.q <- sd(amj$X_00060_00003, na.rm = T)
  mu.q <- mean(amj$X_00060_00003, na.rm = T)
  sigma.c <- sd(amj$X_99133_00003, na.rm = T)
  mu.c <- mean(amj$X_99133_00003, na.rm = T)
  site.record.summary.ls[[i]]$amj.cvccvq <- (sigma.c/mu.c)/(sigma.q/mu.q)
  #r2 amj
  site.record.summary.ls[[i]]$amj.r2 <- summary(lm(log(amj$X_99133_00003)~log(amj$X_00060_00003)))$r.square
  #bsig amj
  site.record.summary.ls[[i]]$amj.bsig <- summary(lm(log(amj$X_99133_00003)~log(amj$X_00060_00003)))$coefficients[2,'Pr(>|t|)']
  
  #jas
  jas <- site.records[[i]][site.records[[i]]$month %in% c(7,8,9),]
  #cq jas
  site.record.summary.ls[[i]]$jas.cq <- coef(lm(log(jas$X_99133_00003)~log(jas$X_00060_00003)))[2] 
  #cvc.cvq jas
  sigma.q <- sd(jas$X_00060_00003, na.rm = T)
  mu.q <- mean(jas$X_00060_00003, na.rm = T)
  sigma.c <- sd(jas$X_99133_00003, na.rm = T)
  mu.c <- mean(jas$X_99133_00003, na.rm = T)
  site.record.summary.ls[[i]]$jas.cvccvq <- (sigma.c/mu.c)/(sigma.q/mu.q)
  #r2 jas
  site.record.summary.ls[[i]]$jas.r2 <- summary(lm(log(jas$X_99133_00003)~log(jas$X_00060_00003)))$r.square
  #bsig jas
  site.record.summary.ls[[i]]$jas.bsig <- summary(lm(log(jas$X_99133_00003)~log(jas$X_00060_00003)))$coefficients[2,'Pr(>|t|)']
  
  #ond
  ond <- site.records[[i]][site.records[[i]]$month %in% c(10,11,12),]
  #cq ond
  site.record.summary.ls[[i]]$ond.cq <- coef(lm(log(ond$X_99133_00003)~log(ond$X_00060_00003)))[2] 
  #cvc.cvq ond
  sigma.q <- sd(ond$X_00060_00003, na.rm = T)
  mu.q <- mean(ond$X_00060_00003, na.rm = T)
  sigma.c <- sd(ond$X_99133_00003, na.rm = T)
  mu.c <- mean(ond$X_99133_00003, na.rm = T)
  site.record.summary.ls[[i]]$ond.cvccvq <- (sigma.c/mu.c)/(sigma.q/mu.q)
  #r2 ond
  site.record.summary.ls[[i]]$ond.r2 <- summary(lm(log(ond$X_99133_00003)~log(ond$X_00060_00003)))$r.square
  #bsig ond
  site.record.summary.ls[[i]]$ond.bsig <- summary(lm(log(ond$X_99133_00003)~log(ond$X_00060_00003)))$coefficients[2,'Pr(>|t|)']
  
}

site.record.summary.ls

#make season a column for all site.records

for(i in 1:5){
  #site.non.events[[i]]$season <- NA
  for(j in 1:nrow(site.records[[i]])){
    if(site.records[[i]][j,'month'] %in% c(1,2,3)){
      site.records[[i]][j,'season'] <- '2-JFM'
    }else if(site.records[[i]][j,'month'] %in% c(4,5,6)){
      site.records[[i]][j,'season'] <- '3-AMJ'
    }else if(site.records[[i]][j,'month'] %in% c(7,8,9)){
      site.records[[i]][j,'season'] <- '4-JAS'
    }else{
      site.records[[i]][j,'season'] <- '1-OND'
    }
  }
}

#make a season-year column identifier
for(i in 1:5){
  site.records[[i]]$season.yr <- paste(site.records[[i]]$season, site.records[[i]]$year, sep = '-')
}
#####
#===================================================================================#

#===================================================================================#
#####CLASSIFYING THE ENTIRE EVENT PERIODS TOGETHER BY SEASON#####

#make a month column
for(i in 1:5){
  site.events.df[[i]]$month <- strptime(site.events.df[[i]]$Date,format = '%Y-%m-%d', tz = 'UTC') %>% 
    format(.,format = '%m') %>% 
    as.numeric()
}

#make a year column
for(i in 1:5){
  site.events.df[[i]]$year <- strptime(site.events.df[[i]]$Date,format = '%Y-%m-%d', tz = 'UTC') %>% 
    format(.,format = '%Y') %>% 
    as.numeric()
}


total.event.summary.df <- data.frame(jfm.cq = NA, amj.cq = NA, jas.cq = NA, ond.cq = NA, jfm.cvccvq = NA, amj.cvccvq = NA, jas.cvccvq = NA, ond.cvccvq = NA, 
                                     jfm.r2 = NA, amj.r2 = NA, jas.r2 = NA, ond.r2 = NA, jfm.bsig = NA, amj.bsig = NA, jas.bsig = NA, ond.bsig = NA)
site.total.event.summary.ls <- list(pan = total.event.summary.df, sc = total.event.summary.df, rf = total.event.summary.df, jeff = total.event.summary.df, vm = total.event.summary.df)

#calculate cq and cvc.cvq by season for non-event periods
for(i in 1:5){
  #jfm
  jfm <- site.events.df[[i]][site.events.df[[i]]$month %in% c(1,2,3),]
  #cq jfm
  site.total.event.summary.ls[[i]]$jfm.cq <- coef(lm(log(jfm$X_99133_00003)~log(jfm$X_00060_00003)))[2] 
  #cvc.cvq jfm
  sigma.q <- sd(jfm$X_00060_00003, na.rm = T)
  mu.q <- mean(jfm$X_00060_00003, na.rm = T)
  sigma.c <- sd(jfm$X_99133_00003, na.rm = T)
  mu.c <- mean(jfm$X_99133_00003, na.rm = T)
  site.total.event.summary.ls[[i]]$jfm.cvccvq <- (sigma.c/mu.c)/(sigma.q/mu.q)
  #r2 jfm
  site.total.event.summary.ls[[i]]$jfm.r2 <- summary(lm(log(jfm$X_99133_00003)~log(jfm$X_00060_00003)))$r.square
  #bsig jfm
  site.total.event.summary.ls[[i]]$jfm.bsig <- summary(lm(log(jfm$X_99133_00003)~log(jfm$X_00060_00003)))$coefficients[2,'Pr(>|t|)']
  
  #amj
  amj <- site.events.df[[i]][site.events.df[[i]]$month %in% c(4,5,6),]
  #cq amj
  site.total.event.summary.ls[[i]]$amj.cq <- coef(lm(log(amj$X_99133_00003)~log(amj$X_00060_00003)))[2] 
  #cvc.cvq amj
  sigma.q <- sd(amj$X_00060_00003, na.rm = T)
  mu.q <- mean(amj$X_00060_00003, na.rm = T)
  sigma.c <- sd(amj$X_99133_00003, na.rm = T)
  mu.c <- mean(amj$X_99133_00003, na.rm = T)
  site.total.event.summary.ls[[i]]$amj.cvccvq <- (sigma.c/mu.c)/(sigma.q/mu.q)
  #r2 amj
  site.total.event.summary.ls[[i]]$amj.r2 <- summary(lm(log(amj$X_99133_00003)~log(amj$X_00060_00003)))$r.square
  #bsig amj
  site.total.event.summary.ls[[i]]$amj.bsig <- summary(lm(log(amj$X_99133_00003)~log(amj$X_00060_00003)))$coefficients[2,'Pr(>|t|)']
  
  #jas
  jas <- site.events.df[[i]][site.events.df[[i]]$month %in% c(7,8,9),]
  #cq jas
  site.total.event.summary.ls[[i]]$jas.cq <- coef(lm(log(jas$X_99133_00003)~log(jas$X_00060_00003)))[2] 
  #cvc.cvq jas
  sigma.q <- sd(jas$X_00060_00003, na.rm = T)
  mu.q <- mean(jas$X_00060_00003, na.rm = T)
  sigma.c <- sd(jas$X_99133_00003, na.rm = T)
  mu.c <- mean(jas$X_99133_00003, na.rm = T)
  site.total.event.summary.ls[[i]]$jas.cvccvq <- (sigma.c/mu.c)/(sigma.q/mu.q)
  #r2 jas
  site.total.event.summary.ls[[i]]$jas.r2 <- summary(lm(log(jas$X_99133_00003)~log(jas$X_00060_00003)))$r.square
  #bsig jas
  site.total.event.summary.ls[[i]]$jas.bsig <- summary(lm(log(jas$X_99133_00003)~log(jas$X_00060_00003)))$coefficients[2,'Pr(>|t|)']
  
  #ond
  ond <- site.events.df[[i]][site.events.df[[i]]$month %in% c(10,11,12),]
  #cq ond
  site.total.event.summary.ls[[i]]$ond.cq <- coef(lm(log(ond$X_99133_00003)~log(ond$X_00060_00003)))[2] 
  #cvc.cvq ond
  sigma.q <- sd(ond$X_00060_00003, na.rm = T)
  mu.q <- mean(ond$X_00060_00003, na.rm = T)
  sigma.c <- sd(ond$X_99133_00003, na.rm = T)
  mu.c <- mean(ond$X_99133_00003, na.rm = T)
  site.total.event.summary.ls[[i]]$ond.cvccvq <- (sigma.c/mu.c)/(sigma.q/mu.q)
  #r2 ond
  site.total.event.summary.ls[[i]]$ond.r2 <- summary(lm(log(ond$X_99133_00003)~log(ond$X_00060_00003)))$r.square
  #bsig ond
  site.total.event.summary.ls[[i]]$ond.bsig <- summary(lm(log(ond$X_99133_00003)~log(ond$X_00060_00003)))$coefficients[2,'Pr(>|t|)']
  
}

lapply(site.events.df, function(x) summary(lm(log(x$X_99133_00003)~log(x$X_00060_00003))))

site.total.event.summary.ls

#make season a column for all site.events.df

for(i in 1:5){
  #site.events.df[[i]]$season <- NA
  for(j in 1:nrow(site.events.df[[i]])){
    if(site.events.df[[i]][j,'month'] %in% c(1,2,3)){
      site.events.df[[i]][j,'season'] <- '2-JFM'
    }else if(site.events.df[[i]][j,'month'] %in% c(4,5,6)){
      site.events.df[[i]][j,'season'] <- '3-AMJ'
    }else if(site.events.df[[i]][j,'month'] %in% c(7,8,9)){
      site.events.df[[i]][j,'season'] <- '4-JAS'
    }else{
      site.events.df[[i]][j,'season'] <- '1-OND'
    }
  }
}


#make a season-year column identifier
for(i in 1:5){
  site.events.df[[i]]$season.yr <- paste(site.events.df[[i]]$season, site.events.df[[i]]$year, sep = '-')
}
#####
#===================================================================================#

#===================================================================================#
#####MAKE A CQ SLOPE SUMMARY FIGURE#####

#data for the jefferson watershed without the low flow period between 7/26/2017-10/19/2017
#---------------------------------------#
cq.full.jeff <- 0.2252684
cq.jeff.JAS <- 0.3019976
cq.baseflow.jeff <- 0.4539984
cq.baseflow.jeff.JAS <- 0.646938
#---------------------------------------#


#data for the jefferson watershed WITH the low flow period between 7/26/2017-10/19/2017
#---------------------------------------#
cq.full.jeff <- 0.76
cq.jeff.JAS <- 1.16
cq.baseflow.jeff <- 1.42
cq.baseflow.jeff.JAS <- 2.13
#---------------------------------------#

cq.full <- lapply(site.records, function(x) as.numeric(coef(lm(log(x$X_99133_00003)~log(x$X_00060_00003)))[2]))

cq.baseflow <- lapply(site.non.events, function(x) as.numeric(coef(lm(log(x$X_99133_00003)~log(x$X_00060_00003)))[2]))

cq.storms <- lapply(site.events.df, function(x) as.numeric(coef(lm(log(x$X_99133_00003)~log(x$X_00060_00003)))[2]))

cols = hcl.colors(12, palette = "Zissou 1", rev = F, alpha = 1)

#With Distribution of events
#exported as 867X1164
par(mfrow = c(5,1), mai = c(0.15,0,0.15,0), omi = c(0.25,1,0.25,1), xpd = F)
for(j in 1:5){
  #Dummy
  hist(site.events.summary.ls[[j]]$cq.slp, breaks = seq(-1.5,1.5,0.1), xlim = c(-1.5,1.5), ylab = 'Number of events', xlab = '', main = wshed.names[j], las = 1, col = 'white', ylim = c(0,25), axes = F)
  #abline(v = c(-0.2,0.2), lty = 2, lwd = 2)
  #abline(v = mean(site.events.summary.ls[[j]]$cq.slp), lwd = 2, lty = 2, col = 'cyan')
  segments(mean(site.events.summary.ls[[j]]$cq.slp),0,mean(site.events.summary.ls[[j]]$cq.slp),30, lwd = 2, lty = 2, col = 'lightpink')
  e.slope <- lm(log(site.events.df[[j]]$X_99133_00003)~log(site.events.df[[j]]$X_00060_00003)) %>%
    coef() %>%
    as.numeric()
  #abline(v = e.slope[2], lwd = 2, lty = 2, col = 'lightpink')
  segments(e.slope[2],0,e.slope[2],30, lwd = 2, lty = 2, col = 'darkgray')
  axis(2,  tck = 0.02, las = 1, labels = c(0,10,20), at = c(0,10,20))
  axis(1, tck = 0.02, labels = T)
  par(new = T)
  #JAS
  hist(site.events.summary.ls[[j]]$cq.slp, breaks = seq(-1.5,1.5,0.1), xlim = c(-1.5,1.5), main = wshed.names[j], las = 1, col = cols[12], ylim = c(0,25), axes = F, xlab = '', ylab = '')
  #AMJ
  par(new = T)
  hist(site.events.summary.ls[[j]][site.events.summary.ls[[j]]$season != '4-JAS',]$cq.slp, breaks = seq(-1.5,1.5,0.1), xlim = c(-1.5,1.5), axes = F, col = cols[8], main = '', xlab = '', ylab = '', ylim = c(0,25))
  #JFM
  par(new = T)
  hist(site.events.summary.ls[[j]][site.events.summary.ls[[j]]$season != '4-JAS' & site.events.summary.ls[[j]]$season != '3-AMJ',]$cq.slp, breaks = seq(-1.5,1.5,0.1), xlim = c(-1.5,1.5), axes = F, col = cols[5], main = '', xlab = '', ylab = '', ylim = c(0,25))
  #OND
  par(new = T)
  hist(site.events.summary.ls[[j]][site.events.summary.ls[[j]]$season == '1-OND',]$cq.slp, breaks = seq(-1.5,1.5,0.1), xlim = c(-1.5,1.5), axes = F, col = cols[1], main = '', xlab = '', ylab = '', ylim = c(0,25))
  #box(lwd = 2)
}
#full record cq
points(cq.full[[j]],50, pch = 21, col = 'black', cex = 3.8, bg = 'darkgray')
#seasonal full record cq slope
points(site.record.summary.ls[[j]][1:4],rep(50,4), bg = cols[c(5,8,12,1)], pch = 21, cex = 3.8, col = 'black')

#full baseflow cq
points(cq.baseflow[[j]],35, pch = 25, col = 'black', cex = 3.8, bg = 'darkgray')
#seasonal baseflow cq slope
points(site.non.event.summary.ls[[j]][1:4],rep(35,4), bg = cols[c(5,8,12,1)], pch = 25, cex = 3.8, col = 'black')

if(j == 4){
  points(c(cq.full.jeff,cq.jeff.JAS), c(50,50), col = c('black',cols[12]), pch = 1, cex = 3.8)
  points(c(cq.baseflow.jeff,cq.baseflow.jeff.JAS), c(35,35), bg = c('black',cols[12]), pch = 6, cex = 3.8)
}


axis(1, labels = T, tck = 0.02)



#saved as 421 X 453
#___BASEFLOW___#
par(mfrow = c(1,1))
#Dummy
plot(seq(-1.5,1.5,0.25),seq(-1.5,1.5,0.25),axes = F, typ = 'n', ylim = c(0.5,5.5), xlim = c(-0.5,1.5), ylab = '', xlab = 'c-Q slope')
#seasonal baseflow cq slope
#ond
ondbf <-lapply(site.non.event.summary.ls, function(x) x[4])%>%
  unlist()
points(ondbf,rep(5,5), bg = alpha(cols[1], 1.0), pch = 25, cex = 2.4, col = 'black')
#jfm
jfmbf <-lapply(site.non.event.summary.ls, function(x) x[1])%>%
  unlist()
points(jfmbf,rep(4,5), bg = alpha(cols[5], 1.0), pch = 25, cex = 2.4, col = 'black')
#amj
amjbf <-lapply(site.non.event.summary.ls, function(x) x[2])%>%
  unlist()
points(amjbf,rep(3,5), bg = alpha(cols[8], 1.0), pch = 25, cex = 2.4, col = 'black')
#jas
jasbf <-lapply(site.non.event.summary.ls, function(x) x[3])%>%
  unlist()
points(jasbf,rep(2,5), bg = alpha(cols[12], 1.0), pch = 25, cex = 2.4, col = 'black')

#full baseflow cq
points(unlist(cq.baseflow),rep(1,5), pch = 25, col = 'black', cex = 2.4, bg = alpha('darkgray', 1.0), main = 'Baseflow')
axis(1, labels = T, tck = 0.02)
axis(2, labels = c('OND','JFM','AMJ','JAS','Annual'), at = seq(5,1), las = 1)
axis(3, labels = F, tck = 0.02)
box()

#___STORMFLOW___# 
#seasonal stormflow cq slope event averaged
#dummy
plot(seq(-1.5,1.5,0.25),seq(-1.5,1.5,0.25),axes = F, typ = 'n', ylim = c(0.5,5.5), xlim = c(-0.5,1.5), ylab = '', xlab = 'c-Q Slope', main= 'Stormflow-Event Averaged')

#OND
points(lapply(site.events.summary.ls, function(x) mean(x[x$season == '1-OND',]$cq.slp)) %>% unlist(),
       rep(5,5), bg = alpha(cols[1], 1.0), pch = 21, cex = 2.4, col = 'black')
#JFM
points(lapply(site.events.summary.ls, function(x) mean(x[x$season == '2-JFM',]$cq.slp)) %>% unlist(),
       rep(4,5), bg = alpha(cols[5], 1.0), pch = 21, cex = 2.4, col = 'black')
#AMJ
points(lapply(site.events.summary.ls, function(x) mean(x[x$season == '3-AMJ',]$cq.slp)) %>% unlist(),
       rep(3,5), bg = alpha(cols[8], 1.0), pch = 21, cex = 2.4, col = 'black')
#JAS
points(lapply(site.events.summary.ls, function(x) mean(x[x$season == '4-JAS',]$cq.slp)) %>% unlist(),
       rep(2,5), bg = alpha(cols[12], 1.0), pch = 21, cex = 2.4, col = 'black')

#points(ssfcq,rep(1, 5), bg = alpha('darkgray', 1.0), pch = 23, cex = 2.5, col = 'black')


#mean storm events
me <- lapply(site.events.summary.ls, function(x) mean(x$cq.slp)) %>%
  unlist()
points(me, rep(1,5), pch = 21, col = 'black', cex = 2.8, bg = alpha('darkgray', 1.0))
axis(1, labels = T, tck = 0.02)
axis(2, labels = c('OND','JFM','AMJ','JAS','Annual'), at = seq(5,1), las = 1)
axis(3, labels = F, tck = 0.02)
box()

#___STORMFLOW___# 
#seasonal stormflow cq slope bulk averaged
#dummy
plot(seq(-1.5,1.5,0.25),seq(-1.5,1.5,0.25),axes = F, typ = 'n', ylim = c(0.5,5.5), xlim = c(-0.5,1.5), ylab = '', xlab = 'c-Q Slope', main = 'Stormflow - Bulk Averaged')

#OND
points(lapply(site.events.df, function(x) as.numeric(coef(lm(log(x[x[,'season'] == '1-OND',]$X_99133_00003)~log(x[x[,'season'] == '1-OND',]$X_00060_00003)))[2])) %>%
         unlist(),
       rep(5,5), bg = alpha(cols[1], 1.0), pch = 21, cex = 2.4, col = 'black')
#JFM
points(lapply(site.events.df, function(x) as.numeric(coef(lm(log(x[x[,'season'] == '2-JFM',]$X_99133_00003)~log(x[x[,'season'] == '2-JFM',]$X_00060_00003)))[2])) %>%
         unlist(),
       rep(4,5), bg = alpha(cols[5], 1.0), pch = 21, cex = 2.4, col = 'black')
#AMJ
points(lapply(site.events.df, function(x) as.numeric(coef(lm(log(x[x[,'season'] == '3-AMJ',]$X_99133_00003)~log(x[x[,'season'] == '3-AMJ',]$X_00060_00003)))[2])) %>%
         unlist(),
       rep(3,5), bg = alpha(cols[8], 1.0), pch = 21, cex = 2.4, col = 'black')
#JAS
points(lapply(site.events.df, function(x) as.numeric(coef(lm(log(x[x[,'season'] == '4-JAS',]$X_99133_00003)~log(x[x[,'season'] == '4-JAS',]$X_00060_00003)))[2])) %>%
         unlist(),
       rep(2,5), bg = alpha(cols[12], 1.0), pch = 21, cex = 2.4, col = 'black')

#full storm cq
fscq <- lapply(site.events.df, function(x) as.numeric(coef(lm(log(x$X_99133_00003)~log(x$X_00060_00003))))[2]) %>%
  unlist()
points(fscq, rep(1,5), pch = 21, col = 'black', cex = 2.8, bg = alpha('darkgray', 1.0))
axis(1, labels = T, tck = 0.02)
axis(2, labels = c('OND','JFM','AMJ','JAS','Annual'), at = seq(5,1), las = 1)
axis(3, labels = F, tck = 0.02)
axis(4, labels = F, tck = 0.02)
box()

#####
#===================================================================================#

#===================================================================================#
#####MAKING A DISTRIBUTION OF EVENTS FOR ALL SITES FIGURE####

par(mar = c(4,4,4,4))
#OND
hist(lapply(site.events.summary.ls, function(x) x$cq.slp) %>%
       unlist(), breaks = seq(-1.5,1.5,0.1), xlim = c(-1.5,1.5), las = 1, col = cols[1], ylim = c(0,80), axes = F, xlab = 'c-Q slope', ylab = 'Number of Events', main = 'c-Q slope distribution for all watersheds')
#JFM
par(new = T)
hist(lapply(site.events.summary.ls, function(x) x[x$season != '1-OND',]$cq.slp) %>%
       unlist(), breaks = seq(-1.5,1.5,0.1), xlim = c(-1.5,1.5), axes = F, col = cols[5], main = '', xlab = '', ylab = '', ylim = c(0,80))
#AMJ
par(new = T)
hist(lapply(site.events.summary.ls, function(x) x[x$season != '1-OND' & x$season != '2-JFM',]$cq.slp) %>%
       unlist(), breaks = seq(-1.5,1.5,0.1), xlim = c(-1.5,1.5), axes = F, col = cols[8], main = '', xlab = '', ylab = '', ylim = c(0,80))
#JAS
par(new = T)
hist(lapply(site.events.summary.ls, function(x) x[x$season == '4-JAS',]$cq.slp) %>%
       unlist(), breaks = seq(-1.5,1.5,0.1), xlim = c(-1.5,1.5), axes = F, col = cols[12], main = '', xlab = '', ylab = '', ylim = c(0,80))
axis(2,  tck = 0.02, las = 1, labels = c(0,20,40,60,80), at = c(0,20,40,60,80))
axis(1, tck = 0.02, labels = T)

#####
#===================================================================================#
