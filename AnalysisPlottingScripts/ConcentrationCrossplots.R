#===================================================================================#
# NOTES: This script creates figures analyzing the seasonal nitrate concentration
# and its relationship to watershed characteristics
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
#####CONCENTRATION CROSSPLOTS FOR EACH SEASON#####

#drainage density and landuse were calculated separately 
drain.d <- c(0.705704462,
             1.109050557,
             0.370089421,
             0.92691941,
             0.701199201)

crop.buffer <- c(0.3192672,
                 0.3935896,
                 0.3107614,
                 0.3508266,
                 0.2813553)

total.crop <- c(0.8745007,
                0.915053,
                0.8516263,
                0.9096986,
                0.884581)

#watershed area
wshed.area <- c(1116,1840,2548,4188,8870)


#making concentration cross plots for each season

#___BASEFLOW VS AREA___# 
#OND
bfond <- lapply(site.non.events, function(x) median(x[x$season == '1-OND',]$X_99133_00003, na.rm = T)) %>%
  unlist()
#JFM
bfjfm <- lapply(site.non.events, function(x) median(x[x$season == '2-JFM',]$X_99133_00003, na.rm = T)) %>%
  unlist()
#AMJ
bfamj <- lapply(site.non.events, function(x) median(x[x$season == '3-AMJ',]$X_99133_00003, na.rm = T)) %>%
  unlist()
#JAS
bfjas <- lapply(site.non.events, function(x) median(x[x$season == '4-JAS',]$X_99133_00003, na.rm = T)) %>%
  unlist()

#define colors
cols = hcl.colors(12, palette = "Zissou 1", rev = F, alpha = 1)
par(mar = c(4,4,4,4))
plot(1,1,axes = F, typ = 'n', ylab = 'Seasonal baseflow [NO3] (mg/L)', xlab = 'Watershed area km2', ylim = c(2,15), xlim = c(0,9500), main = 'Seasonal baseflow [NO3] \n vs Watershed area')
#OND
points(wshed.area, bfond, bg = alpha(cols[1], 1.0), pch = 21, cex = 2.1, col = 'black')
#JFM
points(wshed.area, bfjfm, bg = alpha(cols[5], 1.0), pch = 21, cex = 2.1, col = 'black')
#AMJ
points(wshed.area, bfamj, bg = alpha(cols[8], 1.0), pch = 21, cex = 2.1, col = 'black')
#JAS
points(wshed.area, bfjas, bg = alpha(cols[12], 1.0), pch = 21, cex = 2.1, col = 'black')
#axes
axis(1, labels = T, tck = 0.02)
axis(2, labels = T, tck = 0.02, las = 1)
axis(3, labels = F, tck = 0.02)
axis(4, labels = F, tck = 0.02)
box()
legend('topright', pch = 16, col = cols[c(1,5,8,12)], legend = c('OND','JFM','AMJ','JAS'), bty = 'n', cex = 1.0)


#___STORMFLOW VS DRAINAGE INFRASTRUCUTRE DENSITY___# 
#OND
sfond <- lapply(site.events.df, function(x) median(x[x$season == '1-OND',]$X_99133_00003, na.rm = T)) %>%
  unlist()
#JFM
sfjfm <- lapply(site.events.df, function(x) median(x[x$season == '2-JFM',]$X_99133_00003, na.rm = T)) %>%
  unlist()
#AMJ
sfamj <- lapply(site.events.df, function(x) median(x[x$season == '3-AMJ',]$X_99133_00003, na.rm = T)) %>%
  unlist()
#JAS
sfjas <- lapply(site.events.df, function(x) median(x[x$season == '4-JAS',]$X_99133_00003, na.rm = T)) %>%
  unlist()

plot(1,1,axes = F, typ = 'n', ylab = 'Seasonal Stormflow [NO3]', xlab = 'Drainage density (km/km2)', ylim = c(2,15), xlim = c(0,1.5), main = 'Seasonal stormflow [NO3] \n vs Drainage density')
#OND
points(drain.d, sfond, bg = alpha(cols[1], 1.0), pch = 21, cex = 2.1, col = 'black')
#JFM
points(drain.d, sfjfm, bg = alpha(cols[5], 1.0), pch = 21, cex = 2.1, col = 'black')
#AMJ
points(drain.d, sfamj, bg = alpha(cols[8], 1.0), pch = 21, cex = 2.1, col = 'black')
#JAS
points(drain.d, sfjas, bg = alpha(cols[12], 1.0), pch = 21, cex = 2.1, col = 'black')
#axes
axis(1, labels = T, tck = 0.02)
axis(2, labels = T, tck = 0.02, las = 1)
axis(3, labels = F, tck = 0.02)
axis(4, labels = F, tck = 0.02)
box()
legend('topright', pch = 16, col = cols[c(1,5,8,12)], legend = c('OND','JFM','AMJ','JAS'), bty = 'n', cex = 1.0)

#####
#===================================================================================#