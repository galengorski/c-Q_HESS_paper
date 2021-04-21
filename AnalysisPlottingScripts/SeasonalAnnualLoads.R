#===================================================================================#
# NOTES: This script creates a figure looking at the seasonal nitrate load for each
# watershed. It also calculates annual loads and plot those up with watershed
# characteristics such as drainage density
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
#####MAKE A SEASONAL LOAD FIGURE#####


season.yr <- paste(rep(c('1-OND','2-JFM','3-AMJ','4-JAS'), each = 4),2016:2019, sep = '-')

#load calculations
#------calculate load for the full record-----#
site.records <- lapply(site.records, function(x) cbind(x, 'load.kgday' = x$X_99133_00003*x$X_00060_00003*2.446848))
#number of days with observations per season.yr
tot.seasonal.sum <- lapply(site.records, function(x) aggregate(x$load.kgday, by = list(x$season.yr), FUN = function(y) c(sum(y, na.rm = T), length(na.omit(y)))))
tot.seasonal.sum <- lapply(tot.seasonal.sum, function(x) as.data.frame(x))
#divide by the number of days in that season
tot.sum <- lapply(tot.seasonal.sum, function(z) z$x <- cbind(z$x,z$x[,2]/c(rep(92,4),91,rep(90,3),rep(91,4),rep(92,4))))

#see which seasons had greater than 75% coverage
tot.sum <- lapply(tot.sum, function(x) x <- cbind(x,lapply(x[,3], function(y) if(y>0.75)1 else NA)))

#divide by the % coverage
tot.sum <- lapply(tot.sum, function(x) x <- cbind(x, apply(matrix(c(unlist(x)), ncol = 4),1, function(y) if(y[3]==0)NA else(y[1]*y[4])/y[3])))

#convert to list of matrices
tot.sum.mat <- lapply(tot.sum, function(z) z <- matrix(c(unlist(z)), ncol = 5))

#divide by the total area of each catchment
wshed.area <- c(1116,1840,2548,4188,8870)

for(i in 1:length(tot.sum.mat)){
  tot.sum.mat[[i]] <- cbind(tot.sum.mat[[i]],tot.sum.mat[[i]][,5]/wshed.area[i])
  colnames(tot.sum.mat[[i]]) <- c('tot.load','obs.days','frac.obs','flag','cor.load','cor.load.area')
  row.names(tot.sum.mat[[i]]) <- c(tot.seasonal.sum[[1]]$Group.1)
}


#subset further into yr.season
tot.seasonyr.mat <- lapply(tot.sum.mat, function(x) y <- list(OND = x[1:4,], JFM = x[5:8,], AMJ = x[9:12,], JAS = x[13:16,]))

lapply(tot.seasonyr.mat, function(x) lapply(x, function(y) mean(y[,'cor.load.area'], na.rm = T)))

#----non-events----#
site.non.events <- lapply(site.non.events, function(x) cbind(x, 'load.kgday' = x$X_99133_00003*x$X_00060_00003*2.446848))
#number of days with observations per season.yr
site.non.events.seasonal.sum <- lapply(site.non.events, function(x) aggregate(x$load.kgday, by = list(x$season.yr), FUN = function(y) c(sum(y, na.rm = T), length(na.omit(y)))))
site.non.events.seasonal.sum <- lapply(site.non.events.seasonal.sum, function(x) as.data.frame(x))
#divide by the number of days in that season
ne.sum <- lapply(site.non.events.seasonal.sum, function(z) z$x <- cbind(z$x,z$x[,2]/c(rep(92,4),91,rep(90,3),rep(91,4),rep(92,4))))

#add the percent coverage from the full record, it will be applied to the events and non-event periods
for(i in 1:length(ne.sum)){
  ne.sum[[i]] <- cbind(ne.sum[[i]][,1],tot.sum.mat[[i]][,c('frac.obs','flag')])
}

#convert to matrices
ne.sum.mat <- lapply(ne.sum, function(x) matrix(c(unlist(x)), ncol = 3))

#divide by the % coverage
ne.sum.mat <- lapply(ne.sum.mat, function(x) x <- cbind(x, apply(x, 1, function(y) if(y[2]==0)NA else (y[1]*y[3])/y[2])))

#divide by the total area of each catchment
wshed.area <- c(1116,1840,2548,4188,8870)

for(i in 1:length(ne.sum.mat)){
  ne.sum.mat[[i]] <- cbind(ne.sum.mat[[i]],ne.sum.mat[[i]][,4]/wshed.area[i])
  colnames(ne.sum.mat[[i]]) <- c('tot.load','frac.obs','flag','cor.load','cor.load.area')
  row.names(ne.sum.mat[[i]]) <- c(tot.seasonal.sum[[1]]$Group.1)
}


#subset furher into yr.season
ne.seasonyr.mat <- lapply(ne.sum.mat, function(x) y <- list(OND = x[1:4,], JFM = x[5:8,], AMJ = x[9:12,], JAS = x[13:16,]))

lapply(ne.seasonyr.mat, function(x) lapply(x, function(y) mean(y[,'cor.load.area'], na.rm = T))) %>% unlist()





#----events----#
site.events.df <- lapply(site.events.df, function(x) cbind(x, 'load.kgday' = x$X_99133_00003*x$X_00060_00003*2.446848))
#number of days with observations per season.yr
site.events.seasonal.sum <- lapply(site.events.df, function(x) aggregate(x$load.kgday, by = list(x$season.yr), FUN = function(y) c(sum(y, na.rm = T), length(na.omit(y)))))
site.events.seasonal.sum <- lapply(site.events.seasonal.sum, function(x) as.data.frame(x))

#there are some seasons with no events so we have to put a zero in those seasons making sure to differentiate between zero events
#and periods where no data were measured

#add in the season with zeroes and order based on season year
for(i in 1:length(site.events.seasonal.sum)){
  if(nrow(site.events.seasonal.sum[[i]]) == 16){
    next
  }else{
    missing <- which(season.yr %nin% site.events.seasonal.sum[[i]]$Group.1)
    for(j in 1:length(missing)){
      site.events.seasonal.sum[[i]] <- rbind(site.events.seasonal.sum[[i]], c(season.yr[missing[j]],0,0))
    }
    site.events.seasonal.sum[[i]] <- site.events.seasonal.sum[[i]][order(site.events.seasonal.sum[[i]]$Group.1),]
  }
}

#convert to matrices
e.sum.mat <- lapply(site.events.seasonal.sum, function(z) matrix(c(unlist(z$x)), ncol = 2))


#add the percent coverage from the full record, it will be applied to the events and non-event periods
for(i in 1:length(e.sum.mat)){
  e.sum.mat[[i]] <- cbind(e.sum.mat[[i]], tot.sum.mat[[i]][,'frac.obs'], tot.sum.mat[[i]][,'flag'])
  e.sum.mat[[i]] <- matrix(as.numeric(c(e.sum.mat[[i]])),ncol = 4)
}


#divide by the % coverage
e.sum.mat <- lapply(e.sum.mat, function(t) t <- cbind(t, apply(t, 1, function(y) if(y[3]==0)NA else (y[1]*y[4])/y[3])))

#divide by the total area of each catchment
wshed.area <- c(1116,1840,2548,4188,8870)

for(i in 1:length(e.sum.mat)){
  e.sum.mat[[i]] <- cbind(e.sum.mat[[i]],e.sum.mat[[i]][,5]/wshed.area[i])
  colnames(e.sum.mat[[i]]) <- c('tot.load','n.days', 'frac.obs','flag','cor.load','cor.load.area')
  row.names(e.sum.mat[[i]]) <- c(tot.seasonal.sum[[1]]$Group.1)
}


#subset further into yr.season
e.seasonyr.mat <- lapply(e.sum.mat, function(x) y <- list(OND = x[1:4,], JFM = x[5:8,], AMJ = x[9:12,], JAS = x[13:16,]))

lapply(e.seasonyr.mat, function(x) lapply(x, function(y) mean(y[,'cor.load.area'], na.rm = T))) %>% unlist


ne.loads <- lapply(ne.seasonyr.mat, function(x) lapply(x, function(y) mean(y[,'cor.load.area'], na.rm = T))) %>%
  unlist() %>%
  matrix(.,nrow = 4)

e.loads <- lapply(e.seasonyr.mat, function(x) lapply(x, function(y) mean(y[,'cor.load.area'], na.rm = T))) %>%
  unlist() %>%
  matrix(.,nrow = 4)

#error bars

ne.loads.ranges <- lapply(ne.seasonyr.mat, function(x) lapply(x, function(y) range(y[,'cor.load.area'], na.rm = T))) %>%
  unlist() %>%
  matrix(.,nrow = 2)

e.loads.ranges <- lapply(e.seasonyr.mat, function(x) lapply(x, function(y) range(y[,'cor.load.area'], na.rm = T))) %>%
  unlist() %>%
  matrix(.,nrow = 2)


#non.event seasonal loads
cols = hcl.colors(12, palette = "Zissou 1", rev = F, alpha = 1)
#saved as 650x700
par(mfrow = c(2,1), mai = c(0.35,0.5,0.5,0.35), omi = c(0.5,0.75,0.5,0.5), xpd = T)
barplot(ne.loads[,c(2,4,1,5,3)], beside = T, ylim = c(0, 2.0e3), col = c(cols[c(1,5,8,12)]), axes = F, main = 'Baseflow', names.arg = '')

# #error bars -- ordered by drainage density
arrows(c(1.5,2.5,3.5,4.5), ne.loads.ranges[1,5:8], c(1.5,2.5,3.5,4.5), ne.loads.ranges[2,5:8],code = 3, angle = 90, length = 0.03, lwd = 1.5)
arrows(c(6.5,7.5,8.5,9.5), ne.loads.ranges[1,13:16], c(6.5,7.5,8.5,9.5), ne.loads.ranges[2,13:16],code = 3, angle = 90, length = 0.03, lwd = 1.5)
arrows(c(11.5,12.5,13.5,14.5), ne.loads.ranges[1,1:4], c(11.5,12.5,13.5,14.5), ne.loads.ranges[2,1:4],code = 3, angle = 90, length = 0.03, lwd = 1.5)
arrows(c(16.5,17.5,18.5,19.5), ne.loads.ranges[1,17:20], c(16.5,17.5,18.5,19.5), ne.loads.ranges[2,17:20],code = 3, angle = 90, length = 0.03, lwd = 1.5)
arrows(c(21.5,22.5,23.5,24.5), ne.loads.ranges[1,9:12], c(21.5,22.5,23.5,24.5), ne.loads.ranges[2,9:12],code = 3, angle = 90, length = 0.03, lwd = 1.5)


axis(1, lwd = 2, labels = c('USC','MJF','UPN','DVM','MRF'), at = c(3,8,13,18,23), tck = F)
axis(2, labels = T, lwd = 2, tck = 0.02, las = 1, font = 1)
axis(3, lwd = 2, tck = 0.02, labels = c('','','','',''), at = c(3,8,13,18,23))
axis(4, labels = F, lwd = 2, tck = 0.02)
graphics::box(lwd = 2)

legend('topright', pch = 15, col = cols[c(1,5,8,12)], legend = c('OND','JFM','AMJ','JAS'), bty = 'n', cex = 1.0)

mtext('Cumulative Seasonal Load (kg N/km2)', 2, line = 3.7)
mtext('Baseflow', 3, line = 0, adj = 0, cex = 1.3)


#event seasonal loads
barplot(e.loads[,c(2,4,1,5,3)], beside = T, ylim = c(0, 2.0e03), col = c(cols[c(1,5,8,12)]), axes = F, main = 'Stormflow', names.arg = '')


# #error bars -- ordered by drainage density
arrows(c(1.5,2.5,3.5,4.5), e.loads.ranges[1,5:8], c(1.5,2.5,3.5,4.5), e.loads.ranges[2,5:8],code = 3, angle = 90, length = 0.03, lwd = 1.5)
arrows(c(6.5,7.5,8.5,9.5), e.loads.ranges[1,13:16], c(6.5,7.5,8.5,9.5), e.loads.ranges[2,13:16],code = 3, angle = 90, length = 0.03, lwd = 1.5)
arrows(c(11.5,12.5,13.5,14.5), e.loads.ranges[1,1:4], c(11.5,12.5,13.5,14.5), e.loads.ranges[2,1:4],code = 3, angle = 90, length = 0.03, lwd = 1.5)
arrows(c(16.5,17.5,18.5,19.5), e.loads.ranges[1,17:20], c(16.5,17.5,18.5,19.5), e.loads.ranges[2,17:20],code = 3, angle = 90, length = 0.03, lwd = 1.5)
arrows(c(21.5,22.5,23.5,24.5), e.loads.ranges[1,9:12], c(21.5,22.5,23.5,24.5), e.loads.ranges[2,9:12],code = 3, angle = 90, length = 0.03, lwd = 1.5)


axis(1, lwd = 2, labels = c('USC','MJF','UPN','DVM','MRF'), at = c(3,8,13,18,23), tck = F)
axis(2, labels = T, lwd = 2, tck = 0.02, las = 1, font = 1)
axis(3, lwd = 2, tck = 0.02, labels = c('','','','',''), at = c(3,8,13,18,23))
axis(4, labels = F, lwd = 2, tck = 0.02)
graphics::box(lwd = 2)
mtext('Cumulative Seasonal Load (kg N/km2)', 2, line = 3.7)
mtext('Stormflow', 3, line = 0, adj = 0, cex = 1.3)



#####
#===================================================================================#

#===================================================================================#
#####LOOK AT ANNUAL LOADS#####

#annual loads are not just the sum of seasonal loads because of the threshold for data coverage
#if 0% of a season has observations then 75% of the year can still have observations

#number of days with observations per season.yr
tot.annual.sum <- lapply(site.records, function(x) aggregate(x$load.kgday, by = list(x$year), FUN = function(y) c(sum(y, na.rm = T), length(na.omit(y)))))
tot.annual.sum <- lapply(tot.annual.sum, function(x) as.data.frame(x))
#divide by the number of days in that season
tot.ann.sum <- lapply(tot.annual.sum, function(z) z$x <- cbind(z$x,z$x[,2]/c(366,365,365,365)))

#see which seasons had greater than 75% coverage
tot.ann.sum <- lapply(tot.ann.sum, function(x) x <- cbind(x,lapply(x[,3], function(y) if(y>0.45)1 else NA)))

#divide by the % coverage
tot.ann.sum <- lapply(tot.ann.sum, function(x) x <- cbind(x, apply(matrix(c(unlist(x)), ncol = 4),1, function(y) if(y[3]==0)NA else(y[1]*y[4])/y[3])))

#convert to list of matrices
tot.ann.sum.mat <- lapply(tot.ann.sum, function(z) z <- matrix(c(unlist(z)), ncol = 5))

#divide by the total area of each catchment
wshed.area <- c(1116,1840,2548,4188,8870)

for(i in 1:length(tot.ann.sum.mat)){
  tot.ann.sum.mat[[i]] <- cbind(tot.ann.sum.mat[[i]],tot.ann.sum.mat[[i]][,5]/wshed.area[i])
  colnames(tot.ann.sum.mat[[i]]) <- c('tot.load','obs.days','frac.obs','flag','cor.load','cor.load.area')
  row.names(tot.ann.sum.mat[[i]]) <- c(tot.annual.sum[[1]]$Group.1)
}

#----non-events----#
#site.non.events <- lapply(site.non.events, function(x) cbind(x, 'load.kgday' = x$X_99133_00003*x$X_00060_00003*2.446848))
#number of days with observations per season.yr
site.non.events.annual.sum <- lapply(site.non.events, function(x) aggregate(x$load.kgday, by = list(x$year), FUN = function(y) c(sum(y, na.rm = T), length(na.omit(y)))))
site.non.events.annual.sum <- lapply(site.non.events.annual.sum, function(x) as.data.frame(x))
#divide by the number of days in that season
ne.ann.sum <- lapply(site.non.events.annual.sum, function(z) z$x <- cbind(z$x,z$x[,2]/c(366,365,365,365)))

#add the percent coverage from the full record, it will be applied to the events and non-event periods
for(i in 1:length(ne.ann.sum)){
  ne.ann.sum[[i]] <- cbind(ne.ann.sum[[i]][,1],tot.ann.sum.mat[[i]][,c('frac.obs','flag')])
}

#convert to matrices
ne.ann.sum.mat <- lapply(ne.ann.sum, function(x) matrix(c(unlist(x)), ncol = 3))

#divide by the % coverage
ne.ann.sum.mat <- lapply(ne.ann.sum.mat, function(x) x <- cbind(x, apply(x, 1, function(y) if(y[2]==0)NA else (y[1]*y[3])/y[2])))

#divide by the total area of each catchment
wshed.area <- c(1116,1840,2548,4188,8870)

for(i in 1:length(ne.sum.mat)){
  ne.ann.sum.mat[[i]] <- cbind(ne.ann.sum.mat[[i]],ne.ann.sum.mat[[i]][,4]/wshed.area[i])
  colnames(ne.ann.sum.mat[[i]]) <- c('tot.load','frac.obs','flag','cor.load','cor.load.area')
  row.names(ne.ann.sum.mat[[i]]) <- c(tot.annual.sum[[1]]$Group.1)
}


#----events----#
#number of days with observations per season.yr
site.events.annual.sum <- lapply(site.events.df, function(x) aggregate(x$load.kgday, by = list(x$year), FUN = function(y) c(sum(y, na.rm = T), length(na.omit(y)))))
site.events.annual.sum <- lapply(site.events.annual.sum, function(x) as.data.frame(x))

#there are some seasons with no events so we have to put a zero in those seasons making sure to differentiate between zero events
#and periods where no data were measured

#convert to matrices
e.ann.sum.mat <- lapply(site.events.annual.sum, function(z) matrix(c(unlist(z$x)), ncol = 2))


#add the percent coverage from the full record, it will be applied to the events and non-event periods
for(i in 1:length(e.ann.sum.mat)){
  e.ann.sum.mat[[i]] <- cbind(e.ann.sum.mat[[i]], tot.ann.sum.mat[[i]][,'frac.obs'], tot.ann.sum.mat[[i]][,'flag'])
  e.ann.sum.mat[[i]] <- matrix(as.numeric(c(e.ann.sum.mat[[i]])),ncol = 4)
}

#divide by the % coverage
e.ann.sum.mat <- lapply(e.ann.sum.mat, function(t) t <- cbind(t, apply(t, 1, function(y) if(y[3]==0)NA else (y[1]*y[4])/y[3])))

#divide by the total area of each catchment
wshed.area <- c(1116,1840,2548,4188,8870)

for(i in 1:length(e.ann.sum.mat)){
  e.ann.sum.mat[[i]] <- cbind(e.ann.sum.mat[[i]],e.ann.sum.mat[[i]][,5]/wshed.area[i])
  colnames(e.ann.sum.mat[[i]]) <- c('tot.load','n.days', 'frac.obs','flag','cor.load','cor.load.area')
  row.names(e.ann.sum.mat[[i]]) <- c(tot.annual.sum[[1]]$Group.1)
}

#####
#===================================================================================#

#===================================================================================#
#####SIMPLE CROSSPLOTS WITH DRAINAGE DENSITY#####

#prepare annual loads for plotting
#non.events
ne.avg.load <- lapply(ne.ann.sum.mat, function(x) mean(x[,5])) %>% 
  unlist()
#events
e.avg.load <- lapply(e.ann.sum.mat, function(x) mean(x[,6])) %>% 
  unlist()


#upper and lower limits for annual loads non.events
ne.ranges <- lapply(ne.seasonyr.mat, function(y)
  lapply(y, function(x) x[,5]) %>%
    unlist() %>% 
    matrix(.,nrow = 4) %>% 
    rowSums() %>%
    range(.,na.rm = T))%>%
  unlist() %>%
  matrix(.,nrow = 2)

#upper and lower limits for annual loads events
e.ranges <- lapply(e.seasonyr.mat, function(y)
  lapply(y, function(x) x[,6]) %>%
    unlist() %>% 
    matrix(.,nrow = 4) %>% 
    rowSums() %>%
    range(.,na.rm = T)) %>%
  unlist() %>%
  matrix(.,nrow = 2)



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

col <- colorRamp(c("blue", "red"))(c((unlist(cq.baseflow) %>% sort())/(max(unlist(cq.baseflow))-min(unlist(cq.baseflow))))/1.3654248)
bf <- sort(unlist(cq.baseflow))
col <- colorRamp(c("red", "blue"))( as.vector(bf/bf[5]) )
paste(col[1,],col[2,])
cols <- colorRampPalette(brewer.pal(9,"Blues"))(10)

cols.bf <- c('#0070C0','#BBD9EF','#59A2D6','#6DADDB','#60A6D8')
cols.sf <- c('#A7CEEA','#EEF6FB','#9EC9E7','#D3E6F4','#ADD1EB')

#exported as 716x695
par(mfrow = c(1,1))
plot(drain.d, ne.avg.load, axes = F, cex = 2.6, xlim = c(0.3,1.3), ylim = c(700,3300), ylab = '', xlab = '', typ = 'n')
abline(lm(ne.avg.load~drain.d), lty = 2, lwd = 2, col = 'black')
points(drain.d, ne.avg.load, cex = 2.6, xlim = c(0.3,1.3), ylim = c(700,3300), ylab = '', xlab = '', pch = 25, bg = cols.bf)
abline(lm(e.avg.load~drain.d), lty = 2, lwd = 2, col = 'darkgray')
points(drain.d, e.avg.load, cex = 2.6, xlim = c(0.3,1.3), ylim = c(700,3300), ylab = '', xlab = '', pch = 22, bg = cols.sf)

mtext('Built drainage density (km/km2)', 1, 2)
mtext('Cumulative annual load (kg-N/km2)', 2, 3)
axis(1, labels = T, tck = 0.02)
axis(2, labels = T, tck = 0.02, las = 1)
axis(3, labels = F, tck = 0.02)
axis(4, labels = F, tck = 0.02)
box()

legend('topleft', pch = c(25,22), legend= c("Baseflow", "Stormflow"))

summary(lm(ne.avg.load~drain.d))

#####
#===================================================================================#