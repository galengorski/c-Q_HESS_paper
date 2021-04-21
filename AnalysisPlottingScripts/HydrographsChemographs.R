#===================================================================================#
# NOTES: This script creates a figure visualizing the hydrographs and chemographs   #
# for each watershed with the events shown in red and the non-events shown in blue  # 
# for nitrate and black for discharge                                               #
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
#####DEALING WITH DATA GAPS#####

#in order to plot up gaps in the record we need to find the gaps and split the
#record up and plot each chunk individually so that the gaps will show
site.record.list <- list()
gap.list <- list()

for(j in 1:5){
  temp.comp <- site.records[[j]]
  temp.comp$gap <- NA
  
  #remove any days that have an NA for either flow or nitrate
  
  temp.comp <- temp.comp[!is.na(temp.comp$X_00060_00003)&!is.na(temp.comp$X_99133_00003),]
  #make a column that has the number of days since the last observation, > 1 is a gap
  for (i in 1:(nrow(temp.comp)-1)){
    temp.comp$gap[i] <- as.numeric(difftime(temp.comp$dateTime[i+1],temp.comp$dateTime[i], units = 'days'))
  }
  #subset the gaps
  gaps <- temp.comp[temp.comp$gap > 1,]
  gaps <- gaps[!is.na(gaps$gap),]
  site.record.list[[j]] <- list()
  site.record.list[[j]][[1]] <- temp.comp[temp.comp$dateTime <= gaps$dateTime[1],]
  for(i in 2:nrow(gaps)){
    site.record.list[[j]][[i]] <- temp.comp[temp.comp$dateTime > gaps$dateTime[i-1] & temp.comp$dateTime <= gaps$dateTime[i],]
  }
  site.record.list[[j]][[nrow(gaps)+1]] <- temp.comp[temp.comp$dateTime > gaps$dateTime[nrow(gaps)],]
  gap.list[[j]] <- as.data.frame(matrix(ncol = 3, nrow = nrow(gaps)))
  colnames(gap.list[[j]]) <- c('gap.start', 'gap.end', 'gap.length')
  for(i in 1:nrow(gaps)){
    if(is.na(gaps$dateTime[i])){
      next
    }else{
      gap.list[[j]]$gap.start[i] <- gaps$dateTime[i]
      gap.list[[j]]$gap.end[i] <- temp.comp[which(temp.comp$dateTime == gaps$dateTime[i])+1, ]$dateTime
      #gap.list[[j]]$gap.length[i] <- difftime(gap.list[[j]]$gap.end[i], gap.list[[j]]$gap.start[i])
    }
  }
  
  gap.list[[j]]$gap.start <- as.POSIXct(gap.list[[j]]$gap.start, origin = '1970-01-01')
  gap.list[[j]]$gap.end <- as.POSIXct(gap.list[[j]]$gap.end, origin = '1970-01-01')
  gap.list[[j]]$gap.length <- as.numeric(difftime(gap.list[[j]]$gap.end, gap.list[[j]]$gap.start, units = 'days'))
}

names(gap.list) <- c('pan','sc','rf','jeff','vm')


#####
#===================================================================================#

#===================================================================================#
#####SI FIGURE HYDROGRAPHS AND CHEMOGRAPHS ALTERNATIVE#####
#saved as 850 x 700
#plotted with the same discharge axes
par(mfrow = c(5,1), mai = c(0.35,0,0.35,0), omi = c(0.25,1,0.25,1), xpd = F)
#Panora
#plotting full record discharge
#blank plot
plot(site.records[[1]]$dateTime, site.records[[1]]$X_00060_00003, typ = 'l', lwd = 2, las = 1, ylab = '', xlab = 'x', axes = F, ylim = c(0, 16250), main = 'UPN')
#plotting events discharge
for(i in 1:length(site.events.ls[[1]])){
  temp.event <- site.events.ls[[1]][[i]]
  if(is.null(nrow(temp.event))){
    next
  }else{
    lines(temp.event$Date, temp.event$X_00060_00003, col = 'tomato', lwd = 2)
  }
}
#axes
axis(2, labels = c(0,2500,5000,7500,10000), at = c(0,2500,5000,7500,10000), lwd = 2, tck = 0.02, cex.axis = 1.4, las = 1)
axis.POSIXct(1, at = seq(as.POSIXlt('2016-01-01'),as.POSIXlt('2020-01-01'), by = 'month'), lwd = 2, tck = 0.01, cex.axis = 1.4, labels = F)
axis.POSIXct(1, x = site.records[[1]]$dateTime, labels = T, lwd = 2, tck = 0.02, cex.axis = 1.4)
axis.POSIXct(3, at = seq(as.POSIXlt('2016-01-01'),as.POSIXlt('2020-01-01'), by = 'month'), lwd = 2, tck = 0.01, cex.axis = 1.4, labels = F)
axis.POSIXct(3, x = site.records[[1]]$dateTime, labels = F, lwd = 2, tck = 0.02, cex.axis = 1.4)
mtext('Discharge (cfs)',side = 2, line = 5, cex = 0.8)
mtext('NO3-N (mg/L)', side = 4, line = 3, cex = 0.8, adj = 1, col = 'dodgerblue3')
box(lwd = 2)
#plotting full record NO3
par(new = T)
plot(site.records[[1]]$dateTime, site.records[[1]]$X_99133_00003, typ = 'n', lwd = 2, col = 'dodgerblue3', ylab = '', xlab = 'x', axes = F, ylim = c(-34,27))
#plotting up the chunks so that the gaps show
for(i in 1:length(site.record.list[[1]])){
  points(site.record.list[[1]][[i]]$dateTime, site.record.list[[1]][[i]]$X_99133_00003, lwd = 2, typ = 'l', col = 'dodgerblue3')
}

axis(4, labels = c(0,5,10,15,20,25), at = c(0,5,10,15,20,25), lwd = 2, tck = 0.02, cex.axis = 1.4, las = 1, col = 'dodgerblue3', par(col.lab = 'dodgerblue3'))
#plotting events NO3
for(i in 1:length(site.events.ls[[1]])){
  temp.event <- site.events.ls[[1]][[i]]
  if(is.null(nrow(temp.event))){
    next
  }else{
    lines(temp.event$Date, temp.event$X_99133_00003, col = 'tomato', lwd = 2)
  }
}



#SacCity
#plotting full record discharge
plot(site.records[[1]]$dateTime, site.records[[1]]$X_00060_00003, typ = 'n', lwd = 2, las = 1, ylab = '', xlab = 'x', axes = F, ylim = c(0, 16250), main = 'USC')
points(site.records[[2]]$dateTime, site.records[[2]]$X_00060_00003, typ = 'l', lwd = 2)
#plotting events discharge
for(i in 1:length(site.events.ls[[2]])){
  temp.event <- site.events.ls[[2]][[i]]
  if(is.null(nrow(temp.event))){
    next
  }else{
    lines(temp.event$Date, temp.event$X_00060_00003, col = 'tomato', lwd = 2)
  }
}
#axes
axis(2, labels = c(0,2500,5000,7500,10000), at = c(0,2500,5000,7500,10000), lwd = 2, tck = 0.02, cex.axis = 1.4, las = 1)
axis.POSIXct(1, at = seq(as.POSIXlt('2016-01-01'),as.POSIXlt('2020-01-01'), by = 'month'), lwd = 2, tck = 0.01, cex.axis = 1.4, labels = F)
axis.POSIXct(1, x = site.records[[2]]$dateTime, labels = T, lwd = 2, tck = 0.02, cex.axis = 1.4)
axis.POSIXct(3, at = seq(as.POSIXlt('2016-01-01'),as.POSIXlt('2020-01-01'), by = 'month'), lwd = 2, tck = 0.01, cex.axis = 1.4, labels = F)
axis.POSIXct(3, x = site.records[[2]]$dateTime, labels = F, lwd = 2, tck = 0.02, cex.axis = 1.4)
mtext('Discharge (cfs)',side = 2, line = 5, cex = 0.8)
mtext('NO3-N (mg/L)', side = 4, line = 3, cex = 0.8, adj = 1, col = 'dodgerblue3')
box(lwd = 2)
#plotting full record NO3
par(new = T)
plot(site.records[[1]]$dateTime, site.records[[1]]$X_99133_00003, typ = 'n', lwd = 2, col = 'dodgerblue3', ylab = '', xlab = 'x', axes = F, ylim = c(-34,27))
#plotting up the chunks so that the gaps show
for(i in 1:length(site.record.list[[2]])){
  points(site.record.list[[2]][[i]]$dateTime, site.record.list[[2]][[i]]$X_99133_00003, lwd = 2, typ = 'l', col = 'dodgerblue3')
}
axis(4, labels = c(0,5,10,15,20,25), at = c(0,5,10,15,20,25), lwd = 2, tck = 0.02, cex.axis = 1.4, las = 1, col = 'dodgerblue3', par(col.lab = 'dodgerblue3'))
#plotting events NO3
for(i in 1:length(site.events.ls[[2]])){
  temp.event <- site.events.ls[[2]][[i]]
  if(is.null(nrow(temp.event))){
    next
  }else{
    lines(temp.event$Date, temp.event$X_99133_00003, col = 'tomato', lwd = 2)
  }
}

#Redfield
#plotting full record discharge
plot(site.records[[1]]$dateTime, site.records[[1]]$X_00060_00003, typ = 'n', lwd = 2, las = 1, ylab = '', xlab = 'x', axes = F, ylim = c(0, 32500), main = 'MRF')
points(site.records[[3]]$dateTime, site.records[[3]]$X_00060_00003, typ = 'l', lwd = 2)
#plotting events discharge
for(i in 1:length(site.events.ls[[3]])){
  temp.event <- site.events.ls[[3]][[i]]
  if(is.null(nrow(temp.event))){
    next
  }else{
    lines(temp.event$Date, temp.event$X_00060_00003, col = 'tomato', lwd = 2)
  }
}
#axes
axis(2, labels = c(0,5000,10000,15000,20000), at = c(0,5000,10000,15000,20000), lwd = 2, tck = 0.02, cex.axis = 1.4, las = 1)
axis.POSIXct(1, at = seq(as.POSIXlt('2016-01-01'),as.POSIXlt('2020-01-01'), by = 'month'), lwd = 2, tck = 0.01, cex.axis = 1.4, labels = F)
axis.POSIXct(1, x = site.records[[3]]$dateTime, labels = T, lwd = 2, tck = 0.02, cex.axis = 1.4)
axis.POSIXct(3, at = seq(as.POSIXlt('2016-01-01'),as.POSIXlt('2020-01-01'), by = 'month'), lwd = 2, tck = 0.01, cex.axis = 1.4, labels = F)
axis.POSIXct(3, x = site.records[[3]]$dateTime, labels = F, lwd = 2, tck = 0.02, cex.axis = 1.4)
mtext('Discharge (cfs)',side = 2, line = 5, cex = 0.8)
mtext('NO3-N (mg/L)', side = 4, line = 3, cex = 0.8, adj = 1, col = 'dodgerblue3')
box(lwd = 2)

#plotting full record NO3
par(new = T)
plot(site.records[[1]]$dateTime, site.records[[1]]$X_99133_00003, typ = 'n', lwd = 2, col = 'dodgerblue3', ylab = '', xlab = 'x', axes = F, ylim = c(-34,27))
#plotting up the chunks so that the gaps show
for(i in 1:length(site.record.list[[3]])){
  points(site.record.list[[3]][[i]]$dateTime, site.record.list[[3]][[i]]$X_99133_00003, lwd = 2, typ = 'l', col = 'dodgerblue3')
}

axis(4, labels = c(0,5,10,15,20,25), at = c(0,5,10,15,20,25), lwd = 2, tck = 0.02, cex.axis = 1.4, las = 1, col = 'dodgerblue3', par(col.lab = 'dodgerblue3'))
#plotting events NO3
for(i in 1:length(site.events.ls[[3]])){
  temp.event <- site.events.ls[[3]][[i]]
  if(is.null(nrow(temp.event))){
    next
  }else{
    lines(temp.event$Date, temp.event$X_99133_00003, col = 'tomato', lwd = 2)
  }
}

#Jefferson
#plotting full record discharge
plot(site.records[[1]]$dateTime, site.records[[1]]$X_00060_00003, typ = 'n', lwd = 2, las = 1, ylab = '', xlab = 'x', axes = F, ylim = c(0, 32500), main = 'MJF')
points(site.records[[4]]$dateTime, site.records[[4]]$X_00060_00003, typ = 'l', lwd = 2)
#plotting events discharge
for(i in 1:length(site.events.ls[[4]])){
  temp.event <- site.events.ls[[4]][[i]]
  if(is.null(nrow(temp.event))){
    next
  }else{
    lines(temp.event$Date, temp.event$X_00060_00003, col = 'tomato', lwd = 2)
  }
}
#axes
axis(2, labels = c(0,5000,10000,15000,20000), at = c(0,5000,10000,15000,20000), lwd = 2, tck = 0.02, cex.axis = 1.4, las = 1)
axis.POSIXct(1, at = seq(as.POSIXlt('2016-01-01'),as.POSIXlt('2020-01-01'), by = 'month'), lwd = 2, tck = 0.01, cex.axis = 1.4, labels = F)
axis.POSIXct(1, x = site.records[[4]]$dateTime, labels = T, lwd = 2, tck = 0.02, cex.axis = 1.4)
axis.POSIXct(3, at = seq(as.POSIXlt('2016-01-01'),as.POSIXlt('2020-01-01'), by = 'month'), lwd = 2, tck = 0.01, cex.axis = 1.4, labels = F)
axis.POSIXct(3, x = site.records[[4]]$dateTime, labels = F, lwd = 2, tck = 0.02, cex.axis = 1.4)
mtext('Discharge (cfs)',side = 2, line = 5, cex = 0.8)
mtext('NO3-N (mg/L)', side = 4, line = 3, cex = 0.8, adj = 1, col = 'dodgerblue3')
box(lwd = 2)

#plotting full record NO3
par(new = T)
plot(site.records[[1]]$dateTime, site.records[[1]]$X_99133_00003, typ = 'n', lwd = 2, col = 'dodgerblue3', ylab = '', xlab = 'x', axes = F, ylim = c(-34,27))
#plotting up the chunks so that the gaps show
for(i in 1:length(site.record.list[[4]])){
  points(site.record.list[[4]][[i]]$dateTime, site.record.list[[4]][[i]]$X_99133_00003, lwd = 2, typ = 'l', col = 'dodgerblue3')
}

axis(4, labels = c(0,5,10,15,20,25), at = c(0,5,10,15,20,25), lwd = 2, tck = 0.02, cex.axis = 1.4, las = 1, col = 'dodgerblue3', par(col.lab = 'dodgerblue3'))
#plotting events NO3
for(i in 1:length(site.events.ls[[4]])){
  temp.event <- site.events.ls[[4]][[i]]
  if(is.null(nrow(temp.event))){
    next
  }else{
    lines(temp.event$Date, temp.event$X_99133_00003, col = 'tomato', lwd = 2)
  }
}

#VanMeter
#plotting full record discharge
plot(site.records[[1]]$dateTime, site.records[[1]]$X_00060_00003, typ = 'n', lwd = 2, las = 1, ylab = '', xlab = 'x', axes = F, ylim = c(0, 65000), main = 'DVM')
points(site.records[[5]]$dateTime, site.records[[5]]$X_00060_00003, typ = 'l', lwd = 2)
#plotting events discharge
for(i in 1:length(site.events.ls[[5]])){
  temp.event <- site.events.ls[[5]][[i]]
  if(is.null(nrow(temp.event))){
    next
  }else{
    lines(temp.event$Date, temp.event$X_00060_00003, col = 'tomato', lwd = 2)
  }
}
#axes
axis(2, labels = c(0,10000,20000,30000,40000), at = c(0,10000,20000,30000,40000), lwd = 2, tck = 0.02, cex.axis = 1.4, las = 1)
axis.POSIXct(1, at = seq(as.POSIXlt('2016-01-01'),as.POSIXlt('2020-01-01'), by = 'month'), lwd = 2, tck = 0.01, cex.axis = 1.4, labels = F)
axis.POSIXct(1, x = site.records[[5]]$dateTime, labels = T, lwd = 2, tck = 0.02, cex.axis = 1.4)
axis.POSIXct(3, at = seq(as.POSIXlt('2016-01-01'),as.POSIXlt('2020-01-01'), by = 'month'), lwd = 2, tck = 0.01, cex.axis = 1.4, labels = F)
axis.POSIXct(3, x = site.records[[5]]$dateTime, labels = F, lwd = 2, tck = 0.02, cex.axis = 1.4)
mtext('Discharge (cfs)',side = 2, line = 5, cex = 0.8)
mtext('NO3-N (mg/L)', side = 4, line = 3, cex = 0.8, adj = 1, col = 'dodgerblue3')
box(lwd = 2)

#plotting full record NO3
par(new = T)
plot(site.records[[1]]$dateTime, site.records[[1]]$X_99133_00003, typ = 'n', lwd = 2, col = 'dodgerblue3', ylab = '', xlab = 'x', axes = F, ylim = c(-34,27))
#plotting up the chunks so that the gaps show
for(i in 1:length(site.record.list[[5]])){
  points(site.record.list[[5]][[i]]$dateTime, site.record.list[[5]][[i]]$X_99133_00003, lwd = 2, typ = 'l', col = 'dodgerblue3')
}

axis(4, labels = c(0,5,10,15,20,25), at = c(0,5,10,15,20,25), lwd = 2, tck = 0.02, cex.axis = 1.4, las = 1, col = 'dodgerblue3', par(col.lab = 'dodgerblue3'))
#plotting events NO3
for(i in 1:length(site.events.ls[[5]])){
  temp.event <- site.events.ls[[5]][[i]]
  if(is.null(nrow(temp.event))){
    next
  }else{
    lines(temp.event$Date, temp.event$X_99133_00003, col = 'tomato', lwd = 2)
  }
}
#####
#===================================================================================#

