# ====================================
# Theng et al. 2021 
# meryltheng@gmail.com
# Process simulated movement data
# ====================================

library(dplyr)
library(sp)
library(adehabitatHR)
library(raster)
library(rgeos)
library(foreach)
source("R/0_functions.R")

# Set up for parallel (optional):
library(doParallel) 
( ncore <- detectCores() - 2 ) # 4 cores
cl <- makeCluster(ncore)
cl
clusterSetRNGStream(cl = cl) # set up random number stream 
registerDoParallel(cl) # register as backend for 'foreach'

# Prep
env <- read.csv('Misc/SCRland.csv', header = T) # resource landscape
envm<-t(as.matrix(env))
r <- raster(envm[nrow(envm):1,], xmn=0, xmx=nrow(envm), ymn=0, ymx=ncol(envm))

# Process movement data
n_foragers = 100
run1a = TRUE

foreach(i=1:100, .combine=rbind, .errorhandling='pass', .packages=c(.packages())) %do% { 
  
  if(run1a == TRUE){
  ###### 1a Read and combine F&M movement output ######
  dfF <- read.csv(paste('MovementInput/Females/Tracks',sprintf("%02d", i),'.csv', sep=''), header = T)
  dfM <- read.csv(paste('MovementInput/Males/Tracks',sprintf("%02d", i),'.csv', sep=''), header = T)
  
  # discard burn-in
  dfF %>% 
    group_by(id) %>% 
    filter(row_number() %in% 21601:43200) -> dfF
  dfF <- dfF[, c('id','x','y')]
  
  dfM %>% 
    group_by(id) %>% 
    filter(row_number() %in% 21601:43200) -> dfM
  dfM <- dfM[, c('id','x','y')]
  
  dfF$x<- dfF$x - 50 # correct for empty buffer
  dfF$y<- dfF$y - 50
  dfM$x<- dfM$x - 50
  dfM$y<- dfM$y - 50
  
  dfF <- as.data.frame(dfF)
  dfM <- as.data.frame(dfM)
  
  dfF$sex <- rep('F', nrow(dfF))
  dfM$sex <- rep('M', nrow(dfM))
  
  # first create unique ids for males (to diff against females)
  replace_idM <-  sort(unique(dfM$id)); replace_idF <- sort(unique(dfF$id))
  replace_keyM <- 61:100; replace_keyF <- 1:60
  dfM$key <- rep(NA, nrow(dfM))
  dfF$key <- rep(NA, nrow(dfF))
  k=0
  for (j in unique(dfM$id)){
    k=k+1
    dfM$key[dfM$id == j] <- replace_keyM[k]
  }
  dfM$id <- dfM$key
  dfM <- dfM[,1:4]
  
  k=0
  for (j in unique(dfF$id)){
    k=k+1
    dfF$key[dfF$id == j] <- replace_keyF[k]
  }
  dfF$id <- dfF$key
  dfF <- dfF[,1:4]
  
  df <- rbind(dfF, dfM) # full movement dataset
  
  ###### 2 Convert timesteps to DateTime ######
  n_foragers = 100
  
  # create 12-hour activity days
  datestep <- sapply(seq(1,60,2), function(t){seq(720*(t-1)*60,720*t*60,60)})
  datestep <- as.vector(datestep[-1,])
  df$date <- rep(datestep, n_foragers)
  class(df$date) <- c("POSIXct", "POSIXt")
  attr(df$date, "tzone") <- ""
  
  # save cleaned movement dataset
  write.csv(df, paste("MovementOutput/Movement",
                      sprintf("%02d", i),'.csv', sep=""), row.names = F)
  }
  ###### 1b Obtain true values for movement params ######
  # empty matrix for results
  df <- read.csv(paste("MovementOutput/Movement",sprintf("%02d", i),'.csv', sep=""))
  HRmetrics <- data.frame(matrix(data=NA, nrow=100, ncol=13))
  colnames(HRmetrics) <- c('Replicate', 'id', 'sex',
                           'MCP', "MCPRes",
                           'KDE_hour','KDE_hourRes',
                           'KDE_day', 'KDE_dayRes','90th','80th',
                           'inTrapZone', 'crossBuffer')
  # MCP (all)
  df_all <- df[, c('id','x','y')]
  coordinates(df_all) <- ~x+y
  MCP <- mcp.area(df_all, percent = 95, unin = 'm', unout = 'm2', plotit = F)
  MCPpoly <- mcp(df_all, percent = 95, unin = 'm', unout = 'm2')
  MCP95 <- as.vector(apply(MCP[1,],2, function(x) as.numeric(x)))/10000
  # calc mean resource
  r.vals <- extract(r, MCPpoly)
  r.mean <- lapply(r.vals, FUN=mean)
  MeanRes_MCP <-do.call(rbind,r.mean)
  
  # KDE (hourly)
  df_hour <-
    df %>% 
    group_by(id) %>% 
    filter(row_number() %in% seq(60,21600,60))
  
  df_hour <- df_hour[, c('id','x','y')] # don't mess up original df for later analyses
  coordinates(df_hour) <- ~x+y
  KUD_hour <- kernelUD(df_hour, grid = xy)
  KUD95_hour <- kernel.area(KUD_hour, percent = 95, unin = 'm', unout = 'm2')
  KDE95_hour <- as.vector(apply(KUD95_hour,2, function(x) as.numeric(x))) /10000
  # extract home-range polygons, calculate MeanResource
  homepolys95_hour <- getverticeshr(KUD_hour,percent = 95, grid=xy)
  r.vals <- extract(r, homepolys95_hour)
  r.mean <- lapply(r.vals, FUN=mean)
  MeanRes_KDEhour <-do.call(rbind,r.mean)
  # establish the wide-rangers (90th & 80th percentile)
  HRA90 <- rep(0,100)
  HRA90[which(KDE95_hour[1:60]<quantile(KDE95_hour[1:60],0.9))] <- 1 # females
  HRA90[which(KDE95_hour[61:100]<quantile(KDE95_hour[61:100],0.9))+60] <- 1 # males
  HRA80 <- rep(0,100)
  HRA80[which(KDE95_hour[1:60]<quantile(KDE95_hour[1:60],0.8))] <- 1 # females
  HRA80[which(KDE95_hour[61:100]<quantile(KDE95_hour[61:100],0.8))+60] <- 1 # males
  
  # KDE (daily)
  df_day <-
    df %>% 
    group_by(id) %>% 
    filter(row_number() %in% seq(720,21600,720))
  
  x <- y <- seq(-200, 500, by=1.) # make bigger grid; resolution is the pixel size you desire 
  xy <- expand.grid(x=x,y=y)
  coordinates(xy) <- ~x+y
  gridded(xy) <- TRUE
  
  df_day <- df_day[, c('id','x','y')] # don't mess up original df for later analyses
  coordinates(df_day) <- ~x+y
  KUD_day <- kernelUD(df_day, grid = xy)
  KUD95_day <- kernel.area(KUD_day, percent = 95, unin = 'm', unout = 'm2')
  KDE95_day <- as.vector(apply(KUD95_day,2, function(x) as.numeric(x))) /10000
  homepolys95_day <- getverticeshr(KUD_day,percent = 95, grid = xy)
  r.vals <- extract(r, homepolys95_day)
  r.mean <- lapply(r.vals, FUN=mean)
  MeanRes_KDEday <- do.call(rbind,r.mean)
  
  #### Check for border effects ####
  df %>% 
    filter(x >= 60 & x <= 240, y >= 60 & y <= 240)  %>% 
    count(id) -> TrapZoners
  
  df %>% 
    filter(x <= 0 | x >= 300 | y <= 0 | y >= 300)  %>% 
    count(id) -> BufferCrossers
  
  inTrapZone <- crossBuffer <- rep(0,100)
  TrapZoners <- TrapZoners$id
  BufferCrossers <- BufferCrossers$id
  inTrapZone[TrapZoners] <- 1; crossBuffer[BufferCrossers] <- 1
  
  # fill output
  HRmetrics[,"Replicate"] <- rep(i, 100)
  HRmetrics[,"id"] <- colnames(MCP)
  HRmetrics[,'sex'] <- c(rep('F', 60), rep('M', 40))
  HRmetrics[,'MCP'] <- MCP95 
  HRmetrics[,'MCPRes'] <- MeanRes_MCP
  HRmetrics[,'KDE_hour'] <- KDE95_hour 
  HRmetrics[,'KDE_hourRes'] <- MeanRes_KDEhour
  HRmetrics[,'KDE_day'] <- KDE95_day 
  HRmetrics[,'KDE_dayRes'] <- MeanRes_KDEday
  HRmetrics[,'90th'] <- HRA90
  HRmetrics[,'80th'] <- HRA80
  HRmetrics[,'inTrapZone'] <- inTrapZone
  HRmetrics[,'crossBuffer'] <- crossBuffer
  
  write.csv(HRmetrics, paste("HRmetrics/HR",sprintf("%02d", i),'.csv', sep=""), row.names = F)
}

### Clean up cores if parallel ### 
stopCluster(cl)