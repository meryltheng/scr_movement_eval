# ====================================
# Theng et al. 2021 
# meryltheng@gmail.com
# Detection simulations
# ====================================

library(dplyr)
library(sp)
library(raster)
library(camtrapR)
library(rgeos)
library(foreach)
source("R/0_functions.R")

"-------- Prep detection input --------"
# Make CT array
# Original effort (J=100)
traps.buff <- makeTraps(xlim = c(60,240), ylim = c(60,240), trapspacing = 20, 
                        bufferdist = 0.2) # 10m if 1 unit = 50m

# Double effort (J=196)
#traps.buff <- makeTraps(xlim = c(59,241), ylim = c(59,241), trapspacing = 14, 
#                        bufferdist = 0.2) # 10m if 1 unit = 50m

# Create camera operation matrix
n_traps = length(traps.buff)
camop_no_prob <- matrix(nrow = n_traps, ncol = 30, 1) # all working
rownames(camop_no_prob) <- seq(1:n_traps)
colnames(camop_no_prob) <- as.character(as.Date(seq(1:30), origin = "1969-12-31"))

# Prep CTtable (trap id, x, y, covariate values)
X <- coordinates(traps.buff)
env <- read.csv('Misc/SCRland.csv', header = T) # resource landscape
envm<-t(as.matrix(env))
r <- raster(envm[nrow(envm):1,], xmn=0, xmx=nrow(envm), ymn=0, ymx=ncol(envm))
r.vals <- extract(r, traps.buff) # calculate mean resource value within CT detection zone
r.mean <- lapply(r.vals, FUN=mean)
r.mean<-do.call(rbind,r.mean)
# assemble CTtable
camtraps <- data.frame(overlap = 1:n_traps, x = X[,1], y= X[,2], cov1 = round(scale(r.mean)/0.5)*0.5) # resource cov is scaled and discretised to speed up computation in secr (only the case for detection covs)
colnames(camtraps) <- c('Detector','X','Y','cov1')

"-------- Run detection simulations --------"
foreach(i=1:100, .combine=c, .errorhandling='stop', .packages=c(.packages())) %do% { 
  # read cleaned movement dataset
  df <- read.csv(paste("MovementOutput/Movement",sprintf("%02d", i),'.csv', sep=""), header = T)
  df$date <- as.POSIXct(df$date, format ='%Y-%m-%d %H:%M:%S')
  
  ###### 3 Detection sims ######
  x <- DetectionGenerator(df=df, X = traps.buff)
  
  # clean up
  smurfCH <- as.data.frame(x)# tibble to df
  smurfCH$Species <- "Smurf"
  colnames(smurfCH) <- c('ID','x','y', 'sex','date', 'Detector', 'Species')
  
  # save cleaned detections dataset
  write.csv(smurfCH, paste("Detections/CT",sprintf("%02d", i),'.csv', sep=""), row.names = F)
  
  print(paste('Detections for replicate',sprintf("%02d", i),'complete'))
  
}
 

"-------- Create CaptHist --------"
foreach(i=1:100, .combine=c, .errorhandling='stop', .packages=c(.packages())) %do% {  
  # read detection dataset
  smurfCH <- read.csv(paste("Detections/CT",sprintf("%02d", i),'.csv', sep=""), header = T)
  
  ###### 4 Create capture histories ######
  K = 6 # occ length
  
  # 90th and 80th percentile reference for scenarios B and C
  HRmetrics <- read.csv(paste("HRmetrics/HR",sprintf("%02d", i),'.csv', sep=""))
  
  # A
  sdh100 <- spatialDetectionHistory(recordTableIndividual = smurfCH,
                                    species = "Smurf",
                                    camOp = camop_no_prob,
                                    CTtable = camtraps,
                                    output = "count",
                                    stationCol = "Detector",
                                    stationCovariateCols = "cov1",
                                    speciesCol = "Species",
                                    Xcol = "X",
                                    Ycol = "Y",
                                    individualCol = "ID",
                                    individualCovariateCols = "sex",
                                    recordDateTimeCol = "date",
                                    recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                    occasionLength = K,
                                    day1 = "1970-01-01",
                                    includeEffort = TRUE
  )
  
  save(sdh100, file=paste("CaptHist/CH",sprintf("%02d", i),"A.RData", sep=""))
  # B
  resident_index90 <- HRmetrics[which(HRmetrics$X90th==1),]$id
  smurfCH90 <- smurfCH[smurfCH$ID %in% resident_index90,]
  sdh90 <- spatialDetectionHistory(recordTableIndividual = smurfCH90,
                                   species = "Smurf",
                                   camOp = camop_no_prob,
                                   CTtable = camtraps,
                                   output = "count",
                                   stationCol = "Detector",
                                   stationCovariateCols = "cov1",
                                   speciesCol = "Species",
                                   Xcol = "X",
                                   Ycol = "Y",
                                   individualCol = "ID",
                                   individualCovariateCols = "sex",
                                   recordDateTimeCol = "date",
                                   recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                   occasionLength = K,
                                   day1 = "1970-01-01",
                                   includeEffort = TRUE
  )
  
  save(sdh90, file=paste("CaptHist/CH",sprintf("%02d", i),"B.RData", sep=""))
  
  # C
  resident_index80 <- HRmetrics[which(HRmetrics$X80th==1),]$id
  smurfCH80 <- smurfCH[smurfCH$ID %in% resident_index80,]
  sdh80 <- spatialDetectionHistory(recordTableIndividual = smurfCH80,
                                   species = "Smurf",
                                   camOp = camop_no_prob,
                                   CTtable = camtraps,
                                   output = "count",
                                   stationCol = "Detector",
                                   stationCovariateCols = "cov1",
                                   speciesCol = "Species",
                                   Xcol = "X",
                                   Ycol = "Y",
                                   individualCol = "ID",
                                   individualCovariateCols = "sex",
                                   recordDateTimeCol = "date",
                                   recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                   occasionLength = K,
                                   day1 = "1970-01-01",
                                   includeEffort = TRUE
  )
  
  save(sdh80, file=paste("CaptHist/CH",sprintf("%02d", i),"C.RData", sep=""))
  
  print(paste('CaptHists for replicate',sprintf("%02d", i),'complete'))
  
}

"-------- Summarise CaptHist --------"
CHsumm<- foreach(i=1:100, .combine=rbind, .errorhandling='stop', .packages=c(.packages())) %do% { 
  output <- summarise_CH(replicate = i, scenario = "A", type = 'binary')
  output
}

resultCap <- CHsumm %>%
  group_by(Sex, Scenario) %>%
  summarise(nDet =  mean(nDet, na.rm = T), nCap = mean(nCap, na.rm = T), nTrapVis = mean(nTrapVis, na.rm = T), nOnce = mean(nOnce, na.rm = T), 
            avg.caps = mean(avg.caps, na.rm = T), avg.nlocs = mean(avg.nlocs, na.rm = T), nTrapVis.occ = mean(nTrapVis.occ, na.rm = T))
resultCapSD <- CHsumm %>%
  group_by(Sex, Scenario) %>%
  summarise(nDet =  sd(nDet, na.rm = T), nCap = sd(nCap, na.rm = T), nTrapVis = sd(nTrapVis, na.rm = T), nOnce = sd(nOnce, na.rm = T), 
            avg.caps = sd(avg.caps, na.rm = T), avg.nlocs = sd(avg.nlocs, na.rm = T), nTrapVis.occ = sd(nTrapVis.occ, na.rm = T))

