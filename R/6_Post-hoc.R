# ====================================
# Theng et al. 2021 
# meryltheng@gmail.com
# Post-hoc analyses (for Discussion)
# ====================================
library(tidyr)

"-------- Double trap density results --------"
resultA_DE <- read.csv('Estimates/resultA_DE.csv')

aggregate(p0.cov ~ Model + Sex, data = resultA_DE, mean)
aggregate(p0.cov ~ Model + Sex, data = resultA_DE, sd)
aggregate(alpha2 ~ Model + Sex, data = resultA_DE, mean)
aggregate(alpha2 ~ Model + Sex, data = resultA_DE, sd)

"-------- Occasion-specific centroid distances per individual  --------"
library(rgeos)
library(adehabitatHR)
get_dists <- function(data, xy= xy){
  data <- data[, c('occ','x','y')] 
  coordinates(data) <- ~x+y
  areaPerOcc <- kernelUD(data, grid = xy)
  polyPerOcc <- getverticeshr(areaPerOcc,percent = 95, grid = xy)
  Centroids <- gCentroid(polyPerOcc,byid=TRUE)
  dists <- sqrt(apply(apply(Centroids@coords, 2, diff)^2,1,sum))
  return(dists)
}
# test <- get_dists(df_hour[df_hour$id == 1,])

x <- y <- seq(-200, 500, by=1.) # make bigger grid to prevent error; increases compute time by a lot
xy <- expand.grid(x=x,y=y)
coordinates(xy) <- ~x+y
gridded(xy) <- TRUE

foreach(i=1:100, .combine=c, .errorhandling='stop', .packages=c(.packages())) %do% { 
  # read cleaned movement dataset
  df <- read.csv(paste("MovementOutput/Movement",sprintf("%02d", i),'.csv', sep=""), header = T)
  df_hour <-
    df %>% 
    group_by(id) %>% 
    filter(row_number() %in% seq(60,21600,60)) %>% as.data.frame()
  df_hour$occ <- rep(rep(1:5, each=72), 100) # create occasion index
  
  # apply AC dist function over each individual
  ACdists <- sapply(unique(df_hour$id), function(x)  get_dists(df_hour[df_hour$id == x,],xy= xy))
  
  HRmetrics <- read.csv(paste("HRmetrics/HR",sprintf("%02d", i),'.csv', sep=""))
  ACdists <- t(ACdists); colnames(ACdists) <- c('ACd1','ACd2','ACd3','ACd4')
  HRmetrics <- cbind(HRmetrics, ACdists)
  # write new
  write.csv(HRmetrics, paste("HRmetrics/AC/HR",sprintf("%02d", i),'.csv', sep=""), row.names = F)
}

ACdists <- list()
for (i in 1:100){
  HRmetrics <- read.csv(paste("HRmetrics/AC/HR",sprintf("%02d", i),'.csv', sep=""))
  ACdists[[i]] <- HRmetrics[,c('sex','X90th','X80th','ACd1','ACd2','ACd3','ACd4')]
}

"-------- Rs btw resource and Det/trap at diff time reso  --------"
scenario = 'A'
load('Misc/CHenviron.RData')
nDetMin <- nDetHr <- nDetOcc <- list()
for (i in 1:100){
  x <- read.csv(paste("Detections/CT",sprintf("%02d", i),'.csv', sep="")) 
  data <- x %>% 
    count(Detector) %>%
    complete(Detector = 1:100, fill = list(n = 0)) %>% as.data.frame()
  nDetMin[[i]] <- data[,2]
  
  load(file=paste("CaptHist/CH",sprintf("%02d", i), scenario,".RData", sep="")) # note that CaptHist are in daily resolution
  nDetHr[[i]] <- apply(sdh100, 3, sum)
  y <- sdh100
  y[which(y>1)] <- 1 # convert from count to binary
  nDetOcc[[i]] <- apply(y, 3, sum)
}

nDetMin_c <- do.call(c,nDetMin)
nDetHr_c <- do.call(c,nDetHr)
nDetOcc_c <- do.call(c,nDetOcc)
