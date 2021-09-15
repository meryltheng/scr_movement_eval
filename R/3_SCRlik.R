# ====================================
# Theng et al. 2021 
# meryltheng@gmail.com
# SCR analyses (Likelihood models: standard & RSF)
# ====================================

library(secr)
library(raster)

# Prepare state-space (mask): our ss is the simulation landscape
# for SCR0 (secr)
replicate = 1
scenario = "A"
load(file=paste("CaptHist/CH",sprintf("%02d", replicate), scenario,".RData", sep=""))
mask <- make.mask(traps(sdh100), buffer = 0, spacing = 300/60, nx = 60, ny = 60, type =  "rectangular")
plot(mask, dots = FALSE, mesh = "grey", col = "white", axes = T)
plot(traps(sdh100), detpar = list(pch = 16, cex = 1), add = TRUE)

env <- read.csv('Misc/SCRland.csv', header = T) # resource landscape
envm<-t(as.matrix(env))
r <- raster(envm[nrow(envm):1,], xmn=0, xmx=nrow(envm), ymn=0, ymx=ncol(envm))
coarse.r <- aggregate(r, fact=300/60)
coarse.r [] <- scale(coarse.r [])
g <- as(coarse.r , 'SpatialGridDataFrame')

habmask <- addCovariates(mask, g) # add mask covariate (resource)

# for SCRrsf
gr <- xyFromCell(coarse.r, 1:3600) # coordinates for each pixel
zall <- coarse.r[1:3600] # covs on all cells
traplocs <- traps(sdh100)
ztrap <- coarse.r[cellFromXY(coarse.r, traplocs)] # covs at traplocs

# Prepare for parallel processing
# (Note that we ran the analyses through the HPC)
library(doParallel)
library(parallel)
library(foreach)
# [insert parallel code here]

##### Run SCR analyses ####
resultA <- foreach(i=1:100, .combine=rbind, .errorhandling='stop', .packages=c("secr")) %dopar% { 
  output <- run_scrLik(replicate = i, scenario = "A")
  output
}

write.csv(resultA, file = "Estimates/resultA.csv", row.names = F) # save results

resultB <- foreach(i=1:100, .combine=rbind, .errorhandling='stop', .packages=c("secr")) %dopar% { 
  output <- run_scrLik(replicate = i, scenario = "B")
  output
}

write.csv(resultB, file = "Estimates/resultB.csv", row.names = F)

resultC <- foreach(i=1:100, .combine=rbind, .errorhandling='stop', .packages=c("secr")) %dopar% { 
  output <- run_scrLik(replicate = i, scenario = "C")
  output
}

write.csv(resultC, file = "Estimates/resultC.csv", row.names = F)

# stop parallel
stopCluster(cl)

