# ====================================
# Theng et al. 2021 
# meryltheng@gmail.com
# SCR Markovian transience analyses
# For HPC
# ====================================

'------ For HPC processing ------'
library(foreach)
library(parallel)
library(doParallel)
library(nimble)
library(secr)

# get the input passed from the shell script
args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
  stop("At least one argument must be supplied (input file).\n", call. = FALSE)
} else {
  print(paste0("Arg input:  ", args[1]))
}

# prep input
scenario = "A"
sex = 'F'
sim = as.numeric(args[1])
fName <- paste("input/DoubleEffort/CH", sprintf("%02d", sim),scenario,".RData", sep="")
load(fName)

if (scenario == "A"){
  sdh <- sdh100
  trueN = 100
}

if (scenario == "B"){
  sdh <- sdh90
  trueN = 90
}

if (scenario == "C"){
  sdh <- sdh80
  trueN = 80
}

if (sex == 'F'){
  trueN = trueN * 0.6
} else{
  trueN = trueN * 0.4
}

# split capthist by sex for sex-specific models
y <- subset(sdh, attr(sdh, "covariates") == sex)

y[which(y>1)] <- 1 # convert from count to binary
nind <- dim(y)[1]; print(1:nind)
area = 300*300
K = dim(y)[2] # occ
X <- attr(sdh, "traps") # traplocs
J<-nrow(X)
xlim <- ylim <- c(0,300)
# Data augmentation
M <- max(3*nind, 100)
yaug<-array(0,dim=c(M,dim(y)[2],dim(y)[3]))
for(i in 1:nind){
  yaug[i,,]<- y[i,,]
}
yaug <- apply(yaug,c(3,2),cbind) # transform to fit dumb model

z<-c(rep(1,nind),rep(0,M-nind))

sst<-cbind(runif(M,xlim[1],xlim[2]),runif(M,ylim[1],ylim[2])) # starting values for s

ytot<- apply(y,c(1,2),sum) # summing detections across all traps (per id, per occ)
nind <- nrow(ytot)
ssarr<- array(NA,dim=c(M,2,K))

for(k in 1:K){
  for(i in 1:nind){
    if(sum(ytot[i,])==0) next
    sst[i,1]<- mean( X[ytot[i,]>0,1] )
    sst[i,2]<- mean( X[ytot[i,]>0,2] )
  }
  ssarr[,,k]<- sst
}
#moved.s.st<- array(rnorm(M*2*K,0,.1),dim=c(M,2,K))
#moved.s.st[,,1]<- NA
#s0<- ssarr[,,1]

data <- list(y=yaug)
const <- list (X=X/10,K=K,M=M,J=J,xlim=xlim/10,ylim=ylim/10,area=area/100)

inits <- function(){
  list (p0=runif(1),sigma.det=runif(1,1,2), s=ssarr/10, # s0=s0/10 if phi model
        z=z,sigma.rw=1) # ,phi=0.5,floater=rep(1,M)
}

parameters <-
  c("p0","alpha1","N","D","sigma.det","sigma.rw","psi","alpha0") # ,"s","phi"

code <- nimbleCode({
  p0 ~ dunif(0,1)
  alpha0<- log(p0/(1-p0))
  sigma.det ~ dunif(0, 20)
  alpha1<- 1/(2*sigma.det*sigma.det)
  psi~dunif(0,1) # Data augmentation parameter
  sigma.rw ~ dunif(0,20)
  tau<- 1/(sigma.rw*sigma.rw)
  
  for(i in 1:M){
    z[i] ~ dbern(psi) # Data augmentation parameter
    s[i,1,1]~dunif(xlim[1],xlim[2])
    s[i,2,1]~dunif(ylim[1],ylim[2])
    d[i,1:J,1]<- pow(pow(s[i,1,1]-X[1:J,1],2) + pow(s[i,2,1]-X[1:J,2],2),0.5)
    p[i,1:J,1]<- z[i]*p0*exp(- alpha1*d[i,1:J,1]*d[i,1:J,1])
    for(j in 1:J){
      y[i,j,1] ~ dbin(p[i,j,1],1)
    }
    #occasion specific s:
    for(k in 2:K){
      s[i,1,k]~T(dnorm(s[i,1,k-1],tau),xlim[1],xlim[2]) 
      s[i,2,k]~T(dnorm(s[i,2,k-1],tau),ylim[1],ylim[2])
      d[i,1:J,k]<- pow(pow(s[i,1,k]-X[1:J,1],2) + pow(s[i,2,k]-X[1:J,2],2),0.5)
      p[i,1:J,k]<- z[i]*p0*exp(- alpha1*d[i,1:J,k]*d[i,1:J,k])
      for(j in 1:J){
        y[i,j,k] ~ dbin(p[i,j,k],1) # if size/trails = 1, y needs to be binary?
      }
    }
  }
  N<-sum(z[])
  D<- N/area
})

# Prep cores in HPC cluster
slurm_ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.numeric(slurm_ncores)) {
  cores <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
} else {
  cores = detectCores()
}
cl<-makeCluster(cores)

clusterSetRNGStream(cl = cl) # set up random number stream 
registerDoParallel(cl) # register as backend for 'foreach'

# Run analysis
seeds <- 1:cores

system.time(
  res <- foreach(x = seeds, .packages="nimble",.errorhandling='remove', .inorder=FALSE) %dopar% {
    set.seed(x)
    out <- nimbleMCMC(code, data=data, constants=const,
                      inits=inits, monitors=parameters,
                      niter=22000*5, nburnin=2000*5, thin=1, nchains=1,
                      samplesAsCodaMCMC=TRUE)
  } )

# Save output
jName <- paste("output/rw", sprintf("%02d", sim),scenario,'_',sex,".RData", sep="")
save(res, file = jName)

# stop parallel
stopCluster(cl)

'------ Post-HPC data assembling ------'
# Convert to an mcmcOutput object and look at diagnostics
mclist <- coda::mcmc.list(res)
( mco <- mcmcOutput(mclist) )
diagPlot(mco, parameters)
View(summary(mco))

# Process transience model output
# ----------------------------------
resultA <- foreach(i=1:100, .combine=rbind, .errorhandling='stop', .packages=c(.packages())) %do% { 
  for (sex in c('F','M')){
    replicate = i; scenario = 'A' # edit scenario here
    load(paste("Estimates/transience/rw", sprintf("%02d", replicate),scenario,'_',sex,".RData", sep=""))
    
    if (scenario == "A"){
      trueN = 100
    }
    
    if (scenario == "B"){
      trueN = 90
    }
    
    if (scenario == "C"){
      trueN = 80
    }
    
    # empty matrix for results
    nrows=1
    output <- matrix(data=NA, nrow=nrows, ncol=37)
    colnames(output) <- c('Replicate', 'Scenario', 'Sex',# meta info
                          'N',
                          'Model',  'loglik','nlm.code',
                          'Nhat', 'Nhat.se', 'Nhat.lo', 'Nhat.hi','sigma', 'sigma.se', 'p0', 'p0.se', # param estimates
                          'D.cov', 'D.cov.se', 'p0.cov', 'p0.cov.se', # effect sizes
                          'alpha0', 'alpha0.se', 'alpha2', 'alpha2.se',
                          'Nhat.Rhat','sigma.Rhat','p0.Rhat','alpha0.Rhat','alpha1.Rhat', # rw params/effect sizes
                          'alpha1', 'alpha1.se','alpha1.Rhat', 
                          'sigma.rw', 'sigma.rw.se','sigma.rw.Rhat','psi', 'psi.se','psi.Rhat'
    )
    
    if (sex == 'F'){
      N = trueN * 0.6
    } else{
      N = trueN * 0.4
    }
    
    # meta info
    output[,'Replicate'] <- rep(replicate, nrows)
    output[,'Scenario'] <- rep(scenario, nrows)
    output[,'Sex'] <- rep(sex, nrows)
    output[,'N'] <- rep(N, nrows)
    
    # Convert to an mcmcOutput object and look at diagnostics
    mclist <- coda::mcmc.list(res)
    ( mco <- mcmcOutput(mclist) )
    summMCMC <- summary(mco)
    
    # Model output
    output[1,"Model"] <- "transience0"
    output[1,c('Nhat', 'Nhat.se','Nhat.Rhat')]  <- unlist(summMCMC['N',c(1:2,6)])
    output[1,c('Nhat.lo', 'Nhat.hi')] <- unlist(summMCMC['N',4:5])
    output[1,c('p0', 'p0.se','p0.Rhat')]  <- unlist(summMCMC['p0',c(1:2,6)])
    output[1,c('sigma', 'sigma.se','sigma.Rhat')]  <- unlist(c(summMCMC['sigma.det',1:2]*10,summMCMC['sigma.det',6]))
    output[1,c('sigma.rw', 'sigma.rw.se','sigma.rw.Rhat')]  <- unlist(c(summMCMC['sigma.rw',1:2]*10,summMCMC['sigma.rw',6]))
    output[1,c('alpha0', 'alpha0.se','alpha0.Rhat')]  <- unlist(summMCMC['alpha0',c(1:2,6)])
    output[1,c('alpha1', 'alpha1.se','alpha1.Rhat')]  <- unlist(summMCMC['alpha1',c(1:2,6)])
    output[1,c('psi', 'psi.se','psi.Rhat')]  <- unlist(summMCMC['psi',c(1:2,6)])
    
    if (sex == 'F'){
      outputF <- output
      print(paste(replicate, scenario, sex, 'done'))
    } else{
      output <- rbind(outputF,output)
    }
  }
  output
}

# save results
write.csv(resultA_doub, file='Estimates/resultA.csv', row.names = F)
