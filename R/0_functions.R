# ====================================
# Theng et al. 2021 
# meryltheng@gmail.com
# Functions
# ====================================

"DETECTION SIMULATIONS"
makeTraps <- function(xlim = c(60,240), ylim = c(60,240), trapspacing = 20, bufferdist = 0.2) { # 10m if 1 unit = 50m
  grid <- expand.grid(x=seq(xlim[1],xlim[2],trapspacing),y=seq(ylim[1],ylim[2],trapspacing)) 
  coordinates(grid) <- ~x+y
  grid <- SpatialPoints(grid)
  traps.buff <- gBuffer(grid, width = bufferdist ,quadsegs=round(bufferdist/10,0))
  traps.buff <- disaggregate(traps.buff)
  return(traps.buff)
}
DetectionGenerator <- function(df=NA, X = traps.buff){
  df$overlap <- rep(NA, nrow(df))
  xline <- list()
  for (i in 1:nrow(df)){
    x <- rbind( df[i,c('x','y')], df[i-1,c('x','y')])
    xline[[i]] <- Lines(list(Line(coordinates(x))), ID=paste0(i))
  }
  df.lines <- SpatialLinesDataFrame(SpatialLines(xline, CRS(as.character(NA))), 
                                    df, match.ID = FALSE)
  df %>% 
    add_rownames() %>%
    group_by(id) %>% 
    filter(row_number()==1) %>%
    `[[`("rowname") %>%
    as.numeric() -> index # create index for first line from every id
  
  for(k in 1:length(index)){  # remove first line from every id
    df.lines@lines[[index[k]]]@Lines[[1]] <- NULL
  }
  
  overlap <- over(df.lines, traps.buff, returnList = F)
  df.lines@data[["overlap"]] <- overlap
  overlap.df <- df.lines[!is.na(df.lines@data[["overlap"]]),]
  x=overlap.df
  x <- as.data.frame(x)
  return(x)
}

DetectionGeneratorLoRes <- function(df=NA, X = traps.buff){
  breaks <- seq(60,nrow(df),60)
  xline <- list(); k = 0
  for (i in breaks){
    k = k + 1
    if (i-60 == 0){
      xline[[k]] <- Lines(list(Line(coordinates(df[1:i,c('x','y')]))), ID=paste0(i))
    } else {    
      if (i %in% (seq(21600,21600*100,21600)+60)){
        xline[[k]] <- Lines(list(Line(coordinates(df[(i-59):i,c('x','y')]))), ID=paste0(i))
      } else
        xline[[k]] <- Lines(list(Line(coordinates(df[(i-60):i,c('x','y')]))), ID=paste0(i))}
  }
  
  lineIndex <- c(1,seq(61,nrow(df),60))
  
  df.lines <- SpatialLinesDataFrame(SpatialLines(xline, CRS(as.character(NA))), 
                                    df[lineIndex,], match.ID = FALSE)
  
  index <- seq(1, nrow(df)/60, 21600/60)
  
  overlap <- over(df.lines, traps.buff, returnList = F)
  df.lines@data[["overlap"]] <- overlap
  overlap.df <- df.lines[!is.na(df.lines@data[["overlap"]]),]
  x=overlap.df
  x <- as.data.frame(x)
  return(x)
}

"SUMMARISE CAPTURE HISTORIES"
summarise_CH <- function(replicate = 1, scenario = "A", type = "binary"){
  
  load(file=paste("CaptHist/CH",sprintf("%02d", replicate), scenario,".RData", sep=""))
  
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
  
  if (type == 'binary'){
    sdh[which(sdh>1)] <- 1 # convert from count to binary 
  }
  
  for (sex in c('F','M')){
    # split capthist by sex for sex-specific models
    y <- subset(sdh, attr(sdh, "covariates") == sex)
    
    # empty matrix for results
    output <- matrix(data=NA, nrow=1, ncol=12)
    colnames(output) <- c('Replicate', 'Scenario', 'Sex',# meta info
                          'N', 'nCap', 'nDet', 'nTrapVis', 'nOnce','nRecap','avg.caps', 'avg.nlocs','nTrapVis.occ') # summary stats
    output <- as.data.frame(output)
    
    if (sex == 'F'){
      N = trueN * 0.6
    } else{
      N = trueN * 0.4
    }
    
    # meta info
    output[,'Replicate'] <- replicate
    output[,'Scenario'] <- scenario
    output[,'Sex'] <- sex
    
    # summary stats
    output[,'N'] <- N
    output[,'nCap'] <- dim(y)[1]
    output[,'nDet'] <- sum(y)
    
    yObs <- yCount <- apply(y, c(1,3), sum) # Add up number of captures/detections across occasions
    nCap.ind = rowSums(yCount) # number of independent captures per individual (1day)
    output[,'nOnce'] = length(nCap.ind[nCap.ind==1])
    output[,'nRecap'] = length(nCap.ind[nCap.ind>1])
    output[,'avg.caps'] = mean(nCap.ind) # average captures across captured indivs
    output[,'nTrapVis'] = length(colSums(yCount)[colSums(yCount) > 0]) # no traps visit
    
    yObs[yObs > 0] = 1 # Convert this to binary and count ntraps per individual
    nTraps.ind = rowSums(yObs)
    output[,'avg.nlocs'] = mean(nTraps.ind)
    y[y > 0] = 1 # Convert CH to binary
    yOcc <- apply(y, c(1,2), sum) #Add up number of captures/detections across traps
    output[,'nTrapVis.occ'] = mean(rowMeans(yOcc)) # no traps visited per ind per occ
    
    if (sex == 'F'){
      outputF <- output
    } else{
      output <- rbind(outputF,output)
    }
  }
  return(output)
}

"SCR RSF"
# obtained from Royle et al. 2013
e2dist <- function (x, y){
  i <- sort(rep(1:nrow(y), nrow(x)))
  dvec <- sqrt((x[, 1] - y[i, 1])^2 + (x[, 2] - y[i, 2])^2)
  matrix(dvec, nrow = nrow(x), ncol = nrow(y), byrow = F)
}

spatial.plot <- function(x,y){
  nc<-as.numeric(cut(y,20))
  plot(x,pch=" ")
  points(x,pch=20,col=topo.colors(20)[nc],cex=2)
  ###image.scale(y,col=topo.colors(20))
}

intlik3rsf.v2 <-function(start=NULL,y=y,K=NULL,X=traplocs,ztrap,G,ntel=NULL,zall=NULL,stel=NULL){
  # start = vector of length 4 = starting values
  ###
  ###
  ### ORDER OF STARTING VALUES
  ### (intercept for cloglog-scale detection, log(sigma), coefficient on z(x), log(n0))
  ### where N = n_observed + n0  i.e., n0 = number of uncaptured individuals
  ### to transform back take , e.g., nrow(y)+exp(tmp1$estimate[4]) = Nhat
  ###
  ###
  #y
  #K
  #X
  # ztrap = covariate value at trap locations
  # zall = all covariate values for all nG pixels
  # ntel = nguys x nG matrix of telemetry fixes in each nG pixels
  # stel = home range center of telemetered individuals, IF you wish to estimate it. Not necessar
  
  nG<-nrow(G) # G: centre coordinates of landscape pixels
  D<- e2dist(X,G) # distance matrix from traplocs to G
  
  alpha0<-start[1]
  sigma<- exp(start[2])
  alpha2<- start[3]
  n0<- exp(start[4])
  a0<- 1
  
  if(!is.null(y)){
    loglam<- alpha0 -(1/(2*sigma*sigma))*D*D + alpha2*ztrap # ztrap recycled over nG
    probcap<- 1-exp(-exp(loglam))
    #probcap<- (exp(theta0)/(1+exp(theta0)))*exp(-theta1*D*D)
    Pm<-matrix(NA,nrow=nrow(probcap),ncol=ncol(probcap))
    ymat<-y
    ymat<-rbind(y,rep(0,ncol(y)))
    lik.marg<-rep(NA,nrow(ymat))
    for(i in 1:nrow(ymat)){
      Pm[1:length(Pm)]<- (dbinom(rep(ymat[i,],nG),rep(K,nG),probcap[1:length(Pm)],log=TRUE))
      lik.cond<- exp(colSums(Pm))
      lik.marg[i]<- sum( lik.cond*(1/nG) )
    }
    nv<-c(rep(1,length(lik.marg)-1),n0)
    part1<- lgamma(nrow(y)+n0+1) - lgamma(n0+1)
    part2<- sum(nv*log(lik.marg))
    out<- -1*(part1+ part2)
  }
  else{
    out<-0
  }
  if(!is.null(ntel) & !is.null(stel) ){
    # this is a tough calculation here
    D2<- e2dist(stel,G)^2
    #lam is now nGxnG!
    lam<- t(exp(a0 - (1/(2*sigma*sigma))*t(D2)+ alpha2*zall)) # recycle zall over all ntel guys
    denom<-rowSums(lam)
    probs<- lam/denom # each column is the probs for a guy at column [j]
    
    tel.loglik<- -1*sum( ntel*log(probs) )
    out<- out + tel.loglik
  }
  
  if(!is.null(ntel) & is.null(stel) ){
    # this is a tough calculation here
    D2<- e2dist(G,G)^2
    #lam is now nGxnG!
    lam<- t(exp(a0 - (1/(2*sigma*sigma))*t(D2)+ alpha2*zall)) # recycle zall over all ntel guys
    denom<-rowSums(lam)
    probs<- t(lam/denom) # each column is the probs for a guy at column [j]
    marg<- as.vector(rowSums(exp(ntel%*%log(probs))/nG ))
    
    tel.loglik<- -1*sum(log(marg))
    
    out<- out + tel.loglik
  }
  
  out
}
intlik3rsf.v3 <-function(start=NULL,y=y,K=NULL,X=traplocs,ztrap,G,ntel=NULL,zall=NULL,stel=NULL){
  #this version of the code handles a covariate on log(Density). This is starting value 5
  #start = vector of length 5 = starting values
  #y = nind x ntraps encounter matrix
  #K = how many samples?
  #X = trap locations
  #ztrap = covariate value at trap locations
  #zall = all covariate values for all nG pixels
  #ntel = nguys x nG matrix of telemetry fixes in each nG pixels
  #stel = home range center of telemetered individuals, IF you wish to estimate it. Not necessary
  
  nG<-nrow(G)
  D<- e2dist(X,G)
  
  alpha0<-start[1]
  sigma<- exp(start[2])
  alpha2<- start[3]
  n0<- exp(start[4])
  beta<- start[5]
  a0<- 1
  if(!is.null(zall)){
    psi<- exp(beta*zall)
    psi<-psi/sum(psi)
  }
  else{
    psi<-rep(1/nG,nG)
  }
  if(!is.null(y)){
    loglam<- alpha0 -(1/(2*sigma*sigma))*D*D + alpha2*ztrap # ztrap recycled over nG
    probcap<- 1-exp(-exp(loglam))
    #probcap<- (exp(theta0)/(1+exp(theta0)))*exp(-theta1*D*D)
    Pm<-matrix(NA,nrow=nrow(probcap),ncol=ncol(probcap))
    ymat<-y
    ymat<-rbind(y,rep(0,ncol(y)))
    lik.marg<-rep(NA,nrow(ymat))
    
    for(i in 1:nrow(ymat)){
      Pm[1:length(Pm)]<- (dbinom(rep(ymat[i,],nG),rep(K,nG),probcap[1:length(Pm)],log=TRUE))
      lik.cond<- exp(colSums(Pm))
      lik.marg[i]<- sum( lik.cond*psi ) # this is diff from v2
    }
    nv<-c(rep(1,length(lik.marg)-1),n0)
    part1<- lgamma(nrow(y)+n0+1) - lgamma(n0+1)
    part2<- sum(nv*log(lik.marg))
    out<- -1*(part1+ part2)
  }
  else{
    out<-0
  }
  if(!is.null(ntel) & !is.null(stel) ){
    # this is a tough calculation here
    D2<- e2dist(stel,G)^2
    #lam is now nGxnG!
    lam<- t(exp(a0 - (1/(2*sigma*sigma))*t(D2)+ alpha2*zall)) # recycle zall over all ntel guys
    denom<-rowSums(lam)
    probs<- lam/denom # each column is the probs for a guy at column [j]
    
    tel.loglik<- -1*sum( ntel*log(probs) )
    out<- out + tel.loglik
  }
  
  if(!is.null(ntel) & is.null(stel) ){
    # this is a tough calculation here
    D2<- e2dist(G,G)^2
    #lam is now nGxnG!
    lam<- t(exp(a0 - (1/(2*sigma*sigma))*t(D2)+ alpha2*zall)) # recycle zall over all ntel guys
    denom<-rowSums(lam)
    probs<- t(lam/denom) # each column is the probs for a guy at column [j]
    temp<-exp(ntel%*%log(probs)) # Ntel x nG matrix
    
    marg<- as.vector(rowSums( temp*psi ))
    
    tel.loglik<- -1*sum(log(marg))
    out<- out + tel.loglik
  }
  
  out
}

"RUN LIKELIHOOD MODELS SCR0 and RSF"
run_scrLik <- function(replicate = 1, scenario = "A"){
  
  load(file=paste("input/CH",sprintf("%02d", replicate), scenario,".RData", sep=""))
  
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
  
  
  for (sex in c('F','M')){
    # split capthist by sex for sex-specific models
    y <- subset(sdh, covariates(sdh)$sex == sex)
    y[which(y>1)] <- 1 # convert from count to binary
    detector(traps(y)) <- 'proximity'
    
    # empty matrix for results
    nrows=4
    output <- matrix(data=NA, nrow=nrows, ncol=23)
    colnames(output) <- c('Replicate', 'Scenario', 'Sex',# meta info
                          'N',
                          'Model',  'loglik','nlm.code',
                          'Nhat', 'Nhat.se', 'Nhat.lo', 'Nhat.hi','sigma', 'sigma.se', 'p0', 'p0.se', # param estimates
                          'D.cov', 'D.cov.se', 'p0.cov', 'p0.cov.se', # effect sizes
                          'alpha0', 'alpha0.se', 'alpha2', 'alpha2.se') # rsf params/effect sizes
    
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
    
    ###### 5 SCR ANALYSES ######
    ###### 5a SCR0 ######
    # null model
    output[1,"Model"] <- "basic0"
    M1 <- secr.fit(y, mask = mask, trace = FALSE,)
    output[1,c('Nhat', 'Nhat.se')] <- unlist(region.N(M1)['E.N', 1:2])
    output[1,c('sigma', 'sigma.se')] <- unlist(summary(M1)$predicted['sigma', 2:3])
    output[1,c('p0', 'p0.se')] <- unlist(summary(M1)$predicted['g0', 2:3])
    output[1,"loglik"] <- as.numeric(logLik(M1))
    
    # Cov: g0 ~ z, sigma ~ 1, D ~ z 
    output[2,"Model"] <- "basic1"
    M2 <- secr.fit(y, mask = habmask, trace = FALSE,
                   model = list(g0 ~ cov1, D ~ layer))
    output[2,c('Nhat', 'Nhat.se')] <- unlist(region.N(M2)['E.N', 1:2])
    output[2,c('sigma', 'sigma.se')] <- unlist( exp(coef(M2)['sigma', 1:2]) ) #unlist(summary(M2)$predicted['sigma', 2:3])
    output[2,c('p0', 'p0.se')] <- c(plogis(coef(M2)['g0', 1]),
                                    plogis(coef(M2)['g0', 2])) #unlist(summary(M2)$predicted['g0', 2:3])
    output[2,c('p0.cov', 'p0.cov.se')] <- unlist(coef(M2)['g0.cov1',1:2])
    output[2,c('D.cov', 'D.cov.se')] <- unlist(coef(M2)['D.layer',1:2])
    output[2,"loglik"] <- as.numeric(logLik(M2))
    
    ###### 5b SCRrsf ######
    yCount <- apply(y, c(1,3), sum) # Add up number of captures/detections across occasions
    # basic RSF model
    # with cloglog(p(x)) = covariate
    output[3,"Model"] <- "rsf0"
    M3 <-nlm(intlik3rsf.v2,c(-3,0.28,0,4),y=yCount,K=dim(y)[2],X=traplocs/100,ztrap=ztrap,G=gr/100, hessian = TRUE) # change length units to km
    H <- M3$hessian
    seEsts3 <- sqrt(diag(solve(H)))
    output[3,'Nhat'] <- nrow(yCount)+exp(M3$estimate[4])
    output[3,'Nhat.lo'] <- nrow(yCount)+exp(M3$estimate[4] - 1.96*seEsts3[4])
    output[3,'Nhat.hi'] <- nrow(yCount)+exp(M3$estimate[4] + 1.96*seEsts3[4])
    output[3,c('sigma','sigma.se')] <- c(exp(M3$estimate[2]), exp(seEsts3[2]))*100 # convert back
    output[3,c('alpha0', 'alpha0.se')] <- c(M3[["estimate"]][1], seEsts3[1])
    output[3,c('alpha2', 'alpha2.se')] <- c(M3[["estimate"]][3], seEsts3[3])
    output[3,"loglik"] <- M3[["minimum"]]
    output[3,"nlm.code"] <- M3[["code"]]
    
    # inhomgenous density
    # Fits SCR model with isotropic Gaussian encounter model log(D(x)) = covariate
    # AND with cloglog(p(x)) = covariate
    output[4,"Model"] <- "rsf1"
    M4 <- nlm(intlik3rsf.v3,c(-3,log(3),0,4,0),y=yCount,K=dim(y)[2],X=traplocs/100,ztrap=ztrap,G=gr/100,
              zall=zall,hessian=TRUE)
    H <- M4$hessian
    seEsts4 <- sqrt(diag(solve(H)))
    output[4,'Nhat'] <- nrow(yCount)+exp(M4$estimate[4])
    output[4,'Nhat.lo'] <- nrow(yCount)+exp(M4$estimate[4] - 1.96*seEsts4[4])
    output[4,'Nhat.hi'] <- nrow(yCount)+exp(M4$estimate[4] + 1.96*seEsts4[4])
    output[4,c('sigma','sigma.se')] <- c(exp(M4$estimate[2]), exp(seEsts4[2]))*100 # convert back
    output[4,c('alpha0', 'alpha0.se')] <- c(M4[["estimate"]][1], seEsts4[1])
    output[4,c('alpha2', 'alpha2.se')] <- c(M4[["estimate"]][3], seEsts4[3])
    output[4,c('D.cov', 'D.cov.se')] <- c(M4[["estimate"]][5], seEsts4[5])
    output[4,"loglik"] <- M4[["minimum"]]
    output[4,"nlm.code"] <- M4[["code"]]
    
    if (sex == 'F'){
      outputF <- output
      print(paste(replicate, scenario, sex, 'done'))
    } else{
      output <- rbind(outputF,output)
    }
  }
  return(output)
}

