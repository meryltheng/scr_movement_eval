# ====================================
# Theng et al. 2021 
# meryltheng@gmail.com
# Model evaluation
# ====================================
library(dplyr)

result<-read.csv('Estimates/resultMOD.csv')
result$RB <- (result$Nhat - result$N)/result$N 

# calculate hr estimates
result$sigmaHR <- pi * (result$sigma * sqrt(5.99) )^2 # check how to calc HR from RSF
result$hraH <- result$hraD <- rep(NA,nrow(result))
metric = 'Mean'
for (i in 1:100){
  HRmetrics <- read.csv(paste("HRmetrics/HR",sprintf("%02d", i),'.csv', sep=""))
  
  for (sex in c('F','M')){
    # A
    result$hraH[which(result$Replicate == i & result$Scenario == 'A' & result$Sex == sex)] <- rep(summary(HRmetrics$KDE_hour[HRmetrics$sex == sex])[metric]*10000,5)
    result$hraD[which(result$Replicate == i & result$Scenario == 'A' & result$Sex == sex)] <- rep(summary(HRmetrics$KDE_day[HRmetrics$sex == sex])[metric]*10000,5)
    
    # B
    resident_index90 <- HRmetrics[which(HRmetrics$X90th==1),]$id
    result$hraH[which(result$Replicate == i & result$Scenario == 'B' & result$Sex == sex)] <- rep(summary(HRmetrics$KDE_hour[HRmetrics$sex == sex & HRmetrics$id %in% resident_index90])[metric]*10000,5)
    result$hraD[which(result$Replicate == i & result$Scenario == 'B' & result$Sex == sex)] <- rep(summary(HRmetrics$KDE_day[HRmetrics$sex == sex & HRmetrics$id %in% resident_index90])[metric]*10000,5)
    
    # C
    resident_index80 <- HRmetrics[which(HRmetrics$X80th==1),]$id
    result$hraH[which(result$Replicate == i & result$Scenario == 'C' & result$Sex == sex)] <- rep(summary(HRmetrics$KDE_hour[HRmetrics$sex == sex & HRmetrics$id %in% resident_index80])[metric]*10000,5)
    result$hraD[which(result$Replicate == i & result$Scenario == 'C' & result$Sex == sex)] <- rep(summary(HRmetrics$KDE_day[HRmetrics$sex == sex & HRmetrics$id %in% resident_index80])[metric]*10000,5)
    
  }
}

result$RB.sigmaH <- (result$sigmaHR - result$hraH)/result$hraH
result$RB.sigmaD <- (result$sigmaHR - result$hraD)/result$hraD

"-------- Evaluate --------"

# Abundance N
Nevals <- result  %>%
  mutate(RB = (Nhat - N)/N )  %>% # Relative bias
  mutate(Coverage1 = ifelse(N >= Nhat - 1.96*Nhat.se & N <= Nhat + 1.96*Nhat.se,1,0)) %>% # coverage for basic model
  mutate(Coverage2 = ifelse(N >= Nhat.lo & N <= Nhat.hi,1,0)) %>% # coverage for RSF model is calculated differently (pre-calculated)
  # ----
  group_by(Model, Scenario, Sex) %>%
  mutate(CV = (sd(na.omit(Nhat)))/mean(Nhat, na.rm = T) )%>% 
  summarise(RB_L = quantile(RB, probs=c(0.25, 0.75), na.rm=T)[1], # 50% CI
            RB_U = quantile(RB, probs=c(0.25, 0.75), na.rm=T)[2], # 50% CI
            RB = mean(RB, na.rm = T), 
            CV = unique(CV), 
            Coverage1 = na.omit(length(Coverage1[Coverage1==1]))/na.omit(length(Coverage1)), # coverage for basic model
            Coverage2 = na.omit(length(Coverage2[Coverage2==1]))/na.omit(length(Coverage2)))  %>% # coverage for RSF model
  as.data.frame()

# Sigma 
Sigevals <- result  %>%
  group_by(Model, Scenario, Sex) %>%
  summarise(sigmaMu = mean(sigma, na.rm=T),
            sigmaSD = sd(na.omit(sigma)),
            CV = sd(na.omit(sigma))/mean(sigma, na.rm = T)) %>%
  as.data.frame()


# Sigma (Area used)
SigAevals <- result  %>%
  filter(Model == c('basic0','basic1'))  %>%
  mutate(RB = RB.sigmaD )  %>% # Relative bias
  mutate(Coverage1 = ifelse(hraD >= pi * ((sigma - sigma.se*1.96) * sqrt(5.99) )^2 & hraD <= pi * ((sigma + sigma.se*1.96) * sqrt(5.99) )^2,1,0)) %>% # Coverage
  mutate(Coverage2 = ifelse(hraH >= pi * ((sigma - sigma.se*1.96) * sqrt(5.99) )^2 & hraH <= pi * ((sigma + sigma.se*1.96) * sqrt(5.99) )^2,1,0)) %>% 
group_by(Model, Scenario, Sex) %>%
  mutate(CV = (sd(na.omit(sigmaHR)))/mean(sigmaHR, na.rm = T) )%>% 
  summarise(RB_L = quantile(RB, probs=c(0.25, 0.75), na.rm=T)[1], # 50% CI
            RB_U = quantile(RB, probs=c(0.25, 0.75), na.rm=T)[2], # 50% CI
            RB = mean(RB, na.rm = T), 
            CV = unique(CV),
            Coverage1 = na.omit(length(Coverage1[Coverage1==1]))/na.omit(length(Coverage1)), # coverage for basic model 
            Coverage2 = na.omit(length(Coverage2[Coverage2==1]))/na.omit(length(Coverage2)))  %>% 
  as.data.frame()

# transience 
sigma.rw <- result  %>%
  filter(Model == 'transience0')  %>%
  group_by(Scenario, Sex) %>%
  summarise(sigma = mean(sigma, na.rm=T),
    sigma.rw = mean(sigma.rw, na.rm=T)) %>%
  as.data.frame()

# Rhat
Rhat.mu <- result  %>%
  filter(Model == 'transience0')  %>%
  group_by(Scenario, Sex) %>%
  summarise(Nhat.Rhat = mean(Nhat.Rhat, na.rm=T), sigma.Rhat = mean(sigma.Rhat, na.rm=T)/10, p0.Rhat = mean(p0.Rhat, na.rm=T),
            alpha0.Rhat = mean(alpha0.Rhat, na.rm=T), alpha1.Rhat = mean(alpha1.Rhat, na.rm=T), sigma.rw.Rhat = mean(sigma.rw.Rhat, na.rm=T),
            psi.Rhat = mean(psi.Rhat, na.rm=T))  %>%
  as.data.frame()

# Effect sizes
Effectevals <- result  %>%
  mutate(D.cov.neg = ifelse(D.cov < 0,1,0), 
         p0.cov.neg = ifelse(p0.cov < 0,1,0), 
         alpha2.neg = ifelse(alpha2 < 0,1,0)) %>%
  group_by(Model, Scenario, Sex) %>%
  summarise(D.cov = mean(D.cov, na.rm = T), D.cov.neg = sum(D.cov.neg, na.rm=T),
            p0.cov = mean(p0.cov, na.rm = T), p0.cov.neg = sum(p0.cov.neg, na.rm=T),
            alpha2 = mean(alpha2, na.rm = T), alpha2.neg = sum(alpha2.neg, na.rm=T))  %>%
  as.data.frame()
