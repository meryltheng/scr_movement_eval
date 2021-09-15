# ====================================
# Theng et al. 2021 
# meryltheng@gmail.com
# Plotting
# ====================================

library(ggplot2)
library(scales)

'------ Fig 1 ------'
HRmetrics <- read.csv(paste("HRmetrics/HR",sprintf("%02d", 1),'.csv', sep=""))
tiff(filename="Fig1.tif",height=6500,width=2800,units="px",res=800,compression="lzw")
par(mfrow=c(3,1) , mar = c(4.1, 4.5, 2, 1)) #, oma = c(4,3,3,0))
# Fig. 1a
# females
hist(HRmetrics[HRmetrics$sex=='F','KDE_hour']*10000, breaks = seq(0,40000,500),
     main = '', xlab=expression(paste(italic(A[0.95]), plain(' (hourly)'))), ylab='count', 
     xlim = c(0,50000), ylim = c(0,20),
     col = rgb(0.55,0,0,0.05), border=rgb(1,0,0,0.0),
     cex.axis = 1.5, cex.lab = 1.5)
for (i in 2:100){
  HRmetrics <- read.csv(paste("HRmetrics/HR",sprintf("%02d", i),'.csv', sep=""))
  hist(HRmetrics[HRmetrics$sex=='F','KDE_hour']*10000, breaks = seq(0,40000,500),
       main = '', xlab=' ',ylab=' ', xlim = c(0,50000), ylim = c(0,20),
       col = rgb(0.55,0,0,0.05), border=rgb(1,0,0,0.0), add=T)
}
# males
for (i in 1:100){
  HRmetrics <- read.csv(paste("HRmetrics/HR",sprintf("%02d", i),'.csv', sep=""))
  hist(HRmetrics[HRmetrics$sex=='M','KDE_hour']*10000, breaks = seq(0,50000,500),
       main = '', xlab=' ',ylab=' ', xlim = c(0,50000), ylim = c(0,20),
       col = rgb(0,0,0.55,0.05), border=rgb(1,0,0,0.0), add=T)
}
box(bty="l")
legend("topright", legend=c('Female','Male'), title='Sex',
       col=c('darkred','darkblue'), pch = 15, cex=1.5, bty='n')

# Fig. 1b
ResourceArea <- lapply(1:100, function(x){
  HRmetrics <- read.csv(paste("HRmetrics/HR",sprintf("%02d", x),'.csv', sep=""))
  return(HRmetrics[,c('sex','KDE_hour','KDE_hourRes','KDE_day','KDE_dayRes')])
})
ResourceArea_c <- do.call(rbind, ResourceArea)
#log10(range(ResourceArea_c[,2]))
plot(log10(ResourceArea_c[ResourceArea_c$sex=='F',3]),log10(ResourceArea_c[ResourceArea_c$sex=='F',2]), 
     pch = 16, col = alpha('darkred',.1), cex = 0.8, ylim = c(-1.5,0.68), bty='n',
     ylab=expression(paste(plain('log'[10]),plain('( '),italic(A[0.95]), plain(' (hourly) )'))), 
     xlab=expression(paste(plain('log'[10]),plain('( Mean resource )'))),
     cex.axis = 1.5, cex.lab = 1.5)
points(log10(ResourceArea_c[ResourceArea_c$sex=='M',3]),log10(ResourceArea_c[ResourceArea_c$sex=='M',2]), 
       pch = 16, col = alpha('darkblue',.1), cex = 0.8)
box(bty="l")

# Fig. 1c
library(adehabitatHR)
TAC <- function(n_steps = 21600, burnin = 0,interval = 720, iteration = 1){
  skips <- seq(interval, (n_steps - burnin), by=interval)
  df <- read.csv(paste('MovementOutput/Movement',sprintf("%02d", iteration),'.csv', sep=''), header = T)
  output <- matrix(NA, length(skips), length(unique(df$id)))
  i = 1
  for (q in skips){ 
    df.trunc <-
      df %>% 
      group_by(id) %>% 
      filter(row_number() %in% seq(60,q,60)) # hourly locs
    df.trunc <- df.trunc[, c('id','x','y')]
    coordinates(df.trunc) <- ~x+y
    MCP <- mcp.area(df.trunc, percent = 95, unin = 'm', unout = 'm2', plotit = F)
    # calculate and store home-range estimates
    MCP95v <- as.vector(apply(MCP[1,],2, function(x) as.numeric(x)))
    output[i,] <- MCP95v
    i = i + 1
  }
  return(output)
}

# load('Data/TAC.RData')
x <- sample(1:100, 1) # randomly sample one rep, 77
matplot(TACall[[x]][,1:60], type = "l", main = "", ylab = "TAC", xlab = "Days", col = 'darkred', 
        bty='n', ylim=c(0,21000),
        cex.axis = 1.5, cex.lab = 1.5)
matplot(TACall[[x]][,61:100], type = "l", main = "", ylab = "TAC", xlab = "Days", col = 'darkblue', add=T)
box(bty="l")
dev.off()

'------ Fig 2 ------'
# functions
colfunc <- colorRampPalette(c("white","gray30"))
plotTraject <- function(df = df, plotit = X){
  k = 0
  n_foragers = length(unique(df$id))
  for(i in sort(unique(df$id))){ # lines
    k=k+1
    if(plotit[k] == 1) { # or i
      lines(x=df$x[df$id==i],
            y=df$y[df$id==i],
            lwd = 1.5,
            col = alpha(rainbow(n_foragers)[k],0.25))}
  } 
}

df <- read.csv('MovementInput/Females/Tracks01.csv') # raw movement data
df <- read.csv('MovementInput/Males/Tracks01.csv') # raw movement data
env <- read.csv('Misc/SCRland.csv', header = T) # resource landscape

df$x<- df$x - 50 # correct for empty buffer
df$y<- df$y - 50
df %>% 
  group_by(id) %>% 
  filter(row_number() %in% 21601:43200) -> df

tiff(filename="Fig2.tif",height=2800,width=5200,units="px",res=800,compression="lzw")
par(mfrow=c(1,2), mar = c(2.1,2.1,1,1))
image(1:nrow(env), 1:ncol(env), as.matrix(env), col = colfunc(10), 
      xlim = c(-10,310), ylim = c(-10,310), xlab = '', ylab = '')
plotTraject(df = df, plotit = rep(1,100)) # all
points(X,pch=3)
dev.off()

'------ Fig 3 ------'
library(ggpubr)
scrcols <- c("#dfc27d","#a6611a","#80cdc1","#018571","#f5f5f5")

# Abundance N (RB)
NPlot <- ggviolin(result, x = "Scenario", y = "RB", fill = "Model", palette = scrcols, facet.by = 'Sex',
                  size = 0.3, add = "mean_se", width = 0.5, add.params = list(size=0.1)) + 
  ggtitle('') + xlab('Scenario') + ylab(expression(paste(italic(hat(N))))) + #ylim(-0.6,1.2) +
  geom_hline(yintercept=0, linetype="dashed", color = "black",size=0.2) + 
  #rremove('xlab') + rremove('x.text') + 
  theme(strip.background = element_blank(),
        strip.text.x = element_text(size=12))
tiff(filename="Fig3.tif",height=2500,width=5200,units="px",res=800,compression="lzw")
NPlot
dev.off()

'------ Fig 4 ------'
# Sigma.det (RB)
SigmaPlot <- ggviolin(result, x = "Scenario", y = "sigma", fill = "Model", palette = scrcols, facet.by = 'Sex',
                       size = 0.3, add = "mean_se", width = 0.5, add.params = list(size=0.1)) + 
  ggtitle('') + xlab('Scenario') + ylab(expression(paste(italic(hat(sigma)[det])))) + #ylim(-0.6,1.2) +
  #rremove('legend') + rremove('x.text') +
  rremove('xlab') + 
  theme(strip.background = element_blank(), 
        strip.text.x = element_text(size=11))

SigmahPlot <- ggviolin(result[result$Model== c('basic0','basic1'),], x = "Scenario", y = "RB.sigmaH", fill = "Model", palette = scrcols, facet.by = 'Sex',
                       size = 0.3, add = "mean_se", width = 0.5, add.params = list(size=0.1)) + 
  ggtitle('') + xlab(' ') + ylab(expression(paste(italic(hat(A)[0.95]), plain(' (hourly)')))) + #ylim(-0.6,1.2) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.2) + 
  #rremove('x.text') + 
  rremove('legend') + #rremove('xlab') + 
  theme(strip.background = element_blank(), 
        strip.text.x = element_text(size=11))

SigmadPlot <- ggviolin(result[result$Model== c('basic0','basic1'),], x = "Scenario", y = "RB.sigmaD", fill = "Model", palette = scrcols, facet.by = 'Sex',
                       size = 0.3, add = "mean_se", width = 0.5, add.params = list(size=0.1)) + 
  ggtitle('') + xlab(' ') + ylab(expression(paste(italic(hat(A)[0.95]), plain(' (daily)'))))  + #ylim(-0.6,1.2) +
  geom_hline(yintercept=0, linetype="dashed",color = "black", size=0.2) + 
  rremove('legend')  + #rremove('xlab') + 
  theme(strip.background = element_blank(), 
        strip.text.x = element_text(size=11))

fig4 <- ggarrange(SigmaPlot,
          ggarrange(SigmahPlot, SigmadPlot,ncol = 2, labels = c('(b)','')), # 
          nrow = 2, heights = c(1,1), labels = '(a)', common.legend=TRUE) #
tiff(filename="Fig4.tif",height=3700,width=5200,units="px",res=800,compression="lzw")
annotate_figure(fig4, fig.lab = "Scenario", fig.lab.pos = 'bottom', fig.lab.size = 12, fig.lab.face = 'bold')
dev.off()

'------ Fig 5 ------'
# Effect sizes
# Dcov
Dcov <- ggviolin(result, x = "Scenario", y = "D.cov", fill = "Model", facet.by = 'Sex',#scrcols[c(2,4)], palette = scrcols[c(2,4)]
        add = "mean_se", width = 0.5, add.params = list(size=0.3)) + 
  ggtitle('') + xlab(NULL) + ylab(expression(paste(italic(D),plain('.cov')))) + #ylim(-0.6,1.2) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.2) + 
  scale_fill_manual(values = scrcols,
                    labels = c("basic0","basic1", "rsf0", "rsf1","transience0"), 
                    drop = FALSE) +
  rremove('xlab') + rremove('x.text') +
  theme(strip.background = element_blank(), 
        strip.text.x = element_text(size=11)) + labs(x=NULL)

# p0.cov
p0.cov <- ggviolin(result, x = "Scenario", y = "p0.cov", fill = "Model", palette = scrcols[2], facet.by = 'Sex',
         add = "mean_se", width = 0.3, add.params = list(size=0.3)) + 
  ggtitle('') + xlab(NULL) + ylab(expression(paste(italic(p)[0],plain('.cov')))) + #ylim(-0.6,1.2) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.2) + 
  rremove('xlab') + rremove('x.text') +
  theme(strip.background = element_blank(), 
        strip.text.x = element_blank())

# alpha2
alpha2 <- ggviolin(result, x = "Scenario", y = "alpha2", fill = "Model", palette = scrcols[3:4], facet.by = 'Sex',
         add = "mean_se", width = 0.5, add.params = list(size=0.3)) + 
  ggtitle('') + xlab('Scenario') + ylab(expression(paste(italic(alpha)[2]))) + #ylim(-0.2,0.4) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.2) + 
  theme(strip.background = element_blank(), 
        strip.text.x = element_blank())
fig5 <- ggarrange(Dcov, NULL, p0.cov,NULL, alpha2 ,ncol = 1, nrow = 5, heights = c(1,-0.15,1,-0.15,1),
                  common.legend=TRUE, align = 'hv') # labels = c('(a)','(b)','(c)'),

tiff(filename="Fig5.tif",height=5200,width=4000,units="px",res=800,compression="lzw")
fig5
dev.off()

'------ Fig 6 ------'
# Note that data here is generated in 6_Post-hoc.R
tiff(filename="Fig6.tif",height=6500,width=5200,units="px",res=800,compression="lzw")
par(mfrow=c(3,2), mar = c(1.1, 2.1, 2.1, 1.1), oma = c(4,3,3,0))

# A
sex_ind <- 'F'; sigma <- 9.7022; N = 60; colz <- 'red' # F; sigma is the mean inter-occasion movement estimate across replicates
sex_ind <- 'M'; sigma <- 14.6004; N = 40; colz <- 'blue' # M
# true
hist(do.call(c,ACdists[[1]][ACdists[[1]]$sex == sex_ind, 4:7]), breaks = seq(0,180,5),
     main = '', xlab=' ',ylab=' ',
     col = rgb(0,0,0,0.05), border=rgb(0,0,0,0.0), ylim = c(0,80), cex.axis = 1.5)
for (i in 2:100){
  hist(do.call(c,ACdists[[i]][ACdists[[i]]$sex == sex_ind, 4:7]), breaks = seq(0,180,5),
       main = '', xlab=' ',ylab=' ',
       col = rgb(0,0,0,0.05), border=rgb(0,0,0,0.0),add=T)
}
box(bty="l")
# predicted
tau <- 1/(sigma^2)
s <- cbind(rnorm(4*N,0,sigma),rnorm(4*N,0,sigma))
dists <- apply(s, 1, function(x){
  sqrt(x[1]^2 + x[2]^2)
})
hist(dists, breaks = seq(0,ceiling(max(dists)/5)*5,5), col=rgb(0,0,0,0), border=colz,add=T)

# B
sex_ind <- 'F'; sigma <- 8.6513; N = 54; colz <- 'red' # F
sex_ind <- 'M'; sigma <- 13.1198; N = 36; colz <- 'blue' # M
# true
hist(do.call(c,ACdists[[1]][ACdists[[1]]$X90th==1 & ACdists[[1]]$sex == sex_ind, 4:7]), breaks = seq(0,180,5),
     main = '', xlab=' ',ylab='',
     col = rgb(0,0,0,0.05), border=rgb(0,0,0,0.0), ylim = c(0,80), cex.axis = 1.5)
for (i in 2:100){
  hist(do.call(c,ACdists[[i]][ACdists[[i]]$X90th==1 & ACdists[[i]]$sex == sex_ind, 4:7]), breaks = seq(0,180,5),
       main = '', xlab=' ',ylab='',
       col = rgb(0,0,0,0.05), border=rgb(0,0,0,0.0),add=T)
}
box(bty="l")
# predicted
tau <- 1/(sigma^2)
s <- cbind(rnorm(4*N,0,sigma),rnorm(4*N,0,sigma))
dists <- apply(s, 1, function(x){
  sqrt(x[1]^2 + x[2]^2)
})
hist(dists, breaks = seq(0,ceiling(max(dists)/5)*5,5), col=rgb(0,0,0,0), border=colz,add=T)

# C
sex_ind <- 'F'; sigma <- 8.6513; N = 48; colz <- 'red' # F
sex_ind <- 'M'; sigma <- 13.1198; N = 32; colz <- 'blue' # M
# true
hist(do.call(c,ACdists[[1]][ACdists[[1]]$X80th==1 & ACdists[[1]]$sex == sex_ind, 4:7]), breaks = seq(0,180,5),
     main = '', xlab='',ylab=' ',
     col = rgb(0,0,0,0.05), border=rgb(0,0,0,0.0), ylim = c(0,80), cex.axis = 1.5)
for (i in 2:100){
  hist(do.call(c,ACdists[[i]][ACdists[[i]]$X80th==1 & ACdists[[i]]$sex == sex_ind, 4:7]), breaks = seq(0,180,5),
       main = '', xlab='',ylab=' ',
       col = rgb(0,0,0,0.05), border=rgb(0,0,0,0.0),add=T)
}
box(bty="l")
# predicted
tau <- 1/(sigma^2)
s <- cbind(rnorm(4*N,0,sigma),rnorm(4*N,0,sigma))
dists <- apply(s, 1, function(x){
  sqrt(x[1]^2 + x[2]^2)
})
hist(dists, breaks = seq(0,ceiling(max(dists)/5)*5,5), col=rgb(0,0,0,0), border=colz,add=T)

mtext("Inter-occasion activity centre distance",side=1,line=2,outer=TRUE,cex=1.3)
mtext("Count",side=2,line=1,outer=TRUE,cex=1.3,las=0)
dev.off()

'------ Fig 7 ------'
# plot (Fig. 7)
tiff(filename="Fig7.tif",height=4800,width=5200,units="px",res=800,compression="lzw")
par(mfrow=c(1,1), mar = c(6,4.1,4.1,1), xpd=T)
plot(nDetMin_c,rep(ztrap,100),  pch = 19, cex = 0.5, col = alpha('black', alpha = 0.1), 
     xlab = 'No. detections per trap', ylab = 'Resource (scaled)', frame.plot = FALSE,
     cex.lab=1.2)
points(nDetHr_c,rep(ztrap,100),  pch = 19, cex = 0.5, col = alpha('deepskyblue3', alpha = 0.2))
points(nDetOcc_c,rep(ztrap,100),  pch = 19, cex = 0.5, col = alpha('darkorange3', alpha = 0.1))
legend("top", legend=c("Min", "Day", "Occ"),
       col=c("black", "deepskyblue3", "darkorange3"), pch = 19, cex=1, bty='n', horiz = T,
       y.intersp=0, x.intersp=.5, inset=-0.1)
box(bty="l")
dev.off()

par(mfrow=c(1,1), mar = c(6,4.1,4.1,1))
plot(nDetHr_c,rep(ztrap,100),  pch = 19, cex = 0.5, col = alpha('deepskyblue3', alpha = 0.1), 
     xlab = '', ylab = '', frame.plot = FALSE,cex.axis=1.5)
box(bty="l")
plot(nDetOcc_c,rep(ztrap,100),  pch = 19, cex = 0.5, col = alpha('darkorange3', alpha = 0.1), 
     xlab = '', ylab = '', frame.plot = FALSE,cex.axis=1.5)
box(bty="l")
