rm(list=ls())
#if (!is.null(dev.list())) dev.off()
require(MASS)
require(fitdistrplus)
#setwd("C:/Users/eansar2/Google Drive/Ph.D. Reports/Ranges and Distributions")
#setwd("~/Desktop/Google Drive/Ph.D. Reports/Ranges and Distributions")
data <- read.csv("./Data/Data.csv")

#       plotting and fitting distributions
#--------------------------------------------
#                   Permeability
#-------------------------------------------
opar<-par()
data$Temperature <- (data$Temp...degree.F.-32)/1.8
data$Depth <- (data$Depth..ft.*0.3048)
perm.data<-data[!is.na(data$Permeability..mD.),]
permeability<-perm.data$Permeability..mD.
perm.data$color[perm.data$Source=="Source 1"]<-'red'
perm.data$color[perm.data$Source=="Source 2"]<-'blue'
par(mar=c(4,10,4,2))

fit.perm <- fitdist(perm.data$Permeability..mD.,"lnorm","mme")

plot(fit.perm)

summary(fit.perm)

bc <- boxcox(permeability~1)
lambda <- bc$x[which.max(bc$y)]

# Residuals
## boxCox(permeability~1) # gives lambda = 0.0606
qq.bc <- qqnorm(sort(permeability^(lambda)),type="l",plot.it=FALSE)
qq.none <- qqnorm(sort(permeability),plot.it=FALSE)
qq.log <- qq.log <- qqnorm(sort(log(permeability,10)),
                           type="l",plot.it=FALSE)
y.bc <- (qq.bc$y - mean(qq.bc$y))/sd(qq.bc$y)
y.none <- (qq.none$y - mean(qq.none$y))/sd(qq.none$y)
y.log <- (qq.log$y - mean(qq.log$y))/sd(qq.log$y)
y.lim <- range(y.bc,y.log,y.none)
y.lim[2] <- 4
x.lim <- range(qq.bc$x,qq.log$x,qq.none$x)
plot(1,1,xlim=x.lim,ylim=y.lim,xlab="Theoretical quantile",
     ylab="Standardized observed quantile",
     type="n",xpd=FALSE)
abline(0,1,lty=2)
lines(qq.none$x, y.none, lty=1,col="red")
lines(qq.log$x, y.log, lty=1, col="blue")
lines(qq.bc$x, y.bc, lty=1, col="orange")
legend("bottomright",c("Normal","No transform","Log", "Box-Cox"),
       lty=c(2,1,1,1),col=c("black","red","blue","orange"),bty="n")

#---------------------------------------------
#                   Porosity
#---------------------------------------------
poro.data <- data[!is.na(data$Porosity),]
porosity  <- poro.data$Porosity
poro.data$color[poro.data$Source=="Source 1"]<-'red'
poro.data$color[poro.data$Source=="Source 2"]<-'blue'
par(mar=c(4,10,4,2))




poro.vec <- sort(porosity)
P.vec <- (1:length(poro.vec))/(length(poro.vec) + 1)
poro.ns <- qnorm(P.vec); 
poro.ns <- poro.ns*sd(poro.vec)+mean(poro.vec)
par(mfrow = c(1, 2)) ; 
hist(porosity)
hist(poro.ns)
poro.stat<-list(porosity.sd=sd(poro.ns),porosity.mean=mean(poro.ns))
print(poro.stat)

fit.poro <- fitdist(poro.ns,"norm","mme")
plot(fit.poro)
summary(fit.poro)


#-------------------------------------------------------
#                    Average reservoir temperature
#--------------------------------------------------------

temp.data <- data[!is.na(data$Temp...degree.F.),]
temperature  <- (temp.data$Temp...degree.F.-32)*5/9
temp.data$color[temp.data$Source=="Source 1"]<-'red'
temp.data$color[temp.data$Source=="Source 2"]<-'blue'
par(mar=c(4,10,4,2))



temp.vec <- sort(temperature)
P.vec <- (1:length(temp.vec))/(length(temp.vec) + 1)
temp.ns <- qnorm(P.vec); 
temp.ns <- temp.ns*sd(temp.vec)+mean(temp.vec)
par(mfrow = c(1, 2)) ; 
hist(temperature)
hist(temp.ns)
temp.stat<-list(temperature.sd=sd(temp.ns),
                temperature.mean=mean(temp.ns))
print(temp.stat)

fit.temp <- fitdist(temp.ns,"norm","mme")
plot(fit.temp)
summary(fit.temp)



#-----------------------------------------------
#                       Thicknesss 
#----------------------------------------------
thick.data <- data[!is.na(data$Average.thickness.ft),]
thickness  <- thick.data$Average.thickness.ft*0.3048
thick.data$color[thick.data$Source=="Source 1"]<-'red'
thick.data$color[thick.data$Source=="Source 2"]<-'blue'
par(mar=c(4,10,4,2))




thick.vec <- sort(thickness)
P.vec <- (1:length(thick.vec))/(length(thick.vec) + 1)
thick.ns <- qnorm(P.vec); 
thick.ns <- thick.ns*sd(thick.vec)+mean(thick.vec)
par(mfrow = c(1, 2)) ; 
hist(thickness)
hist(thick.ns)
thick.stat<-list(thickness.sd=sd(thick.ns),
                 thickness.mean=mean(thick.ns))
print(thick.stat)

fit.thick <- fitdist(thick.ns,"norm","mme")
plot(fit.thick)
summary(fit.thick)


#-----------------------------------------------
#                          Length
#-----------------------------------------------
area.data <- data[!is.na(data$Area.SqMi),]
area  <- area.data$Area.SqMi*(1609.34^2);
length <- sqrt(area)
#area.data$color[temp.data$Source=="Source 1"]<-'red'
#area.data$color[temp.data$Source=="Source 2"]<-'blue'
par(mar=c(4,10,4,2),mfrow=c(1,2))


length.vec <- sort(length)
P.vec <- (1:length(length.vec))/(length(length.vec) + 1)
length.ns <- qnorm(P.vec); 
length.ns <- length.ns*sd(length.vec)+mean(length.vec)
par(mfrow = c(1, 2)) ; 
hist(length)
hist(length.ns)
length.stat<-list(length.sd=sd(length.ns),
                  length.mean=mean(length.ns))
print(length.stat)

fit.length <- fitdist(length.ns,"norm","mme")
plot(fit.length)
summary(fit.length)

#-----------------------------------------------
#                 Earth thermal gradient
#-----------------------------------------------
depth.temp = data.frame(depth=(data$Depth..ft.*(0.3048)),
                        temp=((data$Temp...degree.F.-32)*5/9))
depth.temp = depth.temp[!is.na(depth.temp$depth),]
depth.temp = depth.temp[!is.na(depth.temp$temp), ]
depth.temp$temp.20=depth.temp$temp-20;
depth.temp.model = lm(temp.20~depth-1,data=depth.temp)
par(mfrow=c(1,1))

#plotting temp vs. depth
# y.noise is depth
# x is temp
n.plot <- 100
x <- depth.temp$temp
y.noise <- depth.temp$depth
fit.y <- lm (y.noise ~ x -1)
x.plot <- data.frame(x=seq(min(x),max(x),length=n.plot))
pred <- predict(fit.y,newdata=x.plot,se.fit=TRUE,
                interval="confidence")
fit.se <-pred$se.fit
fit.fit <- pred$fit[,"fit"]
conf.lwr <- pred$fit[,"lwr"]
conf.upr <- pred$fit[,"upr"]
pred <- predict(fit.y, newdata=x.plot, se.fit=FALSE,
                interval="prediction")
pred.lwr <- pred[,"lwr"]  #pred[["lwr"]] does not work neigher ...
pred.upr <- pred[,"upr"]
y.lim <- range(y.noise, fit.fit, conf.lwr, conf.upr, 
               pred.lwr,pred.upr) #returns max and min
plot(0,0,xlim=range(x,x.plot),ylim=y.lim,
     ylab="Depth (m)",type="n",axes=FALSE) #plot (0,0)
axis(2)
axis(3)
mtext("Temperature (C)", side=3, line=3) 
x.plot <- x.plot[,1] 
#because length of data frame is the number of its columns.
lines(x.plot,fit.se,lty=1,col=2)  
#class of x.plot can not be data frame. 
polygon(c(x.plot,rev(x.plot)),c(pred.lwr,rev(pred.upr)), 
        col=adjustcolor(7,alpha=1/8),lty=0)
polygon(c(x.plot,rev(x.plot)),c(conf.lwr,rev(conf.upr)), 
        col=adjustcolor(2,alpha=1/8),lty=0)
lines(x.plot,fit.fit,lty=2,col=1)
ix <- sort(x,index=TRUE)$ix
points(x,y.noise,pch=21,col=1,bg="white")
text(145,-3880,sprintf("Temp = %.4f Depth + 20",depth.temp.model$coefficients))
abline(v=135,lty=5,col=adjustcolor(2,alpha=2/8))


#-------------------------------------------------
#    plotting all distributions in one graph
#-----------------------------------------------
histplot <- function(Parameter,Para.name,Color,breaks=6){
par(fig=c(0,0.8,0,0.8))
hist(Parameter,freq=FALSE,xlab=Para.name,col=Color)
rug(jitter(Parameter), col=Color)
par(fig=c(0,0.8,0.55,1),new=TRUE)
boxplot(Parameter,horizontal=TRUE,axes=FALSE)
}

histplot(porosity,'porosity','blue')
histplot(thickness,'thickness','yellow')
histplot(length,'length','green')
histplot(temperature,'temperature','red')


#-------------------------------------------------------
#                     Bootstraping
#-------------------------------------------------------
Tau <- function(formula, data, indices){
  d <- data[indices,]
  fit <- lm(formula,data=d)
  return(fit$coefficients)
  
}

library(boot)
boot.results <- boot(data=depth.temp, statistic=Tau, R=1000, formula=-temp.20~depth-1)
print(boot.results)
boot.ci(boot.results,type=c("bca"))

#lines(density(permeability)$x, density(d.ResDip)$y,
#      col="red", lwd=2)

#-----------------------------------------------
#                Injection temperature
#-----------------------------------------------


#-----------------------------------------------
#               Dip Angle
#-----------------------------------------------
trunc.draws <-function(x, xlim) {
  x[x<xlim[1]]    <- xlim[1]
  x[x>xlim[2]]    <- xlim[2]
  return(x)
}
T_max = 180 #Degree C
source('./functions/Hammersley.R')
Factors <- list('T.avg','Length','Tau')
N=5000
T.mean  <- fit.temp$estimate[[1]]
T.sd    <- fit.temp$estimate[[2]]
T.lim   <- c(T.mean-3*T.sd,T.mean+3*T.sd) 
T.avg.dist = trunc.draws(rnorm(N,T.mean,T.sd),T.lim)
Length.mean  <- fit.length$estimate[[1]]
Length.sd    <- fit.length$estimate[[2]]
Length.lim   <- c(Length.mean-3*Length.sd,Length.mean+3*Length.sd) 
Length.dist = trunc.draws(rnorm(N,Length.mean,Length.sd),Length.lim)
fit.tau <- fitdist(as.vector(boot.results$t),"norm","mme")
Tau.mean   <- fit.tau$estimate[[1]]
Tau.sd     <- fit.tau$estimate[[2]]
Tau.lim    <- c(Tau.mean-3*Tau.sd,Tau.mean+3*Tau.sd) 
Tau.dist = trunc.draws(rnorm(N,Tau.mean,Tau.sd),Tau.lim)

Fact.dist <- data.frame(T.avg=T.avg.dist,Length=Length.dist,Tau=Tau.dist)
HSS.select <- as.data.frame(matrix(rep(NA,N*length(Factors)),
                                   c(N,length(Factors))))
colnames(HSS.select) <- Factors
colnames(Fact.dist)  <-Factors

QuasiMC <- hammersley(prime(50)[1:length(Factors)],N)

  
for (i in 1:length(Factors)){
  HSS.select[,i] <- QuasiMC[,i]*(max(Fact.dist[Factors[[i]]])-
                                   min(Fact.dist[Factors[[i]]]))+ 
                                  min(Fact.dist[Factors[[i]]])
}

Dip.dist = rep(0,N)
Dip.q.Dist <- 2*(T_max-HSS.select[,1])/(HSS.select[,2]*HSS.select[,3])
Dip.dist <- asin(Dip.q.Dist)
#Dip.dist <- Dip.dist[Dip.dist>0]
Dip.dist <- Dip.dist/pi*180

nbreaks = 15 #for hist plot
opar <- par(no.readonly=TRUE,mfrow=c(1,1))
par(fig=c(0, 0.8, 0 ,0.8))
hist(Dip.dist,freq=FALSE,breaks=nbreaks,xlab='Dip angle',col=adjustcolor('burlywood3',alpha=3/8))
#abline(v=mean(Dip.dist),col='red',lwd=3)
rug(jitter(Dip.dist), col='burlywood3')
lines(density(Dip.dist)$x, density(Dip.dist)$y,
      col='burlywood3', lwd=2)
polygon(density(Dip.dist), col=adjustcolor('lightsalmon',alpha=4/8), 
        border='burlywood3')
par(fig=c(0, .8, 0.55, 1), new=TRUE)
boxplot(Dip.dist, horizontal=TRUE, axes=FALSE,col=adjustcolor('burlywood3',alpha=2/8))
points(mean(Dip.dist),pch=21,lwd=20)



#----------------------------------------------
#          Population Tests
#-------------------------------------------
t.test(Permeability..mD.~Source, data=perm.data)
t.test(Porosity~Source, data=poro.data)
t.test(Temperature~Source, data=temp.data)
t.test(Depth~Source, data=temp.data)
summary(aov(Temperature~Depth*Source, data=temp.data))


#-----------------------------------------------
#                   Flow Rate
#----------------------------------------------
g1<-ggplot(data=perm.data,aes(Source,Permeability..mD.))+
  geom_boxplot(aes(fill=Source))+ scale_fill_brewer()+guides(fill=FALSE)
  

g2<-ggplot(data=perm.data,aes(Source,Porosity))+
  geom_boxplot(aes(fill=Source))+ geom_jitter() +guides(fill=FALSE) + scale_fill_brewer(palette='Greens')

g3<-ggplot(data=perm.data,aes(Source,Temperature))+
  geom_boxplot(aes(fill=Source))+ scale_fill_brewer(palette='OrRd')+guides(fill=FALSE)
  
# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


multiplot(g1,g2,g3,cols=3)