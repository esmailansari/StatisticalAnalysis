dip.dist.plot <- function(Dip.dist){
  nbreaks = 15 #for hist plot
  opar <- par(no.readonly=TRUE,mfrow=c(1,1))
  par(fig=c(0, 0.8, 0 ,0.8))
  hist(Dip.dist,freq=FALSE,breaks=nbreaks,xlab='Dip angle',col=adjustcolor('burlywood3',alpha=3/8),main=NULL)
  #abline(v=mean(Dip.dist),col='red',lwd=3)
  rug(jitter(Dip.dist), col='burlywood3')
  lines(density(Dip.dist)$x, density(Dip.dist)$y,
        col='burlywood3', lwd=2)
  polygon(density(Dip.dist), col=adjustcolor('darkgreen',alpha=1/8), 
          border='burlywood3')
  par(fig=c(0, .8, 0.55, 1), new=TRUE)
  boxplot(Dip.dist, horizontal=TRUE, axes=FALSE,col=adjustcolor('burlywood3',alpha=2/8),main=NULL)
  points(mean(Dip.dist),pch=21,lwd=20)
}