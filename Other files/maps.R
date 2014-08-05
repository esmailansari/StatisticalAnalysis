rm(list=ls())
library(maps)
library(ggplot2)
library(ggmap)

mydata <- read.csv("Well-Test summary.csv", header=TRUE)
data<-mydata
mydata <- mydata[mydata$long<0,]
names(mydata)[names(mydata)=="Flow.Rate..BPD."]="Flow_Rate"
states <- map_data("county","louisiana")
#p <- ggplot(states, aes(long, lat)) +  
#  geom_polygon(aes(group=group), colour='black',alpha=0.6)
p <- ggplot()
p <- p + geom_polygon( data=states, aes(x=long, y=lat, group = group),colour="white")
#p <-  p + theme(panel.grid.minor = element_line(colour="white", size=0.5))# + 
  #scale_x_discrete(minor_breaks = seq(-94, -89, 0.2)) +
  #scale_y_discrete(minor_breaks = seq(29, 33, 0.2))
p <- p + geom_jitter(data=mydata,  aes(x=long, y=lat, size = mydata$Permeability..mD., color="coral1"),label=rownames(mydata)) + theme(legend.position="none")
#+ scale_size(name="Total enrollment")
#p <- p + geom_text( data=mydata, aes(x=long, y=lat, label=rownames(mydata), colour="gold2", size=4 ))
print(p)
#gglocator(1)