library(sp)
library(rgeos)
#construct spatial classes and perform geo-processing.
#Read the data and transform them to spatial objects
flowCoord <- flowSiteList[1,c(1,2,3)]
siloCoord <- weatherSiteList[,c(1,3,4)]
merg <- rbind(flowCoord,siloCoord)
sp.merg <- merg
coordinates(sp.merg) <- ~longitude+latitude

#calculate pairwise distances between points
d <- rgeos::gDistance(sp.merg, byid=T)

#Find second shortest distance (closest distance is of point to itself, therefore use second shortest)
min.d <- apply(d, 1, function(x) order(x, decreasing=F)[2])

#Result
newdata <- cbind(merg, merg[min.d,], apply(d, 1, function(x) sort(x, decreasing=F)[2]))
nearestStationNumber<-as.numeric(newdata[1,4])
nearestStationNumber<-as.character(nearestStationNumber)
