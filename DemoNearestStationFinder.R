#Need site and station lists
#Return ID of the neareast SILO station 

library(sp)
library(rgeos)

#Get Site and station lists
flowSiteList <- read.table("siteList.txt",header=T,sep = "",dec = ".")

weatherSiteList <- getWeatherSiteList()


#construct spatial classes and perform geo-processing.
#Read the data and transform them to spatial objects

flowCoord <- flowSiteList[1,c(1,2,3)] #Get coordinate of a site

siloCoord <- weatherSiteList[,c(1,3,4)] #Get coordinate of all SILO site

merg <- rbind(flowCoord,siloCoord) #Merge the 2 above into 1 dataframe

#Turn the merged dataframe into spatia; class
sp.merg <- merg
coordinates(sp.merg) <- ~longitude+latitude

#calculate pairwise distances between points
d <- rgeos::gDistance(sp.merg, byid=T)

#Find second shortest distance (closest distance is of point to itself, therefore use second shortest)
min.d <- apply(d, 1, function(x) order(x, decreasing=F)[2])

#Result
newdata <- cbind(merg, merg[min.d,], apply(d, 1, function(x) sort(x, decreasing=F)[2])) #sort order of station the are near to the site

#Pick station name
nearestStationNumber<-as.numeric(newdata[1,4])#get rid of the space before the station ID
nearestStationNumber<-as.character(nearestStationNumber)#turn back to character class
