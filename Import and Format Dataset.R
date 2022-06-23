library(plyr)

rawdata = read.csv("~/Dropbox/Leman Research/Cyclocross/Data Files/Survey_Data.csv",header=TRUE)
rawdata = rawdata[rawdata$Place!="DNS",]   # Removing competitors that did not show

rawdata$fieldgroup = rep(NA, nrow(rawdata))  # Grouping together fields into field groups
rawdata[rawdata$Field == "Women Junior 15-16", ][, "fieldgroup"] = 1 
rawdata[rawdata$Field == "Women Junior 17-18", ][, "fieldgroup"] = 1 
rawdata[rawdata$Field == "Men Junior 15-16", ][, "fieldgroup"] = 2
rawdata[rawdata$Field == "Men Junior 17-18", ][, "fieldgroup"] = 2
rawdata[rawdata$Field == "Men Elite 23+ Pro/Cat 1/2", ][, "fieldgroup"] = 3
rawdata[rawdata$Field == "Women Elite/U23 Pro/Cat 1/2/3", ][, "fieldgroup"] = 4
rawdata[rawdata$Field == "Women U23", ][, "fieldgroup"] = 4
rawdata[rawdata$Field == "Men U23", ][, "fieldgroup"] = 5

g1=c("MA")
g2=c("MT","TN","VT")
g3=c("ID","ME","NY")
g4=c("NJ","WI","NM","WA","NH")
g5=c("OR","VA","CA","NC","AZ","CO","OH")
g6=c("LA","MN","KY","MD","TX","IL","KS","OK")
g7=c("AL","AR","CT","FL","GA","IA","IN","MI","MO","ND","PA","UT")

rawdata$stategroup = rep(NA,nrow(rawdata))  # Grouping States into state groups (used K-Means)
rawdata[rawdata$State %in% g1,][,"stategroup"] = 1
rawdata[rawdata$State %in% g2,][,"stategroup"] = 2
rawdata[rawdata$State %in% g3,][,"stategroup"] = 3
rawdata[rawdata$State %in% g4,][,"stategroup"] = 4
rawdata[rawdata$State %in% g5,][,"stategroup"] = 5
rawdata[rawdata$State %in% g6,][,"stategroup"] = 6
rawdata[rawdata$State %in% g7,][,"stategroup"] = 7

fitdata = rawdata[rawdata$Survey.Indicator==1,]  # Surveyed data used to fit the model
extrap = rawdata[rawdata$Survey.Indicator==0,c(1,10,11)]   # Non-surveyed data that requires estimates


# These are counts and cost sums used to calculate the posterior cost distribution
extrap.stategroup.count = as.vector(table(extrap$stategroup))
extrap.fieldgroup.count = as.vector(table(extrap$fieldgroup))
fit.stategroup.count = c(as.vector(table(fitdata$stategroup)),0)
fit.fieldgroup.count = as.vector(table(fitdata$fieldgroup))
fit.stategroup.tcost = c(aggregate(fitdata$Surveyed.Travel.Cost,by=list(fitdata$stategroup),FUN=sum)[,2],0)
fit.fieldgroup.ocost = aggregate(fitdata$Surveyed.Other.Cost,by=list(fitdata$fieldgroup),FUN=sum)[,2]

a = table(extrap$fieldgroup,extrap$stategroup)


#################################################################################
#Format Data to Create Plot of Clusters
#################################################################################

temp=fitdata[,c("Surveyed.Travel.Cost","State","stategroup")]
plotdata=aggregate(temp,by=list(temp[,"State"]),FUN="mean")

plot(NA,NA,xlim=c(1,6.1),ylim=c(0,1300),axes=FALSE,xlab="State Cluster",ylab="Travel Cost",main="State Clustering for Travel Costs")
axis(1,seq(1,6,by=1))
axis(2,seq(0,1300,by=100))
points(jitter(plotdata$stategroup[plotdata$stategroup==1],0),plotdata$Surveyed.Travel.Cost[plotdata$stategroup==1],pch=16,col="darkred",cex=.8)
points(jitter(plotdata$stategroup[plotdata$stategroup==2],1.2),plotdata$Surveyed.Travel.Cost[plotdata$stategroup==2],pch=16,col="red3",cex=.8)
points(jitter(plotdata$stategroup[plotdata$stategroup==3],1.2),plotdata$Surveyed.Travel.Cost[plotdata$stategroup==3],pch=16,col="orangered",cex=.8)
points(jitter(plotdata$stategroup[plotdata$stategroup==4],1.2),plotdata$Surveyed.Travel.Cost[plotdata$stategroup==4],pch=16,col="darkorange2",cex=.8)
points(jitter(plotdata$stategroup[plotdata$stategroup==5],1.2),plotdata$Surveyed.Travel.Cost[plotdata$stategroup==5],pch=16,col="orange",cex=.8)
points(jitter(plotdata$stategroup[plotdata$stategroup==6],1.2),plotdata$Surveyed.Travel.Cost[plotdata$stategroup==6],pch=16,col="yellow",cex=.8)

#################################################################################
#Format Data to Create Map
#################################################################################
# NEED TO RE-INSTALL PACKAGES TO MAKE A MAP

# test=data.frame(state.abb)
# colnames(test)="Group.1"
# mapdata=join(test,plotdata,by="Group.1",type="left")[,1:2]
# mapdata=mapdata[mapdata$Group.1!="AK"&mapdata$Group.1!="HI",]
# mapdata[is.na(mapdata$travel),2]=-1
# mapdata

par(mfrow=c(3,2))
for(i in 1:6)
{
qqnorm(fitdata$Surveyed.Travel.Cost[fitdata$stategroup==i],axes=FALSE,main=paste("State Group",i,"Travel Expense"))
axis(1,at=seq(-2,2,by=.5))
axis(2,at=seq(0,1500,by=300))
qqline(fitdata$Surveyed.Travel.Cost[fitdata$stategroup==i])
}

for(i in 1:5)
{
qqnorm(fitdata$Surveyed.Other.Cost[fitdata$fieldgroup==i],axes=FALSE,main=paste("Field",i,"Other Expense"))
axis(1,at=seq(-2,2,by=.5))
axis(2,at=seq(0,1500,by=300))
qqline(fitdata$Surveyed.Other.Cost[fitdata$fieldgroup==i])
}

par(mfrow=c(1,1))
qqnorm(fitdata$Surveyed.Travel.Cost[fitdata$stategroup==5],axes=FALSE,main=paste("State Group 5 Travel Expense"))
axis(1,at=seq(-2,2,by=.5))
axis(2,at=seq(0,1500,by=300))
qqline(fitdata$Surveyed.Travel.Cost[fitdata$stategroup==5])

qqnorm(fitdata$Surveyed.Other.Cost[fitdata$fieldgroup==2],axes=FALSE,main=paste("Men Juniors Other Expense"))
axis(1,at=seq(-2,2,by=.5))
axis(2,at=seq(0,1500,by=300))
qqline(fitdata$Surveyed.Other.Cost[fitdata$fieldgroup==2])
