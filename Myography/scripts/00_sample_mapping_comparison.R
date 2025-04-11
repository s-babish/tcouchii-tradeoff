setwd("C:/Users/sdbab/OneDrive - University of Nevada, Reno/UNR/tcouchii/Myography")

#Libraries -----
library(sf)
library(mapview)
library(webshot)

#load data
#C4P samples:
c4p <- read.csv("Tcouchii_side_hustle/Couchii_C4P/output/Couchii_C4P.csv")

#Tetanus samples:
tetanus <- read.csv("Tcouchii_side_hustle/Couchii_tetanus/output/Couchii_tetanus.csv")

#IC50 samples:
IC50 <- read.csv("IC50.csv")
IC50 <- IC50[!is.na(IC50[10]),]

#Rheobase samples (assuming also in tetanus and c4p):
rheo <- read.csv("200_TPA-sigmoidal-rpt.csv",header=T)

#Plot sample locations -----

#make storage matrix for all of them
allsamps <- matrix(NA,nrow = sum(nrow(IC50),nrow(rheo)),ncol = 4)
colnames(allsamps) <- c("Snake","Latitude","Longitude","Status")
index = 1
for (i in 1:nrow(IC50)) {
  if (IC50[i,2] %in% rheo$Snake) {
    status = "both"
  } else status = "IC_only"
  allsamps[index,] <-c(IC50[i,2], as.numeric(IC50[i,7]), as.numeric(IC50[i,8]), status)
  index = index + 1
}

for (j in 1:nrow(rheo)) {
  if (!(rheo[j,7] %in% IC50$Snake)) {
    allsamps[index,] <- c(rheo[j,7],as.numeric(rheo[j,13]),as.numeric(rheo[j,12]),"Rheo_only")
    index = index + 1
  }
}

#Plot samples
allsampsdf <- as.data.frame(allsamps[!is.na(allsamps[,1]),])
allsampsdf <- transform(allsampsdf, Latitude = as.numeric(Latitude), 
                        Longitude = as.numeric(Longitude), Status = as.factor(Status))
myo_map <- mapview(allsampsdf, xcol="Longitude", ycol="Latitude", zcol=c("Status"),crs=4269, grid=F, legend=T, 
                   map.types=c("CartoDB.Positron","Esri.WorldImagery","OpenTopoMap","Esri.WorldPhysical","USGS.USTopo"),
                   cex=3)
myo_map
#mapshot(myo_map, url="myographysamps.html")

#Make sheet that shows which samples are in which sheets
allIDs <- unique(c(IC50$Snake,rheo$Snake,c4p$Snake,tetanus$Snake))


sample_masterlist <- matrix(NA,nrow=length(allIDs),ncol=5)
colnames(sample_masterlist) <- c("Snake","IC50","Rheobase","C4P","Tetanus")
for (i in 1:nrow(sample_masterlist)) {
  ID <- allIDs[i]
  sample_masterlist[i,] <- c(ID, ID %in% IC50$Snake, ID %in% rheo$Snake, ID %in% 
                           c4p$Snake, ID %in% tetanus$Snake)
}
write.csv(sample_masterlist,"myography_samples.csv")
