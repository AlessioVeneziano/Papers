#################################################################################
############################################ Get occurrence data - Veneziano 2024

### Load packages
library(rredlist)
library(countrycode)
library(rgbif)
library(geodata)
library(CoordinateCleaner)
library(terra)
library(sp)
library(raster)


### Set environment
# The object 'token' is a personal IUCN API key for accessing IUCN data from the
# package 'rredlist'. Request it at https://apiv3.iucnredlist.org/api/v3/token.
# If the link does not work, contact me at veneziano.alessio@gmail.com for help.
# The object 'urban' is the built-environment index from the WorldCover database
# (www.esa-worldcover.org). A list of species is fed to the script to gather data
# of species occurrence coordinates automatically.
token<-"..."
urban<-raster("data/occurrence/WorldCover/WorldCover_built_150s.tif")

spelist<-read.table("data/occurrence/species-list_Veneziano2024.txt",header=F,sep="\t")
  spelist<-spelist$V1


### Iterate procedure
# The following for loop iterates a procedure for download of occurrence data
# from the Global Biodiversity Information Facility (GBIF, www.gbif.org) via the
# package 'rgbif'.
for( i in 1:length(spelist)){
  spe<-spelist[i]
  
  print(paste("Processing:",spe))
  
  ### Retrieve and save IUCN country information
  cou<-rl_occ_country(spe,key=token)$result
  if(length(cou)!=0){
    cou<-cou[cou$origin=="Native",]
  } else {next}
  
  
  ### Get country shapefiles
  cou2<-cou$code
  cou3<-countrycode(cou2,'iso2c','iso3c')
  
  saveas<-"data/occurrence/shapefiles"
  shcou<-gadm(cou3,0,saveas,resolution=2)
  
  
  ### Get occurrence data
  nam<-name_suggest(spe,rank="species")
    nam<-unique(nam$data$canonicalName)
    nam<-na.omit(nam)
  if(length(nam)>1){nam=spe}
  occ<-occ_data(scientificName=nam,limit=10000)
    occ<-as.data.frame(occ$data)
    occ<-occ[!is.na(occ$decimalLongitude) | !is.na(occ$decimalLatitude),]
  
  if(nrow(occ)!=0){
    print(paste(nrow(occ),"found before cleaning!"))
  } else {next}
  
  
  ### First cleaning
  # Cleaning based on type of occurrence and proximity to capitals, country
  # centroids, institutions via the package 'CoordinateCleaner'.
  del<-c("FOSSIL_SPECIMEN","MATERIAL_SAMPLE","UNKNOWN")
  occ<-occ[!occ$basisOfRecord%in%del,]
  
  keep<-c("decimalLongitude","decimalLatitude","country","countryCode","basisOfRecord","year")
  occ<-occ[,keep]
  
  occ<-occ[occ$countryCode%in%cou2,]
  occ<-occ[!duplicated(occ[,1:2]),]
  
  occ<-clean_coordinates(occ,lon="decimalLongitude",lat="decimalLatitude",species=NULL,value="clean",
                         tests=c("capitals","centroids","institutions","zeros"))
  
  if(nrow(occ)!=0){
    print(paste(nrow(occ),"found after cleaning step 1!"))
  } else {next}
  
  
  ### Second cleaning
  # Cleaning based on country shapefiles (kept if within boundaries) and outliers.
  spcou<-as(shcou,"Spatial")
  sp<-SpatialPoints(cbind(occ[,1:2]))
  proj4string(sp)<-slot(spcou, "proj4string")
  keep<-!is.na(over(sp,spcou)$GID_0)
  occ<-occ[keep,]
  
  qu<-quantile(occ[,1])
  lim<-c(qu[2]-1.5*IQR(occ[,1]),qu[4]+1.5*IQR(occ[,1]))
  outlon<-occ[,1]<lim[1] | occ[,1]>lim[2]
  qu<-quantile(occ[,2])
  lim<-c(qu[2]-1.5*IQR(occ[,2]),qu[4]+1.5*IQR(occ[,2]))
  outlat<-occ[,2]<lim[1] | occ[,2]>lim[2]
  occ<-occ[(outlon+outlat)==0,]
  
  if(nrow(occ)!=0){
    print(paste(nrow(occ),"found after cleaning step 2!"))
  } else {next}
  
  
  ### Third cleaning
  # Cleaning based on proximity to urban areas (from the WorldCover database).
  urb<-raster::crop(urban,extent(spcou))
    urb<-raster::mask(urb,spcou)
    urb[urb[]>0.4]<-NA
    urb[!is.na(urb)]<-1
  
  keep<-!is.na(urb[cellFromXY(urb, occ[,1:2])])
  occ<-occ[keep,]
  
  saveas<-paste("data/occurrence/occurrence",sep="",paste(spe,sep="_occurrence",".txt"))
  write.table(occ,saveas,sep="\t",row.names=F,quote=F)
  
  if(nrow(occ)!=0){
    print(paste(nrow(occ),"found after cleaning step 3!"))
  } else {next}
  
  ### Clean environment and proceed to following iteration
  rm(list=ls()[!(ls()%in%c("spelist","token","urban"))])
  gc()
}



################################################################### END OF SCRIPT
#################################################################################


