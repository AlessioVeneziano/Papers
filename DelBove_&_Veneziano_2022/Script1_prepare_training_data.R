##########################################################
###################### Del Bove & Veneziano 2022: SCRIPT 1

# Load functions
source("Functions.R")


# Read Howell dataset (you need an internet connection to read the data)
how<-read.csv("https://web.utk.edu/~auerbach/Howell.csv")


# Compute table showing sex by population and remove populations represented
# by only one sex
tab<-as.matrix(table(how$Population,how$Sex))
del<-rownames(which(tab==0,T))
  how<-how[!how$Population%in%del,]

# Keep a subset of measurements (refer to the main text for definitions)
meas<-c("GOL","BBH","ZYB","AUB","NLH","OBH","OBB","NLB","MDH","OCC")
vars<-c("Sex","Population",meas)
  how<-how[,vars]


# Z-score transformation to remove observations exhibiting 3 or more standard
# deviations from the mean for at least one measurement
rm(del)
scaled<-scale(how[,-c(1,2)])
del<-unique(which(abs(scaled)>3,T)[,1])
  how<-how[-del,]

# Homogenise sex proportion (collecting equal numbers of males and females)
dif<-diff(table(how$Sex))
set.seed(42)
subm<-sample(which(how$Sex=="M"),dif,F)
  how<-how[-subm,]


# Compute Geometric Mean and add to training dataset
how$GM<-apply(how[,3:12],1,gmean)


# Compute Mosimann transformation
how[,3:12]<-t(apply(how[,3:12],1,mosimann))
how[,3:13]<-round(how[,3:13],3)


# Save transformed data
write.csv(how,"how_transformed.csv",quote=F,row.names=F)


############################################# END OF SCRIPT
##########################################################


