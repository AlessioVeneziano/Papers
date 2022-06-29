##########################################################
###################### Del Bove & Veneziano 2022: SCRIPT 3

# This script shows how to use the best model to estimate sex on the data
# measured on new crania.


# Load packages and functions
source("Functions.R")
library(h2o)


# Set memory and h2o environment
memory.limit(1e+09)

h2o.shutdown()
h2o.init(nthreads=-1, max_mem_size="4g")


# Read the model and the data to be estimated (it is assumed that the data comes
# in the same format as the training data - same measurements, same field names),
# except for the Mosimann transformation (which is performed below); here, just
# to explain how the process works, we use some example data from the Howell
# dataset
mod<-h2o.loadModel("BestModel_hidden11_L2reg1e-04")
newdat<-read.csv("example_new_data.csv")


# Transform data (add Geometric Mean and apply Mosimann transformation)
newdat$GM<-apply(newdat,1,gmean)

newdat[,-ncol(newdat)]<-t(apply(newdat[,-ncol(newdat)],1,mosimann))
newdat<-round(newdat,3)


# For consistency with the package 'h2o', the data is transformed into an object
# of class 'H2OFrame'
newdat<-as.h2o(newdat)


# Predict the sex of the individual crania measured in 'newdat'; the model produces
# a best estimate at a probability threshold that can be different from 0.5, but
# having trained the model with exactly equal number of females and males, we
# decide to assume a prior equal probability for each sex (thus we re-compute the
# estimated sex based on a 0.5 threshold)
pred<-h2o.predict(mod,newdat)
  pred<-as.data.frame(pred)
    pred$predict<-ifelse(pred$'F'>0.5,"F","M")



############################################# END OF SCRIPT
##########################################################




