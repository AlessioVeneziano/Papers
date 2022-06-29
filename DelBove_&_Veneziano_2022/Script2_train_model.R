##########################################################
###################### Del Bove & Veneziano 2022: SCRIPT 2

# Load packages
library(h2o)


# Set memory and h2o environment (the memory allocated and the specifics of
# 'h2o.init' can be changed depending on the specifics of your system); the 
# function 'h2o.shutdown' removes any h2o environment previously initialised;
# the launch of a new h2o environment can take a few seconds
memory.limit(1e+09)

h2o.shutdown()
h2o.init(nthreads=-1, max_mem_size="4g")


# Read dataset (this is the dataset prepared in the Script 1)
how<-read.csv("how_transformed.csv")


# Prepare training and validation sets from the whole training dataset (for
# consistency with the h2o package, the data is transformed into an object of
# class'H2OFrame')
dtv<-how[,-2]
dtv$Sex<-as.factor(dtv$Sex)
  dtv<-as.h2o(dtv)

# Define outcome (sex) and variables
xs<-colnames(dtv)[-1]
ys<-colnames(dtv)[1]


# Prepare parameters for model tuning using a 'brute force' approach: generate
# a dataframe with all the combinations of the values to test for each parameter
# to tune
hlist<-list(3,7,11,15,19,23,27,c(3,3),c(7,7),c(11,11),c(15,15),c(19,19),c(23,23),c(27,27))
l2list<-list(0,1e-4,5e-4,1e-3,5e-3,1e-2,5e-2)
  pars<-expand.grid(hlist,l2list)

# Model tuning: models with different combinations of parameters are fitted on
# the training set and their performance is assessed on the validation set
# (a more computationally efficient tuning can be attained using the function
# 'h2o.grid', but due to system limitations we preferred a for loop - it can
# take up to a few hours)
stats<-matrix(NA,nrow(pars),40)
for(i in 1:nrow(pars)){
  print(i)
  
mod<-h2o.deeplearning(x=xs,y=ys,training_frame=dtv,epochs=1000,stopping_metric="logloss",
                      loss="CrossEntropy",distribution="bernoulli",standardize=T,
                      adaptive_rate=F,rate=0.0005,momentum_start=0.5,momentum_ramp=1e6,
                      momentum_stable=0.99,stopping_rounds=20,stopping_tolerance=0.0001,
                      input_dropout_ratio=0,reproducible=T,seed=42,nfolds=10,
                      fold_assignment="Stratified",activation="RectifierWithDropout",
                      fast_mode=F,hidden=pars[i,1][[1]],l2=pars[i,2][[1]])

  statcv<-mod@model$cross_validation_metrics_summary[,1:2]
  stats[i,]<-c(t(statcv))
  
  rm("mod")
  if(i%%10==0){print("Cleaning in progress...") ; gc()}
}

# Put together the results obtained from the tuning ('mean' and 'sd' refer to the
# statistics calculated for each model tuned and computed over the cross-validated
# results) and save the table
stats<-cbind.data.frame(pars,stats)
  colnam<-c(t(cbind(rownames(statcv),rownames(statcv))))
  colnam<-paste(colnam,c("mean","sd"),sep="_")
    colnames(stats)<-c("hlayer","l2reg",colnam)

saveRDS(stats,"tunedModel_stats.RData")


# Choosing best model: the choice can be based on multiple considerations (refer
# to the main text of the paper for the criteria used); here, for example, we
# choose the model maximising the Area Under Curve (AUC) of a Receiver Operating
# Characteristic curve (ROC); the parameters of the best model are extracted
bmPars<-stats[which.max(stats$auc_mean),1:2]
  par1<-unlist(bmPars[1])
    names(par1)<-NULL
  par2<-unlist(bmPars[2])
    names(par2)<-NULL


# Train best model using the whole training data (no validation set provided - 
# this is the final model that can be used to estimate sex to new observations)
mod<-h2o.deeplearning(x=xs,y=ys,training_frame=dtv,epochs=1000,stopping_metric="logloss",
                      loss="CrossEntropy",distribution="bernoulli",standardize=T,
                      adaptive_rate=F,rate=0.0005,momentum_start=0.5,momentum_ramp=1e6,
                      momentum_stable=0.99,stopping_rounds=20,stopping_tolerance=0.00001,
                      input_dropout_ratio=0,reproducible=T,seed=42,activation="RectifierWithDropout",
                      variable_importances=T,fast_mode=F,export_weights_and_biases=T,
                      hidden=par1,l2=par2)


# Save best model for later use
h2o.saveModel(mod,getwd(),filename="BestModel_hidden11_l2reg1e-04")



############################################# END OF SCRIPT
##########################################################



