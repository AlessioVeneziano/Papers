##########################################################
########################## Veneziano et al 2026: FUNCTIONS


# adds transparency to colours
colAlpha<-function(colour,alpha=0.3){
  acol<-c(col2rgb(colour),T)/255
  return(rgb(acol[1],acol[2],acol[3],alpha))
}

# transforms the output of lm and related model objects to show explicit slope and
# intercept for additional classes within a factor (it only works for two classes)
modelEst2groups<-function(model,rounded=3){
  coefs<-coef(model)
  vc<-vcov(model)
  
  intercept_A<-coefs[1]
  slope_A<-coefs[2]
  se_intercept_A<-sqrt(vc[1,1])
  se_slope_A<-sqrt(vc[2,2])
  
  intercept_B<-coefs[1]+coefs[3]
  slope_B<-coefs[2]+coefs[4]
  
  se_intercept_B<-sqrt(vc[1,1]+vc[3,3]+2*vc[1,3])
  se_slope_B<-sqrt(vc[2,2]+vc[4,4]+2*vc[2,4])
  
  tab<-data.frame(
    group=c("A","B"),
    intercept=c(intercept_A, intercept_B),
    se_intercept=c(se_intercept_A,se_intercept_B),
    slope=c(slope_A,slope_B),
    se_slope=c(se_slope_A,se_slope_B)
  )
  tab[,-1]<-round(tab[,-1],rounded)
  
  return(tab)
}


################################################################## END OF SCRIPT
################################################################################


