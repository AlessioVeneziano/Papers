##########################################################
##################### Shamaei-Tousi et al 2022: FUNCTIONS

procDist<-function(m1,m2){
  sqrt(sum((m1-m2)^2))
}

splitAsym<-function(input,orig,mirr){
  O<-input[orig,]
  M<-input[mirr,]
  
  tot<-apply((O-M)^2,1,function(x) (sum(x))) #sqrt
  tot<-sum(tot)
  
  da<-(sum((colMeans(O)-colMeans(M))^2)) #sqrt
  fl<-tot-da
  
  out<-list(tot=tot,da=da,fl=fl)
  return(out)
}

groupHull<-function(x,y,group,cols=NULL,alpha=0.2,...){
  if(missing(cols)){cols<-palette()}
  
  xg<-split(x,group)
  yg<-split(y,group)
  
  for(i in 1:length(xg)){
    pos<-chull(xg[[i]],yg[[i]])
    pos<-c(pos,pos[1])
    
    acol<-col2rgb(cols[i])
    acol<-acol/255
    acol<-rgb(acol[1],acol[2],acol[3],alpha=alpha,maxColorValue=1)
    polygon(xg[[i]][pos],yg[[i]][pos],col=acol,...)
  }
}

revertPCA<-function(scores,loadings){
  scores%*%t(loadings)
}

dist2plane<-function(plane,target,absolute=F){
  x<-plane[,1]
  y<-plane[,2]
  z<-plane[,3]
  
  r<-lm(z~x+y)
  cc<-coef(r)[c(2,3,1)]
  
  num<-cc[1]*target[1] + cc[2]*target[2] - target[3] +cc[3]
  den<-sqrt(sum(c(cc[1:2],-1)^2))
  d<-num/den
  
  if(absolute==T){d<-abs(d)}
  
  return(distance=d)
}



############################################# END OF SCRIPT
##########################################################