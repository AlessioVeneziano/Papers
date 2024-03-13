################################################################################
########## Analysis: Size Dimorphism vs seasonality in primates - Veneziano 2024

### Load packages and functions
# The function 'plot.pcaCor' plots variable loadings onto the Principal Component
# scores. The function 'predict.kmeans' assigns new observations to the clusters
# identified by k-means algorithm based on their proximity to the centroids of
# said clusters.
library(terra)
library(ape)
library(RRphylo)
library(phytools)

plot.pcaCor<-function(pca,dat,n1=1,n2=2,len.factor=0.3,label=NULL,col.label="black",pos.text=4,off.text=0.5,...){
  f<-attributes(pca)$class
  if(f=="princomp"){
    scr<-pca$scores
  } else if (f=="prcomp"){
    scr<-pca$x
  } else (stop("Use 'prcomp' or 'princomp' outputs only!"))
  
  cors<-cor(scr,dat)
  
  x<-pca$loadings[,n1]*len.factor*sqrt(nrow(dat)-1)
  y<-pca$loadings[,n2]*len.factor*sqrt(nrow(dat)-1)
  arrows(0,0,x,y,...)
  
  if(missing(label)){
    text(x,y,labels=colnames(cors),col=col.label)
  } else {
    text(x,y,labels=label,col=col.label,pos=pos.text,offset=off.text)
  }
}
predict.kmeans<-function(km,newdata){
  apply(newdata,1,function(a) which.min(colSums((t(km$centers)-a)^2)))
}


### Read body weight, occurrence, climate and phylogenetic data
bsize<-read.table("data/body size/primates-body-weight_Veneziano2024.txt",header=T,sep="\t")

n.occ<-list.files("data/occurrence/occurrence",full.names=T)
  occ<-lapply(as.list(n.occ),read.table,sep="\t",header=T)

clim<-list.files("data/geo",full.names=T)
  clim<-rast(clim)

phy<-read.tree("data/consensus-tree_Veneziano2024")
  phy<-multi2di(phy)
  phy$tip.label<-sub("_"," ",phy$tip.label)


### Set useful factors and check sample size by family
spe<-bsize$Species
fam<-bsize$Family

table(fam)


### Get seasonality climatic variables at each occurrence site for each species 
amb<-lapply(occ,function(x) extract(clim,x[,1:2]))
  amb<-lapply(amb,unique)
    names(amb)<-spe
  amb<-do.call(rbind,amb)
  amb$species<-sub("\\..*","",rownames(amb))
    rownames(amb)<-NULL
  amb<-amb[,c(5,2:4)]

tablong<-table(amb$species)


### Prepare size variables
# The object 'sexd' represents sexual size dimorphism and it is computed as the
# logarithm of the male-to-female ratio.
mw<-bsize$MaleWg
  names(mw)<-bsize$Species

fw<-bsize$FemaleWg
  names(fw)<-bsize$Species

sexd<-log(mw/fw)
  names(sexd)<-spe


### PCA weighted on species sample size
samb<-scale(amb[,-1],center=T,scale=T)
wt<-rep(1/tablong,tablong)
wtcov<-cov.wt(samb,wt=wt,center=T)
pca<-princomp(samb,covmat=wtcov)

summary(pca)
pca$loadings
cor(pca$scores,samb)


### Identify species in highly seasonal environments
# Two clusters are identified based on PC1 and PC2. It's important to highlight
# that this method works because all the PC loadings represent an increase in
# seasonality along common directions (visible from the plots in the following
# section).
spelong<-rep(spe,c(tablong))
pcmedian<-aggregate(pca$scores[,1:2],list(spelong),median)

med<-apply(pca$scores[,1:2],2,median)
d<-sqrt(colSums((t(pca$scores[,1:2])-med)^2))
kdat<-scale(cbind(d,pca$scores[,1:2]))

set.seed(42)
km<-kmeans(kdat,2,iter.max=1e+05)

med<-apply(pcmedian[,-1],2,median)
d<-sqrt(colSums((t(pcmedian[,-1])-med)^2))
kdat<-scale(cbind(d,pcmedian[,-1]))
fac<-predict.kmeans(km,kdat)

facvar<-factor(fac,labels=c("seasonal","aseasonal"))
  names(facvar)<-spe


### Useful parameters for graphics
# The colour 'grey80' (light grey) is assigned to the seasonal group, while
# the colour 'grey60' (darker) is assigned to the aseasonal group.
col2<-c("grey80","grey60")
kmcol<-ifelse(km$cl==1,"grey80","grey60")
faccol<-ifelse(fac==1,"grey80","grey60")
facpch<-ifelse(fac==1,23,21)


### Plot PCA and seasonal species (Figure 1a)
plot(pca$scores[,1:2],pch=16,cex=0.3,col=kmcol,xlim=c(-2,8),ylim=c(-4,4),asp=1,
     xlab="PC1 (77.80% of Variance)",ylab="PC2 (18.65% of Variance)")
points(pcmedian[,-1],cex=1,pch=facpch,bg=faccol)
abline(h=0,v=0,lty=2,col="black")
plot.pcaCor(pca,samb,n1=1,n2=2,col="black",length=0.1,lwd=1,len.factor=0.02,
            label=c("BIO15","BIO4","BIO7"),col.label="black")


### Test differences in body size and sexual dimorphism vs seasonality
# This analysis is absent from the paper because it shows a result similar to
# other findings already present in literature. It shows that no difference in
# body weight dimorphism, or female/male body weight is present between seasonal
# and aseasonal species (as identified in this analysis).
phylANOVA(phy,facvar,sexd,nsim=1e+04)
phylANOVA(phy,facvar,log(fw),nsim=1e+04)
phylANOVA(phy,facvar,log(mw),nsim=1e+04)


### Compute rates via RRphylo
rr.d<-RRphylo(tree=phy,y=sexd)
rr.f<-RRphylo(tree=phy,y=log(fw))
rr.m<-RRphylo(tree=phy,y=log(mw))


### Test for rate shifts between seasonal and non-seasonal species via RRphylo
sh.d<-search.shift(rr.d,status.type="sparse",state=facvar,nrep=1e+05)
sh.f<-search.shift(rr.f,status.type="sparse",state=facvar,nrep=1e+05)
sh.m<-search.shift(rr.m,status.type="sparse",state=facvar,nrep=1e+05)

rbind(sh.d$state,
      c(sh.f$state$rate,1-sh.f$state$p),
      c(sh.m$state$rate,1-sh.m$state$p))


### Plot rate differences between seasonal and aseasonal groups (Figure 1b)
# The first few lines of the following code collect tip rates from the overall
# rates calculated by RRphylo (which computes also node rates).
keep<-grep("[A-Za-z]",rownames(rr.d$rates))
d.rate<-abs(rr.d$rates)[keep,]
  d.rate<-d.rate[match(spe,names(d.rate))]
f.rate<-abs(rr.f$rates)[keep,]
  f.rate<-f.rate[match(spe,names(f.rate))]
m.rate<-abs(rr.m$rates)[keep,]
  m.rate<-m.rate[match(spe,names(m.rate))]

boxfac<-paste(facvar,set=rep(1:3,each=length(facvar)))
  boxfac<-factor(boxfac,levels=c("seasonal 1","aseasonal 1",
                                 "seasonal 2","aseasonal 2",
                                 "seasonal 3","aseasonal 3"))
boxdf<-data.frame(rates=c(d.rate,f.rate,m.rate),boxfac)

boxplot(rates~boxfac,data=boxdf,at=c(1,2,4,5,7,8),outline=F,col=col2,
        xaxt="n",xlab="",ylab="Absolute rate")
  axis(1,las=1,at=c(1.5,4.5,7.5),labels=c("SSD","Female log-weight","Male log-weight"),tick=F)


### Compare female and male body size rate scaling in seasonal and aseasonal groups
# The regression tests for different intercepts and different slopes. Both Ordinary
# Least Regression (OLS) and Phylogenetic Generalized Least-Squares (PGLS) are
# performed.
df<-data.frame(x=f.rate,y=m.rate,fac=facvar)
mod<-lm(y~x*fac,data=df)
summary(mod)

bm<-corBrownian(1,phy,form=~spe|facvar)
mod1<-gls(m.rate~f.rate*facvar,correlation=bm)
summary(mod1)
mod1$coef[2]+mod1$coef[4]


### Plot rate scaling  (Figure 1c)
# This plot is based on the parameters estimated in the OLS because similar to
# those obtained via PGLS.
plot(df$x,df$y,pch=facpch,bg=faccol,xlab="Rate of female log-weight",
     ylab="Rate of male log-weight",xaxt="n",yaxt="n")
  axis(1,at=seq(0,0.3,0.05))
  axis(2,at=seq(0,0.3,0.05))
  abline(a=mod$coef[1],b=mod$coef[2],col=col2[1])
  abline(a=mod$coef[1]+mod$coef[3],b=mod$coef[2]+mod$coef[4],col=col2[2])
  abline(a=0,b=1,col="black",lty=2)



################################################################## END OF SCRIPT
################################################################################


