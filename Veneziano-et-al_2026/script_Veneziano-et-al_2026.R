################################################################################
#######               Veneziano, Alfieri & Riga, 2026                    #######
#####          Testing the effect of environmental seasonality             ##### 
###          on sexual size dimorphism and its phylogenetic rate             ###
##                   in catarrhines and platyrrhines                          ##
#                                                                              #
################################################################################

### Load packages and bespoke functions
library(ape)
library(RRphylo)
library(phytools)
library(nlme)

source("functions_Veneziano-et-al_2026.R")


### Read weight, climate and phylogeny data (and modify phylogenetic tip names to match those in "dat")
dat<-read.csv("data_Veneziano-et-al_2026.csv",header=T)

phy<-read.tree("consensus-tree-primates_Veneziano-et-al_2026")
  phy$tip.label<-sub("_"," ",phy$tip.label)


### Prepare factors for species and infraorder
spe<-dat$species
gro<-dat$infraorder


### Prepare variables for analysis
fw<-dat$Fweight # female weight
  fw<-log(fw)
  names(fw)<-spe
mw<-dat$Mweight # male weight
  mw<-log(mw)
  names(mw)<-spe
sexd<-log(dat$Mweight/dat$Fweight) # sexual size dimorphism
  names(sexd)<-spe

pseas<-dat$PrecSeas # precipitation seasonality
  names(pseas)<-spe
tseas<-dat$TempSeas # temperature seasonality
  names(tseas)<-spe


### Set colours and folder to save plots
cols<-c("orange","forestgreen")
colsA<-tapply(cols,1:2,colAlpha,alpha=0.5)

colsG<-cols[as.numeric(as.factor(gro))]
colsGA<-colsA[as.numeric(as.factor(gro))]

dir.create("images_Veneziano-et-al_2026")


### Compute phylogenetic rates for sexual size dimorphism, female and male log weights
rr.d<-RRphylo(tree=phy,y=sexd)
rr.f<-RRphylo(tree=phy,y=fw)
rr.m<-RRphylo(tree=phy,y=mw)

keep<-grep("[A-Za-z]",rownames(rr.d$rates)) # exclude rates computed at internal nodes
d.rate<-abs(rr.d$rates)[keep,]
  d.rate<-d.rate[match(spe,names(d.rate))]
f.rate<-abs(rr.f$rates)[keep,]
  f.rate<-f.rate[match(spe,names(f.rate))]
m.rate<-abs(rr.m$rates)[keep,]
  m.rate<-m.rate[match(spe,names(m.rate))]

ld.rate<-log(d.rate)
lf.rate<-log(f.rate)
lm.rate<-log(m.rate)


### Analysis 1: seasonality vs dimorphism, seasonality vs female-male weight, linear correlation
bm<-corPagel(1,phy,form=~spe,fixed=F)

mod.Td<-gls(sexd~tseas*gro,correlation=bm) # PGLS models for SSD, female and male weights
mod.Pd<-gls(sexd~pseas*gro,correlation=bm)
mod.Tf<-gls(fw~tseas*gro,correlation=bm)
mod.Pf<-gls(fw~pseas*gro,correlation=bm)
mod.Tm<-gls(mw~tseas*gro,correlation=bm)
mod.Pm<-gls(mw~pseas*gro,correlation=bm)

summary(mod.Td)
summary(mod.Pd)
summary(mod.Tf)
summary(mod.Pf)
summary(mod.Tm)
summary(mod.Pm)
  
modelEst2groups(mod.Td)
modelEst2groups(mod.Pd)
modelEst2groups(mod.Tf)
modelEst2groups(mod.Pf)
modelEst2groups(mod.Tm)
modelEst2groups(mod.Pm)

{
  pos<-gro=="Catarrhini"
  r2.cat<-c(cor(sexd[pos],predict(mod.Td)[pos])^2,
            cor(sexd[pos],predict(mod.Pd)[pos])^2,
            cor(fw[pos],predict(mod.Tf)[pos])^2,
            cor(fw[pos],predict(mod.Pf)[pos])^2,
            cor(mw[pos],predict(mod.Tm)[pos])^2,
            cor(mw[pos],predict(mod.Pm)[pos])^2)
  r2.pla<-c(cor(sexd[!pos],predict(mod.Td)[!pos])^2,
            cor(sexd[!pos],predict(mod.Pd)[!pos])^2,
            cor(fw[!pos],predict(mod.Tf)[!pos])^2,
            cor(fw[!pos],predict(mod.Pf)[!pos])^2,
            cor(mw[!pos],predict(mod.Tm)[!pos])^2,
            cor(mw[!pos],predict(mod.Pm)[!pos])^2)
  r2<-cbind(r2.cat,r2.pla)
  r2<-round(r2,3)
  rownames(r2)<-paste(rep(c("sexd","fw","mw"),each=2),rep(c("temp","prec"),3))
} # Pseudo R-squared for the PGLS models

{
  dpred<-data.frame(tseas=seq(3.3,6.5,length=5),pseas=seq(3,4.9,length=5),
                    gro=sort(rep(c("Catarrhini","Platyrrhini"),5)))
  pred.Td<-predict(mod.Td,dpred)
  pred.Pd<-predict(mod.Pd,dpred)
  
  fpred<-data.frame(tseas=seq(3.3,6.5,length=5),pseas=seq(3,4.9,length=5),
                    gro=sort(rep(c("Catarrhini","Platyrrhini"),5)))
  pred.Tf<-predict(mod.Tf,fpred)
  pred.Pf<-predict(mod.Pf,fpred)
  
  mpred<-data.frame(tseas=seq(3.3,6.5,length=5),pseas=seq(3,4.9,length=5),
                    gro=sort(rep(c("Catarrhini","Platyrrhini"),5)))
  pred.Tm<-predict(mod.Tm,mpred)
  pred.Pm<-predict(mod.Pm,mpred)
} # Prediction for fitting lines in Figure 1
{
  pdf("images_Veneziano-et-al_2026/Figure1.pdf",width=5.5,height=7.5)
  par(mfrow=c(3,2))
  
  axd<-seq(min(sexd),max(sexd),length=5)
  axf<-seq(min(fw),max(fw),length=5)
  axm<-seq(min(mw),max(mw),length=5)
  
  axt<-seq(min(tseas),max(tseas),length=5)
  axp<-seq(min(pseas),max(pseas),length=5)
  
  plot(sexd~tseas,col=colsGA,pch=16,xlab="Temperature seasonality",ylab="SSD",xaxt="n",yaxt="n")
    axis(1,axt,round(axt,1))
    axis(2,axd,round(axd,2),las=2)
    lines(dpred$tseas[1:5],pred.Td[1:5],col=cols[1],lwd=2)
    lines(dpred$tseas[6:10],pred.Td[6:10],col=cols[2],lwd=2)
  plot(sexd~pseas,col=colsGA,pch=16,xlab="Precipitation seasonality",ylab="SSD",xaxt="n",yaxt="n")
    axis(1,axp,round(axp,1))
    axis(2,axd,round(axd,2),las=2)
    lines(dpred$pseas[1:5],pred.Pd[1:5],col=cols[1],lwd=2)
    lines(dpred$pseas[6:10],pred.Pd[6:10],col=cols[2],lwd=2)
  plot(fw~tseas,col=colsGA,pch=16,xlab="Temperature seasonality",ylab="Female Size (log g)",xaxt="n",yaxt="n")
    axis(1,axt,round(axt,1))
    axis(2,axf,round(axf,1),las=2)
    lines(fpred$tseas[1:5],pred.Tf[1:5],col=cols[1],lwd=2)
    lines(fpred$tseas[6:10],pred.Tf[6:10],col=cols[2],lwd=2)
  plot(fw~pseas,col=colsGA,pch=16,xlab="Precipitation seasonality",ylab="Female Size (log g)",xaxt="n",yaxt="n")
    axis(1,axp,round(axp,1))
    axis(2,axf,round(axf,1),las=2)
    lines(fpred$pseas[1:5],pred.Pf[1:5],col=cols[1],lwd=2)
    lines(fpred$pseas[6:10],pred.Pf[6:10],col=cols[2],lwd=2)
  plot(mw~tseas,col=colsGA,pch=16,xlab="Temperature seasonality",ylab="Male Size (log g)",xaxt="n",yaxt="n")
    axis(1,axt,round(axt,1))
    axis(2,axm,round(axm,1),las=2)
    lines(mpred$tseas[1:5],pred.Tm[1:5],col=cols[1],lwd=2)
    lines(mpred$tseas[6:10],pred.Tm[6:10],col=cols[2],lwd=2)
  plot(mw~pseas,col=colsGA,pch=16,xlab="Precipitation seasonality",ylab="Male Size (log g)",xaxt="n",yaxt="n")
    axis(1,axp,round(axp,1))
    axis(2,axm,round(axm,1),las=2)
    lines(mpred$pseas[1:5],pred.Pm[1:5],col=cols[1],lwd=2)
    lines(mpred$pseas[6:10],pred.Pm[6:10],col=cols[2],lwd=2)
  dev.off()
} # Figure 1


### Analysis 2: seasonality vs rate of dimorphism, seasonality vs rate of female-male weight, linear
bm<-corPagel(1,phy,form=~spe,fixed=F)

mod.Tdr<-gls(ld.rate~tseas*gro,correlation=bm) # PGLS models for rate of SSD, female and male weights
mod.Pdr<-gls(ld.rate~pseas*gro,correlation=bm)
mod.Tfr<-gls(lf.rate~tseas*gro,correlation=bm)
mod.Pfr<-gls(lf.rate~pseas*gro,correlation=bm)
mod.Tmr<-gls(lm.rate~tseas*gro,correlation=bm)
mod.Pmr<-gls(lm.rate~pseas*gro,correlation=bm)

summary(mod.Tdr)
summary(mod.Pdr)
summary(mod.Tfr)
summary(mod.Pfr)
summary(mod.Tmr)
summary(mod.Pmr)

modelEst2groups(mod.Tdr)
modelEst2groups(mod.Pdr)
modelEst2groups(mod.Tfr)
modelEst2groups(mod.Pfr)
modelEst2groups(mod.Tmr)
modelEst2groups(mod.Pmr)

{
  pos<-gro=="Catarrhini"
  r2.cat<-c(cor(ld.rate[pos],predict(mod.Tdr)[pos])^2,
            cor(ld.rate[pos],predict(mod.Pdr)[pos])^2,
            cor(lf.rate[pos],predict(mod.Tfr)[pos])^2,
            cor(lf.rate[pos],predict(mod.Pfr)[pos])^2,
            cor(lm.rate[pos],predict(mod.Tmr)[pos])^2,
            cor(lm.rate[pos],predict(mod.Pmr)[pos])^2)
  r2.pla<-c(cor(ld.rate[!pos],predict(mod.Tdr)[!pos])^2,
            cor(ld.rate[!pos],predict(mod.Pdr)[!pos])^2,
            cor(lf.rate[!pos],predict(mod.Tfr)[!pos])^2,
            cor(lf.rate[!pos],predict(mod.Pfr)[!pos])^2,
            cor(lm.rate[!pos],predict(mod.Tmr)[!pos])^2,
            cor(lm.rate[!pos],predict(mod.Pmr)[!pos])^2)
  r2.rate<-cbind(r2.cat,r2.pla)
    r2.rate<-round(r2.rate,3)
  rownames(r2.rate)<-paste(rep(c("sexd","fw","mw"),each=2),rep(c("temp","prec"),3))
} # Pseudo R-squared for the PGLS models

{
  dpred<-data.frame(tseas=seq(3.3,6.5,length=5),pseas=seq(3,4.9,length=5),
                    gro=sort(rep(c("Catarrhini","Platyrrhini"),5)))
  pred.Tdr<-predict(mod.Tdr,dpred)
  pred.Pdr<-predict(mod.Pdr,dpred)
  
  fpred<-data.frame(tseas=seq(3.3,6.5,length=5),pseas=seq(3,4.9,length=5),
                    gro=sort(rep(c("Catarrhini","Platyrrhini"),5)))
  pred.Tfr<-predict(mod.Tfr,fpred)
  pred.Pfr<-predict(mod.Pfr,fpred)
  
  mpred<-data.frame(tseas=seq(3.3,6.5,length=5),pseas=seq(3,4.9,length=5),
                    gro=sort(rep(c("Catarrhini","Platyrrhini"),5)))
  pred.Tmr<-predict(mod.Tmr,mpred)
  pred.Pmr<-predict(mod.Pmr,mpred)
} # Prediction for fitting lines in Figure 2
{
  pdf("images_Veneziano-et-al_2026/Figure2.pdf",width=5.5,height=7.5)
  par(mfrow=c(3,2))
  
  axd<-seq(min(ld.rate),max(ld.rate),length=5)
  axf<-seq(min(lf.rate),max(lf.rate),length=5)
  axm<-seq(min(lm.rate),max(lm.rate),length=5)
  
  axt<-seq(min(tseas),max(tseas),length=5)
  axp<-seq(min(pseas),max(pseas),length=5)
  
  plot(ld.rate~tseas,col=colsGA,pch=16,xlab="Temperature seasonality",ylab="Phylo-rate SSD",xaxt="n",yaxt="n")
    axis(1,axt,round(axt,1))
    axis(2,axd,round(axd,1),las=2)
    lines(dpred$tseas[1:5],pred.Tdr[1:5],col=cols[1],lwd=2)
    lines(dpred$tseas[6:10],pred.Tdr[6:10],col=cols[2],lwd=2)
  plot(ld.rate~pseas,col=colsGA,pch=16,xlab="Precipitation seasonality",ylab="Phylo-rate SSD",xaxt="n",yaxt="n")
    axis(1,axp,round(axp,1))
    axis(2,axd,round(axd,1),las=2)
    lines(dpred$pseas[1:5],pred.Pdr[1:5],col=cols[1],lwd=2)
    lines(dpred$pseas[6:10],pred.Pdr[6:10],col=cols[2],lwd=2)
  plot(lf.rate~tseas,col=colsGA,pch=16,xlab="Temperature seasonality",ylab="Phylo-rate Female Size",xaxt="n",yaxt="n")
    axis(1,axt,round(axt,1))
    axis(2,axf,round(axf,1),las=2)
    lines(fpred$tseas[1:5],pred.Tfr[1:5],col=cols[1],lwd=2)
    lines(fpred$tseas[6:10],pred.Tfr[6:10],col=cols[2],lwd=2)
  plot(lf.rate~pseas,col=colsGA,pch=16,xlab="Precipitation seasonality",ylab="Phylo-rate Female Size",xaxt="n",yaxt="n")
    axis(1,axp,round(axp,1))
    axis(2,axf,round(axf,1),las=2)
    lines(fpred$pseas[1:5],pred.Pfr[1:5],col=cols[1],lwd=2)
    lines(fpred$pseas[6:10],pred.Pfr[6:10],col=cols[2],lwd=2)
  plot(lm.rate~tseas,col=colsGA,pch=16,xlab="Temperature seasonality",ylab="Phylo-rate Male Size",xaxt="n",yaxt="n")
    axis(1,axt,round(axt,1))
    axis(2,axm,round(axm,1),las=2)
    lines(mpred$tseas[1:5],pred.Tmr[1:5],col=cols[1],lwd=2)
    lines(mpred$tseas[6:10],pred.Tmr[6:10],col=cols[2],lwd=2)
  plot(lm.rate~pseas,col=colsGA,pch=16,xlab="Precipitation seasonality",ylab="Phylo-rate Male Size",xaxt="n",yaxt="n")
    axis(1,axp,round(axp,1))
    axis(2,axm,round(axm,1),las=2)
    lines(mpred$pseas[1:5],pred.Pmr[1:5],col=cols[1],lwd=2)
    lines(mpred$pseas[6:10],pred.Pmr[6:10],col=cols[2],lwd=2)
  dev.off()
} # Figure 2


### Define seasonal categories (i.e. find species in "high" seasonal conditions)
zp<-pseas/sd(pseas) # compute z-scores of climate variables without centering
zt<-tseas/sd(tseas)
pt<-sqrt(zp^2 + zt^2) # euclidean distance of z-scores from zero seasonality

thr<-quantile(pt,probs=c(0.50,0.70,0.90))

seas1<-ifelse(pt>=thr[1],"high","low") # categorical seasonality at sequential percentile thresholds
seas2<-ifelse(pt>=thr[2],"high","low")
seas3<-ifelse(pt>=thr[3],"high","low")

table(seas1,gro) # number of species at high/low seasonality by infraorder
table(seas2,gro)
table(seas3,gro)

spe[seas1=="high"] # species at high seasonality
spe[seas2=="high"]
spe[seas3=="high"]

{
  pdf("images_Veneziano-et-al_2026/Figure3a.pdf",width=7.5,height=5.0)
  plot(NA,xlim=c(7,15),ylim=c(-0.2,1.0),xlab="Distance from zero seasonality",ylab="SSD",yaxt="n")
  polygon(c(thr[1],thr[1],18,18),c(-2,2,2,-2),col=colAlpha("grey60",0.3),border=NA)
  polygon(c(thr[2],thr[2],18,18),c(-2,2,2,-2),col=colAlpha("grey60",0.3),border=NA)
  polygon(c(thr[3],thr[3],18,18),c(-2,2,2,-2),col=colAlpha("grey60",0.3),border=NA)
    points(pt,sexd,pch=21,bg=colsG)
  axis(2,seq(-0.2,1.0,length=5),las=2)
    text(thr[1]+0.2,-0.16,labels="PC50",col="grey30",cex=0.8,srt=90)
    text(thr[2]+0.2,-0.16,labels="PC70",col="grey30",cex=0.8,srt=90)
    text(thr[3]+0.2,-0.16,labels="PC90",col="grey30",cex=0.8,srt=90)
  dev.off()
} # Figure 3a (top)

{
  nam<-paste(substr(phy$tip.label,1,1),".",sub(".* ","",phy$tip.label),sep="")
  phy.n<-phy
    phy.n$tip.label<-nam
  sexd.n<-sexd
    names(sexd.n)<-paste(substr(names(sexd.n),1,1),".",sub(".* ","",names(sexd.n)),sep="")
  RR<-RRphylo(tree=phy.n,y=sexd.n)
  
  pdf("images_Veneziano-et-al_2026/Figure4.pdf",width=7,height=7)
  pRR<-plotRR(RR,y=sexd.n)
  colramp<-colorRampPalette(c("cyan","blue"))
  pRR$plotRRrates(tree.args=list(type="fan",edge.width=2,show.tip.label=T,cex=0.6,no.margin=T,label.offset=0.5),
                  colorbar.args=NULL,color.pal=colramp)
  dev.off()
} # Figure 4


### Analysis 3: seasonality vs dimorphism, categorical
set.seed(42)
pad1<-phylANOVA(phy,seas1,sexd,nsim=1e+04) # Phylogenetic ANOVA 50th percentile
set.seed(42)
paf1<-phylANOVA(phy,seas1,fw,nsim=1e+04)
set.seed(42)
pam1<-phylANOVA(phy,seas1,mw,nsim=1e+04)

set.seed(42)
pad2<-phylANOVA(phy,seas2,sexd,nsim=1e+04) # Phylogenetic ANOVA 70th percentile
set.seed(42)
paf2<-phylANOVA(phy,seas2,fw,nsim=1e+04)
set.seed(42)
pam2<-phylANOVA(phy,seas2,mw,nsim=1e+04)
  
set.seed(43)
pad3<-phylANOVA(phy,seas3,sexd,nsim=1e+04) # Phylogenetic ANOVA 90th percentile
set.seed(43)
paf3<-phylANOVA(phy,seas3,fw,nsim=1e+04)
set.seed(43)
pam3<-phylANOVA(phy,seas3,mw,nsim=1e+04)

{
  res<-rbind(round(c(pad1$'F',pad1$Pf),3),round(c(pad2$'F',pad2$Pf),3),round(c(pad3$'F',pad3$Pf),3),
             round(c(paf1$'F',paf1$Pf),3),round(c(paf2$'F',paf2$Pf),3),round(c(paf3$'F',paf3$Pf),3),
             round(c(pam1$'F',pam1$Pf),3),round(c(pam2$'F',pam2$Pf),3),round(c(pam3$'F',pam3$Pf),3))
  rownames(res)<-c(paste("sexd",c("50th","70th","90th"),sep="_"),
                   paste("female",c("50th","70th","90th"),sep="_"),
                   paste("male",c("50th","70th","90th"),sep="_"))
  res[,2]<-res[,2]*3 # Bonferroni correction
} # table of results of phylogenetic ANOVA


### Analysis 4: seasonality vs rate of dimorphism, categorical (threshold)
set.seed(42)
sh1<-search.shift(rr.d,status.type="sparse",state=seas1,nrep=1e+05) # Rate difference test, 50th percentile
set.seed(42)
shf1<-search.shift(rr.f,status.type="sparse",state=seas1,nrep=1e+05)
set.seed(42)
shm1<-search.shift(rr.m,status.type="sparse",state=seas1,nrep=1e+05)

set.seed(42)
sh2<-search.shift(rr.d,status.type="sparse",state=seas2,nrep=1e+05) # Rate difference test, 70th percentile
set.seed(42)
shf2<-search.shift(rr.f,status.type="sparse",state=seas2,nrep=1e+05)
set.seed(42)
shm2<-search.shift(rr.m,status.type="sparse",state=seas2,nrep=1e+05)

set.seed(42)
sh3<-search.shift(rr.d,status.type="sparse",state=seas3,nrep=1e+05) # Rate difference test, 90th percentile
set.seed(42)
shf3<-search.shift(rr.f,status.type="sparse",state=seas3,nrep=1e+05)
set.seed(42)
shm3<-search.shift(rr.m,status.type="sparse",state=seas3,nrep=1e+05)

{
  res<-rbind(round(rbind(sh1$state,sh2$state,sh3$state),3),
             round(rbind(shf1$state,shf2$state,shf3$state),3),
             round(rbind(shm1$state,shm2$state,shm3$state),3))
  rownames(res)<-c(paste("sexd",c("50th","70th","90th"),sep="_"),
                   paste("female",c("50th","70th","90th"),sep="_"),
                   paste("male",c("50th","70th","90th"),sep="_"))
  # change direction of difference and of p-value (see ?search.shift)
  for(i in 1:nrow(res)){
    if(sign(res[i,1])>=0){
      res[i,1]<-res[i,1]*-1
      res[i,2]<-1-res[i,2]
    }
  }
  res[,2]<-res[,2]*3 # Bonferroni correction
} # table of results of rate difference test


### Repeat analyses 3 and 4 for Catarrhini
pos<-gro=="Catarrhini"

sexd.cat<-sexd[pos]
fw.cat<-fw[pos]
mw.cat<-mw[pos]
spe.cat<-spe[pos]

pseas.cat<-pseas[pos]
tseas.cat<-tseas[pos]

phy.cat<-keep.tip(phy,spe[pos])

{
  zp.cat<-pseas.cat/sd(pseas.cat)
  zt.cat<-tseas.cat/sd(tseas.cat)
  pt.cat<-sqrt(zp.cat^2 + zt.cat^2)
  
  thr.cat<-quantile(pt.cat,probs=c(0.50,0.70,0.90))
  
  seas1.cat<-ifelse(pt.cat>=thr.cat[1],"high","low")
  seas2.cat<-ifelse(pt.cat>=thr.cat[2],"high","low")
  seas3.cat<-ifelse(pt.cat>=thr.cat[3],"high","low")
  
  table(seas1.cat)
  table(seas2.cat)
  table(seas3.cat)
  
  spe.cat[seas1.cat=="high"]
  spe.cat[seas2.cat=="high"]
  spe.cat[seas3.cat=="high"]
} # Define seasonal categories for Catarrhini
{
  pdf("images_Veneziano-et-al_2026/Figure3b.pdf",width=7.5,height=5.0)
  plot(NA,xlim=c(6,14),ylim=c(-0.2,1.0),xlab="Distance from zero seasonality",ylab="SSD",yaxt="n")
  polygon(c(thr.cat[1],thr.cat[1],18,18),c(-2,2,2,-2),col=colAlpha("grey60",0.3),border=NA)
  polygon(c(thr.cat[2],thr.cat[2],18,18),c(-2,2,2,-2),col=colAlpha("grey60",0.3),border=NA)
  polygon(c(thr.cat[3],thr.cat[3],18,18),c(-2,2,2,-2),col=colAlpha("grey60",0.3),border=NA)
  points(pt.cat,sexd.cat,pch=21,bg="orange")
  axis(2,seq(-0.2,1.0,length=5),las=2)
  text(thr.cat[1]+0.2,-0.16,labels="PC50",col="grey30",cex=0.8,srt=90)
  text(thr.cat[2]+0.2,-0.16,labels="PC70",col="grey30",cex=0.8,srt=90)
  text(thr.cat[3]+0.2,-0.16,labels="PC90",col="grey30",cex=0.8,srt=90)
  dev.off()
} # Figure 3b (centre)

{
  rr.d.cat<-RRphylo(tree=phy.cat,y=sexd.cat)
  rr.f.cat<-RRphylo(tree=phy.cat,y=fw.cat)
  rr.m.cat<-RRphylo(tree=phy.cat,y=mw.cat)
} # Compute phylogenetic rates

{
  set.seed(42)
  pad1.cat<-phylANOVA(phy.cat,seas1.cat,sexd.cat,nsim=1e+04)
  set.seed(42)
  paf1.cat<-phylANOVA(phy.cat,seas1.cat,fw.cat,nsim=1e+04)
  set.seed(42)
  pam1.cat<-phylANOVA(phy.cat,seas1.cat,mw.cat,nsim=1e+04)

  set.seed(42)
  pad2.cat<-phylANOVA(phy.cat,seas2.cat,sexd.cat,nsim=1e+04)
  set.seed(42)
  paf2.cat<-phylANOVA(phy.cat,seas2.cat,fw.cat,nsim=1e+04)
  set.seed(42)
  pam2.cat<-phylANOVA(phy.cat,seas2.cat,mw.cat,nsim=1e+04)

  set.seed(43)
  pad3.cat<-phylANOVA(phy.cat,seas3.cat,sexd.cat,nsim=1e+04)
  set.seed(43)
  paf3.cat<-phylANOVA(phy.cat,seas3.cat,fw.cat,nsim=1e+04)
  set.seed(43)
  pam3.cat<-phylANOVA(phy.cat,seas3.cat,mw.cat,nsim=1e+04)

  res<-rbind(round(c(pad1.cat$'F',pad1.cat$Pf),3),round(c(pad2.cat$'F',pad2.cat$Pf),3),round(c(pad3.cat$'F',pad3.cat$Pf),3),
             round(c(paf1.cat$'F',paf1.cat$Pf),3),round(c(paf2.cat$'F',paf2.cat$Pf),3),round(c(paf3.cat$'F',paf3.cat$Pf),3),
             round(c(pam1.cat$'F',pam1.cat$Pf),3),round(c(pam2.cat$'F',pam2.cat$Pf),3),round(c(pam3.cat$'F',pam3.cat$Pf),3))
  rownames(res)<-c(paste("sexd",c("50th","70th","90th"),sep="_"),
                   paste("female",c("50th","70th","90th"),sep="_"),
                   paste("male",c("50th","70th","90th"),sep="_"))
  res[,2]<-res[,2]*3 # Bonferroni correction
} # Analysis 3

{
  set.seed(42)
  sh1.cat<-search.shift(rr.d.cat,status.type="sparse",state=seas1.cat,nrep=1e+05)
  set.seed(42)
  shf1.cat<-search.shift(rr.f.cat,status.type="sparse",state=seas1.cat,nrep=1e+05)
  set.seed(42)
  shm1.cat<-search.shift(rr.m.cat,status.type="sparse",state=seas1.cat,nrep=1e+05)
  
  set.seed(42)
  sh2.cat<-search.shift(rr.d.cat,status.type="sparse",state=seas2.cat,nrep=1e+05)
  set.seed(42)
  shf2.cat<-search.shift(rr.f.cat,status.type="sparse",state=seas2.cat,nrep=1e+05)
  set.seed(42)
  shm2.cat<-search.shift(rr.m.cat,status.type="sparse",state=seas2.cat,nrep=1e+05)
  
  set.seed(42)
  sh3.cat<-search.shift(rr.d.cat,status.type="sparse",state=seas3.cat,nrep=1e+05)
  set.seed(42)
  shf3.cat<-search.shift(rr.f.cat,status.type="sparse",state=seas3.cat,nrep=1e+05)
  set.seed(42)
  shm3.cat<-search.shift(rr.m.cat,status.type="sparse",state=seas3.cat,nrep=1e+05)
  
  round(rbind(sh1.cat$state,shf1.cat$state,shm1.cat$state),3) # when rate.difference > 0, p = 1-p.status.diff
  round(rbind(sh2.cat$state,shf2.cat$state,shm2.cat$state),3) # when rate.difference > 0, p = 1-p.status.diff
  round(rbind(sh3.cat$state,shf3.cat$state,shm3.cat$state),3) # when rate.difference > 0, p = 1-p.status.diff
  
  res<-rbind(round(rbind(sh1.cat$state,sh2.cat$state,sh3.cat$state),3),
             round(rbind(shf1.cat$state,shf2.cat$state,shf3.cat$state),3),
             round(rbind(shm1.cat$state,shm2.cat$state,shm3.cat$state),3))
  rownames(res)<-c(paste("sexd",c("50th","70th","90th"),sep="_"),
                   paste("female",c("50th","70th","90th"),sep="_"),
                   paste("male",c("50th","70th","90th"),sep="_"))
  # change direction of difference and of p-value (see ?search.shift)
  for(i in 1:nrow(res)){
    if(sign(res[i,1])>=0){
      res[i,1]<-res[i,1]*-1
      res[i,2]<-1-res[i,2]
    }
  }
  res[,2]<-res[,2]*3 # Bonferroni correction
} # Analysis 4


### Repeat analyses 3 and 4 for Platyrrhini
pos<-gro=="Platyrrhini"

sexd.pla<-sexd[pos]
fw.pla<-fw[pos]
mw.pla<-mw[pos]
spe.pla<-spe[pos]

pseas.pla<-pseas[pos]
tseas.pla<-tseas[pos]

phy.pla<-keep.tip(phy,spe[pos])

{
  zp.pla<-pseas.pla/sd(pseas.pla)
  zt.pla<-tseas.pla/sd(tseas.pla)
  pt.pla<-sqrt(zp.pla^2 + zt.pla^2)
  
  thr.pla<-quantile(pt.pla,probs=c(0.50,0.70,0.90))
  
  seas1.pla<-ifelse(pt.pla>=thr.pla[1],"high","low")
  seas2.pla<-ifelse(pt.pla>=thr.pla[2],"high","low")
  seas3.pla<-ifelse(pt.pla>=thr.pla[3],"high","low")
  
  table(seas1.pla)
  table(seas2.pla)
  table(seas3.pla)
  
  spe.pla[seas1.pla=="high"]
  spe.pla[seas2.pla=="high"]
  spe.pla[seas3.pla=="high"]
} # Define seasonal categories for Platyrrhini
{
  pdf("images_Veneziano-et-al_2026/Figure3c.pdf",width=7.5,height=5.0)
  plot(NA,xlim=c(10,16),ylim=c(-0.2,1.0),xlab="Distance from zero seasonality",ylab="SSD",yaxt="n")
  polygon(c(thr.pla[1],thr.pla[1],18,18),c(-2,2,2,-2),col=colAlpha("grey60",0.3),border=NA)
  polygon(c(thr.pla[2],thr.pla[2],18,18),c(-2,2,2,-2),col=colAlpha("grey60",0.3),border=NA)
  polygon(c(thr.pla[3],thr.pla[3],18,18),c(-2,2,2,-2),col=colAlpha("grey60",0.3),border=NA)
  points(pt.pla,sexd.pla,pch=21,bg="forestgreen")
  axis(2,seq(-0.2,1.0,length=5),las=2)
  text(thr.pla[1]+0.2,-0.16,labels="PC50",col="grey30",cex=0.8,srt=90)
  text(thr.pla[2]+0.2,-0.16,labels="PC70",col="grey30",cex=0.8,srt=90)
  text(thr.pla[3]+0.2,-0.16,labels="PC90",col="grey30",cex=0.8,srt=90)
  dev.off()
} # Figure 3c (bottom)

{
  rr.d.pla<-RRphylo(tree=phy.pla,y=sexd.pla)
  rr.f.pla<-RRphylo(tree=phy.pla,y=fw.pla)
  rr.m.pla<-RRphylo(tree=phy.pla,y=mw.pla)
} # Compute phylogenetic rates

{
  set.seed(42)
  pad1.pla<-phylANOVA(phy.pla,seas1.pla,sexd.pla,nsim=1e+04)
  set.seed(42)
  paf1.pla<-phylANOVA(phy.pla,seas1.pla,fw.pla,nsim=1e+04)
  set.seed(42)
  pam1.pla<-phylANOVA(phy.pla,seas1.pla,mw.pla,nsim=1e+04)

  set.seed(42)
  pad2.pla<-phylANOVA(phy.pla,seas2.pla,sexd.pla,nsim=1e+04)
  set.seed(42)
  paf2.pla<-phylANOVA(phy.pla,seas2.pla,fw.pla,nsim=1e+04)
  set.seed(42)
  pam2.pla<-phylANOVA(phy.pla,seas2.pla,mw.pla,nsim=1e+04)

  set.seed(43)
  pad3.pla<-phylANOVA(phy.pla,seas3.pla,sexd.pla,nsim=1e+04)
  set.seed(43)
  paf3.pla<-phylANOVA(phy.pla,seas3.pla,fw.pla,nsim=1e+04)
  set.seed(43)
  pam3.pla<-phylANOVA(phy.pla,seas3.pla,mw.pla,nsim=1e+04)

  res<-rbind(round(c(pad1.pla$'F',pad1.pla$Pf),3),round(c(pad2.pla$'F',pad2.pla$Pf),3),round(c(pad3.pla$'F',pad3.pla$Pf),3),
             round(c(paf1.pla$'F',paf1.pla$Pf),3),round(c(paf2.pla$'F',paf2.pla$Pf),3),round(c(paf3.pla$'F',paf3.pla$Pf),3),
             round(c(pam1.pla$'F',pam1.pla$Pf),3),round(c(pam2.pla$'F',pam2.pla$Pf),3),round(c(pam3.pla$'F',pam3.pla$Pf),3))
  rownames(res)<-c(paste("sexd",c("50th","70th","90th"),sep="_"),
                   paste("female",c("50th","70th","90th"),sep="_"),
                   paste("male",c("50th","70th","90th"),sep="_"))
  res[,2]<-res[,2]*3 # Bonferroni correction
} # Analysis 3

{
  set.seed(42)
  sh1.pla<-search.shift(rr.d.pla,status.type="sparse",state=seas1.pla,nrep=1e+05)
  set.seed(42)
  shf1.pla<-search.shift(rr.f.pla,status.type="sparse",state=seas1.pla,nrep=1e+05)
  set.seed(42)
  shm1.pla<-search.shift(rr.m.pla,status.type="sparse",state=seas1.pla,nrep=1e+05)
  
  set.seed(42)
  sh2.pla<-search.shift(rr.d.pla,status.type="sparse",state=seas2.pla,nrep=1e+05)
  set.seed(42)
  shf2.pla<-search.shift(rr.f.pla,status.type="sparse",state=seas2.pla,nrep=1e+05)
  set.seed(42)
  shm2.pla<-search.shift(rr.m.pla,status.type="sparse",state=seas2.pla,nrep=1e+05)
  
  set.seed(42)
  sh3.pla<-search.shift(rr.d.pla,status.type="sparse",state=seas3.pla,nrep=1e+05)
  set.seed(42)
  shf3.pla<-search.shift(rr.f.pla,status.type="sparse",state=seas3.pla,nrep=1e+05)
  set.seed(42)
  shm3.pla<-search.shift(rr.m.pla,status.type="sparse",state=seas3.pla,nrep=1e+05)
  
  round(rbind(sh1.pla$state,shf1.pla$state,shm1.pla$state),3) # when rate.difference > 0, p = 1-p.status.diff
  round(rbind(sh2.pla$state,shf2.pla$state,shm2.pla$state),3) # when rate.difference > 0, p = 1-p.status.diff
  round(rbind(sh3.pla$state,shf3.pla$state,shm3.pla$state),3) # when rate.difference > 0, p = 1-p.status.diff
  
  res<-rbind(round(rbind(sh1.pla$state,sh2.pla$state,sh3.pla$state),3),
             round(rbind(shf1.pla$state,shf2.pla$state,shf3.pla$state),3),
             round(rbind(shm1.pla$state,shm2.pla$state,shm3.pla$state),3))
  rownames(res)<-c(paste("sexd",c("50th","70th","90th"),sep="_"),
                   paste("female",c("50th","70th","90th"),sep="_"),
                   paste("male",c("50th","70th","90th"),sep="_"))
  # change direction of difference and of p-value (see ?search.shift)
  for(i in 1:nrow(res)){
    if(sign(res[i,1])>=0){
      res[i,1]<-res[i,1]*-1
      res[i,2]<-1-res[i,2]
    }
  }
  res[,2]<-res[,2]*3 # Bonferroni correction
} # Analysis 4



################################################################## END OF SCRIPT
################################################################################


