library("igraph")
library("ggplot2")
library("reshape2")


#Necessary functions
erfun<-function(x){
  return(2 * pnorm(x * sqrt(2)) - 1)
}

alpha.gt<-function(a,b,omega,t){
  if(is.na(a)==TRUE | is.na(b)==TRUE){return(0)}
  alpha<-exp(-((a-b)^2/(omega)^2))
  denom<-(pi/2)*(erfun(t/omega)-erfun(-t/omega))
  alpha<-alpha/denom
  if(abs(a-b)>t){alpha<-0}
  return(alpha)           
}

network.sk<-function(m,geno,omega,t){
  A<-matrix(nrow=length(geno),ncol=length(geno))
  for(i1 in 1:length(geno)) {
    for(i2 in 1:length(geno)){
      A[i1,i2]<-alpha.gt(geno[i1],geno[i2],omega,t)
    }
  }
  m<-na.omit(m)
  mat<-matrix(ncol=nrow(m),nrow=nrow(m))
  for(j1 in 1:nrow(m)){
    for(j2 in (1:nrow(m))[-i]){
      mat1<-outer(m[j1,],m[j2,])
      mat[j1,j2]<-sum(mat1*A)
      
    }
  }
  list<-is.na(mat)
  mat[list]<-0
  diag(mat)<-0
  mat<-mat/max(mat)
  return(mat)
}


mod<-function(a,geno,omega,t){
  a<-na.omit(a)
  if(length(a)>0){
    net<-network.sk(a,geno,omega,t)
    G<-graph_from_adjacency_matrix(net,mode="undirected",weighted=TRUE,diag=FALSE)
    fgreedy<-fastgreedy.community(G,merges=TRUE, modularity=TRUE)
    nclump<-max(fgreedy$membership)
    modmax<-1-(1/nclump)
    return(max(fgreedy$modularity)/modmax)}
  else{return(NA)}
}

#Modularity-remove the species with extreme values from the calculation
mod2<-function(a,geno,omega,t){
  a<-na.omit(a)
  means<-apply(a,1,function(i1) sum(i1*geno))
  a1<-a[which(abs(means)<4.5),]
  if(length(a1)>0){
    net<-network.sk(a1,geno,omega,t)
    G<-graph_from_adjacency_matrix(net,mode="undirected",weighted=TRUE,diag=FALSE)
    fgreedy<-fastgreedy.community(G,merges=TRUE, modularity=TRUE)
    nclump<-max(fgreedy$membership)
    modmax<-1-(1/nclump)
    return(max(fgreedy$modularity)/modmax)}
  else{return(NA)}
}


#Mean Nearest neighbot distance- A measure of deviation from overdispersion pattern
mnnd<-function(a){
  a<-na.omit(a)
  if(length(a)>1){
    a<-a[order(a,decreasing=TRUE)]
    nnd<-c(length=length(a))
    nnd[1]<-a[1]-a[2]
    nnd[length(a)]<-a[length(a)-1]-a[length(a)]
    for(i in 2:(length(a)-1)){
      nnd[i]<-min((a[i-1]-a[i]),(a[i]-a[i+1]))
    }
    ranges<-a[1]-a[length(a)]
    mnnd<-sum(nnd)/length(nnd)
    mmax<-(ranges/(length(a)-1))
    return(mnnd/mmax)
  }else{return(NA)}
}

mnnd2<-function(a){
  a<-na.omit(a)
  a1<-a[which(abs(a)<4.5)]
  if(length(a1)>1){
    a1<-a1[order(a1,decreasing=TRUE)]
    nnd<-c(length=length(a1))
    nnd[1]<-a1[1]-a1[2]
    nnd[length(a1)]<-a1[length(a1)-1]-a1[length(a1)]
    for(i in 2:(length(a1)-1)){
      nnd[i]<-min((a1[i-1]-a1[i]),(a1[i]-a1[i+1]))
    }
    ranges<-a1[1]-a1[length(a1)]
    mnnd<-sum(nnd)/length(nnd)
    mmax<-(ranges/(length(a1)-1))
    return(mnnd/mmax)
  } else{return(NA)}
}


nsp<-20    #No. of species
n<-20     #No. of loci
geno<-seq(-5,5,l=2*n+1) #Vector of genotypes
nt<-length(geno)   # No. of genotypes
reps<-20 #Number of replicates for initial conditions

samples<-c(1,5*(1:10),100,200,500,ceiling(seq(1000,5000,length.out=25)))

#Load raw data
tot<-readRDS("SK_GT2020_r.Rds")

param<-matrix(ncol=3)
res<-matrix(ncol=nt*length(samples))
pops<-res
for(i in 1:length(tot)){
  temp.p<-tot[[i]][[1]]
  temp.r<-tot[[i]][[2]]
  temp.pop<-tot[[i]][[3]]
  param<-rbind(param,matrix(unlist(temp.p),ncol=3,byrow=TRUE))
  res<-rbind(res,do.call(rbind,temp.r))
  pops<-rbind(pops,do.call(rbind,temp.pop)[,-1])
}


param<-(param[-1,])
res<-res[-1,]
pops<-pops[-1,]


#Calculate trait means with time
meandat<-matrix(nrow=nrow(param),ncol=nsp*length(samples))
tps<-rep(1:length(samples),each=nt)
tsteps<-rep(1:length(samples),each=nsp)

for(i in 1:nrow(meandat)){
  for(j in 1:length(samples)){
  mat<-res[(nsp*(i-1)+1):(nsp*i),which(tps==j)]
  meandat[i,which(tsteps==j)]<-apply(mat,1,function(i3) sum(geno*i3))
  }
}


#Calculate species richness
tps<-rep(1:length(samples),each=nt)
tsteps<-rep(1:length(samples),each=nsp)
richness<-matrix(nrow=nrow(param),ncol=length(samples))
for(i in 1:nrow(richness)){
  for(j in 1:length(samples)){
  samp<-pops[(nsp*(i-1)+1):(nsp*i),which(tps==j)]
  richness[i,j]<-sum(rowSums(samp)>1)
  }
}

#Calculate modularity
tps<-rep(1:length(samples),each=nt)

modat2<-matrix(nrow=nrow(param),ncol=length(samples))
for(i in 1:nrow(param)){
  for(j in 1:length(samples)){
    mat<-res[(nsp*(i-1)+1):(nsp*i),which(tps==j)]
    mat<-mat/rowSums(mat)
    modat2[i,j]<-mod2(mat,geno,param[i,2],param[i,1])
  }
}

modat<-matrix(nrow=nrow(param),ncol=length(samples))
for(i in 1:nrow(param)){
  for(j in 1:length(samples)){
    mat<-pops[(nsp*(i-1)+1):(nsp*i),which(tps==j)]
    mat<-mat/rowSums(mat)
    modat[i,j]<-mod(mat,geno,param[i,2],param[i,1])
  }
}


#Calculate MNND
tsteps<-rep(1:length(samples),each=nsp)
mnnds<-matrix(nrow=nrow(meandat),ncol=length(samples))
for(i in 1:nrow(meandat)){
  for(j in 1:length(samples)){
    samp<-meandat[i,which(tsteps==j)]
    samp<-samp[which(samp!=0)]
    mnnds[i,j]<-mnnd(samp)
  }
}

mnnds2<-matrix(nrow=nrow(meandat),ncol=length(samples))
for(i in 1:nrow(meandat)){
  for(j in 1:length(samples)){
    samp<-meandat[i,which(tsteps==j)]
    samp<-samp[which(samp!=0)]
    mnnds2[i,j]<-mnnd2(samp)
  }
}

#Plot modularity vs. time
modtimes<-modat
modtimes<-cbind(param,modtimes)
colnames(modtimes)<-c("t","Omega","Tvar",samples)
modtimes<-as.data.frame(modtimes)
modtimes$t<-as.factor(modtimes$t)
modtimes$Tvar<-as.factor(modtimes$Tvar)
modtimes$Omega<-as.factor(modtimes$Omega)

plot1dat<-melt(modtimes,id=c("t","Tvar","Omega"))
plot1dat$variable<-as.numeric(plot1dat$variable)
dumb<-function(i){return(samples[i])}
plot1dat$variable<-sapply(plot1dat$variable,dumb)

plot1dat<-na.omit(plot1dat)



plot1datm<-aggregate(plot1dat$value,list(plot1dat$t,plot1dat$Tvar,plot1dat$Omega,plot1dat$variable),mean,na.rm=TRUE)
plot1dats<-aggregate(plot1dat$value,list(plot1dat$t,plot1dat$Tvar,plot1dat$Omega,plot1dat$variable),sd,na.rm=TRUE)

plot1dat1<-cbind(plot1datm,plot1datm$x+plot1dats$x,plot1datm$x-plot1dats$x)

colnames(plot1dat1)<-c("t","Intra","omega","Time","mean","upper","lower")

plot1dat1<-plot1dat1[which(plot1dat1$t%in%c(0.5,100) & plot1dat1$Intra%in%c(0.2,0.9)),]

plot1dat1<-plot1dat1[which(plot1dat1$t%in%c(0.5,100)),]

t.labs<-cbind(levels(plot1dat1$t),c("Strong Threshold","Intermediate threshold",
                                    "Weak Threshold","No Threshold"))
idx.t<-match(plot1dat1$t,t.labs[,1])
plot1dat1$t1<-t.labs[idx.t,2]

var.labs<-cbind(levels(plot1dat1$Intra),c("Weak Competition",
                                          "Intermediate Competition","Strong Competition"))
idx.var<-match(plot1dat1$Intra,var.labs[,1])
plot1dat1$Intra1<-var.labs[idx.var,2]

o.labs<-cbind(levels(plot1dat1$omega),c("Narrow","Intermediate","Wide"))

idx.o<-match(plot1dat1$omega,o.labs[,1])
plot1dat1$omega1<-o.labs[idx.o,2]


plot1dat1$Intra1<-factor(plot1dat1$Intra1,levels=c("Weak Competition","Intermediate Competition",
                                                   "Strong Competition"))

plot1dat1$t1<-factor(plot1dat1$t1,levels=c("No Threshold","Weak Threshold",
                                           "Intermediate threshold","Strong Threshold"))

plot1dat1$omega1<-factor(plot1dat1$omega1,levels=c("Narrow","Intermediate","Wide"))


plot1<-ggplot(plot1dat1,aes(x=Time,y=mean, group=omega1,colour=omega1)) +geom_line(size=1.5) +ylim(0,1)+
  geom_ribbon(aes(ymin=lower, ymax=upper, x=Time, fill=omega1),alpha=0.2)+ylab("Modularity")+
  facet_grid(Intra1~t1,labeller=labeller(Intra1=label_wrap_gen(10),t1=label_wrap_gen(10)),switch="y")+
  labs(fill="Competition Width",color="Competition Width")+
  theme(
    strip.text.x = element_text(
      size = 11, face = "bold"
    ),
    strip.text.y = element_text(
      size = 10, face = "bold"
    ),
    legend.title = element_text(color = "black", size = 15),
    legend.key.size = unit(1.2, "cm"),
    legend.text=element_text(size=15),
    axis.title=element_text(size=15,face="bold")
  )
plot1

  


tiff('SKGT2020_short_Mod2.tiff', units="in", width=10, height=5, res=300)
plot1
dev.off()






#plot MNND vs. time
mnnds1<-cbind(param,mnnds)
colnames(mnnds1)<-c("t","Omega","Intra",samples)
mnnds1<-as.data.frame(mnnds1)
mnnds1$t<-as.factor(mnnds1$t)
mnnds1$Intra<-as.factor(mnnds1$Intra)
mnnds1$Omega<-as.factor(mnnds1$Omega)


plot2dat<-melt(mnnds1,id=c("t", "Omega","Intra"))
plot2dat$variable<-as.numeric(plot2dat$variable)
dumb<-function(i){return(samples[i])}
plot2dat$variable<-sapply(plot2dat$variable,dumb)

plot2dat<-na.omit(plot2dat)

plot2datm<-aggregate(plot2dat$value,list(plot2dat$t,plot2dat$Omega,plot2dat$Intra,plot2dat$variable),mean,na.rm=TRUE)
plot2dats<-aggregate(plot2dat$value,list(plot2dat$t,plot2dat$Omega,plot2dat$Intra,plot2dat$variable),sd,na.rm=TRUE)
plot2dat1<-cbind(plot2datm,plot2datm$x+plot2dats$x,plot2datm$x-plot2dats$x)
colnames(plot2dat1)<-c("t","omega","Intra","Time","mean","upper","lower")

plot2dat1<-plot2dat1[which(plot2dat1$t%in%c(0.5,100) & plot2dat1$Intra%in%c(0.2,0.9)),]

t.labs<-cbind(levels(plot2dat1$t),c("Strong Threshold","Intermediate threshold",
                                    "Weak Threshold","No Threshold"))
idx.t<-match(plot2dat1$t,t.labs[,1])
plot2dat1$t1<-t.labs[idx.t,2]

var.labs<-cbind(levels(plot2dat1$Intra),c("Weak Competition",
                                          "Intermediate Competition","Strong Competition"))
idx.var<-match(plot2dat1$Intra,var.labs[,1])
plot2dat1$Intra1<-var.labs[idx.var,2]

o.labs<-cbind(levels(plot2dat1$omega),c("Narrow","Intermediate","Wide"))

idx.o<-match(plot2dat1$omega,o.labs[,1])
plot2dat1$omega1<-o.labs[idx.o,2]


plot2dat1$Intra1<-factor(plot2dat1$Intra1,levels=c("Weak Competition","Intermediate Competition",
                                                   "Strong Competition"))

plot2dat1$t1<-factor(plot2dat1$t1,levels=c("No Threshold","Weak Threshold",
                                           "Intermediate threshold","Strong Threshold"))

plot2dat1$omega1<-factor(plot2dat1$omega1,levels=c("Narrow","Intermediate","Wide"))

plot2<-ggplot(plot2dat1,aes(x=Time,y=mean, group=omega1,colour=omega1)) +geom_line(size=1.5) + ylim(0,1) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=Time, fill=omega1),alpha=0.2)+ylab("Mean nearest neighbor distance")+
  facet_grid(Intra1~t1,labeller=labeller(Intra1=label_wrap_gen(10),t1=label_wrap_gen(10)),switch="y")+
  labs(fill="Competition Width",color="Competition Width")+
  theme(
    strip.text.x = element_text(
      size = 11, face = "bold"
    ),
    strip.text.y = element_text(
      size = 10, face = "bold"
    ),
    legend.title = element_text(color = "black", size = 15),
    legend.key.size = unit(1.2, "cm"),
    legend.text=element_text(size=15),
    axis.title=element_text(size=15,face="bold")
  )
plot2



tiff('SKGT2020_short_MNND.tiff', units="in", width=10, height=5, res=300)
plot2
dev.off()


################################################
#Plot richness vs. time
rich2<-cbind(param,richness)
colnames(rich2)<-c("t","Omega","inter",samples)
rich2<-as.data.frame(rich2)
rich2$t<-as.factor(rich2$t)
rich2$inter<-as.factor(rich2$inter)
rich2$Omega<-as.factor(rich2$Omega)

plot3dat<-melt(rich2,id=c("t","Omega","inter"))
plot3dat$variable<-as.numeric(plot3dat$variable)
dumb<-function(i){return(samples[i])}
plot3dat$variable<-sapply(plot3dat$variable,dumb)

plot3datm<-aggregate(plot3dat$value,list(plot3dat$t,plot3dat$Omega,plot3dat$inter,plot3dat$variable),mean,na.rm=TRUE)
plot3dats<-aggregate(plot3dat$value,list(plot3dat$t,plot3dat$Omega,plot3dat$inter,plot3dat$variable),sd,na.rm=TRUE)
plot3dat1<-cbind(plot3datm,plot3datm$x+plot3dats$x,plot3datm$x-plot3dats$x)
colnames(plot3dat1)<-c("t","omega","Intra","Time","mean","upper","lower")

plot3dat1<-plot3dat1[which(plot3dat1$t%in%c(0.5,100) & plot3dat1$Intra%in%c(0.2,0.9) & plot3dat1$Time<5001),]

t.labs<-cbind(levels(plot3dat1$t),c("Strong Threshold","Intermediate threshold",
                                    "Weak Threshold","No Threshold"))
idx.t<-match(plot3dat1$t,t.labs[,1])
plot3dat1$t1<-t.labs[idx.t,2]

var.labs<-cbind(levels(plot3dat1$Intra),c("Weak Competition",
                                          "Intermediate Competition","Strong Competition"))
idx.var<-match(plot3dat1$Intra,var.labs[,1])
plot3dat1$Intra1<-var.labs[idx.var,2]

o.labs<-cbind(levels(plot3dat1$omega),c("Narrow","Intermediate","Wide"))

idx.o<-match(plot3dat1$omega,o.labs[,1])
plot3dat1$omega1<-o.labs[idx.o,2]


plot3dat1$Intra1<-factor(plot3dat1$Intra1,levels=c("Weak Competition","Intermediate Competition",
                                                   "Strong Competition"))

plot3dat1$t1<-factor(plot3dat1$t1,levels=c("No Threshold","Weak Threshold",
                                           "Intermediate threshold","Strong Threshold"))

plot3dat1$omega1<-factor(plot3dat1$omega1,levels=c("Narrow","Intermediate","Wide"))



plot3<-ggplot(plot3dat1,aes(x=Time,y=mean, group=omega1,colour=omega1)) +geom_line(size=1.5) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, x=Time, fill=omega1),alpha=0.2)+ylab("Species richness")+
  facet_grid(Intra1~t1,labeller=labeller(Intra1=label_wrap_gen(10),t1=label_wrap_gen(10)),switch="y")+
  labs(fill="Competition Width",color="Competition Width")+
  theme(
    strip.text.x = element_text(
      size = 11, face = "bold"
    ),
    strip.text.y = element_text(
      size = 10, face = "bold"
    ),
    legend.title = element_text(color = "black", size = 15),
    legend.key.size = unit(1.2, "cm"),
    legend.text=element_text(size=15),
    axis.title=element_text(size=15,face="bold")
  )
plot3



tiff('SKGT2020_Short_Rich.tiff', units="in", width=10, height=5, res=300)
plot3
dev.off()



#Compare initial vs. final trait distributions
initdat<-res[,1:nt]
findat<-res[,(ncol(res)+1-nt):ncol(res)]
shap.ini<-t(apply(initdat,1,shaptest,geno=geno))
shap.fin<-t(apply(findat,1,shaptest,geno=geno))

plot(density(shap.fin[,1]),xlim=c(0,1),ylim=c(0,20),xlab="Shapiro test statistic",main="")
lines(density(shap.ini[,1]),lty=2)

sum(which(shap.fin[,2]>0.05))

param1<-param[rep((1:nrow(param)),each=nsp),]
iniplot<-cbind(param1,shap.fin[,1])
colnames(iniplot)<-c("t","omega","inter","shap")
iniplot<-as.data.frame(iniplot)
iniplot$t<-as.factor(iniplot$t)
iniplot$omega<-as.factor(iniplot$omega)
iniplot$inter<-as.factor(iniplot$inter)

inibox<-ggplot(iniplot,aes(x=t,y=shap,fill=omega))+geom_boxplot(position=position_dodge(1))+
  facet_grid(~inter)                                                        
inibox


iniplot.m<-aggregate(iniplot$shap,list(iniplot$t,iniplot$omega,iniplot$inter),mean,na.rm=TRUE)
iniplot.s<-aggregate(iniplot$shap,list(iniplot$t,iniplot$omega,iniplot$inter),sd,na.rm=TRUE)
iniplot1<-cbind(iniplot.m,iniplot.m$x+iniplot.s$x,iniplot.m$x-iniplot.s$x)
colnames(iniplot1)<-c("t","omega","inter","mean","upper","lower")

var.labs<-cbind(levels(iniplot1$inter),c("Weak Competition",
                                         "Intermediate Competition","Strong Competition"))
idx.var<-match(iniplot1$inter,var.labs[,1])
iniplot1$inter1<-var.labs[idx.var,2]

iniplot1$inter1<-factor(iniplot1$inter1,levels=c("Weak Competition","Intermediate Competition",
                                                 "Strong Competition"))

plot.ini<-ggplot(iniplot1,aes(x=t,y=mean,group=omega,fill=omega))+ xlab("Threshold")+ylab("Shapiro Statistic")+
  geom_bar(stat="identity",position=position_dodge())+geom_errorbar(aes(ymin=lower,ymax=upper),width=0.2,
                                                                    position=position_dodge(.9))+ facet_grid(~inter1,labeller=labeller(Intra1=label_wrap_gen(10)))

plot.ini

#Plot the final modularity
modat.d<-modat[,length(samples)]-modat[,1]
modfin<-cbind(param,modat.d)
colnames(modfin)<-c("t","omega","inter","change")
modfin<-as.data.frame(modfin)
modfin$t<-as.factor(modfin$t)
modfin$inter<-as.factor(modfin$inter)
modfin$omega<-as.factor(modfin$omega)

modfin.m<-aggregate(modfin$change,list(modfin$t,modfin$omega,modfin$inter),mean,na.rm=TRUE)
modfin.s<-aggregate(modfin$change,list(modfin$t,modfin$omega,modfin$inter),sd,na.rm=TRUE)

plot12dat<-cbind(modfin.m,modfin.s$x)
colnames(plot12dat)<-c("t","omega","inter","mean","sd")

var.labs<-cbind(levels(plot12dat$inter),c("Weak Competition",
                                          "Intermediate Competition","Strong Competition"))
idx.var<-match(plot12dat$inter,var.labs[,1])
plot12dat$inter1<-var.labs[idx.var,2]
plot12dat$inter1<-factor(plot12dat$inter1,levels=c("Weak Competition","Intermediate Competition",
                                                   "Strong Competition"))

plot12<-ggplot(plot12dat,aes(x=t,y=mean,fill=omega))+ylab("Change in Modularity")+xlab("Threshold")+ylim(-1,1)+
  geom_bar(stat="identity",position=position_dodge())+geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.2,
                                                                    position=position_dodge(.9))+ facet_grid(~inter1,labeller=labeller(Intra1=label_wrap_gen(10)))
#+  theme_classic()
plot12
tiff('SKGT2020_final_finmod.tiff', units="in", width=10, height=5, res=300)
plot12
dev.off()