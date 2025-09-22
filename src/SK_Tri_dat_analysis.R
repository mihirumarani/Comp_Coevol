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

alpha.tri<-function(a,b,t){
  if(is.na(a)==TRUE | is.na(b)==TRUE){return(0)}
  if(abs(a-b)>=t){return(alpha<-0)}
  slope<-1/t
  if(((a-b)>-t) & ((a-b)<=0)){alpha<-(1+(slope*(a-b)))}
  if(((a-b)>0) & ((a-b)<t)){alpha<-(1-(slope*(a-b)))}
  return(alpha)
}

network.sk<-function(m,geno,t){
  A<-matrix(nrow=length(geno),ncol=length(geno))
  for(i1 in 1:length(geno)) {
    for(i2 in 1:length(geno)){
      A[i1,i2]<-alpha.tri(geno[i1],geno[i2],t)
    }
  }

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


mod<-function(a,geno,t){
  a<-na.omit(a)
  if(length(a)>1){
    net<-network.sk(a,geno,t)
    G<-graph_from_adjacency_matrix(net,mode="undirected",weighted=TRUE,diag=FALSE)
    fgreedy<-fastgreedy.community(G,merges=TRUE, modularity=TRUE)
    return(max(fgreedy$modularity))}
  else{return(NA)}
}



#Nearest neighbot distance
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
  }
  else{return(NA)}
}


#Parameter values
nsp<-20    #No. of species
n<-20     #No. of loci
geno<-seq(-5,5,l=2*n+1) #Vector of genotypes
nt<-length(geno)   # No. of genotypes
reps<-20 #Number of replicates for initial conditions

samples<-c(1,5*(1:10),100,200,500,ceiling(seq(1000,5000,length.out=25)))

#Load raw data
tot<-readRDS("SK_Tri2020.Rds")
param<-matrix(ncol=2)
res<-matrix(ncol=nt*length(samples))
pops<-res
for(i in 1:length(tot)){
  temp.p<-tot[[i]][[1]]
  temp.r<-tot[[i]][[2]]
  temp.pop<-tot[[i]][[3]]
  param<-rbind(param,matrix(unlist(temp.p),ncol=2,byrow=TRUE))
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

modat<-matrix(nrow=nrow(param),ncol=length(samples))
for(i in 1:nrow(param)){
  for(j in 1:length(samples)){
    mat<-res[(nsp*(i-1)+1):(nsp*i),which(tps==j)]
    modat[i,j]<-mod(mat,geno,param[i,1])
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



#Plot modularity vs. time
modtimes<-modat
modtimes<-cbind(param,modtimes)
colnames(modtimes)<-c("t","Tvar",1:length(samples))
modtimes<-as.data.frame(modtimes)
modtimes$t<-as.factor(modtimes$t)
modtimes$Tvar<-as.factor(modtimes$Tvar)

plot1dat<-melt(modtimes,id=c("t","Tvar"))
plot1dat$t<-as.factor(plot1dat$t)
plot1dat$Tvar<-as.factor(plot1dat$Tvar)
plot1dat$variable<-as.numeric(plot1dat$variable)
dumb<-function(i){return(samples[i])}
plot1dat$variable<-sapply(plot1dat$variable,dumb)

plot1datm<-aggregate(plot1dat$value,list(plot1dat$t,plot1dat$Tvar,plot1dat$variable),mean,na.rm=TRUE)

plot1dats<-aggregate(plot1dat$value,list(plot1dat$t,plot1dat$Tvar,plot1dat$variable),sd,na.rm=TRUE)

plot1dat1<-cbind(plot1datm,plot1datm$x+plot1dats$x,plot1datm$x-plot1dats$x)
colnames(plot1dat1)<-c("t","Intra","Time","mean","upper","lower")


t.labs<-cbind(levels(plot1dat1$t),c("Strong Threshold","Intermediate threshold",
                                    "Weak Threshold","No Threshold"))
idx.t<-match(plot1dat1$t,t.labs[,1])
plot1dat1$t1<-t.labs[idx.t,2]

var.labs<-cbind(levels(plot1dat1$Intra),c("Weak Competition",
                                          "Intermediate Competition","Strong Competition"))
idx.var<-match(plot1dat1$Intra,var.labs[,1])
plot1dat1$Intra1<-var.labs[idx.var,2]

plot1dat1$Intra1<-factor(plot1dat1$Intra1,levels=c("Weak Competition","Intermediate Competition",
                                                   "Strong Competition"))

plot1dat1$t1<-factor(plot1dat1$t1,levels=c("No Threshold","Weak Threshold",
                                           "Intermediate threshold","Strong Threshold"))



plot1<-ggplot(plot1dat1,aes(x=Time,y=mean)) +geom_line(size=1.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=Time),alpha=0.2)+ylab("Modularity") + 
  facet_grid(Intra1~t1,labeller=labeller(Intra1=label_wrap_gen(10),t1=label_wrap_gen(10)),switch="y")+
  theme(
    strip.text.x = element_text(
      size = 15, face = "bold"
    ),
    strip.text.y = element_text(
      size = 12, face = "bold"
    ),
    legend.key.size = unit(1.5, "cm"),
    legend.text=element_text(size=15),
    axis.title=element_text(size=20,face="bold")
  )
plot1

tiff('SKTri2020_final_Mod.tiff', units="in", width=10, height=5, res=300)
plot1
dev.off()


#plot MNND vs. time
mnnds1<-cbind(param,mnnds)
colnames(mnnds1)<-c("t","inter",samples)
mnnds1<-as.data.frame(mnnds1)
mnnds1$t<-as.factor(mnnds1$t)
mnnds1$inter<-as.factor(mnnds1$inter)

plot2dat<-melt(mnnds1,id=c("t","inter"))
plot2dat$variable<-as.numeric(plot2dat$variable)
dumb<-function(i){return(samples[i])}
plot2dat$variable<-sapply(plot2dat$variable,dumb)
plot2dat<-na.omit(plot2dat)

plot2datm<-aggregate(plot2dat$value,list(plot2dat$t,plot2dat$inter,plot2dat$variable),mean,na.rm=TRUE)
plot2dats<-aggregate(plot2dat$value,list(plot2dat$t,plot2dat$inter,plot2dat$variable),sd,na.rm=TRUE)
plot2dat1<-cbind(plot2datm,plot2datm$x+plot2dats$x,plot2datm$x-plot2dats$x)
colnames(plot2dat1)<-c("t","inter","Time","mean","upper","lower")


t.labs<-cbind(levels(plot2dat1$t),c("Strong Threshold","Intermediate threshold",
                                    "Weak Threshold","No Threshold"))
idx.t<-match(plot2dat1$t,t.labs[,1])
plot2dat1$t1<-t.labs[idx.t,2]

var.labs<-cbind(levels(plot2dat$inter),c("Weak Competition",
                                         "Intermediate Competition","Strong Competition"))
idx.var<-match(plot2dat1$inter,var.labs[,1])
plot2dat1$inter1<-var.labs[idx.var,2]

plot2dat1$inter1<-factor(plot2dat1$inter1,levels=c("Weak Competition","Intermediate Competition",
                                                   "Strong Competition"))

plot2dat1$t1<-factor(plot2dat1$t1,levels=c("No Threshold","Weak Threshold",
                                           "Intermediate threshold","Strong Threshold"))



plot2<-ggplot(plot2dat1,aes(x=Time,y=mean)) +geom_line(size=1.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=Time),alpha=0.2)+ylab("Mean nearest neighbor distance")+
  facet_grid(inter1~t1,labeller=labeller(Intra1=label_wrap_gen(10),t1=label_wrap_gen(10)),switch="y")+
  theme(
    strip.text.x = element_text(
      size = 13, face = "bold"
    ),
    strip.text.y = element_text(
      size = 13, face = "bold"
    ),
    legend.key.size = unit(1.5, "cm"),
    legend.text=element_text(size=15),
    axis.title=element_text(size=20,face="bold")
  )
plot2

tiff('SKTri2020_final_MNND.tiff', units="in", width=10, height=5, res=300)
plot2
dev.off()


################################################
#Plot richness vs. time

rich2<-cbind(param,richness)
colnames(rich2)<-c("t","inter",samples)
rich2<-as.data.frame(rich2)
rich2$t<-as.factor(rich2$t)
rich2$inter<-as.factor(rich2$inter)

plot3dat<-melt(rich2,id=c("t","inter"))
plot3dat$variable<-as.numeric(plot3dat$variable)
dumb<-function(i){return(samples[i])}
plot3dat$variable<-sapply(plot3dat$variable,dumb)

plot3dat<-na.omit(plot3dat)

plot3datm<-aggregate(plot3dat$value,list(plot3dat$t,plot3dat$inter,plot3dat$variable),mean,na.rm=TRUE)
plot3dats<-aggregate(plot3dat$value,list(plot3dat$t,plot3dat$inter,plot3dat$variable),sd,na.rm=TRUE)
plot3dat1<-cbind(plot3datm,plot3datm$x+plot3dats$x,plot3datm$x-plot3dats$x)
colnames(plot3dat1)<-c("t","inter","Time","mean","upper","lower")


t.labs<-cbind(levels(plot3dat1$t),c("Strong Threshold","Intermediate threshold",
                                    "Weak Threshold","No Threshold"))
idx.t<-match(plot3dat1$t,t.labs[,1])
plot3dat1$t1<-t.labs[idx.t,2]

var.labs<-cbind(levels(plot3dat$inter),c("Weak Competition",
                                         "Intermediate Competition","Strong Competition"))
idx.var<-match(plot3dat1$inter,var.labs[,1])
plot3dat1$inter1<-var.labs[idx.var,2]

plot3dat1$inter1<-factor(plot3dat1$inter1,levels=c("Weak Competition","Intermediate Competition",
                                                   "Strong Competition"))

plot3dat1$t1<-factor(plot3dat1$t1,levels=c("No Threshold","Weak Threshold",
                                           "Intermediate threshold","Strong Threshold"))

plot3<-ggplot(plot3dat1,aes(x=Time,y=mean)) +geom_line(size=1.5) + 
  geom_ribbon(aes(ymin=lower, ymax=upper, x=Time),alpha=0.2)+ylab("Species richness")+
  facet_grid(inter1~t1,labeller=labeller(Intra1=label_wrap_gen(10),t1=label_wrap_gen(10)),switch="y")+
  theme(
    strip.text.x = element_text(
      size = 11, face = "bold"
    ),
    strip.text.y = element_text(
      size = 10, face = "bold"
    ),
    legend.key.size = unit(1.5, "cm"),
    legend.text=element_text(size=15),
    axis.title=element_text(size=15,face="bold")
  )
plot3


tiff('SKTri2020_final_Rich.tiff', units="in", width=10, height=5, res=300)
plot3
dev.off()






