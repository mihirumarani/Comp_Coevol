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
  if(abs(a-b)>t){return(alpha<-0)}
  slope<-1/t
  if(((a-b)>-t) & ((a-b)<=0)){alpha<-(1+(slope*(a-b)))}
  if(((a-b)>0) & ((a-b)<t)){alpha<-(1-(slope*(a-b)))}
  return(alpha)
}

network<-function(a,t){
  a<-na.omit(a)
  mat<-matrix(ncol=length(a),nrow=length(a))
  samples<-1:length(a)
  for(i in samples){
    for(j in samples){
      
      mat[i,j]<-alpha.tri(a[i],a[j],t)
      
    }
  }
  list<-is.na(mat)
  mat[list]<-0
  diag(mat)<-0
  mat<-mat/max(mat)
  return(mat)
}

mods<-function(a,t){
  a<-na.omit(a)
  if(length(a)>1){
    net<-network(a,omega,t)
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
nsp<-10    #No. of species
n<-20     #No. of loci
geno<-seq(-5,5,l=2*n+1) #Vector of genotypes
nt<-length(geno)   # No. of genotypes
reps<-20 #Number of replicates for initial conditions

#Parameters to vary
ts<-c(0.5,1,1.5,100)                           #Threshold to competition kernel
a1s<-c(0.2,0.5,0.9)                            #Strength of interspecific competition relative to intraspecific comp.

samples<-c(1,5*(1:10),100,200,500,ceiling(seq(1000,5000,length.out=25)))

#Load raw data
tot<-readRDS("SK_Tri1020_55tr.Rds")
param<-matrix(ncol=2)
res<-matrix(ncol=nt*length(samples))
for(i in 1:length(tot)){
  temp.p<-tot[[i]][[1]]
  temp.r<-tot[[i]][[2]]
  temp.pop<-tot[[i]][[3]]
  param<-rbind(param,matrix(unlist(temp.p),ncol=2,byrow=TRUE))
  res<-rbind(res,do.call(rbind,temp.r))
}
param<-(param[-1,])
res<-res[-1,]

#Calculate trait means with time
meandat<-matrix(nrow=nrow(param),ncol=nsp*length(samples))
tps<-rep(1:length(samples),each=nt)

for(i in 1:nrow(meandat)){
  mat<-res[(nsp*(i-1)+1):(nsp*i),]
  mat.l<-lapply(1:length(samples),function(j) mat[,which(tps==j)])
  means<-c()
  for(k in 1:length(mat.l)){
    means<-c(means,apply(mat.l[[k]],1,function(x) sum(x*geno)))
  }
  meandat[i,]<-means
}

#Calculate trait variances with time
vardat<-matrix(nrow=length(param),ncol=nsp*length(samples))
tps<-rep(1:length(samples),each=nt)
for(i in 1:nrow(vardat)){
  mat<-res[(nsp*(i-1)+1):(nsp*i),]
  mat.l<-lapply(1:length(samples),function(j) mat[,which(tps==j)])
  vars<-c()
  for(j in 1:length(mat.l)){
   mean<-apply(mat.l[[j]],1,function(x) sum(x*geno))
   for(k in 1:nrow(mat.l[[j]])){
     temp<-mat.l[[j]][k,]
     var<-sum(temp*((geno-mean[k])^2))
     vars<-c(vars,var)
   }
  }
  vardat[i,]<-vars
}


#Calculate modularity
tsteps<-rep(1:length(samples),each=nsp)
modat<-matrix(nrow=nrow(meandat),ncol=length(samples))
for(i in 1:nrow(meandat)){
  for(j in 1:length(samples)){
    modat[i,j]<-mods(meandat[i,which(tsteps==j)],param[i])
  }
}



#Calculate MNND
tsteps<-rep(1:length(samples),each=nsp)
mnnds<-matrix(nrow=nrow(meandat),ncol=length(samples))
for(i in 1:nrow(meandat)){
  for(j in 1:length(samples)){
    samp<-meandat[i,which(tsteps==j)]
    samp<-samp[-which(samp==0)]
    mnnds[i,j]<-mnnd(samp)
  }
}



#Plot modularity vs. time
modtimes<-modat
modtimes<-cbind(param,modtimes)
colnames(modtimes)<-c("t","Tvar","Omega",1:length(samples))
modtimes<-as.data.frame(modtimes)
modtimes$t<-as.factor(modtimes$t)
modtimes$Tvar<-as.factor(modtimes$Tvar)
modtimes$Omega<-as.factor(modtimes$Omega)

plot1dat<-melt(modtimes,id=c("t","Tvar","Omega"))
plot1dat$t<-as.factor(plot1dat$t)
plot1dat$Tvar<-as.factor(plot1dat$Tvar)
plot1dat$Omega<-as.factor(plot1dat$Omega)
plot1dat$variable<-as.numeric(plot1dat$variable)
dumb<-function(i){return(samples[i])}
plot1dat$variable<-sapply(plot1dat$variable,dumb)

plot1datm<-aggregate(plot1dat$value,list(plot1dat$t,plot1dat$Tvar,plot1dat$Omega,plot1dat$variable),mean,na.rm=TRUE)

plot1dats<-aggregate(plot1dat$value,list(plot1dat$t,plot1dat$Tvar,plot1dat$Omega,plot1dat$variable),sd,na.rm=TRUE)

plot1dat1<-cbind(plot1datm,plot1datm$x+plot1dats$x,plot1datm$x-plot1dats$x)
colnames(plot1dat1)<-c("t","Intra_var","omega","Time","mean","upper","lower")

plot1dat1$t<-factor(plot1dat1$t,levels=c("100","1.5","1","0.5"))
plot1dat1$Intra_var<-factor(plot1dat1$Intra_var,levels=c("1","2"))
plot1dat1$omega<-factor(plot1dat1$omega,levels=c("0.1", "0.5"))


var.labs<-c("Low ITV","High ITV")
names(var.labs)<-c("1","2")


t.labs<-c("Threshold=0.5","Threshold=1","Threshold=1.5","No Threshold")
names(t.labs)<-c("0.5","1","1.5","100")

plot1<-ggplot(plot1dat1,aes(x=Time,y=mean, group=omega,colour=omega)) +geom_line(size=1.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=Time, fill=omega),alpha=0.2)+ylab("Modularity") + 
  facet_grid(Intra_var~t,labeller=labeller(Intra_var=var.labs,t=t.labs),switch="y")+
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


plot2dat1$t<-factor(plot2dat1$t,levels=c("100","5","2.5","1"))
plot2dat1$inter<-factor(plot2dat1$inter,levels=c("0.2","0.5","0.9"))


t.labs<-c("Thr=1","Thr=2.5","Thr=5","No Thr")
names(t.labs)<-c("1","2.5","5","100")

i.labs<-c("Weak competition","Intermediate competition","Strong competition")
names(i.labs)<-c("0.2","0.5","0.9")


plot2<-ggplot(plot2dat1,aes(x=Time,y=mean)) +geom_line(size=1.5,col="red") +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=Time,colour="red",fill="red"),alpha=0.2)+ylab("Mean nearest neighbor distance")+
  facet_grid(t~inter,labeller=labeller(t=t.labs,inter=i.labs),switch="y")+
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


################################################
#Plot extinctions

tps<-rep(1:length(samples),each=nt)
exdat<-matrix(nrow=nrow(param),ncol=nsp*length(samples))

for(i in 1:nrow(exdat)){
  mat<-res[(nsp*(i-1)+1):(nsp*i),]
  mat.l<-lapply(1:length(samples),function(j) mat[,which(tps==j)])
  means<-c()
  for(k in 1:length(mat.l)){
    means<-c(means,apply(mat.l[[k]],1,function(x) if(sum(x)==0){return(0)} else{return(1)}))
  }
  exdat[i,]<-means
}

exdat1<-matrix(nrow=nrow(exdat),ncol=length(samples))
tsteps<-rep(1:length(samples),each=nsp)
for(i in 1:nrow(exdat)){
  for(j in 1:length(samples)){
    exdat1[i,j]<-sum(exdat[i,which(tsteps==j)])
  }
}


exdat2<-cbind(param,exdat1)
colnames(exdat2)<-c("t","inter",samples)
exdat2<-as.data.frame(exdat2)
exdat2$t<-as.factor(exdat2$t)
exdat2$inter<-as.factor(exdat2$inter)

plot3dat<-melt(exdat2,id=c("t","inter"))
plot3dat$variable<-as.numeric(plot3dat$variable)
dumb<-function(i){return(samples[i])}
plot3dat$variable<-sapply(plot3dat$variable,dumb)

plot3datm<-aggregate(plot3dat$value,list(plot3dat$t,plot3dat$inter,plot3dat$variable),mean,na.rm=TRUE)
plot3dats<-aggregate(plot3dat$value,list(plot3dat$t,plot3dat$inter,plot3dat$variable),sd,na.rm=TRUE)
plot3dat1<-cbind(plot3datm,plot3datm$x+plot3dats$x,plot3datm$x-plot3dats$x)
colnames(plot3dat1)<-c("t","inter","Time","mean","upper","lower")


plot3dat1$t<-factor(plot3dat1$t,levels=c("100","5","2.5","1"))
plot3dat1$inter<-factor(plot3dat1$inter,levels=c("0.2","0.5","0.9"))


t.labs<-c("Thr=1","Thr=2.5","Thr=5","No Thr")
names(t.labs)<-c("1","2.5","5","100")

i.labs<-c("Weak competition","Intermediate competition","Strong competition")
names(i.labs)<-c("0.2","0.5","0.9")


plot3<-ggplot(plot3dat1,aes(x=Time,y=mean)) +geom_line(size=1.5,col="red") +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=Time,colour="red",fill="red"),alpha=0.2)+ylab("Mean nearest neighbor distance")+
  facet_grid(t~inter,labeller=labeller(t=t.labs,inter=i.labs),switch="y")+
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
plot3





