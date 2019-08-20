library("igraph")
library("ggplot2")
library("reshape")


#Necessary functions

network<-function(a,omega,t){
  a<-na.omit(a)
  mat<-matrix(ncol=length(a),nrow=length(a))
  samples<-1:length(a)
  for(i in samples){
    for(j in samples){
      
      mat[i,j]<-alpha1(a[i],a[j],omega,t)
      
    }
  }
  list<-is.na(mat)
  mat[list]<-0
  diag(mat)<-0
  mat<-mat/max(mat)
  return(mat)
}

mod<-function(a,omega,t){
  a<-na.omit(a)
  net<-network(a,omega,t)
  G<-graph_from_adjacency_matrix(net,mode="undirected",weighted=TRUE,diag=FALSE)
  fgreedy<-fastgreedy.community(G,merges=TRUE, modularity=TRUE)
  return(max(fgreedy$modularity))
}




#Raw data:
paramdat<-read.csv("paramdatsfull.csv")
mnnds<-read.csv("mnndsfull.csv")
popdat<-read.csv("popdatsfull.csv")
trdat<-read.csv("trdatsfull.csv")



trdat<-as.matrix(trdat[,-1])
mnnds<-mnnds<-as.matrix(mnnds[,-1])
popdat<-as.matrix(popdat[,-1])




dumb2<-function(a){
  i<-a[1]
  j<-a[2]
  if(i+j==2){return(1)}
  if(i+j==4){return(3)}
  if(i+j==6){return(2)}
}
newvar<-apply(paramdat[,3:4],1,dumb2)
paramdat<-cbind(paramdat[,2],newvar,paramdat[,5])


dumb2<-function(a){
  if(a==0.5){return(1)}
  if(a==2.0){return(2)}
}


newvar<-apply(paramdat[,3:4],1,dumb2)
paramdat<-cbind(paramdat[,2],newvar,paramdat[,5])

index<-cbind(c(1,2,3),c("low","high","wide"))



modsamples<-ceiling(c(seq(1,2000,length.out=20),seq(2100,49999,length.out=30)))


trlist<-list()
seq1<-1+20*(0:49)
for(i in 1:nrow(trdat)){
  res2<-lapply(seq1,function(j,vec) vec[j:(j+19)], vec=trdat[i,])
  res2<-do.call(rbind,res2)
  trlist[[i]]<-res2
}

poplist<-list()
seq1<-1+20*(0:49)
for(i in 1:nrow(popdat)){
  res2<-lapply(seq1,function(j,vec) vec[j:(j+19)], vec=popdat[i,])
  res2<-do.call(rbind,res2)
  poplist[[i]]<-res2
}


#Calculate modularity
modat<-matrix(nrow=nrow(trdat),ncol=50)
for(i in 1:length(trlist)){
  for(j in 1:50){
    modat[i,j]<-mod(trlist[[i]][j,],paramdat[i,3],paramdat[i,1])
  }
}



#Species richness vs. time
richness<-matrix(nrow=nrow(popdat),ncol=50)
for(i in 1:length(poplist)){
  for(j in 1:50){
    richness[i,j]<-length(which(poplist[[i]][j,]<=0))
  }
}
richness<-n-richness
richdat<-cbind(paramdat,richness)
colnames(richdat)<-c("t","Tvar","Omega",1:50)
richdat<-as.data.frame(richdat)
richdat$t<-as.factor(richdat$t)
richdat$Tvar<-as.factor(richdat$Tvar)
richdat$Omega<-as.factor(richdat$Omega)
richdat2<-melt(richdat,id=c("t","Tvar", "Omega"))
richdat2$variable<-as.numeric(richdat2$variable)
dumb<-function(i){return(modsamples[i])}
richdat2$variable<-sapply(richdat2$variable,dumb)
dumb3<-function(i){return(index[i,2])}
richdat2$Tvar<-sapply(richdat2$Tvar,dumb3)
richdat.m<-aggregate(richdat2$value,list(richdat2$t,richdat2$Tvar,richdat2$Omega,richdat2$variable),mean,na.rm=TRUE)
richdat.s<-aggregate(richdat2$value,list(richdat2$t,richdat2$Tvar,richdat2$Omega,richdat2$variable),sd,na.rm=TRUE)
richdata<-cbind(richdat.m,richdat.m$x+richdat.s$x,richdat.m$x-richdat.s$x)
colnames(richdata)<-c("t","Tvar","omega","time","mean","upper","lower")
richdata$lower[which(richdata$lower<=0)]<-0

richdata$Tvar<-factor(richdata$Tvar,levels=c("low","high","wide"),labels=c("U[0.5,1,5]","U[2.5,3.5]","U[0.5,3.5]"))

richdata$t<-factor(richdata$t,labels=c("5","7","10","No threshold"))

richnplot<-ggplot(richdata,aes(x=time,y=mean,color=omega))+geom_line(size=1.5)+ylab("Species Richness")+
 facet_grid(Tvar~t,labeller=label_both)
richplot

geom_ribbon(aes(ymin=lower, ymax=upper, x=time, fill=omega),alpha=0.2)
#Plot connectance vs. time


conndat<-matrix(nrow=nrow(trdat),ncol=50)
times<-1+(0:49)*20
for(i in 1:nrow(trdat)){
  for(j in 1:50){
    dat<-na.omit(trdat[i,(times[j]:(times[j]+19))])
    net<-network(dat,paramdat[i,3],paramdat[i,1])
    n<-nrow(net)
    conn<-length(which(net!=0))/(n*(n-1))
    conndat[i,j]<-conn
  }
}




conndat<-cbind(paramdat,conndat)
colnames(conndat)<-c("t","Tvar","Omega",1:50)
conndat<-as.data.frame(conndat)
conndat$t<-as.factor(conndat$t)
conndat$Tvar<-as.factor(conndat$Tvar)
conndat$Omega<-as.factor(conndat$Omega)
conndat2<-melt(conndat,id=c("t","Tvar", "Omega"))
conndat2$variable<-as.numeric(conndat2$variable)
dumb<-function(i){return(modsamples[i])}
conndat2$variable<-sapply(conndat2$variable,dumb)
dumb3<-function(i){return(index[i,2])}
conndat2$Tvar<-sapply(conndat2$Tvar,dumb3)
conndat.m<-aggregate(conndat2$value,list(conndat2$t,conndat2$Tvar,conndat2$Omega,conndat2$variable),mean,na.rm=TRUE)
conndat.s<-aggregate(conndat2$value,list(conndat2$t,conndat2$Tvar,conndat2$Omega,conndat2$variable),sd,na.rm=TRUE)
conndata<-cbind(conndat.m,conndat.m$x+conndat.s$x,conndat.m$x-conndat.s$x)
colnames(conndata)<-c("t","Tvar","omega","time","mean","upper","lower")
conndata$lower[which(conndata$lower<=0)]<-0

conndata$Tvar<-factor(conndata$Tvar,levels=c("low","high","wide"),labels=c("U[0.5,1,5]","U[2.5,3.5]","U[0.5,3.5]"))

conndata$t<-factor(conndata$t,labels=c("5","7","10","No threshold"))
  



connplot<-ggplot(conndata,aes(x=time,y=mean,color=omega))+geom_line(size=1.5)+ylab("connectance")+
  geom_ribbon(aes(ymin=lower, ymax=upper, x=time, fill=omega),alpha=0.2)+facet_grid(Tvar~t,labeller=label_both)
tiff('plot11.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')

connplot  
dev.off()

#Modularity index vs. time and how it differs across param space

modtimes<-modat
modtimes<-cbind(paramdat,modtimes)
colnames(modtimes)<-c("t","Tvar","Omega",1:50)
modtimes<-as.data.frame(modtimes)
modtimes$t<-as.factor(modtimes$t)
modtimes$Tvar<-as.factor(modtimes$Tvar)
modtimes$Omega<-as.factor(modtimes$Omega)


#Plot no. 1: Modularity vs. time 
modtimes1<-modtimes
plot1dat<-melt(modtimes1,id=c("t","Tvar","Omega"))
plot1dat$variable<-as.numeric(plot1dat$variable)
dumb<-function(i){return(modsamples[i])}
plot1dat$variable<-sapply(plot1dat$variable,dumb)
dumb3<-function(i){return(index[i,2])}
plot1dat$Tvar<-sapply(plot1dat$Tvar,dumb3)

plot1datm<-aggregate(plot1dat$value,list(plot1dat$t,plot1dat$Tvar,plot1dat$Omega,plot1dat$variable),mean,na.rm=TRUE)

plot1dats<-aggregate(plot1dat$value,list(plot1dat$t,plot1dat$Tvar,plot1dat$Omega,plot1dat$variable),sd,na.rm=TRUE)

plot1dat<-cbind(plot1datm,plot1datm$x+plot1dats$x,plot1datm$x-plot1dats$x)
colnames(plot1dat)<-c("t","Tvar","omega","time","mean","upper","lower")


plot1dat$Tvar<-factor(plot1dat$Tvar,levels=c("low","high","wide"),labels=c("U[0.5,1,5]","U[2.5,3.5]","U[0.5,3.5]"))

plot1dat$t<-factor(plot1dat$t,labels=c("5","7","10","No threshold"))


plot1<-ggplot(plot1dat,aes(x=time,y=mean, group=omega,colour=omega)) +geom_line(size=1.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=time, fill=omega),alpha=0.2)+ylab("Modularity") + facet_grid(Tvar~t,labeller=label_both)

tiff('plot12.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')
plot1
dev.off()


#Plot data with Mod vs. time for threshold data wrt omega and Tvar.

#resample data points where initial modularity for different omegas is similar
#Filter the values with means and sd closer to overall mean and sd


initmean<-mean(modtimes$`1`)
initsd<-sd(modtimes$`1`)
range<-c(initmean-(0.5*initsd),initmean+(0.5*initsd))
dat05<-which(modtimes$Omega==0.05 & findInterval(modtimes$`1`,range)==1)
dat10<-which(modtimes$Omega==0.1 & findInterval(modtimes$`1`,range)==1)
dat25<-which(modtimes$Omega==0.25 & findInterval(modtimes$`1`,range)==1)
dat50<-which(modtimes$Omega==0.5 & findInterval(modtimes$`1`,range)==1)
dat100<-which(modtimes$Omega==1 & findInterval(modtimes$`1`,range)==1)


omegadat<-rbind(modtimes[dat05,],modtimes[dat10,],modtimes[dat25,],modtimes[dat50,],modtimes[dat100,])

omegadat1<-melt(omegadat,id=c("t","Tvar", "Omega"))
omegadat1$variable<-as.numeric(omegadat1$variable)

dumb<-function(i){return(modsamples[i])}
omegadat1$variable<-sapply(omegadat1$variable,dumb)
dumb3<-function(i){return(index[i,2])}
omegadat1$Tvar<-sapply(omegadat1$Tvar,dumb3)
omegameans.t<-aggregate(omegadat1$value,list(omegadat1$t,omegadat1$Omega,omegadat1$Tvar,omegadat1$variable),mean,na.rm=TRUE)
omegasd.t<-aggregate(omegadat1$value,list(omegadat1$t,omegadat1$Omega,omegadat1$Tvar,omegadat1$variable),sd,na.rm=TRUE)

lower.o<-omegameans.t$x-omegasd.t$x
lower.o[which(lower.o<=0)]<-0


omegadats<-cbind(omegameans.t,omegameans.t$x+omegasd.t$x,lower.o)
colnames(omegadats)<-c("t","Omega","Tvar","time","mean","upper","lower")



pomega<-ggplot(omegadats,aes(x=time,y=mean,colour=Omega)) +geom_line(size=1.5)  +
  geom_ribbon(aes(ymin=lower, ymax=upper,fill=Omega),alpha=0.2)+ facet_wrap(~t+Tvar)

pomega


#Fractions of simulations showing increased modularity

modchange<-modtimes$`50`-modtimes$`1`
df<-cbind(modtimes[,1:5],modchange=modchange)

df.t<-df[which(df$Threshold!=200),]
df.nt<-df[which(df$Threshold==200),]
colnames(df.t)<-c("Threshold","variance1","variance 2","omega","h","Change.in.Modularity")
colnames(df.nt)<-c("Threshold","variance1","variance 2","omega","h","Change.in.Modularity")


stat2<-glm(df$modchange~df$t + df$Tvar1 + df$Omega+ df$h, family=gaussian)
summary(stat2)

boxplot(df.t$Change.in.Modularity~df.t$omega,xlab="Omega",ylab="Change in Modularity")
boxplot(df.t$Change.in.Modularity~df.t$h,xlab="Heritability",ylab="Change in Modularity")
boxplot(df.t$Change.in.Modularity~df.t$Threshold,xlab="Threshold",ylab="Change in Modularity")
boxplot(df.t$Change.in.Modularity~df.t$variance1,xlab="Intraspecific variance",ylab="Change in Modularity")

boxplot(df.nt$Change.in.Modularity~df$omega,xlab="Omega",ylab="Change in Modularity")
boxplot(df.nt$Change.in.Modularity~df$h,xlab="Heritability",ylab="Change in Modularity")
boxplot(df.nt$Change.in.Modularity~df$Threshold,xlab="Threshold",ylab="Change in Modularity")
boxplot(df.nt$Change.in.Modularity~df$variance1,xlab="Intraspecific variance",ylab="Change in Modularity")



#MNND plots


mnnds<-cbind(paramdat,mnnds)
colnames(mnnds)<-c("t","Tvar","Omega",modsamples)
mnnds<-as.data.frame(mnnds)
mnnds$t<-as.factor(mnnds$t)
mnnds$Tvar<-as.factor(mnnds$Tvar)
mnnds$Omega<-as.factor(mnnds$Omega)


mnnds1<-mnnds


plot2dat<-melt(mnnds1,id=c("t","Tvar", "Omega"))
plot2dat$variable<-as.numeric(plot2dat$variable)
dumb<-function(i){return(modsamples[i])}
plot2dat$variable<-sapply(plot2dat$variable,dumb)
dumb3<-function(i){return(index[i,2])}
plot2dat$Tvar<-sapply(plot2dat$Tvar,dumb3)
plot2datm<-aggregate(plot2dat$value,list(plot2dat$t,plot2dat$Tvar,plot2dat$Omega,plot2dat$variable),mean,na.rm=TRUE)
plot2dats<-aggregate(plot2dat$value,list(plot2dat$t,plot2dat$Tvar,plot2dat$Omega,plot2dat$variable),sd,na.rm=TRUE)
plot2dat<-cbind(plot2datm,plot2datm$x+plot2dats$x,plot2datm$x-plot2dats$x)
colnames(plot2dat)<-c("t","Tvar","omega","time","mean","upper","lower")


plot2dat$Tvar<-factor(plot2dat$Tvar,levels=c("low","high","wide"),labels=c("U[0.5,1,5]","U[2.5,3.5]","U[0.5,3.5]"))

plot2dat$t<-factor(plot2dat$t,labels=c("5","7","10","No threshold"))

plot2dat1<-transform(plot2dat,Tvar=factor(Tvar,levels=c("low","high","wide")))


plot2<-ggplot(plot2dat,aes(x=time,y=mean, group=omega,colour=omega)) +geom_line(size=1.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=time, fill=omega),alpha=0.2)+ylab("Mean nearest neighbor distance")+facet_grid(Tvar~t,labeller=label_both)

tiff('plot13.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')
plot2
dev.off()






##Null model: Keep the initial demographic parameters the same and run them through the simulation.


paramdat1<-read.csv("paramdatsnull.csv")
mnnds1<-read.csv("mnndsnull.csv")
popdat1<-read.csv("popdatsnull.csv")
trdat1<-read.csv("trdatsnull.csv")
trshift1<-read.csv("trshiftnull.csv")

trdat1<-as.matrix(trdat1[,-1])
mnnds1<-mnnds<-as.matrix(mnnds1[,-1])
popdat1<-as.matrix(popdat1[,-1])
trshift1<-as.matrix(trshift1[,-1])



trlist1<-list()
for(i in 1:nrow(trdat1)){
  res2<-lapply(seq1,function(j,vec) vec[j:(j+19)], vec=trdat1[i,])
  res2<-do.call(rbind,res2)
  trlist1[[i]]<-res2
}

#Calculate modularity
modat1<-matrix(nrow=nrow(trdat1),ncol=50)
for(i in 1:length(trlist1)){
  for(j in 1:50){
    modat1[i,j]<-mod(trlist1[[i]][j,],paramdat1[i,3],paramdat1[i,1])
  }
}


uniques<-unique(trdat[,1])
uniques1<-unique(trdat1[,1])

newdat<-cbind(paramdat,trdat[,1],modat)
newdat<-newdat[which(newdat[,4]%in%uniques1),]
newdat<-newdat[,-4]

newdat1<-cbind(paramdat1,modat1)

colnames(newdat)<-c("t","Tvar","omega",1:50)
colnames(newdat1)<-c("t","Tvar","omega",1:50)

newdat<-as.data.frame(newdat)
newdat$t<-as.factor(newdat$t)
newdat$Tvar<-as.factor(newdat$Tvar)
newdat$omega<-as.factor(newdat$omega)
newdat1<-as.data.frame(newdat1)
newdat1$t<-as.factor(newdat1$t)
newdat1$Tvar<-as.factor(newdat1$Tvar)
newdat1$omega<-as.factor(newdat1$omega)

newdat<-melt(newdat,id=c("t","Tvar","omega"))
newdat$variable<-as.numeric(newdat$variable)
newdat1<-melt(newdat1, id=c("t","Tvar","omega"))
newdat1$variable<-as.numeric(newdat1$variable)



newdat2<-cbind(newdat,newdat1[,5])
newdat2<-na.omit(newdat2)
colnames(newdat2)<-c("t","Tvar","omega","time","mod","nullmod")
newdat2<-melt(newdat2,id=c("t","Tvar","omega","time"))

newdat2.m<-aggregate(newdat2$value,list(newdat2$t,newdat2$Tvar,newdat2$omega,newdat2$time,newdat2$variable),mean,na.rm=TRUE)
newdat2.s<-aggregate(newdat2$value,list(newdat2$t,newdat2$Tvar,newdat2$omega,newdat2$time,newdat2$variable),sd,na.rm=TRUE)
newdat3<-cbind(newdat2.m,newdat.m$x+newdat2.s$x,newdat.m$x-newdat2.s$x)
colnames(newdat3)<-c("t","Tvar","omega","time","label","mean","upper","lower")
dumb3<-function(i){return(index[i,2])}
newdat3$Tvar<-sapply(newdat3$Tvar,dumb3)

newdat4<-newdat3[which(newdat3$t==200),]

nplot<-ggplot(newdat4,aes(x=time,y=mean,group=label,colour=label))+geom_line(size=1.5)+
  geom_ribbon(aes(ymin=lower, ymax=upper, x=time, fill=label),alpha=0.2) +facet_wrap(~Tvar+omega)

nplot

#####################

newdat<-cbind(paramdat,trdat[,1],mnnds)
newdat<-newdat[which(newdat[,4]%in%uniques),]
newdat<-newdat[,-4]

newdat1<-cbind(paramdat1,mnnds1)


colnames(newdat)<-c("t","Tvar","omega",1:50)
colnames(newdat1)<-c("t","Tvar","omega",1:50)

newdat<-as.data.frame(newdat)
newdat$t<-as.factor(newdat$t)
newdat$Tvar<-as.factor(newdat$Tvar)
newdat$omega<-as.factor(newdat$omega)
newdat1<-as.data.frame(newdat1)
newdat1$t<-as.factor(newdat1$t)
newdat1$Tvar<-as.factor(newdat1$Tvar)
newdat1$omega<-as.factor(newdat1$omega)

newdat<-melt(newdat,id=c("t","Tvar","omega"))
newdat$variable<-as.numeric(newdat$variable)
newdat1<-melt(newdat1, id=c("t","Tvar","omega"))
newdat1$variable<-as.numeric(newdat1$variable)


newdat2<-cbind(newdat,newdat1[,5])
newdat2<-na.omit(newdat2)
colnames(newdat2)<-c("t","Tvar","omega","time","mod","nullmod")
newdat2<-melt(newdat2,id=c("t","Tvar","omega","time"))

newdat2.m<-aggregate(newdat2$value,list(newdat2$t,newdat2$Tvar,newdat2$omega,newdat2$time,newdat2$variable),mean,na.rm=TRUE)
newdat2.s<-aggregate(newdat2$value,list(newdat2$t,newdat2$Tvar,newdat2$omega,newdat2$time,newdat2$variable),sd,na.rm=TRUE)
newdat3<-cbind(newdat2.m,newdat.m$x+newdat2.s$x,newdat.m$x-newdat2.s$x)
colnames(newdat3)<-c("t","Tvar","omega","time","label","mean","upper","lower")
dumb3<-function(i){return(index[i,2])}
newdat3$Tvar<-sapply(newdat3$Tvar,dumb3)

newdat4<-newdat3[which(newdat3$t==5),]

nplot<-ggplot(newdat4,aes(x=time,y=mean,group=label,colour=label))+geom_line(size=1.5)+
  geom_ribbon(aes(ymin=lower, ymax=upper, x=time, fill=label),alpha=0.2) +facet_wrap(~Tvar+omega)

nplot




