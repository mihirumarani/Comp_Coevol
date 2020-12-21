library("igraph")
library("ggplot2")
library("reshape2")


#Necessary functions
alpha1<-function(a,b,omega,t){
  if(is.na(a)==TRUE | is.na(b)==TRUE){return(0)}
  alpha<-exp(-((a-b)^2/(omega)^2))
  if(abs(a-b)>=t){alpha<-0}
  return(alpha)           
}

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
  if(length(a)>0){
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



#Parameters
nsp<-20
h<-rep(0.5,nsp)
omegas<-c(0.1,0.5)                             #Width of competition kernel
ts<-c(1,2.5,5,100)                           #Threshold to competition kernel
a0<-1
a1s<-c(0.2,0.5,0.9)                            #Strength of interspecific competition relative to intraspecific comp.
omegas<-c(0.1,0.5)
tvar<-runif(nsp,0.5,2.5) 

reps<-10

samples<-ceiling(c(seq(1,2000,length.out=10),seq(2100,29999,length.out=30)))

#Raw data:


tot<-readRDS("QG20.Rds")
param<-matrix(ncol=3)
trdat<-matrix(ncol=length(samples)*nsp)
popdat<-matrix(ncol=length(samples)*nsp)
for(i in 1:length(tot)){
  temp.p<-tot[[i]][[1]]
  temp.pop<-tot[[i]][[2]]
  temp.tr<-tot[[i]][[3]]
  param<-rbind(param,matrix(unlist(temp.p),ncol=3,byrow=TRUE))
  popdat<-rbind(popdat,matrix(unlist(temp.pop),ncol=length(samples)*nsp,byrow=TRUE))
  trdat<-rbind(trdat,matrix(unlist(temp.tr),ncol=length(samples)*nsp,byrow=TRUE))
}

param<-param[-1,]    #1st column- Threshold, 2nd column- omega, 3rd column-a1
popdat<-popdat[-1,]
trdat<-trdat[-1,]



#Calculate modularity
tsteps<-rep(1:length(samples),each=nsp)
modat<-matrix(nrow=nrow(trdat),ncol=length(samples))
for(i in 1:nrow(trdat)){
  for(j in 1:length(samples)){
    modat[i,j]<-mods(trdat[i,which(tsteps==j)],param[i,3],param[i,1])
  }
}





#Calculate MNND
tsteps<-rep(1:length(samples),each=nsp)
mnnds<-matrix(nrow=nrow(trdat),ncol=length(samples))
for(i in 1:nrow(trdat)){
  for(j in 1:length(samples)){
    samp<-trdat[i,which(tsteps==j)]
    samp<-samp[-which(samp==0)]
    mnnds[i,j]<-mnnd(samp)
  }
}





#Modularity index vs. time and how it differs across param space

modtimes<-modat

#plot heatmaps
modtimes<-as.data.frame(modtimes)
modmeans1<-as.data.frame(cbind(param,modmeans))
colnames(modmeans1)<-c("t","Tvar","omega","means")
modmeans1$t<-as.factor(modmeans1$t)
modmeans1$Tvar<-as.factor(modmeans1$Tvar)
modmeans1$omega<-as.factor(modmeans1$omega)
ggplot(modmeans1,aes(x=t,y=Tvar,fill=means))+geom_tile()+scale_fill_gradient(low="white", high="blue")

modvars<-apply(modtimes,1,var)
modvars1<-as.data.frame(cbind(paramdat,modvars))
colnames(modvars1)<-c("t","Tvar","omega","vars")
modvars1$t<-as.factor(modvars1$t)
modvars1$Tvar<-as.factor(modvars1$Tvar)
modvars1$omega<-as.factor(modvars1$omega)
ggplot(modvars1,aes(x=t,y=Tvar,fill=vars))+geom_tile()+scale_fill_gradient(low="white", high="blue")




modtimes<-cbind(param,modtimes)
modtimes1<-modtimes
modtimes1<-as.data.frame(modtimes1)
colnames(modtimes1)<-c("t","Tvar","Omega",1:length(samples))



modtimes1<-as.data.frame(modtimes1)
modtimes1$t<-as.factor(modtimes1$t)
modtimes1$Tvar<-as.factor(modtimes1$Tvar)
modtimes1$Omega<-as.factor(modtimes1$Omega)



#Plot no. 1: Modularity vs. time 
plot1dat<-melt(modtimes1,id=c("t","Tvar","Omega"))
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



plot1dat1$t<-factor(plot1dat1$t,levels=c("200","7","5","3"))
plot1dat1$Intra_var<-factor(plot1dat1$Intra_var,levels=c("1","2"))
plot1dat1$omega<-factor(plot1dat1$omega,levels=c("0.1", "0.5", "1"))


var.labs<-c("Low ITV","High ITV")
names(var.labs)<-c("1","2")

t.labs<-c("Threshold=3","Threshold=5","Threshold=7","No Threshold")
names(t.labs)<-c("3","5","7","200")


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

tiff('plot.mod2.tiff', units="in", width=10, height=5, res=300)
plot1
dev.off()

#, compression = 'lzw'
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
#Calculate mnnds from trait data
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

plot2datm<-aggregate(plot2dat$value,list(plot2dat$t,plot2dat$Intra,plot2dat$Omega,plot2dat$variable),mean,na.rm=TRUE)
plot2dats<-aggregate(plot2dat$value,list(plot2dat$t,plot2dat$Intra,plot2dat$Omega,plot2dat$variable),sd,na.rm=TRUE)
plot2dat1<-cbind(plot2datm,plot2datm$x+plot2dats$x,plot2datm$x-plot2dats$x)
colnames(plot2dat1)<-c("t","omega","Intra","Time","mean","upper","lower")


plot2dat1$Intra<-factor(plot2dat1$Intra,levels=c("0.2","0.5","0.9"))
plot2dat1$t<-factor(plot2dat1$t,levels=c("100","1.5","1","0.5"))
plot2dat1$omega<-factor(plot2dat1$omega,levels=c("0.1", "0.5"))

t.labs<-c("Thr=0.5","Thr=1","Thr=1.5","No Thr")
names(t.labs)<-c("0.5","1","1.5","100")

var.labs<-c("Weak Competition","Intermediate Competition","Strong Competition")
names(var.labs)<-c("0.2","0.5","0.9")


plot2<-ggplot(plot2dat1,aes(x=Time,y=mean, group=omega,colour=omega)) +geom_line(size=1.5) +
  geom_ribbon(aes(ymin=lower, ymax=upper, x=Time, fill=omega),alpha=0.2)+ylab("Mean nearest neighbor distance")+
  facet_grid(Intra~t,labeller=labeller(Intra_var=var.labs,t=t.labs),switch="y")+
  theme(
    strip.text.x = element_text(
      size = 13, face = "bold"
    ),
    strip.text.y = element_text(
      size = 12, face = "bold"
    ),
    legend.key.size = unit(1.5, "cm"),
    legend.text=element_text(size=15),
    axis.title=element_text(size=20,face="bold")
  )
plot2


tiff('plot.mnnd2.tiff', units="in", width=10, height=5, res=300)
plot2
dev.off()


#, compression = 'lzw'



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


#Build a convergence vs. divergence plot
conv<-matrix(nrow=length(trlist),ncol=49)
for(i in 1:length(trlist)){
  dat<-trlist[[i]]
  dat<-dat[,order(dat[1,])]
  diffs<-t(apply(dat,1,diff))
  diffs2<-(apply(diffs,2,diff))
  conv[i,]<-apply(diffs2,1,function(x) sum(x<0))
}
div<-matrix(nrow=nrow(conv),ncol=ncol(conv))
div<-20-conv
condiv1<-cbind(paramdat,conv[,1],div[,1])

condiv<-cbind(paramdat,conv)


colnames(condiv)<-c("t","Tvar","Omega",1:49)
condiv<-as.data.frame(condiv)
condiv$t<-as.factor(condiv$t)
condiv$Tvar<-as.factor(condiv$Tvar)
condiv$Omega<-as.factor(condiv$Omega)


#Plot no. 1: Modularity vs. time 
condiv1<-condiv
plot1dat<-melt(condiv1,id=c("t","Tvar","Omega"))
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
  geom_ribbon(aes(ymin=lower, ymax=upper, x=time, fill=omega),alpha=0.2)+ylab("convergences") + facet_grid(Tvar~t,labeller=label_both)
plot1

######################################################
#Plot convergence vs divergece matrix
#Function to calculate pairwise trait difference
pairdiff<-function(a){
  mat<-matrix(nrow=length(a),ncol=length(a))
  for(i in 1:length(a)){
    mat[,i]<-a[i]-a
  }
  return(mat)
}


conv<-matrix(nrow=50,ncol=60)
div<-matrix(nrow=50,ncol=60)
for(i in 1:length(trlist)){
  for(j in 1:50){
    traits<-trlist[[i]][j,]
    premat<-pairdiff(traits)
    newtr<-traits+trchangelist[[i]][j,]
    postmat<-pairdiff(newtr)
    conv[j,i]<-(length(which((abs(postmat)-abs(premat))<0))-20)/2
    div[j,i]<-(length(which((abs(postmat)-abs(premat))>0))-20)/2
  }
}
conv<-t(conv)
div<-t(div)
conv<-cbind(paramdat,conv)
div<-cbind(div,paramdat)


dumb2<-function(a){
  i<-a[1]
  j<-a[2]
  if(i+j==2){return(1)}
  if(i+j==4){return(3)}
  if(i+j==6){return(2)}
}
newvar<-apply(paramdat[,3:4],1,dumb2)
paramdat<-cbind(paramdat[,2],newvar,paramdat[,5])
conv<-cbind(paramdat,conv)
conv<-conv[,-4]

colnames(conv)<-c("t","tvar","omega",1:50)
conv<-as.data.frame(conv)
conv$t<-as.factor(conv$t)
conv$tvar<-as.factor(conv$tvar)
conv$omega<-as.factor(conv$omega)
conv1<-melt(conv,id=c("t","tvar","omega"))
conv1$variable<-as.numeric(conv1$variable)
dumb<-function(i){return(modsamples[i])}
conv1$variable<-sapply(conv1$variable,dumb)
dumb3<-function(i){return(index[i,2])}
conv1$tvar<-sapply(conv1$tvar,dumb3)
conv1$tvar<-factor(conv1$tvar,levels=c("low","high","wide"),labels=c("U[0.5,1,5]","U[2.5,3.5]","U[0.5,3.5]"))

conv1$t<-factor(conv1$t,labels=c("5","7","10","No threshold"))
colnames(conv1)<-c("t","tvar","omega","time","convergences")
converge<-ggplot(conv1,aes(x=time,y=convergences,group=omega,colour=omega))+geom_line(size=1.5)+facet_grid(tvar~t)
converge


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
richdat2$Tvar<-as.factor(richdat2$Tvar)
richdat.m<-aggregate(richdat2$value,list(richdat2$t,richdat2$Tvar,richdat2$Omega,richdat2$variable),mean,na.rm=TRUE)
richdat.s<-aggregate(richdat2$value,list(richdat2$t,richdat2$Tvar,richdat2$Omega,richdat2$variable),sd,na.rm=TRUE)
richdata<-cbind(richdat.m,richdat.m$x+richdat.s$x,richdat.m$x-richdat.s$x)
colnames(richdata)<-c("t","Tvar","omega","time","mean","upper","lower")
richdata$lower[which(richdata$lower<=0)]<-0

richdata$Tvar<-factor(richdata$Tvar,levels=c("Low","Wide","High"),labels=c("1","2","3"))

richdata$t<-factor(richdata$t,labels=c("5","7","10","No threshold"))

richplot<-ggplot(richdata,aes(x=time,y=mean,color=omega))+geom_line(size=1.5)+ylab("Species Richness")+
  facet_grid(Tvar~t,labeller=label_both)
richplot

#geom_ribbon(aes(ymin=lower, ymax=upper, x=time, fill=omega),alpha=0.2)


#Plot connectance vs. time


conndat<-matrix(nrow=nrow(trdat),ncol=50)
times<-1+(0:49)*20
for(i in 1:nrow(trdat)){
  for(j in 1:50){
    dat<-na.omit(trdat[i,(times[j]:(times[j]+19))])
    if(length(dat)>0){
      net<-network(dat,paramdat[i,3],paramdat[i,1])
      n<-nrow(net)
      conn<-length(which(net!=0))/(n*(n-1))
    }
    else{conn<-NA}
    conndat[i,j]<-conn
  }
}


#Plot heatmaps

connmean<-apply(conndat,1,mean)
connmean1<-as.data.frame(cbind(paramdat,connmean))
colnames(connmean1)<-c("t","Tvar","omega","mean")
connmean1$t<-as.factor(connmean1$t)
connmean1$Tvar<-as.factor(connmean1$Tvar)
connmean1$omega<-as.factor(connmean1$omega)

ggplot(connmean1,aes(t,omega,fill=mean))+geom_tile()+scale_fill_gradient(low="white", high="blue") 



connvar<-apply(conndat,1,var)
connvar1<-as.data.frame(cbind(paramdat,connvar))
colnames(connvar1)<-c("t","Tvar","omega","variance")
connvar1$t<-as.factor(connvar1$t)
connvar1$Tvar<-as.factor(connvar1$Tvar)
connvar1$omega<-as.factor(connvar1$omega)

ggplot(connvar1,aes(t,Tvar,fill=variance))+geom_tile()+scale_fill_gradient(low="white", high="blue")






conndat<-cbind(paramdat,conndat)
conndat<-as.data.frame(conndat)
colnames(conndat)<-c("t","Tvar","Omega",1:50)


conndat1<-conndat[which(conndat$t%in%c(200,5)),]
conndat1<-conndat1[which(conndat1$Tvar%in%c(1,2)),]
conndat1<-conndat[which(conndat1$Omega%in%c(0.05,0.25,1)),]
conndat1<-as.data.frame(conndat1)



conndat1<-as.data.frame(conndat)
conndat1$t<-as.factor(conndat1$t)
conndat1$Tvar<-as.factor(conndat1$Tvar)
conndat1$Omega<-as.factor(conndat1$Omega)
conndat2<-melt(conndat1,id=c("t","Tvar", "Omega"))
conndat2$variable<-as.numeric(conndat2$variable)
dumb<-function(i){return(modsamples[i])}
conndat2$variable<-sapply(conndat2$variable,dumb)
#dumb3<-function(i){return(index[i,2])}
#conndat2$Tvar<-sapply(conndat2$Tvar,dumb3)
conndat.m<-aggregate(conndat2$value,list(conndat2$t,conndat2$Tvar,conndat2$Omega,conndat2$variable),mean,na.rm=TRUE)
conndat.s<-aggregate(conndat2$value,list(conndat2$t,conndat2$Tvar,conndat2$Omega,conndat2$variable),sd,na.rm=TRUE)
conndata<-cbind(conndat.m,conndat.m$x+conndat.s$x,conndat.m$x-conndat.s$x)
colnames(conndata)<-c("t","Intra_var","omega","Time","mean","upper","lower")
conndata$lower[which(conndata$lower<=0)]<-0


conndata$Intra_var<-factor(conndata$Intra_var,levels=c( "1", "3", "2"))

conndata$t<-factor(conndata$t,levels=c("200","10","7","5","3"))
conndata$omega<-factor(conndata$omega,levels=c("0.05", "0.25", "1"),labels=c("Weak Competition","Intermediate Competition","Strong Competition"))

o.labs<-c("Weak Competition","Intermediate Competition","Strong Competition")
names(o.labs)<-c("0.05","0.25","1")

var.labs<-c("Low ITV","High ITV","Wide range of ITV")
names(var.labs)<-c("1","3","2")


t.labs<-c("Threshold=3","Threshold=5","Threshold=7","Threshold=10","No Threshold")
names(t.labs)<-c("3","5","7","10","200")


connplot<-ggplot(conndata,aes(x=Time,y=mean,color=omega))+geom_line(size=1.5)+ylab("Connectance")+
  geom_ribbon(aes(ymin=lower, ymax=upper, x=Time, fill=omega),alpha=0.2)+
  facet_grid(Intra_var~t,labeller=labeller(Intra_var=var.labs,t=t.labs),switch="y")+
  theme(
    strip.text.x = element_text(
      size = 12, face = "bold"
    ),
    strip.text.y = element_text(
      size = 10, face = "bold"
    ),
    legend.key.size = unit(1.5, "cm"),
    legend.text=element_text(size=12,face="bold"),
    axis.title=element_text(size=17,face="bold")
  )

connplot

tiff('plot.conn1.tiff', units="in", width=10, height=5, res=300)
connplot  
dev.off()







