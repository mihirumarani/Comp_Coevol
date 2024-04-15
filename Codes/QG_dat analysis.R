library("igraph")
library("ggplot2")
library("reshape2")

mnnddat=NULL

setwd="C:/Users/mihir"

files=list.files(pattern="compSK")

for(i in 1:length(files)){
  
  dat=read.csv(files[i])%>%as_tibble()
  
  mnnddat=bind_rows(mnnddat,
                    dat%>%
                      group_by(nloci,reps,omega,t,a1,time)%>%
                      res=mnnd(trmean)%>%
                      ungroup())
                    
  
}



#Necessary functions
alpha1<-function(a,b,omega,t){
  if(is.na(a)==TRUE | is.na(b)==TRUE){return(0)}
  alpha<-exp(-((a-b)^2/(omega)^2))
  if(abs(a-b)>=t){alpha<-0}
  return(alpha)           
}

erfun<-function(x){
  return(2 * pnorm(x * sqrt(2)) - 1)
}


al.ana<-function(t1,t2,v1,v2,omega,t){
  if(any(is.na(c(t1,t2)))){
    return(0)
  } else{
    
    x1<-0.5*omega/sqrt(((omega^2)+(2*(v1^2))+(2*(v2^2))))
    x2<-exp(-((t1-t2)^2)/((omega^2)+(2*(v1^2))+(2*(v2^2))))
    
    a<-((omega^2)+(2*(v1^2))+(2*(v2^2)))/(2*(v1^2)*((omega^2)+(2*(v2^2))))
    
    b<-(((omega^2)*t1)+(2*(v2^2)*t1)+(2*(v1^2)*t2))/((v1^2)*((omega^2)+(2*(v2^2))))
    
    a1<-omega/(v2*sqrt((2*(omega^2))+(4*(v2^2))))
    
    b1<-(((omega^2)*t)+(2*(v2^2)*t)-((omega^2)*t2))/(v2*omega*sqrt((2*(omega^2))+(4*(v2^2))))
    
    b2<-(-((omega^2)*t)-(2*(v2^2)*t)-((omega^2)*t2))/(v2*omega*sqrt((2*(omega^2))+(4*(v2^2))))
    
    
    
    c1<-((2*a*b1)+(a1*b))/(2*sqrt((a^2)+(a*(a1^2))))
    c2<-((2*a*b2)+(a1*b))/(2*sqrt((a^2)+(a*(a1^2))))
    
    int<-x1*x2*(erfun(c1)-erfun(c2))
    
    return(int)
  }
}


#beta function to see the total trait change from pairwise interaction


be.ana<-function(t1,t2,v1,v2,omega,t){
  if(any(is.na(c(t1,t2)))){
    return(0)
  } else{
    
    
    x1<-omega*(v1^2)*(t2-t1)/(((omega^2)+(2*(v1^2))+(2*(v2^2)))^1.5)
    x2<-exp(-((t1-t2)^2)/((omega^2)+(2*(v1^2))+(2*(v2^2))))
    
    x3<-sqrt(2)*(v1^2)*(omega^2)/(sqrt((v1^2)+(v2^2))*((omega^2)+(2*(v1^2))+(2*(v2^2))))
    
    
    
    # x1<-0.5*omega/(v1*sqrt(2*pi*((omega^2)+(2*(v2^2)))))
    #x2<-exp(-(((omega^2)*t1^2)+(2*(v2^2)*(t1^2))+(2*(v1^2)*(t2^2)))/(2*(v1^2)*((omega^2)+(2*(v2^2)))))
    
    a<-((omega^2)+(2*(v1^2))+(2*(v2^2)))/(2*(v1^2)*((omega^2)+(2*(v2^2))))
    
    b<-(((omega^2)*t1)+(2*(v2^2)*t1)+(2*(v1^2)*t2))/((v1^2)*((omega^2)+(2*(v2^2))))
    
    a1<-omega/(v2*sqrt((2*(omega^2))+(4*(v2^2))))
    
    b1<-(((omega^2)*t)+(2*(v2^2)*t)-((omega^2)*t2))/(v2*omega*sqrt((2*(omega^2))+(4*(v2^2))))
    
    b2<-(-((omega^2)*t)-(2*(v2^2)*t)-((omega^2)*t2))/(v2*omega*sqrt((2*(omega^2))+(4*(v2^2))))
    
    x4<-exp((-(((t1^2)/(2*(v1^2)))+((t2^2)/((omega^2)+(2*(v2^2))))))+((b^2)/((4*a)+(4*(a1^2)))))
    #c1<-sqrt(pi/a)
    #c2<-exp((b^2)/(4*a))
    c1<-((2*a*b1)+(a1*b))/(2*sqrt((a^2)+(a*(a1^2))))
    c2<-((2*a*b2)+(a1*b))/(2*sqrt((a^2)+(a*(a1^2))))
    
    #int1<-t1*c1*c2*erfun(c3)
    #int2<-t1*c1*c2*erfun(c4)
    
    int1<-x1*x2*(erfun(c1)-erfun(c2))
    
    #c5<-sqrt(pi)*b/(2*a*sqrt(a))
    #c6<-a1/(a*sqrt(a+(a1^2)))
    c3<-(-(4*a*(b1^2))-(4*a1*b*b1))/((4*a)+(4*(a1^2)))
    c4<-(-(4*a*(b2^2))-(4*a1*b*b2))/((4*a)+(4*(a1^2)))
    
    int2<-x3*x4*(exp(c3)-exp(c4))
    return(int1+int2)
  }
}

network<-function(means,vars,omega,t){
  means<-na.omit(means)
  mat<-matrix(ncol=length(means),nrow=length(means))
  samples<-1:length(means)
  for(i in samples){
    for(j in samples){
      
      mat[i,j]<-al.ana(means[i],means[j],vars[i],vars[j],omega,t)/
        al.ana(means[i],means[i],vars[i],vars[j],omega,t)
      
    }
  }
  lists<-is.na(mat)
  mat[lists]<-0
  diag(mat)<-0
  mat<-mat/max(mat)
  return(mat)
}



mod<-function(a,tvar,omega,t){
  tvar<-tvar[!is.na(a)]
  a<-na.omit(a)
  if(length(a)>0){
    net<-network(a,tvar,omega,t)
    G<-graph_from_adjacency_matrix(net,mode="undirected",weighted=TRUE,diag=FALSE)
    fgreedy<-fastgreedy.community(G,merges=TRUE, modularity=TRUE)
    return(c(max(fgreedy$modularity),length(unique(fgreedy$membership))))}
  else{return(NA)}
}



#Nearest neighbot distance
mnnd<-function(a){
  a<-na.omit(a)
  if(length(a)>1){
    a<-a[order(a,decreasing=TRUE)]
    nnd<-c(length=length(a))
    vector(length=length(a))
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

reps<-20

samples<-ceiling(c(seq(1,1000,length.out=10),seq(2000,29999,length.out=20)))


#Raw data:
tot<-readRDS("QG20_f.Rds")
tvar<-tot[[length(tot)]]
tot<-tot[-length(tot)]
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


#Calculate species richness
tsteps<-rep(1:length(samples),each=nsp)
richness<-matrix(nrow=nrow(popdat),ncol=length(samples))
for(i in 1:length(samples)){
  samp<-popdat[,which(tsteps==i)]
  richness[,i]<-apply(samp,1,function(j) length(which(j>0)))
}


#Calculate modularity
tsteps<-rep(1:length(samples),each=nsp)
modat<-matrix(nrow=nrow(trdat),ncol=length(samples))

for(i in 1:nrow(trdat)){
  for(j in 1:length(samples)){
    modat[i,j]<-mod(trdat[i,which(tsteps==j)],tvar,param[i,2],param[i,1])[1]
  }
}

modat.m<-matrix(nrow=nrow(trdat),ncol=length(samples))

for(i in 1:nrow(trdat)){
  for(j in 1:length(samples)){
    modat.m[i,j]<-mod(trdat[i,which(tsteps==j)],tvar,param[i,2],param[i,1])[2]
  }
}

#Calculate MNND
tsteps<-rep(1:length(samples),each=nsp)
mnnds<-matrix(nrow=nrow(trdat),ncol=length(samples))
for(i in 1:nrow(trdat)){
  for(j in 1:length(samples)){
    samp<-trdat[i,which(tsteps==j)]
    samp<-samp[which(samp!=0)]
    mnnds[i,j]<-mnnd(samp)
  }
}



#Plot no. 1: Modularity vs. time 
modtimes1<-modat
modtimes1<-cbind(param,modtimes1)
modtimes1<-as.data.frame(modtimes1)
colnames(modtimes1)<-c("t","Omega","Tvar",1:length(samples))
modtimes1$t<-as.factor(modtimes1$t)
modtimes1$Tvar<-as.factor(modtimes1$Tvar)
modtimes1$Omega<-as.factor(modtimes1$Omega)

plot1dat<-melt(modtimes1,id=c("t","Tvar","Omega"))
plot1dat$variable<-as.numeric(plot1dat$variable)
dumb<-function(i){return(samples[i])}
plot1dat$variable<-sapply(plot1dat$variable,dumb)


plot1datm<-aggregate(plot1dat$value,list(plot1dat$t,plot1dat$Tvar,plot1dat$Omega,plot1dat$variable),mean,na.rm=TRUE)
plot1dats<-aggregate(plot1dat$value,list(plot1dat$t,plot1dat$Tvar,plot1dat$Omega,plot1dat$variable),sd,na.rm=TRUE)
plot1dat1<-cbind(plot1datm,plot1datm$x+plot1dats$x,plot1datm$x-plot1dats$x)
colnames(plot1dat1)<-c("t","Intra","omega","Time","mean","upper","lower")

plot1dat1<-plot1dat1[which(plot1dat1$t%in%c(0.5,100) & plot1dat1$Intra%in%c(0.2,0.9) & plot1dat1$Time<5001),]

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

plot1<-ggplot(plot1dat1,aes(x=Time,y=mean,colour=omega1)) +geom_line(size=1.5) + ylim(0,1)+
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



tiff('QG20_short_Mod.tiff', units="in", width=10, height=5, res=300)
plot1
dev.off()

#, compression = 'lzw'


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

plot2dat<-na.omit(plot2dat)

plot2datm<-aggregate(plot2dat$value,list(plot2dat$t,plot2dat$Omega,plot2dat$Intra,plot2dat$variable),mean,na.rm=TRUE)
plot2dats<-aggregate(plot2dat$value,list(plot2dat$t,plot2dat$Omega,plot2dat$Intra,plot2dat$variable),sd,na.rm=TRUE)
plot2dat1<-cbind(plot2datm,plot2datm$x+plot2dats$x,plot2datm$x-plot2dats$x)
colnames(plot2dat1)<-c("t","omega","Intra","Time","mean","upper","lower")

plot2dat1<-plot2dat1[which(plot2dat1$t%in%c(0.5,100) & plot2dat1$Intra%in%c(0.2,0.9) & plot2dat1$Time<5001),]


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


tiff('QG20_short_MNND.tiff', units="in", width=10, height=5, res=300)
plot2
dev.off()


#, compression = 'lzw'

#########################
#Plot species richness

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


tiff('QG20_short_Rich.tiff', units="in", width=10, height=5, res=300)
plot3
dev.off()

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







