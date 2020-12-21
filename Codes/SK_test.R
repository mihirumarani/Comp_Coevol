require("tensor")
require(pracma)
library("igraph")
library("ggplot2")
library("gganimate")
theme_set(theme_bw())
library("reshape2")
library("transformr")
library("png")
library("grid")
library(gifski)


#Comeptitive kernel

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


# Model parameters
nsp<-10    #No. of species
n<-5      # No. of loci
geno<-seq(-1,1,l=2*n+1) #Vector of genotypes
nt<-length(geno)   # No. of genotypes
omega<-0.5
t<-0.5
a0<-1
a1<-0.2

#Parameters to vary



#Initial growth rates
r<-abs(rnorm(nsp,0,0.1))


#Timepoints to record
samples<-c(1,5*(1:10),100,200,300,500,ceiling(seq(501,1000,length.out=7)))


  #Initial densities of genotypes across species.
  N<-matrix(runif(nsp*nt),nrow=nsp,ncol=nt)
  Ng0<-N/rowSums(N)
  Np0<-runif(nsp,500,1500)*N/rowSums(N)
  
    
#SK model of inheritance: Pre-calculate the probabilities of offspring having a genotype u if parents' genotypes are v and w.
#The output is a 3-d array 
haplR <- array(0, dim=rep(n+1,3))
for(i in 0:n) for(j in 0:i) for(k in 0:min(n,(i+j))){ 
  haplR[1+i,1+j,1+k] <- sum(dhyper(max(0, i+j-n):min(i, j), i, n-i, j) *
                              dbinom(k-(max(0, i+j-n):min(i,j)), i+j -
                                       2*(max(0, i+j-n):min(i, j)), prob=0.5))
}
for (k in 0:n) {
  haplR[,,1+k] <- haplR[,,1+k] + t(haplR[,,1+k])
  diag(haplR[,,1+k]) <- diag(haplR[,,1+k])/2
}
indexsum.haplR <- matrix(0, 2*n+1, 2*n+1)
for(k in 0:n){
  for(i in 0:n) indexsum.haplR[1+i,1+k] <- haplR[1+i,1,1+k]
  for(j in 0:n) indexsum.haplR[1+j+n,1+k] <- haplR[1+n,1+j,1+k]
}

R <- array(dim=rep(nt, 3))
for (i in 0:(2*n)) for (j in 0:(2*n)) for (q in 0:(2*n)) {
  R[1+i,1+j,1+q] <- sum(indexsum.haplR[1+i,1+(0:q)] *
                          indexsum.haplR[1+j,1+q-(0:q)])
}

#Pre-calculate coefficients of competition between pairs of genotypes
A<-array(dim=c(nsp,nt,nsp,nt))
    
    for(x in 1:nsp){
      for(y in 1:nt){
        for(z in 1:nsp){
          for(xy in 1:nt)
            if(x==z){A[x,y,z,xy]<-a0*alpha.gt(geno[xy],geno[y],omega,t)
            } else{A[x,y,z,xy]<-a1*alpha.gt(geno[xy],geno[y],omega,t)}
        }
      }
    }


#Simluation
Ngen<-Ng0
Np<-Np0

freqs<-matrix(nrow=nsp)
Npop<-matrix(nrow=nsp)


for(m in 1:1000){
  
  if(m%in%samples){
    freqs<-cbind(freqs,Ngen)
    Npop<-cbind(Npop,Np)
  }
  
  #Reproduction process
  #Np <- sapply(seq(dim(R)[3]),
   #            function(i) rowSums(tensor::tensor(Ngen,R,2,1)[,,i]*Np))
  ext.sp<-which(rowSums(Np)!=0)
  exti.sp<-which(rowSums(Np)==0)
  Ngen[ext.sp,]<-Np[ext.sp,]/rowSums(Np)[ext.sp]
  Ngen[exti.sp,]<-0
  if(sum(is.nan(rowSums(Ngen)))>0) break
  
  Ngen[Ngen<(10^(-8))] <- 0
  #Selection process
  Np <- Np*(1+r*(1-(tensor::tensor(A, Ngen, c(3,4), c(1,2)))))
  
  Np[Np<1]<-0
}
freqs<-freqs[,-1]
Npop<-Npop[,-1]



########################3
#Movie time! How the trait distributions change with time...

tps<-rep(1:length(samples),each=nt)
freqs1<-matrix(ncol=nt)
for(i in 1:length(samples)){
  freqs1<-rbind(freqs1,freqs[,which(tps==i)])
}
freqs1<-freqs1[-1,]
colnames(freqs1)<-geno
freqs2<-melt(freqs1)
freqs2<-cbind(rep(rep(samples,each=nsp),nt),rep(rep(1:nsp,length(samples)),nt),freqs2[,-1])
colnames(freqs2)<-c("time","species","genotype","frequency")

freqs2<-as.data.frame(freqs2)
freqs2$species<-as.factor(freqs2$species)
freqs2$time<-as.integer(freqs2$time)


p<-ggplot(subset(freqs2,time==1000), aes(x=genotype,y=frequency)) + geom_line(size=1.5) +  ylim(0,0.5)+
 facet_wrap(vars(species))
p
  transition_time(time)
  #facet_wrap(vars(species),nrow=2,ncol=5)+  ylim(0,0.5)
animate(p)

plot(freqs[1,which(tps==6)]~geno,type="l",ylim=c(0,0.5))  
for(i in 2:nsp){lines(freqs[i,which(tps==6)]~geno,col=i)}



#####################################################
#Data analysis
#####################################################
#Nearest neighbot distance
mnnd<-function(a){
  a<-na.omit(a)
  if(length(a)>0){
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

# Model parameters
nsp<-10    #No. of species
as<-rbind(c(0,0),c(1,0),c(1,0.5),c(1,1))
reps.table<-cbind(c(10,20,50),c(10,5,1))
omega<-0.5
t<-1
reps<-10
omega<-0.5
t<-1

#Parameters to vary

ns<-c(10,50,100)
as<-rbind(c(0,0),c(1,0),c(1,0.5),c(1,1))


#Timepoints to record
samples<-c(1,100,200,300,500,ceiling(seq(501,1000,length.out=7)))

#Load data
tot<-readRDS("SKtest.Rds")
param<-matrix(ncol=3)
for(i in 1:length(tot)){
  temp.p<-tot[[i]][[1]]
  param<-rbind(param,matrix(unlist(temp.p),ncol=3,byrow=TRUE))
}
param<-as.matrix(param[-1,])

trdat<-list()
popdat<-list()
for(i in 1:length(tot)){
  temp.1<-tot[[i]][[2]]
  temp.2<-tot[[i]][[3]]
  for(j in 1:length(temp.1)){
    trdat<-c(trdat,list(temp.1[[j]]))
    popdat<-c(popdat,list(temp.2[[j]]))
  }
}

##########
#plot progressions of trait distributions under no selection
traits<-list()
popl<-list()
for(i in 1:length(ns)){
  n<-ns[i]
  geno<-seq(-1,1,l=2*n+1)
  nt<-length(geno)
  ind<-which(param[,1]==ns[i] & param[,2]==1 & param[,3]==0)
  tr<-matrix(ncol=nt*length(samples))
  for(j in 1:length(ind)){
    tr<-rbind(tr,trdat[[ind[j]]])
  }
  tr<-tr[-1,]
  traits<-c(traits,list(tr))
}

#Plot sample of changes in trait distribution for different number of loci
for(i in 1:3){
  n<-ns[i]
  geno<-seq(-1,1,l=2*n+1)
  nt<-length(geno)
  tps<-rep(1:length(samples),each=nt)
  tr<-traits[[i]]
  temp<-matrix(nrow=nrow(tr),ncol=length(samples))
  for(j in 1:length(samples)){
    
  }
}

#calculate the deviation from normality vs. time for all the parameter combinations, use the Shapiro test
shapstat<-matrix(ncol=1+length(samples))
for(i in 1:length(traits)){
  n<-ns[i]
  geno<-seq(-1,1,l=2*n+1)
  nt<-length(geno)
  tps<-rep(1:length(samples),each=nt)
  tr<-traits[[i]]
  temp<-matrix(nrow=nsp*reps,ncol=length(samples))
  for(j in 1:length(samples)){
    a<-tr[,which(tps==j)]
    temp[,j]<-apply(a,1,function(x) if(sum(x!=0)) shapiro.test(sample(geno,100,replace=TRUE,prob=x))$statistic else 0)
  }
  temp<-cbind(rep(n,nsp*reps),temp)
shapstat<-rbind(shapstat,temp)
}

shapstat<-shapstat[-1,]
colnames(shapstat)<-c("n",samples)
shapstat<-as.data.frame(shapstat)
shapstat$n<-as.factor(shapstat$n)
#shapstat<-shapstat[-which(rowSums(shapstat)==0),]
shapstat1<-melt(shapstat,id=c("n"))
plot(shapstat1$value~shapstat1$Var2)
shap<-shapstat1
colnames(shap)<-c("loci","time","normality")

shap.mean<-aggregate(shap$normality,by=list(shap$loci,shap$time),mean)
shap.sd<-aggregate(shap$normality,by=list(shap$loci,shap$time),sd)
shapdat<-cbind(shap.mean,shap.mean$x+shap.sd$x,shap.mean$x-shap.sd$x)
colnames(shapdat)<-c("loci","time","normality","upper","lower")
shapdat$loci<-factor(shapdat$loci,levels=c("10","50","100"))

shaplot<-ggplot(shapdat,aes(x=time,y=normality,group=loci,colour=loci))+geom_line(size=1.5)+
  geom_ribbon(aes(ymin=lower,ymax=upper,x=time,fill=loci),alpha=0.2)+ylab("Normality")
shaplot
