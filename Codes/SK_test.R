require("tensor")
require(pracma)
library("igraph")
library("ggplot2")
library("reshape2")
library(truncnorm)



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


alpha.tri<-function(a,b,t){
  if(is.na(a)==TRUE | is.na(b)==TRUE){return(0)}
  if(abs(a-b)>=t){return(alpha<-0)}
  slope<-1/t
  if(((a-b)>-t) & ((a-b)<=0)){alpha<-(1+(slope*(a-b)))}
  if(((a-b)>0) & ((a-b)<t)){alpha<-(1-(slope*(a-b)))}
  return(alpha)
}


# Model parameters
nsp<-20   #No. of species
n<-5   # No. of loci
geno<-seq(-5,5,l=2*n+1) #Vector of genotypes
nt<-length(geno)   # No. of genotypes
omega<-2
t<-0.5
a0<-1
a1<-0.5

    
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

R <- array(dim=rep(1+2*n, 3))
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
            if(x==z){A[x,y,z,xy]<-1
              #a0*alpha.gt(geno[xy],geno[y],omega,t)
            } else{A[x,y,z,xy]<-a1*alpha.gt(geno[xy],geno[y],omega,t)}
        }
      }
    }


#Simluation

#Initial densities of genotypes across species.

N0<-t(replicate(nsp,sample(c(1,runif((nt-1),0.05,0.15)),nt,replace=FALSE)))
N0<-N0/rowSums(N0)

N.ini<-rep(1000,nsp)
K<-2*N.ini

r1<-0.25-(pnorm(-5.5,geno,1)+1-pnorm(5.5,geno,1))

r<-outer(rep(0.1,nsp),r1)

freqs<-matrix(nrow=nsp)
Npop<-matrix(nrow=nsp)


#Timepoints to record

Ngen<-N0

Np<-Ngen*N.ini

#plot of initial trait distributions
plot(Np[1,]~geno,type="l",ylim=c(0,1000))
for(i in 2:nsp){
  lines(Np[i,]~geno,col=i)
}

for(m in 1:10000){
  
  #Find the adjust the values for extinct phenotypes
  ext.sp<-which(rowSums(Np)!=0)
  exti.sp<-which(rowSums(Np)==0)
  Ngen[ext.sp,]<-Np[ext.sp,]/rowSums(Np)[ext.sp]
  Ngen[exti.sp,]<-0
  Ngen[Ngen<(10^(-8))] <- 0
  
  if(m%in%samples){
    Ngen1<-Ngen
    Ngen1[which(rowSums(Ngen1)==0),]<-NA
    freqs<-cbind(freqs,Ngen1)
    Npop<-cbind(Npop,Np)
  }
  
  #Reproduction process
  Ngen <- sapply(seq(dim(R)[3]),
                 function(i) rowSums(tensor::tensor(Ngen,R,2,1)[,,i]*Ngen))
  Np<-Ngen*rowSums(Np)
  
  #Selection process
  Np <- Np*(1+r*(1-(tensor::tensor(A, Np, c(3,4), c(1,2))/K)))
  Np[Np<1]<-0


}

freqs<-freqs[,-1]
#freqs<-cbind(Ng0,freqs)
Npop<-Npop[,-1]
#Npop<-cbind(Np0,Npop)

#Plot trajectory of trait means
tps<-rep(1:length(samples),each=nt)

means<-matrix(nrow=nsp,ncol=length(samples))
pops<-matrix(nrow=nsp,ncol=length(samples))
for(i in 1:length(samples)){
  samp<-freqs[,which(tps==i)]
  samp2<-Npop[,which(tps==i)]
  means[,i]<-apply(samp,1,function(j) sum(geno*j)) 
  pops[,i]<-rowSums(samp2)
}

plot(means[1,]~samples,type="l",ylim=c(-5,5),ylab="Trait means",xlab="Time",cex.lab=1.5,mgp=c(1.8,0.5,0))
for(i in 2:nrow(means)){
  lines(means[i,]~samples)
}

pop1<-apply(pops,2,function(j) length(which(j!=0)))


plot(pops[1,]~samples,type="l",ylim=c(0,3500))
for(i in 2:nsp){
  lines(pops[i,]~samples)
}

plot(pops[,length(samples)]~means[,length(samples)],ylab="population size",xlab="Mean trait",ylim=c(0,max(K)))
points(K~means[,length(samples)],pch=3)


#plot the trait trajectory with the line widths indicating the population sizes
par(mgp=c(axis.title.position, axis.label.position, axis.line.position))


tiff('SK_sample2.tiff', units="in", width=10, height=7, res=300)

plot(means[1,1:2]~samples[1:2],type="l",lwd=log(pops[1,2]+1)/2,xlim=c(0,max(samples)),ylim=c(-5,5),ylab="Trait mean",xlab="Time",
     mgp=c(1.5,0.3,0),cex.lab=1.5,font.lab=1.5)

for(i in 2:nsp){
  lines(means[i,1:2]~samples[1:2],lwd=log(pops[i,2]+1)/2)
}
for(j in 2:(length(samples)-1)){
  for(k in 1:nsp){
    lines(means[k,j:(j+1)]~samples[j:(j+1)],lwd=log(pops[k,(j+1)]+1)/2)
  }
}
dev.off()


#############################
#Beta analysis
#function to evaluate beta
beta.num<-function(a,b,omega,t){
  mean.a<-mean(a)
  int<-0
  for(i in 1:10000){
    a1<-sample(a,1)
    b1<-sample(b,1)
    int<-int+((a1-mean.a)*alpha.gt(a1,b1,omega,t))
  }
  int<-int/10000
  return(int)
}


#Choose converging pairs
means1<-means[order(means[,1]),]
diffs<-apply(means1,2,diff)
lms<-c()
for(i in 1:nrow(diffs)){
  lms<-c(lms,lm(diffs[i,]~samples)$coefficients[2])
}

pair<-c(which.min(lms),which.min(lms)+1)



beta.times1<-matrix(nrow=nrow(freqs),ncol=length(samples))

beta.times2<-matrix(nrow=nrow(freqs),ncol=length(samples))

freqs2<-freqs[order(means[,1]),]


tps<-rep(1:length(samples),each=nt)

for(i in 1:length(samples)){
  
  a1<-sample(geno,500,replace=TRUE,prob=freqs2[pair[1],which(tps==i)])
  a2<-sample(geno,500,replace=TRUE,prob=freqs2[pair[2],which(tps==i)])

  for(j in 1:nrow(freqs)){
    beta.times1[j,i]<-beta.num(a1,sample(geno,500,replace=TRUE,prob=freqs2[j,which(tps==i)]),omega,t)
    beta.times2[j,i]<-beta.num(a2,sample(geno,500,replace=TRUE,prob=freqs2[j,which(tps==i)]),omega,t)
  }
}
  

  



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

#######################





