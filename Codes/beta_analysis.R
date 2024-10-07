#Functions
library(truncnorm)
library(tidyverse)
library(patchwork)


theme_set(theme_bw())
theme_update(panel.grid = element_blank(),
             strip.background = element_rect(fill = 'orange'))

cbpalette <-
  c(
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#C5E772",
    "#4F2CAA"
  )


#Functions for competitive kernels
alpha_gt=function(x,y,omega,t){
  
  al=0
  if(abs(x-y)<t){
    
    al=exp(-((x-y)^2)/(omega^2))
  }
  
  return(al)
}


alpha_tri=function(x,y,slope,t){
  
  al=0
  if(abs(x-y)<t){
    
    al=1-(slope*abs(x-y))
  }
  
  return(max(0,al))
  
}

#function to get the probability of a getting a phenotype under SK model
qgprob=function(n){
  
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
  
  
  return(R)
  
  
}


single_sim1=function(time,r,K1,K2,a1,A,R,Ng0,Npop){
  
  Np0=Ng0*Npop
  Ngen=Ng0
  Np=Np0
  
  nsp=nrow(Ng0)
  nt=ncol(Ng0)
  
  dat=array(dim=c(time+1,nsp,nt))
  
  dat[1,,]=Np
  
  #Start the simulation
  for(m in 2:(time+1)){
  
    #Determine the extinct species
    Np[which(rowSums(Np) < 10 ),]=0
    Ngen[which(rowSums(Np)==0),]=0
    
    if(all(rowSums(Np)==0)){
      break
    }else{
      
      newgen=matrix(0,nsp,nt)
      Np1=newgen
      
      #Reproduction
      for(i in which(rowSums(Ngen)!=0)){
        
        probs=Ngen[i,]%*%t(Ngen[i,])
        
        for(j in 1:nt){
          
          newgen[i,j]=sum(probs*R[,,j])
        }
      }
      
      newp=newgen*rowSums(Np)
      
      #Selection
      
      if(nrow(A)>1){
        
        for(i1 in 1:nrow(newp)){
          
          #Intraspecific density dependence: Logistic
          rdash=r[i1]*(1-(sum(newp[i1,])/K2))
          
          #Interspecific comp: Lotka-Volterra
          
          for(i2 in 1:nt){
            
            comps=a1*sum(A[i2,]%*%t(newp[-i1,]))
            
            Np[i1,i2]= newp[i1,i2]+ (newp[i1,i2]*rdash*(1-(comps/K1)))
            
          }
          
        }
        
        Np[which(Np < 1)]=0
        Ngen=Np/rowSums(Np)
        Ngen[is.nan(Ngen)]=0
      }
      
    }
    
    dat[m,,]=Np
        
  }
  
  return(dat)

}

getsum=function(dat){
  
  pars=dim(dat)
  
  pops=matrix(0,pars[1],pars[2])
  trmeans=pops
  
  for(i in 1:pars[2]){
    for(j in 1:pars[1]){
    
    pops[j,i]=sum(dat[j,i,])
    trmeans[j,i]=sum(res[j,i,]*geno)/sum(res[j,i,])
    }
  }
  
  l=list()
  
  l[[1]]=pops
  l[[2]]=trmeans
  
  return(l)
  
  
}

################################################################################ 
#Plot beta values as functions of competition kernels and trait distributions
################################################################################ 

n=10
geno=seq(-1,1,length.out=2*n+1)

nt=length(geno)

#Parameters
omega=1
slope=0.75
t=1
a1=0.25

#Set the mean trait values for spp no. 2
means=seq(-1.15,1.15,0.005)


#Pre-calculate coefficients of competitoin between pairs of phenotypes

Ag=matrix(0,nt,nt)
At=Ag

for(i1 in 1:nt){
  for(i2 in 1:nt){
    
    Ag[i1,i2]=alpha_gt(geno[i1],geno[i2],omega,t)
    At[i1,i2]=alpha_tri(geno[i1],geno[i2],slope,t)
    
  }
}


#Case 1: Uniformly distributed traits + truncated gaussian comp. kernel

#Assign probabilities for all phenotypes of spp 1 and assume that mean trait value
#for sp 1 is zero.

N1=dunif(geno,-0.2,0.2)
N1=N1/sum(N1)

betas1=vector(length=length(means))

for(i in 1:length(means)){
  
  #Create a trait distribution around the mean trait value for spp. 2
  N2=dunif(geno,means[i]-0.2,means[i]+0.2)
  
  N2=N2/sum(N2)
  
  N1n=vapply(1:length(N1), function(x) sum(Ag[x,]*N2) , 1)
  
  newfreq=N1*(1-N1n)
  
  newfreq= newfreq/sum(newfreq)
  
  betas1[i]=sum(newfreq*geno)
  
}



#Case 2: Uniformly distributed traits + triangular comp. kernel

N1=dunif(geno,-0.2,0.2)
N1=N1/sum(N1)

betas2=vector(length=length(means))

for(i in 1:length(means)){
  
  #Create a trait distribution around the mean trait value for spp. 2
  N2=dunif(geno,means[i]-0.2,means[i]+0.2)
  
  N2=N2/sum(N2)
  
  N1n=vapply(1:length(N1), function(x) sum(At[x,]*N2) , 1)
  
  newfreq=N1*(1-N1n)
  
  newfreq= newfreq/sum(newfreq)
  
  betas2[i]=sum(newfreq*geno)
  
}


#Case 3: Gaussian distributed traits + truncated gaussian comp. kernel

N1=dtruncnorm(geno,-1,1,mean=0,sd=0.2)
N1=N1/sum(N1)

betas3=vector(length=length(means))

for(i in 1:length(means)){
  
  #Create a trait distribution around the mean trait value for spp. 2
  N2=dtruncnorm(geno,-1,1,mean=means[i],sd=0.2)
  
  N2=N2/sum(N2)
  
  N1n=vapply(1:length(N1), function(x) sum(Ag[x,]*N2) , 1)
  
  newfreq=N1*(1-N1n)
  
  newfreq= newfreq/sum(newfreq)
  
  betas3[i]=sum(newfreq*geno)
  
}


#Case 4: Gaussian distributed traits + triangular comp. kernel
N1=dtruncnorm(geno,-1,1,mean=0,sd=0.2)
N1=N1/sum(N1)

betas4=vector(length=length(means))

for(i in 1:length(means)){
  
  #Create a trait distribution around the mean trait value for spp. 2
  N2=dtruncnorm(geno,-1,1,mean=means[i],sd=0.2)
  
  N2=N2/sum(N2)
  
  N1n=vapply(1:length(N1), function(x) sum(At[x,]*N2) , 1)
  
  newfreq=N1*(1-N1n)
  
  newfreq= newfreq/sum(newfreq)
  
  betas4[i]=sum(newfreq*geno)
  
}

res=data.frame(mean=means,beta1=betas1,beta2=betas2,
               beta3=betas4,beta4=betas4)

p1=res%>%ggplot(aes(x=means,y=beta1))+
  geom_line()+
  ggtitle("Uniform distribution + Gaussian Kernel")

p2=res%>%ggplot(aes(x=means,y=beta2))+
  geom_line()+
  ggtitle("Uniform distribution + Triangular Kernel")

p3=res%>%ggplot(aes(x=means,y=beta3))+
  geom_line()+
  ggtitle("Gaussian distribution + Gaussian Kernel")

p4=res%>%ggplot(aes(x=means,y=beta4))+
  geom_line()+
  ggtitle("Gaussian distribution + Triangular Kernel")


(p1|p2)/(p3|p4)




#############################################################################
##Long term simulations for multispecies competition
#############################################################################
#Sample initial conditions and param values

nsp=10
n=5

#Traits range between -1 and 1
geno=seq(-1,1,length.out=(2*n+1))

nt=length(geno)

#Params
#kernel related
omega=1
slope=0.7
t=1.0
a1=0.25

#Demography related
r=abs(runif(nsp,0.2,0.3))

K1=1500000
K2=1500


#pre-calculate coefficients of competition between pairs of phenotypes

Ag=matrix(0,nt,nt)
At=Ag

for(i1 in 1:nt){
  for(i2 in 1:nt){
    Ag[i1,i2]=alpha_gt(geno[i1],geno[i2],omega,t)
    At[i1,i2]=alpha_tri(geno[i1],geno[i2],slope,t)
  }
}

R=qgprob(n)

#Starting population
randm=runif(nsp,-0.6,0.6)

N_uni=matrix(0,nsp,nt)

N_gau=N_uni

for(i in 1:nsp){
  
  N_uni[i,]=dunif(geno,randm[i]-0.4,randm[i]+0.4)
  
  N_gau[i,]=dtruncnorm(geno,-1,1,randm[i],0.3)
  
}

N_uni=N_uni/rowSums(N_uni)
N_gau=N_gau/rowSums(N_gau)


#Case 1: Uniformly distributed traits + truncated gaussian comp. kernel

res=single_sim1(5000,r,K2,K2,a1,Ag,R,N_uni,1000)
dat=getsum(res)

pops=dat[[1]]%>%as_tibble()

names(pops)=as.character(seq(nsp))
  
pplot=pops%>%
  mutate(time=1:n())%>%
  pivot_longer(1:nsp,names_to = "species",values_to = "N")%>%
  ggplot(aes(x=time,y=N,col=species))+
  geom_line()


trmeans=dat[[2]]%>%as_tibble()
names(trmeans)=as.character(seq(nsp))

tplot=trmeans%>%
  mutate(time=1:n())%>%
  pivot_longer(1:nsp,names_to = "species",values_to = "trait")%>%
  ggplot(aes(x=time,y=trait,col=species))+
  geom_line()

pplot|tplot



########################################################################
#Data analysis
#Load files



path1="C:/Users/mihir/Documents/Comp_Coevol/Codes"



files=list.files(path1)[grep("compdat",list.files(path1))]

dat=tibble()

for(i in 1:length(files)){
  
  dat=dat%>%
      bind_rows(
        read.csv(paste0(path1,"/",files[1]))
      )
  
}

#Variables- nloci,a1s,K1s,kernel,traits

dat%>%
  group_by(kernel,time)%>%
  summarize(means=mean(mnnd),
            sd=sd(mnnd))%>%
  ungroup()%>%
  ggplot(aes(x=time,y=means,col=kernel,fill=kernel))+
  geom_line()+
  geom_point()+
  geom_ribbon(aes(ymin=means-sd,ymax=means+sd),alpha=0.2)

dat%>%
  group_by(a1s,time)%>%
  summarize(means=mean(mnnd),
            sd=sd(mnnd))%>%
  ungroup()%>%
  mutate(a1s=as.factor(a1s))%>%
  ggplot(aes(x=time,y=means,col=a1s,fill=a1s))+
  geom_line()+
  geom_point()+
  geom_ribbon(aes(ymin=means-sd,ymax=means+sd),alpha=0.2)


###############################################################################
#OLD CODE
###############################################################################
#Error function
erfun<-function(x){
  return(2 * pnorm(x * sqrt(2)) - 1)
}


#Determine competitive kernels

#Gaussian kernel
alpha.g<-function(a,b,omega){
  if(is.na(a)==TRUE | is.na(b)==TRUE){return(0)}
  alpha<-exp(-((a-b)^2/(omega)^2))
  return(alpha)           
}

#Gaussian kernel with threshold (Truncated Gaussian). Very high value of threshold leads to a regular Gaussian kernel.

alpha.gt<-function(a,b,omega,t){
  if(is.na(a)==TRUE | is.na(b)==TRUE){return(0)}
  alpha<-exp(-((a-b)^2/(omega)^2))
  if(abs(a-b)>t){alpha<-0}
  return(alpha)           
}

#Triangle kernel with threshold t
alpha.tri<-function(a,b,slope,t){
  if(is.na(a)==TRUE | is.na(b)==TRUE){return(0)}
  if(abs(a-b)>t){return(alpha<-0)}
  if(((a-b)>-t) & ((a-b)<=0)){alpha<-(1+(slope*(a-b)))}
  if(((a-b)>0) & ((a-b)<t)){alpha<-(1-(slope*(a-b)))}
  return(alpha)
}


#Plot kernels 
a<-0
bs<-seq(-1,1,length.out=100)
alphas<-matrix(ncol=3,nrow=100)
omega<-0.5
t<-0.5
slope=2
for(i in 1:100){
  b<-bs[i]
  alphas[i,1]<-alpha.g(a,b,omega)
  alphas[i,2]<-alpha.gt(a,b,omega,t)
  alphas[i,3]<-alpha.tri(a,b,slope,t)
}
plot(alphas[,1]~bs,type="l",ylim=c(-0.2,1.3))
lines(alphas[,2]~bs,col=2)
lines(alphas[,3]~bs,col=3)




#This function gives beta functions for three kernel functions for a given pair of trait distributions
beta.gt<-function(a,b,omega,t){
  mean.a<-mean(a)
  int<-0
  for(i in 1:100000){
    a1<-sample(a,1)
    b1<-sample(b,1)
    int<-int+((a1-mean.a)*alpha.gt(a1,b1,omega,t))
  }
  int<-int/10000
  return(int)
}


beta.tri<-function(a,b,slope,t){
  mean.a<-mean(a)
  int<-0
  for(i in 1:100000){
    a1<-sample(a,1)
    b1<-sample(b,1)
    int<-int+((a1-mean.a)*alpha.tri(a1,b1,slope,t))
  }
  int<-int/10000
  return(int)
}


dat=NULL

pars=expand.grid(s1=c(0.1,0.5,1,5,10),s2=c(0.1,0.5,1,5,10))

for(i in 1:nrow(pars)){
  
  dat=bind_rows(dat,tibble(par1=pars[i,1],par2=pars[i,2],ind=i,vals=rbeta(1000,pars[i,1],pars[i,2])))
  
}

dat%>%
ggplot(aes(vals))+
geom_density()+
facet_wrap(vars(ind),scales="free")  


#Scenario 1: Flat trait distributions 
means=seq(-4,4,length.out=50)
a=runif(1000,-1,1)

res1=NULL

for(i in 1:length(means)){
  
  b=runif(1000,means[i]-1,means[i]+1)
  
  res1=bind_rows(res1,tibble(mean=means[i],
                       beta_gt=beta.gt(a,b,0.5,0.5),
                       beta_tri=beta.tri(a,b,2,0.5)))
  
}



#Scenario 2: Normally distributed traits

means=seq(-4,4,length.out=50)
a=rtruncnorm(1000,-5,5,0,1)

res=NULL

for(i in 1:length(means)){
  
  b=rtruncnorm(1000,-5,5,means[i],1)
  
  res=bind_rows(res,tibble(mean=means[i],
                           beta_gt=beta.gt(a,b,0.5,0.5),
                           beta_tri=beta.tri(a,b,2,0.5)))
  
}

res%>%
  ggplot()+
  geom_line(aes(mean,beta_gt,col="red"))+
  geom_line(aes(mean,beta_tri,col="green"))+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)


results=bind_rows(
  res1%>%
    pivot_longer(cols=beta_gt:beta_tri,names_to = "kernel",values_to = "beta")%>%
    mutate(traits="flat"),
  res%>%
    pivot_longer(cols=beta_gt:beta_tri,names_to = "kernel",values_to = "beta")%>%
    mutate(traits="normal"))

results%>%
  ggplot(aes(mean,beta,col=kernel))+
  geom_line()+
  geom_hline(yintercept=0)+
  geom_vline(xintercept=0)+
  facet_wrap(~traits)

#plot beta function for normal trait distributions (truncated) and truncated Gaussian kernel
m=seq(-5,5,length.out=200)
sdval=1
omega=1
t=2

int=0
a=rtruncnorm(10^4,-5,5,m[i],sdval)
b=rtruncnorm(10^4,-5,5,0,sdval)

betas=NULL

for (i in 1:length(m)){
  
  int=0
  
for (j in 1:length(a)){
  
  a1<-runif(1,-5,5)
  b1<-runif(1,m[i]-2.5,m[i]+2.5)
  int<-int+(a1*alpha.tri(a1,b1,t))  
  
  
}
  
  betas=c(betas,int/length(a))
  
}


plot(-betas~m,type="l")





#plot beta function for different trait distribution shapes. Functions for trait distributions:
#1. Normal, 2. Uniform 3. Multimodal
betas<-c()

t<-1
omega<-0.5

t1<-rnorm(100,0,1)
t1<-runif(100,-5,5)
t1<-c(rnorm(50,0,1)+rnorm(50,3,1))


a.mean<-mean(t1)
bs<-matrix(nrow=100,ncol=100)
for(i in 1:nrow(bs)){
  bs[i,]<-c(rnorm(50,sample((-5:5),1),1),rnorm(50,sample((-5:5),1),1))
}

b.means<-apply(bs,1,mean)

for(i in 1:length(b.means)){
 # b.mean<-b.means[i]
  #t2<-rtruncnorm(500,-1,1,b.mean,0.5)
  #t2<-runif(100,b.mean-5,b.mean+5)
  #t2<-rnorm(100,b.mean,1)
  t2<-bs[i]
  
  beta.g<-c(beta.g,beta.num(t1,t2)[1])
  beta.gt<-c(beta.gt,beta.num(t1,t2)[2])
  beta.tri<-c(beta.tri,beta.num(t1,t2)[3])
}
plot(beta.gt~b.means,ylim=c(-0.2,0.2),cex.lab=1.5,ylab="",
     xlab=expression(mu[i]*-mu[j]),main="Truncated Gaussian kernel")
title(ylab=expression(beta[ij]),line=1.9, cex.lab=1.5)
abline(h=0,lty=2)
abline(v=0,lty=2)
legend("topleft",bty="n",legend=c(expression('z'[i]*'~ U(-5,5)'),expression('z'[j]*'~ U('*mu[j]*'- 5,'*mu[j]*'+ 5)')))

legend("topright",bty="n",legend=c(expression('z'[i]*'~ N('*mu[i]*',0.5)'),expression('z'[j]*'~ N('*mu[j]*',0.5)')))
legend("topleft",bty="n",legend=c("Both species have bimodal trait distributions"))




t1<-rnorm(1000,-5,1)
t2<-t1[t1>-5]
t2<-t2[t2<5]

ss<-seq(-5,0,by=0.3)
dff<-c()
for(i in 1:length(ss)){
  s1<-rnorm(1000,ss[i],1)
  s2<-s1[s1>-5]
  s2<-s2[s2<5]
  dff<-c(dff,(beta.num(t1,s1,omega,t)[2]-beta.num(t2,s2,omega,t)[2]))
}

plot(dff~ss,type="l")




##############
#Draw a generic beta function
b<-seq(-8,5,length.out=200)
bet1<-vector(length=200)
for(i in 1:length(b)){
  a1<-runif(500,-5,5)
  b1<-runif(500,b[i],b[i]+3)
  bet1[i]<--beta.num(a1,b1,0.5,0.5)
  
}
b0<-b+1.5
plot(bet1~b0,type="l")
#,axes=FALSE,ann=FALSE)

p1<-plot(bet1~b,type="l",axes=FALSE,ann=FALSE)

tiff('beta_num.tiff', units="in", width=10, height=5, res=300)
plot(bet1~b,type="l",axes=FALSE,ann=FALSE)
dev.off()


######################################################################
# New analysis
#Estimate the conditions under which species may exhibit convergent evolution
#Represent this probability in terms of parameters of competition kernel, and trait 
#distributions







