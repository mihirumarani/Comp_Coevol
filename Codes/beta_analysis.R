#Functions
library(truncnorm)
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
alpha.tri<-function(a,b,t){
  if(is.na(a)==TRUE | is.na(b)==TRUE){return(0)}
  if(abs(a-b)>t){return(alpha<-0)}
  slope<-1/t
  if(((a-b)>-t) & ((a-b)<=0)){alpha<-(1+(slope*(a-b)))}
  if(((a-b)>0) & ((a-b)<t)){alpha<-(1-(slope*(a-b)))}
  return(alpha)
}


#Plot kernels 
a<-0
bs<-seq(-5,5,length.out=100)
alphas<-matrix(ncol=3,nrow=100)
omega<-0.1
t<-0.25
for(i in 1:100){
  b<-bs[i]
  alphas[i,1]<-alpha.g(a,b,omega)
  alphas[i,2]<-alpha.gt(a,b,omega,t)
  alphas[i,3]<-alpha.tri(a,b,t)
}
plot(alphas[,1]~bs,type="l",ylim=c(-0.2,1.3))
lines(alphas[,2]~bs,col=2)
lines(alphas[,3]~bs,col=3)


#This function gives beta functions for three kernel functions for a given pair of trait distributions
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
