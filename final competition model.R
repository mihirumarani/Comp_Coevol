library(igraph)
library(reshape2)


#Functions
#functions for network parameters
network01<-function(a,omega,t){
  a<-na.omit(a)
  mat<-matrix(ncol=length(a),nrow=length(a))
  samples<-1:length(a)
  for(i in samples){
    for(j in samples){
      if(alpha1(a[i],a[j],omega,t)>0){
        mat[i,j]<-1
      }
      else{mat[i,j]<-0}
    }
  }
  list<-is.na(mat)
  mat[list]<-0
  diag(mat)<-0
  return(mat)
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



#Calculate modularity
mods<-function(a,omega,t){
  a<-na.omit(a)
  net<-network01(a,omega,t)
  G<-as.undirected(graph_from_adjacency_matrix(net))
  fgreedy<-fastgreedy.community(G,merges=TRUE, modularity=TRUE)
  return(max(fgreedy$modularity))
}

mod2<-function(a,omega,t){
  a<-na.omit(a)
  net<-network(a,omega,t)
  G<-graph_from_adjacency_matrix(net,mode="undirected",weighted=TRUE,diag=FALSE)
  fgreedy<-fastgreedy.community(G,merges=TRUE, modularity=TRUE)
  return(max(fgreedy$modularity))
}



#Nearest neighbot distance
mnnd<-function(a){
  a<-na.omit(a)
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

#Nearest neighbor st.deviation
sdnnr<-function(a){
  a<-a[order(a,decreasing=TRUE)]
  nnd<-c(length=length(a))
  nnd[1]<-a[1]-a[2]
  nnd[length(a)]<-a[length(a)-1]-a[length(a)]
  for(i in 2:(length(a)-1)){
    nnd[i]<-min((a[i-1]-a[i]),(a[i]-a[i+1]))
  }
  ranges<-a[1]-a[length(a)]
  sdnnr<-sqrt(sum((nnd-mean(nnd))^2)/(length(a)-1))
  return(sdnnr/ranges)
}

#Determine alpha values

#Symmetric function
alpha1<-function(a,b,omega,t){
  if(is.na(a)==TRUE | is.na(b)==TRUE){return(0)}
  alpha<-exp(-((a-b)^2/(omega)^2))
  if(abs(a-b)>t){alpha<-0}
  return(alpha)           
}

#Asymmetric function
alpha2<-function(a,b,o){
  alpha<-1/(1+exp(-(a-b)/o))
  return(alpha)
}



#al function to calculate the total rate change from pairwise interaction
#Analytical function taken from Barabasi and D'Andrea (2016)
al<-function(t1,t2,v1,v2,omega,t){
  if(is.na(t1-t2)){al<-0}
  else{
  if(abs(t1-t2)>t){al<-0}
  else{
  den<-sqrt(2*v1^2+2*v2^2+omega^2)
  al<-(omega/den)*exp(-((t1-t2)^2)/(den^2))
  }
  }
  return(al)
  }
  

#Numerical solution to the function from Barabasi and D'Andrea
al2<-function(t1,t2,v1,v2,omega,t){
  a<-0
  for(i in 1:10^5){
    a<-a+alpha1(rnorm(1,t1,v1),rnorm(1,t2,v2),omega,t)
  }
  
  return(a/10^5)
}




#beta function to see the total trait change from pairwise interaction taken from Barabasi and D'Andrea (2016)
be<-function(t1,t2,v1,v2,omega,t){
  if(is.na(t1-t2)){be<-0}
  else{
  if(abs(t1-t2)>t){be<-0}
  else{
  den<-sqrt(2*(v1^2)+2*(v2^2)+omega^2)
  be<-((2*omega*(v1^2)*(t2-t1))/(den^3))*exp(-((t1-t2)^2)/(den^2))
  }
  }
  return(be)
}

#Numerical solution of beta function
be2<-function(t1,t2,v1,v2,omega,t){
  a<-0
  for(i in 1:10^5){
    t<-rnorm(1,t1,v1)
    a<-a+((t-t1)*alpha1(t,rnorm(1,t2,v2),omega,t))
  }
  return(a/10^5)
}

#Pop change in a generation for a network

pops<-function(N,K,r,tmean,tvar,omega,t){
  pop1<-c()
  N[which(N<=0)]<-0
  for(i in 1:length(N)){
    int<-0
    for(j in (1:length(N))[-i]){
      a<-al(tmean[i],tmean[j],tvar[i],tvar[j],omega,t)*N[j]
      int<-int+a
    }
    pop1[i]<-N[i]+N[i]*(r[i]*(1-((N[i]+int)/K[i])))
  }
  pop1[which(pop1<=0)]<-0
  return(pop1)
}


#Trait change in a generation for a network

trait<-function(N,h,r,K,tmean,tvar,omega,t){
  tmean[which(N<=0)]<-NA
  N[which(N<=0)]<-0
  traits<-c()
  for(i in 1:length(N)){
    int<-0
    for(j in (1:length(N))[-i]){
      a<-be(tmean[i],tmean[j],tvar[i],tvar[j],omega,t)*N[j]
      int<-int+a
    }
    traits[i]<-tmean[i] - ((h[i]^2)*r[i]*int/K[i])
  }
  return(traits)
}


#initiate simulation

n<-20

#Heredity values
h<-rep(0.5,n)

#Breadth of competition
omega<-0.25

#Trait values
#Need to specify means and one variance value
tmean2<-runif(n,-8.530317,  8.312025)
tmean1<-rnorm(n,0,5)
tmean<-tmean1
tmean<-tmean[order(tmean,decreasing = FALSE)]
tvar<-abs(rnorm(n,3,1))


#Demographic parameters
#N<-500+1000*abs(tmean/max(tmean))
#N<-1500-1000*abs(tmean/max(tmean))

#N<-rep(abs(rnorm(1,1000,200)),n)
N<-abs(rnorm(n,1000,200))
K<-2*N
#r<-rep(abs(rnorm(1,1,0.5)),n)
r<-abs(rnorm(n,1,0.5))
t<-5



#Simulation
 
  data1<-N
  data2<-tmean
  Nnew<-N
  tmeannew<-tmean
  

  for(i in 1:10000){
    tmeannew<-trait(Nnew,h,r,K,tmeannew,tvar,omega,t)
    Nnew1<-pops(Nnew,K,r,tmeannew,tvar,omega,t)
    Nnew<-Nnew1
    data1<-rbind(data1,Nnew1)
    data2<-rbind(data2,tmeannew)
  }

#plots
#Population trend

plot(data1[,1],type="l",ylim=c(min(data1,na.rm=TRUE),max(data1,na.rm=TRUE)),col=cols[1],ylab="N",xlab="Time",lwd=3)
for(i in 2:ncol(data1)){
  lines(data1[,i],col=cols[i],lwd=3)
}


#Trait evolution 
#tiff("PlotN4.tiff", width m= 5, height = 5, units = 'in', res = 300)
finalcols<-cols
color<-grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
cols1<-sample(color,20,replace=FALSE)
cols<-cols1


#tiff('plot5.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')

plot(data2[,1],type="l",ylab="trait value",ylim=c(min(data2,na.rm=TRUE),max(data2,na.rm=TRUE)+2),xlab="time",col=cols[1],lwd=3)
for(i in 2:ncol(data2)){
  lines(data2[,i],col=cols[i],lwd=3)
}
#dev.off()


#Plots traits + abundance (circle size)
 
data3<-data2[c(1,5000,10000),] 
#tiff('plot6.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')


plot(data3[,1]~c(1,5000,10000),type="o",cex=(data1[c(1,5000,10000),1]/500),lwd=3,ylim=c(min(data3),max(data3)),xlab="time",ylab="Trait values",col=cols[1])
for(i in 2:n){
  lines(data3[,i]~c(1,5000,10000),type="o",pch=1,cex=(data1[c(1,5000,10000),i]/500),lwd=3,col=cols[i])
}
#dev.off()

plot(rep(3,n)~data2[1000,],cex=(data1[1001,]/500),ylim=c(-2,5),col=1:n,lwd=2,xlim=c(-30,30))
points((rep(-1,n))~data2[1,],cex=(data1[1,]/500),col=1:n,lwd=2)
points((rep(1,n))~data2[5000,],cex=(data1[500,]/500),col=1:n,lwd=2)

#Plot trends in grand mean and grand variance of all species
plot(apply(data2,1,mean),ylim=c(-10,10),xlab="time",ylab="Grand mean")
plot(apply(data2,1,var),xlab="time",ylab="Grand Variance")


net1<-network(data3[1,],0.25,5)
g<-graph.adjacency(net1,weighted=TRUE,mode="undirected")
cl<-fastgreedy.community((g))
weights <- ifelse(crossing(cl, g), 1, 100)
l<-layout_with_kk(g,weights=weights)
V(g)$size <- data1[1,]/50
V(g)$label<-NA
V(g)$frame.color <- "white"
V(g)$color <- cols
E(g)$arrow.mode <- 0
#tiff('plot7.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')
plot(g,edge.width=E(g)$weight,layout=l)
#dev.off()


net1<-network(data3[2,],0.25,5)
g<-graph.adjacency(net1,weighted=TRUE,mode="undirected")
V(g)$size <- data1[5000,]/80
V(g)$frame.color <- "white"
V(g)$color <- cols
V(g)$label<-NA
E(g)$arrow.mode <- 0
tiff('plot8.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')
plot(g,edge.width=E(g)$weight/50,layout=l)
dev.off()

net1<-network(data3[3,],0.25,5)
g<-graph.adjacency(net1,weighted=TRUE,mode="undirected")
V(g)$size <- data1[10000,]/80
V(g)$frame.color <- "white"
V(g)$color <- cols
V(g)$label<-NA
E(g)$arrow.mode <- 0
tiff('plot10.tiff', units="in", width=5, height=4, res=300, compression = 'lzw')
plot(g,edge.width=E(g)$weight/50,layout=l)
dev.off()





#############################
#Run the trait dynamic without population dynamic
###Alternative model of trait evolution
#For the interaction component, trait change is linearly proportional to  trait matching function with thresholds
#Inputs: z= vector of all trait values, S=coefficient of rate change
traitevol<-function(z,omega,t){
  zdash<-c()
  for(i in 1:length(z)){
    x<-z-z[i]
    x[which(abs(x)>=t)]<-0
    y<-sum(-sign(x)*exp(-abs(x)/omega))
    zdash<-c(zdash,z[i]+y)
  }
  return(zdash)
}



#Simulate the process over multiple time steps
n<-20
traits<-tmean
primer<-traits
newdat<-traits
for(i in 1:10000){
  newdat<-traitevol(newdat,0.25,5)
  traits<-rbind(traits,newdat)
}

plot(traits[,1],type="l",ylim=c(min(na.omit(traits)),max(na.omit(traits))))
for(i in 2:20){
  lines(traits[,i],col=i)
}

dat<-traits[modsamples,]
modnull<-apply(dat,1,mods,omega=0.25,t=5)
mdnull<-apply(dat,1,mnnd)
plot(modnull~modsamples,ylim=c(0,1),type="l")
lines(mod1~modsamples,col=2)
lines(mod2~modsamples,col=3)

plot(mdnull~modsamples,ylim=c(0,1),type="l")
lines(md1~modsamples,col=2)
lines(md2~modsamples,col=3)



mat0<-jacob(data1[1,],data2[1,])
mat0<-mat0[order(mat0[,1],decreasing=TRUE),]
plot(abs(mat0[,2])~mat0[,1],xlab="Trait value", ylab="Per capita effect of competition",type="l")
points(rep(0,length(data2[1,]))~data2[1,],col=2,cex=(data1[1,]/500))

mat1<-jacob(data1[101,],data2[101,])
mat1<-mat1[order(mat1[,1],decreasing=TRUE),]
plot(abs(mat1[,2])~mat1[,1],xlab="Trait value", ylab="Per capita effect of competition",type="l",ylim=c(-50,max(mat1[,2])))
points(rep(0,length(data2[101,]))~data2[101,],col=2,cex=(data1[101,]/500))

mat2<-jacob(data1[501,],data2[501,])
mat2<-mat2[order(mat2[,1],decreasing=TRUE),]
plot(abs(mat2[,2])~mat2[,1],xlab="Trait value", ylab="Per capita effect of competition",type="l",ylim=c(-50,max(mat2[,2])))
points(rep(0,length(data2[501,]))~data2[501,],col=2,cex=(data1[501,]/500))

mat3<-jacob(data1[801,],data2[801,])
mat3<-mat3[order(mat3[,1],decreasing=TRUE),]
plot(abs(mat3[,2])~mat3[,1],xlab="Trait value", ylab="Per capita effect of competition",type="l",ylim=c(-50,max(mat3[,2])))
points(rep(0,length(data2[101,]))~data2[101,],col=2,cex=(data1[101,]/500))

mat4<-jacob(data1[1001,],data2[1001,])
mat4<-mat4[order(mat4[,1],decreasing=TRUE),]
plot(abs(mat4[,2])~mat4[,1],xlab="Trait value", ylab="Per capita effect of competition",type="l",ylim=c(-50,max(mat4[,2])))
points(rep(0,length(data2[1001,]))~data2[1001,],col=2,cex=(data1[1001,]/500))


###########################################


#Fitness landscape of an individual with trait t and variance tvar when competing with 20 species.
fits<-function(traitm,traitv,N,tmean,tvar,omega,t){
  values<-c()
  for(i in 1:length(traitm)){
    int<-0
    for(j in 1:length(N)){
      int<-int+N[j]*al(traitm[i],tmean[j],traitv,tvar[j],omega,t)
    }
    values<-c(values,int)
  }
  return(values)
}

#partial derivative of dN/dt wrt trait mean of each species 
jacob.t<-function(N,K,tmean,tvar,omega,t){
  values<-c()
  for(i in 1:length(tmean)){
    int<-0
    for(j in 1:length(N)){
      int<-int+N[j]*(al(tmean[i],tmean[j],tvar[i],tvar[j],omega,t)*(2*(tmean[j]-tmean[i])/(2*(tvar[i]^2)+2*(tvar[j]^2)+(omega^2))))
    }
    values<-c(values,-N[i]*int/K[i])
  }
  return(values)
}

#partial derivative of dN/dt wrt population size of each species 
jacob.n<-function(N,K,r,tmean,tvar,omega,t){
  values<-c()
  for(i in 1:length(tmean)){
    int<-0
    for(j in 1:length(N)){
      int<-int+N[j]*(al(tmean[i],tmean[j],tvar[i],tvar[j],omega,t))
    }
    int<-r[i]-((2*N[i]+int)/K[i])
    values<-c(values,int)
  }
  return(values)
}



#initiate simulation

n<-20

#Heredity values
h<-rep(0.5,n)

#Breadth of competition
omega<-0.25

#reference trait means and variance
reft<-seq(-30,30,length.out=50)
refv<-3

#Trait values
#Need to specify means and one variance value
tmean<-rnorm(n,0,5)
tmean<-tmean[order(tmean,decreasing = FALSE)]
tvar<-abs(rnorm(n,3,1))


#Demographic parameters
N<-abs(rnorm(n,1000,200))
K<-2*N
r<-abs(rnorm(n,1,0.5))
t<-5



#Simulation

data1<-N
data2<-tmean
data3<-c(length=length(reft))
data4<-c(length=length(tmean))
data5<-c(length=length(tmean))

Nnew<-N
tmeannew<-tmean

for(i in 1:10000){
  tmeannew<-trait(Nnew,h,r,K,tmeannew,tvar,omega,t)
  Nnew1<-pops(Nnew,K,r,tmeannew,tvar,omega,t)
  Nnew<-Nnew1
  data1<-rbind(data1,Nnew1)
  data2<-rbind(data2,tmeannew)
  data3<-rbind(data3,fits(reft,refv,Nnew,tmeannew,tvar,omega,t))
  data4<-rbind(data4,jacob.t(Nnew,K,tmeannew,tvar,omega,t))
  data5<-rbind(data5,jacob.n(Nnew,K,r,tmeannew,tvar,omega,t))
}

intervals<-c(2,500*(1:20))


mds<-apply(data2,1,mnnd)

data3<-max(data3)-data3




for(i in 1:length(intervals)){
  plot(data3[i,]~reft,type="l",ylim=c(0,max(data3[i,])),xlim=c(min(data2[i,]),max(data2[i,])))
  points(rep(0,20)~data2[i,],col="red",cex=data1[i,]/500,lwd=2)
}

plot(data3[2,]~reft,type="l",ylim=c(0,max(data3[2,])),xlim=c(-25,25),ylab="fitness", xlab="trait of incipient species")
points(rep(0,20)~data2[2,],col="red",cex=data1[2,]/500,lwd=2)

plot(data3[500,]~reft,type="l",ylim=c(0,max(data3[500,])),xlim=c(-25,25),ylab="fitness", xlab="trait of incipient species")
points(rep(0,20)~data2[500,],col="red",cex=data1[500,]/500,lwd=2)

plot(data3[1000,]~reft,type="l",ylim=c(0,max(data3[1000,])),xlim=c(-25,25),ylab="fitness", xlab="trait of incipient species")
points(rep(0,20)~data2[1000,],col="red",cex=data1[1000,]/500,lwd=2)

plot(data3[2000,]~reft,type="l",ylim=c(0,max(data3[2000,])),xlim=c(-25,25),ylab="fitness", xlab="trait of incipient species")
points(rep(0,20)~data2[2000,],col="red",cex=data1[2000,]/500,lwd=2)

plot(data3[4000,]~reft,type="l",ylim=c(0,max(data3[4000,])),xlim=c(-25,25),ylab="fitness", xlab="trait of incipient species")
points(rep(0,20)~data2[4000,],col="red",cex=data1[4000,]/500,lwd=2)

plot(data3[6000,]~reft,type="l",ylim=c(0,max(data3[6000,])),xlim=c(-25,25),ylab="fitness", xlab="trait of incipient species")
points(rep(0,20)~data2[6000,],col="red",cex=data1[6000,]/500,lwd=2)


plot(data3[8000,]~reft,type="l",ylim=c(0,max(data3[8000,])),xlim=c(-25,25),ylab="fitness", xlab="trait of incipient species")
points(rep(0,20)~data2[8000,],col="red",cex=data1[8000,]/500,lwd=2)


plot(data3[10000,]~reft,type="l",ylim=c(0,max(data3[10000,])),xlim=c(-25,25),ylab="fitness", xlab="trait of incipient species")
points(rep(0,20)~data2[10000,],col="red",cex=data1[10000,]/500,lwd=2)



plot(data2[,1],type="l",ylab="trait value",ylim=c(min(data2,na.rm=TRUE),max(data2,na.rm=TRUE)+2),xlab="time",col=cols[1],lwd=3)
for(i in 2:ncol(data2)){
  lines(data2[,i],col=cols[i],lwd=3)
}
abline(v=5000,lwd=3)
#Plot Jacobian

plot(data4[2,]~data2[2,],lwd=4,cex=data1[2,]/250,ylab="Change in fitness wrt trait increase",xlab="Trait mean",xlim=c(-25,25))
abline(h=0)

plot(data4[100,]~data2[100,],lwd=4,cex=data1[100,]/250,ylab="Change in fitness wrt trait increase",xlab="Trait mean",xlim=c(-25,25))
abline(h=0)

plot(data4[500,]~data2[500,],lwd=4,cex=data1[500,]/250,col=cols,ylab="Change in fitness wrt trait increase",xlab="Trait mean",xlim=c(-25,25),ylim=c(-30,30))
abline(h=0)

plot(data4[1000,]~data2[1000,],lwd=4,cex=data1[1000,]/250,col=cols,ylab="Change in fitness wrt trait increase",xlab="Trait mean",xlim=c(-25,25),ylim=c(-30,30))
abline(h=0)

plot(data4[2000,]~data2[2000,],lwd=4,cex=data1[2000,]/250,col=cols,ylab="Change in fitness wrt trait increase",xlab="Trait mean",xlim=c(-25,25),ylim=c(-30,30))
abline(h=0)

plot(data4[5000,]~data2[5000,],lwd=4,cex=data1[5000,]/250,col=cols,ylab="Change in fitness wrt trait increase",xlab="Trait mean",xlim=c(-25,25),ylim=c(-30,30))
abline(h=0)

plot(data4[10000,]~data2[10000,],lwd=4,cex=data1[10000,]/250,col=cols,ylab="Change in fitness wrt trait increase",xlab="Trait mean",xlim=c(-25,25),ylim=c(-30,30))
abline(h=0)





#Create sets of trait distributions with different MNND values and measure fitness change wrt trait and population sizes
trset<-matrix(rnorm(2000,0,10),nrow=100)

popset<-matrix(abs(rnorm(2000,1000,500)),nrow=100)

tset<-c(5,7,10,200)

tvarset<-rbind(c(0.5,1.5),c(2.5,3.5),c(0.5,3.5))

omegaset<-c(0.1,0.25,0.5,1,2)

reft<-seq(-30,30,length.out=50)
refv<-2

dataset<-vector(length=50)

  for(j in 1:length(tset)){
    for(k in 1:nrow(tvarset)){
      vars<-runif(20,tvarset[k,1],tvarset[k,2])
      for(l in 1:length(omegaset)){
        for(i in 1:100){
          dataset<-rbind(dataset,fits(reft,refv,popset[i,],trset[i,],vars,omegaset[l],tset[j]))
        }
      }
    }
  }
dataset<-dataset[-1,]
dataset<-max(dataset)-dataset


localMaxima <- function(x) {
  # Use -Inf instead if x is numeric (non-integer)
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

totalmax<-function(x){
  return(length(localMaxima(x)))
}

mnnds<-apply(trset,1,mnnd)
mds<-rep(mnnds,60)
p<-plot_ly(x=reft,y=mnnds,z=~dataset[1:100,]) %>% add_surface()
p


optimas<-apply(dataset,1,totalmax)


#plot local fitness optimas in a biotic environment of 20 competing species as a function of mnnd.
plot(optimas~mds)


