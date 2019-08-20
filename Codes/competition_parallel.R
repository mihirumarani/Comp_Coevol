library(igraph)


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
  max(fgreedy$modularity)
  return(max(fgreedy$modularity))
}



mods2<-function(a,omega,t){
  a<-na.omit(a)
  net<-network(a,omega,t)
  G<-as.undirected(graph.adjacency(net,mode="undirected",weighted=TRUE,diag=FALSE))
  fgreedy<-fastgreedy.community(G,merges=TRUE, modularity=TRUE)
  max(fgreedy$modularity)
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
  

al2<-function(t1,t2,v1,v2){
  a<-0
  for(i in 1:10^5){
    a<-a+alpha1(rnorm(1,t1,v1),rnorm(1,t2,v2))
  }
  
  return(a/10^5)
}


#beta function to see the total trait change from pairwise interaction
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

  #No of species +param values
n<-50
h<-rep(0.5,n)
omegas<-c(0.05,0.1,0.25,0.5,1)
tvars<-rbind(c(0.5,1.5),c(2.5,3.5),c(0.5,3.5))
ts<-c(5,7,10,200)


#############################################
#Start the cluster and split the loop on 8 cores
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

#Create a set of 100 initializations for tmean, N, r and K
tmeanset<-matrix(rnorm(50*n,0,5),nrow=100,ncol=n)
#Nset<-matrix(rep(abs(rnorm(50,1000,200)),n),ncol=n)
Nset<-matrix(abs(rnorm(50*n,1000,200)),nrow=100,ncol=n)
Kset<-2*Nset
#rset<-matrix(rep(abs(rnorm(50,0,0.5)),n),ncol=n)
rset<-matrix(abs(rnorm(50*n,0,0.5)),nrow=100,ncol=n)


modsamples<-ceiling(c(seq(1,2000,length.out=20),seq(2100,49999,length.out=30)))
numcores<-detectCores()
clust<-makeCluster(numcores)
clusterExport(clust,c("trait","pops","be","al","alpha1","n","h","omegas","tvars",
                      "ts","mods2","modsamples","tmeanset",
                      "Nset","Kset","rset","network","mnnd","sdnnr"))
registerDoParallel(clust)


multiResultClass <- function(result1=NULL,result2=NULL,result3=NULL)
{
  me <- list(
    result1 = result1,
    result2 = result2,
    result3 = result3
  )
  
  ## Set the name for the class
  class(me) <- append(class(me),"multiResultClass")
  return(me)
}



dat<-foreach(i=1:50,.packages="igraph",.combine=rbind)%dopar%
{
  
  paramdat<-c(length=4)
  #modat<-c(length=50)
  trdat<-c(length=1000)
  popdat<-c(length=1000)
  for(p in 1:length(ts)){
    t<-ts[p]
  for(k in 1:5){
    omega<-omegas[k]
      for(l in 1:3){
          tvard<-tvars[l,]
          tvar<-runif(n,tvard[1],tvard[2])
          tmean<-tmeanset[i,]
          N<-Nset[i,]
          K<-Kset[i,]
          r<-rset[i,]
          Nnew<-N
          tmeannew<-tmean
          mod1<-c()
          popdat1<-c()
          trdat1<-c()
          for(m in 1:50000){
            tmeannew1<-trait(Nnew,h,r,K,tmeannew,tvar,omega,t)
            Nnew<-pops(Nnew,K,r,tmeannew1,tvar,omega,t)
            if(m%in%modsamples){
             # mod1<-c(mod1,mods2(tmeannew,omega,t))
             # mds1<-c(mds1,mnnd(tmeannew))
              popdat1<-c(popdat1,Nnew)
              trdat1<-c(trdat1,tmeannew)
            }
            tmeannew<-tmeannew1
          }
          paramdat<-rbind(paramdat,c(t,tvard,omega))
         #modat<-rbind(modat,mod1)
         # mds<-rbind(mds,mds1)
          popdat<-rbind(popdat,popdat1)
          trdat<-rbind(trdat,trdat1)
        }
    }
  }
  result <- multiResultClass()
  result$result1<-paramdat
  result$result2<-trdat
  result$result3<-popdat
  
  return(result)
}

stopCluster(clust)

paramdats<-dat[,1]
paramdats1<-do.call(rbind,paramdats)
paramdats<-paramdats1[-which(paramdats1[,1]==4),]
paramdats<-as.matrix(paramdats)
write.csv(paramdats,file="paramdatsfull50.csv")

trdats<-dat[,2]
trdats1<-do.call(rbind,trdats)
trdats<-trdats1[-which(trdats1[,1]==1000),]
trdats<-as.matrix(trdats)
write.csv(trdats,file="trdatsfull50.csv")

popdats<-dat[,3]
popdats1<-do.call(rbind,popdats)
popdats<-popdats1[-which(popdats1[,1]==1000),]
popdats<-as.matrix(popdats)
write.csv(popdats,file="popdatsfull.csv")

mds<-dat[,4]
mds1<-do.call(rbind,mds)
mds<-mds1[-which(mds1[,1]==50),]
mds<-as.matrix(mds)
write.csv(mds,file="mnndsfull.csv")


###############################




####Clean data###Final set should have #omegasX#hsX100 rows
rawdat<-read.csv("motherlode1.csv")
rawdat<-rawdat[-1,]
rawdat<-rawdat[,-1]
rawdat1<-rawdat[-which(rawdat[,1]==82),]

####column info: 1=omega, 2=h, 3:22=initial pop., 23:42=final pop, 43:62=initial traits, 63:82=final traits

#plot results: Central tendancies initial vs. final
initmean<-apply(rawdat1[,43:62],1,mean)
finalmean<-apply(rawdat1[,63:82],1,mean)
initvar<-apply(rawdat1[,43:62],1,var)
finalvar<-apply(rawdat1[,63:82],1,var)
initmnnd<-apply(rawdat1[,43:62],1,mnnd)
finalmnnd<-apply(rawdat1[,63:82],1,mnnd)
initsdnnr<-apply(rawdat1[,43:62],1,sdnnr)
finalsdnnr<-apply(rawdat1[,63:82],1,sdnnr)


meansdat<-data.frame(means=c(rep("Initial",length(initmean)),rep("Final",length(finalmean))),replicates=c(initmean,finalmean))
boxplot(meansdat$replicates~meansdat$means)
varsdat<-data.frame(Variance=c(rep("Initial",length(initvar)),rep("Final",length(finalvar))),replicates=c(initvar,finalvar))
boxplot(varsdat$replicates~varsdat$Variance)
mnnddat<-data.frame(Mean_Neighbor_dist.=c(rep("Initial",length(initmnnd)),rep("Final",length(finalmnnd))),replicates=c(initmnnd,finalmnnd))
boxplot(mnnddat$replicates~mnnddat$Mean_Neighbor_dist.)



#Calculate network measures: Modularity, nestedness?, connectivity :This is highly iffy!
net1<-network(data2[1,])
diag(net1)<-0
graph1<-graph_from_adjacency_matrix(net1)
graph1<-as.undirected(graph1)
V(graph1)$size<-20
E(graph1)$arrow.size<-0.2

plot(graph1)
plot(density(na.omit(data2[1,])))
cluster_edge_betweenness(graph1)





#Population distribution vs. traits values/ranks
traitdat<-rawdat[,23:42]
for(i in 1:nrow(traitdat)){
  traitdat[i,]<-traitdat[1,][order(traitdat[1,],decreasing=TRUE)]
}


#Heredity values
h<-rep(0.5,n)

#Breadth of competition
omega<-0.1

#Trait values
#Need to specify means and one variance value
tmean<-rnorm(n,0,5)

tvar<-abs(rnorm(n,3.0))


#Demographic parameters
N<-abs(rnorm(n,1000,200))
K<-2*N
r<-abs(rnorm(n,1,0.5))



#Simulation
 
  data1<-N
  data2<-tmean
  Nnew<-N
  tmeannew<-tmean
  
  #net1<-network(tmeannew)
  #modmat[j,1]<-mods1(net1)
  
 # ranges<-c(-20,20)
  #binsize<-3
  #bins<-ceiling((ranges[2]-ranges[1])/binsize)
  #memberships1<-matrix(nrow=1001,ncol=n)
  #memberships<-matrix(nrow=1001,ncol=ceiling((ranges[2]-ranges[1])/3))
  for(i in 1:10000){
    tmeannew<-trait(Nnew,h,r,K,tmeannew,tvar,omega)
    Nnew1<-pops(Nnew,K,r,tmeannew,tvar,omega)
    
    Nnew<-Nnew1
    data1<-rbind(data1,Nnew1)
    data2<-rbind(data2,tmeannew)
    
   
    #Track membership
   # membs<-c()
    #a<-tmeannew
    #for(k in 1:length(a)){
     # if(is.na(a[k])==TRUE){membs<-c(membs,0)}
      #res<-ceiling((a[k]-ranges[1])/binsize)
      #membs<-c(membs,res)
    #}
    #memberships1[i,]<-membs
    #freq<-c()
    
    #for(l in 1:bins){
     # freq<-c(freq,sum(membs==l))
    #}
    #memberships[i,]<-freq
  }

#plots
#Population trend

plot(data1[,1],type="l",ylim=c(min(data1,na.rm=TRUE),max(data1,na.rm=TRUE)),ylab="N",xlab="Time")
for(i in 2:ncol(data1)){
  lines(data1[,i],col=i)
}

#Trait data
#tiff("PlotN4.tiff", width m= 5, height = 5, units = 'in', res = 300)

plot(data2[,1],type="l",ylab="trait value",ylim=c(min(data2,na.rm=TRUE),max(data2,na.rm=TRUE)+2),xlab="time")
for(i in 2:ncol(data2)){
  lines(data2[,i],col=i)
}
 legend("topleft",c("a=5","omega=2"))
dev.off()


#Plots traits + abundance (circle size)

plot((rep(3,n))~data2[1001,],cex=(data1[1001,]/500),ylim=c(-2,5),col=1:n,lwd=2)
points((rep(-1,n))~data2[1,],cex=(data1[1,]/500),col=1:n,lwd=2)
points((rep(1,n))~data2[500,],cex=(data1[500,]/500),col=1:n,lwd=2)

plot(apply(data2,1,mean),ylim=c(-10,10),xlab="time",ylab="Grand mean")
plot(apply(data2,1,var),xlab="time",ylab="Grand Variance")




#Build a jacabian 
jacob<-function(N,t){
  tr<-runif(100,-15,15)
  effect<-c()
  v1<-3
  v2<-3
  for(i in 1:100){
    a<-0
    for(j in 1:length(N)){
      a<-a+al(tr[i],t[j],v1,v2,omega)*N[j]
    }
  effect<-c(effect,a)
  }
  return(cbind(tr,effect))
}


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





#Plot changes in memberships
plot(memberships[,1],type="l",ylim=c(0,20))
for(i in 2:ncol(memberships)){
  lines(memberships[,i],col=i)
}
#rates
rates<-data1[2:1001,1]/data1[1:1000,1]
plot(rates,type="l")

x<-1:nrow(data3)
s<-seq(length(x)-1)
cols<-grey.colors(nrow(data3),0.3,0.9)
plot(x,data3[,1],xlab="time steps",ylab="trait value",ylim=c(min(data3),max(data3)),type="l")
segments(x[s],data3[,1][s],x[s+1],data3[,1][s+1],lwd=popdat[1:(nrow(popdat)-1),1]/500,col=cols[1])
for(i in 2:nrow(data3)){
  segments(x[s],data3[,i][s],x[s+1],data3[,i][s+1],lwd=popdat[1:(nrow(popdat)-1),i]/500,col=cols[i])
}






lines(means,lty=8, lwd=3)
plot(means)
plot(vars)
##############################################
#For connectivity
#Plot a network based on traits
network<-function(a){
  a<-na.omit(a)
  mat<-matrix(ncol=length(a),nrow=length(a))
  samples<-1:length(a)
  for(i in samples){
    for(j in samples){
      mat[i,j]<-alpha1(a[i],a[j],1)
    }
  }
  list<-is.na(mat)
  mat[list]<-0
  mat[which(mat>0)]<-1
  diag(mat)<-0
  return(mat)
}


mynet<-network(tmean)
diag(mynet)<-0
mynet[which(mynet>0)]<-1

mygraph<-graph_from_adjacency_matrix(mynet)
plot(mygraph, edge.arrow.size=.2, edge.color="orange",
     vertex.color="orange", vertex.frame.color="#ffffff",
      vertex.label.color="black")

#Plot initial and final networks
net1<-network(data2[1,])
diag(net1)<-0
graph1<-graph_from_adjacency_matrix(net1)
graph1<-as.undirected(graph1)
V(graph1)$size<-20
E(graph1)$arrow.size<-0.2

plot(graph1)
plot(density(na.omit(data2[1,])))
cluster_edge_betweenness(graph1)
spectralOptimization(net1)

#Final Network
net2<-network(na.omit(data2[nrow(data2),]))
diag(net2)<-0
graph2<-graph_from_adjacency_matrix(net2)
#Plot the network
V(graph2)$size<-20
E(graph2)$arrow.size<-0.2
plot(graph2)
#plot trait distribution of extant species
plot(density(na.omit(data2[1001,])))
#Modules and memberships of networks
membership(cluster_edge_betweenness(as.undirected(graph2)))
spectralOptimization(net2)

g <- sample_pa(100, m = 2, directed = FALSE)


g <- make_full_graph(10) %du% make_full_graph(10)
g <- add_edges(g, c(1,11))
eb <- cluster_edge_betweenness(g)
eb


plot(graph2)

#Plot the evolving network
#Step 1: Find the modules and memberships for each step
#Step 2: Assign colors according to membership
samp<-c(1,50*(1:floor(nrow(data2)/50)))


M<-vector(length=n)
for(i in samp){
  net<-network(data2[i,])
  diag(net)<-0
  net[which(net>0)]<-1
  g<-graph_from_adjacency_matrix(net)
  g<-as.undirected(g,mode="collapse")
  member<-edge.betweenness.community(g)$membership
  M<-rbind(M,member)
}

M<-M[-1,]

  

tr<-data2[samp,]

#plot trait data from sampled time points
x<-1:nrow(tr)
s<-seq(length(x)-1)
plot(x,tr[,1],xlab="time",ylab="trait value",ylim=c(min(data2),max(data2)),type="l"))
segments(x[s],tr[,1][s],x[s+1],tr[,1][s+1],col=M[1:(nrow(tr)-1),1])
for(i in 2:nrow(tr)){
  segments(x[s],tr[,i][s],x[s+1],tr[,i][s+1],col=M[1:(nrow(tr)-1),i])
}



  nz<-function(a){
  return(length(which(a!=0)))
}

con<-matrix(ncol=n,nrow=nrow(data2))


for(i in 1:nrow(data2)){
  mat1<-network(data2[i,])
  con[i,]<-apply(mat1,2,nz)
}

plot(con[,1],type="l")
for(i in 2:n){
  lines(con[,i],col=i)
}




 


