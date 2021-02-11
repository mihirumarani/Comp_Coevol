library(truncnorm)
library(igraph)

#Determine alpha values

erfun<-function(x){
  return(2 * pnorm(x * sqrt(2)) - 1)
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
  #mat<-mat/max(mat)
  return(mat)
}



mod<-function(a,tvar,omega,t){
  tvar<-tvar[!is.na(a)]
  a<-na.omit(a)
  tvar<-tvar
  if(length(a)>0){
    net<-network(a,tvar,omega,t)
    G<-graph_from_adjacency_matrix(net,mode="undirected",weighted=TRUE,diag=FALSE)
    fgreedy<-fastgreedy.community(G,merges=TRUE, modularity=TRUE)
    return(c(max(fgreedy$modularity),length(unique(fgreedy$membership))))}
  else{return(NA)}
}


al.ana<-function(t1,t2,v1,v2,omega,t){
  
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


#beta function to see the total trait change from pairwise interaction


be.ana<-function(t1,t2,v1,v2,omega,t){
    
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


#Pop change in a generation for a network
pops<-function(N,K,r,tmean,tvar,omega,t,a0,a1){
  alive<-which(N>0)
  pop1<-vector(length=length(N))
  for(i in alive){
    int<-0
    for(j in alive){
      if(i==j){a<-0
      #a0*al.ana(tmean[i],tmean[j],tvar[i],tvar[j],omega,t)*Ndash[j]
      } else{a<-a1*al.ana(tmean[i],tmean[j],tvar[i],tvar[j],omega,t)*N[j]}
      int<-int+a
    }
    r1<-0.75-(pnorm(-5,tmean[i],0.5)+1-pnorm(5,tmean[i],0.5))
    pop1[i]<-N[i]+(N[i]*r[i]*r1*(1-((N[i]+int)/K[i])))
  }
  pop1[which(pop1<=0)]<-0
  return(pop1)
}

a<-seq(-5,5,length.out=100)
b<-0.75-(pnorm(-5,tmean[i],0.5)+1-pnorm(5,tmean[i],0.5))

#Trait change in a generation for a network

trait<-function(N,K,h,r,tmean,tvar,omega,t,a0,a1){
  alive<-which(N>0)
  dead<-which(N<=0)
  traits<-vector(length=length(N))
  for(i in alive){
    int<-0
    Ndash<-N/K[i]
    for(j in alive){
      if(i==j){a<-0
      #a0*be.ana(tmean[i],tmean[j],tvar[i],tvar[j],omega,t)*Ndash[j]
      } else{a<-a1*be.ana(tmean[i],tmean[j],tvar[i],tvar[j],omega,t)*Ndash[j]}
      int<-int+a
    }
    r1<-1-(pnorm(-5,tmean[i],tvar[i])+1-pnorm(5,tmean[i],tvar[i]))
    traits[i]<-tmean[i] - ((h[i]^2)*r1*r[i]*int)
  }
  traits[dead]<-NA
  return(traits)
}


#No of species +param values
nsp<-20
h<-rep(0.5,nsp)                             
t<-2                           #Threshold to competition kernel
a0<-1
a1<-0.5                           #Strength of interspecific competition relative to intraspecific comp.
omega<-1
tvar<-rep(1,nsp) 

reps<-10

tmean<-means[,1]
N<-rep(1000,nsp)
N<-runif(nsp,500,1500)
K<-2*N
r<-runif(nsp,0.02,0.05)
r<-rep(1/30,nsp)



samples<-ceiling(c(seq(1,100,length.out=10),seq(500,9999,length.out=20)))


#Start simulation
        Nnew<-N
        tmeannew<-tmean
        Ndat<-matrix(nrow=nsp,ncol=length(samples))
        trdat<-matrix(nrow=nsp,ncol=length(samples))
        for(m in 1:10000){
          tmeannew1<-trait(Nnew,K,h,r,tmeannew,tvar,omega,t,a0,a1)
          Nnew<-pops(Nnew,K,r,tmeannew1,tvar,omega,t,a0,a1)
          Nnew[Nnew<1]<-0
          if(m%in%samples){
            Ndat[,which(samples==m)]<-Nnew
            trdat[,which(samples==m)]<-tmeannew1
          }
          tmeannew<-tmeannew1
        }
        

###################
# Plot results
        
#Plot trait evolution trajectories
plot(trdat[1,]~samples,type="l",ylim=c(-5,5),ylab="Trait means",xlab="Time",cex.lab=1.5,mgp=c(1.8,0.5,0))
for(i in 2:nsp){
  lines(trdat[i,]~samples)
}        

plot(Ndat[1,]~samples,type="l",ylim=c(0,3000))        
for(i in 2:nsp){
  lines(Ndat[i,]~samples)
}        

means<-trdat
pops<-Ndat
#plot the trait trajectory with the line widths indicating the population sizes
par(mgp=c(axis.title.position, axis.label.position, axis.line.position))


tiff('GT_sample2.tiff', units="in", width=10, height=7, res=300)

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

modat<-c() 
for(i in 1:length(samples)){
  modat<-c(modat,mod(trdat[,i],tvar,omega,t)[1])
}
plot(modat~samples)

modat1<-c() 
for(i in 1:length(samples)){
  modat1<-c(modat1,mod(trdat[,i],tvar,omega,t)[2])
}
plot(modat1~samples)

