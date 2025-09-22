

#Determine alpha values

erfun<-function(x){
  return(2 * pnorm(x * sqrt(2)) - 1)
}



al.ana<-function(t1,t2,v1,v2,omega,t){
#  if(any(is.na(c(t1,t2)))){
 #  return(0)
  #} else{
  
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
 # }
}


#beta function to see the total trait change from pairwise interaction


be.ana<-function(t1,t2,v1,v2,omega,t){
  #if(any(is.na(c(t1,t2)))){
   #  return(0)
  #  } else{
  
  
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
#  }
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
    r1<-r[i]*(1-(max(dnorm(5.5,tmean[i],tvar[i]),dnorm(-5.5,tmean[i],tvar[i]))/dnorm(tmean[i],tmean[i],tvar[i])))
    pop1[i]<-N[i]+(N[i]*r1*(1-((N[i]+int)/K[i])))
  }
  pop1[which(pop1<=0)]<-0
  return(pop1)
}




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
    r1<-r[i]*(1-(max(dnorm(5.5,tmean[i],tvar[i]),dnorm(-5.5,tmean[i],tvar[i]))/dnorm(tmean[i],tmean[i],tvar[i])))
    traits[i]<-tmean[i] - ((h[i]^2)*r1*int)
  }
  traits[dead]<-NA
  return(traits)
}

  #No of species +param values
nsp<-10
h<-rep(0.5,nsp)
omegas<-c(0.5,1,2)                             #Width of competition kernel
ts<-c(1,2.5,5,100)                           #Threshold to competition kernel
a0<-1
a1s<-c(0.2,0.5,0.9)                            #Strength of interspecific competition relative to intraspecific comp.
tvar<-runif(nsp,0.5,2.5) 

reps<-1

#############################################
#Start the cluster and split the loop on 8 cores
library(parallel)
library(foreach)
library(iterators)
library(doParallel)
library(truncnorm)

#Create a set of 100 initializations for tmean, N, r and K
tmeanset<-matrix(rtruncnorm(reps*nsp,-5,5,mean=0,sd=10),ncol=nsp)  
#Nset<-matrix(rep(abs(rnorm(50,1000,200)),n),ncol=n)
Nset<-matrix(runif(reps*nsp,500,1500),nrow=reps,ncol=nsp)
#rset<-matrix(rep(abs(rnorm(50,0,0.5)),n),ncol=n)
rset<-matrix(abs(rnorm(reps*nsp,0,0.1)),ncol=nsp)
Kset<-2*Nset


modsamples<-ceiling(c(seq(1,2000,length.out=10),seq(2100,29999,length.out=20)))

numcores<-detectCores()
clust<-makeCluster(min(numcores,24))
clusterExport(clust,c("reps","trait","pops","be.ana","al.ana","erfun","nsp","h","omegas","tvar",
                      "ts","modsamples","tmeanset","Nset","rset","Kset","a0","a1s"))
registerDoParallel(clust)


dat<-foreach(i=1:reps)%dopar%
{ 
   tmean<-tmeanset[i,]
   N<-Nset[i,]
   r<-rset[i,]
   K<-Kset[i,]
      
   Nnew<-N
   tmeannew<-tmean
 
  paramdat<-list()
  trdat<-list()
  popdat<-list()
  
  for(j in 1:length(omegas)){
    omega<-omegas[j]
    for(k in 1:length(ts)){
      t<-ts[k]
      for(l in 1:length(a1s)){
          a1<-a1s[l]
          
  

          popdat1<-c()
          trdat1<-c()
            for(m in 1:30000){
              tmeannew1<-trait(Nnew,K,h,r,tmeannew,tvar,omega,t,a0,a1)
              Nnew<-pops(Nnew,K,r,tmeannew1,tvar,omega,t,a0,a1)
              Nnew[Nnew<1]<-0
              if(m%in%modsamples){
                
                popdat1<-c(popdat1,Nnew)
                trdat1<-c(trdat1,tmeannew)
              }
            tmeannew<-tmeannew1
          }
          
          paramdat<-c(paramdat,list(c(t,omega,a1)))
          popdat<-c(popdat,list(popdat1))
          trdat<-c(trdat,list(trdat1))
        }
  }
  }
  mlode<-list(p=paramdat,pop=popdat,tr=trdat)
  return(mlode)
}

stopCluster(clust)
dat<-c(dat,list(tvar))
saveRDS(dat,file="QG20_55.Rds")



###############################




