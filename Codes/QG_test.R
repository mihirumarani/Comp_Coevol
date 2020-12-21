library(truncnorm)

#Determine alpha values

erfun<-function(x){
  return(2 * pnorm(x * sqrt(2)) - 1)
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

pops<-function(N,r,tmean,tvar,omega,t,a0,a1){
  pop1<-c()
  N[which(N<=0)]<-0
  for(i in 1:length(N)){
    int<-0
    for(j in 1:length(N)){
      Ndash<-N/sum(N)
      if(i==j){a<-a0*al.ana(tmean[i],tmean[j],tvar[i],tvar[j],omega,t)*Ndash[j]
      } else{a<-a1*al.ana(tmean[i],tmean[j],tvar[i],tvar[j],omega,t)*Ndash[j]}
      int<-int+a
    }
    pop1[i]<-N[i]+N[i]*(r[i]*(1-int))
  }
  pop1[which(pop1<=0)]<-0
  return(pop1)
}



#Trait change in a generation for a network

trait<-function(N,h,r,tmean,tvar,omega,t,a0,a1){
  tmean[which(N<=0)]<-NA
  N[which(N<=0)]<-0
  traits<-c()
  for(i in 1:length(N)){
    int<-0
    for(j in 1:length(N)){
      Ndash<-N/sum(N)
      if(i==j){a<-a0*be.ana(tmean[i],tmean[j],tvar[i],tvar[j],omega,t)*Ndash[j]
      } else{a<-a1*be.ana(tmean[i],tmean[j],tvar[i],tvar[j],omega,t)*Ndash[j]}
      int<-int+a
    }
    traits[i]<-tmean[i] - ((h[i]^2)*r[i]*int)
  }
  return(traits)
}

#No of species +param values
nsp<-20
h<-rep(0.5,nsp)                             
t<-2                           #Threshold to competition kernel
a0<-1
a1<-0.5                           #Strength of interspecific competition relative to intraspecific comp.
omega<-0.5
tvar<-runif(nsp,0.5,2.5) 

reps<-10

tmean<-rtruncnorm(nsp,a=-5,b=5,0,10)
N<-runif(nsp,500,1500)
r<-abs(rnorm(nsp,0,0.01))



samples<-ceiling(c(seq(1,1000,length.out=10),seq(1500,29999,length.out=10)))


#Start simulation
        Nnew<-N
        tmeannew<-tmean
        Ndat<-N
        trdat<-tmean
        for(m in 1:30000){
          tmeannew1<-trait(Nnew,h,r,tmeannew,tvar,omega,t,a0,a1)
          Nnew<-pops(Nnew,r,tmeannew1,tvar,omega,t,a0,a1)
          if(m%in%samples){
            Ndat<-cbind(Ndat,Nnew)
            trdat<-cbind(trdat,tmeannew1)
          }
          tmeannew<-tmeannew1
        }
        
trdat<-trdat[,-1]
Ndat<-Ndat[,-1]
###################
# Plot results
        
#Plot trait evolution trajectories
plot(trdat[1,]~samples,type="l",ylim=c(-10,10))
for(i in 2:nsp){
  lines(trdat[i,]~samples,col=i)
}        
        
        