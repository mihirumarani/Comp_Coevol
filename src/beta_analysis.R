#Functions
library(truncnorm)
library(tidyverse)
library(patchwork)
library(plotly)

#Regression trees
library(rsample)     # data splitting 
library(rpart)       # performing regression trees
library(rpart.plot)  # plotting regression trees
library(ipred)       # bagging
library(caret)       # bagging


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

alpha.gt=function(x,y,omega,t){
  a1=0
  if(abs(x-y)<t){
    a1=exp(-((x-y)^2)/(omega^2))
  }
  return(a1)
}

mnnd<-function(a){
  a<-na.omit(a)
  if(length(a)>1 & length(unique(a)>1)){
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

#################################
#Analyse 3 species convergence test
####################################

#Plot hypothetical extreme cases

dat=read.csv("D:/project_files/comp_coevol/data/3sp_test7.csv")%>%as_tibble()
dat_gau=read.csv("D:/project_files/comp_coevol/data/3sp_test_gaus.csv")%>%
  as_tibble()%>%
  filter((t>(-1-sp1)) & (t>(sp3-1)))
  

dat_uni=dat%>%
  filter(dist=="Uniform")%>%
  filter((t>(-0.8-sp1)) & (t>(sp3-0.8)))

dat=dat_gau%>%
  bind_rows(dat_uni)%>%
  filter(sp1 !=0 & sp3!=0)


m3=rpart(formula=diff1~.,
         data=dat%>%select(-diff2),
         method='anova',
         xval=10,
         )

dat5=dat%>%filter(loci==5)
dat7=dat%>%filter(loci==7)
dat10=dat%>%filter(loci==10)
dat20=dat%>%filter(loci==20)
dat3=dat%>%filter(loci==3)

t.labs=str_wrap(c("Kernel Threshold: 0.1","Kernel Threshold: 0.3",
         "Kernel Threshold: 0.5","Kernel Threshold: 0.7"),width=18)
names(t.labs)=c("0.1","0.3","0.5","0.7")
omega.labs=str_wrap(c("Kernel width: 0.05","Kernel width: 0.25",
             "Kernel width: 0.45","Kernel width: 0.65",
             "Kernel width: 0.85"),width=10)
names(omega.labs)=c("0.05","0.25","0.45","0.65","0.85")


p20g=dat20%>%
      filter(dist=="Gaussian")%>%
      mutate(diff12=sign(diff1))%>%
      ggplot(aes(sp3,sp1,fill=diff12))+
      geom_tile()+
      facet_grid(omega~t,labeller=labeller(t=t.labs,omega=omega.labs))+
      scale_fill_gradient(low="red",high="white")+
      ggtitle("Gaussian traits")+
      labs(x="Mean trait of species 3",
           y="Mean trait of species 1")+
      theme(legend.position="none")
      

p20u=dat20%>%
  filter(dist=="Uniform")%>%
  mutate(diff12=sign(diff1))%>%
  ggplot(aes(sp3,sp1,fill=diff12))+
  geom_tile()+
  facet_grid(omega~t,labeller=labeller(t=t.labs,omega=omega.labs))+
  scale_fill_gradient(low="red",high="white")+
  ggtitle("Uniform traits")+
  labs(x="Mean trait of species 3",
       y="Mean trait of species 1")+
  theme(legend.position="none")

p20g|p20u

p10g=dat10%>%
  filter(dist=="Gaussian")%>%
  mutate(diff12=sign(diff1))%>%
  ggplot(aes(sp3,sp1,fill=diff12))+
  geom_tile()+
  facet_grid(omega~t,labeller=labeller(t=t.labs,omega=omega.labs))+
  scale_fill_gradient(low="red",high="white")+
  ggtitle("Gaussian traits")+
  labs(x="Mean trait of species 3",
       y="Mean trait of species 1")+
  theme(legend.position="none")


p10u=dat10%>%
  filter(dist=="Uniform")%>%
  mutate(diff12=sign(diff1))%>% 
  ggplot(aes(sp3,sp1,fill=diff12))+
  geom_tile()+
  facet_grid(omega~t,labeller=labeller(t=t.labs,omega=omega.labs))+
  scale_fill_gradient(low="red",high="white")+
  ggtitle("Uniform traits")+
  labs(x="Mean trait of species 3",
       y="Mean trait of species 1")+
  theme(legend.position="none")

p10g|p10u

"Caption: Result of single instance competition between three species where the
mean trait of species 2 was fixed at zero while means of species 1 and 3 were allowed
to vary ((-1,0) for species 1 and (0,1) for species 3).
Red colour indicates that mean traits species 1 and 2 converged while the white color
shows the cases of divergence.
Two panels show whether the traits of each species followed Gaussian or Uniform distribitions.
Horizontal (competition trait threshold) and vertical (kernel width) facets show
different values of parameters of the competition kernel.
"


        
#dat2=read.csv(paste0(loc,"3sp_test5.csv"))%>%as_tibble()

#dat=dat1%>%bind_rows(dat2)
# Columns 'sp1', 'sp2' and 'sp3' show the initial trait values of 3 species.
# Column "diff1 and diff2" show the total divergence between two species pairs (1-2 and 2-3) at 
#the end of the simulation. (traits of species are ordered to begin with
#sp1 has the lowest trait value)

#Do the convergence events differ systematically with model parameters?

p11=dat%>%
    mutate(tr=-sp1/sp3)%>%
    filter(dist=="Gaussian")%>%
    mutate(loci=as.factor(loci))%>%
    ggplot(aes(tr,diff1,col=loci))+
    geom_line()+
    facet_grid(omega~t,labeller=label_both)+
    geom_hline(yintercept=0)+
    ggtitle("Gaussian traits")+
    labs(x="Ratio of gaps between sp. pairs",
         y="Trait divergence between sp 1 and 2")

p12=dat%>%
  mutate(tr=-sp1/sp3)%>%
  filter(dist=="Gaussian")%>%
  mutate(loci=as.factor(loci))%>%
  ggplot(aes(tr,diff1,col=loci))+
  geom_line()+
  facet_grid(omega~t,labeller=label_both)+
  geom_hline(yintercept=0)+
  ggtitle("Uniform traits")+
  labs(x="Ratio of gaps between sp. pairs",
       y="Trait divergence between sp 1 and 2")

p11|p12


#Figure captions:
#  Trait divergence between pairs of species vs. the ratio of gaps between two species pairs.
#   (sp2-sp1)/(sp3-sp2). Lower value of x-axis indicates that the spp. 1 and 2 are much closer
#   while sp. 3 is farther on trait axis.Higher x-axis value indicates vice-versa. Y axis values
#   below zero indicate convergence.



#Draw a heatmap
h11=dat%>%
  drop_na()%>%
  group_by(loci,dist,omega,t)%>%
  summarize(conv1=sum(diff1<0)/length(diff1),
            conv2=sum(diff2)/length(diff2))%>%
  ungroup()%>%
  ggplot(aes(omega,t,fill=conv1))+
  geom_tile()+
  facet_grid(loci~dist)+
  scale_fill_gradientn(colours=c("white","red"),
                       values=c(0,1))

h12=dat%>%
  drop_na()%>%
  group_by(loci,dist,omega,t)%>%
  summarize(conv1=sum(diff1<0)/length(diff1),
            conv2=sum(diff2)/length(diff2))%>%
  ungroup()%>%
  ggplot(aes(omega,t,fill=conv2))+
  geom_tile()+
  facet_grid(loci~dist)+
  scale_fill_gradientn(colours=c("white","red"),
                       values=c(0,1))



#Do convergent events depend on the how the species are located on trait axis?

sp1dat1=dat%>%
  group_by(sp1,loci,dist,omega,t)%>%
  summarize(conv1=sum(diff1<0)/length(diff1))%>%
  ungroup()%>%
  filter(dist=="Gaussian")%>%
  mutate(loci=as.factor(loci))%>%
  ggplot(aes(sp1,conv1,col=loci))+
  geom_line()+
  facet_grid(omega~t,labeller=label_both)+
  geom_hline(yintercept=0)+
  ggtitle("Gaussian traits")+
  labs(x="Species 1 trait",
       y="Convergence probability for sp1 and sp2")

sp1dat2=dat%>%
  group_by(sp1,loci,dist,omega,t)%>%
  summarize(conv1=sum(diff1<0)/length(diff1))%>%
  ungroup()%>%
  filter(dist=="Uniform")%>%
  mutate(loci=as.factor(loci))%>%
  ggplot(aes(sp1,conv1,col=loci))+
  geom_line()+
  facet_grid(omega~t,labeller=label_both)+
  geom_hline(yintercept=0)+
  ggtitle("Uniform traits")+
  labs(x="Species 1 trait",
       y="Convergence probability for sp1 and sp2")

sp1dat1|sp1dat2


sp3dat1=dat%>%
  group_by(sp3,loci,dist,omega,t)%>%
  summarize(conv1=sum(diff1<0)/length(diff1))%>%
  ungroup()%>%
  filter(dist=="Gaussian")%>%
  mutate(loci=as.factor(loci))%>%
  ggplot(aes(sp3,conv1,col=loci))+
  geom_line()+
  facet_grid(omega~t,labeller=label_both)+
  geom_hline(yintercept=0)+
  ggtitle("Gaussian traits")+
  labs(x="Species 3 trait",
       y="Convergence probability for sp1 and sp2")

sp3dat2=dat%>%
  group_by(sp3,loci,dist,omega,t)%>%
  summarize(conv1=sum(diff1<0)/length(diff1))%>%
  ungroup()%>%
  filter(dist=="Uniform")%>%
  mutate(loci=as.factor(loci))%>%
  ggplot(aes(sp3,conv1,col=loci))+
  geom_line()+
  facet_grid(omega~t,labeller=label_both)+
  geom_hline(yintercept=0)+
  ggtitle("Uniform traits")+
  labs(x="Species 3 trait",
       y="Convergence probability for sp1 and sp2")

sp3dat1|sp3dat2


spgap1=dat%>%
  group_by(sp1,sp3,dist)%>%
  summarize(conv1=sum(diff1<0)/length(diff1),
            conv2=sum(diff2<0)/length(diff2))%>%
  ungroup()%>%
  mutate(sp1=-sp1)%>%
  mutate(sp1=as.factor(sp1),
         sp3=as.factor(sp3))%>%
  ggplot(aes(sp1,sp3,fill=conv1))+
  geom_tile()+
  facet_wrap(vars(dist))+
  scale_fill_gradientn(colours=c("white","red"),
                       values=c(0,1))+
  labs(x="species 1 trait",
       y="species 3 trait")


spgap2=dat%>%
  group_by(sp1,sp3,dist)%>%
  summarize(conv1=sum(diff1<0)/length(diff1),
            conv2=sum(diff2<0)/length(diff2))%>%
  ungroup()%>%
  mutate(sp1=-sp1)%>%
  mutate(sp1=as.factor(sp1),
         sp3=as.factor(sp3))%>%
  ggplot(aes(sp1,sp3,fill=conv2))+
  geom_tile()+
  facet_wrap(vars(dist))+
  scale_fill_gradientn(colours=c("white","red"),
                       values=c(0,1))

spgap1|spgap2

###################################################################
#Plot 20 sp simulation results
#####################################################################

#################################################
#Find proportions of pairs of species converging 
##################################################


convdat=read.csv("D:/project_files/comp_coevol/data/convdatfull.csv")%>%
  as_tibble()

convdat%>%
  filter(t==1 & omega==0.5 & dist=="Gaussian")%>%
  mutate(a=as.factor(a))%>%
  group_by(nsp,loci,dist,a,t,time)%>%
  summarize(mean=mean(conv),
            sd=sd(conv))%>%
  ungroup()%>%
  ggplot(aes(time,mean,col=a,fill=a))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd),alpha=0.2)+
  facet_grid(nsp~loci,labeller=label_both)

times=unique(convdat$time)
nsps=c(5,10,20)
pars=expand.grid(nsp=nsps,time=times)
result2=tibble()

varset=tibble(var=c("loci","dist","a","omega","t","rep"))

for (i in 1:nrow(pars)){
  
  dat=convdat%>%
    filter(nsp==pars[i,1],time==pars[i,2])
  
  m=rpart(
    formula=conv~.,
    data=dat,
    method='anova',
    xval=10
  )
  
  res=m$variable.importance

  vars=names(res)
  
  rs=tibble(var=vars,imp=res)
  
  if(nrow(rs)>0){
    rs=rs%>%full_join(varset)
  }else{
    rs=tibble(var=varset$var,imp=0)
  }

  erdat=m$cptable
  
  error=erdat[nrow(erdat),3]
  
  result2=result2%>%
    bind_rows(
      tibble(
        nsp=pars[i,1],
        time=pars[i,2],
        vars=rs$var,
        imp=rs$imp,
        error=error
      )
    )
  
}

result2%>%
  mutate(nsp=as.factor(nsp),
          rsq=1-error)%>%
  ggplot(aes(time,rsq,col=nsp))+
  geom_line()


result2%>%
  replace_na(list(imp=0))%>%
  ggplot(aes(time,imp,col=vars))+
  geom_line()+
  facet_wrap(vars(nsp))



  convdat%>%
  group_by(nsp,loci,dist,a,omega,t)%>%
  summarize(mean=mean(conv),
            sd=sd(conv))%>%
  ungroup()%>%
  rename("species"=nsp)%>%
  mutate(omega=as.factor(omega))%>%
  filter(loci==5 & dist=="Gaussian")%>%
  ggplot(aes(t,mean,col=omega,fill=omega))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd),alpha=0.2)+
  facet_grid(a~species,labeller=label_both)+
  labs(x="Competition trait threshold",
       y="Average proportion of convergences")
  
  #Draw a decision tree
  
  m1= rpart(
    formula = conv ~ .,
    data    = convdat,
    method  = "anova"
  )
  
  rpart.plot(m1)
  
  plotcp(m1)
  
  convdat10=read.csv("D:/project_files/comp_coevol/data/convdat10.csv")%>%
    as_tibble()
  
  m10=rpart(
    formula = conv ~ .,
    data    = convdat10,
    method  = "anova"
  )
  
  rpart.plot(m10)
  
  plotcp(m10)

############################################################
#MNND values

mnndat=read.csv("D:/project_files/comp_coevol/data/mnndat.csv")%>%
    as_tibble()
  
mndat2=mnndat%>%
  filter(time==max(time))%>%
  select(-time)

mnndat=mndat2%>%
    group_by(nsp,loci,dist,a,omega,t)%>%
    summarize(mean=mean(mnnds),
            sd=sd(mnnds))%>%
    ungroup()

mnndat5=mnndat%>%filter(nsp==5)
mnndat10=mnndat%>%filter(nsp==10)
mnndat20=mnndat%>%filter(nsp==20)

mn5g=mnndat10%>%
  filter(dist=="Gaussian")%>%
  mutate(loci=as.factor(loci))%>%
  ggplot(aes(a,mean,col=loci,fill=loci))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd),alpha=0.2)+
  facet_grid(omega~t,labeller=label_both)+
  labs(x="Relative strength of interspecific competition",
       y="MNND")+
  ggtitle("Gaussian traits")
  
mn5u=mnndat10%>%
    filter(dist=="Uniform")%>%
    mutate(loci=as.factor(loci))%>%
    ggplot(aes(a,mean,col=loci,fill=loci))+
    geom_line()+
    geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd),alpha=0.2)+
    facet_grid(omega~t,labeller=label_both)+
    labs(x="Relative strength of interspecific competition",
         y="MNND")+
  ggtitle("Uniform traits")
  
  mn5g|mn5u
  
  
m2=rpart(
  formula = mnnds ~ .,
  data    = mndat2,
  method  = "anova"
)

rpart.plot(m2)

plotcp(m2)
  
mndat=mnndat%>%
  filter(time==max(time))%>%
  group_by(nsp,loci,dist,a,omega,t)%>%
  summarize(mean=mean(mnnds),
            sd=sd(mnnds))%>%
  ungroup()

mndat%>%
  ggplot(aes(omega,t,fill=mean))+
  geom_tile()+
  facet_grid(loci~a)+
  scale_fill_gradientn(colours=c("white","red"),
                       values=c(0,1))
  

mplot1=mnndat%>%
  rename("species"=nsp)%>%
  mutate(omega=as.factor(omega))%>%
  filter(loci==10 & dist=="Uniform" & time==1000)%>%
  ggplot(aes(t,mean,col=omega,fill=omega))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd),alpha=0.2)+
  facet_grid(a~species,labeller=label_both)+
  labs(x="Trait threshold to competition",
       y="Average MNND value")+
  ggtitle("Trait dispersion pattern")

mplot2=mnndat%>%
  rename("species"=nsp)%>%
  mutate(loci=as.factor(loci))%>%
  filter(omega==.5 & t==0.5 & time==1000)%>%
  ggplot(aes(a,mean,col=loci,fill=loci))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd),alpha=0.2)+
  facet_grid(dist~species,labeller=label_both)+
  labs(x="Strength of interspecific competition",
       y="Average MNND value")+
  ggtitle("Trait dispersion pattern")



############################################################
#Extinctions
extdat=read.csv("D:/project_files/comp_coevol/data/extdat.csv")%>%
  as_tibble()

ext.fin=extdat%>%
  filter(time == max(time))%>%
  select(-time)

ex.m=rpart(
    formula = extinct ~ .,
    data    = ext.fin,
    method  = "anova"
  )
  
rpart.plot(ex.m)
  
plotcp(ex.m)

  


eplot2=extdat%>%
  mutate(loci=as.factor(loci))%>%
  filter(omega==.5 & t==0.5 & time==1000)%>%
  group_by(nsp,loci,dist,a)%>%
  summarize(mean=mean(live),
            sd=sd(live))%>%
  ungroup()%>%
  rename("species"=nsp)%>%
  ggplot(aes(a,mean,col=loci,fill=loci))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd),alpha=0.2)+
  facet_grid(dist~species,labeller=label_both)+
  labs(x="Strength of interspecific competition",
       y="Average No. of extant species")+
  ggtitle("extinction pattern")

eplot1=extdat%>%
  mutate(omega=as.factor(omega))%>%
  filter(loci==10 & dist=="Uniform")%>%
  group_by(nsp,loci,dist,omega,t,a)%>%
  summarize(mean=mean(live),
            sd=sd(live))%>%
  ungroup()%>%
  rename("species"=nsp)%>%
  ggplot(aes(t,mean,col=omega,fill=omega))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd),alpha=0.2)+
  facet_grid(a~species,labeller=label_both)+
  labs(x="Trait threshold to competition",
       y="Average No. of extant species")+
  ggtitle("extinction pattern")

mplot1|eplot1

mplot2|eplot2

################################################################################
#There are too many parameters to get clear signals from plots
################################################################################

#Regression trees
library(rsample)     # data splitting 
library(rpart)       # performing regression trees
library(rpart.plot)  # plotting regression trees
library(ipred)       # bagging
library(caret)       # bagging

convdat=read.csv("D:/project_files/comp_coevol/data/convdat.csv")%>%
  as_tibble()

convdat100=read.csv("D:/project_files/comp_coevol/data/convdat100.csv")%>%
  as_tibble()

convdat100%>%
  group_by(nsp,loci,dist,a,omega,t)%>%
  summarize(mean=mean(conv),
            sd=sd(conv))%>%
  ungroup()%>%
  rename("species"=nsp)%>%
  mutate(omega=as.factor(omega))%>%
  filter(loci==5 & dist=="Gaussian")%>%
  ggplot(aes(a,mean,col=omega,fill=omega))+
  geom_line()+
  geom_ribbon(aes(ymin=mean-sd,ymax=mean+sd),alpha=0.2)+
  facet_grid(t~species,labeller=label_both)+
  labs(x="Strength of interspecific interaction",
       y="Average proportion of convergences")


convdat=convdat%>%
  mutate(
         dist=as.factor(dist),
        )




m1 <- rpart(
  formula = conv ~ .,
  data    = convdat,
  method  = "anova"
)

rpart.plot(m1)

plotcp(m1)


m100=rpart(
  formula=conv~.,
  data=convdat100,
  method='anova'
)

rpart.plot(m100)



convdat2=convdat%>%
  group_by(nsp,loci,dist,omega,t,a)%>%
  summarize(conv=mean(conv))%>%
  ungroup()%>%
  filter(nsp==5)%>%
  select(-nsp)

models=list()

pars=expand(convdat,nsp,rep,a)

for (i in 1:nrow(pars)){
  
  dat1=convdat%>%
    filter(rep==pars$rep[i] & nsp== pars$nsp[i] & a==pars$a[i])%>%
    select(dist,loci,omega,t,conv)
  
  models[[i]]=rpart(
    formula=conv~.,
    data=dat1,
    method="anova"
  )
}

get_n_l=function(x){
  err=min(x$cptable[,5])
  fr=x$frame
  nodes=fr[which(fr$var != "<leaf>"),1]
  grps=fr[which(fr$var == "<leaf>"),c(4,5)]
  
  return(list(nodes=nodes,grps=grps))
}

get_freq=function(x){
  res=tibble(pars=c("dist","loci","omega","t"),
             times=0)
  res$times[1]=sum(x == "dist")
  res$times[2]=sum(x == "loci")
  res$times[3]=sum(x == "omega")
  res$times[4]=sum(x == "t")
  
  return(res)
}

res.nodes=tibble()
res.grps=tibble()

for(i in 1:nrow(pars)){
  
  res=get_n_l(models[[i]])
  
  res.nodes=res.nodes%>%
    bind_rows(tibble(
      nsp=pars$nsp[i],
      a=pars$a[i],
      rep=pars$rep[i],
      nodes=res$nodes
    ))
  
  res.grps=res.grps%>%
    bind_rows(tibble(
      nsp=pars$nsp[i],
      a=pars$a[i],
      rep=pars$rep[i],
      means=res$grps[,2],
      devs=res$grps[,1]
    ))

}

res.nodes=res.nodes%>%
  group_by(nsp,a,rep)%>%
  reframe(get_freq(nodes))

res.nodes%>%
  mutate(rep=as.factor(rep))%>%
  ggplot(aes(pars,times,fill=rep))+
  geom_bar(stat="identity")+
  facet_grid(nsp~a,labeller=label_both)

res.nodes%>%
  group_by(nsp,a,pars)%>%
  summarize(mean=mean(times),
            sd=sd(times))%>%
  ungroup()%>%
  mutate(a=as.factor(a))%>%
  ggplot(aes(pars,mean))+
  geom_point()+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd))+
  facet_grid(nsp~a,labeller=label_both)


m1 <- rpart(
  formula = conv ~ .,
  data    = conv_train,
  method  = "anova"
)

rpart.plot(m1)

plotcp(m1)
