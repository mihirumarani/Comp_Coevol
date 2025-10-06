#Functions
library(truncnorm)
library(tidyverse)
library(patchwork)
library(plotly)


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

dat=read.csv("D:/project_files/comp_coevol/data/3sp_test6.csv")%>%as_tibble()
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


#Convergence likelihoods are similar for Gaussian vs. uniformly distributed traits!

###################################################################
#Plot 20 sp simulation results
#####################################################################

loc5="D:/project_files/comp_coevol/data/5sp/"
loc10="D:/project_files/comp_coevol/data/10sp/"
loc20="D:/project_files/comp_coevol/data/20sp/"

dat5=tibble()
files5=list.files(loc5)
for(i in 1:length(files5)){
  dat5=dat5%>%
        bind_rows(
        (read.csv(paste0(loc5,files5[i]))%>%
           as_tibble()))
}
dat5=dat5%>%mutate(nsp=5)

dat10=tibble()
files10=list.files(loc10)
for(i in 1:length(files10)){
  dat10=dat10%>%
    bind_rows(nsps
      (read.csv(paste0(loc10,files10[i]))%>%
         as_tibble()))
}
dat10=dat10%>%mutate(nsp=10)

dat20=tibble()
files20=list.files(loc20)
for(i in 1:length(files20)){
  dat20=dat20%>%
    bind_rows(
      (read.csv(paste0(loc20,files20[i]))%>%
         as_tibble()))
}
dat20=dat20%>%mutate(nsp=20)


dat=dat5%>%bind_rows(dat10)%>%bind_rows(dat20)

#################################################
#Find proportions of pairs of species converging 
##################################################

##########
##USE THIS CODE ONCE TO CREATE THE CONVDAT FILE
##########

#Code to calculate no. of convergences between pairs of adjacent species.

nsps=c(5,10,20)

dat=dat%>%
  mutate(trmean1=ifelse(pop>1,trmean,NA))

convdat=tibble()

for (j in nsps){
  
  dat0=dat%>%filter(nsp ==j)
  pars=dat0%>%count(loci,dist,a,omega,t,rep)%>%as.data.frame()
  
  for (i in 1:nrow(pars)){
    
    dat1=dat0%>%
      filter(loci==pars$loci[i],
             dist==pars$dist[i],
             a==pars$a[i],
             omega==pars$omega[i],
             t==pars$t[i],
             rep==pars$rep[i])
  
    ord1=order(dat1%>%
                 filter(time==0)%>%
                 pull(trmean))
  
    cps=dat1%>%
        #filter(time %in% c(0,max(time)))%>%
        group_by(time)%>%
        slice(ord1)%>%
        ungroup()%>%
      group_by(time)%>%
      reframe(td=diff(trmean))%>%
      ungroup()%>%
      group_by(time)%>%
      mutate(sppair=1:length(td))%>%
      ungroup()%>%
      group_by(sppair)%>%
      reframe(conv=diff(td))%>%
      #summarize(conv=sum(diff(td)<0))%>%
      ungroup()%>%
      summarize(sum(conv<0)/length(sppair))%>%pull()
    
    convdat=convdat%>%
      bind_rows(tibble(
        nsp=j,
        loci=pars$loci[i],
        dist=pars$dist[i],
        a=pars$a[i],
        omega=pars$omega[i],
        t=pars$t[i],
        rep=pars$rep[i],
        conv=cps))
    }
}

write.csv(convdat,"D:/project_files/comp_coevol/data/convdat.csv")

################################################################################


convdat=read.csv("D:/project_files/comp_coevol/data/convdat.csv")%>%
  as_tibble()

  convdat%>%
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
  

############################################################
#MNND values

mnndat=dat%>%
        group_by(nsp,loci,dist,a,omega,t,rep,time)%>%
        summarize(mnnds=mnnd(trmean))%>%
        ungroup()%>%
        group_by(nsp,loci,dist,a,omega,t,time)%>%
        summarize(mean=mean(mnnds),
                  sd=sd(mnnds))%>%
        ungroup()

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
extdat=dat%>%
  group_by(nsp,rep,loci,dist,omega,t,a,time)%>%
  summarize(live=sum(pop>1)/length(pop))%>%
  ungroup()%>%
  filter(time==1000)



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


convdat=convdat%>%
  mutate(
         dist=as.factor(dist),
        )

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
