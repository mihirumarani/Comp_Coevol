using CSV, DataFrames, Random, LinearAlgebra, Distances, Distributions, SpecialFunctions, Plots

#Necessary functions

function alpha_gt(x::Float64, y::Float64, omega::Float64, t::Float64)
    
    al=0
    if(abs(x-y)<t) 
         al=exp(-((x-y)^2)/(omega^2))
    end
    return al
end

function alpha_tri(x::Float64, y::Float64, omega::Float64, t::Float64)
    
    alpha=0.0
    
    if ( abs(x-y) < t) 

        slope=0.5/omega
        
        alpha=1- (slope*abs(x-y))
        
    end
        
    return max(0,alpha)
    
end
    
function mnnds(a::Vector{Float64})

    a=filter!(!isnan,a)
    if(length(a)>1)
        a=sort(a,rev=true)
        nnd=zeros(Float64,length(a))
        nnd[1]=a[2]-a[1]
        nnd[length(nnd)]=a[length(a)]-a[length(a)-1]
        for i in 2:((length(a)-1))
            nnd[i]=min((a[i-1]-a[i]),(a[i]-a[i+1]))
        end
        ranges=a[1]-a[end]
        mnnd=sum(nnd)/length(nnd)
        mmax=ranges/(length(a)-1)

        return mnnd/mmax
    
    else return 0.0
    end
    
end
    
function qgprob(n::Int64)

    #All possible phenotypes
    pheno= collect(1:(2*n+1)) ./ (2*n+1)
    nt=length(pheno)

    G=zeros(Float64,n+1,n+1,n+1)

    for i in 0:n, j in 0:i, k in max(0,(i+j-n)):min(n,(i+j))
                m=collect(0:min(j,k,i+j-k))
                G[1+i,1+j,1+k]=sum(pdf.(Hypergeometric(i,n-i,j),m).*pdf.(Binomial.(i+j .- (2 .* m)),k .- m))
    end

    for k in 0:n
        G[:,:,1+k]=G[:,:,1+k]+transpose(G[:,:,1+k])
        for i1 in 0:n
            G[i1+1,i1+1,k+1] /= 2
        end
    end

    ind_haplR=zeros(Float64,2*n+1, 2*n+1)

    for k in 0:n
        for i in 0:n
             ind_haplR[1+i,1+k] = G[1+i,1,1+k]
            for j in 0:n
                ind_haplR[1+j+n,1+k]=G[1+n,1+j,1+k]
            end
        end
    end

    R=zeros(Float64,nt,nt,nt)

    for i in 0:(2*n), j in 0:(2*n), q in 0:(2*n)
         R[1+i,1+j,1+q]= sum(ind_haplR[1+i,1 .+ (0:q)] .*
                             ind_haplR[1+j,1+q .- (0:q)])
    end

    return R
end

function single_sim(r,K,a1,A,R,Ng0,Npop)
    
        Np0=Ng0 .*Npop
        Ngen=deepcopy(Ng0)
        Np=deepcopy(Np0)
    
        dat=zeros(Float64,1000,nsp,nt)
        
        #Start the simulation
        for m in 1:1000
    
            #Determine the extinct species
            Np[findall(sum(Np,dims=2) .< 10),:] .= 0
            Ngen[findall(sum(Np,dims=2) .==0),:] .= 0

            if all(sum(Np,dims=2) ==0) 
                break
            else

                newgen=zeros(Float64,nsp,nt)

                #Reproduction event
                for i in findall(!iszero,sum(eachcol(Ngen)))

                    probs=Ngen[i,:]*Ngen[i,:]'
                
                    newgen[i,:]=[sum(probs.*R[:,:,j]) for j in 1:nt]
                
                    newgen[i,:] ./= sum(newgen[i,:])
                    
                end

                newp=newgen .* sum(Np,dims=2)

                #Selection event
            
                if size(A)[1]>1

                    for i in 1:size(newp)[1]

                        comps=[sum((a1*A[j,:]) .* newp[1:end .!=i,:]') + sum(A[j,:] .* newp[i,:]) for j in 1:nt]
                    
                        Np[i,:] += newp[i,:] .* r[i] .* (1 .-(comps ./K[i]))

                    end
                
                end

                Np[findall(Np .<1)] .= 0
                Ngen= Np ./ sum(Np,dims=2)
                Ngen[isnan.(Ngen)].=0
            end
            
        dat[m,:,:]=Np
            
        end
    
    return dat
    
end

function single_sim1(steps,r,K1,K2,a1,A,R,Ng0,Npop)

    nsp=size(Ng0)[1]
    nt=size(Ng0)[2]
        
    Np0=Ng0 .*Npop
    Ngen=deepcopy(Ng0)
    Np=deepcopy(Np0)
    
    dat=zeros(Float64,steps+1,nsp,nt)
    dat[1,:,:]=Np
        
    #Start the simulation
    for m in 2:(steps+1)
    
        #Determine the extinct species
        Np[findall(sum(Np,dims=2) .< 10),:] .= 0
        Ngen[findall(sum(Np,dims=2) .==0),:] .= 0

        if all(sum(Np,dims=2) ==0) 
            break
        else
            newgen=zeros(Float64,nsp,nt)

            #Reproduction event
            for i in findall(!iszero,sum(eachcol(Ngen)))
                
                probs=Ngen[i,:]*Ngen[i,:]'
                
                newgen[i,:]=[sum(probs.*R[:,:,x1]) for x1 in 1:nt]
                
                newgen[i,:] ./= sum(newgen[i,:])
                    
            end

            newgen .*= sum(Np,dims=2)

            #Selection event
            
            if size(A)[1]>1

                for i1 in 1:size(newgen)[1]
                
                    #Reduction in growth rate due to intraspecific competition (not trait-dependent!). The carrying capacity
                    #here is assumed to be the same for all species.

                    rdash=r[i1] * (1 -(sum(newgen[i1,:])/K1[i1]))
                    
                    #Impact of interspecific competition

                    comps=[(a1.*sum((A[x2,:]) .* newgen[1:end .!=i1,:]')) for x2 in 1:nt]
                    
                    Np[i1,:] = newgen[i1,:] + (newgen[i1,:] .* rdash .* (1 .-(comps ./K2)))
                    
                end
                
            end

            Np[findall(Np .<1)] .= 0
            Ngen= Np ./ sum(Np,dims=2)
            Ngen[isnan.(Ngen)].=0
        end
            
        dat[m,:,:]=Np
            
    end

    result=DataFrame()

    for i in 1:nsp

        append!(result,DataFrame(species=i,
                                time=1:steps,
                                N=[sum(dat[x,i,:]) for x in 1:steps],
                                trmean=[sum(dat[x,i,:] .* geno)/sum(res[x,i,:])  for x in 1:tstep]))
    end

    return result
    
end

function single_sim12(steps,r,K1,K2,a1,A,R,Ng0,Npop,samples)

    nsp=size(Ng0)[1]
    nt=size(Ng0)[2]
    geno=collect(range(-1.0,stop=1.0,length=nt))
        
    Np0=Ng0 .*Npop
    Ngen=deepcopy(Ng0)
    Np=deepcopy(Np0)

    result=DataFrame()
        
    #Start the simulation
    for m in 1:steps
    
        #Determine the extinct species
        Np[findall(sum(Np,dims=2) .< 10),:] .= 0
        Ngen[findall(sum(Np,dims=2) .==0),:] .= 0

        if all(sum(Np,dims=2) ==0) 
            break
        else
            newgen=zeros(Float64,nsp,nt)

            #Reproduction event
            for i in findall(!iszero,sum(eachcol(Ngen)))
                
                probs=Ngen[i,:]*Ngen[i,:]'
                
                newgen[i,:]=[sum(probs.*R[:,:,x1]) for x1 in 1:nt]
                
                newgen[i,:] ./= sum(newgen[i,:])
                    
            end

            newgen .*= sum(Np,dims=2)

            #Selection event
            
            if size(A)[1]>1

                for i1 in 1:size(newgen)[1]
                
                    #Reduction in growth rate due to intraspecific competition (not trait-dependent!). The carrying capacity
                    #here is assumed to be the same for all species.

                    rdash=r[i1] * (1 -(sum(newgen[i1,:])/K1[i1]))
                    
                    #Impact of interspecific competition

                    comps=[(a1.*sum((A[x2,:]) .* newgen[1:end .!=i1,:]')) for x2 in 1:nt]
                    
                    Np[i1,:] = newgen[i1,:] + (newgen[i1,:] .* rdash .* (1 .-(comps ./K2)))
                    
                end
                
            end

            Np[findall(Np .<1)] .= 0
            Ngen= Np ./ sum(Np,dims=2)
            Ngen[isnan.(Ngen)].=0
        end

        if m ∈ samples

            append!(result,
                    DataFrame(time=m,
                       sp=1:nsp,
                       pop=pop=vec(sum(Np,dims=2)),
                       trmean=[sum(Ngen[x,:] .* geno) for x in 1:size(Ngen)[1]]))

        end
            
    end

    return result
    
end


function single_sim2(steps,r,K,a1,genos,A,R,Ng0,Npop,samples)

    nsp=size(Ng0)[1]
    nt=size(Ng0)[2]
        
    Np0=Ng0 .*Npop
    Ngen=deepcopy(Ng0)
    Np=deepcopy(Np0)

    res=zeros(Float64,length(samples)+1,nsp,nt)
    res[1,:,:]=Np
    ct=1

    #Start the simulation
    for m in 1:steps
    
        #Determine the extinct species
        Np[findall(sum(Np,dims=2) .< 10),:] .= 0
        Ngen[findall(sum(Np,dims=2) .==0),:] .= 0

        if all(sum(Np,dims=2) ==0) 
            break
        else
            newgen=zeros(Float64,nsp,nt)

            #Reproduction event
            for i in findall(!iszero,sum(eachcol(Ngen)))
                
                probs=Ngen[i,:]*Ngen[i,:]'
                
                newgen[i,:]=[sum(probs.*R[:,:,x1]) for x1 in 1:nt]
                
                newgen[i,:] ./= sum(newgen[i,:])
                    
            end

            newgen .*= sum(Np,dims=2)

            #Selection event
            
            if size(A)[1]>1

                for i1 in 1:size(newgen)[1]
                
                    #Reduction in growth rate due to intraspecific competition (not trait-dependent!). The carrying capacity
                    #here is assumed to be the same for all species.
                    
                    #Impact of interspecific competition

                    comps=[(a1.*sum(A[i1,x2,:,:] .* newgen)) for x2 in 1:nt]

                    #Apply LV equation to each phenotype of each species
                    Np[i1,:] = newgen[i1,:] + (newgen[i1,:] .* r[i1] .* (1 .- (((a1 .* comps) .+ sum(newgen[i1,:]))/K[i1]))) 
                    #Np[i1,:] = newgen[i1,:] + (newgen[i1,:] .* r[i1] .* (1 .-comps .- ((sum(newgen[i1,:])/K[i1]))))
                    
                end
                
            end

            Np[findall(Np .<1)] .= 0
            Ngen= Np ./ sum(Np,dims=2)
            Ngen[isnan.(Ngen)].=0
        end

        if m ∈ samples

            res[ct+1,:,:]=Np
            ct += 1
        end
        
    end

    return res
    
end

function single_sim3(steps,r,K,a1,genos,A,R,Ng0,Npop,samples)

    nsp=size(Ng0)[1]
    nt=size(Ng0)[2]
        
    Np0=Ng0 .*Npop
    Ngen=deepcopy(Ng0)
    Np=deepcopy(Np0)

    result=DataFrame()
    append!(result,
            DataFrame(time=0,
                      sp=1:nsp,
                      pop=pop=vec(sum(Np,dims=2)),
                      trmean=[sum(Ngen[x,:] .* genos[:,x]) for x in 1:size(Ngen)[1]]))
        
    #Start the simulation
    for m in 1:steps
    
        #Determine the extinct species
        Np[findall(sum(Np,dims=2) .< 10),:] .= 0
        Ngen[findall(sum(Np,dims=2) .==0),:] .= 0

        if all(sum(Np,dims=2) ==0) 
            break
        else
            newgen=zeros(Float64,nsp,nt)

            #Reproduction event
            for i in findall(!iszero,sum(eachcol(Ngen)))
                
                probs=Ngen[i,:]*Ngen[i,:]'
                
                newgen[i,:]=[sum(probs.*R[:,:,x1]) for x1 in 1:nt]
                
                newgen[i,:] ./= sum(newgen[i,:])
                    
            end

            newgen .*= sum(Np,dims=2)

            #Selection event
            
            if size(A)[1]>1

                for i1 in 1:size(newgen)[1]
                
                    #Reduction in growth rate due to intraspecific competition (not trait-dependent!). The carrying capacity
                    #here is assumed to be the same for all species.
                    
                    #Impact of interspecific competition

                    comps=[(a1.*sum(A[i1,x2,:,:] .* newgen)) for x2 in 1:nt]

                    #Apply LV equation to each phenotype of each species
                    Np[i1,:] = newgen[i1,:] + (newgen[i1,:] .* r[i1] .* (1 .- (((a1 .* comps) .+ sum(newgen[i1,:]))/K[i1]))) 
                    #Np[i1,:] = newgen[i1,:] + (newgen[i1,:] .* r[i1] .* (1 .-comps .- ((sum(newgen[i1,:])/K[i1]))))
                    
                end
                
            end

            Np[findall(Np .<1)] .= 0
            Ngen= Np ./ sum(Np,dims=2)
            Ngen[isnan.(Ngen)].=0
        end

        if m ∈ samples

            append!(result,
                    DataFrame(time=m,
                       sp=1:nsp,
                       pop=pop=vec(sum(Np,dims=2)),
                       trmean=[sum(Ngen[x,:] .* genos[:,x]) for x in 1:size(Ngen)[1]]))

        end
            
    end

    return result
    
end

function single_sim_nc(steps,R,Ng0)

    nsp=size(Ng0)[1]
    nt=size(Ng0)[2]
    Ngen=deepcopy(Ng0)
    
    #result=DataFrame()
    result=zeros(Float64,steps+1,nsp,nt)
    result[1,:,:]=Ngen 
    
    #Start the simulation
    for m in 1:steps
    
        #Determine extinct species
        #Ngen[findall(sum(Np,dims=2) .==0),:] .= 0

        if all(sum(Ngen,dims=2) ==0) 
            break
        else
            newgen=zeros(Float64,nsp,nt)

            #Reproduction event
            for i in findall(!iszero,sum(eachcol(Ngen)))
                
                prob=Ngen[i,:]*Ngen[i,:]'
                
                newgen[i,:]=[sum(prob*R[:,:,x1]) for x1 in 1:nt]
                
                newgen[i,:] ./= sum(newgen[i,:])
                    
            end

            Ngen=newgen
            #Ngen[findall(Ngen .<0)] .= 0
            #Ngen[isnan.(Ngen)].=0
            
        end

        result[m+1,:,:]=Ngen
    end

    return result
    
end

function SK_proc(steps,R,Np)

    nt=length(Np)
           
    Np ./= sum(Np)
    
    res=zeros(Float64,steps+1,nt)
    res[1,:]=Np

    for m in 1:steps

        #Reproduction event   
        prob1=Np*Np'
        newp=[sum(prob1*R[:,:,x1]) for x1 in 1:nt]
        Np = newp ./ sum(newp)
        Np[findall(Np .<0.001)] .= 0
        Np[isnan.(Np)].=0
        res[m+1,:]=Np
    end
    return res
end

#Set working directory

dir="C:\\Users\\mihir\\Documents\\Comp_Coevol\\data\\"

#############################################################################
#Try out comp kernel shapes
############################################################################
#comp. kernel parameters
omegas=collect(0.05:0.1:1.0)
t=0.5

xs=collect(-1:0.05:1)

kernelp=plot()

for i in omegas
    ys=alpha_gt.(xs,0.0,i,t)
    plot!(kernelp,xs,ys)
end

kernelp

############################################################################
#Demonstrate when convergence occurs
############################################################################

# 3 species case

#Premise: We first fix the trait distributions of species A. We then vary mean traits of species B and C
#and mark when species A and B converge.

#Different comp. kernel functions do not seem to have qualitative effect on evolutionary dynamics. So,
# at this point, we will just use truncated Gaussian kernels.

#Parameters to vary: 1. #loci, 2. kernel omega 3. comp. threshold 3. Trait distributions
nsp=3  #No. of species
ns=[3,5,7,10,20] #No. of loci
#kernel=["Gaussian","Triangle"]
omegas=collect(0.05:0.2:1.0)
ts=collect(0.1:0.2:0.8)
dists=["Gaussian","Uniform"]
a1s=[0.0,0.00001,0.0001,0.001,0.01]
reps=50
t1s=-0.6:0.05:0.0
t2s=0.0:0.05:0.6
trsets=collect(Iterators.product(t1s,t2s))

res=DataFrame()

for i1 in ns, i2 in dists, i3 in omegas, i4 in ts, i5 in 1:length(trsets)

    nt=2*i1 + 1
    trmeans=[trsets[i5][1],0.0,trsets[i5][2]]
    
    #Assigning trait distributions around given means for 3 species.
    genos=zeros(Float64,nsp,nt)
    for i in 1:nsp
        genos[i,:]=collect(range(trmeans[i]-0.5,trmeans[i]+0.5,length=nt))
    end

    Ns=zeros(Float64,nsp,nt)

    if i2=="Gaussian"
        [Ns[x,:]=pdf.(truncated(Normal(trmeans[x],0.2),-1.0,1.0),genos[x,:]) for x in 1:nsp]
    elseif i2=="Uniform"
        [Ns[x,:]=pdf.(Uniform(trmeans[x]-0.4,trmeans[x]+0.4),genos[x,:]) for x in 1:nsp]
    end
    Ns= Ns ./ sum(Ns,dims=2)

    A=zeros(Float64,nsp,nt,nsp,nt)
    for j1 in 1:nsp, j2 in 1:nt, j3 in 1:nsp, j4 in 1:nt
        if(j1 != j3)
            A[j1,j2,j3,j4]=alpha_gt(genos[j1,j2],genos[j3,j4],i3,i4)
        end
    end

    comps=zeros(Float64,nsp,nt)
    for j5 in 1:nsp, j6 in 1:nt
        comp=0.0
        for j7 in 1:nsp, j8 in 1:nt
            comp = comp + A[j5,j6,j7,j8]*Ns[j5,j6]*Ns[j7,j8]
        end
        comps[j5,j6]=comp
    end
    comps=comps ./ [maximum(comps[x,:]) for x in 1:nsp]      
    Ns=Ns .* (1 .- comps)
    newmeans=[sum(genos[x,:] .* Ns[x,:])/sum(Ns[x,:]) for x in 1:nsp]

    diffs=diff(newmeans)-diff(trmeans)

    append!(res,DataFrame(rep=i5,loci=i1,dist=i2,omega=i3,t=i4,sp1=trmeans[1],
                        sp3=trmeans[3],diff1=diffs[1],diff2=diffs[2]))
end

CSV.write("3sp_test6.csv",res)


####################################
#Trial run

nsp=5
n=5
omega=0.2
t=0.5
R=qgprob(n)
nt=2*n + 1
a1=0.1

#Set up trait distributions of species populations
#Each population has #nt discrete phenotypes centered around a specific trait value. Under no selection pressure, populations' trait means converge to this value.
#Important note: This assumption DOES NOT imply stabilizing selection, just an outcome of the quantitative genetic architecture of the trait. To further stress this
# point, sets of populations of #nsp species will be assigned randomized mean trait values.

tmeans=rand(Uniform(-0.6,0.6),nsp)
#create all possible phenotypes for all species
genos=zeros(Float64,nt,nsp)
for i in 1:nsp
    genos[:,i]=collect(range(tmeans[i]-0.5,tmeans[i]+0.5,length=nt))
end

#Create alpha matrix 
A=zeros(Float64,nsp,nt,nsp,nt)
for j1 in 1:nsp, j2 in 1:nt, j3 in 1:nsp, j4 in 1:nt
    if(j1 != j3)
        A[j1,j2,j3,j4]=alpha_gt(genos[j2,j1],genos[j4,j3],omega,t)
    end
end

Npop=fill(10000.0,nsp)
K=fill(20000,nsp)
r=rand(Uniform(0.0,0.1),nsp)

#Assign initial trait distributions to species
Ng0=zeros(Float64,nsp,nt)
#Gaussian
[Ng0[x,:]=pdf.(truncated(Normal(tmeans[x],0.2),-1.0,1.0),genos[:,x]) for x in 1:nsp]

#Uniform
[Ng0[x,:]=pdf.(Uniform(tmeans[x]-0.5,tmeans[x]+0.5),genos[:,x]) for x in 1:nsp]

Ng0=Ng0 ./ sum(Ng0,dims=2)

initd=plot()
for i in 1:nsp
    plot!(initd,genos[:,i],Ng0[i,:])
end
initd


steps=5000 #keep 1000 at minimum
samples=[collect(1:9); collect(10:10:100); collect(100:100:steps)]

res0=single_sim3(steps,r,K,a1,genos,A,R,Ng0,Npop,samples)

trdat0=zeros(Float64,size(res0)[1],size(res0)[2])
popdat0=zeros(Float64,size(res0)[1],size(res0)[2])
for i in 1:size(res0)[1]
    [popdat0[i,x1]=sum(res0[i,x1,:]) for x1 in 1:size(res0)[2]]
    [trdat0[i,x2] = sum(res0[i,x2,:] .* genos[:,x2])/sum(res0[i,x2,:]) for x2 in 1:size(res0)[2]]
end

plot([0;samples],trdat0,
    marker=:circle,
    markersize=popdat0 ./2000)


anim= @animate for i in 1:size(res0)[1]
        plot(res0[i,:,:])
    end

gif(anim,fps=2)



#################################################################################
#Run 20 spp simulations
#################################################################################

#Parameters to vary
reps=1:5
nsp=20
ns=[5,7,10,20,30,50] #No. of loci
kernels=["Gaussian","Triangle"]
omegas=[0.1,0.2,0.5,1.0]
ts=[0.1,0.2,0.5,1.0]
dists=["Gaussian","Uniform"]
a1s=[0.1,0.25,0.5,1.0]

n=5
inits=[rand(Uniform(-1.0,1.0),nsp) for i in reps]
parset=collect(Iterators.product(reps,omegas,ts,a1s,dists))

#Fix demographic parameters
r= abs.(rand(Normal(0.1,0.05),nsp))
Npop=fill(1000.0,nsp)
K=fill(2000.0,nsp)


geno=collect(range(-1.0,stop=1.0,length=2*n +1))
nt=length(geno)
K2=fill(K2/nt,nt)
a1=0.1
A=zeros(Float64,nt,nt)

    for j1 in 1:nt, j2 in 1:nt
        A[j1,j2]=alpha_gt(geno[j1],geno[j2],omega,t)
    end

R=qgprob(n)
Ng0=zeros(Float64,nsp,nt)
rands=rand(Uniform(-0.6,0.6),nsp)
[Ng0[i,:]=pdf.(truncated(Normal(rands[i],0.2),-1.0,1.0),geno) for i in 1:nsp]
Ng0=Ng0 ./ sum(Ng0,dims=2)


tsteps=5000
res1=getsum(tsteps,r,K1,K21,0.1,A,R,Ng0,Npop)
pop1=res1[1]
trdat1=res1[2]

df1=stack(DataFrame(trdat1,:auto))
df1.time=repeat(0:tsteps,outer=nsp)
df1.species=repeat(1:nsp,inner=1+tsteps)
df1=df1[:,[:time,:species,:value]]


plot(df1.time,df1.value,
    group=df1.species)



