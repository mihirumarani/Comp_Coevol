using Distributed

addprocs(11)

@everywhere begin

    using CSV, DataFrames, Random, LinearAlgebra, Distances, Distributions, SpecialFunctions


    cd1="D:\\Project_files\\comp_coevol\\data\\20sp\\"

    #Load functions
    function alpha_gt(x::Float64, y::Float64, omega::Float64, t::Float64)
    
        al=0
        if(abs(x-y)<t) 
             al=exp(-((x-y)^2)/(omega^2))
             #denom=(pi/2)*(erf(t/omega)-erf(-t/omega))
             #al/=denom
        end
        return al
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

    function Amat(genos,nsp,nt,omega,t)

        A=zeros(Float64,nsp,nt,nsp,nt)

        for j1 in 1:nsp, j2 in 1:nt, j3 in 1:nsp, j4 in 1:nt
            if(j1 != j3)
                A[j1,j2,j3,j4]= alpha_gt(genos[j2,j1],genos[j4,j3],omega,t)
            end
        end
        return A
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

    #Parameters to vary

    nsp=20
    #kernels=["Gaussian","Triangle"]
    omegas=[0.1,0.2,0.5,1.0]
    ts=[0.1,0.2,0.5,1.0]
    dists=["Gaussian","Uniform"]
    a1s=[0.1,0.25,0.5,1.0]

    #Demographic params
    r=rand(Uniform(0.0,0.1),nsp)
    Npop=fill(1000.0,nsp)
    K=fill(2000.0,nsp)

    #create diff files for the following param combos
    ns=[10] #No. of loci
    reps=1:50
    inits=[rand(Uniform(-0.8,0.8),nsp) for i in reps] #Initial mean trait values for species

    pars=collect(Iterators.product(ns,reps))

    steps=5000 #keep 1000 at minimum
    samples=[collect(1:9); collect(10:10:100); collect(100:100:steps)]

    #Final boss to run
    function compsim(par::Tuple{Int64, Int64})

        loci=par[1]
        nt=2*loci + 1

        #Set up possible phenotypes for all species
        tmeans=inits[par[2]]
        genos=zeros(Float64,nt,nsp)
        for i in 1:nsp
            genos[:,i]=collect(range(tmeans[i]-0.5,tmeans[i]+0.5,length=nt))
        end
        
        #Mating prob. matrix
        R=qgprob(loci)

        result=DataFrame()

        for i1 in omegas,i2 in ts, i3 in dists, i4 in a1s 

            #Create alpha matrix 
            A=Amat(genos,nsp,nt,i1,i2)

            #Assign initial trait distributions to species
            Ng0=zeros(Float64,nsp,nt)

            if i3=="Gaussian"
                [Ng0[x,:]=pdf.(truncated(Normal(tmeans[x],0.2),tmeans[x]-0.5,tmeans[x]+0.5),genos[:,x]) for x in 1:nsp]
            elseif i3=="Uniform"
                [Ng0[x,:]=pdf.(Uniform(tmeans[x]-0.5,tmeans[x]+0.5),genos[:,x]) for x in 1:nsp]
            end
            Ng0=Ng0 ./ sum(Ng0,dims=2)

            res=single_sim3(steps,r,K,i4,genos,A,R,Ng0,Npop,samples)
            res2=unstack(res,:time,:sp,:trmean)

            res.dist .= i3
            res.a .= i4
            res.omega .= i1
            res.t .= i2
            append!(result,res)
        end

        result.loci .= loci
        result.rep .= par[2]

        CSV.write(string(cd1,loci,"_",par[2],".csv"),result)
    end    

end

pmap(compsim,pars)
 

#######################################################################################
#For the following code, restart the julia kernel. No need to parallelization.

using CSV, DataFrames, Random, LinearAlgebra, Distances, Distributions, SpecialFunctions

loc5="D:\\project_files\\comp_coevol\\data\\5sp\\"
loc10="D:\\project_files\\comp_coevol\\data\\10sp\\"
loc20="D:\\project_files\\comp_coevol\\data\\20sp\\"

dat5=DataFrame()
files5=readdir(loc5)
for i in files5
    append!(dat5,
            CSV.read(string(loc5,i),DataFrame))
end
dat5.nsp .= 5


dat10=DataFrame()
files10=readdir(loc10)
for i in files10
    append!(dat10,
            CSV.read(string(loc10,i),DataFrame))
end
dat10.nsp .= 10

dat20=DataFrame()
files20=readdir(loc20)
for i in files20
    append!(dat20,
            CSV.read(string(loc20,i),DataFrame))
end
dat20.nsp .= 20

dat=append!(dat5,dat10,dat20)

nsps=[5,10,20]

trmean1=[dat.pop[i] .>1 ? dat.trmean[i] : missing for i in 1:size(dat)[1]]
dat.trmean1 = trmean1

convdat=DataFrame()

for i in nsps

    dat0=dat[dat.nsp .== i, :]
    pars=unique(dat0[:,5:10])

    for j in 1:nrow(pars)

        dat1=dat0[
            dat0.loci .== pars.loci[j] .&&
            dat0.dist .== pars.dist[j] .&&
            dat0.a .== pars.a[j] .&&
            dat0.omega .== pars.omega[j] .&&
            dat0.t .== pars.t[j] .&&
            dat0.rep .== pars.rep[j], :]

        reord=sortperm(dat1[dat1.time .==0,:].trmean1)
        times=unique(dat1.time)
        tsteps=length(times)
        dat1= dat1[(sort(repeat(i .* (1:tsteps),i)) .- i) + repeat(reord,tsteps),[:time,:sp,:trmean1]]

        diffdat=DataFrame()
        
        for k in times
            dat2=dat1[dat1.time .== k,:]
            td=collect(skipmissing(diff(dat2.trmean1)))

            if length(td)>0
                append!(diffdat,
                DataFrame(time=k,
                sppair=1:length(td),
                td=td))
            end
        end

        sppairs=1:(i-1)
        diffdat2=DataFrame()
        for l in sppairs
            dat2=diffdat[diffdat.sppair .== l,:]
            td=diff(dat2.td)
            append!(diffdat2,
            DataFrame(sppair=l,
            times=1:length(td),
            td=td))
        end
        conv=sum(diffdat2.td .< 0)/nrow(diffdat2)

        append!(convdat,
        DataFrame(nsp=i,
                loci=pars.loci[j],
                dist=pars.dist[j],
                a=pars.a[j],
                omega=pars.omega[j],
                t=pars.t[j],
                rep=pars.rep[j],
                conv=conv))
    end
end

CSV.write(string("D:\\project_files\\comp_coevol\\data\\","convdat.csv"),convdat)



########################################################################################
#Recalculate the convergence data only for the first 100 time steps

using CSV, DataFrames, Random, LinearAlgebra, Distances, Distributions, SpecialFunctions

loc5="D:\\project_files\\comp_coevol\\data\\5sp\\"
loc10="D:\\project_files\\comp_coevol\\data\\10sp\\"
loc20="D:\\project_files\\comp_coevol\\data\\20sp\\"

dat5=DataFrame()
files5=readdir(loc5)
for i in files5
    append!(dat5,
            CSV.read(string(loc5,i),DataFrame))
end
dat5.nsp .= 5

dat10=DataFrame()
files10=readdir(loc10)
for i in files10
    append!(dat10,
            CSV.read(string(loc10,i),DataFrame))
end
dat10.nsp .= 10

dat20=DataFrame()
files20=readdir(loc20)
for i in files20
    append!(dat20,
            CSV.read(string(loc20,i),DataFrame))
end
dat20.nsp .= 20

dat=append!(dat5,dat10,dat20)


nsps=[5,10,20]

trmean1=[dat.pop[i] .>1 ? dat.trmean[i] : missing for i in 1:size(dat)[1]]
dat.trmean1 = trmean1

convdat=DataFrame()

for i in nsps

    dat0=dat[dat.nsp .== i, :]
    pars=unique(dat0[:,5:10])

    for j in 1:nrow(pars)

        dat1=dat0[
            dat0.loci .== pars.loci[j] .&&
            dat0.dist .== pars.dist[j] .&&
            dat0.a .== pars.a[j] .&&
            dat0.omega .== pars.omega[j] .&&
            dat0.t .== pars.t[j] .&&
            dat0.rep .== pars.rep[j], :]

        reord=sortperm(dat1[dat1.time .==0,:].trmean1)
        times=unique(dat1.time)
        tsteps=length(times)
        dat1= dat1[(sort(repeat(i .* (1:tsteps),i)) .- i) + repeat(reord,tsteps),[:time,:sp,:trmean1]]

        dat2=Matrix(unstack(dat1,:time,:trmean1)[:,Not(1)])

        etimes=[findfirst(ismissing.(dat2[i,:])) for i in 1:size(dat2)[1]]

        etimes=replace(etimes,nothing=>tsteps+1)

        etimes=DataFrame(sp=1:length(etimes),ext=etimes)

        sort!(etimes,order(:ext))

        convs=DataFrame()
        
        if minimum(etimes.ext) < tsteps

            samp=deepcopy(dat2)
            counter=1            
            
            for k in 1:findfirst(etimes.ext .==tsteps+1)
                
                samp=dat2[:,counter: etimes.ext[k]-1]
                
                if size(samp)[2] >=2

                    if counter !=1
                        samp=samp[Not(etimes.sp[1:k-1]),:]
                    end

                    diffdat1=Matrix(undef,size(samp)[1]-1,size(samp)[2])
                    [diffdat1[:,i1]=diff(samp[:,i1]) for i1 in 1:size(samp)[2]]
                    diffdat2=Matrix(undef,size(diffdat1)[1],size(diffdat1)[2]-1)
                    [diffdat2[i2,:] =diff(diffdat1[i2,:]) for i2 in 1:size(diffdat1)[1]]
                    conv=[sum(diffdat2[:,i3] .<0)/size(diffdat2)[1] for i3 in 1:size(diffdat2)[2]]

                    append!(convs,
                    DataFrame(time=times[counter:(etimes.ext[k]-2)],
                            conv=conv))
                end

                counter=etimes.ext[k]
            end
        else

            diffdat1=Matrix(undef,size(dat2)[1]-1,size(dat2)[2])
            [diffdat1[:,i4]=diff(dat2[:,i4]) for i4 in 1:size(dat2)[2]]
            diffdat2=Matrix(undef,size(diffdat1)[1],size(diffdat1)[2]-1)
            [diffdat2[i5,:] =diff(diffdat1[i5,:]) for i5 in 1:size(diffdat1)[1]]
            conv=[sum(diffdat2[:,i6] .<0)/size(diffdat2)[1] for i6 in 1:size(diffdat2)[2]]

            append!(convs,
                DataFrame(time=times[1:(length(times)-1)],
                        conv=conv))
        end

        
        append!(convdat,
        DataFrame(nsp=i,
                loci=pars.loci[j],
                dist=pars.dist[j],
                a=pars.a[j],
                omega=pars.omega[j],
                t=pars.t[j],
                rep=pars.rep[j],
                time=convs.time,
                conv=convs.conv))
    end
end

CSV.write(string("D:\\project_files\\comp_coevol\\data\\","convdatfull.csv"),convdat)


####################################################################################
#Calculate MNND and extinction values

using CSV, DataFrames, Random, LinearAlgebra, Distances, Distributions, SpecialFunctions

function mnnds(a)

    a=a[Not(ismissing.(a))]

    if length(unique(a)) > 1
        
        a=sort(a,rev=true)
        nnd=zeros(Float64,length(a))
        nnd[1]=a[1]-a[2]
        nnd[length(nnd)]=a[length(a)-1]-a[length(a)]
        [nnd[i] = min(a[i-1]-a[i],a[i]-a[i+1]) for i in 2:(length(a)-1)]
        ranges=a[1]-a[length(a)]
        mnnd=sum(nnd)/length(nnd)
        mmax=ranges/(length(a)-1)
        return mnnd/mmax

    else
        return missing
    end
end

loc5="D:\\project_files\\comp_coevol\\data\\5sp\\"
loc10="D:\\project_files\\comp_coevol\\data\\10sp\\"
loc20="D:\\project_files\\comp_coevol\\data\\20sp\\"

dat5=DataFrame()
files5=readdir(loc5)
for i in files5
    append!(dat5,
            CSV.read(string(loc5,i),DataFrame))
end
dat5.nsp .= 5


dat10=DataFrame()
files10=readdir(loc10)
for i in files10
    append!(dat10,
            CSV.read(string(loc10,i),DataFrame))
end
dat10.nsp .= 10

dat20=DataFrame()
files20=readdir(loc20)
for i in files20
    append!(dat20,
            CSV.read(string(loc20,i),DataFrame))
end
dat20.nsp .= 20

dat=append!(dat5,dat10,dat20)

trmean1=[dat.pop[i] .>1 ? dat.trmean[i] : missing for i in 1:size(dat)[1]]
dat.trmean1 = trmean1

mnndat=combine(groupby(dat,[:nsp,:loci,:dist,:a,:omega,:t,:rep,:time,]), :trmean1 => (x -> mnnds(x)) => :mnnds)

CSV.write(string("D:\\project_files\\comp_coevol\\data\\","mnndat.csv"),mnndat)


extdat=combine(groupby(dat,[:nsp,:loci,:dist,:a,:omega,:t,:rep,:time,]), :pop => (x -> sum(x .< 1)) => :extinct)

CSV.write(string("D:\\project_files\\comp_coevol\\data\\","extdat.csv"),extdat)

        

    

        