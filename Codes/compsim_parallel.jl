using Distributed

addprocs(11)

@everywhere begin

    using CSV, DataFrames, Random, LinearAlgebra, Distances, Distributions, SpecialFunctions


    cd1="C:\\Users\\mihir\\Documents\\Comp_Coevol\\data"

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
    omegas=[0.1,0.2,0.5,0.75,1.0]
    ts=[0.1,0.2,0.5,1.0,2.0]
    dist=["Gaussian","Uniform"]
    a1s=[0.0,0.1,0.5,1]

    #create diff files for the following param combos
    ns=[3,5,7,10,20] #No. of loci
    reps=1:5
    inits=[rand(Uniform(-0.8,0.8),nsp) for i in reps] #Initial trait distributions for species

    pars=collect(Iterators.product(ns,reps))

    steps=2000 #keep 1000 at minimum
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

        #Demographic params
        r=rand(Uniform(0.0,0.1),nsp)
        Npop=fill(1000.0,nsp)
        K=fill(2000.0,nsp)

        result=DataFrame()

        for i1 in omegas,i2 in ts, i3 in dist, i4 in a1s 

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
