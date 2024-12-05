using CSV, DataFrames, Random, LinearAlgebra, Distances, Distributions, SpecialFunctions, Plots, GLMakie

#Makie tutorial: https://www.youtube.com/watch?v=odpoatozNz8&ab_channel=doggodotjl
#Make interactive plots to demonstrate the context under which convergence occurs

#Necessary function

function alpha_gt(x::Float64, y::Float64, omega::Float64, t::Float64)
    
    al=0
    if(abs(x-y)<t) 
         al=exp(-((x-y)^2)/(omega^2))
         #denom=(pi/2)*(erf(t/omega)-erf(-t/omega))
         #al/=denom
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

function single_sim1(time,r,K1,K2,a1,A,R,Ng0,Npop)
    
        nsp=size(Ng0)[1]
        nt=size(Ng0)[2]
        
        Np0=Ng0 .*Npop
        Ngen=deepcopy(Ng0)
        Np=deepcopy(Np0)
    
        dat=zeros(Float64,time+1,nsp,nt)

        dat[1,:,:]=Np
        
        #Start the simulation
        for m in 2:(time+1)
    
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

                        rdash=r[i1] * (1 -(sum(newgen[i1,:])/K2))
                    
                        #Impact of interspecific competition

                        comps=[(a1.*sum((A[x2,:]) .* newgen[1:end .!=i1,:]')) for x2 in 1:nt]
                    
                        Np[i1,:] = newgen[i1,:] + (newgen[i1,:] .* rdash .* (1 .-(comps ./K1)))
                    
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


function getsum(time,r,K1,K2,a1,A,R,Ng0,Npop)
        
    res=single_sim1(time,r,K1,K2,a1,A,R,Ng0,Npop)
    
    tstep=size(res)[1]
    nsp=size(res)[2]
    nt=size(res)[3]
            
    geno=collect(range(-1.0,stop=1.0,length=nt))
    
    pops=zeros(Float64,tstep,nsp)
    trmeans=zeros(Float64,tstep,nsp)

    for i in 1:nsp

    pops[:,i]=[sum(res[x,i,:]) for x in 1:tstep]
    trmeans[:,i]=[sum(res[x,i,:] .* geno)/sum(res[x,i,:])  for x in 1:tstep]

    end

    return([pops,trmeans])
    
end


#############################################################################
#Try out comp kernel shapes
############################################################################
GLMakie.activate!()
#comp. kernel parameters
omega=.7
t=0.5

xs=-1:0.05:1
ys1=[alpha_gt(xs[i],0.0,omega,t) for i in 1:length(xs)]
ys2=[alpha_tri(xs[i],0.0,omega,t) for i in 1:length(xs)]

Plots.plot(xs,ys1)
Plots.plot!(xs,ys2)


fig = Figure(size = (3840, 2160))

ax1 = fig[1, 1] = Axis(fig,
    # borders
    aspect = 1,
    # title
    title = "Competition kernel shape",
    titlegap = 48, titlesize = 60,
    # x-axis
    xautolimitmargin = (0, 0), xgridwidth = 2, xticklabelsize = 36,
    xticks = LinearTicks(20), xticksize = 18,
    # y-axis
    yautolimitmargin = (0, 0), ygridwidth = 2, yticklabelpad = 14,
    yticklabelsize = 36, yticks = LinearTicks(20), yticksize = 18
)

lsgrid= SliderGrid(fig,
        (label= "Omega", range=0.1:0.25:1.5,startvalue=0.75),
        (label= "Threshold", range=0.2:0.1:1.5,startvalue=0.7),
        width = 350,
        tellheight = false
    )