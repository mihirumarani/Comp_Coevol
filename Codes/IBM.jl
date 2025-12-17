#IBM simulations

#Quantitative genetic model:
# N independent bialleilic loci contribute additively to the phenotypic value. 
# For each locus i, one allele "0" contributes 0 to the phenotype while the other allele,
# "1" contribute some non-zero amount to the phenotype.
#  Loci may or may not contribute equally to the phentypic value.

using CSV, DataFrames, Random, LinearAlgebra, Distances, Distributions, SpecialFunctions, Plots, StatsPlots,
FreqTables,Combinatorics


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

function pheno_from_geno(a::Matrix{Tuple{Float64, Float64}},
                        eff::Vector{Float64},
                        tmin::Float64)
        
        return(a|> 
            x-> sum.(x)|>
            x-> tmin .+ map(y-> sum(x[y,:] .* eff), 1:size(a)[1]) |>
            x-> round.(x,digits=2))
end
        

#Set up parameter space

#Competition kernel
omega=0.5  #kernel width
t=10.0     #trait threshold

nsp= 3 #no. of species
Np=repeat([100],nsp) #Population sizes
tstep=100

nloci=5   #No. of loci    

#Effect size of each locus: try 4 different combinations
eff= ones(Float64,nloci,4)
eff[:,2] += [1,0,0,0,0]  #One locus slightly more important than others
eff[:,3] += [4,0,0,0,0]  #One locus contributes much more than others
eff[:,4] += [1,1,0,0,0]  #Two loci more important than other 3.
eff = eff ./ sum(eff,dims=1)

#Initial trait means
tmins=[-0.5,0.0,0.2]

#nitiate the populations 
copy1=sample([0.0,0.1],nsp*nloci*Np[1],replace=true)
copy2=sample([0.0,0.1],nsp*nloci*Np[1],replace=true)
gametes=tuple.(copy1,copy2)
probs=rand(Uniform(0.3,0.7),3)

init=[reshape(tuple.(wsample([0.0,0.1],[probs[1],1-probs[1]],nloci*Np[1],replace=true),
            wsample([0.0,0.1],[probs[1],1-probs[1]],nloci*Np[1],replace=true)),(Np[1],nloci)),
    reshape(tuple.(wsample([0.0,0.1],[probs[2],1-probs[2]],nloci*Np[2],replace=true),
            wsample([0.0,0.1],[probs[2],1-probs[2]],nloci*Np[2],replace=true)),(Np[2],nloci)),
    reshape(tuple.(wsample([0.0,0.1],[probs[3],1-probs[3]],nloci*Np[3],replace=true),
            wsample([0.0,0.1],[probs[3],1-probs[3]],nloci*Np[3],replace=true)),(Np[3],nloci))]


#Store result here
res=fill(DataFrame(eff=[],time=[],pheno=[]),nsp)

tsteps=100


for i in 1:size(eff)[2]

    pop=deepcopy(init)

     ##Put the initial population data into the result file
    [res[k]=append!(res[k],
                    DataFrame(eff=i,time=0,pheno=pheno_from_geno(pop[k],eff[:,i],tmins[k]))) for k in 1:nsp]

    #Find all possible phenotypes for each species and construct the matrix
    #of competition effect between each pair of interspecific phenotypes.
    gtypes= vec(collect(Iterators.product(fill([0,0.1,0.2],nloci)...)))
    pts=[gtypes |> 
        x->map(y -> eff[:,i] .* x[y], 1:length(x)) |>
        x->tmins[i0] .+ sum.(x) |>
        x->round.(x,digits=2) |>
        unique |> sort 
        for i0 in 1:nsp]  

    nt=length(pts[1])

    A=zeros(Float64,nsp,nt,nsp,nt)
    for j1 in 1:nsp, j2 in 1:nt, j3 in 1:nsp, j4 in 1:nt
        if(j1 != j3)
            A[j1,j2,j3,j4]=alpha_gt(pts[j1][j2],pts[j3][j4],omega,t)
        end
    end

    for i1 in 1:tsteps

        changes=[ceil(Int,0.2*size(pop[x])[1]) for x in 1:nsp]

        #Reproduction process: Random mating, add 0.2 times more individuals
        for i2 in 1:nsp

            #Separate the gametes
            c1=[first.(pop[i2][:,i3]) for i3 in 1:nloci]
            c2=[last.(pop[i2][:,i4]) for i4 in 1:nloci]

            #Add 20% individuals by randomly sampling gametes and combining them.

            nxgen=[tuple.(sample(c1[i5],changes[i2],replace=true),
                        sample(c2[i5],changes[i2],replace=true))   for i5 in 1:nloci]
            
            pop[i2]=vcat(pop[i2],reshape(collect(Iterators.flatten(nxgen)),(:,nloci)))
        end

        #Selection process: calculate reduced fitness due to competition
        pop1=[sum.(pop[x]) for x  in 1:nsp]

        popfreqs=[pop1[i6] |> 
                x-> map(y-> sum(eff[:,i] .* x[y,:]), 1:size(x)[1]) |> 
                x-> round.(tmins[i6] .+ x,digits=2)
                for i6 in 1:nsp]
                
        pfreqs=[popfreqs[i7] |> 
                freqtable |>
                x-> DataFrame(phenos=names(x)[1],freq=x) |>
                x-> [x;DataFrame(phenos=setdiff(pts[i7],x.phenos),freq=0)] |>
                x-> sort(x,:phenos)
                for i7 in 1:nsp]

        [pfreqs[i8][:,:prob] = pfreqs[i8].freq ./ sum(pfreqs[i8].freq) for i8 in 1:nsp]

        freqmat=zeros(Float64,nsp,nt)
        [freqmat[i9,:]= pfreqs[i9].prob for i9 in 1:nsp]

        for i10 in 1:nsp
            comps= [(sum(A[i10,i11,:,:] .* freqmat)) for i11 in 1:nt]
            d1=Dict(pts[i10][i12]=>comps[i12] for i12 in 1:length(pts[i10]))
            prob1=[get(d1,popfreqs[i10][i13],0.0) for i13 in 1:length(popfreqs[i10])]
            pop[i10]=pop[i10][setdiff(1:size(pop[i10])[1],
                                wsample(1:length(prob1),prob1,changes[i10],replace=false)),:]
        end

        ##Put the population data into the result file
        [res[k]=append!(res[k],
                        DataFrame(eff=i,time=i1,pheno=pheno_from_geno(pop[k],eff[:,i],tmins[k]))) for k in 1:nsp]

    end
end

res1=res[1]
res2=res[2]
res3=res[3]

result=DataFrame()

append!(result,DataFrame(sp=1,eff=res1.eff,time=res1.time,pheno=res1.pheno))
append!(result,DataFrame(sp=2,eff=res2.eff,time=res2.time,pheno=res2.pheno))
append!(result,DataFrame(sp=3,eff=res3.eff,time=res3.time,pheno=res3.pheno))

CSV.write("D:\\project_files\\comp_coevol\\data\\ibm1.csv",result)

anim= @animate for i in 0:tsteps

    dat1=res1[res1.eff .==1,:]
    dat2=res2[res2.eff .==1,:]
    dat3=res3[res3.eff .==1,:]

    density(dat1[dat1.time .==i,:].pheno)
    density!(dat2[dat2.time .==i,:].pheno)
    density!(dat3[dat3.time .==i,:].pheno)
    
    end

gif(anim,fps=2)