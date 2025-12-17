#IBM simulations

#Quantitative genetic model:
# N independent bialleilic loci contribute additively to the phenotypic value. 
# For each locus i, one allele "0" contributes 0 to the phenotype while the other allele,
# "1" contribute some non-zero amount to the phenotype.
#  Loci may or may not contribute equally to the phentypic value.

using CSV, DataFrames, Random, LinearAlgebra, Distances, Distributions, SpecialFunctions, Plots, StatsPlots,
FreqTables,Combinatorics, AlgebraOfGraphics,CairoMakie


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

#Function to create multiloci bialleilic diploid populations (for n species) with randomly assigned alleles (with fixed frequencies)
#Each individual is described as a vector of tuples (containing the presence/absence of A/a allele) of length equalling number of loci. 
function get_genpop(nsp::Int64,
                    Np::Vector{Int64},   #population sizes for each species.
                    nloci::Int64,
                    probs::Vector{Float64}  #Frequency of A allele for each species. 
                    )

    init=[reshape(tuple.(wsample([0.0,0.1],[probs[i],1-probs[i]],nloci*Np[i],replace=true),
            wsample([0.0,0.1],[probs[i],1-probs[i]],nloci*Np[i],replace=true)),(Np[1],nloci)) 
            for i in 1:nsp]

    return init
end

#Function to get phenotypes of individuals described as explicit bialleilic genotypes.
function pheno_from_geno(a::Matrix{Tuple{Float64, Float64}},
                        eff::Vector{Float64},
                        tmin::Float64)
        
        return(a|> 
            x-> sum.(x)|>
            x-> tmin .+ map(y-> sum(x[y,:] .* eff), 1:size(a)[1]) |>
            x-> round.(x,digits=2))
end

#Function to simulate population and trait dynamics of competing species
function simIBM(init::Vector{Matrix{Tuple{Float64,Float64}}},
                eff::Vector{Float64},
                tmins::Vector{Float64},
                omega::Float64,
                t::Float64,
                tsteps::Int64)

    #Store result here
    res=DataFrame()

    pop=deepcopy(init)
    
    nloci=size(init[1])[2]
    eff ./= sum(eff)

    ##Put the initial population data into the result file
    for k in 1:nsp
        append!(res,
                DataFrame(time=0,sp=k,pheno=pheno_from_geno(pop[k],eff,tmins[k])))
    end


    #Find all possible phenotypes for each species and construct the matrix
    #of competition effect between each pair of interspecific phenotypes.
    gtypes= vec(collect(Iterators.product(fill([0,0.1,0.2],nloci)...)))
    pts=[gtypes |> 
        x->map(y -> eff .* x[y], 1:length(x)) |>
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
        pop1=[sum.(pop[x]) for x  in 1:nsp]    #Contribution from each locus (column) for each individual (row)

        popfreqs=[pop1[i6] |> 
                x-> map(y-> sum(eff .* x[y,:]), 1:size(x)[1]) |> 
                x-> round.(tmins[i6] .+ x,digits=2)
                for i6 in 1:nsp]               #Phenotype of each individual of each species                
                
        pfreqs=[popfreqs[i7] |> 
                freqtable |>
                x-> DataFrame(phenos=names(x)[1],freq=x) |>
                x-> [x;DataFrame(phenos=setdiff(pts[i7],x.phenos),freq=0)] |>
                x-> sort(x,:phenos)
                for i7 in 1:nsp]                        

        [pfreqs[i8][:,:prob] = pfreqs[i8].freq ./ sum(pfreqs[i8].freq) for i8 in 1:nsp]  #Relative frequcies of each phenotype for each species

        freqmat=zeros(Float64,nsp,nt)
        [freqmat[i9,:]= pfreqs[i9].prob for i9 in 1:nsp]

        for i10 in 1:nsp    #For each species...
            comps= [(sum(A[i10,i11,:,:] .* freqmat)) for i11 in 1:nt]                           #Compute fitness reduction for each phenotype
            d1=Dict(pts[i10][i12]=>comps[i12] for i12 in 1:length(pts[i10]))                    #Create directory key-> phenotype, value- Fitness loss- use it to assign fitness to each individual
            prob1=[get(d1,popfreqs[i10][i13],0.0) for i13 in 1:length(popfreqs[i10])]
            pop[i10]=pop[i10][setdiff(1:size(pop[i10])[1],
                                wsample(1:length(prob1),prob1 .+ 0.001,changes[i10],replace=false)),:]
        end

        ##Put the population data into the result file
        for k in 1:nsp
            append!(res,
                    DataFrame(time=i1,sp=k,pheno=pheno_from_geno(pop[k],eff,tmins[k])))
        end

    end

    return res
end

#Set up parameter space

Random.seed!(11)

#Competition kernel
omega=0.5  #kernel width
t=10.0     #trait threshold

nsp= 3 #no. of species
Np=repeat([100],nsp) #Population sizes

nloci=5   #No. of loci
prob=fill(0.2,nsp)  # probability of alleles A across loci and populations    

#Effect size of each locus: try 4 different combinations
effs= ones(Float64,nloci,4)
effs[:,2] += [1,0,0,0,0]  #One locus slightly more important than others
effs[:,3] += [4,0,0,0,0]  #One locus contributes much more than others
effs[:,4] += [1,1,0,0,0]  #Two loci more important than other 3.
effs ./=    sum(effs,dims=1)

#Initial trait means
tmins=[-0.5,0.0,0.2]

tsteps=100

init=get_genpop(nsp,Np,nloci,prob)

################################################################################################
#Trial runs

omega=1.0
t=10.0

tmins=[-0.05,0.0,0.2]
eff=effs[:,1]

res=simIBM(init,eff[:,1],tmins,omega,t,tsteps)

anim1=@animate for i in 0:tsteps

    @df res[res.time .==i,:] StatsPlots.density(
        :pheno,
        group=:sp,
        xlabel="Trait value",
        ylabel="Frequency",
        size=(500,500)
    )

end
    
gif(anim1,fps=2)


meanres=res |>
        x-> groupby(x,[:time,:sp]) |>
        x-> combine(x,:pheno => mean,renamecols=false)

draw(data(meanres)*
    mapping(:time,
            :pheno =>L"Mean Population trait",
            color=:sp => x-> nonnumeric(x))*
    visual(Lines))
        
##################################################################################################
#Try the simulations for different initial trait means 

#Competition kernel
omega=0.5  #kernel width
t=0.5    #trait threshold

#Probabilities of alleles
prob=fill(0.4,nsp) 

t1s=collect(-0.9:0.1:-0.1)
t2s=collect(0.1:0.1:0.9)

tdat=collect(Iterators.product(t1s,t2s))

eff=fill(1.0,nloci)

resdat=DataFrame()

for i in 1:length(tdat)

    tmins=[tdat[i][1],tdat[i][2],0.0]

    res=simIBM(init,eff,tmins,omega,t,tsteps) |>
        x-> combine(groupby(x,[:sp,:time]), :pheno => mean, renamecols=false) 

    res.t1 .= tdat[i][1]
    res.t2 .= tdat[i][2]
    append!(resdat,res)

end

resdat1=resdat |>
filter(:sp => x-> x!= 3) |>
x-> combine(groupby(x,[:t1,:t2,:sp]), :pheno => diff => [:trdiff])

resdat1.tdiff = collect(Iterators.flatten(fill(collect(Iterators.flatten(fill(collect(1:100),2))),length(tdat))))


resdat1 |> 
x-> combine(groupby(x,[:t1,:t2,:tdiff]), :trdiff => diff => [:trdiff])


p1=data(effdat)*
    mapping(:time,
            :pheno=> L"Mean Population trait",
            color=:sp => x->nonnumeric(x),
            layout=:eff => x-> nonnumeric(x))*
    visual(Lines)

draw(p1)



