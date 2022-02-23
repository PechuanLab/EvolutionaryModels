last_id = 1

# First founder mutation
mutation0 = Mutation(last_id,0.8)
last_id = last_id + 1
mutation1 = Mutation(last_id,0.4)

# Initial Size
N_0 = 1
# Initial clone
clone = Clone([mutation0,mutation1],N_0)
# Initial Population
population = Population([clone])
push!(population.clones,clone)
divide(clone)

#=
Packages
=#
 using Distributions
 using Random

#=
Functions
=#

"""
Divides a clonal lineage according to its fitness
...
# Arguments
- `clone::Clone`: Clonal lineage selected for division.
- `dt::Float64=0.1`: Time interval step.
- `B0::Float64=0.2`: Initial Birth rate.
- `D0::Float64=0.2`: Initial Death rate.
...
"""
function Divide(clone::Clone,dt::Float64,B0::Float64,D0::Float64)

    # Get attributes
    mu_vec=clone.Mutations
    clone_size=clone.N
    fitness_effects=[x.Fitness for x in mu_vec ]
    
    # Update clone size
    F=sum(fitness_effects) #Additive, it could be multiplicative
    B=B0*(1+F)             #Here, it becomes multplicative for some reason
    pvar1=Poisson(clone_size*B*dt) 
    number_of_births=rand(pvar1,1)
    pvar2=Poisson(clone_size*D0*dt)
    number_of_deaths=rand(pvar2,1)
    clone.N=clone.N+number_of_births[1]+number_of_deaths[1]

    return clone
end 


"""
Samples the fitness effect for a new mutation
...
# Arguments
- `w::Array{Float64}`: Fitness Landscape, 
    #neutral benefit, benefidial, deleterious# 
    sn = 0 ; sb = 0.05 ; sd = 0 ; w = [sn,sb,sd]
- `θ::Array{Float64}`: Probability of sampling a mutation of a type Landscape, #neutral benefit, benefidial, deleterious# 
    pn = 1/3 ; pb = 2/3 ; pd = 0 ; θ = [pn,pb,pd]
...
"""
function MutationFitness(w::Array{Float64},θ::Array{Float64}) 

    # Sample mutation fitnesss
    FitnessLandscape = w
    DFE = Categorical(θ) 
    FitnessEffect = FitnessLandscape[rand(DFE,1)][1]

    return FitnessEffect
end 

"""
Mutates a clonal lineage generating a new one
...
# Arguments
- `clone::Clone`: Clonal lineage selected for division.
- `dt::Float64=0.1`: Time interval step.
- `μ::Float64=0.2`: Mutation rate
- `last_id::Int`: Global bookeeping of mutation id
- `w::Array{Float64}`: Fitness Landscape, w = [sn,sb,sd],sn = 0 #neutral benefit
    sb = 0.05 #beneficial , sd = 0  #deleterious
- `θ::Array{Float64}`: Probability of sampling a mutation of a type Landscape, θ = [pn,pb,pd],pn = 1/3 #neutral benefit
    pb = 2/3 #beneficial , pd = 0  #deleterious
...
"""
function Mutate(clone::Clone,dt::Float64,μ::Float64,last_id::Int,w::Array{Float64},θ::Array{Float64})
    
    # Get the number of mutations that happen on dt interval
    muta=clone.Mutations
    pvar3=Poisson(clone.N*μ*dt)
    num_muta=rand(pvar3,1)[1]

    # Add each of these mutations and independet
    new_clones = []
    for i in 1:num_muta
        last_id=last_id+1
        new_fit=MutationFitness(w,θ)
        clone_mut = deepcopy(clone)
        push!(clone_mut.Mutations,Mutation(last_id,new_fit))
        clone_mut.N = 1
        push!(new_clones, clone_mut)
    end
    new_clones = convert(Array{Clone},new_clones)  
    return new_clones,last_id
end 

"""
Returns the total population size
...
# Arguments
- `population::Population`: Current population
...
"""
function PopSize(population::Population)

    pop_size = sum([x.N for x in population.clones])
    
    return(pop_size)

end

#for i in new_clones push!(population.clones,i) end

