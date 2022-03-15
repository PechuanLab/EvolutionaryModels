#=
Packages
=#
 using Distributions
 using Random
 using DataFrames

#=
Functions
=#

"""
Unwraps a population on its fields
...
# Arguments
- `population::Population`: Population.
- `fieldname::String`: Field to unwrap.
...
"""
function unwrapper_pop(population::Population, fieldname::String)
    x = []
    for i=1:length(population.clones)
        push!(x,getfield(population.clones[i],Symbol(fieldname)) )
    end
    return x
end

"""
Unwraps a mutatin on its fields
...
# Arguments
- `mutations::Mutations`: Population.
- `fieldname::String`: Field to unwrap.
...
"""
function unwrapper_mut(mutations::Array{Mutation}, fieldname::String)
    x = []
    for i=1:length(mutations)
        push!(x,getfield(mutations[i],Symbol(fieldname)) )
    end
    return x
end


"""
Divides a clonal lineage according to its fitness
...
# Arguments
- `clone::Clone`: Clonal lineage selected for division.
- `dt::Float64=0.1`: Time interval step.
- `B0::Float64=0.2`: Initial Birth rate.
- `D0::Float64=0.1`: Initial Death rate.
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
    clone.N=clone.N+number_of_births[1]-number_of_deaths[1]

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
function Mutate(clone::Clone,dt::Float64,μ::Float64,last_id::Int64,w::Array{Float64},θ::Array{Float64})
    
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

"""
Models sequencing and other sampling aspects
...
# Arguments
- `population::populationulation: plug in populationulation 
...
"""
function Sequence(population::Population,mean_depth::Float64, last_id::Int64)

    # Initialize DF
    CloneDF0 = DataFrame()
    CloneDF0[!,:ID] = 1:last_id
    CloneDF0[!,:Fitness] = zeros(last_id)
    CloneDF0[!,:N] = zeros(last_id)
    
    # Unwrapp
    clone_sizes = unwrapper_pop(population,"N")
    mutations = unwrapper_pop(population,"Mutations")


    #For each clone unwrap
    # Grow the seed with mutatioons
    for i = 1:length(mutations)

        CloneDFnext = DataFrame()
        clone_mutations = mutations[i]
        CloneDFnext[!,:ID] = unwrapper_mut(clone_mutations,"ID")
        CloneDFnext[!,:Fitness] = unwrapper_mut(clone_mutations,"Fitness")
        CloneDFnext[!,:N] = clone_sizes[i]*ones(length(clone_mutations))

        CloneDF0[!,:N][CloneDFnext[!,:ID]] = CloneDF0[!,:N][CloneDFnext[!,:ID]] + CloneDFnext[!,:N]
        CloneDF0[!,:Fitness][CloneDFnext[!,:ID]] = CloneDFnext[!,:Fitness]

    end

    # Let's add VAF (Variant Allele Frequency)
    CloneDF0[!,:VAF] = CloneDF0[!,:N]*1/2*1/PopSize(population)
    # Let's add the expected reads
    depth = round(Int64,mean_depth*exp.(rand(Normal(0,0.3), 1))[1])
    expected_reads = depth*CloneDF0[!,:VAF]
    reads_array =rand.(Poisson.(expected_reads),1)
    CloneDF0[!,:Reads] =  [reads_array[i][1] for i in 1:length(reads_array)]
    CloneDF0[!,:measuredVAF] = CloneDF0[!,:Reads]/depth

return CloneDF0

end  


"""
Population dynamics sampling from a Poission distribution given a time interval
...
# Arguments
- `population::populationulation: plug in populationulation 
...
"""

function PoissonDynamics(t,dt,B0,D0,population::Population,last_id::Int64,μ,w,θ)
          
          # Time update
            t = t+dt
            pop_size = PopSize(population) 
            new_clones0 = []

            # Main loop
            for i in 1:length(population.clones)
                # Divide
                Divide(population.clones[i],dt,B0,D0)
                # Avoid negative numbers
                if population.clones[i].N < 1
                   population.clones[i].N = 0
                end
                # Mutate
                new_clones,last_id = Mutate(population.clones[i],dt,dt,μ,last_id,w,θ)
                if isempty(new_clones) 
                    continue
                end
                for i = 1:length(new_clones)
                    #push !
                    push!(new_clones0,new_clones[i])
                end
            end
            #Clean-up
            if isempty(new_clones0)        
            else
                for i in 1:length(new_clones0)
                push!(population.clones,new_clones0[i])
                end
            end

            # Clean-up
            clone_sizes = unwrapper_pop(population,"N")
            # Remove a lineage if has negative or zero cell numbers
            deleteat!(population.clones, findall(x->x<1,clone_sizes))
return population, last_id
end


#=
Main Simulation for only one patient:

# Popuation Dynamics Parameters
    eq_population_size = 10^5 #number of stem cells (HSCs)
    lifespan = 500 #measured in cell divisions (e.g. 100 = 100 cell divisions)
    dt=0.1 # measured in units of the overall cell division rate
    μ=3*10^(-6) #mutation rate
    lam=5 # HSC divisions per year (symmetric and asymmetric)
    last_id=1
    mean_depth=3000 #sequencing depth

# Fitness Landscape
    sn = 0 ; sb = 0.05 ; sd = 0 ; w = [sn,sb,sd]
    pn = 1/3 ; pb = 2/3 ; pd = 0 ; θ = [pn,pb,pd]
=#
function EvolutionaryDynamics(eq_population_size::Int64,dt::Float64,lifespan::Int,μ::Float64,w::Array{Float64},θ::Array{Float64})


    # Initialize the population, ab initio we can consider this the HSC compartment size
    last_id = 1
    mutation0 = Mutation(last_id,0.0) # Reference fitness founder
    N_0 = eq_population_size # Ne for the initial condition
    # Initial clone
    clone0 = Clone([mutation0],N_0)
    # Initial Population
    population = Population([clone0])
    # t
    t = 0
    mean_depth = 300.0
    # Initial growth phase
    pop_size = PopSize(population)
    
    Ndts = Int(floor(lifespan/dt))

    #Initialize mut_histories 
    mut_histories = Sequence(population,mean_depth,last_id)
    mut_histories[!,:T] = [0]
    # Remove the while as it is dangerous
    for i in 1:Ndts
        if pop_size < eq_population_size
            B0 = 1.0
            D0 = 0.0
            population, last_id = PoissonDynamics(t,dt,B0,D0,population,last_id,μ,w,θ)
            pop_size = PopSize(population)
            # mut_histories 
            mut_histories0 = Sequence(population,mean_depth,last_id)
            mut_histories0[!,:T] = i*ones(nrow(mut_histories0))  
            mut_histories = vcat(mut_histories,mut_histories0)    
        else 
            B0 = 0.2
            D0 = 0.2
            population, last_id = PoissonDynamics(t,dt,B0,D0,population,last_id,μ,w,θ)
            pop_size = PopSize(population)
            # mut_histories 
            mut_histories0 = Sequence(population,mean_depth,last_id)
            mut_histories0[!,:T] = i*ones(nrow(mut_histories0))  
            mut_histories = vcat(mut_histories,mut_histories0)  
        end

    end

return population, last_id, mut_histories
end
