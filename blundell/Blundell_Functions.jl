# First founder mutation
mutation0 = Mutation(1,0.8)
mutation1 = Mutation(2,0.4)

# Initial Size
N_0 = 1
# Initial clone
clone = Clone([mutation0,mutation1],N_0)
# Initial Population
Population0 = Population([clone0])
push!(Population0.clones,clone0)
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
function divide(clone::Clone,dt::Float64,B0::Float64,D0::Float64)

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


function mutate(clone_info,dt,u,last_id,DFE)
    muta=clone_info.Mutations
    csize=clone_info.N
    pvar3=Poisson(csize*u*dt)
    num_muta=rand(pvar3,1)

    for i in range(start=1,stop=num_muta,step=1)
        last_id=last_id+1
        new_fit=mutation_fitness(DFE)
        push!(clone_info.Mutations,Mutation(last_id,new_fit))
    end
    
    return (clone_info.Mutations,last_id)
end 
