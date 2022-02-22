#=
Functions
=#

function InitPop(barcodes::Int,s_coef::Array{Float64,1},n0::Array{Float64,1})
	#initializes population
    lin1 = Lineage(1,s_coef[1],n0[1])
    population = fill(lin1,1)
	for i in 2:barcodes
		 lin2 = Lineage(i,s_coef[i],n0[i])
		 push!(population,lin2) 
	end
	population = Population(population)
return population
end

# First founder mutation
mutation0 = mutation(1,0.8)
# Initial Size
N_0 = 1
# Initial clone
clone0 = Clone([mutation0],N_0)
# Initial Population
Population0 = Population([clone0])
push!(Population0.clones,clone0)





function divide(clone_info,dt,B0,D0)

    mu_vec=clone_info.Mutations
    clone_size=clone_info.N
    fitness_effects=[x.Fitness for x in mu_vec ]

    F=sum(fitness_effects)
    B=B0*(1+F)
    pvar1=Poisson(clone_size*B*dt) 
    number_of_births=rand(pvar1,1)
    pvar2=Poisson(clone_size*D0*dt)
    number_of_deaths=rand(pvar2,1)
    clone_info.N=clone_info.N+number_of_births[1]+number_of_deaths[1]

    return clone_info
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
