#= 
External Libraries 
=#
using Catalyst

#= 
Population Dynamics Functions:
=#

################### Population Initialization

"""
Initializes an Asexual Population with N Barcoded Lineages

# Examples

```julia-repl
julia> InitBarcodedPop(PureBirth,(1,0.5,100),2)
Population(Lineages[PureBirth(1, 0.5, 100.0), PureBirth(1, 0.5, 100.0)])
```

...
# Arguments
- `LinType::DataType`: Lineage Type
- `TypeParams::Tuple`: Assoaciated compatible type paramters.
- `Barcodes::Int`: Total Number of Barcoded Lineages 
...

"""
function InitBarcodedPop(LinType::DataType,TypeParams::Tuple,Barcodes::Int)
	# First initialize the population
    lin1 = LinType(TypeParams...)
    population = fill(lin1,1)
    for i in 2:Barcodes
		 lin2 = LinType(TypeParams...)
		 push!(population,lin2) 
	end
	population = Population(population)
return population
end

################### Lineage Dynamics

function StochGrowth(CRN,ReactionParams::Tuple,)
	# CRNT full

return Lineage





function CulturePasss(poblacio::Population,TimeCulture)
	if TimeCulture == 0
	    poblacio = poblacio
	else
		# Exponential assumption
		fitness_vec = unwrapper(poblacio,"fitness")
		N_vec = unwrapper(poblacio,"N")
		mass_action = N_vec.*exp.(TimeCulture*fitness_vec)
		sdev_action = ((N_vec.*exp.(2*TimeCulture*fitness_vec)).*(ones(length(fitness_vec))-exp.(-TimeCulture*fitness_vec))).^(1/2)
		for i in 1:length(poblacio.clones)
			mu = mass_action[i]
			sigma = sdev_action[i]
			d = Normal(mu,sigma)
			x = rand(d, 1)
			if x[1] < 0
				x[1] = 0
			end
	    	poblacio.clones[i].N = x[1]
		end
	end
	return poblacio
end
