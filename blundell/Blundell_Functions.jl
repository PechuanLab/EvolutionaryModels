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