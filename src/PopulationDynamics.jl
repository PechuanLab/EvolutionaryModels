#= 
External Libraries 
=#
using Catalyst
using DifferentialEquations
using Distributions
using FiniteStateProjection

#= 
Population Dynamics Functions:
=#

#################################################### Population Initialization #####################################################

"""
Initializes an Asexual Population with N Barcoded Lineages

# Examples

```julia-repl
julia> InitBarcodedPop(AsexualClone,(1,PureBirthParam([0.5]),100),2)
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

#################################################### Lineage Dynamics ###############################################################

"""
Progates a clonal lineage according to a compatible single type stochastic 
model

# Examples

```julia-repl
julia> CloneStochGrowth(pure_birth,(0.0,10),AsexualClone(1,PureBirthParam([0.5]),100))

```

...
# Arguments
- `CRN`: A Chemical Reaction Network
- `interval::Tuple`: Time interval of growth
- `clone::Lineages`: Clonal lineage to act on 
...

"""
function CloneStochGrowth(CRN,interval::Tuple,clone::Lineages)
	# Parameters
	N_0 = [clone.N]
	ReactionParams = clone.Parameters.ParameterValues
	# Gillespie Samples according to the reaction network
	prob = DiscreteProblem(CRN, N_0,interval, ReactionParams)
	jump_prob = JumpProblem(pure_birth, prob, Direct())
	solution = solve(jump_prob,SSAStepper())
	clone.N = solution.u[end][1]
return clone
end


"""
Progates a clonal lineage according to a compatible single type stochastic 
model with the Langevin Chemical Equation

# Examples

```julia-repl
julia> CloneStochGrowth(pure_birth,(0.0,10),AsexualClone(1,PureBirthParam([0.5]),100))

```

...
# Arguments
- `CRN`: A Chemical Reaction Network
- `interval::Tuple`: Time interval of growth
- `clone::Lineages`: Clonal lineage to act on 
...

"""
function CloneLangevinGrowth(CRN,interval::Tuple,clone::Lineages)
	# Parameters
	N_0 = [clone.N]
	ReactionParams = clone.Parameters.ParameterValues
	# Gillespie Samples according to the reaction network
	sprob = SDEProblem(CRN, N_0,interval, ReactionParams)
	ssol  = solve(sprob, EM(), dt=.01)
	clone.N = solution.u[end][1]
return clone
end

"""
Progates a clonal lineage according to a compatible single type stochastic 
model in the ode limit

# Examples

```julia-repl
julia> CloneODEGrowth(pure_birth,(0.0,10),AsexualClone(1,PureBirthParam([0.5]),100))

```

...
# Arguments
- `CRN`: A Chemical Reaction Network
- `interval::Tuple`: Time interval of growth
- `clone::Lineages`: Clonal lineage to act on 
...

"""
function CloneODEGrowth(CRN,interval::Tuple,clone::Lineages)
	# Parameters
	N_0 = [clone.N]
	ReactionParams = clone.Parameters.ParameterValues
	# Gillespie Samples according to the reaction network
	oprob = ODEProblem(CRN, N_0, interval,ReactionParams)
	osol  = solve(oprob, Tsit5())
	clone.N = osol.u[end][1]
return clone
end


function CloneYuleGrowth(clone::Lineages,t::Float64)
	# Parameters
	r = clone.N
	ReactionParams = clone.Parameters.ParameterValues
	p = exp(-ReactionParams[1] * t)
	# Sample from the negative binomial
	d = NegativeBinomial(r,p)
	clone.N = rand(d,1)[1]
return clone
end

function CloneYuleGrowthNormal(clone::Lineages,t::Float64)
	# Parameters
	ReactionParams = clone.Parameters.ParameterValues
	p = exp(-ReactionParams[1] * t)
	μ = clone.N*(p)^-1
	σ = (clone.N*(1-p)*p^2)^(1/2)
	# Sample from the negative binomial
	d = Normal(μ,σ)
	clone.N = floor(rand(d,1)[1])
return clone
end

d1 = NegativeBinomial(r,p)
d2 = Normal(μ,σ)
vec1 = rand(d1,1000)
vec2 = rand(d1,1000)

using EarthMoversDistance

histogram1 = rand(d1,1000)
histogram2 = rand(d1,1000)
earthmovers(histogram1, histogram2, (x, y) -> abs(x - y)) # custom ground distance function


