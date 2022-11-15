#= 
External Libraries 
=#
using Catalyst
using DifferentialEquations
using Distributions
using StatsBase


#= 
Accessory Functions:
=#

"""
Extracts and array of fields from a Population

# Examples

```julia-repl
julia> unwrapper(poblacio,"N")
```
...
# Arguments
...

"""
function unwrapper(poblacio::Population, fieldname)
    x = fill(0.0,length(poblacio.lineage))
    for i=1:length(poblacio.lineage)
        x[i]=getfield(poblacio.lineage[i],Symbol(fieldname)) 
    end
    return x
end

"""
Vector of vectors
# Examples

```julia-repl
julia> VectorofVecs(rand(d,10))
```

...
...

"""
function VectorofVecs(s::Vector{Float64})
	sels = []
	for i in 1:length(s)
		push!(sels,[s[i]])
	end
return sels
end

#= 
Population Dynamics Functions:
=#

#################################################### Population Initialization #####################################################

"""
Initializes an Asexual Population with N Barcoded Lineages

# Examples

```julia-repl
        Barcodes = collect(1:10)
        DynParamsVec = PureBirthParam.(VectorofVecs(rand(d,10)))
        n0 = 100*collect(1:10)
julia> InitBarcodedPop(AsexualClone,tuple.(Barcodes,DynParamsVec,n0))
```

...
# Arguments
- `LinType::DataType`: Lineage Type
- `TypeParamsVector::Vector{Tuple}`: Associated vector with tuples of compatible type paramters.
...

"""
function InitBarcodedPop(LinType::DataType,TypeParamsVector)
	# First initialize the population
    lin1 = LinType(TypeParamsVector[1]...)
    population = fill(lin1,1)
    for i in 2:length(TypeParamsVector)
		 lin2 = LinType(TypeParamsVector[i]...)
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

"""
Progates a clonal lineage according to the exact solution of the Yule Process

# Examples

```julia-repl
julia> CloneYuleGrowth(AsexualClone(1,PureBirthParam([0.5]),100),100.0)

```

...
# Arguments
- `clone::Lineages`: Clonal lineage to act on 
...

"""
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


"""
Progates a clonal lineage according to an approximation of the Yule Process

# Examples

```julia-repl
julia> CloneYuleGrowthNormal.(population.lineage ,20.0)

```

...
# Arguments
- `clone::Lineages`: Clonal lineage to act on 
-
...

"""
function CloneYuleGrowthNormal(clone::Lineages,t::Float64)
	# Safeguard
	if t == 0
		 clone = clone
	end
	# Parameters
	ReactionParams = clone.Parameters.ParameterValues
	p = exp(-ReactionParams[1] * t)
	μ = clone.N*(p)^-1
	σ = (clone.N*(1-p)*p^2)^(1/2)
	# Sample from the negative binomial
	d = Normal(μ,σ)
	x = floor(rand(d,1)[1])
	# Reflecting boundary for safeguard
	if x[1] < 0
	 x[1] = 0
	end
	clone.N = x[1]
return clone
end

#################################################### Transfer Dynamics ###############################################################

"""
Progates a clonal lineage according to a compatible single type stochastic 
model

# Examples

```julia-repl
julia> PoissonTransfer(clone,10,1000)

```

...
# Arguments
- `CRN`: A Chemical Reaction Network
- `interval::Tuple`: Time interval of growth
- `clone::Lineages`: Clonal lineage to act on 
...

"""

function PoissonTransfer(clone::Lineages,Ntransfer::Float64,Ntot::Float64)
	# This works if clones are independent
	λ = clone.N/Ntot*Ntransfer
	d = Poisson(λ)
	N = rand(d)
	clone.N = N
	return clone
end


"""
Progates a clonal lineage according to a compatible single type stochastic 
model

# Examples

```julia-repl
julia> MultinomialTransfer(poblacio,10)

```

...
# Arguments
...

"""
function MultinomialTransfer(poblacio::Population,Ntransfer::Int)
	#function body
	N_vec = unwrapper(poblacio,"N")
	barcodes = length(poblacio.lineage)
	total=sum(N_vec)
	Linage_freq=N_vec/total
	d=Categorical(Linage_freq)
	R = rand(d,Ntransfer)
	NewPopN = counts(R,barcodes)
	for i in 1:length(poblacio.lineage)
	    poblacio.lineage[i].N  = NewPopN[i] 
	end
	return poblacio
end


#################################################### Population Dynamics ###############################################################


"""
Simulates Simple Asexual Population being transferred in batch culture.

# Examples

```julia-repl
Barcodes = collect(1:10)
TypeParamsVector = rand(d,10)
n0 = 100*collect(1:10)
LinType = AsexualClone
NPasses = 10
CulturePass = CloneYuleGrowthNormal
TimeCulture =  5.0*collect(1:NPasses)
Ntransfer = 10000*collect(1:NPasses)
TransferFunction = PoissonTransfer
@time BatchCulture(LinType,Barcodes,n0,TypeParamsVector,NPasses,CulturePass,TimeCulture,TransferFunction,Ntransfer)
```

...
# Arguments
...

"""
function BatchCulture(LinType::DataType,Barcodes,n0,TypeParamsVector,NPasses::Int,CulturePass::Function,TimeCulture,TransferFunction::Function,Ntransfer)
	# Prepare the fitness
	DynParamsVec = PureBirthParam.(VectorofVecs(TypeParamsVector))
    # Initial Population
	poblacio =  InitBarcodedPop(LinType,tuple.(Barcodes,DynParamsVec,n0))	
	# Things we will record
	BarCodeMat = zeros(Float64,NPasses,length(poblacio.lineage))
	# Culture cycles
	for i in 1:NPasses
		CulturePass.(poblacio.lineage ,TimeCulture[i])
	    Nvec = unwrapper(poblacio,"N")
	    BarCodeMat[i,:] = Nvec
	    TransferFunction.(poblacio.lineage,Ntransfer[i],sum(Nvec)))
	end
	return BarCodeMat
end
