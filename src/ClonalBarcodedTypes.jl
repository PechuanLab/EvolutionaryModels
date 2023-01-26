############### Type Hierarchy Top Element
abstract type Lineages end
abstract type Populations end
abstract type DynamicsParameters end

############### DynamicsParameters Types
mutable struct PureBirthParam <: DynamicsParameters
    ParameterValues::Vector{Float64}
end

############### Population Types
mutable struct Population <: Populations
    lineage::Array{Lineages}
end

############### Lineage Types
# Lineages without mutation
mutable struct AsexualClone <: Lineages
    barcode::Int
    Parameters::DynamicsParameters
    N::Float64
end

# Mutating lineages
mutable struct MutClone <: Lineages
    barcode::Int
    ParametersParent::DynamicsParameters
    NParent::Float64
    ParametersMutant::DynamicsParameters
    NMutant::Float64
    τ::Float64
    mutates::Bool
end

############# Specific Types for a particular model

#=
Blundell model
=#

struct Mutation <: Lineages
    ID::Int
    Fitness::Float64
end

mutable struct BlundellClone <: Lineages
    Mutations::Array{Mutation}
    N::Int
end

