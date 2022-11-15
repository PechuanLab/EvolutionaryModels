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

