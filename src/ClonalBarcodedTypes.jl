############### Type Hierarchy Top Element
abstract type Lineages end
abstract type Populations end

############### Population Types
mutable struct Population <: Populations
    lineage::Array{Lineages}
end

############### Lineage Types
# Requires a compatible reset_pop function for each Population type for inference

# Lineages without mutation
mutable struct PureBirth <: Lineages
    barcode::Int
    Kb::Float64
    N::Float64
end


