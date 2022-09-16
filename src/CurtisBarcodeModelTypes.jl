############### Top Types
abstract type Lineages end
abstract type Populations end

############### Lineage Types
# Lineage without mutation
mutable struct Lineage <: Lineages
    barcode::Int
    fitness::Float64
    N::Float64
end

# Lineage with mutation at a particular passage
mutable struct LineageMutPass <: Lineages
    barcode::Int
    fitness::Float64
    N::Float64
    PassMut::Int
    fitnessMut::Float64
    NMut::Float64
end
# Lineage with a mutation at a particular time
mutable struct LineageMutTime <: Lineages
    barcode::Int
    fitness::Float64
    N::Float64
    τ::Float64
    fitnessMut::Float64
    NMut::Float64
    mutates::Bool
end
# Lineage with timedepent fitness
mutable struct LineageW <: Lineages
	barcode::Int
    w_0::Float64
    N::Float64
    K::Float64
    r::Float64
end

mutable struct LineageW2 <: Lineages
    barcode::Int
    w_0::Float64
    N::Float64
    K::Float64
    r::Float64
    M::Float64
end
# Population type
mutable struct Population <: Populations
    llinatges::Array{Lineages}
end