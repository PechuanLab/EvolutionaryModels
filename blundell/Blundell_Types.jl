abstract type Lineages end

struct Mutation <: Lineages
    ID::Int
    Fitness::Float64
end

mutable struct Clone <: Lineages
    Mutations::Array{Mutation}
    N::Int
end

mutable struct Population <: Lineages
    clones::Array{Clone}
end

