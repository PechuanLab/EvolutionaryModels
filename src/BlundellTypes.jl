export Mutation
export Clone

struct Mutation <: Lineages
    ID::Int
    Fitness::Float64
end

mutable struct Clone <: Lineages
    Mutations::Array{Mutation}
    N::Int
end



