# Export
export Population
# Highlevel types for any evolutionary simulation
abstract type Lineages end
abstract type Populations end

mutable struct Population <: Populations
    clones::Array{Lineages}
end