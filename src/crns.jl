using Catalyst
using DifferentialEquations

export pure_immigration
export pure_death
export pure_birth
export birth_death
export birth_death_immigration

#=
BID Chemical Reaction Network Representation
=#

# CRN
pure_immigration = @reaction_network begin
    K_i, 0 --> A
end K_i

pure_death = @reaction_network begin
    K_d, A --> 0
end K_d

pure_birth = @reaction_network begin
    K_b, A --> 2A
end K_b

birth_death = @reaction_network begin
    K_b, A --> 2A
    K_d, A --> 0
end K_b K_d

birth_death_immigration = @reaction_network begin
    K_b, A --> 2A
    K_d, A --> 0
    K_i, 0 --> A
end K_b K_d K_i

#=
General Simple Models
=#

SimpleEquilibrium = @reaction_network begin
    K_f, A --> B
    K_b, B --> A
end K_f K_b
