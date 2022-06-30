# Dependencies
using BigCombinatorics
using Distributions
using CSV
using DataFrames
# using Catalyst
using DifferentialEquations
using Plots
using DSP

#=
 Analytical solutions for the BID family Proability distribution of having n individuals at time t
=#

"""
Pure immigration time dependent probability distribution
0 --> A
...
# Arguments
- `K_i::Float64`: Immigration kinetic constant integer numbers.
- `t::Float64`: Time from t0 = 0
- `n::Int`: Number of individuals of type A
- `n_0::Int`: Initial number of individuals of type A
...
"""
function PureImmigration(K_i::Float64, t::Float64, n::Int, n_0::Int)
    # PureBirth Immigration process poisson 
    # t: Time
    # n_0: Starting Population
    # n: final Population
    # K_i: Immigration rate
    p = (K_i * t)^n / factorial(big(n)) * exp(-K_i * t) + n_0
    return p
end

function PureDeath(K_d, t, N, N_0)
    # Pure Death Process Solution: binomial
    # t: Time in days
    # n: Starting Population
    # k: final Population
    # K_b: Birth rate
    p = binomial(N_0, N) * exp(-K_d * N * t) * (1 - exp(-K_d * t))^(N_0 - N)
    return p
end


function PureBirth(K_b, t, N, N_0)
    # PureBirth Process negative binomial
    # t: Time in days
    # n: Starting Population
    # k: final Population
    # K_b: Birth rate
    p = binomial(N - 1, N - N_0) * exp(-K_b * N_0 * t) * (1 - exp(-K_b * t))^(N - N_0)
    return p
end


function P_BD_n01(σ, τ, n)
    # Probability distribution for the BID with m = 0 and n_0 =1
    # Define the quantities of section III.C
    a = ((1 - σ) * (1 - exp(-τ))) / ((1 + σ) - (1 - σ) * exp(-τ))
    b = (4 * σ^2 * exp(-τ)) / ((1 + σ) - (1 - σ) * exp(-τ))^2
    c = ((1 + σ) * (1 - exp(-τ))) / ((1 + σ) - (1 - σ) * exp(-τ))
    # Probability
    p = δ(0, n) * a + (!δ(0, n)) * b * c^(n - 1)
    return p
end

function P_BD(σ, τ, n, n_0)
    # Convolve to get exact Birth-Death process 
    q = P_BD_n01.(σ, τ, n)
    p = conv(q, q)
    for i = 2:n_0
        #convolve with itself
        p = conv(q, p)
    end
    return p
end

function P_BID_n0(n, t, K_d, K_b, K_i, t_0)
    #=
    Birth Death Immigration for n0 
    Some default paramters
    n=10
    n_0=0
    t=100
    K_d=0
    K_b=0
    K_i=0.1
    t_0=0
    =#
    #Non-dimensional form
    τ = (K_b - K_d) * (t - t_0)
    σ = (K_b - K_d) / (K_b + K_d)
    m = K_i / (K_b - K_d)
    # collect on smaller terms
    c = (2 * m * σ) / (1 + σ)
    a = ((1 + σ - (1 - σ) * exp(-τ)) / (2 * σ * exp(-τ)))
    b = ((1 + σ) * (1 - exp(-τ))) / (1 + σ - (1 - σ) * exp(-τ))
    # Probability n_0 = 0  with immigration
    p = a^(-c) * 1 / factorial(big(n)) * Poch(c, n) * b^n
    return p
end



function P_BID(n, n_0, t, K_d, K_b, K_i, t_0)
    # Full BID Solution
    #Non-dimensional form
    τ = (K_b - K_d) * (t - t_0)
    σ = (K_b - K_d) / (K_b + K_d)
    m = K_i / (K_b - K_d)
    # sum
    q = P_BID_n0.(n, t, K_d, K_b, K_i, t_0)
    q = convert.(Float64, q)
    p = P_BD(σ, τ, n, n_0)
    p1 = conv(p, q)
    return p1
end


# Solution
function Simple_Equilibrium(K_f, K_b, N, N_0)
    # Pure Death Process Solution: binomial
    # t: Time in days
    # n: Starting Population
    # k: final Population
    # K_b: Birth rate
    p = binomial(N, N_0) * (K_f / (K_f + K_b))^N_0 * (K_b / (K_f + K_b))^(N - N_0)
    return p
end
