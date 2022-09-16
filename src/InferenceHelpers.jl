# Dependencies
using KissABC
using Distributions
using ClusterManagers
using DataFrames
using Setfield
using StatsBase
using CSV
using Distances
using Plots

##################################### Cost Functions (need to work on how to make them universal)
function cost((s_RG0,s_RG1,s_RG2,s_RG3,s_RG4,s_RG5,s_RG6,s_RG7,s_RG8,s_RG9,s_RG10,S_Remainder))
    x = sim((s_RG0,s_RG1,s_RG2,s_RG3,s_RG4,s_RG5,s_RG6,s_RG7,s_RG8,s_RG9,s_RG10,S_Remainder))
    y = tdata
    x = Matrix(select!(x, Not(:TimePoint)))
    d = (sum(colwise(Euclidean(), x, y)))^(1/2)
    d
end

function cost((s_RG1,S_Remainder))
    x = sim((s_RG1,S_Remainder))
    y = tdata
    x = Matrix(select!(x, Not(:TimePoint)))
    d = (sum(colwise(Euclidean(), x, y)))^(1/2)
    d
end
# Cost Function
function cost((s_RG0,s_RG1,S_Remainder,r_0,r_1,r_2,K_1,K_2,K_3))
    x = sim((s_RG0,s_RG1,S_Remainder,r_0,r_1,r_2,K_1,K_2,K_3))
    y = tdata
    x = Matrix(select!(x, Not(:TimePoint)))
    d = (sum(colwise(Euclidean(), x, y)))^(1/2)
    d
end
# Cost Function
function cost((s_RG0,s_RG1,S_Remainder,r_0,r_1,r_2,K_1,K_2,K_3,M1,M2,M3))
    x = sim((s_RG0,s_RG1,S_Remainder,r_0,r_1,r_2,K_1,K_2,K_3,M1,M2,M3))
    y = tdata
    x = Matrix(select!(x, Not(:TimePoint)))
    d = (sum(colwise(Euclidean(), x, y)))^(1/2)
    d
end