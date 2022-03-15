# EvolutionaryModels

## Overview

This package contains a Julia implementation of useful population dynamics models. We also plan to implement some inference tools to fit the models data generated in experimental or sampling designs. 

## Example

The following example recreates the results presented in Blundell et al. 2020.

```julia
using EvolutionaryModels

# Popuation Dynamics Parameters
eq_population_size = 10^5 #number of stem cells (HSCs)
lifespan = 200 #measured in cell divisions (e.g. 100 = 100 cell divisions)
dt=0.1 # measured in units of the overall cell division rate
μ=3*10^(-6) #mutation rate
lam=5 # HSC divisions per year (symmetric and asymmetric)
last_id=1
mean_depth=3000.0 #sequencing depth

# Fitness Landscape
sn = 0 ; sb = 0.05 ; sd = 0 ; sk = 10*sb ; w = [sn,sb,sd,sk]
pn = 1/4 ; pb = 2/3 ; pd = 0 ; pk = 1-1/4-2/3 ; θ = [pn,pb,pd,pk]


lifespan/dt

@time population,last_id,mut_histories = EvolutionaryModels.EvolutionaryDynamics(eq_population_size,dt,lifespan,μ,w,θ)

using Plots

histogram(mut_histories[!,:VAF],nbins=100)
histogram(df[!,:Fitness])
histogram(df[!,:measuredVAF][.!(df[!,:measuredVAF] .== 0)], nbins=100) 
plot(df[2:end,:VAF],df[2:end,:Fitness],seriestype = :scatter)
plot(mut_histories[!,:T],mut_histories[!,:VAF],seriestype = :scatter)

using StatsPlots, RDatasets
gr()
mut_histories1 = filter(:ID => !=(1),mut_histories)
mut_histories1 = filter(:VAF => !=(0),mut_histories1)

mut_histories1[!,:VAF] = log.(mut_histories1[!,:VAF])
@df mut_histories1 scatter(
    :T,
    :VAF,
    group = :ID
)

```
