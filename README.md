# EvolutionaryModels

<p align="center">
  <img width="700"  src="https://github.roche.com/pechuanj/EvolutionaryModels/blob/master/harveys.png.001.png">
</p>


## Installation

```julia
add https://github.roche.com/pechuanj/EvolutionaryModels
```

## Overview

This package contains a Julia implementation of useful population dynamics models. We also plan to implement some inference tools to fit the models data generated in experimental or sampling designs. 

## Example 1

The following example recreates the results presented in Blundell et al. 2020.

```julia
using EvolutionaryModels

############## Evolutionary Simulation
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


# How many time epochs ?
lifespan/dt

# Main function
@time population,last_id,mut_histories = EvolutionaryModels.EvolutionaryDynamics(eq_population_size,dt,lifespan,μ,w,θ)

############# Plottinng the results
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

## Example 2

This package contains the companion code for the paper "Comparing exact and approximate solutions to a stochastic process model of clone size distribution". Mostly encodes the useful analytical expressions and provides wrappers to compare them with stochastic simulations generated using other Julia packages. Here are the lattice of limiting cases for the Birth-Death-Immigration process.

```julia
# Libraries
using BID
using DifferentialEquations

##################### Birth-Death Process Starting from n_0 = 1
K_b = (5.5)
K_d = (4.5)
N_0 = 1
t = 1
t_0 = 0

prob = DiscreteProblem(birth_death, [N_0], (0.0,t), (K_b,K_d))
jump_prob = JumpProblem(birth_death, prob, Direct())
solutions = BID.SamplesGillespie(jump_prob,SSAStepper(),100)
trajectories = BID.solutionstoDF(solutions)
BID.plotCRNGillespie(solutions)


# Let's compare
nmostres = 2000
FinalMostra = SampleFinal(nmostres,jump_prob,SSAStepper())
Nmin = minimum(FinalMostra[1])
Nmax = maximum(FinalMostra[1])
orange = Nmin:1:Nmax
τ = (K_b-K_d)*(t-t_0)
σ = (K_b-K_d)/(K_b+K_d)
DistroUs = P_BD_n01.(σ,τ,orange)
histogram(FinalMostra[1],normalize = true,bins=100)
#plot!(orange,Distro)
plot!(orange,DistroUs)


##################### Birth-Death Process Starting from n_0 = 20
K_b = (5.5)
K_d = (4.5)
N_0 = 20
t = 1
t_0 = 0

prob = DiscreteProblem(birth_death, [N_0], (0.0,t), (K_b,K_d))
jump_prob = JumpProblem(birth_death, prob, Direct())
solutions = SamplesGillespie(jump_prob,SSAStepper(),100)
trajectories = solutionstoDF(solutions)
plotCRNGillespie(solutions)


# Let's compare
nmostres = 2000
FinalMostra = SampleFinal(nmostres,jump_prob,SSAStepper())
Nmin = minimum(FinalMostra[1])
Nmax = maximum(FinalMostra[1])
orange = Nmin:1:Nmax
τ = (K_b-K_d)*(t-t_0)
σ = (K_b-K_d)/(K_b+K_d)
DistroUs = P_BD(σ,τ,orange,N_0)
histogram(FinalMostra[1],normalize = true,bins=100)
#plot!(orange,Distro)
plot!(orange,DistroUs[1:length(orange)])

##################### Full Birth-Death-Immigration Process Starting from n_0 = 2
K_b = (5.5)
K_d = (4.5)
K_i = (1)
N_0 = 2
t = 0.4
t_0 = 0

prob = DiscreteProblem(BID, [N_0], (0.0,t), (K_b,K_d,K_i))
jump_prob = JumpProblem(BID, prob, Direct())
solutions = SamplesGillespie(jump_prob,SSAStepper(),100)
trajectories = solutionstoDF(solutions)
plotCRNGillespie(solutions)


# Let's compare
nmostres = 2000
FinalMostra = SampleFinal(nmostres,jump_prob,SSAStepper())
Nmin = minimum(FinalMostra[1])
Nmax = maximum(FinalMostra[1])
orange = Nmin:1:Nmax
τ = (K_b-K_d)*(t-t_0)
σ = (K_b-K_d)/(K_b+K_d)
m = K_i/(K_b-K_d)
DistroUs = P_BID(orange,N_0,t,K_d,K_b,K_i,t_0)
histogram(FinalMostra[1],normalize = true,bins=20)
#plot!(orange,Distro)
plot!(orange,DistroUs[1:length(orange)])


##################### Pure Immigration
K_i = (2.8)
N_0 = 0
t = 18

prob = DiscreteProblem(pure_immigration, [N_0], (0.0,t), K_i)
jump_prob = JumpProblem(pure_immigration, prob, Direct())
solutions = SamplesGillespie(jump_prob,SSAStepper(),100)
trajectories = solutionstoDF(solutions)
plotCRNGillespie(solutions)



# Let's compare
nmostres = 2000
FinalMostra = SampleFinal(nmostres,jump_prob,SSAStepper())
Nmin = minimum(FinalMostra[1])
Nmax = maximum(FinalMostra[1])
orange = Nmin:1:Nmax
Distro = PureImmigration.(K_i,t,orange,N_0)
DistroUs = P_BID.(orange,N_0,t,0.000000001,0.00001,K_i,0)
histogram(FinalMostra[1],normalize = true)
#plot!(orange,Distro)
plot!(orange,DistroUs)


```

