# EvolutionaryModels

<p align="center">
  <img width="700"  src="https://github.com/PechuanLab/EvolutionaryModels/blob/main/harveys.png.001.png">
</p>


## Installation

```julia
add https://github.com/PechuanLab/EvolutionaryModels.git
```

## Overview

This package contains a Julia implementation of useful population dynamics models I have played with. We also plan to implement some inference tools to fit the models data generated in experimental or sampling designs. 


## Simulating an Evolutionary Model

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
histogram(mut_histories[!,:Fitness])
histogram(mut_histories[!,:measuredVAF][.!(mut_histories[!,:measuredVAF] .== 0)], nbins=100) 
plot(mut_histories[2:end,:VAF],mut_histories[2:end,:Fitness],seriestype = :scatter)
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


## Fitting an Evolutionary Model to Data

Once an evolutionaray model is defined, we can make use of other Julia packages to perform inferences on data. As an example, we demonstrate the parameter inference on data coming from a barcoded cell line propagated for several passages on batch culture conditions.

```julia
# Libraries
using EvolutionaryModels 
using KissABC
using Distributions
using ClusterManagers
using DataFrames
using Setfield
using StatsBase
using CSV
using Distances
using Plots

######### Prepare data to be fitted
data = CSV.read("../data/Top2_Population_Counts.csv",DataFrame;drop=[1])
Tdata = EvolutionaryModels.TargetData(data)
# Timepoints sampled
TimePoint = Tdata[1]
# Target data to fit
tdata = Tdata[2]
# Number of barcodes
barcodes = Tdata[3]
# Initial condition of each barcoded Lineage
n0 = Tdata[4]

######### TimeSeries Parameters
parameters = CSV.read("../data/Complete_MasterSizes.csv",DataFrame;drop=[1])
Parametres = EvolutionaryModels.Params(parameters)
# Times between culutres
TimeCultures = Parametres[1]
# Cells transferred
Ntransferes = Parametres[2]
Passes = Parametres[3]
# Global Parameters
NPasses = Parametres[4]

# Define the prior
prior=Factored(Uniform(0,1),# s_RG0
              Uniform(0,1),  #s_RG1
              Uniform(0,1))  #S_Remainder

# Perform the inference
@time ressmc = smc(prior,EvolutionaryModels.costtest, nparticles=20, epstol=10,verbose=true,parallel=true)

```


## Authors

This package was developed by:

- Ximo Pechuan Jorge ([@jopejor](https://github.com/jopejor))
- Devra van Egeren ([@dvegeren](https://github.com/dvegeren))


