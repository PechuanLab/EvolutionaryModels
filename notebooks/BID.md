## Analytical Solutions to the BID 

This section contains a piece of code part of a larger set required to reproduce the results of the paper "Comparing exact and approximate solutions to a stochastic process model of clone size distribution". Mostly encodes the useful analytical expressions and provides wrappers to compare them with stochastic simulations generated using other Julia packages. Here are the lattice of limiting cases for the Birth-Death-Immigration process as depicted in the paper.

```julia
# Libraries
using EvolutionaryModels
using DifferentialEquations
using Plots


##################### Pure Immigration
K_i = (2.8)
N_0 = 0
t = 18

prob = DiscreteProblem(pure_immigration, [N_0], (0.0,t), K_i)
jump_prob = JumpProblem(pure_immigration, prob, Direct())
solutions = EvolutionaryModels.SamplesGillespie(jump_prob,SSAStepper(),100)
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


##################### Birth-Death Process Starting from n_0 = 1
K_b = (5.5)
K_d = (4.5)
N_0 = 1
t = 1
t_0 = 0

prob = DiscreteProblem(birth_death, [N_0], (0.0,t), (K_b,K_d))
jump_prob = JumpProblem(birth_death, prob, Direct())
solutions = EvolutionaryModels.SamplesGillespie(jump_prob,SSAStepper(),100)
trajectories = EvolutionaryModels.solutionstoDF(solutions)
EvolutionaryModels.plotCRNGillespie(solutions,"Birth-Death")


# Let's compare
nmostres = 2000
FinalMostra = EvolutionaryModels.SampleFinal(nmostres,jump_prob,SSAStepper(),1)
Nmin = minimum(FinalMostra[1])
Nmax = maximum(FinalMostra[1])
orange = Nmin:1:Nmax
τ = (K_b-K_d)*(t-t_0)
σ = (K_b-K_d)/(K_b+K_d)
DistroUs = EvolutionaryModels.P_BD_n01.(σ,τ,orange)
histogram(FinalMostra[1],normalize = true,bins=100, color = "#BEDAA5")
#plot!(orange,Distro)
plot!(orange,DistroUs, color = "#1F654C",w=3)


##################### Birth-Death Process Starting from n_0 = 20
K_b = (5.5)
K_d = (4.5)
N_0 = 20
t = 1
t_0 = 0

prob = DiscreteProblem(birth_death, [N_0], (0.0,t), (K_b,K_d))
jump_prob = JumpProblem(birth_death, prob, Direct())
solutions = EvolutionaryModels.SamplesGillespie(jump_prob,SSAStepper(),100)
trajectories = EvolutionaryModels.solutionstoDF(solutions)
EvolutionaryModels.plotCRNGillespie(solutions)


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
solutions = EvolutionaryModels.SamplesGillespie(jump_prob,SSAStepper(),100)
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


```