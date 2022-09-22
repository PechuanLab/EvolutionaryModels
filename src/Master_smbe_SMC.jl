#=
Parameter scan for the barcode dynamics SMC inference
by Ximo Pechuan i Jorge 25/03/2021
=#

######### Load the model and accessory functions
#ARGoS = ["6077_PC1_R2__Top2_Population_Counts.csv","6077_PC1_R2__Top2_Complete_MasterSizes.csv"]
include("./Master_smbe_Types.jl") # Types
include("./Master_smbe_Functions.jl")

######### Prepare data to be fitted
data = CSV.read("../data/",DataFrame;drop=[1])
Tdata = TargetData(data)
# Comment for full data
#data=data[[x in ["T0","T1","T2","T3"] for x in data.TimePoint],:]
# Timepoints sampled
TimePoint = Tdata[1]
# Target data to fit
tdata = Tdata[2]
# Number of barcodes
barcodes = Tdata[3]
# Initial condition of each barcoded Lineage
n0 = Tdata[4]

######### TimeSeries Parameters
parameters = CSV.read(ARGS[2],DataFrame;drop=[1])
# Comment for full data
#parameters=parameters[[x in ["T0","T1","T2","T3"] for x in parameters.TimePoint],:]
Parametres = Params(parameters)
# Times between culutres
TimeCultures = Parametres[1]
# Cells transferred
Ntransferes = Parametres[2]
Passes = Parametres[3]
# Global Parameters
NPasses = Parametres[4]
#= 
Master Inference
=#
######### Multidimensional Uniform Prior on the model paramters
prior=Factored(Uniform(0,1), 
              Uniform(0,1),
              Uniform(0,1),
              Uniform(0,1),
              Uniform(0,1),
              Uniform(0,1),
              Uniform(0,1),
              Uniform(0,1),
              Uniform(0,1),
              Uniform(0,1),
              Uniform(0,1),
              Uniform(0,1))

prior2=Factored(
              Uniform(0,1),
              Uniform(0,1))

priorW=Factored(Uniform(0,1),# s_RG0
              Uniform(0,1),  #s_RG1
              Uniform(0,1),  #S_Remainder
              Uniform(0,0.5),   #r_0
              Uniform(0,0.5),   #r_1
              Uniform(0,0.5),   #r_2
              Uniform(1,2),   #K1
              Uniform(1,2),   #K2
              Uniform(1,2),   #K3
           )

priorM=Factored(Uniform(0,1),# s_RG0
              Uniform(0,1),  #s_RG1
              Uniform(0,1),  #S_Remainder
              Uniform(0,0.5),   #r_0
              Uniform(0,0.5),   #r_1
              Uniform(0,0.5),   #r_2
              Uniform(1,2),   #K1
              Uniform(1,2),   #K2
              Uniform(1,2),   #K3
              Uniform(0,1000), #M1
              Uniform(0,1000), #M2
              Uniform(0,1000) #M3
           )
 
#approx_density = ApproxKernelizedPosterior(prior,cost,750) 
# ABC parameter inference
# Sequential Monte Carlo
ressmc = smc(prior2, cost, nparticles=20, epstol=10,verbose=true,parallel=true)
#param_names = ["s_RG0", "s_RG1", "S_Remainder","r_0","r_1","r_2","K1","K2","K3"]
#param_names = ["s_RG0", "s_RG1", "S_Remainder","r_0","r_1","r_2","K1","K2","K3","M1","M2","M3"]
param_names = ["s_RG1", "S_Remainder"]
#Save the result
ABCResult = SaveParams(ressmc,param_names)
CSV.write(ARGS[3],ABCResult)
