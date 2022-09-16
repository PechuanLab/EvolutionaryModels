#=
Model Plots scan for the barcode dynamics SMC inference
by Ximo Pechuan i Jorge 25/03/2021
=#

######### Load the model and accessory functions
#ARGoS = ["6077_PC1_R2__Top2_Population_Counts.csv","6077_PC1_R2__Top2_Complete_MasterSizes.csv"]
include("./Master_smbe_Types.jl") # Types
include("./Master_smbe_Functions.jl")

######### Prepare data to be fitted
data = CSV.read(ARGS[1],DataFrame;drop=[1])
Tdata = TargetData(data)
# Comment for full data
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
Parametres = Params(parameters)
# Times between culutres
TimeCultures = Parametres[1]
# Cells transferred
Ntransferes = Parametres[2]
Passes = Parametres[3]
# Global Parameters
NPasses = Parametres[4]

#Read the parameter file
Nreplicates = 10
s = CSV.read(ARGS[3],DataFrame;)
mean_s = mean.(eachcol(s))
# sampled

for i = 1:Nreplicates
  simulated  = sim(mean_s) 
  name = "R$i"*ARGS[4]
  CSV.write(name,simulated)
end
