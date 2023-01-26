#=
How to perform smc on simulated dataset.
=#
reference = TypeParamsVector

# Define the fixed string
fixed_string = "P"
# Define the vector of numbers
numbers = 1:length(reference)
# Use the `map` function to concatenate the fixed string to each number
result = map(x -> fixed_string * string(x), numbers)
# Cost function
InferenceParameterVector = Tuple(reference)

function SimulationCostFunc(InferenceParameterVector)
	params = [i for i in InferenceParameterVector]
    x = BatchCulture(LinType,Barcodes,n0,params,NPasses,CulturePass,TimeCulture,TransferFunction,Ntransfer)
    y = tdata
    d = (sum(colwise(Euclidean(), x, y)))^(1/2)
    d
end

prior=Factored(
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
ressmc = smc(prior, SimulationCostFunc, nparticles=100, epstol=0.5,verbose=true,parallel=true)
# Generate the simulated data
tdata = BatchCulture(LinType,Barcodes,n0,TypeParamsVector,NPasses,CulturePass,TimeCulture,TransferFunction,Ntransfer)

