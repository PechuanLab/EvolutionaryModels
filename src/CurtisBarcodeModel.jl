######################### Required libraries
using KissABC
using Distributions
using ClusterManagers
using DataFrames
using Setfield
using StatsBase
using CSV
using Distances
using Plots

#########################
#= 
Ancillary Functions  s
=#

function SaveParams(ressmc,param_names)
	#Save result
	ABCResult = DataFrame()
	for i = 1:length(param_names)
		colname = param_names[i]
		ABCResult[!,colname]= ressmc[1][i].particles  
	end
	return ABCResult
end

function Params(parameters::DataFrame)
	# Prepares Paramters 
	# Times between culutres
	TimeCultures = parameters[:,:CultureTime]
	# Cells transferred
	Ntransferes = parameters[:,:CellsTransfer]
	Passes = parameters[:,:TimePoint]
	# Global Parameters
	NPasses = size(Passes,1) 
	return TimeCultures,Ntransferes,Passes,NPasses
end

function TargetData(data::DataFrame)
	# Prepares target Data 
	TimePoint = data[:,"TimePoint"]  
	# Target data to fit
	tdata = data[:,1:size(data,2)-1]
	tdata = Matrix(tdata)
	barcodes = size(data,2)-1
	# Initial condition of each barcoded Lineage
	n0 = [tdata[1,j] for j=1:barcodes]
	return TimePoint,tdata,barcodes,n0
end

function unwrapper(poblacio::Population, fieldname)
    x = fill(0.0,length(poblacio.clones))
    for i=1:length(poblacio.clones)
        x[i]=getfield(poblacio.clones[i],Symbol(fieldname)) 
    end
    return x
end

function measurement_Sample(BarCode_mat,TimePoint,Passes,data)
	# Sample according to the experimental obsvervations
	df_BarcodeMat = DataFrame(BarCode_mat,:auto)
    df_BarcodeMat.TimePoint = Passes
	rename!(df_BarcodeMat,names(data)) 
    subset_BarcodeMat = df_BarcodeMat[[x in TimePoint for x in df_BarcodeMat.TimePoint],:]
	return subset_BarcodeMat
end 

function AverageFitSIG(t,w_0,K,r)
	# Sigmoid average fitness
	# k = K*w_0
	w = K*w_0*w_0*exp(r*t)/(K*w_0+w_0*(exp(r*t)-1))
	return w
end

function AverageFitSIG(t,w_0,K,r,M)
	# Sigmoid average fitness
	# k = K*w_0
	w = w_0+(K*w_0-w_0)/(K*w_0+(exp(-r*(t-M))))
	return w
end


######################################################### Initialize Population

function InitPop(barcodes::Int,s_coef::Array{Float64,1},n0::Array{Float64,1},passmut,fitnessMut,nmut)
	#initializes population
    lin1 = LineageMutPass(1,s_coef[1],n0[1],passmut[1],fitnessMut[1],nmut[1])
    population = fill(lin1,1)
	for i in 2:barcodes
		 lin2 = LineageMutPass(i,s_coef[i],n0[i],passmut[i],fitnessMut[i],nmut[i])
		 push!(population,lin2) 
	end
	population = Population(population)
return population
end

function InitPop(barcodes::Int,s_coef::Array{Float64,1},n0::Array{Float64,1},passmut,fitnessMut,nmut,mutates)
	#initializes population
    lin1 = LineageMutTime(1,s_coef[1],n0[1],passmut[1],fitnessMut[1],nmut[1],mutates[1])
    population = fill(lin1,1)
	for i in 2:barcodes
		 lin2 = LineageMutTime(i,s_coef[i],n0[i],passmut[i],fitnessMut[i],nmut[i],mutates[i])
		 push!(population,lin2) 
	end
	population = Population(population)
return population
end

function InitPop(barcodes::Int,s_coef::Array{Float64,1},n0::Array{Float64,1})
	#initializes population
    lin1 = Lineage(1,s_coef[1],n0[1])
    population = fill(lin1,1)
	for i in 2:barcodes
		 lin2 = Lineage(i,s_coef[i],n0[i])
		 push!(population,lin2) 
	end
	population = Population(population)
return population
end

function InitPop(barcodes::Int,s_coef::Array{Float64,1},n0::Array{Float64,1},K::Array{Float64,1},r::Array{Float64,1})
	#initializes population
    lin1 = LineageW(1,s_coef[1],n0[1],K[1],r[1])
    population = fill(lin1,1)
	for i in 2:barcodes
		 lin2 = LineageW(i,s_coef[i],n0[i],K[i],r[i])
		 push!(population,lin2) 
	end
	population = Population(population)
return population
end

function InitPop(barcodes::Int,s_coef::Array{Float64,1},n0::Array{Float64,1},K::Array{Float64,1},r::Array{Float64,1},M::Array{Float64,1})
	#initializes population
    lin1 = LineageW2(1,s_coef[1],n0[1],K[1],r[1],M[1])
    population = fill(lin1,1)
	for i in 2:barcodes
		 lin2 = LineageW2(i,s_coef[i],n0[i],K[i],r[i],M[i])
		 push!(population,lin2) 
	end
	population = Population(population)
return population
end

########################################################## Cell Culture Growth

function CulturePasss(poblacio::Population,TimeCulture)
	if TimeCulture == 0
	    poblacio = poblacio
	else
		# Exponential assumption
		fitness_vec = unwrapper(poblacio,"fitness")
		N_vec = unwrapper(poblacio,"N")
		mass_action = N_vec.*exp.(TimeCulture*fitness_vec)
		sdev_action = ((N_vec.*exp.(2*TimeCulture*fitness_vec)).*(ones(length(fitness_vec))-exp.(-TimeCulture*fitness_vec))).^(1/2)
		for i in 1:length(poblacio.clones)
			mu = mass_action[i]
			sigma = sdev_action[i]
			d = Normal(mu,sigma)
			x = rand(d, 1)
			if x[1] < 0
				x[1] = 0
			end
	    	poblacio.clones[i].N = x[1]
		end
	end
	return poblacio
end

# Passage mut
function CulturePasss(poblacio::Population,TimeCulture,NPass)
	# Normal approximation Yule
	if TimeCulture == 0
	    poblacio = poblacio
	else
		# Obtain parameters
		fitness_vec = unwrapper(poblacio,"fitness")
		N_vec = unwrapper(poblacio,"N")
		PassMut_vec = unwrapper(poblacio,"PassMut")
		# Mutate
		for i in 1:length(poblacio.clones)
			# Initiate the mutation 
			if PassMut_vec[i] == NPass
			   # For the Yule process we are fine
			   poblacio.clones[i].NMut = 1
			end
		end
		# For clonal
		mass_action = N_vec.*exp.(TimeCulture*fitness_vec)
		sdev_action = ((N_vec.*exp.(2*TimeCulture*fitness_vec)).*(ones(length(fitness_vec))-exp.(-TimeCulture*fitness_vec))).^(1/2)
		for i in 1:length(poblacio.clones)
			mu = mass_action[i]
			sigma = sdev_action[i]
			d = Normal(mu,sigma)
			x = rand(d, 1)
			if x[1] < 0
				x[1] = 0
			end
	    	poblacio.clones[i].N = x[1]
		end
		# For sub-clonal
		# Get mutants
		NMut_vec = unwrapper(poblacio,"NMut")
        fitnessMut_vec = unwrapper(poblacio,"fitnessMut")
        # Yule
		mass_action_mut = NMut_vec.*exp.(TimeCulture*fitnessMut_vec)
		sdev_action_mut = ((NMut_vec.*exp.(2*TimeCulture*fitnessMut_vec)).*(ones(length(fitnessMut_vec))-exp.(-TimeCulture*fitnessMut_vec))).^(1/2)
		for i in 1:length(poblacio.clones)
			mu = mass_action_mut[i]
			sigma = sdev_action_mut[i]
			d = Normal(mu,sigma)
			x = rand(d, 1)
			if x[1] < 0
				x[1] = 0
			end
	    	poblacio.clones[i].NMut = x[1]
		end

	end
	return poblacio
end
# time mut 

function CulturePasss(poblacio::Population,TimeCulture,T1,T2,mutated)
	# Normal approximation Yule
	if TimeCulture == 0
	    poblacio = poblacio
	else
		# Obtain parameters
		fitness_vec = unwrapper(poblacio,"fitness")
		N_vec = unwrapper(poblacio,"N")
		# Get mutants
		NMut_vec = unwrapper(poblacio,"NMut")
        fitnessMut_vec = unwrapper(poblacio,"fitnessMut")
		PassMut_vec = unwrapper(poblacio,"τ")
		Mutates_vec = unwrapper(poblacio,"mutates")
		# For clonal
		mass_action = N_vec.*exp.(TimeCulture*fitness_vec)
		sdev_action = ((N_vec.*exp.(2*TimeCulture*fitness_vec)).*(ones(length(fitness_vec))-exp.(-TimeCulture*fitness_vec))).^(1/2)
		for i in 1:length(poblacio.clones)
			mu = mass_action[i]
			sigma = sdev_action[i]
			d = Normal(mu,sigma)
			x = rand(d, 1)
			if x[1] < 0
				x[1] = 0
			end
	    	poblacio.clones[i].N = x[1]
		end
		# For sub-clonal
        # Yule
		mass_action_mut = NMut_vec.*exp.(TimeCulture*fitnessMut_vec)
		sdev_action_mut = ((NMut_vec.*exp.(2*TimeCulture*fitnessMut_vec)).*(ones(length(fitnessMut_vec))-exp.(-TimeCulture*fitnessMut_vec))).^(1/2)
		for i in 1:length(poblacio.clones)
			mu = mass_action_mut[i]
			sigma = sdev_action_mut[i]
			d = Normal(mu,sigma)
			x = rand(d, 1)
			if x[1] < 0
				x[1] = 0
			end
	    	poblacio.clones[i].NMut = x[1]
		end
		# Mutate
		for i in 1:length(poblacio.clones)
			# Initiate the mutation 
			if (PassMut_vec[i] >= T1 && mutated[i] == false && Mutates_vec[i] == true && PassMut_vec[i] <= T2)
			   # For the Yule process we are fine
			   ΔT = T2 - PassMut_vec[i]
			   mass_action = exp(ΔT*fitness_vec[i])
		       sdev_action = ((exp(2*ΔT*fitnessMut_vec[i]))*(1-exp(-ΔT*fitnessMut_vec[i])))^(1/2)
		       d = Normal(mass_action,sdev_action)
			   x = rand(d, 1)
					if x[1] < 0
						x[1] = 0
					end
		       poblacio.clones[i].NMut = x[1]+1
		       mutated[i] = true
			end
		end
	end
	return poblacio
end

# Mut fitness W
function CulturePasss(poblacio::Population,TimeCulture,T1)
	if TimeCulture == 0
	    poblacio = poblacio
	else
		# Exponential assumption
		w_vec = unwrapper(poblacio,"w_0")
		k_vec = unwrapper(poblacio,"K")
		r_vec = unwrapper(poblacio,"r")
		N_vec = unwrapper(poblacio,"N")
		fitness_vec = AverageFitSIG.(T1,w_vec,k_vec,r_vec) 
		mass_action = N_vec.*exp.(TimeCulture*fitness_vec)
		sdev_action = ((N_vec.*exp.(2*TimeCulture*fitness_vec)).*(ones(length(fitness_vec))-exp.(-TimeCulture*fitness_vec))).^(1/2)
		for i in 1:length(poblacio.clones)
			mu = mass_action[i]
			sigma = sdev_action[i]
			d = Normal(mu,sigma)
			x = rand(d, 1)
			if x[1] < 0
				x[1] = 0
			end
	    	poblacio.clones[i].N = x[1]
		end
	end
	return poblacio
end

# Mut fitness W2
function CulturePasss(poblacio::Population,TimeCulture,T1)
	if TimeCulture == 0
	    poblacio = poblacio
	else
		# Exponential assumption
		w_vec = unwrapper(poblacio,"w_0")
		k_vec = unwrapper(poblacio,"K")
		r_vec = unwrapper(poblacio,"r")
		N_vec = unwrapper(poblacio,"N")
		M_vec = unwrapper(poblacio,"M")
		fitness_vec = AverageFitSIG.(T1,w_vec,k_vec,r_vec,M_vec) 
		mass_action = N_vec.*exp.(TimeCulture*fitness_vec)
		sdev_action = ((N_vec.*exp.(2*TimeCulture*fitness_vec)).*(ones(length(fitness_vec))-exp.(-TimeCulture*fitness_vec))).^(1/2)
		for i in 1:length(poblacio.clones)
			mu = mass_action[i]
			sigma = sdev_action[i]
			d = Normal(mu,sigma)
			x = rand(d, 1)
			if x[1] < 0
				x[1] = 0
			end
	    	poblacio.clones[i].N = x[1]
		end
	end
	return poblacio
end



####################################################################### Subculturing Function
function TransferMut(poblacio::Population,Ntransfer::Int)
	# Multinomial sampling
	N_vec = unwrapper(poblacio,"N")
	NMut_vec = unwrapper(poblacio,"NMut")
	NTotal_vec= vcat(N_vec,NMut_vec)
	total=sum(NTotal_vec)
	Linage_freq=NTotal_vec/total
	d=Categorical(Linage_freq)
	R = rand(d,Ntransfer)
	NewPopN = counts(R,2*barcodes)
	# Exact
	for i in 1:length(poblacio.clones)
	    poblacio.clones[i].N  = NewPopN[i]
	    poblacio.clones[i].NMut  = NewPopN[i+barcodes] 
	end
	return poblacio
end

function Transfer(poblacio::Population,Ntransfer::Int,barcodes::Int)
	#function body
	N_vec = unwrapper(poblacio,"N")
	total=sum(N_vec)
	Linage_freq=N_vec/total
	d=Categorical(Linage_freq)
	R = rand(d,Ntransfer)
	NewPopN = counts(R,barcodes)
	for i in 1:length(poblacio.clones)
	    poblacio.clones[i].N  = NewPopN[i] 
	end
	return poblacio
end



############################## Main Function ##################

#=
barcodes = 3
s_coef = [0.3,0.25,0.2]
n0 = [10000.0,10000.0,10000.0]
TimeCulture = 16
NPasses = 10
passmut = [0,4,0]
fitnessMut = [0,0.5,0]
nmut = [0,0,0]
Ntransferes = Int.(100000*ones(10))
K = [0.2,0.4,0.1]
r = [0.04,0.05,0.06]
=#

## Basic function
function SimPop(barcodes::Int,s_coef::Array{Float64,1},
	n0::Array{Int64,1},TimeCulture,NPasses::Int)
	# Things we will record
	BarCode_mat = zeros(Float64,NPasses,barcodes)
	# Initial Population
	poblacio = InitPop(barcodes,s_coef,n0)
	# Culture cycles
	for i in 1:NPasses
	    poblacio = CulturePasss(poblacio::Population,TimeCulture)
	    N_vec = unwrapper(poblacio,"N")
	    BarCode_mat[i,:] = N_vec
	    poblacio = Transfer(poblacio::Population,Ntransfer::Int,barcodes::Int)
	end
	return BarCode_mat
end

function SimPop(barcodes,s_coef,n0,TimeCultures,NPasses,passmut,fitnessMut,nmut,Ntransferes)
	# Pass
	NPass = 1
	# Things we will record
	BarCode_mat = zeros(Float64,NPasses,barcodes)
	BarCode_mat_mut = zeros(Float64,NPasses,barcodes)
	# Initial Population
	poblacio = InitPop(barcodes,s_coef,n0,passmut,fitnessMut,nmut)
	# Culture cycles
	for i in 1:NPasses
		TimeCulture = TimeCultures[i]
		Ntransfer = Ntransferes[i]
	    poblacio = CulturePasss(poblacio::Population,TimeCulture,NPass)
	    N_vec = unwrapper(poblacio,"N")
	    N_mut = unwrapper(poblacio,"NMut")
	    BarCode_mat[i,:] = N_vec
	    BarCode_mat_mut[i,:] = N_mut
	    poblacio = Transfer(poblacio::Population,Ntransfer)
	    NPass = NPass+1
	end
	return BarCode_mat,BarCode_mat_mut
end


function SimPop(barcodes,s_coef,n0,TimeCultures,NPasses,passmut,fitnessMut,nmut,Ntransferes,mutates,mutated)
	# Pass
	T1 = 0.0
	# Things we will record
	BarCode_mat = zeros(Float64,NPasses,barcodes)
	BarCode_mat_mut = zeros(Float64,NPasses,barcodes)
	# Initial Population
	poblacio = InitPop(barcodes,s_coef,n0,passmut,fitnessMut,nmut,mutates)
	# Culture cycles
	for i in 1:NPasses
		TimeCulture = TimeCultures[i]
		Ntransfer = Ntransferes[i]
		T2 = T1+TimeCulture
	    poblacio = CulturePasss(poblacio::Population,TimeCulture,T1,T2,mutated)
	    N_vec = unwrapper(poblacio,"N")
	    N_mut = unwrapper(poblacio,"NMut")
	    BarCode_mat[i,:] = N_vec
	    BarCode_mat_mut[i,:] = N_mut
	    poblacio = Transfer(poblacio::Population,Ntransfer)
	    T1 = T1+TimeCulture
	end
	return BarCode_mat,BarCode_mat_mut
end


function SimPop(barcodes,s_coef,n0,K,r,TimeCultures,NPasses,Ntransferes)
	# Pass
	T1 = 0.0
	# Things we will record
	BarCode_mat = zeros(Float64,NPasses,barcodes)
	BarCode_mat_mut = zeros(Float64,NPasses,barcodes)
	# Initial Population
	poblacio = InitPop(barcodes,s_coef,n0,K,r)
	# Culture cycles
	for i in 1:NPasses
		TimeCulture = TimeCultures[i]
		Ntransfer = Ntransferes[i]
	    poblacio = CulturePasss(poblacio::Population,TimeCulture,T1)
	    N_vec = unwrapper(poblacio,"N")
	    BarCode_mat[i,:] = N_vec
	    poblacio = Transfer(poblacio::Population,Ntransfer)
	    T1 = T1+TimeCulture
	end
	return BarCode_mat
end

function SimPop(barcodes,s_coef,n0,K,r,M,TimeCultures,NPasses,Ntransferes)
	# Pass
	T1 = 0.0
	# Things we will record
	BarCode_mat = zeros(Float64,NPasses,barcodes)
	BarCode_mat_mut = zeros(Float64,NPasses,barcodes)
	# Initial Population
	poblacio = InitPop(barcodes,s_coef,n0,K,r,M)
	# Culture cycles
	for i in 1:NPasses
		TimeCulture = TimeCultures[i]
		Ntransfer = Ntransferes[i]
	    poblacio = CulturePasss(poblacio::Population,TimeCulture,T1)
	    N_vec = unwrapper(poblacio,"N")
	    BarCode_mat[i,:] = N_vec
	    poblacio = Transfer(poblacio::Population,Ntransfer)
	    T1 = T1+TimeCulture
	end
	return BarCode_mat
end

################################# Simulacion
function sim((s_RG0,s_RG1,s_RG2,s_RG3,s_RG4,s_RG5,s_RG6,s_RG7,s_RG8,s_RG9,s_RG10,S_Remainder)) 
  # Selective coefficients     
  s_coef = [s_RG0,s_RG1,s_RG2,s_RG3,s_RG4,s_RG5,s_RG6,s_RG7,s_RG8,s_RG9,s_RG10,S_Remainder]
  # Things we will record
  BarCode_mat = zeros(Float64,NPasses,barcodes)
  # Initial Population
  poblacio = InitPop(barcodes,s_coef,n0)
  # Culture cycles
  for i in 1:NPasses
      TimeCulture = TimeCultures[i]
      Ntransfer = Ntransferes[i]
      poblacio = CulturePasss(poblacio::Population,TimeCulture)
      N_vec = unwrapper(poblacio,"N")
      BarCode_mat[i,:] = N_vec
      poblacio = Transfer(poblacio::Population,Ntransfer::Int,barcodes::Int)
  end
  subset_BarcodeMat = measurement_Sample(BarCode_mat,TimePoint,Passes,data)
  return subset_BarcodeMat
end

function simtest((s_RG0,s_RG1,S_Remainder),NPasses,barcodes,TimeCultures,Ntransferes,n0,TimePoint) 
  # Fix
  # Selective coefficients     
  s_coef = [s_RG0,s_RG1,S_Remainder]
  # Things we will record
  BarCode_mat = zeros(Float64,NPasses,barcodes)
  # Initial Population
  poblacio = InitPop(barcodes,s_coef,n0)
  # Culture cycles
  for i in 1:NPasses
      TimeCulture = TimeCultures[i]
      Ntransfer = Ntransferes[i]
      poblacio = CulturePasss(poblacio::Population,TimeCulture)
      N_vec = unwrapper(poblacio,"N")
      BarCode_mat[i,:] = N_vec
      poblacio = Transfer(poblacio::Population,Ntransfer::Int,barcodes::Int)
  end
  subset_BarcodeMat = measurement_Sample(BarCode_mat,TimePoint,Passes,data)
  return subset_BarcodeMat
end

function sim((s_RG0,s_RG1,S_Remainder,passmut_1,s_mut)) 
  # Selective coefficients     
  s_coef = [s_RG0,s_RG1,S_Remainder]
  # Passage for mutation to appear
  passmut = [0,passmut_1,0]
  # Fitness mutation
  fitnessMut = [0,s_mut,0]
  # Initial mutants 
  nmut = [0,0,0]
  # Simulate
  simulacio = SimPop(barcodes,s_coef,n0,TimeCultures,NPasses,passmut,fitnessMut,nmut,Ntransferes)
  # Look at the added matrix
  BarCode_mat = simulacio[1]+simulacio[2]
  # Subset to compare with data
  subset_BarcodeMat = measurement_Sample(BarCode_mat,TimePoint,Passes,data)
  return subset_BarcodeMat
end
#= 
s_RG0 = 0.4
s_RG1 = 0.5
S_Remainder =0.1
r_0 = 0.2
r_1 = 0.4
r_2 = 0.5
K_1 = 0.5
K_2 = 0.5
K_3 = 0.5
=#
function sim((s_RG0,s_RG1,S_Remainder,r_0,r_1,r_2,K_1,K_2,K_3)) 
  # Selective coefficients     
  s_coef = [s_RG0,s_RG1,S_Remainder]
  # fitness change rate
  r = [r_0,r_1,r_2]
  # Maximum fitness
  K = [K_1,K_2,K_3]
  # Simulate
  simulacio = SimPop(barcodes,s_coef,n0,K,r,TimeCultures,NPasses,Ntransferes)
  # Look at the added matrix
  BarCode_mat = simulacio
  # Subset to compare with data
  subset_BarcodeMat = measurement_Sample(BarCode_mat,TimePoint,Passes,data)
  return subset_BarcodeMat
end

function sim((s_RG0,s_RG1,S_Remainder,r_0,r_1,r_2,K_1,K_2,K_3,M1,M2,M3)) 
  # Selective coefficients     
  s_coef = [s_RG0,s_RG1,S_Remainder]
  # fitness change rate
  r = [r_0,r_1,r_2]
  # Maximum fitness
  K = [K_1,K_2,K_3]
  M = [M1,M2,M3]
  # Simulate
  simulacio = SimPop(barcodes,s_coef,n0,K,r,M,TimeCultures,NPasses,Ntransferes)
  # Look at the added matrix
  BarCode_mat = simulacio
  # Subset to compare with data
  subset_BarcodeMat = measurement_Sample(BarCode_mat,TimePoint,Passes,data)
  return subset_BarcodeMat
end




