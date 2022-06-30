using Plots

#=
 Sampling and plotting accessory functions
=#

function plotCRNGillespie(sol_array, titleEq)
    trajectories = plot(
        sol_array[1],
        xlabel = "Time",
        ylabel = "Number of Cells",
        thickness_scaling = 1.5,
        legend = false,
        title = titleEq,
    )
    # Plots
    for i = 2:length(sol_array)
        plot!(
            trajectories,
            sol_array[i],
            xlabel = "Time",
            thickness_scaling = 1.5,
            legend = false,
        )
    end

    return trajectories
end


function SamplesGillespie(jump_prob, solver, nsamples)
    solutions = []
    # get samples from Gillespie Solutions
    for i = 1:nsamples
        sol = solve(jump_prob, solver)
        push!(solutions, sol)
    end
    return solutions
end

"""
Prepares a data frame containing the trajectories for plotting
...
# Arguments
- `sol_array::Array{Any}`: Array of solutions comimng from SamplesGillespie()
...
"""
function solutionstoDF(sol_array::Array{Any})
    df = DataFrame()
    df[!, :Time] = sol_array[1].t
    df[!, :Time] = sol_array[1].u
    df[!, :Replicate] = fill("R1", size(df, 1))
    # Plots
    for i = 2:length(sol_array)
        df1 = DataFrame()
        df1[!, :Time] = sol_array[1].t
        df1[!, :Time] = sol_array[1].u
        df1[!, :Replicate] = fill("R$i", size(df1, 1))
        df = [df; df1]
    end
    return df
end

function SampleFinal(nmostres::Int, problemata::JumpProblem, solver, species::Int)
    Finals = []
    for i = 1:nmostres
        #for body
        sol = solve(problemata, solver)
        append!(Finals, sol[end][species])
    end
    mitja = mean(Finals)
    desvest = (std(Finals))^2
    return Finals, mitja, desvest
end
