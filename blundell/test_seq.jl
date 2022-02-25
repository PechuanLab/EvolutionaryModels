
#=
sequence function 
=#

"""

...
# Arguments
- `pop::Population: plug in Population 
...
"""

function sequence(pop::Population,mean_depth::Float64)

	N=pop_size(pop)

	variant_dict=Dict()
	lov_list=[]
    lov_dic=Dict()
    

for set in pop.clones

    #println(set)
    #println(set.Mutations)
    #println(set.N)

    clone1=deepcopy(set.Mutations)
    size=set.N

    if length(clone1)>1 

        deleteat!(clone1, 1)

        for  mutation in clone1
            mutat=deepcopy(mutation)
            mID=mutat.ID

            if mID in keys(variant_dict)

                true_freq=variant_dict[mID][2]+0.5*(size/N)
                variant_dict[mID][2]=true_freq

            else 
                variant_dict[mID]=["true_freq",0.5*(size/N)]

            end
        end

    end

end 

for set in variant_dict
    m=set[1]
    true_freq=set[2][2]

    x = rand(Normal(0,0.3), 1)
    x2=exp(x[1])
    x3=mean_depth*x2
    depth1=round(Int64,x3) 

    expected_reads1=true_freq*depth1
    reads1=rand(Poisson(expected_reads1),1)[1]

    y = rand(Normal(0,0.3), 1)
    y2=exp(y[1])
    y3=mean_depth*y2
    depth2=round(Int64,y3) 

    expected_reads2=true_freq*depth2
    reads2=rand(Poisson(expected_reads2),1)[1]


    if reads1<4
        reads1=0
    end 
    if reads2<4
        reads2=0
    end 
    freq1=reads1/depth1
    freq2=reads2/depth2

    lov_dic["mutation_ID"]=m
    lov_dic["true_freq"]=true_freq
    lov_dic["freq1"]=freq1
    lov_dic["reads1"]=reads1
    lov_dic["depth1"]=depth1 
    lov_dic["freq2"]=freq2
    lov_dic["reads2"]=reads2
    lov_dic["depth2"]=depth2 

    push!(lov_list,lov_dic)

end


return lov_list

end  