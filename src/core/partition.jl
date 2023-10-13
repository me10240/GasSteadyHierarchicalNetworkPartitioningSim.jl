
function create_partition(filepath::AbstractString)::Dict{String, Any}


        data = JSON.parsefile(filepath)
        partition = Dict{String, Any}()
        partition["level"] = Dict{Int64, Any}()
        
        for i = 1:4
                partition["level"][i] =  Vector{Int64}()
        end
        push!(partition["level"][1], 1)



        num_partition = length(data)
        

        partition["subnetworks"] = Vector{Dict{String, Any}}()

        for i = 1: num_partition
                push!(partition["subnetworks"], Dict{String, Any}())
        end
        #
        for i = 1:num_partition
                @assert  length(data[string(i)]) > 1
                partition["subnetworks"][i]["node_list"] = data[string(i)]
                partition["subnetworks"][i]["interface"] = Dict{Int, Vector}()
        end
        #
        for sn = 1: num_partition
                for j in partition["level"][1]
                        if sn in  partition["level"][1]
                                continue
                        end
                        common_nodes = intersect(partition["subnetworks"][sn]["node_list"], partition["subnetworks"][j]["node_list"])
                
                        if length(common_nodes) != 0    
                                partition["subnetworks"][sn]["interface"][j] = common_nodes
                                push!(partition["level"][2], sn)  #these have interface with sn1
                        end
                end
        end
        #
        if length(union(partition["level"][1], partition["level"][2])) == num_partition
                println("All non-slack subnetworks are 2nd level...partitioning has radial layout")
                return partition
        end

        for sn = 1 : num_partition
                for j in partition["level"][2]
                        if sn in   union(partition["level"][1], partition["level"][2])
                                continue
                        end
                        common_nodes = intersect(partition["subnetworks"][sn]["node_list"], partition["subnetworks"][j]["node_list"])

                        if length(common_nodes) != 0
                                partition["subnetworks"][sn]["interface"][j] = common_nodes
                                push!(partition["level"][3], sn)  #these have interface with sn1
                        end                      
                end
        end 

        if length(union(partition["level"][1], partition["level"][2], partition["level"][3])) == num_partition
                println("At least one non-slack subnetwork is 3rd level...partitioning has chain layout")
                return partition
        end

        for sn = 1 : num_partition
                for j in partition["level"][3]
                        if sn in  union(partition["level"][1], partition["level"][2], partition["level"][3])
                                continue
                        end
                        common_nodes = intersect(partition["subnetworks"][sn]["node_list"], partition["subnetworks"][j]["node_list"])

                        if length(common_nodes) != 0
                                partition["subnetworks"][sn]["interface"][j] = common_nodes
                                push!(partition["level"][4], sn)  #these have interface with sn1
                        end                      
                end
        end

        if length(union(partition["level"][1], partition["level"][2], partition["level"][3], partition["level"][4])) == num_partition 
                println("At least one non-slack subnetwork is 4th level...partitioning has chain layout")
                return partition
        else

                error("This partition is not acceptable since at least one non-slack subnetwork is 5th level.")
        end 

        return partition
            

end



function create_interface_transfers!(ssp_array::Vector{SteadySimulator},  partition::Dict{String, Any}; val::Real=0.0)

    for level = 1:3
        for sn_id_i in partition["level"][level + 1]
            for sn_id_j in partition["level"][level]
                for node_id in get(partition["subnetworks"][sn_id_i]["interface"], sn_id_j, [])
                    if haskey(ssp_array[sn_id_j].ref[:node][node_id], "transfer") == false
                        ssp_array[sn_id_j].ref[:node][node_id]["transfer"] = Dict{Int, Float64}()
                    end
                    ssp_array[sn_id_j].ref[:node][node_id]["transfer"][sn_id_i] = val
                end
            end 
        end
    end 

    return

end

function interface_array(ssp_array::Vector{SteadySimulator}, partition::Dict{String,Any})::Vector{Float64}

    var = Vector{Float64}()
    
    for level = 1:3
        for sn_id_i in partition["level"][level + 1]
            for sn_id_j in partition["level"][level]
                for node_id in get(partition["subnetworks"][sn_id_i]["interface"], sn_id_j, [])
                    push!(var, ssp_array[sn_id_j].ref[:node][node_id]["transfer"][sn_id_i])
                end
            end
        end
    end


    return var
end

function designate_interface_nodes_as_slack!(ssp_array::Vector{SteadySimulator}, partition::Dict{String, Any})

    for level = 1:3
        for sn_id_i in partition["level"][level + 1]
            for sn_id_j in partition["level"][level]
                for node_id in get(partition["subnetworks"][sn_id_i]["interface"], sn_id_j, [])
                    ssp_array[sn_id_i].ref[:node][node_id]["is_slack"] = 1
                    ssp_array[sn_id_i].ref[:node][node_id]["withdrawal"] = NaN
                end
            end 
        end
    end

    return 

end


function update_interface_potentials_of_nbrs!(ssp_array::Vector{SteadySimulator}, partition::Dict{String, Any}, level::Int64)

    for sn_id_i in partition["level"][level + 1]
        for sn_id_j in partition["level"][level]
            for node_id in get(partition["subnetworks"][sn_id_i]["interface"], sn_id_j, [])
                ssp_array[sn_id_i].ref[:node][node_id]["potential"] =  ssp_array[sn_id_j].ref[:node][node_id]["potential"]
            end
        end
    end
    return
end

function update_interface_transfers!(ssp_array::Vector{SteadySimulator},  partition::Dict{String, Any})
   
    for level = 1:3
        for sn_id_i in partition["level"][level + 1]
            for sn_id_j in partition["level"][level]
                for node_id in get(partition["subnetworks"][sn_id_i]["interface"], sn_id_j, [])
                    ssp_array[sn_id_j].ref[:node][node_id]["transfer"][sn_id_i] = 0 * ssp_array[sn_id_j].ref[:node][node_id]["transfer"][sn_id_i] - 1.0 * ssp_array[sn_id_i].ref[:node][node_id]["withdrawal"]
                end
            end 
        end
    end 
    return
end

function combine_subnetwork_solutions(ss::SteadySimulator, ssp_array::Vector{SteadySimulator})::Vector{Float64}
    ndofs = length(ref(ss, :dof))
    x_dofs = zeros(Float64, ndofs) 
    dofs_updated = 0

    for sn_id = 1 : length(ssp_array)
        for i = 1:length(ref(ssp_array[sn_id], :dof))
            comp, id = ssp_array[sn_id].ref[:dof][i]
            if comp == :node
                x_dofs[ss.ref[comp][id]["dof"]] = ssp_array[sn_id].ref[comp][id]["potential"]
            elseif comp == :pipe
                x_dofs[ss.ref[comp][id]["dof"]] = ssp_array[sn_id].ref[comp][id]["flow"]
            elseif comp == :compressor
                x_dofs[ss.ref[comp][id]["dof"]] = ssp_array[sn_id].ref[comp][id]["flow"]
            end
        end
    end

    return x_dofs
end
