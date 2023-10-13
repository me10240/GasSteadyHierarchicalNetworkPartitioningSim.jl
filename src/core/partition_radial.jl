
function create_partition(filepath::AbstractString)::Dict{String, Any}


        data = JSON.parsefile(filepath)
        partition = Dict{String, Any}()
        # partition["zeroth_level_subnetworks"] =  Vector{Int64}()
        # push!(partition["zeroth_level_subnetworks"], 1)
        # partition["first_level_subnetworks"] =  Vector{Int64}()
        # partition["second_level_subnetworks"] =  Vector{Int64}()

        num_partition = length(data)
        

        partition["subnetworks"] = Vector{Dict{String, Any}}()

        for i = 1: num_partition
                push!(partition["subnetworks"], Dict{String, Any}())
        end
        #
        for i = 1:num_partition
                partition["subnetworks"][i]["node_list"] = data[string(i)]
                partition["subnetworks"][i]["interface"] = Dict{Int, Vector}()
        end
        #
        for j = 2: num_partition
                common_nodes = intersect(partition["subnetworks"][1]["node_list"], partition["subnetworks"][j]["node_list"])
                
                if length(common_nodes) != 0    
                        partition["subnetworks"][1]["interface"][j] = common_nodes
                        # push!(partition["first_level_subnetworks"], j)  #these have interface with sn1
                end
        end
        #
        # if length(partition["first_level_subnetworks"]) == num_partition - 1
        #         println("All non-slack subnetworks are first level...partitioning has radial layout")
        #         return partition
        # end

        # for sn = 2 : num_partition
        #         for j in partition["first_level_subnetworks"]
        #                 if sn in   partition["first_level_subnetworks"]
        #                         continue
        #                 end
        #                 common_nodes = intersect(partition["subnetworks"][sn]["node_list"], partition["subnetworks"][j]["node_list"])

        #                 if length(common_nodes) != 0
        #                         partition["subnetworks"][sn]["interface"][j] = common_nodes
        #                         push!(partition["second_level_subnetworks"], sn)  #these have interface with sn1
        #                 end                      
        #         end
        # end 

        # if length(partition["first_level_subnetworks"]) + length(partition["second_level_subnetworks"])== num_partition - 1
        #         println("At least one non-slack subnetwork is second level...partitioning has chain layout")
        #         return partition
        # else

        #         error("This partition is not acceptable since at least one non-slack subnetwork is third level.")
        # end

        return partition
            

end


function create_interface_transfers!(ssp::SteadySimulator,  partition::Vector{Dict{String, Any}})

    for sn_id = 2 : length(partition) 
        for node_id in partition[1]["interface"][sn_id]
            if haskey(ssp.ref[:node][node_id], "transfer") == false
                ssp.ref[:node][node_id]["transfer"] = Dict{Int, Float64}()
            end
            ssp.ref[:node][node_id]["transfer"][sn_id] = 0.0 

        end
    end
    return

end

function interface_array(ssp::SteadySimulator, partition::Vector{Dict{String,Any}})::Vector{Float64}

    var = Vector{Float64}()
    for sn_id = 2 : length(partition) 
        for node_id in partition[1]["interface"][sn_id]
            push!(var, ssp.ref[:node][node_id]["transfer"][sn_id])
        end
    end

return var
end

function designate_interface_nodes_as_slack!(ssp_array::Vector{SteadySimulator}, partition::Vector{Dict{String, Any}})

    for sn_id = 2 : length(partition)
        for node_id in partition[1]["interface"][sn_id]
            ssp_array[sn_id].ref[:node][node_id]["is_slack"] =  1
            ssp_array[sn_id].ref[:node][node_id]["withdrawal"] =  NaN 
        end
    end

end

function update_interface_potentials_of_first_level_neighbours!(ssp_array::Vector{SteadySimulator}, partition::Vector{Dict{String, Any}})

    for sn_id = 2 : length(partition)
        for node_id in partition[1]["interface"][sn_id]
            ssp_array[sn_id].ref[:node][node_id]["potential"] =  ssp_array[1].ref[:node][node_id]["potential"]
        end
    end

end

function update_interface_transfers!(ssp_array::Vector{SteadySimulator},  partition::Vector{Dict{String, Any}})
    for sn_id = 2 : length(partition) 
        for node_id in partition[1]["interface"][sn_id]
            ssp_array[1].ref[:node][node_id]["transfer"][sn_id] = 0.5 * ssp_array[1].ref[:node][node_id]["transfer"][sn_id]   - 0.5 * ssp_array[sn_id].ref[:node][node_id]["withdrawal"]
        end
    end
    return
end

function combine_subnetwork_solutions!(ss::SteadySimulator, ssp_array::Vector{SteadySimulator})
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
