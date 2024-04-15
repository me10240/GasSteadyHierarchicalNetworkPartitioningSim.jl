
function create_partition(filepath::AbstractString)::Dict{Any, Any}

    data = JSON.parsefile(filepath)
    partition = Dict{Any, Any}()
    num_partitions = data["num_partitions"]
    partition["num_partitions"] = num_partitions
    partition["num_interfaces"] = length(data["interface_nodes"])
    partition["interface_nodes"] = Vector{Int64}(data["interface_nodes"])
    partition["slack_network_ids"] = Vector{Int64}(data["slack_network_ids"])
    partition["slack_nodes"]  = Vector{Int64}(data["slack_nodes"])


    for i = 1: num_partitions

        partition[i] = Dict{String, Any}()
        partition[i]["node_list"] = Vector{Int64}(data[string(i)])
        partition[i]["interface"] = Vector{Int64}()
        partition[i]["bdry"] = Dict{Int, Vector{Int64}}()
        partition[i]["transfer"] = Dict{Int, Float64}()

        for j in data["interface_nodes"]
            if j in partition[i]["node_list"]
                push!(partition[i]["interface"], j)
                partition[i]["transfer"][j] = 0.0
            end
        end

    end
    
    order = Vector{Int64}()
    push!(order, 1)


    stage = 1
    partition["level"] = Dict{Int, Vector{Int64}}()
    partition["level"][stage] =  Vector{Int64}()
    push!(partition["level"][stage], 1) 

    while length(order) < num_partitions
        var = Vector{Int64}()
        for j in order
            for i in setdiff(collect(1:num_partitions), order)
                if length(intersect(partition[i]["node_list"], partition[j]["node_list"])) != 0
                    push!(var, i)
                    continue
                end
            end
        end
        stage = stage + 1
        partition["level"][stage] =  Vector{Int64}()
        for item in var
            push!(order, item)
            push!(partition["level"][stage], item)
        end
    end
    partition["num_level"] = stage

    for n = 1 : partition["num_level"]-1
        for sn_i in partition["level"][n]
            for sn_j in partition["level"][n+1]
                bdry = intersect(partition[sn_i]["node_list"], partition[sn_j]["node_list"])
                if length(bdry) == 1
                    partition[sn_i]["bdry"][sn_j] = Vector{Int64}(bdry)
                end
            end
        end
    end


    num_edges = partition["num_partitions"] + partition["num_interfaces"] - 1  # since tree

    global_to_local = Dict{Int64, Tuple{String, Int64}}()
    vertex_to_global_map = Dict{Any, Vector{Int64}}()

    for j in data["interface_nodes"]
        vertex_to_global_map[j] = Vector{Int64}()
    end


    edge_index = 1
    for i = 1: partition["num_partitions"] 
        vertex_to_global_map["N-$i"] = Vector{Int64}()

        interface_nodes = partition[i]["interface"]
        for j in interface_nodes
            global_to_local[edge_index] = ("N-$i", j)
            push!(vertex_to_global_map["N-$i"], edge_index)
            push!(vertex_to_global_map[j], edge_index)
            edge_index += 1
        end

    end
    partition["global_to_local"] = global_to_local
    partition["vertex_to_global_map"] = vertex_to_global_map
    
    return partition
end


function designate_interface_nodes_as_slack!(ssp_array::Vector{SteadySimulator}, partition::Dict{Any, Any})

    for level = 1 : partition["num_level"] - 1
        for sn_id_i in partition["level"][level]
            for sn_id_j in partition["level"][level+1]
                for node_id in get(partition[sn_id_i]["bdry"], sn_id_j, [])
                    ssp_array[sn_id_j].ref[:node][node_id]["is_slack"] = 1
                    ssp_array[sn_id_j].ref[:node][node_id]["withdrawal"] = NaN
                end
            end 
        end
    end

    return 

end

function update_interface_potentials_of_nbrs!(ssp_array::Vector{SteadySimulator}, partition::Dict{Any, Any}, level::Int64)

    if level == partition["num_level"]
        return
    end
    for sn_id_i in partition["level"][level]
        for sn_id_j in partition["level"][level + 1]
            for node_id in get(partition[sn_id_i]["bdry"], sn_id_j, [])
                ssp_array[sn_id_j].ref[:node][node_id]["potential"] =  ssp_array[sn_id_i].ref[:node][node_id]["potential"]
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
            else
                x_dofs[ss.ref[comp][id]["dof"]] = ssp_array[sn_id].ref[comp][id]["flow"]
            end
        end
    end

    return x_dofs
end
