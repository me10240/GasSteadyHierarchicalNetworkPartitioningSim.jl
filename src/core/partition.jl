function read_partition_file(filepath::AbstractString)::Dict{Any, Any}
    data = JSON.parsefile(filepath)
    return load_partition_data(data)
end

function create_graph(ss::SteadySimulator)::Tuple{SimpleGraph, Dict{Any, Int64}, Vector{Any}}
    V = ref(ss, :node) |> collect 
    new_node_from_old = Dict(V[i][1] => i for i in range(1, length(V)))
    old_node_from_new =[V[i][1] for i in range(1, length(V)) ]
    edge_component_universe = Vector{Symbol}([:pipe, :compressor, :control_valve, :valve, :short_pipe, :resistor, :loss_resistor])
    E = []
    for comp in edge_component_universe
        append!(E, get(ref(ss), comp, []))
    end
    g = SimpleGraph(length(V)) 

    for (i, edge) in E 
        fr = edge["fr_node"]
        to = edge["to_node"]
        if has_edge(g, new_node_from_old[fr], new_node_from_old[to])
            @warn "MORE THAN ONE EDGE BETWEEN NODES $(fr) AND $(to)"
            continue
        end
        add_edge!(g, new_node_from_old[fr], new_node_from_old[to])
    end 
    num_connected_components = length(connected_components(g))
    @assert (num_connected_components == 1) "Disconnected network ! $num_connected_components connected components found!"
    return g, new_node_from_old, old_node_from_new
end

function choose_articulation_point(g::SimpleGraph, vertex_list::Vector{Int64}, slack_nodes::Vector{Int64}, allow_slack_node_partitioning::Bool)::Union{Int64, Nothing}
    g_sub, vmap = induced_subgraph(g, vertex_list)
    if isempty(articulation(g_sub))
        return nothing
    end
    
    P = sort(articulation(g_sub), by = x->degree(g_sub, x), rev=true) #max degree first
    if allow_slack_node_partitioning == true
        return vmap[P[1]]
    end
    
    if allow_slack_node_partitioning == false
        P_new = [i for i in P if vmap[i]  ∉ slack_nodes]
        if isempty(P_new)
            return nothing
        else
            return vmap[P_new[1]]
        end
    end
    
end

function find_connected_components(g::SimpleGraph, vertex_list::Vector{Int64}, cut_vertex::Int64)::Vector{Vector}
    vertex_list1 = setdiff(vertex_list, [cut_vertex])
    g_cut, vmap = induced_subgraph(g, vertex_list1)
    C_list = connected_components(g_cut)
    for k = 1: length(C_list)
        C_list[k] = [vmap[v] for v in C_list[k]]
        push!(C_list[k], cut_vertex)
    end
    return C_list
end

function test_partitions_put_slack_network_first!(completed_sns::Vector{Vector{Int64}}, interface_nodes::Vector{Int64}, slack_nodes::Vector{Int64})::Vector{Int64}

    sort!(completed_sns, by = C->length(intersect(C, slack_nodes)), rev=true) #slack networks first
    slack_network_ids = []

    # identify slack networks
    for slack_id in slack_nodes
        for i = 1:length(completed_sns)
            if slack_id in completed_sns[i] && !(i  in slack_network_ids)
                push!(slack_network_ids, i)
            end
        end
    end

    return slack_network_ids

end

function all_interface_nodes_are_slack(interface_nodes::Vector{Int64}, slack_nodes::Vector{Int64})::Bool
    
    for id in  interface_nodes
        if !(id in slack_nodes)
            return false
            break
        end
    end
    return true
end



function create_partition_with_cut_points(ss::SteadySimulator; max_nodes::Int64 = 10, allow_slack_node_partitioning::Bool=false, write_to_file = false, filepath = "partition-expt.json")::Dict{String, Any}

    data = Dict{String,Any}(
        "num_partitions" => 0, 
        "slack_nodes" => [i for (i, node) in ref(ss, :node) if node["is_slack"] == 1]
    )

    completed_sns = Vector{Vector{Int64}}()
    other_sns = Vector{Vector{Int64}}()
    interface_nodes = Vector{Int64}()
    g, new_node_from_old, old_node_from_new = create_graph(ss);
    slack_nodes = [new_node_from_old[v] for v in data["slack_nodes"]]
    vertex_list = collect(vertices(g))
    count = 1
    while true
        c  = choose_articulation_point(g, vertex_list, slack_nodes, allow_slack_node_partitioning);
        if isnothing(c)
            push!(completed_sns, vertex_list)
        else
            push!(interface_nodes, c)
            C_list = find_connected_components(g, vertex_list, c)
            append!(other_sns, C_list)
        end
        sort!(other_sns, by = C->length(C), rev=true) #longest first

        if isempty(other_sns)
            @info "Cannot partition further!"
            break
        end

        if length(other_sns[1]) < max_nodes
            @info "Partitioning complete!" 
            append!(completed_sns, other_sns)
            break
        end
        
        vertex_list = splice!(other_sns, 1)
        count += 1
    end

    # put slack network first
    slack_network_ids = test_partitions_put_slack_network_first!(completed_sns, interface_nodes, slack_nodes)

    if all_interface_nodes_are_slack(interface_nodes, slack_nodes)
        @info "All interface nodes are slack nodes, no need for flow solve on block cut tree."
        data["block_cut_tree_solve"] = false
    else
        if length(intersect(completed_sns[1], slack_nodes)) == length(slack_nodes) #
            @info "At least one partition contains all slack nodes!"
            data["block_cut_tree_solve"] = true
        else
            @error "No subnetwork contains all slack nodes, invalid partition!"
            return Dict{String, Any}()
        
        end
    end

    data["slack_network_ids"] = Vector{Int64}(slack_network_ids)

    #back to old numbering
    for k = 1 : length(completed_sns)
        data[string(k)] = Vector{Any}
        data[string(k)] = [old_node_from_new[v] for v in completed_sns[k]]
    end

    data["num_partitions"] = length(completed_sns)
    interface_nodes = [old_node_from_new[v] for v in interface_nodes]
    data["interface_nodes"] = interface_nodes

    if write_to_file == true
        open(filepath, "w") do f 
            JSON.print(f, data)
        end
        @info "Partition data saved to $filepath"
    end 
    return data
end


function load_partition_data(data::Dict{String, Any})::Dict{Any, Any}

    partition = Dict{Any, Any}()
    num_partitions = data["num_partitions"]
    partition["num_partitions"] = num_partitions
    partition["num_interfaces"] = length(data["interface_nodes"])
    partition["interface_nodes"] = Vector{Int64}(data["interface_nodes"])
    partition["slack_network_ids"] = Vector{Int64}(data["slack_network_ids"])
    partition["slack_nodes"]  = Vector{Int64}(data["slack_nodes"])
    partition["block_cut_tree_solve"] = data["block_cut_tree_solve"]

    @assert num_partitions > 1
    
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

    first_slack_network = partition["slack_network_ids"][1]
    if partition["block_cut_tree_solve"] == true
        @assert length(intersect(partition["slack_nodes"], partition[first_slack_network]["node_list"])) == length(partition["slack_nodes"])
    end

    
    order = Vector{Int64}(partition["slack_network_ids"])


    stage = 1
    partition["level"] = Dict{Int, Vector{Int64}}()
    partition["level"][stage] =  Vector{Int64}(partition["slack_network_ids"])

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
                if ssp_array[sn_id_j].ref[:is_pressure_node][node_id] == true
                    if isnan(ssp_array[sn_id_i].ref[:node][node_id]["pressure"])
                        println("subnetwork_id: ", sn_id_i, " node_id: ", node_id)
                        @warn("NaN found !")
                    end
                    ssp_array[sn_id_j].ref[:node][node_id]["pressure"] =  ssp_array[sn_id_i].ref[:node][node_id]["pressure"]
                else
                    if isnan(ssp_array[sn_id_i].ref[:node][node_id]["potential"])
                        println("subnetwork_id: ", sn_id_i, " node_id: ", node_id)
                        @warn("NaN found !")
                    end
                    ssp_array[sn_id_j].ref[:node][node_id]["potential"] =  ssp_array[sn_id_i].ref[:node][node_id]["potential"]
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
                if ss.ref[:is_pressure_node][id] == true
                    x_dofs[ss.ref[comp][id]["dof"]] = ssp_array[sn_id].ref[comp][id]["pressure"]
                else
                    x_dofs[ss.ref[comp][id]["dof"]] = ssp_array[sn_id].ref[comp][id]["potential"]
                end
            else
                x_dofs[ss.ref[comp][id]["dof"]] = ssp_array[sn_id].ref[comp][id]["flow"]
            end
        end
    end

    return x_dofs
end
