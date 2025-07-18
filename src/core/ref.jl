function _add_components_to_ref!(ref::Dict{Symbol,Any}, data::Dict{String,Any}, bc::Dict{Symbol,Any})

    ref[:control_vars] = Dict{Int64, Any}()
    num_control_vars = 0

    for (i, node) in get(data, "nodes", [])
        name = :node
        (!haskey(ref, name)) && (ref[name] = Dict())
       
        id = parse(Int64, i)
        ref[name][id] = Dict()
        @assert id == (haskey(node, "id") ? node["id"] : node["node_id"])
        
        ref[name][id]["id"] = id
        ref[name][id]["is_slack"] = node["slack_bool"]
        ref[name][id]["pressure"] = NaN
        ref[name][id]["density"] = NaN 
        ref[name][id]["withdrawal"] = NaN
        ref[name][id]["potential"] = NaN

    

        if haskey(bc[:node], id)
            num_control_vars += 1
        end
        

            

        if ref[:node][id]["is_slack"] == 1 #slack
            ref[:control_vars][num_control_vars] = (:node, id, "pressure") ##
            ref[:node][id]["withdrawal"] = NaN
            ref[:node][id]["pressure"] = bc[:node][id]["val"]
        else
            ref[:node][id]["withdrawal"] = haskey(bc[:node], id) ? bc[:node][id]["val"] : 0.0 #added
            ref[:control_vars][num_control_vars] = (:node, id, "injection") ##
        end

    end

    for (i, pipe) in get(data, "pipes", [])
        name = :pipe
        (!haskey(ref, name)) && (ref[name] = Dict())
        id = parse(Int64, i)
        ref[name][id] = Dict()
        @assert id == (haskey(pipe, "id") ? pipe["id"] : pipe["pipe_id"])

        ref[name][id]["id"] = id
        ref[name][id]["fr_node"] = (haskey(pipe, "fr_node") ? pipe["fr_node"] : pipe["from_node"])
        ref[name][id]["to_node"] = pipe["to_node"]
        ref[name][id]["diameter"] = pipe["diameter"]
        ref[name][id]["area"] = pipe["area"]
        ref[name][id]["length"] = pipe["length"]
        ref[name][id]["friction_factor"] = pipe["friction_factor"]
        ref[name][id]["flow"] = NaN
    end

    for (i, compressor) in get(data, "compressors", [])
        name = :compressor
        (!haskey(ref, name)) && (ref[name] = Dict())
        id = parse(Int64, i)
        ref[name][id] = Dict()
        @assert id == (haskey(compressor, "id") ? compressor["id"] : compressor["comp_id"])

        ref[name][id]["id"] = id
        ref[name][id]["to_node"] = compressor["to_node"]
        ref[name][id]["fr_node"] = (haskey(compressor, "fr_node") ? compressor["fr_node"] : compressor["from_node"])
        ref[name][id]["control_type"] = unknown_control
        ref[name][id]["c_ratio"] = bc[name][id]["val"]
        ref[name][id]["discharge_pressure"] = NaN
        ref[name][id]["flow"] = NaN
        num_control_vars += 1
        ref[:control_vars][num_control_vars] = (:compressor, id, "c_ratio") ##
    end

    for (i, control_valve) in get(data, "control_valves", [])
        (parse(Int64, i) in bc[:control_valve_status][:off]) && (continue)
        name = :control_valve
        (!haskey(ref, name)) && (ref[name] = Dict())
        id = parse(Int64, i)
        ref[name][id] = Dict()
        @assert id == (haskey(control_valve, "id") ? control_valve["id"] : control_valve["control_valve_id"])

        ref[name][id]["id"] = id
        ref[name][id]["to_node"] = control_valve["to_node"]
        ref[name][id]["fr_node"] = (haskey(control_valve, "fr_node") ? control_valve["fr_node"] : control_valve["from_node"])
        ref[name][id]["control_type"] = unknown_control
        ref[name][id]["c_ratio"] = bc[name][id]["val"]
        ref[name][id]["discharge_pressure"] = NaN
        ref[name][id]["flow"] = NaN
    end 

    for (i, valve) in get(data, "valves", []) 
        (parse(Int64, i) in bc[:valve_status][:off]) && (continue)
        name = :valve 
        (!haskey(ref, name)) && (ref[name] = Dict()) 
        id = parse(Int64, i)
        ref[name][id] = Dict()
        @assert id == (haskey(valve, "id") ? valve["id"] : valve["valve_id"])

        ref[name][id]["id"] = id
        ref[name][id]["to_node"] = valve["to_node"]
        ref[name][id]["fr_node"] = (haskey(valve, "fr_node") ? valve["fr_node"] : valve["from_node"])
        ref[name][id]["flow"] = NaN
    end 

    for (i, resistor) in get(data, "resistors", [])
        name = :resistor
        (!haskey(ref, name)) && (ref[name] = Dict())
        id = parse(Int64, i)
        ref[name][id] = Dict()
        @assert id == (haskey(resistor, "id") ? resistor["id"] : resistor["resistor_id"])

        ref[name][id]["id"] = id
        ref[name][id]["to_node"] = resistor["to_node"]
        ref[name][id]["fr_node"] = (haskey(resistor, "fr_node") ? resistor["fr_node"] : resistor["from_node"])
        ref[name][id]["drag"] = resistor["drag"]
        ref[name][id]["diameter"] = resistor["drag"]
        ref[name][id]["flow"] = NaN
    end 

    for (i, loss_resistor) in get(data, "loss_resistors", [])
        name = :loss_resistor
        (!haskey(ref, name)) && (ref[name] = Dict())
        id = parse(Int64, i)
        ref[name][id] = Dict()
        @assert id == (haskey(loss_resistor, "id") ? loss_resistor["id"] : loss_resistor["loss_resistor_id"])

        ref[name][id]["id"] = id
        ref[name][id]["fr_node"] = (haskey(loss_resistor, "fr_node") ? loss_resistor["fr_node"] : loss_resistor["from_node"])
        ref[name][id]["to_node"] = loss_resistor["to_node"]
        ref[name][id]["pressure_drop"] = loss_resistor["p_loss"]
        ref[name][id]["flow"] = NaN
    end 

    for (i, pipe) in get(data, "short_pipes", [])
        name = :short_pipe
        (!haskey(ref, name)) && (ref[name] = Dict())
        id = parse(Int64, i)
        ref[name][id] = Dict()
         @assert id == (haskey(pipe, "id") ? pipe["id"] : pipe["short_pipe_id"])

        ref[name][id]["id"] = id
        ref[name][id]["fr_node"] = (haskey(pipe, "fr_node") ? pipe["fr_node"] : pipe["from_node"])
        ref[name][id]["to_node"] = pipe["to_node"]
        ref[name][id]["flow"] = NaN
    end 

    ref[:total_control_vars] = num_control_vars 
    return
end

function _add_index_info!(ref::Dict{Symbol, Any}, data::Dict{String, Any})
    dofid = 1
    ref[:dof] = Dict{Int64, Any}()
    
    for (i, node) in ref[:node]
        node["dof"] = dofid
        ref[:dof][dofid] = (:node, i)
        dofid += 1
    end

    edge_component_universe = Vector{Symbol}([:pipe, :compressor, :control_valve, :valve, :short_pipe, :resistor, :loss_resistor])
    edge_component_list = Vector{Symbol}()

    for key in edge_component_universe
        if key in keys(ref)
            push!(edge_component_list, key)
        end
    end

    for comp in edge_component_list
        for (i, comp_handle) in ref[comp]
            comp_handle["dof"] = dofid
            ref[:dof][dofid] = (comp, i)
            dofid += 1
        end
    end

    
    return
end

function _add_incident_dofs_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
    ref[:incoming_dofs] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )
    ref[:outgoing_dofs] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )

    edge_component_universe = Vector{Symbol}([:pipe, :compressor, :control_valve, :valve, :short_pipe, :resistor, :loss_resistor])
    edge_component_list = Vector{Symbol}()

    for key in edge_component_universe
        if key in keys(ref)
            push!(edge_component_list, key)
        end
    end

    for comp in edge_component_list
        for (_, comp_handle) in get(ref, comp, [])
            push!(ref[:incoming_dofs][comp_handle["to_node"]], comp_handle["dof"])
            push!(ref[:outgoing_dofs][comp_handle["fr_node"]], comp_handle["dof"])
        end
    end


    return
end

function _add_pipe_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
    ref[:incoming_pipes] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )
    ref[:outgoing_pipes] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )

    for (id, pipe) in ref[:pipe]
        push!(ref[:incoming_pipes][pipe["to_node"]], id)
        push!(ref[:outgoing_pipes][pipe["fr_node"]], id)
    end
    return
end

function _add_compressor_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
    ref[:incoming_compressors] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )
    ref[:outgoing_compressors] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )

    for (id, compressor) in get(ref, :compressor, [])
        push!(ref[:incoming_compressors][compressor["to_node"]], id)
        push!(ref[:outgoing_compressors][compressor["fr_node"]], id)
    end
    return
end

function _add_control_valve_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
    ref[:incoming_control_valves] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )
    ref[:outgoing_control_valves] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )

    for (id, control_valve) in get(ref, :control_valve, [])
        push!(ref[:incoming_control_valves][control_valve["to_node"]], id)
        push!(ref[:outgoing_control_valves][control_valve["fr_node"]], id)
    end
    return
end

function _add_valve_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
    ref[:incoming_valves] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )
    ref[:outgoing_valves] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )
    
    for (id, valve) in get(ref, :valve, [])
        push!(ref[:incoming_valves][valve["to_node"]], id)
        push!(ref[:outgoing_valves][valve["fr_node"]], id)
    end
    return
end

function _add_resistor_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
    ref[:incoming_resistors] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )
    ref[:outgoing_resistors] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )
    
    for (id, resistor) in get(ref, :resistor, [])
        push!(ref[:incoming_resistors][resistor["to_node"]], id)
        push!(ref[:outgoing_resistors][resistor["fr_node"]], id)
    end
    return
end

function _add_loss_resistor_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
    ref[:incoming_loss_resistors] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )
    ref[:outgoing_loss_resistors] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )
    
    for (id, loss_resistor) in get(ref, :loss_resistor, [])
        push!(ref[:incoming_loss_resistors][loss_resistor["to_node"]], id)
        push!(ref[:outgoing_loss_resistors][loss_resistor["fr_node"]], id)
    end
    return
end

function _add_short_pipe_info_at_nodes!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
    ref[:incoming_short_pipes] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )
    ref[:outgoing_short_pipes] = Dict{Int64, Vector{Int64}}(
        i => [] for i in keys(ref[:node])
    )

    for (id, short_pipe) in get(ref, :short_pipe, [])
        push!(ref[:incoming_short_pipes][short_pipe["to_node"]], id)
        push!(ref[:outgoing_short_pipes][short_pipe["fr_node"]], id)
    end
    return
end

function _add_pressure_node_flag!(ref::Dict{Symbol,Any}, data::Dict{String,Any})
    ref[:is_pressure_node] = Dict{Int64, Bool}(
        i => false for i in keys(ref[:node])
    )

    for (_, compressor) in get(ref, :compressor, [])
        ref[:is_pressure_node][compressor["fr_node"]] = true 
        ref[:is_pressure_node][compressor["to_node"]] = true
    end
    return 
end 

function _update_node_flag!(ref::Dict{Symbol,Any})
    for i in keys(ref[:is_pressure_node])
        ref[:is_pressure_node][i] = false
    end 
    return
end 


function build_ref(data::Dict{String,Any}, bc::Dict{Symbol,Any};
    ref_extensions=[])::Dict{Symbol,Any}

    ref = Dict{Symbol,Any}()
    _add_components_to_ref!(ref, data, bc)

    for extension in ref_extensions
        extension(ref, data)
    end

    return ref
end

function _add_components_to_subnetwork_ref!(ref::Dict{Symbol,Any},  full_ref::Dict{Symbol,Any}, node_list::Vector)
    
    ref[:node] = Dict()

    for id in node_list
        ref[:node][id] = deepcopy( full_ref[:node][id] )
    end

    edge_component_universe = Vector{Symbol}([:pipe, :compressor, :control_valve, :valve, :short_pipe, :resistor, :loss_resistor])
    edge_component_list = Vector{Symbol}()

    for key in edge_component_universe
        if key in keys(full_ref)
            ref[key] = Dict()
            push!(edge_component_list, key)
        end
    end

    for comp in edge_component_list
        for (id, comp_handle) in full_ref[comp]
            if (full_ref[comp][id]["fr_node"] in node_list) && (full_ref[comp][id]["to_node"] in node_list)
                ref[comp][id] = deepcopy( full_ref[comp][id])
            end
        end
    end
    return

end


function build_subnetwork_ref(full_ref::Dict{Symbol,Any}, node_list::Vector; ref_extensions=[])::Dict{Symbol,Any}

    ref = Dict{Symbol,Any}()
     _add_components_to_subnetwork_ref!(ref,  full_ref, node_list)
    
    # set pressure_node flag from full_ref so that node flag does not change across subnetworks
    ref[:is_pressure_node] = Dict{Int64, Bool}(
        i => full_ref[:is_pressure_node][i] for i in keys(ref[:node])
    )
 
    var = Dict{String, Any}()
    for extension in ref_extensions
        extension(ref, var)
    end

    return ref
end