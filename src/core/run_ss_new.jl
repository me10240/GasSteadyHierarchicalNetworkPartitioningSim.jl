
using LinearAlgebra
using NLSolversBase


function prepare_for_nonlin_solve!(ss::SteadySimulator)::OnceDifferentiable
    

    residual_fun! = (r_dof, x_dof) -> assemble_residual!(ss, x_dof, r_dof)
    Jacobian_fun! = (J_dof, x_dof) -> assemble_mat!(ss, x_dof, J_dof)
    n = length(ref(ss, :dof))
    J0 = spzeros(n, n)
    assemble_mat!(ss, rand(n), J0)
    df = OnceDifferentiable(residual_fun!, Jacobian_fun!, rand(n), rand(n), J0)

    return df
end

function solve_on_network!(ss::SteadySimulator, df::OnceDifferentiable; x_guess::Vector=Vector{Float64}(), 
    method::Symbol=:newton,
    iteration_limit::Int64=2000, 
    show_trace_flag::Bool=false,
    kwargs...)::SolverReturn

    
    if isempty(x_guess)
        n = length(ref(ss, :dof))
        x_guess = ones(n)
    end

    time = @elapsed soln = nlsolve(df, x_guess; method = method, iterations = iteration_limit, show_trace=show_trace_flag, kwargs...)


    convergence_state = converged(soln)

    if convergence_state == false
        return SolverReturn(nl_solve_failure, 
            soln.iterations, 
            soln.residual_norm, 
            time, soln.zero, 
            Int[], Int[], Int[])
    end

    sol_return = update_solution_fields_in_ref!(ss, soln.zero)



    return SolverReturn(unique_physical_solution, 
        soln.iterations, 
        soln.residual_norm, 
        time, soln.zero, 
        sol_return[:compressors_with_neg_flow], 
        sol_return[:nodes_with_neg_potential],
        sol_return[:nodes_with_pressure_not_in_domain])
end



function run_partitioned_ss(filepath::AbstractString, ss::SteadySimulator; eos::Symbol=:ideal, show_trace_flag::Bool=false, iteration_limit::Int=200, method::Symbol=:newton, cond_number::Bool=true)::Tuple{Vector{Float64}, Vector{Float64}} 

    partition = create_partition(filepath)

    # assume if multiple  slack nodes, then one (first) partition must contain all of them. Others can contain one or more.
    num_partition = partition["num_partitions"]

    cond_number_array = Vector{Float64}()
    ssp_array = Vector{SteadySimulator}()
    for i = 1 : num_partition
            push!(ssp_array, initialize_simulator_subnetwork(ss, partition[i]["node_list"], eos))
    end

    println("Initialized steady state simulator...")

    designate_interface_nodes_as_slack!(ssp_array, partition)

    # at this point solve flow problem on block tree for exact interface transfers assuming single slack
    flow_solve_on_block_cut_tree!(ssp_array, partition)

    println("Propagating  subnetwork solution ...")
    
    for level = 1: partition["num_level"]
        println("Solving level $level subnetworks")
        for sn_id in partition["level"][level]
            @show sn_id

            if length(ssp_array[sn_id].ref[:node]) == 2
                solver1 = solve_edge!(ssp_array[sn_id])
                # println(sn_id, " ", solver1.iterations, "  ", solver1.residual_norm)
                # continue
            end

            df = prepare_for_nonlin_solve!(ssp_array[sn_id])
            if cond_number == true
                push!(cond_number_array, cond(gradient(df), 1) )
            end
            solver = solve_on_network!(ssp_array[sn_id], df, show_trace_flag=show_trace_flag, iteration_limit=iteration_limit, method=method)
            println(sn_id, " ", solver.iterations, " ", solver.residual_norm)
            # if sn_id in [6, 16, 10, 11, 20, 3]
            # if length(ssp_array[sn_id].ref[:node]) == 2
            #     println(solver1.solution, "\n", solver.solution, "\n", solver1.solution - solver.solution)
            # end


            # if length(ssp_array[sn_id].ref[:node]) == 2
            #     println(norm(x_dof - solver.solution))
            # end
        end
        update_interface_potentials_of_nbrs!(ssp_array, partition, level)
    end

    x_dof = combine_subnetwork_solutions(ss, ssp_array)

    println("Completed")


    return x_dof, cond_number_array

end


function _sum_q(ssp_array::Vector{SteadySimulator},  partition::Dict{Any, Any}, vertex_id::Int64)::Float64
    c = 0.0
    for i in partition[vertex_id]["node_list"]
        if ssp_array[vertex_id].ref[:node][i]["is_slack"] != 1
            c = c + ssp_array[vertex_id].ref[:node][i]["withdrawal"]
        end
    end
    return c
end

function _assemble_system(ssp_array::Vector{SteadySimulator},  partition::Dict{Any, Any})::Tuple{Matrix{Float64}, Vector{Float64}}
    
    num_edges = partition["num_partitions"] + partition["num_interfaces"] - 1  # since it is a tree

    #convention - flow is towards interface (so flow is always withdrawal from network)
    A = zeros(num_edges, num_edges)
    b = zeros(num_edges)
    for i = 1 : partition["num_partitions"]
        if i == 1 #slack network
            continue
        end
        node_i = "N-$i"
        for edge_dof in partition["vertex_to_global_map"][node_i]
            eq_no = i - 1
            A[eq_no, edge_dof] = -1
            b[eq_no] = _sum_q(ssp_array, partition, i) #withdrawal
        end
    end

    for j = 1: length(partition["interface_nodes"])
        node_j = partition["interface_nodes"][j]
        for edge_dof in partition["vertex_to_global_map"][node_j]
            eq_no = j + partition["num_partitions"] - 1
            A[eq_no, edge_dof] = 1 
            b[eq_no] = 0.0
        end
    end

    return A, b
    
end

function _assemble_transfers_into_network!(ssp_array::Vector{SteadySimulator},  partition::Dict{Any, Any}, f::Vector{Float64})
    num_edges = partition["num_partitions"] + partition["num_interfaces"] - 1  # since it is a tree
    for i = 1 : num_edges
        (n1, node_id) = partition["global_to_local"][i]
        n_index= parse(Int64, split(n1,'-')[2])
        if ssp_array[n_index].ref[:node][node_id]["is_slack"] != 1 # slack withdrawal will be NaN
            ssp_array[n_index].ref[:node][node_id]["transfer"] = f[i] # f is a withdrawal from network
        end
    end
    return
end


function flow_solve_on_block_cut_tree!(ssp_array::Vector{SteadySimulator},  partition::Dict{Any, Any})

    A, b = _assemble_system(ssp_array, partition)
    f = A \ b # net inflow in terms of A = withdrawal (outflow) in b
    _assemble_transfers_into_network!(ssp_array, partition, f)

    return
end

function solve_edge!(ss::SteadySimulator)::SolverReturn

    x_dof = zeros(3)
    if length( keys( get( ref(ss), :pipe, [] ) ) ) == 1
        p_key = collect(keys(ref(ss, :pipe)))[1]
        residual = _solve_pipe!(ss, p_key, x_dof)
    end

    if length( keys( get( ref(ss), :compressor, [] ) ) ) == 1
        c_key = collect(keys(ref(ss, :compressor)))[1]
        residual = _solve_compressor!(ss, c_key, x_dof)
    end
    
    sol_return = update_solution_fields_in_ref!(ss, x_dof)
        
    return SolverReturn(unique_physical_solution, 
    0.0, 
    residual, 
    0.0, x_dof, 
    sol_return[:compressors_with_neg_flow], 
    sol_return[:nodes_with_neg_potential],
    sol_return[:nodes_with_pressure_not_in_domain])
end

function _solve_pipe!(ss::SteadySimulator, p_key::Int64, x_dof::Vector{Float64})::Float64

    to_node = ref(ss, :pipe, p_key, "to_node") 
    fr_node = ref(ss, :pipe, p_key, "fr_node")

    if ref(ss, :node, to_node, "is_slack") == 1
        val = get(ref(ss, :node,  fr_node), "transfer", 0) +  ref(ss, :node,  fr_node, "withdrawal")
        flow = -val
    else
        val = get(ref(ss, :node,  to_node), "transfer", 0) + ref(ss, :node,  to_node, "withdrawal")
        flow = val
    end

    x_dof[ref(ss, :pipe, p_key, "dof")] = flow
    
    c = nominal_values(ss, :mach_num)^2 / nominal_values(ss, :euler_num) 
    pipe = ref(ss, :pipe, p_key)
    resistance = pipe["friction_factor"] * pipe["length"] * c / (2 * pipe["diameter"] * pipe["area"]^2)

    if ref(ss, :node, to_node, "is_slack") == 1
        if ref(ss, :is_pressure_node, to_node) == true
            x_dof[ref(ss, :node, to_node, "dof")] = ref(ss, :node, to_node, "pressure")
            pi_to = get_potential(ss, x_dof[ref(ss, :node, to_node, "dof")]) 
        else
            x_dof[ref(ss, :node, to_node, "dof")] = ref(ss, :node, to_node, "potential")
            pi_to =  x_dof[ref(ss, :node, to_node, "dof")]
        end
        pi_fr = pi_to + flow * abs(flow) * resistance
        x_dof[ref(ss, :node, fr_node, "dof")] = ref(ss, :is_pressure_node, fr_node) ? invert_positive_potential(ss, pi_fr) :  pi_fr
    else
        if ref(ss, :is_pressure_node, fr_node) == true
            x_dof[ref(ss, :node, fr_node, "dof")] = ref(ss, :node, fr_node, "pressure")
            pi_fr = get_potential(ss, x_dof[ref(ss, :node, fr_node, "dof")]) 
        else
            x_dof[ref(ss, :node, fr_node, "dof")] = ref(ss, :node, fr_node, "potential")
            pi_fr =  x_dof[ref(ss, :node, fr_node, "dof")]
        end
        pi_to = pi_fr - flow * abs(flow) * resistance
        x_dof[ref(ss, :node, to_node, "dof")] = ref(ss, :is_pressure_node, to_node) ? invert_positive_potential(ss, pi_to) :  pi_to
    end
    residual = abs(pi_fr - pi_to - flow * abs(flow) * resistance)

    return residual
end

function _solve_compressor!(ss::SteadySimulator, c_key::Int64, x_dof::Vector{Float64})::Float64

    to_node = ref(ss, :compressor, c_key, "to_node") 
    fr_node = ref(ss, :compressor, c_key, "fr_node") 

    if ref(ss, :node, to_node, "is_slack") == 1
        val = get(ref(ss, :node,  fr_node), "transfer", 0) + ref(ss, :node,  fr_node, "withdrawal")
        flow = -val
    else
        val = get(ref(ss, :node,  to_node), "transfer", 0) + ref(ss, :node,  to_node, "withdrawal")
        flow = val
    end
    x_dof[ref(ss, :compressor, c_key, "dof")] = flow

    comp = ref(ss, :compressor, c_key)
    cmpr_val = comp["c_ratio"] 

    c0, c1, c2, c3 = ss.potential_ratio_approx
    cmpr_val_expr =  c0  +  c1 * cmpr_val + c2 * cmpr_val^2 + c3 * cmpr_val^3

    is_pressure_eq = ref(ss, :is_pressure_node, fr_node) || ref(ss, :is_pressure_node, to_node)
    if is_pressure_eq == true
        val = cmpr_val
        if ref(ss, :node, to_node, "is_slack") == 1
            x_dof[ref(ss, :node, to_node, "dof")] = ref(ss, :node, to_node, "pressure")
            x_dof[ref(ss, :node, fr_node, "dof")] = x_dof[ref(ss, :node, to_node, "dof")] / val
        else
            x_dof[ref(ss, :node, fr_node, "dof")] = ref(ss, :node, fr_node, "pressure")
            x_dof[ref(ss, :node, to_node, "dof")] = x_dof[ref(ss, :node, fr_node, "dof")] *  val
        end
        
    else
        val = cmpr_val_expr
        if ref(ss, :node, to_node, "is_slack") == 1
            x_dof[ref(ss, :node, to_node, "dof")] = ref(ss, :node, to_node, "potential")
            x_dof[ref(ss, :node, fr_node, "dof")] = x_dof[ref(ss, :node, to_node, "dof")] / val
        else
            x_dof[ref(ss, :node, fr_node, "dof")] = ref(ss, :node, fr_node, "potential")
            x_dof[ref(ss, :node, to_node, "dof")] = x_dof[ref(ss, :node, fr_node, "dof")] *  val
        end
    end

    residual = abs(x_dof[ref(ss, :node, fr_node, "dof")] * val - x_dof[ref(ss, :node, to_node, "dof")])

    return residual
end