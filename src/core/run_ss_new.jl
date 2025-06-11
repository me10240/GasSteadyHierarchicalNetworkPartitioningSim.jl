
using LinearAlgebra
using NLSolversBase


function interface_residual_jacobian!(ssp_array::Vector{SteadySimulator}, df_array::Vector{OnceDifferentiable}, partition::Dict{Any, Any}, x_dof::AbstractArray, r_dof::AbstractArray, J_dof::AbstractArray)

    
    update_interface_slack_dofs!(ssp_array, partition, x_dof)
    for id = 1 : partition["num_partitions"]
        @assert length(ssp_array[id].ref[:node]) > 2
        solver = solve_on_network!(ssp_array[id], df_array[id], show_trace_flag=false, iteration_limit=100, method=:newton, sensitivity_nodes= partition[id]["interface"])
        update_transfers!(ssp_array[id], partition[id])
        update_transfer_sensitivities!(solver.sensitivity_mat, partition[id])
    end

    for node_id  in partition["interface_nodes"]
        r_dof[partition["vertex_to_global"][node_id]] = -partition["interface_withdrawals"][node_id]
    end


    for i = 1 : partition["num_partitions"]
        for node_j in partition[i]["interface"]
            r_dof[partition["vertex_to_global"][node_j]] += partition[i]["transfer"][node_j]
            for node_k in partition[i]["interface"]
                J_dof[partition["vertex_to_global"][node_j], partition["vertex_to_global"][node_k]] += partition[i]["transfer_sensitivity_mat"][node_j][node_k]
            end
        end
    end

    return
end

function prepare_for_partition_interface_solve!(ssp_array::Vector{SteadySimulator}, df_array::Vector{OnceDifferentiable}, partition::Dict{Any, Any})::OnceDifferentiable
    
    my_fun! = (r_dof, J_dof, x_dof) -> interface_residual_jacobian!(ssp_array, df_array, partition, x_dof, r_dof, J_dof)

    n = partition["num_interfaces"]
    r_dof = zeros(n)
    J_dof = zeros(n, n)

    residual_fun! = (r_dof, x_dof) -> my_fun!(r_dof, J_dof, x_dof)
    Jacobian_fun! = (J_dof, x_dof) -> my_fun!(r_dof, J_dof, x_dof)
    
    my_fun!(r_dof, J_dof, rand(n))
    df = OnceDifferentiable(residual_fun!, Jacobian_fun!, rand(n), rand(n), J_dof)

    return df
end

function solve_on_partition_interface!(partition::Dict{Any, Any}, df::OnceDifferentiable; x_guess::Vector=Vector{Float64}(), 
    method::Symbol=:newton,
    iteration_limit::Int64=2000, 
    show_trace_flag::Bool=false,
    kwargs...)::SolverReturnPartitionInterface

    n = partition["num_interfaces"]
    if isempty(x_guess)
        x_guess = ones(n)
    end

    time = @elapsed soln = nlsolve(df, x_guess; method = method, iterations = iteration_limit, show_trace=show_trace_flag, kwargs...)

    convergence_state = converged(soln)

    if convergence_state == false
        return SolverReturnPartitionInterface(nl_solve_failure, 
            soln.iterations, 
            soln.residual_norm, 
            time, soln.zero)
    end

    return SolverReturnPartitionInterface(success, 
        soln.iterations, 
        soln.residual_norm, 
        time, soln.zero)
end

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
    sensitivity_nodes::Vector{Int64}=Vector{Int64}(),
    kwargs...)::SolverReturn

    
    n = length(ref(ss, :dof))
    if isempty(x_guess)
        x_guess = ones(n)
    end

    time = @elapsed soln = nlsolve(df, x_guess; method = method, iterations = iteration_limit, show_trace=show_trace_flag, kwargs...)


    convergence_state = converged(soln)

    if convergence_state == false
        return SolverReturn(nl_solve_failure, 
            soln.iterations, 
            soln.residual_norm, 
            time, soln.zero, Dict{Int64, Any}(),
            Int[], Int[], Int[])
    end

    sol_return = update_solution_fields_in_ref!(ss, soln.zero)

    sensitivity_dict = Dict{Int64, Any}()
    if !isempty(sensitivity_nodes)
        Jsol = spzeros(n, n)
        assemble_mat!(ss, soln.zero, Jsol) # Jacobian at the solution
        RHS = zeros(n, length(sensitivity_nodes))
        # create dict to hold sensitivity info
        initialize_sensitivity_dict!(sensitivity_dict, sensitivity_nodes)
        form_sensitivity_RHS!(ss, RHS, sensitivity_nodes)
        calculate_sensitivities!(Jsol, RHS)
        update_withdrawal_sensitivities!(ss, RHS, sensitivity_dict, sensitivity_nodes)
    end

    return SolverReturn(success, 
        soln.iterations, 
        soln.residual_norm, 
        time, soln.zero, sensitivity_dict,
        sol_return[:compressors_with_neg_flow], 
        sol_return[:nodes_with_neg_potential],
        sol_return[:nodes_with_pressure_not_in_domain])
end



function run_partitioned_ss(filepath::AbstractString, ss::SteadySimulator; eos::Symbol=:ideal, show_trace_flag::Bool=false, iteration_limit::Int=200, method::Symbol=:newton)::Vector{Float64} 

    partition = create_partition(filepath)
    set_interface_withdrawals!(ss, partition)


    num_partition = partition["num_partitions"]

    ssp_array = Vector{SteadySimulator}()
    for i = 1 : num_partition
            push!(ssp_array, initialize_simulator_subnetwork(ss, partition[i]["node_list"], eos))
    end

    println("Initialized steady state simulator...")

    set_interface_nodes_as_slack!(ssp_array, partition)

    df_array = Vector{OnceDifferentiable}()
    for i = 1 : num_partition
        push!(df_array, prepare_for_nonlin_solve!(ssp_array[i]))
    end
    # x_dof = ones(partition["num_interfaces"])
    # x_dof = [2.837016258123215, 2.5775965089419297] # these are exact values
    # x_dof = [2.83, 2.57]
    x_dof = [1.0, 2.0]
    
    df = prepare_for_partition_interface_solve!(ssp_array, df_array, partition)
    solver = solve_on_partition_interface!(partition, df; x_guess=x_dof, 
    method=method, iteration_limit=iteration_limit, 
    show_trace_flag=show_trace_flag)
    println(solver)

    x_dof = combine_subnetwork_solutions(ss, ssp_array)

    println("Completed")

    return x_dof
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
    # change this  bcos here pressure/pot is given, not flow/transfer
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