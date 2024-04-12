
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
        x_guess = rand(n)
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



function run_partitioned_ss(filepath::AbstractString, ss::SteadySimulator; show_trace_flag::Bool=false)::Vector{Float64}

    partition = create_partition(filepath)

    num_partition = partition["num_partitions"]
    
    ssp_array = Vector{SteadySimulator}()
    for i = 1 : num_partition
            push!(ssp_array, initialize_simulator_subnetwork(ss, partition[i]["node_list"]))
    end
    println("Initialized steady state simulator...")

    designate_interface_nodes_as_slack!(ssp_array, partition)

    # at this point solve flow problem on block tree for exact interface transfers assuming single slack
    flow_solve_on_block_cut_tree!(ssp_array, partition)


    df_array = Vector{OnceDifferentiable}()  # datatype OnceDifferentiable not found here
    for i = 1: num_partition
            push!(df_array, prepare_for_nonlin_solve!(ssp_array[i]))
    end
    println("Initialized nonlinear solve...")

    println("Propagating  subnetwork solution ...")
    
    for level = 1: partition["num_level"]
        for sn_id in partition["level"][level]
            solver = solve_on_network!(ssp_array[sn_id], df_array[sn_id], show_trace_flag=show_trace_flag)
        end
        update_interface_potentials_of_nbrs!(ssp_array, partition, level)
    end

    x_dof = combine_subnetwork_solutions(ss, ssp_array)

    println("Completed")


    return x_dof

end

# function _create_initial_guess_dof!(ss::SteadySimulator)::Array
#     ndofs = length(ref(ss, :dof))
#     x_guess = 0.5 * ones(Float64, ndofs) 
#     dofs_updated = 0

#     components = [:node, :pipe, :compressor, 
#         :control_valve, :valve, 
#         :resistor, :loss_resistor, :short_pipe]

#     for component in components 
#         for (i, val) in get(ss.initial_guess, component, [])
#             x_guess[ref(ss, component, i, "dof")] = val 
#             dofs_updated += 1
#         end 
#     end 
#     return x_guess
# end

function sum_q(ssp_array::Vector{SteadySimulator},  partition::Dict{Any, Any}, vertex_id::Int64)::Float64
    c = 0.0
    for i in partition[vertex_id]["node_list"]
        if ssp_array[vertex_id].ref[:node][i]["is_slack"] != 1
            c = c + ssp_array[vertex_id].ref[:node][i]["withdrawal"]
        end
    end
    return c
end


function flow_solve_on_block_cut_tree!(ssp_array::Vector{SteadySimulator},  partition::Dict{Any, Any})

    num_edges = partition["num_partitions"] + partition["num_interfaces"] - 1  # since it is a tree

    #convention - flow is towards interface (so flow is always withdrawal from network)
    A = zeros(num_edges, num_edges)
    b = zeros(num_edges)
    eq_no = 1
    for i in keys(partition["vertex_to_global_map"])
        if i == 1 #slack network
            continue
        end
        for edge_dof in partition["vertex_to_global_map"][i]
            if i in partition["interface_nodes"]
                A[eq_no, edge_dof] = 1 
                b[eq_no] = 0.0
            else
                A[eq_no, edge_dof] = -1
                b[eq_no] = sum_q(ssp_array, partition, i) #withdrawal
            end
        end
        eq_no += 1
    end
    f = A \ b # net inflow in terms of A = withdrawal (outflow) in b

    for i = 1 : num_edges
        (n1, node_id) = partition["global_to_local"][i]
        if ssp_array[n1].ref[:node][node_id]["is_slack"] != 1
            ssp_array[n1].ref[:node][node_id]["transfer"] = f[i] # f is a withdrawal from network
        end
    end

    return
end