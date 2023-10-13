
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

    # if length(partition["level"][4]) != 0
    #         error("This partitioning is not acceptable")
    # end

    num_partition = length(partition["subnetworks"])

    
    ssp_array = Vector{SteadySimulator}()
    for i = 1 : num_partition
            push!(ssp_array, initialize_simulator_subnetwork(ss, partition["subnetworks"][i]["node_list"]))
    end
    println("Initialized steady state simulator...")


    create_interface_transfers!(ssp_array,  partition)
    designate_interface_nodes_as_slack!(ssp_array, partition)


    df_array = Vector{OnceDifferentiable}()  # datatype OnceDifferentiable not found here
    for i = 1: num_partition
            push!(df_array, prepare_for_nonlin_solve!(ssp_array[i]))
    end
    println("Initialized nonlinear solve...")




    println("Starting  iterations...")

    for i = 1: 200

            old = interface_array(ssp_array, partition)

            for sn_id in partition["level"][1]
                    solver1 = solve_on_network!(ssp_array[sn_id], df_array[sn_id], show_trace_flag=show_trace_flag)
            end
            update_interface_potentials_of_nbrs!(ssp_array, partition, 1)

            for sn_id in partition["level"][2]
                    solver2 = solve_on_network!(ssp_array[sn_id], df_array[sn_id], show_trace_flag=show_trace_flag)
            end
            update_interface_potentials_of_nbrs!(ssp_array, partition, 2)

            for sn_id in partition["level"][3]
                    solver3 = solve_on_network!(ssp_array[sn_id], df_array[sn_id], show_trace_flag=show_trace_flag)
            end
            update_interface_potentials_of_nbrs!(ssp_array, partition, 3)

            for sn_id in partition["level"][4]
                    solver4 = solve_on_network!(ssp_array[sn_id], df_array[sn_id], show_trace_flag=show_trace_flag)
            end


            update_interface_transfers!(ssp_array,  partition)

            new = interface_array(ssp_array, partition)
            

            err = norm(new - old)

            if err < 1e-4
                    println(i, " iterations, converged ")
                    break
            end

            
    end

    x_dof = combine_subnetwork_solutions(ss, ssp_array)

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

    