

@testset "test GasLib-40-multiple-slacks ideal run" begin
    file = "./data/GasLib-40-multiple-slacks/"
    exact_sol = GasSteadySim._parse_json("data/GasLib-40-multiple-slacks/exact_sol_ideal.json")

    eos_var = :ideal
    ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
    partition_data_or_file = create_partition(ss; num_partitions=2, write_to_file=false, filepath="")
    # @info partition_data_or_file["interface_nodes"]
    x_guess = Vector{Float64}()
    for i in partition_data_or_file["interface_nodes"]
        pval = exact_sol["nodal_pressure"][string(i)] / ss.nominal_values[:pressure]
        push!(x_guess, get_potential(ss, pval))
    end
    @show x_guess
    # x_guess = [0.9985186255278212, 2.2509512091126638]
    # x_guess = [0.5501218168950411, 1.2339331568082517]
    # x_guess = [0.55, 1.234]
    

    x_dof = run_partitioned_ss(partition_data_or_file, ss, eos=eos_var, show_trace_flag=true, iteration_limit=100, method=:trust_region, x_guess=x_guess);
    # @show solver_return.iterations
    # _check_correctness(ss.sol, exact_sol)


end

# @testset "test GasLib-40-multiple-slacks simple CNGA run" begin
#     file = "./data/GasLib-40-multiple-slacks/"
#     eos_var = :simple_cnga
#     ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
#     # solver_return = run_simulator!(ss)
#     partition_data_or_file = create_partition(ss; num_partitions=2, write_to_file=false, filepath="")
#     x_dof = run_partitioned_ss(partition_data_or_file, ss, eos=eos_var, show_trace_flag=true, iteration_limit=100, method=:newton);
#     exact_sol = GasSteadySim._parse_json("data/GasLib-40-multiple-slacks/exact_sol_simple_cnga.json")
#     _check_correctness(ss.sol, exact_sol)
# end


# @testset "test GasLib-40-multiple-slacks full CNGA run" begin
#     file = "./data/GasLib-40-multiple-slacks/"
#     eos_var = :full_cnga
#     ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
#     # solver_return = run_simulator!(ss)
#     partition_data_or_file = create_partition(ss; num_partitions=2, write_to_file=false, filepath="")
#     x_dof = run_partitioned_ss(partition_data_or_file, ss, eos=eos_var, show_trace_flag=true, iteration_limit=100, method=:newton);
#     exact_sol = GasSteadySim._parse_json("data/GasLib-40-multiple-slacks/exact_sol_full_cnga.json")
#     _check_correctness(ss.sol, exact_sol)
# end
