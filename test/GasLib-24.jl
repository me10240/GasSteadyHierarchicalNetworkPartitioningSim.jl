
@testset "test GasLib-24 ideal run" begin
    file = "./data/GasLib-24/"
    eos_var = :ideal
    ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
    filepath = file * "partition-test-script.json"
    x_dof1, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    exact_sol = GasSteadyHierarchicalNetworkPartitioningSim._parse_json("data/GasLib-24/exact_sol_ideal.json")
    _check_correctness(ss.sol, exact_sol)
    partition = create_partition_with_cut_points(ss; max_nodes= 15, allow_slack_node_partitioning=true, write_to_file = false, filepath = file * "partition-expt.json")
    x_dof2, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    @test isapprox(norm(x_dof1 - x_dof2, Inf), 0, atol = 5e-6)


end

@testset "test GasLib-24 simple CNGA run" begin
    file = "./data/GasLib-24/"
    eos_var = :simple_cnga
    ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
    filepath = file * "partition-test-script.json"

    x_dof1, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    exact_sol = GasSteadyHierarchicalNetworkPartitioningSim._parse_json("data/GasLib-24/exact_sol_simple_cnga.json")
    _check_correctness(ss.sol, exact_sol)
    partition = create_partition_with_cut_points(ss; max_nodes= 15, allow_slack_node_partitioning=true, write_to_file = false, filepath = file * "partition-expt.json")
    x_dof2, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    @test isapprox(norm(x_dof1 - x_dof2, Inf), 0, atol = 5e-6)
end


@testset "test GasLib-24 full CNGA run" begin
    file = "./data/GasLib-24/"
    eos_var = :full_cnga
    ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
    filepath = file * "partition-test-script.json"

    x_dof1, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    exact_sol = GasSteadyHierarchicalNetworkPartitioningSim._parse_json("data/GasLib-24/exact_sol_full_cnga.json")
    _check_correctness(ss.sol, exact_sol)
    partition = create_partition_with_cut_points(ss; max_nodes= 15, allow_slack_node_partitioning=true, write_to_file = false, filepath = file * "partition-expt.json")
    x_dof2, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    @test isapprox(norm(x_dof1 - x_dof2, Inf), 0, atol = 5e-6)
end
