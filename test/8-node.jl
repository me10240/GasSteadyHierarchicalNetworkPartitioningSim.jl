

@testset "test 8 node ideal run" begin
    file = "./data/8-node/"
    eos_var = :ideal
    ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
    filepath = file * "partition-test-script.json"
    x_dof1, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    @test ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow) ≈ -275.00 atol = 1e-2
    exact_sol = GasSteadyHierarchicalNetworkPartitioningSim._parse_json("data/8-node/exact_sol_ideal.json")
    _check_correctness(ss.sol, exact_sol)
    partition = create_partition_with_cut_points(ss; max_nodes= 4, allow_slack_node_partitioning=true, write_to_file = false, filepath = file * "partition-expt.json")
    x_dof2, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    @test isapprox(norm(x_dof1 - x_dof2, Inf), 0, atol = 1e-7)


end

@testset "test 8 node simple CNGA run" begin
    file = "./data/8-node/"
    eos_var = :simple_cnga
    ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
    filepath = file * "partition-test-script.json"
    x_dof1, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    @test ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow) ≈ -275.00 atol = 1e-2
    exact_sol = GasSteadyHierarchicalNetworkPartitioningSim._parse_json("data/8-node/exact_sol_simple_cnga.json")
    _check_correctness(ss.sol, exact_sol)
    partition = create_partition_with_cut_points(ss; max_nodes= 4, allow_slack_node_partitioning=true, write_to_file = false, filepath = file * "partition-expt.json")
    x_dof2, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    @test isapprox(norm(x_dof1 - x_dof2, Inf), 0, atol = 1e-7)
end


@testset "test 8 node full CNGA run" begin
    file = "./data/8-node/"
    eos_var = :full_cnga
    ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
    filepath = file * "partition-test-script.json"
    x_dof1, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    @test ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow) ≈ -275.00 atol = 1e-2
    exact_sol = GasSteadyHierarchicalNetworkPartitioningSim._parse_json("data/8-node/exact_sol_full_cnga.json")
    _check_correctness(ss.sol, exact_sol)
    partition = create_partition_with_cut_points(ss; max_nodes= 4, allow_slack_node_partitioning=true, write_to_file = false, filepath = file * "partition-expt.json")
    x_dof2, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    @test isapprox(norm(x_dof1 - x_dof2, Inf), 0, atol = 1e-7)
end
