

@testset "test 8 node ideal run" begin
    file = "./data/8-node/"
    eos_var = :ideal
    ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
    filepath = file * "partition-test-script.json"
    _, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    @test ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow) ≈ -275.00 atol = 1e-2
    exact_sol = GasSteadyHierarchicalNetworkPartitioningSim._parse_json("data/8-node/exact_sol_ideal.json")
    _check_correctness(ss.sol, exact_sol)


end

@testset "test 8 node simple CNGA run" begin
    file = "./data/8-node/"
    eos_var = :simple_cnga
    ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
    filepath = file * "partition-test-script.json"
    _, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    @test ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow) ≈ -275.00 atol = 1e-2
    exact_sol = GasSteadyHierarchicalNetworkPartitioningSim._parse_json("data/8-node/exact_sol_simple_cnga.json")
    _check_correctness(ss.sol, exact_sol)
end


@testset "test 8 node full CNGA run" begin
    file = "./data/8-node/"
    eos_var = :full_cnga
    ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
    filepath = file * "partition-test-script.json"
    _, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    @test ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow) ≈ -275.00 atol = 1e-2
    exact_sol = GasSteadyHierarchicalNetworkPartitioningSim._parse_json("data/8-node/exact_sol_full_cnga.json")
    _check_correctness(ss.sol, exact_sol)
end
