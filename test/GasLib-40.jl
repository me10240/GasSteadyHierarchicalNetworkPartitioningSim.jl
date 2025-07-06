
@testset "test GasLib-40 ideal run" begin
    file = "./data/GasLib-40/"
    eos_var = :ideal
    ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
    filepath = file * "partition-test-script.json"
    _, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    exact_sol = GasSteadyHierarchicalNetworkPartitioningSim._parse_json("data/GasLib-40/exact_sol_ideal.json")
    _check_correctness(ss.sol, exact_sol)


end

@testset "test GasLib-40 simple CNGA run" begin
    file = "./data/GasLib-40/"
    eos_var = :simple_cnga
    ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
    filepath = file * "partition-test-script.json"

    _, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    exact_sol = GasSteadyHierarchicalNetworkPartitioningSim._parse_json("data/GasLib-40/exact_sol_simple_cnga.json")
    _check_correctness(ss.sol, exact_sol)
end


@testset "test GasLib-40 full CNGA run" begin
    file = "./data/GasLib-40/"
    eos_var = :full_cnga
    ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
    filepath = file * "partition-test-script.json"

    _, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    exact_sol = GasSteadyHierarchicalNetworkPartitioningSim._parse_json("data/GasLib-40/exact_sol_full_cnga.json")
    _check_correctness(ss.sol, exact_sol)
end
