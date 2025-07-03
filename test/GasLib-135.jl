
@testset "test GasLib-135 ideal run" begin
    file = "./data/GasLib-135/"
    eos_var = :ideal
    ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
    filepath = file * "partition-test-script.json"
    _, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    exact_sol = GasSteadySim._parse_json("data/GasLib-135/exact_sol_ideal.json")
    _check_correctness(ss.sol, exact_sol)


end

@testset "test GasLib-135 simple CNGA run" begin
    file = "./data/GasLib-135/"
    eos_var = :simple_cnga
    ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
    filepath = file * "partition-test-script.json"

    _, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    exact_sol = GasSteadySim._parse_json("data/GasLib-135/exact_sol_simple_cnga.json")
    _check_correctness(ss.sol, exact_sol)
end


@testset "test GasLib-135 full CNGA run" begin
    file = "./data/GasLib-135/"
    eos_var = :full_cnga
    ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
    filepath = file * "partition-test-script.json"

    _, _ = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
    exact_sol = GasSteadySim._parse_json("data/GasLib-135/exact_sol_full_cnga.json")
    _check_correctness(ss.sol, exact_sol)
end
