

using GasSteadySim
using JSON
using LinearAlgebra
using NLSolversBase
using PyCall

# file = "./data/8-node/"
# file = "./data/GasLib-40/"
file = "./data/Texas7k_Gas/"
run_pycall = false

if run_pycall == true
    println(pwd())
    pushfirst!(pyimport("sys")."path", "")
    partition_module = pyimport("network_partition_script")
    partition_module.run_script(file, loglevel="debug", allow_slack_node_partitioning = false, num_max=10, round_max=2, plotting_flag=true) #bool True in python is true in julia
end

eos_var = :ideal
t11 = @elapsed ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="") 
t12 = @elapsed df = prepare_for_nonlin_solve!(ss)
t13 = @elapsed solver = solve_on_network!(ss, df, show_trace_flag=true, iteration_limit=2000, method=:newton)
# solver = solve_on_network!(ss, df, x_guess=x_dof, iteration_limit=1)
println(solver.iterations, " ", solver.residual_norm)

partition_data_or_file = create_partition(ss; num_partitions=5, write_to_file=true, filepath=file *"partition_data.json")
# partition_data_or_file = file * "partition_data.json"
# partition_data_or_file = file * "partition-test-script-dummy.json"




t21 = @elapsed x_dof = run_partitioned_ss(partition_data_or_file, ss, eos=eos_var, show_trace_flag=true, iteration_limit=100, method=:trust_region, x_guess = []);

if isnothing(x_dof) == false
    value!(df, x_dof)
    # println(value(df))
    println(norm(value(df), Inf))
    println("x-x*:", norm(x_dof - solver.solution))

    solver = solve_on_network!(ss, df, x_guess=x_dof, method=:trust_region)
    # solver = solve_on_network!(ss, df, x_guess=x_dof, iteration_limit=1)
    println(solver.iterations, " ", solver.residual_norm)
end

# Q1 - if all subnetwork residuals are small, then is over all residual with xdof large at interface nodes only ? check 
# Q2 - condition number of Jacobian in sensitivity linear system
# Q3 - is lufact doing sparse solve ? if not ensure sparse solve.. maybe by krylov methods ?











