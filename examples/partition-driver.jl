

using GasSteadySim
using JSON
using LinearAlgebra
using NLSolversBase
using PyCall
println(pwd())
pushfirst!(pyimport("sys")."path", "")


file = "./data/8-node/"
# file = "./data/Texas7k_Gas/"
partition_module = pyimport("network_partition_script")
partition_module.run_script(file, allow_slack_node_partitioning = false, num_max=2, round_max=10, plotting_flag=true) #bool True in python is true in julia


eos_var = :ideal
t11 = @elapsed ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="") 
t12 = @elapsed df = prepare_for_nonlin_solve!(ss)
t13 = @elapsed solver = solve_on_network!(ss, df, iteration_limit=30,  show_trace_flag=true)
# solver = solve_on_network!(ss, df, x_guess=x_dof, iteration_limit=1)
println(solver.iterations, " ", solver.residual_norm)



filepath = file * "partition-test-script.json"

t21 = @elapsed x_dof, cond_number_array = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false)
push!(cond_number_array, cond(gradient(df), 1))
println("Condition numbers: \n", cond_number_array)
t22 = @elapsed var = value!(df, x_dof)
println(norm(var))
solver = solve_on_network!(ss, df, x_guess=x_dof)
# solver = solve_on_network!(ss, df, x_guess=x_dof, iteration_limit=1)
println(solver.iterations, " ", solver.residual_norm)



time_full = t11 + t12 + t13
time_partition = t21 + t22

println(time_full, "\t", time_partition)
# println(norm(x_dof - solver.solution))


# populate_solution!(ss)


# time comparison - total time for single network vs total time for




