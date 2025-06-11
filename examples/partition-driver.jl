

using GasSteadySim
using JSON
using LinearAlgebra
using NLSolversBase
using PyCall

# file = "./data/8-node/"
file = "./data/GasLib-40/"
# file = "./data/Texas7k_Gas/"
run_pycall = false

if run_pycall == true
    println(pwd())
    pushfirst!(pyimport("sys")."path", "")
    partition_module = pyimport("network_partition_script")
    partition_module.run_script(file, loglevel="debug", allow_slack_node_partitioning = false, num_max=10, round_max=2, plotting_flag=true) #bool True in python is true in julia
end

# why is partition solve failing at second level ? check what's happening
eos_var = :ideal
t11 = @elapsed ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="") 
t12 = @elapsed df = prepare_for_nonlin_solve!(ss)
t13 = @elapsed solver = solve_on_network!(ss, df, show_trace_flag=true, iteration_limit=2000, method=:trust_region)
# solver = solve_on_network!(ss, df, x_guess=x_dof, iteration_limit=1)
println(solver.iterations, " ", solver.residual_norm)



filepath = file * "partition-test-script-new.json"

t21 = @elapsed x_dof = run_partitioned_ss(filepath, ss, eos=eos_var, show_trace_flag=true, iteration_limit=100, method=:trust_region);
# t22 = @elapsed var = value!(df, x_dof)
# println(norm(var))
# solver = solve_on_network!(ss, df, x_guess=x_dof, method=:trust_region)
# solver = solve_on_network!(ss, df, x_guess=x_dof, iteration_limit=1)
# println(solver.iterations, " ", solver.residual_norm)



# time_full = t11 + t12 + t13
# time_partition = t21 + t22

# println(time_full, "\t", time_partition)
# println(norm(x_dof - solver.solution))


# populate_solution!(ss)


# time comparison - total time for single network vs total time for




