

using GasSteadyHierarchicalNetworkPartitioningSim
using JSON
using LinearAlgebra
using NLSolversBase
using PyCall
import GasSteadySim as GSS

# file = "./data/8-node/"
file = "./data/GasLib-24/"
# file = "./data/GasLib-135/"

# file = "./data/GasLib-40/"
# file = "./data/Texas7k_Gas/"
# file = "../test/data/GasLib-40/"
run_pycall = false

if run_pycall == true
    println(pwd())
    pushfirst!(pyimport("sys")."path", "")
    partition_module = pyimport("network_partition_script")
    partition_module.run_script(file, loglevel="info", allow_slack_node_partitioning = false, num_max=2, round_max=5, plotting_flag=true) #bool True in python is true in julia
end

# why is partition solve failing at second level ? check what's happening
eos_var = :ideal
ss = GSS.initialize_simulator(file, eos=eos_var, initial_guess_filename="") 
ss_copy = deepcopy(ss)
# t12 = @elapsed df = prepare_for_nonlin_solve!(ss)
# t13 = @elapsed solver = solve_on_network!(ss, df, show_trace_flag=true, iteration_limit=2000, method=:trust_region)
#============== Save solution data for use =============================#
# populate_solution!(ss)
# filename = file * "exact_sol_$eos_var.json"
# open(filename, "w") do f 
#         JSON.print(f, ss.sol, 2)
# end
#=======================================================================#
# solver = solve_on_network!(ss, df, x_guess=x_dof, iteration_limit=1)
# println(solver.iterations, " ", solver.residual_norm)



filepath = file * "partition-test-script.json"

x_dof, cond_number_array = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=true, iteration_limit=2000, method=:trust_region)
df = prepare_for_nonlin_solve!(ss_copy)
push!(cond_number_array, cond(gradient(df), 1))
println("Condition numbers: \n", cond_number_array)
var = value!(df, x_dof)
println(norm(var))
solver = solve_on_network!(ss_copy, df, x_guess=x_dof, method=:trust_region)
# solver = solve_on_network!(ss, df, x_guess=x_dof, iteration_limit=1)
println(solver.iterations, " ", solver.residual_norm)








