

using GasSteadyHierarchicalNetworkPartitioningSim
using Graphs
using JSON
using LinearAlgebra
using NLSolversBase
using PyCall

file = "./data/8-node/"
# file = "./data/GasLib-40-multiple-slacks-2/"
# file = "./data/GasLib-135/"
# file = "./data/GasLib-24/"
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

eos_var = :simple_cnga
ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="") 
# df = prepare_for_nonlin_solve!(ss)
# solver = solve_on_network!(ss, df, show_trace_flag=true, iteration_limit=2000, method=:trust_region)
#============== Save solution data for use =============================#
# populate_solution!(ss)
# filename = file * "exact_sol_$eos_var.json"
# open(filename, "w") do f 
#         JSON.print(f, ss.sol, 2)
# end
#=======================================================================#


partition_data_or_file = create_partition_with_cut_points(ss; max_nodes= 4, allow_slack_node_partitioning=true, write_to_file = true, filepath = file * "partition-expt.json")
# partition_data_or_file = file * "partition-test-script.json"

out = run_partitioned_ss(partition_data_or_file, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
if !isnothing(out)
    x_dof, cond_number_array = out
    df = prepare_for_nonlin_solve!(ss)
    @info "Checking obtained solution by plugging into full system..."
    value!(df, x_dof)
    @info("Residual (inf norm) is : $(norm(value(df), Inf))")
    @info "Running NR with obtained solution as initial guess..."
    solver = solve_on_network!(ss, df, x_guess=x_dof, method=:trust_region)
    @info("x-x* (inf norm): $(norm(x_dof - solver.solution, Inf))")
    @info("Solver iterations and residual norm: $(solver.iterations), $(solver.residual_norm)")
end


    



