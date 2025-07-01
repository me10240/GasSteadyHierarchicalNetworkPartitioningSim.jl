

using GasSteadySim
using JSON
using LinearAlgebra
using NLSolversBase

# file = "./data/8-node/"
# file = "./data/GasLib-40/"
file = "./data/GasLib-40-multiple-slacks-new/"
# file = "../test/data/GasLib-40-multiple-slacks/"


# file = "./data/Texas7k_Gas/"


eos_var = :ideal
t11 = @elapsed ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="") 
# t12 = @elapsed df = prepare_for_nonlin_solve!(ss)
# t13 = @elapsed solver = solve_on_network!(ss, df, show_trace_flag=true, iteration_limit=2000, method=:newton)

partition_data_or_file = create_partition(ss; num_partitions=2, write_to_file=true, filepath=file *"partition_data.json")
# partition_data_or_file = file * "partition_data.json"
# partition_data_or_file = file * "partition-test-script-dummy.json"
@show partition_data_or_file["interface_nodes"]
# x_guess = Vector{Float64}()
# for i in partition_data_or_file["interface_nodes"]
#         val = ss.ref[:node][i]["potential"]
#         push!(x_guess, val)
# end
# x_guess = [0.9985186255278212, 2.2509512091126638]
# x_guess = [0.5501218168950411, 1.2339331568082517]
x_guess = [0.9985, 2.2509]


t21 = @elapsed x_dof = run_partitioned_ss(partition_data_or_file, ss, eos=eos_var, show_trace_flag=true, iteration_limit=100, method=:trust_region, x_guess=x_guess);

if isnothing(x_dof) == false
    df = prepare_for_nonlin_solve!(ss)
    @info "Checking obtained solution by plugging into full system..."
    value!(df, x_dof)
    @info("Residual (inf norm) is : $(norm(value(df), Inf))")
    @info "Running NR with obtained solution as initial guess..."
    solver = solve_on_network!(ss, df, x_guess=x_dof, method=:trust_region)
    @info("x-x* (inf norm): $(norm(x_dof - solver.solution, Inf))")
    @info("Solver iterations and residual norm: $(solver.iterations), $(solver.residual_norm)")
end


# Q2 - condition number of Jacobian in sensitivity linear system











