


# splitting can be radial --  meaning sn_1 (with slack)  intersects sn_j for all j > 1 , but 
# for any distinct j, k sn_j intersects sn_k is null whenever, neither j, k =1.




using GasSteadySim
using JSON
using LinearAlgebra
using NLSolversBase


file = "./data/GasLib-40-multiple-slacks-2/"


eos_var = :ideal
ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="") 
df = prepare_for_nonlin_solve!(ss)

filepath = file * "partition-test-script.json"

x_dof = run_partitioned_ss(filepath, ss)


var = value!(df, x_dof)
println(norm(var))

solver = solve_on_network!(ss, df, x_guess=x_dof)
# solver = solve_on_network!(ss, df, x_guess=x_dof, iteration_limit=1)

println(solver.iterations, " ", solver.residual_norm)
# println(norm(x_dof - solver.solution))


# populate_solution!(ss)




