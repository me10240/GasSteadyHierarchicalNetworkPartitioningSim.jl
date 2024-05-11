

using GasSteadySim
using JSON
using LinearAlgebra
using NLSolversBase


file = "./data/GasLib-582/"


eos_var = :ideal
ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="") 
df = prepare_for_nonlin_solve!(ss)

filepath = file * "partition-test-script.json"

x_dof = run_partitioned_ss(filepath, ss, eos=eos_var)


var = value!(df, x_dof)
println(norm(var))

solver = solve_on_network!(ss, df, x_guess=x_dof)
# solver = solve_on_network!(ss, df, x_guess=x_dof, iteration_limit=1)

println(solver.iterations, " ", solver.residual_norm)
# println(norm(x_dof - solver.solution))


# populate_solution!(ss)




