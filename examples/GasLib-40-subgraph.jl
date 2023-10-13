using GasSteadySim

file = "./data/GasLib-40-subgraph/"
eos_var = :ideal
ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")
solver_return = run_simulator!(ss)

println(solver_return.status)
