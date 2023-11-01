




using GasSteadySim
using JSON
using LinearAlgebra
using NLSolversBase


# file = "./data/GasLib-135/"
# file = "./data/GasLib-134/"
# file = "./data/GasLib-40-split/"
# file = "./data/GasLib-24/"
file = "./data/GasLib-11/"

# file = "./data/8-node/"


function create_cnga_solution_pressure_formulation(ss_p::SteadySimulator, ss_pot::SteadySimulator)::Vector{Float64}
    ndofs = length(ref(ss_p, :dof))
    x_dofs = zeros(Float64, ndofs) 
    dofs_updated = 0

    for i = 1:length(ref(ss_pot, :dof))
        comp, id = ss_pot.ref[:dof][i]
        if comp == :node
            val = ref(ss_p, :is_pressure_node, id) ? ss_pot.ref[comp][id]["pressure"] : ss_pot.ref[comp][id]["potential"]
            x_dofs[ss_p.ref[comp][id]["dof"]] = val
        elseif comp == :pipe
            x_dofs[ss_p.ref[comp][id]["dof"]] = ss_pot.ref[comp][id]["flow"]
        elseif comp == :compressor
            x_dofs[ss_p.ref[comp][id]["dof"]] = ss_pot.ref[comp][id]["flow"]
        end
    end

    return x_dofs
end


eos_var = :simple_cnga
ss_pot = initialize_simulator(file, eos=eos_var, potential_formulation_flag=true, potential_ratio_approx=[0.0, 0.0, 0.9, 0.1], initial_guess_filename="") 
df_pot = prepare_for_nonlin_solve!(ss_pot)
solver_pi = solve_on_network!(ss_pot, df_pot)

ss_p = initialize_simulator(file, eos=eos_var, initial_guess_filename="") 
df_p = prepare_for_nonlin_solve!(ss_p)

x_dof = create_cnga_solution_pressure_formulation(ss_p, ss_pot)

num_compressors = length(ss_p.ref[:compressor])
var = value!(df_p, x_dof)

println("Pipe + Node Residual (inf norm): ", norm(var[1:end-num_compressors], Inf), "\nCompressor  Residual (inf norm): ", norm(var[end-num_compressors+1:end], Inf))

# println(findall(x->abs(x)>1e-3, var)," ", var, " ", norm(var), " ", norm(var, Inf))

solver_p = solve_on_network!(ss_p, df_p, x_guess= x_dof, show_trace_flag=true)
# solver_p = solve_on_network!(ss, df, x_guess=x_dof, iteration_limit=1)

println(solver_p.iterations, " ", solver_p.residual_norm)
println("Absolute error of approx soln (inf norm): ", norm(x_dof - solver_p.solution, Inf), "\nRelative error of approximate soln (inf norm): ", norm(x_dof - solver_p.solution, Inf)/norm(solver_p.solution, Inf))



# for i = 1:length(ss_p.ref[:node])
#     v1 = ss_p.ref[:node][i]["potential"]
#     v2 = ss_pot.ref[:node][i]["potential"]
#     println(v1, " ", v2, " ", abs(v1-v2))
# end
# println("nodes finished")
# for i = 1:length(ss_p.ref[:pipe])
#     v1 = ss_p.ref[:pipe][i]["flow"]
#     v2 = ss_pot.ref[:pipe][i]["flow"]
#     println(v1, " ", v2, " ", abs(v1-v2))
# end
# println("pipes finished")

# for i = 1:length(ss_p.ref[:compressor])
#     v1 = ss_p.ref[:compressor][i]["flow"]
#     v2 = ss_pot.ref[:compressor][i]["flow"]
#     println(v1, " ", v2, " ", abs(v1-v2))
# end
# println("compressors finished")


