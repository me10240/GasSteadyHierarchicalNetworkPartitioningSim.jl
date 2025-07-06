


# splitting can be radial --  meaning sn_1 (with slack)  intersects sn_j for all j > 1 , but 
# for any distinct j, k sn_j intersects sn_k is null whenever, neither j, k =1.




using GasSteadyHierarchicalNetworkPartitioningSim
using JSON
using LinearAlgebra
using NLSolversBase


file = "./data/GasLib-40-split/"


eos_var = :ideal
ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="") 
df = prepare_for_nonlin_solve!(ss)

filepath = "./data/GasLib-40-split/partition-7/partition7.json"

x_dof = run_partitioned_ss(filepath, ss)



# partition = create_partition(filepath)

# if length(partition["level"][4]) != 0
#         error("This partitioning is not acceptable")
# end
# num_partition = length(partition["subnetworks"])



# ssp_array = Vector{SteadySimulator}()
# for i = 1 : num_partition
#         push!(ssp_array, initialize_simulator_subnetwork(ss, partition["subnetworks"][i]["node_list"]))
# end
# println("Initialized steady state simulator...")


# create_interface_transfers!(ssp_array,  partition, val=0)
# designate_interface_nodes_as_slack!(ssp_array, partition)


# df_array = Vector{OnceDifferentiable}()  
# for i = 1: num_partition
#         push!(df_array, prepare_for_nonlin_solve!(ssp_array[i]))
# end
# println("Initialized nonlinear solve...")


# # println(ssp_array[2].ref[:pipe])
# # global solver1, solver2
# # println("Starting  iterations...")

# for i = 1: 20

#         println(i)
#         old = interface_array(ssp_array, partition)
#         # println(old)

#         for sn_id in partition["level"][1]
#                 solver1 = solve_on_network!(ssp_array[sn_id], df_array[sn_id])
#                 # println(solver1.solution)
#                 for (i, _) in ssp_array[1].ref[:node]
#                         # println(ssp_array[1].ref[:node][i])
#                 end
#                 println(ssp_array[1].ref[:node][21]["potential"], "\t",  ssp_array[1].ref[:node][27]["potential"])
#         end
#         println(" a")
#         update_interface_potentials_of_nbrs!(ssp_array, partition, 1)
#         # println(ssp_array[2].ref[:node][21]["potential"], ssp_array[2].ref[:node][27]["potential"])

#         for sn_id in partition["level"][2]
#                 solver2 = solve_on_network!(ssp_array[sn_id], df_array[sn_id])
#                 # println(solver2.solution)
#                 # for (i, _) in ssp_array[2].ref[:pipe]
#                 #         println(ssp_array[2].ref[:pipe][i])
#                 # end
#         end
#         println("b")

#         # update_interface_potentials_of_nbrs!(ssp_array, partition, 2)

#         # for sn_id in partition["level"][3]
#         #         solver3 = solve_on_network!(ssp_array[sn_id], df_array[sn_id])
#         # end
#         # update_interface_potentials_of_nbrs!(ssp_array, partition, 3)

#         # for sn_id in partition["level"][4]
#         #         solver4 = solve_on_network!(ssp_array[sn_id], df_array[sn_id])
#         # end


#         update_interface_transfers!(ssp_array,  partition)
#         println("c")


#         new = interface_array(ssp_array, partition)
#         println(new)

#         err = norm(new - old)
#         println(err)

#         if err < 1e-4
#                 println(i, " iterations, converged ")
#                 break
#         end

        
# end

# x_dof = combine_subnetwork_solutions(ss, ssp_array)


var = value!(df, x_dof)
println(norm(var))

solver = solve_on_network!(ss, df, x_guess=x_dof)
# solver = solve_on_network!(ss, df, x_guess=x_dof, iteration_limit=1)

println(solver.iterations, " ", solver.residual_norm)
# println(norm(x_dof - solver.solution))


# populate_solution!(ss)




