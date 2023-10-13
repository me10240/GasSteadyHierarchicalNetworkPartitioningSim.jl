#issues 

#does it matter if slack node is split ? complications since injection at slack not known ?

# need to supply nominal values in input file.
# make sure they are same set of nominal values
# right now nominal pressure is set from slack alone

# instead of accessing bc from boundary_condition each time, what if we update nodes with withdrawal (given val or 0) or pressure (if slack)
# and compressors with c_ratio so that there is no need for get_nodal_control etc
# this way, can directly overwrite bc values and have seamless simulation  
# only argument against it is that we actually solve for slack pressures through an equation
# easy for nodes, compressors...how to do this for loss resistors, valves etc ?


# suppose we partition network at a vertex which has non-zero injection/withdrawal (q) --- need to absorb the q into first network, 
# since vertex will become a slack node for second network. 
# need to store the given q without modification since we need it for update in iteration
# how to do this cleanly without causing confusion -- store q in withdrawal, add key transfer that is updated


# can we create full ref and then supply partition network  as set of vertex ids ? Then simply delete each vertex and edge that 
# involves vertices not in set, and also prune bc data struc. thus starting from two copies, ss1.ref and ss2.ref can be obrtained. 
# this will have advantage that non-dim etc have already been done correctly on data before partitioning.
# is this a good way to do this instead of using Plasmo extensively ?
# problems -- will two copies of network ref take up too much memory ?

# idea of vertex splitting is to imagine the concerned vertices  divides into two and  connected by transfers
# topology of vertex split, how complicated can connections get ?


# having set of partitioning vertices will allow easy iteration update and checking for convergence. How ? By taking intersection,
# can find vertices belonging to both sets.
# How will this work if we have more than two subnetworks -- gasLib40 can be split into 3 networks with single common vertex, say X.
# N1 with slack and assumed q at split vertex X gives slack pressure for X in N2, N3. slack injections for X in N2, N3 is used to update 
# q for X in N1 and start another round of iteration. 
# so at each common node, find number of networks it belongs to, say n. need n-1 transfers vto update ?

# after convergence, combine solution from both networks to get solution for whole network
# if complete ref exists, can incorporate subnetwork solutions

# theory: can we prove ideal gas system has generalized solution (ie, potential, flow)  by starting with existence  when edge reln is linear.


# theory ; this is like neumanntodirichlet rather than dirichlet-to=neumann
# cannot do dirichletneumann because  then first network multi-slack, second network no slack
# understand convergence of dirichlet-neumann...can we use it to prove neumannndirichlet and specifically for this network
# does dirichletneumann work for more than 2 networks





using GasSteadySim
using JSON

file = "./data/GasLib-40-split/"
eos_var = :ideal
ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="")


# ss1 = deepcopy(ss)
# node_list = [9, 12, 18, 22, 23, 25, 33, 34, 35]

# for key_level_1 in keys(ss1.ref)
#         for key2 in keys(ss1.ref[key_level_1])
#                 if typeof(key2) == Int in node_list

#         end
# end

# ss2 = deepcopy(ss)




# slack is node 38
# N2 = [9, 12, 18, 22, 23, 25, 33, 34, 35]
# N1 = [1-40]\N2 U [33, 18]
# Only 33, 18 common to both N1, N2






file1 = "./data/GasLib-40-split/N1/"
file2 = "./data/GasLib-40-split/N2/"

eos_var = :ideal
ss1 = initialize_simulator(file1, eos=eos_var, initial_guess_filename="") 
println("n1 done")


# add transfer fields to network 1
# since there is only one network connected to N1, no ambiguity here
# but consider if N1, N2, N3 share vertex 4, or suppose N1, N2 share 4, and N1 and N3 share 7 etc..in such cases probably need 
# ss1.ref["transfer"][N2] will be a dictionary with nodeid and vals


ss1.ref[:node][18]["transfer"] = Dict{String, Any}("n1"=>0) 
ss1.ref[:node][33]["transfer"] = Dict{String, Any}("n1"=>0) 



ss2 = initialize_simulator(file2, eos=eos_var, initial_guess_filename="") 
println("n2 done")

df1 = prepare_for_nonlin_solve!(ss1)
df2 = prepare_for_nonlin_solve!(ss2)
df = prepare_for_nonlin_solve!(ss)






function get_interface_dofs(ss1)
        return [ss1.ref[:node][18]["transfer"]["n1"], ss1.ref[:node][33]["transfer"]["n1"]]
end

function norm(a)
        return sqrt(sum( a .* a))
end

function _combine_solutions!(ss::SteadySimulator, ss1::SteadySimulator, ss2::SteadySimulator)
    ndofs = length(ref(ss, :dof))
    x_dofs = zeros(Float64, ndofs) 
    dofs_updated = 0


    for i = 1:length(ref(ss1, :dof))
        comp, id = ss1.ref[:dof][i]
        if comp == :node
                x_dofs[ss.ref[comp][id]["dof"]] = ss1.ref[comp][id]["potential"]
        elseif comp == :pipe
                x_dofs[ss.ref[comp][id]["dof"]] = ss1.ref[comp][id]["flow"]
        elseif comp == :compressor
                x_dofs[ss.ref[comp][id]["dof"]] = ss1.ref[comp][id]["flow"]
        end
    end

    for i = 1:length(ref(ss2, :dof))
        comp, id = ss2.ref[:dof][i]
        if comp == :node
                x_dofs[ss.ref[comp][id]["dof"]] = ss2.ref[comp][id]["potential"]
        elseif comp == :pipe
                x_dofs[ss.ref[comp][id]["dof"]] = ss2.ref[comp][id]["flow"]
        elseif comp == :compressor
                x_dofs[ss.ref[comp][id]["dof"]] = ss2.ref[comp][id]["flow"]
        end
    end

    return x_dofs
end



for i = 1: 20

        # take values from ss2.ref and update ss1.ref
        
        solver1 = solve_on_network!(ss1, df1)
        old = get_interface_dofs(ss1)

        println(solver1.status)
        println(ss1.ref[:node][5]["potential"], "\t", ss1.ref[:node][16]["potential"]) 
        # println(ss1.ref[:node][4]["withdrawal"], "\t", ss1.ref[:node][7]["withdrawal"])  #debug this
        # println(ss1.ref[:node][4]["transfer"], "\t", ss1.ref[:node][7]["transfer"])  #debug this

        println(i, "a")


        # BC for N2 comes from N1 solution
        ss2.ref[:node][18]["potential"] = ss1.ref[:node][18]["potential"]  
        ss2.ref[:node][33]["potential"] = ss1.ref[:node][33]["potential"]
        println(ss2.ref[:node][18]["potential"], "\t", ss2.ref[:node][18]["potential"]) 


        solver2 = solve_on_network!(ss2, df2)
        println(i, "b")

        println(solver2.status)
        # println(ss2.ref[:node][7]["potential"])

        #Update BC for N1
        ss1.ref[:node][18]["transfer"]["n1"] = 0.5 * ss1.ref[:node][18]["transfer"]["n1"] -0.5 * ss2.ref[:node][18]["withdrawal"]  
        ss1.ref[:node][33]["transfer"]["n1"] = 0.5 * ss1.ref[:node][33]["transfer"]["n1"] -0.5 * ss2.ref[:node][33]["withdrawal"]
        new = get_interface_dofs(ss1)

        err = norm(new - old)

        if err < 1e-3
                println("converged")
                println(solver1.solution, solver2.solution)
                break
        end
        # @show ss1.ref[:node][5]["potential"], ss1.ref[:node][16]["potential"], ss1.ref[:node][20]["potential"]
        # println(ss1.ref[:node][4]["withdrawal"], "\t", ss1.ref[:node][7]["withdrawal"])  #debug this

        
end

x_dof = _combine_solutions!(ss, ss1, ss2)
solver = solve_on_network!(ss, df, x_guess=x_dof, iteration_limit=1)
println(norm(x_dof - solver.solution))






# println(solver_return.status)
# println(ref(ss, :node, 1, "withdrawal") * nominal_values(ss, :mass_flow))






 

