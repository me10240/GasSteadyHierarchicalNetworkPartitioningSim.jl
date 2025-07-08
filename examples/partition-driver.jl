

using GasSteadyHierarchicalNetworkPartitioningSim
using JSON
using LinearAlgebra
using NLSolversBase
using PyCall

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
ss = initialize_simulator(file, eos=eos_var, initial_guess_filename="") 
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

x_dof, cond_number_array = run_partitioned_ss(filepath, ss, eos=eos_var, cond_number=false, show_trace_flag=false, iteration_limit=2000, method=:trust_region)
df = prepare_for_nonlin_solve!(ss_copy)
push!(cond_number_array, cond(gradient(df), 1))
println("Condition numbers: \n", cond_number_array)
var = value!(df, x_dof)
println(norm(var))
solver = solve_on_network!(ss_copy, df, x_guess=x_dof, method=:trust_region)
# solver = solve_on_network!(ss, df, x_guess=x_dof, iteration_limit=1)
println(solver.iterations, " ", solver.residual_norm)


# choose articulation point pt, say c (of max degree)
# find connected comps C_i  of SG\c, add c to interface_nodes
# C_i = C_i U c, arrange C_i in decreasing size. 
# if any C_i biconnected, remove and push to diff list
# check if remaining  #C_i < N
# if not look at subgraph induced by C_i and its articulation point

function choose_articulation_point
function create_partition(ss::SteadySimulator; 
    num_partitions=4, write_to_file = false, filepath = "partition-expt.json")::Dict{String, Any}

    V = ref(ss, :node) |> collect 
    node_map = Dict(V[i][1] => i for i in range(1, length(V)))
    E = [ref(ss, :pipe)..., ref(ss, :compressor)...] 
    g = SimpleGraph(length(V)) 

    for (i, edge) in E 
        fr = edge["fr_node"]
        to = edge["to_node"]
        add_edge!(g, node_map[fr], node_map[to])
    end 

    # ## for vertex separator: 
    # # parts = Metis.separator(g)

    # ## for partition: 
    # parts = Metis.partition(g, num_partitions) #, alg=:RECURSIVE)

    data = Dict{Any,Any}(
        "num_partitions" => 0, 
        "interface_nodes" => [], 
        "slack_nodes" => [i for (i, node) in ref(ss, :node) if node["is_slack"] == 1]
    )
    
    cut_edge_list = []
    for edge in edges(g) 
        src = edge.src 
        dst = edge.dst 
        (parts[src] == parts[dst]) && (continue)
        # record edges that are being broken
        push!(cut_edge_list, edge)    
    end

    conn_comps = test_partition_consistency(g, cut_edge_list)
    true_num_partitions = length(conn_comps)
    
    if  true_num_partitions == 1
        @error "Given vertex separators do not partition network!"
        return Dict{String, Any}()
    end

    partitions = Dict( i => [] for i = 1 : true_num_partitions) 

    for i in 1 : true_num_partitions
        partitions[i] = [ V[k][1] for k in conn_comps[i] ]
    end 

    vertex_list = []
    for edge in cut_edge_list
        src = edge.src
        dst = edge.dst
        push!(vertex_list, [src, dst])
    end
    interface_seq = collect(Iterators.product(vertex_list...))[:]

    selected_seq  = nothing
    for seq in interface_seq
        if allunique(seq) == false
            continue
        end
        num_edges = test_vertex_sequence(g, collect(seq))
        if num_edges == 0
            selected_seq = seq
            break
        end
    end
    if isa(selected_seq, Nothing)
        @error "No valid interface sequence found!"
        return Dict{String, Any}()
    end

    for i in selected_seq
        push!(data["interface_nodes"], V[i][1])
        for edge in cut_edge_list
            if i in [edge.src, edge.dst]
                v  = (i == edge.src) ? edge.dst : edge.src
                for k = 1 :  true_num_partitions
                    # if partition has other vertex of the edge, add i
                    (V[v][1] in partitions[k] && V[i][1] ∉ partitions[k]) && push!(partitions[k], V[i][1])
                end
            end
        end
    end
                
        
    for (p, partition) in partitions 
        data[string(p)] = partition
    end
    data["num_partitions"] = true_num_partitions 

    for i = 1 : true_num_partitions
        if length(data[string(i)]) ==  2
                @error "At least one partition  has only 2 nodes - invalid partition!"
                return Dict{String, Any}()
        end
    end

    if write_to_file == true
        open(filepath, "w") do f 
            JSON.print(f, data, 2)
        end
        @info "Partition data saved to $filepath"
    end 

    return data
end 


