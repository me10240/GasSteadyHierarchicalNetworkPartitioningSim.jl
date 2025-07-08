module GasSteadyHierarchicalNetworkPartitioningSim
# import GasSteadySim as GSS
import GasSteadySim: initialize_simulator, SteadySimulator, update_solution_fields_in_ref!, populate_solution!, get_potential, invert_positive_potential, run_simulator!, assemble_mat!, assemble_residual!
import JSON
using NLsolve
using SparseArrays
using LineSearches
using LinearAlgebra

# include("io/json.jl")
# include("io/data_utils.jl")

# include("unit_conversion/unit_convertor_utils.jl")
# include("unit_conversion/to_si.jl")
# include("unit_conversion/to_english.jl")
# include("unit_conversion/to_pu.jl")
# include("unit_conversion/unit_convertors.jl")

# include("core/eos.jl")
include("core/types.jl")
include("core/ref.jl")
include("core/ig.jl")
# include("core/bc.jl")
# include("core/sol.jl")
include("core/initialize_ss.jl")
# include("core/assemble.jl")
include("core/run_ss_new.jl")  # new
include("core/partition.jl") # new
# include("core/output.jl")
# include("io/writer.jl")
include("core/export.jl")

end # module
