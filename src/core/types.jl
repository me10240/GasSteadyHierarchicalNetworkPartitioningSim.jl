

struct SteadySimulator
    data::Dict{String,Any}
    ref::Dict{Symbol,Any}
    sol::Dict{String,Any}
    nominal_values::Dict{Symbol,Any}
    params::Dict{Symbol,Any}
    initial_guess::Dict{Symbol,Any}
    boundary_conditions::Dict{Symbol,Any}
    potential_ratio_approx::Vector{Float64}
    pu_eos_coeffs::Function
    pu_pressure_to_pu_density::Function
    pu_density_to_pu_pressure::Function
end

ref(ss::SteadySimulator) = ss.ref
ref(ss::SteadySimulator, key::Symbol) = ss.ref[key]
ref(ss::SteadySimulator, key::Symbol, id::Int64) = ss.ref[key][id]
ref(ss::SteadySimulator, key::Symbol, id::Int64, field) = ss.ref[key][id][field]

params(ss::SteadySimulator) = ss.params
params(ss::SteadySimulator, key::Symbol) = ss.params[key]

nominal_values(ss::SteadySimulator) = ss.nominal_values
nominal_values(ss::SteadySimulator, key::Symbol) = ss.nominal_values[key]

initial_pipe_mass_flow(ss::SteadySimulator, id::Int64) = 
    ss.initial_guess[:pipe][id]

initial_compressor_flow(ss::SteadySimulator, id::Int64) = 
    ss.initial_guess[:compressor][id]

initial_nodal_pressure(ss::SteadySimulator, id::Int64) = 
    ss.initial_guess[:node][id]

initial_control_valve_flow(ss::SteadySimulator, id::Int64) = 
    ss.initial_guess[:control_valve][id]

initial_valve_flow(ss::SteadySimulator, id::Int64) = 
    ss.initial_guess[:valve][id]

initial_resistor_flow(ss::SteadySimulator, id::Int64) = 
    ss.initial_guess[:resistor][id]

initial_loss_resistor_flow(ss::SteadySimulator, id::Int64) = 
    ss.initial_guess[:loss_resistor][id]

initial_short_pipe_flow(ss::SteadySimulator, id::Int64) = 
    ss.initial_guess[:short_pipe][id]


get_eos_coeffs(ss::SteadySimulator) = ss.pu_eos_coeffs(nominal_values(ss), params(ss))
get_pressure(ss::SteadySimulator, density) = ss.pu_density_to_pu_pressure(density, nominal_values(ss), params(ss))
get_density(ss::SteadySimulator, pressure) = ss.pu_pressure_to_pu_density(pressure, nominal_values(ss), params(ss))

function get_potential(ss::SteadySimulator, pressure)
    b1, b2 = get_eos_coeffs(ss)
    return (b1/2) * pressure^2 + (b2/3) * pressure^3
end 

function get_potential_derivative(ss::SteadySimulator, pressure) 
    b1, b2 = get_eos_coeffs(ss)
    return b1 * pressure + b2 * pressure^2
end 


@enum CONTROL_TYPE begin
    c_ratio_control = 0
    discharge_pressure_control = 1
    flow_control = 2
    pressure_control = 3
    unknown_control = 100
end

@enum SOLVER_STATUS begin 
    unique_physical_solution = 0 
    nl_solve_failure = 1 
    unique_unphysical_solution = 2 
    unphysical_solution = 3
end

struct SolverReturn 
    status::SOLVER_STATUS
    iterations::Int 
    residual_norm::Float64 
    time::Float64 
    solution::Vector{Float64}
    negative_flow_in_compressors::Vector{Int64}
    negative_potentials_in_nodes::Vector{Int64}
    pressure_domain_not_satisfied_in_nodes::Vector{Int64}
end 




