

function _check_correctness(sol::Dict{String,Any}, exact_sol::Dict{String,Any})

    entries = ["compressor_flow", "pipe_flow", "nodal_pressure"]  # can add more quantities here if required
    for key in entries
        for i in keys(sol[key])
            # @show sol[key][i] - exact_sol[key][string(i)]
            @test isapprox(sol[key][i], exact_sol[key][string(i)]; rtol = 1e-4)
        end
    end
    return
end

""" helper function to check all residuals (other than compressors, control_valves) are zero """ 
function _check_residuals(ss::SteadySimulator, r::Vector{Float64}, tol = 1e-2)
    components = [:node, :pipe, :valve, :resistor]
    for component in components 
        for (_, comp) in get(ref(ss), component, [])
            eqn_no = comp["dof"]
            @test r[eqn_no] ≈ 0.0 atol = tol 
        end 
    end 
end 