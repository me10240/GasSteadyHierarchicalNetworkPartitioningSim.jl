

function _check_correctness(sol::Dict{String,Any}, exact_sol::Dict{String,Any})

    entries = ["compressor_flow", "pipe_flow", "nodal_pressure", "control_valve_flow", "valve_flow", "resistor_flow", "loss_resistor_flow", "short_pipe_flow"]  # can add more quantities here if required
    for key in entries
        if haskey(sol, key)
            for i in keys(sol[key])
                # @show sol[key][i] - exact_sol[key][string(i)]
                @test isapprox(sol[key][i], exact_sol[key][string(i)]; rtol = 5e-3)
            end
        end
    end
    return
end