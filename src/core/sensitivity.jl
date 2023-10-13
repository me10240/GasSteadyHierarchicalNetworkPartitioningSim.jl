
function initialize_for_sensitivity_computations(ss::SteadySimulator)::AbstractArray

	ndofs = length(ref(ss, :dof))
    RHS = zeros(ndofs, ref(ss, :total_control_vars))
    return RHS 

end

function assemble_rhs_corresponding_to_single_control_var!(ss::SteadySimulator,  i::Int64,  RHS::AbstractArray, x_dof::AbstractArray)

	component, index, ctrl_type = ref(ss, :control_vars, i)

	if component == :node && ctrl_type == "injection"

		RHS[ref(ss, :node, index, "dof"), i] = -1.0  # want -dr/dq 

	elseif component == :node && ctrl_type == "pressure"
		
		RHS[ref(ss, :node, index, "dof"), i] = 1.0  # want -dr/dp_slack

	elseif component == :compressor && ctrl_type == "c_ratio"

        fr_node = ref(ss, :compressor, index, "fr_node")

        RHS[ref(ss, :compressor, index, "dof"), i] = -1 * x_dof[ref(ss, :node, fr_node, "dof")] # want -d(alpha *p)/d alpha

	end
	return
end

function assemble_rhs_with_all_control_vars!(ss::SteadySimulator,  RHS::AbstractArray, x_dof::AbstractArray)
	for i = 1 : ref(ss, :total_control_vars)
        assemble_rhs_corresponding_to_single_control_var!(ss, i,  RHS, x_dof)
    end
    return
end

function calculate_slack_injection_sensitivities!(ss::SteadySimulator, RHS::AbstractArray)

	slack_node_indices = Vector{Int64}()

	@inbounds for (id, -) in ref(ss, :node)

		if ref(ss, :node, id, "is_slack") == 1

			push!(slack_node_indices, id)
		end

	end


	for id in slack_node_indices

		out_edge = ref(ss, :outgoing_dofs, id)
	    in_edge = ref(ss, :incoming_dofs, id)

		for i = 1 : ref(ss, :total_control_vars)

	        r = 0.0
	        r -= sum(RHS[e, i] for e in out_edge; init=0.0) 
	        r += sum(RHS[e, i] for e in in_edge; init=0.0)
	        RHS[ref(ss, :node, id, "dof"), i] = -r# q + inflow -outflow = 0
		
		end
	end

	return
end

function calculate_sensitivities!(ss::SteadySimulator, J::AbstractArray, RHS::AbstractArray)

	lufact = lu(J)
	for i = 1 : ref(ss, :total_control_vars)
		x = lufact \ RHS[:, i]
		RHS[:, i] .= x
	end

	calculate_slack_injection_sensitivities!(ss, RHS)

	return
end


function convert_sensitivities_to_given_units(ss::SteadySimulator, RHS::AbstractArray)::AbstractArray

	ndofs = length(ref(ss, :dof))
	Sensitivity_matrix = zeros(ndofs, ref(ss, :total_control_vars))

	for i = 1: length(ref(ss, :dof))

		x_scaling = 0.0
		param_scaling = 0.0
		sym, local_id = ref(ss, :dof, i)

		if sym == :node

			if ref(ss, :node, local_id, "is_slack") == 1
				x_scaling = ss.nominal_values[:mass_flow]
			else
				x_scaling = (ref(ss, :is_pressure_node, local_id)) ? ss.nominal_values[:pressure] : get_potential(ss, ss.nominal_values[:pressure])
			end 

		elseif sym in [:pipe, :compressor]

			x_scaling = ss.nominal_values[:mass_flow]

		else

			@error "Not implemented anything other than nodal pressures/injections, pipe flows and compressor ratios"

		end

		for col = 1 : ref(ss, :total_control_vars)
			component, index, - = ref(ss, :control_vars, col)

			if component == :node 
				if ref(ss, :node, index, "is_slack") == 1
					param_scaling = (ref(ss, :is_pressure_node, local_id)) ? ss.nominal_values[:pressure] : get_potential(ss, ss.nominal_values[:pressure])
				else
					param_scaling = ss.nominal_values[:mass_flow]
				end 

			elseif component == :compressor
				param_scaling = 1.0
			end

			Sensitivity_matrix[i, col] = (x_scaling / param_scaling) * RHS[i, col]

		end
	end

	return Sensitivity_matrix

end