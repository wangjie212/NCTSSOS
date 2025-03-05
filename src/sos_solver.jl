function dualize(method::MomentMethod, model::Model)
	dual_model = Model()	


    variables = [
        begin
            G_dim = sdp_constraint |> (x -> getfield(x, :side_dimension)) ∘ JuMP.shape ∘ constraint_object
            @variable(dual_model, [1:G_dim, 1:G_dim] in PSDCone())
        end for sdp_constraint in model[:sdp_constraint]
    ]

	@variable(dual_model, b)
    @objective(dual_model, Max, b)

	total_basis = keys(get_total_basis2var_dict(method))

	
    return model
end

