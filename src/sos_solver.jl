function get_C_α_j(basis::VE, localizing_mtx::VectorConstraint) where {VE<:AbstractVector{<:JuMP.AbstractVariableRef}}

	dim = localizing_mtx |> (x -> getfield(x, :side_dimension)) ∘ JuMP.shape 

    Is, Js, Vs = [[
        begin
            I = i <= 2 ? Int64[] : Float64[]
            sizehint!(I, max(5, dim))
            I
        end for _ in eachindex(basis)
    ] for i in 1:3]

	for (i,cur_expr) in enumerate(getfield(localizing_mtx, :func))
        row_idx, col_idx = mod1(i, dim), min(i ÷ dim + 1, dim)
        for (α, coeff) in getfield(cur_expr, :terms)
            α_idx = findfirst(isequal(α), basis)
			push!(Is[α_idx], row_idx)
			push!(Js[α_idx], col_idx)
			push!(Vs[α_idx], coeff)
		end
	end

    return [
        sparse(Is[i], Js[i], Vs[i], dim, dim)
        for i in eachindex(basis)
    ]
end

function dualize(model::Model)
	dual_model = Model()	


    dual_variables = [
        begin
            G_dim = sdp_constraint |> (x -> getfield(x, :side_dimension)) ∘ JuMP.shape ∘ constraint_object
            @variable(dual_model, [1:G_dim, 1:G_dim] in PSDCone())
        end for sdp_constraint in model[:mtx_constraints]
    ]

	@variable(dual_model, b)
    @objective(dual_model, Max, b)


    primal_objective_terms = model |> (x -> getfield(x, :terms)) ∘ objective_function


    f_α_constraints = [AffExpr(haskey(primal_objective_terms, α) ? primal_objective_terms[α] : 0.0) for α in all_variables(model)]

	f_α_constraints[1] -= b

    [
        begin
            C_α_j = get_C_α_j(all_variables(model), constraint_object(sdp_constraint));
            f_α_constraints .-= (map(x -> LinearAlgebra.tr(x * dual_variables[i]), C_α_j))
        end for (i, sdp_constraint) in enumerate(model[:mtx_constraints])
    ]

    @constraint(dual_model, f_α_constraints in MOI.Zeros(length(all_variables(model))))
    return dual_model
end

