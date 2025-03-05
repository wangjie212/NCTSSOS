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
        row_idx, col_idx = mod1(i, dim), div(i, dim, RoundUp)
        for (α, coeff) in getfield(cur_expr, :terms)
            α_idx = findfirst(isequal(α), basis)
            isnothing(α_idx) && continue
			push!(Is[α_idx], row_idx)
			push!(Js[α_idx], col_idx)
			push!(Vs[α_idx], coeff * (row_idx == col_idx ? 1.0 : 2.0))
		end
	end

    return [
        sparse(Is[i], Js[i], Vs[i], dim, dim)
        for i in eachindex(basis)
    ]
end

init_dual_variables(model::Model,dual_model::Model)= [
		begin
			G_dim = sdp_constraint |> (x -> getfield(x, :side_dimension)) ∘ JuMP.shape ∘ constraint_object
			@variable(dual_model, [1:G_dim, 1:G_dim] in PSDCone())
		end
		for sdp_constraint in model[:mtx_constraints]
	]

function dualize(model::Model, total_basis2var_dict::Dict{M,VRef}) where {M<:AbstractMonomialLike,VRef<:AbstractVariableRef}

	dual_model = Model()	

    dual_variables = init_dual_variables(model,dual_model)

    dual_model[:dual_variables] = dual_variables

	@variable(dual_model, b)
    @objective(dual_model, Max, b)

    primal_objective_terms = model |> (x -> getfield(x, :terms)) ∘ objective_function

    symmetric_basis = sort(unique([symmetric_canonicalize(basis) for basis in keys(total_basis2var_dict)]))

	symmetric_variables = getindex.(Ref(total_basis2var_dict), symmetric_basis)

    f_α_constraints = [AffExpr(haskey(primal_objective_terms, α) ? primal_objective_terms[α] : 0.0) for α in symmetric_variables]

	f_α_constraints[1] -= b

    
    for (i, sdp_constraint) in enumerate(model[:mtx_constraints])
        C_α_j = get_C_α_j(symmetric_variables, constraint_object(sdp_constraint))
        f_α_constraints .-= (map(x -> LinearAlgebra.tr(x * dual_variables[i]), C_α_j))
    end

    @constraint(dual_model, f_α_constraints in MOI.Zeros(length(symmetric_basis)))
    return dual_model
end

