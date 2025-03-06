struct SOSProblem{T} <: OptimizationProblem
    model::GenericModel{T}
end


function get_C_α_j(
    basis::VE, localizing_mtx::VectorConstraint
) where {VE<:AbstractVector{<:JuMP.AbstractVariableRef}}
    dim = ((x -> getfield(x, :side_dimension)) ∘ JuMP.shape)(localizing_mtx)

    Is, Js, Vs = [
        [
            begin
                I = i <= 2 ? Int64[] : Float64[]
                sizehint!(I, max(5, dim))
                I
            end for _ in eachindex(basis)
        ] for i in 1:3
    ]

    for (i, cur_expr) in enumerate(getfield(localizing_mtx, :func))
        row_idx, col_idx = mod1(i, dim), div(i, dim, RoundUp)
        for (α, coeff) in getfield(cur_expr, :terms)
            α_idx = findfirst(isequal(α), basis)
            push!(Is[α_idx], row_idx)
            push!(Js[α_idx], col_idx)
            push!(Vs[α_idx], coeff)
        end
    end

    return [sparse(Is[i], Js[i], Vs[i], dim, dim) for i in eachindex(basis)]
end

function init_dual_variables!(model::GenericModel{T}, dual_model::GenericModel{T}) where {T}
    return [
        begin
            G_dim = sdp_constraint |> ((x -> getfield(x, :side_dimension)) ∘ JuMP.shape ∘ constraint_object)
            @variable(dual_model, [1:G_dim, 1:G_dim] in PSDCone())
        end for sdp_constraint in model[:mtx_constraints]
    ]
end

function sos_dualize(
   moment_problem::MomentProblem{C,T} 
) where {C,T}
    dual_model = GenericModel{T}()

    dual_variables = init_dual_variables!(moment_problem.model, dual_model)

    dual_model[:dual_variables] = dual_variables

    @variable(dual_model, b)
    @objective(dual_model, Max, b)

    primal_objective_terms = ((x -> getfield(x, :terms)) ∘ objective_function)(moment_problem.model)

    symmetric_basis = sort(
        unique([symmetric_canonicalize(basis) for basis in keys(moment_problem.monomap)])
    )

    symmetric_variables = getindex.(Ref(moment_problem.monomap), symmetric_basis)

    f_α_constraints = [
        AffExpr(haskey(primal_objective_terms, α) ? primal_objective_terms[α] : 0.0) for
        α in symmetric_variables
    ]

    f_α_constraints[1] -= b

    for (i, sdp_constraint) in enumerate(moment_problem.model[:mtx_constraints])
        C_α_j = get_C_α_j(
            getindex.(Ref(moment_problem.monomap), collect(keys(moment_problem.monomap))),
            constraint_object(sdp_constraint),
        )
        for (k, α) in enumerate(keys(moment_problem.monomap))
            l = findfirst(
                isequal(moment_problem.monomap[symmetric_canonicalize(α)]),
                symmetric_variables,
            )
            f_α_constraints[l] -= LinearAlgebra.tr(C_α_j[k] * dual_variables[i])
        end
    end

    @constraint(dual_model, f_α_constraints in MOI.Zeros(length(symmetric_basis)))
    return SOSProblem(dual_model)
end
