struct SOSProblem{T} <: OptimizationProblem
    model::GenericModel{T}
end

# Decompose the matrix into the form sum_j C_αj * g_j
# j: index of the constraint
# α: the monomial (JuMP variable)
function get_Cαj(basis::Vector{GenericVariableRef{T}}, localizing_mtx::VectorConstraint{F,S,Shape}) where {T,F,S,Shape}
    dim = JuMP.shape(localizing_mtx).side_dimension
    cis = CartesianIndices((dim, dim))
    nbasis = length(basis)

    # Is, Js, Vs: for storing sparse repr. of C_αj,
    # each element corresponds to a monomial in basis.
    Is, Js, Vs = [Int[] for _ in 1:nbasis], [Int[] for _ in 1:nbasis], [T[] for _ in 1:nbasis]

    for (ci, cur_expr) in zip(cis, localizing_mtx.func)
        for (α, coeff) in cur_expr.terms
            α_idx = findfirst(==(α), basis)
            push!(Is[α_idx], ci.I[1])
            push!(Js[α_idx], ci.I[2])
            push!(Vs[α_idx], coeff)
        end
    end

    return [sparse(Is[i], Js[i], Vs[i], dim, dim) for i in eachindex(basis)]
end

function sos_dualize(
   moment_problem::MomentProblem{C,T} 
) where {C,T}
    dual_model = GenericModel{T}()

    dual_variables = [
        begin
            G_dim = sdp_constraint |> ((x -> getfield(x, :side_dimension)) ∘ JuMP.shape ∘ constraint_object)
            @variable(dual_model, [1:G_dim, 1:G_dim] in PSDCone())
        end for sdp_constraint in moment_problem.model[:mtx_constraints]
    ]

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
        C_α_j = get_Cαj(
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
