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

function sos_dualize(moment_problem::MomentProblem{C,T}) where {C,T}
    dual_model = GenericModel{T}()

    # Initialize Gj as variables
    dual_variables = map(moment_problem.constraints) do cons
        G_dim = JuMP.shape(constraint_object(cons)).side_dimension
        @variable(dual_model, [1:G_dim, 1:G_dim] in PSDCone())
    end

    # b: to bound the minimum value of the primal problem
    @variable(dual_model, b)
    @objective(dual_model, Max, b)

    primal_objective_terms = objective_function(moment_problem.model).terms

    # NOTE: objective is Symmetric, hence when comparing polynomials, we need to canonicalize them first
    # TODO: fix this for trace
    symmetric_basis = sort(unique!([symmetric_canonicalize(basis) for basis in keys(moment_problem.monomap)]))

    # JuMP variables corresponding to symmetric_basis
    symmetric_variables = getindex.(Ref(moment_problem.monomap), symmetric_basis)

    # specify constraints
    fα_constraints = [AffExpr(get(primal_objective_terms, α, zero(T))) for α in symmetric_variables]
    fα_constraints[1] -= b   # constant term in the primal objective
    for (i, sdp_constraint) in enumerate(moment_problem.constraints)
        Cαj = get_Cαj(collect(values(moment_problem.monomap)), constraint_object(sdp_constraint))
        for (k, α) in enumerate(keys(moment_problem.monomap))
            l = findfirst(==(moment_problem.monomap[symmetric_canonicalize(α)]), symmetric_variables)
            fα_constraints[l] -= LinearAlgebra.tr(Cαj[k] * dual_variables[i])
        end
    end
    @constraint(dual_model, fα_constraints .== 0)

    return SOSProblem(dual_model)
end
