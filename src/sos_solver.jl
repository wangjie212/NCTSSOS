struct SOSProblem{T} <: OptimizationProblem
    model::GenericModel{T}
end

# Decompose the matrix into the form sum_j C_αj * g_j
# j: index of the constraint
# α: the monomial (JuMP variable)
function get_Cαj(basis::Vector{GenericVariableRef{T}}, localizing_mtx::VectorConstraint{F,S,Shape}) where {T,F,S,Shape}
    dim = get_dim(localizing_mtx)
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
    dual_variables = map(constraint_object.(moment_problem.constraints)) do cons
        # FIXME: HEDIOUS
        G_dim = get_dim(cons)
        @variable(dual_model, [1:G_dim, 1:G_dim] in ((cons.set isa MOI.Zeros) ? SymmetricMatrixSpace() : PSDCone()))
    end

    # b: to bound the minimum value of the primal problem
    @variable(dual_model, b)
    @objective(dual_model, Max, b)

    primal_objective_terms = objective_function(moment_problem.model).terms

    # NOTE: objective is Symmetric, hence when comparing polynomials, we need to canonicalize them first
    # TODO: fix this for trace
    unsymmetrized_basis = sort(collect(keys(moment_problem.monomap)))

    symmetric_basis = sort(unique!([moment_problem.reduce_func(symmetric_canonicalize(basis)) for basis in unsymmetrized_basis]))

    # JuMP variables corresponding to symmetric_basis
    symmetric_variables = getindex.(Ref(moment_problem.monomap), symmetric_basis)

    # specify constraints
    fα_constraints = [AffExpr(get(primal_objective_terms, α, zero(T))) for α in symmetric_variables]

    symmetrized_α2cons_dict = Dict(zip(unsymmetrized_basis, map(x -> searchsortedfirst(symmetric_basis, moment_problem.reduce_func(symmetric_canonicalize(x))), unsymmetrized_basis)))

    unsymmetrized_basis_vals = getindex.(Ref(moment_problem.monomap), unsymmetrized_basis)

    add_to_expression!(fα_constraints[1], -1.0, b)

    for (i, sdp_constraint) in enumerate(moment_problem.constraints)
        Cαj = get_Cαj(unsymmetrized_basis_vals, constraint_object(sdp_constraint))
        for (k, α) in enumerate(unsymmetrized_basis)
            for jj in 1:size(Cαj[k], 2)
                for ii in nzrange(Cαj[k], jj)
                    add_to_expression!(fα_constraints[symmetrized_α2cons_dict[α]], -Cαj[k].nzval[ii], dual_variables[i][jj, Cαj[k].rowval[ii]])
                end
            end
        end
    end
    @constraint(dual_model, fα_constraints .== 0)

    return SOSProblem(dual_model)
end
