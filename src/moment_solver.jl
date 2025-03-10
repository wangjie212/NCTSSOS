# C: is the varaibles commuting
# T: type of the coefficients
# order: order of the moment problem
# monomap: map from monomials in DynamicPolynomials to variables in JuMP
struct MomentProblem{C,T,CR<:ConstraintRef} <: OptimizationProblem
    order::Int
    model::GenericModel{T}
    constraints::Vector{CR}
    monomap::Dict{Monomial{C},GenericVariableRef{T}}  # TODO: maybe refactor.
end

function substitute_variables(poly::Polynomial{C,T}, monomap::Dict{Monomial{C},GenericVariableRef{T}}) where {C,T}
    return mapreduce(x -> coefficient(x) * monomap[monomial(x)], +, terms(poly))
end

function moment_relax(pop::PolynomialOptimizationProblem{C,T}, order::Int) where {C,T}
    objective = symmetric_canonicalize(pop.objective)

    # NOTE: objective and constraints may have integer coefficients, but popular JuMP solvers does not support integer coefficients
    # left type here to support BigFloat model for higher precision
    model = GenericModel{T}()

    # get the basis of the moment matrix, then sort it
    moment_basis = get_basis(pop.variables, order)
    # total_basis that is the union of support of all constraints and objective and moment matrix and localizing matrices
    total_basis = sort(unique!([neat_dot(row_idx, col_idx) for row_idx in moment_basis for col_idx in moment_basis]))

    # map the monomials to JuMP variables, the first variable must be 1
    @variable(model, y[1:length(total_basis)])
    @constraint(model, first(y) == 1)
    monomap = Dict(zip(total_basis, y))

    # add the constraints
    constraint_matrices = [
        constrain_moment_matrix!(
            model,
            poly,
            get_basis(pop.variables, order - ceil(Int, maxdegree(poly) / 2)),
            monomap,
        ) for poly in [one(objective), pop.constraints...]
    ]

    # store the constraints for future reference

    @objective(model, Min, substitute_variables(objective, monomap))

    return MomentProblem(order, model, constraint_matrices, monomap)
end


function constrain_moment_matrix!(
    model::GenericModel{T},
    poly::Polynomial{C,T},
    local_basis::Vector{Monomial{C}},
    monomap::Dict{Monomial{C},GenericVariableRef{T}},
) where {C,T}
    moment_mtx = [
        substitute_variables(poly * neat_dot(row_idx, col_idx), monomap) for
        row_idx in local_basis, col_idx in local_basis
    ]
    return @constraint(model, moment_mtx in PSDCone())
end
