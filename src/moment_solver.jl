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
    mapreduce(x -> (coefficient(x) * monomap[monomial(x)]), +, terms(poly))
end

# clique_alg: algorithm for clique decomposition
# cliques_sub_mtx_col_basis: each clique, each obj/constraint, each ts_clique, each basis needed to index moment matrix
function moment_relax(pop::PolynomialOptimizationProblem{C,T}, order::Int, cliques::Vector{Vector{PolyVar{C}}},
    cliques_cons::Vector{Vector{Int}}, cliques_mtcs_bases::Vector{Vector{Vector{Vector{Monomial{C}}}}}) where {C,T}
    objective = symmetric_canonicalize(pop.objective)

    # NOTE: objective and constraints may have integer coefficients, but popular JuMP solvers does not support integer coefficients
    # left type here to support BigFloat model for higher precision
    model = GenericModel{T}()

    # the union of clique_total_basis
    total_basis = sorted_union(map(zip(cliques_cons, cliques_mtcs_bases)) do (cons, mtcs_bases)
        union(vec(reduce(vcat, [
            map(monomials(poly)) do m
                neat_dot(rol_idx, m * col_idx)
            end
            for (poly, bases) in zip([one(objective); pop.constraints[cons]], mtcs_bases) for basis in bases for rol_idx in basis for col_idx in basis
        ])))
    end...)


    # map the monomials to JuMP variables, the first variable must be 1
    @variable(model, y[1:length(total_basis)])
    @constraint(model, first(y) == 1)
    monomap = Dict(zip(total_basis, y))

    # add the constraints
    constraint_matrices = vec(reduce(vcat, [
        mapreduce(vcat, ts_cliques) do ts_sub_basis
            constrain_moment_matrix!(
                model,
                poly,
                ts_sub_basis,
                monomap,
            )
        end for (i, clique_vars) in enumerate(cliques) for (ts_cliques, poly) in zip(cliques_mtcs_bases[i], [one(objective), pop.constraints[cliques_cons[i]]...])
    ]))

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
        substitute_variables(sum([coef * neat_dot(row_idx, mono * col_idx) for (coef, mono) in zip(coefficients(poly), monomials(poly))]), monomap) for
        row_idx in local_basis, col_idx in local_basis
    ]
    return @constraint(model, moment_mtx in PSDCone())
end
