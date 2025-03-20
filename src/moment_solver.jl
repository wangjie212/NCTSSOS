# C: is the varaibles commuting
# T: type of the coefficients
# monomap: map from monomials in DynamicPolynomials to variables in JuMP
struct MomentProblem{C,T,CR<:ConstraintRef} <: OptimizationProblem
    model::GenericModel{T}
    constraints::Vector{CR}
    monomap::Dict{Monomial{C},GenericVariableRef{T}}  # TODO: maybe refactor.
end

function substitute_variables(poly::Polynomial{C,T}, monomap::Dict{Monomial{C},GenericVariableRef{T}}) where {C,T}
    mapreduce(x -> (coefficient(x) * monomap[monomial(x)]), +, terms(poly))
end

# order: order of the moment problem
# clique_alg: algorithm for clique decomposition
# cliques_sub_mtx_col_basis: each clique, each obj/constraint, each ts_clique, each basis needed to index moment matrix
function moment_relax(pop::PolyOpt{C,T}, cliques_cons::Vector{Vector{Int}}, cliques_term_sparsities::Vector{Vector{TermSparsity{C}}}) where {C,T}

    # NOTE: objective and constraints may have integer coefficients, but popular JuMP solvers does not support integer coefficients
    # left type here to support BigFloat model for higher precision
    model = GenericModel{T}()

    # the union of clique_total_basis
    total_basis = sorted_union(map(zip(cliques_cons, cliques_term_sparsities)) do (cons_idx, term_sparsities)
        union(vec(reduce(vcat, [
            map(monomials(poly)) do m
                neat_dot(rol_idx, m * col_idx)
            end
            for (poly, term_sparsity) in zip([one(pop.objective); pop.constraints[cons_idx]], term_sparsities) for basis in term_sparsity.block_bases for rol_idx in basis for col_idx in basis
        ])))
    end...)

    # map the monomials to JuMP variables, the first variable must be 1
    @variable(model, y[1:length(total_basis)])
    @constraint(model, first(y) == 1)
    monomap = Dict(zip(total_basis, y))

    constraint_matrices =
        mapreduce(vcat, zip(cliques_term_sparsities, cliques_cons)) do (term_sparsities, cons_idx)
            mapreduce(vcat, zip(term_sparsities, [one(pop.objective), pop.constraints[cons_idx]...])) do (term_sparsity, poly)
                map(term_sparsity.block_bases) do ts_sub_basis
                    constrain_moment_matrix!(
                        model,
                        poly,
                        ts_sub_basis,
                        monomap,
                    )
                end
            end
        end

    @objective(model, Min, substitute_variables(symmetric_canonicalize(pop.objective), monomap))

    return MomentProblem(model, constraint_matrices, monomap)
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
