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

# outputs: a vector of polyvars
function clique_decomp(pop::PolynomialOptimizationProblem, clique_alg::EliminationAlgorithm, order::Int)
    return map(x -> pop.variables[x], collect(Vector{Int}, cliquetree(get_correlative_graph(pop.variables, [pop.objective, pop.constraints...], order), alg=clique_alg)[2]))
end

# clique_alg: algorithm for clique decomposition
function moment_relax(pop::PolynomialOptimizationProblem{C,T}, order::Int, cliques::Vector{Vector{PolyVar{C}}}) where {C,T}
    objective = symmetric_canonicalize(pop.objective)

    clique_cons, ignored_cons = assign_constraint(cliques, pop.constraints)

    # NOTE: objective and constraints may have integer coefficients, but popular JuMP solvers does not support integer coefficients
    # left type here to support BigFloat model for higher precision
    model = GenericModel{T}()

    # clique_total_basis that is the union of support of all constraints and objective and moment matrix and localizing matrices in a clique
    clique_total_basis = map(cliques) do clique
        # get the basis of the moment matrix in a clique, then sort it
        clique_moment_basis = get_basis(clique, order)
        sorted_unique([neat_dot(row_idx, col_idx) for row_idx in clique_moment_basis for col_idx in clique_moment_basis])
    end

    # the union of clique_total_basis
    total_basis = sorted_union(clique_total_basis...)

    # map the monomials to JuMP variables, the first variable must be 1
    @variable(model, y[1:length(total_basis)])
    @constraint(model, first(y) == 1)
    monomap = Dict(zip(total_basis, y))

    # add the constraints
    constraint_matrices = vec([
        constrain_moment_matrix!(
            model,
            poly,
            get_basis(clique_vars, order - ceil(Int, maxdegree(poly) / 2)),
            monomap,
        ) for (i, clique_vars) in enumerate(cliques) for poly in [one(objective), pop.constraints[clique_cons[i]]...]
    ])

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
        substitute_variables(sum([coef * neat_dot(row_idx, mono * col_idx) for (coef, mono) in zip(coefficients(poly), monomials(poly))]), monomap) for
        row_idx in local_basis, col_idx in local_basis
    ]
    return @constraint(model, moment_mtx in PSDCone())
end

# ordered_vars: variables in the order to be appeared in graph
# polys: objective + constraints order is important
# order: order of the moment problem
function get_correlative_graph(ordered_vars::Vector{PolyVar{C}}, obj::Polynomial{C,T}, cons::Vector{Polynomial{C,T}}, order::Int) where {C,T}
    nvars = length(ordered_vars)
    G = SimpleGraph(nvars)

    add_clique!(G, map(v -> findfirst(==(v), ordered_vars), unique!(effective_variables(obj))))

    for poly in enumerate(cons)
        # for clearer logic, I didn't combine the two branches
        if order == ceil(Int, maxdegree(poly) // 2)
            # if objective or order too large, each term forms a clique
            map(monomials(poly)) do mono
                add_clique!(G, map(v -> findfirst(==(v), ordered_vars), unique!(effective_variables(mono))))
            end
        else
            # NOTE: if is constraint and order not too large, all variables in the constraint forms a clique
            # this ensures each "small" constraint is in a clique ?
            add_clique!(G, map(v -> findfirst(==(v), ordered_vars), unique!(effective_variables(poly))))
        end
    end
    return G
end

function assign_constraint(cliques::Vector{Vector{PolyVar{C}}}, cons::Vector{Polynomial{C,T}}) where {C,T}
    # assign each constraint to a clique
    # there might be constraints that are not captured by any single clique,
    # NOTE: we ignore this constraint. This should only occur at lower order of relaxation.

    # clique_cons: vector of vector of constraints index, each belong to a clique
    clique_cons = map(cliques) do clique
        findall(g -> issubset(unique!(effective_variables(g)), clique), cons)
    end
    return clique_cons, setdiff(1:length(cons), union(clique_cons...))
end
