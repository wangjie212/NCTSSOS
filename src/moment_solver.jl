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
    # mapreduce(x -> (haskey(monomap, monomial(x)) ? (coefficient(x) * monomap[monomial(x)]) : (haskey(monomap, star(monomial(x))) ? coefficient(x) * monomap[star(monomial(x))] : AffExpr(zero(T)))), +, terms(poly))
    mapreduce(x -> (coefficient(x) * monomap[monomial(x)]), +, terms(poly))
end

# outputs: a vector of polyvars
function clique_decomp(variables::Vector{PolyVar{C}}, objective::Polynomial{C,T}, cons::Vector{Polynomial{C,T}}, clique_alg::EliminationAlgorithm, order::Int) where {C,T}
    label, tree = cliquetree(get_correlative_graph(variables, objective, cons, order), alg=clique_alg)
    return map(x -> sort(variables[label[x]]), collect(Vector{Int}, tree))
end

# TODO: it is possible to merge the logic of clique_decomp and block_decomp, is using Vector{Union{PolyVar{C},Monomial{C}}} ok?
# outputs: a vector of monomials
function block_decomp(G::SimpleGraph, basis::Vector{Monomial{C}}, clique_alg::EliminationAlgorithm) where {C}
    label, tree = cliquetree(G, alg=clique_alg)
    return map(x -> sort(basis[label[x]]), collect(Vector{Int}, tree))
end


# clique_alg: algorithm for clique decomposition
# cliques_sub_mtx_col_basis: each clique, each obj/constraint, each ts_clique, each basis needed to index moment matrix
function moment_relax(pop::PolynomialOptimizationProblem{C,T}, order::Int, cliques::Vector{Vector{PolyVar{C}}},
    cliques_cons::Vector{Vector{Int}}, cliques_total_basis::Vector{Vector{Monomial{C}}}, cliques_sub_mtx_col_basis::Vector{Vector{Vector{Vector{Monomial{C}}}}}) where {C,T}
    objective = symmetric_canonicalize(pop.objective)

    # NOTE: objective and constraints may have integer coefficients, but popular JuMP solvers does not support integer coefficients
    # left type here to support BigFloat model for higher precision
    model = GenericModel{T}()

    # the union of clique_total_basis
    total_basis = sorted_union(sorted_union.(cliques_total_basis)...)

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
        end for (i, clique_vars) in enumerate(cliques) for (ts_cliques, poly) in zip(cliques_sub_mtx_col_basis[i], [one(objective), pop.constraints[cliques_cons[i]]...])
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

# ordered_vars: variables in the order to be appeared in graph
# polys: objective + constraints order is important
# order: order of the moment problem
function get_correlative_graph(ordered_vars::Vector{PolyVar{C}}, obj::Polynomial{C,T}, cons::Vector{Polynomial{C,T}}, order::Int) where {C,T}
    nvars = length(ordered_vars)
    G = SimpleGraph(nvars)

    # find index of all unique variables in polynomial/monomial p
    vmap(p) = map(v -> findfirst(==(v), ordered_vars), unique!(effective_variables(p)))

    map(mono -> add_clique!(G, vmap(mono)), monomials(obj))

    for poly in cons
        # for clearer logic, I didn't combine the two branches
        if order == ceil(Int, maxdegree(poly) // 2)
            # if objective or order too large, each term forms a clique
            map(mono -> add_clique!(G, vmap(mono)), monomials(poly))
        else
            # NOTE: if is constraint and order not too large, all variables in the constraint forms a clique
            # this ensures each "small" constraint is in a clique ?
            add_clique!(G, vmap(poly))
        end
    end
    return G
end

#   porting nccpop.jl's  get_graph
#   constructs the graph according to (7.5) and (7.14) together
#   total_support: support of objective, constraint and their corresponding term sparsity graph in previous iteration (7.14)
#   total_basis: basis used to index the moment matrix
function get_term_sparsity_graph(cons_support::Vector{Monomial{C}}, total_support::Vector{Monomial{C}}, total_basis::Vector{Monomial{C}}) where {C}
    nterms = length(total_basis)
    G = SimpleGraph(nterms)
    for i in 1:nterms, j in i+1:nterms
        map(c_supp -> neat_dot(total_basis[i], c_supp * total_basis[j]) in total_support && add_edge!(G, i, j), cons_support)
    end
    return G
end

#   supp(G,g): monomials that are either v^† g_supp v where v is a vertex in G, or β^† g_supp γ where {β,γ} is an edge in G following (10,4)
#   given term sparsity graph G, which terms needs to be considered as a variable for describing the localizing/moment matrix with respect to g
function term_sparsity_graph_supp(G::SimpleGraph, basis::Vector{Monomial{C}}, g::Polynomial) where {C}
    # following (10.4) in Sparse Polynomial Optimization: Theory and Practise
    mapreducegsupp(a, b) = (mapreduce(g_supp -> neat_dot(a, g_supp * b), vcat, monomials(g)))
    return union([mapreducegsupp(basis[v], basis[v]) for v in vertices(G)], [mapreducegsupp(basis[e.src], basis[e.dst]) for e in edges(G)])
end


# total_basis: basis of the moment matrix for a clique
# prev_localizing_mtx_basis: support of previous iteration of term sparse graph
#
# localizing_mtx_basis: Vector of Vector of Monomial, each is the basis of the localizing matrix
function iterate_term_sparse_basis(total_basis::Vector{Monomial{C}}, prev_localizing_mtx_basis::Vector{Vector{Monomial{C}}}, elim_alg::EliminationAlgorithm) where {C}

    Fs = [SimpleGraph(length(total_basis))]

    Gs = block_decomp.(Fs, Ref(total_basis), Ref(elim_alg))

    prev_localizing_mtx_basis == localizing_mtx_basis && @info "Term Sparse Graph has stabilized"
    return localizing_mtx_basis
end

# calling nccpop.jl's get_blocks
# tsupp : total support for objective and constraint
# basis : vector of total_basis for each clique's moment matrix
function get_blocks(clique::Vector{PolyVar{C}}, clique_cons::Vector{Int}, obj::Polynomial{C,T}, cons::Vector{Polynomial{C,T}}) where {C,T}

end

# function support_extension(G::SimpleGraph, basis::Vector{Monomial{C}}, g::Polynomial) where {C}
#     G_nxt = SimpleGraph(nv(G))
#     G_supp = term_sparsity_graph_supp(G, basis, g)

#     return G_nxt
# end

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
