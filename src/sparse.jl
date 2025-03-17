# cliques: grouping of variables into sets of variables
# cliques_cons: groups constraints according to cliques,
# constraints in a clique only has support on corresponding variables
# discarded_cons: constraints that are not in any clique
# cliques_idcs_bases: within each clique, the vectors of monomials used to index moment/localizing matrices
struct CorrelativeSparsity{C}
    cliques::Vector{Vector{PolyVar{C}}}
    cliques_cons::Vector{Vector{Int}}
    discarded_cons::Vector{Int}
    cliques_idcs_bases::Vector{Vector{Vector{Monomial{C}}}}
end

# ordered_vars: variables in the order to be appeared in graph
# polys: objective + constraints order is important
# order: order of the moment problem
function get_correlative_graph(ordered_vars::Vector{PolyVar{C}}, obj::Polynomial{C,T}, cons::Vector{Polynomial{C,T}}, order::Int) where {C,T}
    # NOTE: Ordering in DynamicPolynomials is funky
    @assert issorted(ordered_vars, rev=true) "Variables must be sorted"

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

function clique_decomp(G::SimpleGraph, clique_alg::EliminationAlgorithm)
    label, tree = cliquetree(G, alg=clique_alg)
    return map(x -> label[x], collect(Vector{Int}, tree))
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

function correlative_sparsity(variables::Vector{PolyVar{C}}, objective::Polynomial{C,T}, cons::Vector{Polynomial{C,T}}, order::Int, elim_algo::Union{Nothing,EliminationAlgorithm}) where {C,T}
    cliques = isnothing(elim_algo) ? [variables] : map(x -> variables[x], clique_decomp(get_correlative_graph(variables, objective, cons, order), elim_algo))

    cliques_cons, discarded_cons = assign_constraint(cliques, cons)

    # get the operators needed to index columns of moment/localizing mtx in each clique
    # depending on the clique's varaibles each is slightly different
    cliques_idx_basis = map(zip(cliques, cliques_cons)) do (clique, clique_cons)
        # get the basis of the moment matrix in a clique, then sort it
        [[get_basis(clique, order)]; get_basis.(Ref(clique), order .- ceil.(Int, maxdegree.(cons[clique_cons]) / 2))]
    end

    return CorrelativeSparsity{C}(cliques, cliques_cons, discarded_cons, cliques_idx_basis)
end


# term_sparse_graph_supp: support of the current term sparsity graph for an obj/cons
# block_bases: the bases of the moment/localizing matrix in each clique of term sparse graph
struct TermSparsity{C}
    term_sparse_graph_supp::Vector{Monomial{C}}
    block_bases::Vector{Vector{Monomial{C}}}
end

# porting nccpop.jl's  get_graph
# constructs the graph according to (7.5) and (7.14) together
# activated_supp: support of objective, constraint and their corresponding term sparsity graph in previous iteration (7.14)
# basis: basis used to index the moment matrix
function get_term_sparsity_graph(cons_support::Vector{Monomial{C}}, activated_supp::Vector{Monomial{C}}, basis::Vector{Monomial{C}}) where {C}
    nterms = length(basis)
    G = SimpleGraph(nterms)
    for i in 1:nterms, j in i+1:nterms
        for supp in cons_support
            if symmetric_canonicalize(neat_dot(basis[i], supp * basis[j])) in activated_supp
                add_edge!(G, i, j)
                continue
            end
        end
    end
    return G
end

# supp(G,g): monomials that are either v^† g_supp v where v is a vertex in G, or β^† g_supp γ where {β,γ} is an edge in G following (10,4)
# given term sparsity graph G, which terms needs to be considered as a variable for describing the localizing/moment matrix with respect to g
function term_sparsity_graph_supp(G::SimpleGraph, basis::Vector{Monomial{C}}, g::Polynomial) where {C}
    # following (10.4) in Sparse Polynomial Optimization: Theory and Practise
    # NOTE: Do I need to symmetric canonicalize it?
    gsupp(a, b) = map(g_supp -> neat_dot(a, g_supp * b), monomials(g))
    return union([gsupp(basis[v], basis[v]) for v in vertices(G)]..., [gsupp(basis[e.src], basis[e.dst]) for e in edges(G)]...)
end


# returns: F (the chordal graph), blocks in basis
function iterate_term_sparse_supp(activated_supp::Vector{Monomial{C}}, poly::Polynomial, basis::Vector{Monomial{C}}, elim_algo::EliminationAlgorithm) where {C}
    F = get_term_sparsity_graph(collect(monomials(poly)), activated_supp, basis)
    blocks = clique_decomp(F, elim_algo)
    map(block -> add_clique!(F, block), blocks)
    return TermSparsity(term_sparsity_graph_supp(F, basis, poly), map(x -> basis[x], blocks))
end
