struct PolyOptResult{C,T}
    objective::T # support for high precision solution
    cliques_cons::CorrelativeSparsity{C}
    cliques_term_sparsities::Vector{Vector{TermSparsity{C}}}
    # should contain objective and moment matrix or gram matrix for manually checking what happened
end

# implement NoElimination on clique_decomp
struct NoElimination <: EliminationAlgorithm end

# FIXME: how do I properly create a CliqueTree for complete graph?
function cliquetree(graph, alg::NoElimination, snd::SupernodeType)
    return cliquetree(complete_graph(nv(graph)), BFS(), snd)
end

@kwdef struct SolverConfig
    optimizer
    mom_order::Int
    cs_algo::EliminationAlgorithm = NoElimination()
    ts_algo::EliminationAlgorithm = NoElimination()
end

# TODO: add user interface
# learn from OMEinsum.jl
# consider adding Solver Interface
# consider obtaining enough information on Moment matrix etc to check if problem solved correctly
# prev_ans::Union{Nothing,PolyOptResult{C,T}}=nothing
function cs_nctssos(pop::PolynomialOptimizationProblem{C,T}, solver_config::SolverConfig) where {C,T}
    @assert solver_config.mom_order >= maximum([ceil(Int, maxdegree(poly) / 2) for poly in [pop.objective; pop.constraints]]) "Moment Matrix Order Must be large enough to include all monomials in problem"

    corr_sparsity = correlative_sparsity(pop.variables, pop.objective, pop.constraints, solver_config.mom_order, solver_config.cs_algo)

    cliques_objective = [reduce(+, [issubset(effective_variables(mono), clique) ? coef * mono : zero(mono) for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))]) for clique in corr_sparsity.cliques]

    # prepare the support for each term sparse localizing moment
    initial_activated_supp = [sorted_union(symmetric_canonicalize.(monomials(obj_part)), mapreduce(a -> monomials(a), vcat, pop.constraints[cons_idx]; init=Monomial{C}[]), [neat_dot(b, b) for b in idcs_bases[1]])
                              for (obj_part, cons_idx, idcs_bases) in zip(cliques_objective, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)]

    cliques_term_sparsities = map(zip(initial_activated_supp, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)) do (activated_supp, cons_idx, idcs_bases)
        [iterate_term_sparse_supp(activated_supp, poly, basis, solver_config.ts_algo) for (poly, basis) in zip([one(pop.objective); pop.constraints[cons_idx]], idcs_bases)]
    end

    moment_problem = moment_relax(pop, corr_sparsity.cliques_cons, cliques_term_sparsities)
    sos_problem = sos_dualize(moment_problem)
    set_optimizer(sos_problem.model, solver_config.optimizer)

    optimize!(sos_problem.model)
    return PolyOptResult(objective_value(sos_problem.model), corr_sparsity, cliques_term_sparsities)
end

function cs_nctssos_higher(pop::PolynomialOptimizationProblem{C,T}, prev_res::PolyOptResult, solver_config::SolverConfig) where {C,T}
    cliques_term_sparsities = [sorted_union([poly_term_sparsity.term_sparse_graph_supp for poly_term_sparsity in term_sparsities])
                               for term_sparsities in prev_ans.cliques_term_sparsities]

    moment_problem = moment_relax(pop, prev_res.corr_sparsity.cliques_cons, cliques_term_sparsities)
    sos_problem = sos_dualize(moment_problem)
    set_optimizer(sos_problem.model, solver_config.optimizer)
    optimize!(sos_problem.model)
    return PolyOptResult(objective_value(sos_problem.model), prev_res.corr_sparsity, cliques_term_sparsities)
end
