struct PolyOptResult{C,T}
    objective::T # support for high precision solution
    cliques_cons::CorrelativeSparsity{C}
    cliques_term_sparsities::Vector{Vector{TermSparsity{C}}}
    # should contain objective and moment matrix or gram matrix for manually checking what happened
end

# implement NoElimination on clique_decomp
struct NoElimination <: EliminationAlgorithm end

struct SolverConfig
    optimizer
    mom_order::Int
    cs_algo::EliminationAlgorithm
    ts_algo::EliminationAlgorithm
end

# TODO: add user interface
# learn from OMEinsum.jl
# consider adding Solver Interface
# consider obtaining enough information on Moment matrix etc to check if problem solved correctly
# prev_ans::Union{Nothing,PolyOptResult{C,T}}=nothing
function cs_nctssos(pop::PolynomialOptimizationProblem{C,T}, solver_config::SolverConfig) where {C,T}
    @assert mom_order >= maximum([ceil(Int, maxdegree(poly) / 2) for poly in [pop.objective; pop.constraints]]) "Moment Matrix Order Must be large enough to include all monomials in problem"
    @assert (!isnothing(prev_ans) && !isnothing(ts_algo)) "No Chordal Decomposition Algorithm provided for next iteration of Term Sparsity"

    isnothing(prev_ans) ? corr_sparsity = correlative_sparsity(pop.variables, pop.objective, pop.constraints, mom_order, cs_algo) : corr_sparsity = prev_ans.correlative_sparsity

    if isnothing(ts_algo)
        cliques_term_sparsities = [
            [TermSparsity(Vector{Monomial{C}}(), [basis]) for basis in idx_basis]
            for idx_basis in corr_sparsity.cliques_idcs_bases
        ]
    else
        cliques_objective = [reduce(+, [issubset(effective_variables(mono), clique) ? coef * mono : zero(mono) for (coef, mono) in zip(coefficients(pop.objective), monomials(pop.objective))]) for clique in corr_sparsity.cliques]

        # prepare the support for each term sparse localizing moment
        initial_activated_supp =  if isnothing(prev_ans)
            [sorted_union([poly_term_sparsity.term_sparse_graph_supp for poly_term_sparsity in term_sparsities])
                for term_sparsities in prev_ans.cliques_term_sparsities]
        else
            [sorted_union(symmetric_canonicalize.(monomials(obj_part)), mapreduce(a -> monomials(a), vcat, pop.constraints[cons_idx];init = Monomial{C}[]), [neat_dot(b, b) for b in idcs_bases[1]])
                for (obj_part, cons_idx, idcs_bases) in zip(cliques_objective, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)]
        end

        cliques_term_sparsities = map(zip(initial_activated_supp, corr_sparsity.cliques_cons, corr_sparsity.cliques_idcs_bases)) do (activated_supp, cons_idx, idcs_bases)
            [iterate_term_sparse_supp(activated_supp, poly, basis, ts_algo) for (poly, basis) in zip([one(pop.objective); pop.constraints[cons_idx]], idcs_bases)]
        end
    end

    moment_problem = moment_relax(pop, corr_sparsity.cliques_cons, cliques_term_sparsities)
    sos_problem = sos_dualize(moment_problem)
    set_optimizer(sos_problem.model, optimizer)
    optimize!(sos_problem.model)
    return PolyOptResult(objective_value(sos_problem.model), corr_sparsity, cliques_term_sparsities)
end
